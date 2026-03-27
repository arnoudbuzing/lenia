//! Lenia continuous cellular automaton — Rust core.
//!
//! Compiled as a cdylib for use from Wolfram Language via LibraryFunctionLoad.
//! Exports LibraryLink-compatible functions using the standard C ABI.

use num_complex::Complex64;
use rustfft::FftPlanner;
use std::os::raw::c_void;

// ---------------------------------------------------------------------------
// LibraryLink type definitions (matching WolframLibrary.h)
// ---------------------------------------------------------------------------

/// mint is ptrdiff_t (i64 on 64-bit platforms)
type MInt = i64;
/// mreal is double
type MReal = f64;
/// MTensor is an opaque pointer
type MTensor = *mut c_void;

/// MArgument is a union — on 64-bit it's 8 bytes (pointer-sized).
/// We treat it as a raw pointer since all variants are pointer types.
type MArgument = *mut c_void;

/// WolframLibraryData is a pointer to a struct of function pointers.
/// We define just enough of the struct to access MTensor callbacks.
///
/// IMPORTANT: The field offsets must match WolframLibrary.h exactly.
/// We use a raw byte array and manually compute offsets.
type WolframLibraryData = *mut c_void;

// Offsets into the st_WolframLibraryData struct for Wolfram 15.0
// Counted from the actual WolframLibrary.h  (each slot is 8 bytes / one pointer)
//   0: UTF8String_disown
//   1: MTensor_new
//   2: MTensor_free
//   3: MTensor_clone
//   4: MTensor_shareCount
//   5: MTensor_disown
//   6: MTensor_disownAll
//   7: MTensor_setInteger
//   8: MTensor_setReal
//   9: MTensor_setComplex
//  10: MTensor_setMTensor
//  11: MTensor_getInteger
//  12: MTensor_getReal
//  13: MTensor_getComplex
//  14: MTensor_getMTensor
//  15: MTensor_getRank
//  16: MTensor_getDimensions
//  17: MTensor_getType
//  18: MTensor_getFlattenedLength
//  19: MTensor_getIntegerData
//  20: MTensor_getRealData
//  21: MTensor_getComplexData

const MTENSOR_NEW_OFFSET: usize = 1;
const MTENSOR_GETRANK_OFFSET: usize = 15;
const MTENSOR_GETDIMENSIONS_OFFSET: usize = 16;
const MTENSOR_GETREALDATA_OFFSET: usize = 20;

// MType_Real = 3 (from WolframLibrary.h)
const MTYPE_REAL: MInt = 3;
// LIBRARY_NO_ERROR = 0
const LIBRARY_NO_ERROR: i32 = 0;
// WolframLibraryVersion = 8 (from WolframLibrary.h for v15)
const WOLFRAM_LIBRARY_VERSION: MInt = 8;

/// Read a function pointer from the WolframLibraryData struct at a given slot offset.
unsafe fn get_callback<T>(lib_data: WolframLibraryData, slot: usize) -> T {
    let base = lib_data as *const *const c_void;
    let ptr = base.add(slot);
    std::mem::transmute_copy(&*ptr)
}

// Callback function types
type MTensorNewFn = unsafe extern "C" fn(MInt, MInt, *const MInt, *mut MTensor) -> i32;
type MTensorGetRankFn = unsafe extern "C" fn(MTensor) -> MInt;
type MTensorGetDimensionsFn = unsafe extern "C" fn(MTensor) -> *const MInt;
type MTensorGetRealDataFn = unsafe extern "C" fn(MTensor) -> *mut MReal;

/// Read MTensor from MArgument (MArgument is union of pointers; tensor field = *MTensor)
unsafe fn margument_get_tensor(arg: MArgument) -> MTensor {
    // MArgument.tensor is a *MTensor, and for the tensor variant,
    // the union contains a pointer to the MTensor.
    // So: *(MArgument as *mut MTensor)
    *(arg as *mut MTensor)
}

/// Write MTensor to MArgument
unsafe fn margument_set_tensor(arg: MArgument, tensor: MTensor) {
    *(arg as *mut MTensor) = tensor;
}

/// Read MInt from MArgument
unsafe fn margument_get_integer(arg: MArgument) -> MInt {
    *(arg as *mut MInt)
}

/// Read MReal from MArgument
unsafe fn margument_get_real(arg: MArgument) -> MReal {
    *(arg as *mut MReal)
}

// ---------------------------------------------------------------------------
// 2D FFT helpers
// ---------------------------------------------------------------------------

fn fft2d(data: &mut [Complex64], rows: usize, cols: usize, inverse: bool) {
    let mut planner = FftPlanner::<f64>::new();

    // --- rows ---
    let row_plan = if inverse {
        planner.plan_fft_inverse(cols)
    } else {
        planner.plan_fft_forward(cols)
    };
    let mut scratch = vec![Complex64::new(0.0, 0.0); row_plan.get_inplace_scratch_len()];
    for i in 0..rows {
        let start = i * cols;
        row_plan.process_with_scratch(&mut data[start..start + cols], &mut scratch);
    }

    // --- columns ---
    let col_plan = if inverse {
        planner.plan_fft_inverse(rows)
    } else {
        planner.plan_fft_forward(rows)
    };
    let mut scratch = vec![Complex64::new(0.0, 0.0); col_plan.get_inplace_scratch_len()];
    let mut col_buf = vec![Complex64::new(0.0, 0.0); rows];
    for j in 0..cols {
        for i in 0..rows {
            col_buf[i] = data[i * cols + j];
        }
        col_plan.process_with_scratch(&mut col_buf, &mut scratch);
        for i in 0..rows {
            data[i * cols + j] = col_buf[i];
        }
    }
}

// ---------------------------------------------------------------------------
// Kernel construction
// ---------------------------------------------------------------------------

fn build_kernel(rows: usize, cols: usize, radius: usize) -> Vec<f64> {
    let mid_r = (rows + 1) / 2 - 1;
    let mid_c = (cols + 1) / 2 - 1;
    let r_inv = 1.0 / radius as f64;

    let n = rows * cols;
    let mut kernel = vec![0.0f64; n];
    let mut total = 0.0f64;

    for i in 0..rows {
        for j in 0..cols {
            let di = i as f64 - mid_r as f64;
            let dj = j as f64 - mid_c as f64;
            let r = (di * di + dj * dj).sqrt() * r_inv;
            if r > 0.0 && r < 1.0 {
                let val = (4.0 - 1.0 / (r * (1.0 - r))).exp();
                kernel[i * cols + j] = val;
                total += val;
            }
        }
    }

    if total > 0.0 {
        let inv_total = 1.0 / total;
        for v in kernel.iter_mut() {
            *v *= inv_total;
        }
    }

    // FFT-shift
    let shift_r = mid_r;
    let shift_c = mid_c;
    let mut shifted = vec![0.0f64; n];
    for i in 0..rows {
        for j in 0..cols {
            let si = (i + shift_r) % rows;
            let sj = (j + shift_c) % cols;
            shifted[i * cols + j] = kernel[si * cols + sj];
        }
    }

    shifted
}

// ---------------------------------------------------------------------------
// Growth function
// ---------------------------------------------------------------------------

#[inline]
fn growth(n: f64, mu: f64, sigma: f64) -> f64 {
    let d = n - mu;
    2.0 * (-d * d / (2.0 * sigma * sigma)).exp() - 1.0
}

// ---------------------------------------------------------------------------
// Single Lenia step
// ---------------------------------------------------------------------------

fn lenia_step_internal(
    grid: &mut [f64],
    kernel_fft: &[Complex64],
    rows: usize,
    cols: usize,
    mu: f64,
    sigma: f64,
    dt: f64,
    buf: &mut Vec<Complex64>,
) {
    let n = rows * cols;
    let norm = 1.0 / n as f64;

    for i in 0..n {
        buf[i] = Complex64::new(grid[i], 0.0);
    }

    fft2d(buf, rows, cols, false);

    for i in 0..n {
        buf[i] *= kernel_fft[i];
    }

    fft2d(buf, rows, cols, true);

    for i in 0..n {
        let potential = buf[i].re * norm;
        let g = growth(potential, mu, sigma);
        grid[i] = (grid[i] + dt * g).clamp(0.0, 1.0);
    }
}

// =========================================================================
// LibraryLink entry points
// =========================================================================

/// Required: return the library version.
#[no_mangle]
pub extern "C" fn WolframLibrary_getVersion() -> MInt {
    WOLFRAM_LIBRARY_VERSION
}

/// Required: library initialization.
#[no_mangle]
pub extern "C" fn WolframLibrary_initialize(_lib_data: WolframLibraryData) -> i32 {
    LIBRARY_NO_ERROR
}

/// lenia_simulate(grid, steps, radius, mu, sigma, dt, returnHistory) -> result
///
/// Wolfram Language call:
///   LibraryFunctionLoad[lib, "lenia_simulate",
///     {{Real, 2, "Constant"}, Integer, Integer, Real, Real, Real, Integer},
///     {Real, 2}]    (* non-history *)
///   -- or for history: returns {Real, 3}
///
/// Args[0]: MTensor grid (input, rank 2, Real)
/// Args[1]: Integer steps
/// Args[2]: Integer radius
/// Args[3]: Real mu
/// Args[4]: Real sigma
/// Args[5]: Real dt
/// Args[6]: Integer returnHistory (0 or 1)
/// Res: MTensor (rank 2 or rank 3)
#[no_mangle]
pub unsafe extern "C" fn lenia_simulate(
    lib_data: WolframLibraryData,
    _argc: MInt,
    args: *mut MArgument,
    res: MArgument,
) -> i32 {
    // Extract callbacks from WolframLibraryData
    let mt_new: MTensorNewFn = get_callback(lib_data, MTENSOR_NEW_OFFSET);
    let mt_get_rank: MTensorGetRankFn = get_callback(lib_data, MTENSOR_GETRANK_OFFSET);
    let mt_get_dims: MTensorGetDimensionsFn = get_callback(lib_data, MTENSOR_GETDIMENSIONS_OFFSET);
    let mt_get_real_data: MTensorGetRealDataFn = get_callback(lib_data, MTENSOR_GETREALDATA_OFFSET);

    // Parse arguments
    let grid_tensor = margument_get_tensor(*args.add(0));
    let steps = margument_get_integer(*args.add(1)) as usize;
    let radius = margument_get_integer(*args.add(2)) as usize;
    let mu = margument_get_real(*args.add(3));
    let sigma = margument_get_real(*args.add(4));
    let dt = margument_get_real(*args.add(5));
    let return_history = margument_get_integer(*args.add(6)) != 0;

    // Get grid dimensions
    let rank = mt_get_rank(grid_tensor);
    if rank != 2 {
        return 1; // LIBRARY_DIMENSION_ERROR
    }
    let dims_ptr = mt_get_dims(grid_tensor);
    let rows = *dims_ptr as usize;
    let cols = *dims_ptr.add(1) as usize;
    let n = rows * cols;

    // Get input data
    let grid_data = mt_get_real_data(grid_tensor);
    let mut grid = vec![0.0f64; n];
    std::ptr::copy_nonoverlapping(grid_data, grid.as_mut_ptr(), n);

    // Build kernel & FFT
    let kernel = build_kernel(rows, cols, radius);
    let mut kernel_fft: Vec<Complex64> = kernel.iter().map(|&v| Complex64::new(v, 0.0)).collect();
    fft2d(&mut kernel_fft, rows, cols, false);

    // Allocate result tensor
    let result_tensor: MTensor;
    let result_data: *mut MReal;

    if return_history {
        // Rank-3 tensor: (steps+1) x rows x cols
        let out_dims: [MInt; 3] = [(steps + 1) as MInt, rows as MInt, cols as MInt];
        let mut tensor: MTensor = std::ptr::null_mut();
        let err = mt_new(MTYPE_REAL, 3, out_dims.as_ptr(), &mut tensor);
        if err != 0 {
            return err;
        }
        result_tensor = tensor;
        result_data = mt_get_real_data(result_tensor);

        // Write initial state
        std::ptr::copy_nonoverlapping(grid.as_ptr(), result_data, n);
    } else {
        // Rank-2 tensor: rows x cols
        let out_dims: [MInt; 2] = [rows as MInt, cols as MInt];
        let mut tensor: MTensor = std::ptr::null_mut();
        let err = mt_new(MTYPE_REAL, 2, out_dims.as_ptr(), &mut tensor);
        if err != 0 {
            return err;
        }
        result_tensor = tensor;
        result_data = mt_get_real_data(result_tensor);
    }

    // Simulation loop
    let mut buf = vec![Complex64::new(0.0, 0.0); n];
    for step in 0..steps {
        lenia_step_internal(&mut grid, &kernel_fft, rows, cols, mu, sigma, dt, &mut buf);

        if return_history {
            let offset = (step + 1) * n;
            std::ptr::copy_nonoverlapping(grid.as_ptr(), result_data.add(offset), n);
        }
    }

    // Write final result (non-history mode)
    if !return_history {
        std::ptr::copy_nonoverlapping(grid.as_ptr(), result_data, n);
    }

    // Set result
    margument_set_tensor(res, result_tensor);

    LIBRARY_NO_ERROR
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kernel_normalisation() {
        let k = build_kernel(32, 32, 5);
        let sum: f64 = k.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "Kernel sum = {sum}, expected 1.0");
    }

    #[test]
    fn growth_at_mu() {
        let g = growth(0.15, 0.15, 0.015);
        assert!((g - 1.0).abs() < 1e-10, "Growth at mu = {g}, expected 1.0");
    }

    #[test]
    fn growth_far_from_mu() {
        let g = growth(1.0, 0.15, 0.015);
        assert!((g - (-1.0)).abs() < 1e-10, "Growth far = {g}, expected -1.0");
    }
}

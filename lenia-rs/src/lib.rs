//! Lenia continuous cellular automaton — Rust core.
//!
//! Compiled as a cdylib for use from Wolfram Language via LibraryFunctionLoad.
//! Exports LibraryLink-compatible functions using the standard C ABI.

use num_complex::Complex64;
use rustfft::{Fft, FftPlanner};
use std::os::raw::c_void;
use std::sync::Arc;

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

fn fft2d(
    data: &mut [Complex64],
    rows: usize,
    cols: usize,
    inverse: bool,
    planner: &mut FftPlanner<f64>,
) {
    if inverse {
        let row_plan = planner.plan_fft_inverse(cols);
        let mut scratch = vec![Complex64::new(0.0, 0.0); row_plan.get_inplace_scratch_len()];
        for i in 0..rows {
            let start = i * cols;
            row_plan.process_with_scratch(&mut data[start..start + cols], &mut scratch);
        }

        let col_plan = planner.plan_fft_inverse(rows);
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
        
        let scale = 1.0 / (rows * cols) as f64;
        for val in data.iter_mut() {
            *val *= scale;
        }
    } else {
        let row_plan = planner.plan_fft_forward(cols);
        let mut scratch = vec![Complex64::new(0.0, 0.0); row_plan.get_inplace_scratch_len()];
        for i in 0..rows {
            let start = i * cols;
            row_plan.process_with_scratch(&mut data[start..start + cols], &mut scratch);
        }

        let col_plan = planner.plan_fft_forward(rows);
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
// Lenia Step Implementation
// ---------------------------------------------------------------------------

fn lenia_step_internal(
    grid: &mut [f64],
    kernel_fft: &[Complex64],
    rows: usize,
    cols: usize,
    mu: f64,
    sigma: f64,
    dt: f64,
    buf: &mut [Complex64],
    row_plan: &Arc<dyn Fft<f64>>,
    col_plan: &Arc<dyn Fft<f64>>,
    inv_row_plan: &Arc<dyn Fft<f64>>,
    inv_col_plan: &Arc<dyn Fft<f64>>,
) {
    let n = rows * cols;
    let mut col_buf = vec![Complex64::new(0.0, 0.0); rows];

    // 1. Grid FFT
    for i in 0..n {
        buf[i] = Complex64::new(grid[i], 0.0);
    }
    
    // Manual reuse of buffers for rows/cols in lenia_step to avoid reallocation
    let mut scratch_row = vec![Complex64::new(0.0, 0.0); row_plan.get_inplace_scratch_len()];
    for i in 0..rows {
        row_plan.process_with_scratch(&mut buf[i * cols..(i + 1) * cols], &mut scratch_row);
    }

    let mut scratch_col = vec![Complex64::new(0.0, 0.0); col_plan.get_inplace_scratch_len()];
    for j in 0..cols {
        for i in 0..rows {
            col_buf[i] = buf[i * cols + j];
        }
        col_plan.process_with_scratch(&mut col_buf, &mut scratch_col);
        for i in 0..rows {
            buf[i * cols + j] = col_buf[i];
        }
    }

    // 2. Pointwise multiplication
    for i in 0..n {
        buf[i] *= kernel_fft[i];
    }

    // 3. Inverse FFT
    let mut inv_scratch_row = vec![Complex64::new(0.0, 0.0); inv_row_plan.get_inplace_scratch_len()];
    for i in 0..rows {
        inv_row_plan.process_with_scratch(&mut buf[i * cols..(i + 1) * cols], &mut inv_scratch_row);
    }

    let mut inv_scratch_col = vec![Complex64::new(0.0, 0.0); inv_col_plan.get_inplace_scratch_len()];
    for j in 0..cols {
        for i in 0..rows {
            col_buf[i] = buf[i * cols + j];
        }
        inv_col_plan.process_with_scratch(&mut col_buf, &mut inv_scratch_col);
        for i in 0..rows {
            buf[i * cols + j] = col_buf[i];
        }
    }

    // 4. Growth and Update
    let scale = 1.0 / n as f64;
    for i in 0..n {
        let potential = buf[i].re * scale;
        let growth = 2.0 * (-(potential - mu).powi(2) / (2.0 * sigma.powi(2))).exp() - 1.0;
        grid[i] = (grid[i] + dt * growth).clamp(0.0, 1.0);
    }
}

// =========================================================================
// LibraryLink entry points
// =========================================================================

// ---------------------------------------------------------------------------
// Simulation logic helper
// ---------------------------------------------------------------------------

unsafe fn run_lenia_core(
    lib_data: WolframLibraryData,
    grid_tensor: MTensor,
    kernel_tensor: MTensor,
    steps: usize,
    mu: f64,
    sigma: f64,
    dt: f64,
    return_history: bool,
) -> Result<MTensor, i32> {
    // Extract callbacks
    let mt_new: MTensorNewFn = get_callback(lib_data, MTENSOR_NEW_OFFSET);
    let mt_get_rank: MTensorGetRankFn = get_callback(lib_data, MTENSOR_GETRANK_OFFSET);
    let mt_get_dims: MTensorGetDimensionsFn = get_callback(lib_data, MTENSOR_GETDIMENSIONS_OFFSET);
    let mt_get_real_data: MTensorGetRealDataFn = get_callback(lib_data, MTENSOR_GETREALDATA_OFFSET);

    // Get grid dimensions
    let rank = mt_get_rank(grid_tensor);
    if rank != 2 { return Err(3); }
    let dims_ptr = mt_get_dims(grid_tensor);
    let rows = *dims_ptr as usize;
    let cols = *dims_ptr.add(1) as usize;
    let n = rows * cols;

    // Get grid data
    let grid_data = mt_get_real_data(grid_tensor);
    if grid_data.is_null() { return Err(5); }
    let mut grid = vec![0.0f64; n];
    std::ptr::copy_nonoverlapping(grid_data, grid.as_mut_ptr(), n);

    // Get kernel data and compute its FFT
    let kernel_data = mt_get_real_data(kernel_tensor);
    if kernel_data.is_null() { return Err(5); }
    let mut kernel_fft: Vec<Complex64> = (0..n)
        .map(|i| Complex64::new(*kernel_data.add(i), 0.0))
        .collect();
    
    // Efficiently compute kernel FFT
    let mut planner = FftPlanner::new();
    fft2d(&mut kernel_fft, rows, cols, false, &mut planner);

    // Pre-calculate FFT plans for the simulation steps
    let row_plan = planner.plan_fft_forward(cols);
    let col_plan = planner.plan_fft_forward(rows);
    let inv_row_plan = planner.plan_fft_inverse(cols);
    let inv_col_plan = planner.plan_fft_inverse(rows);

    // Allocate result
    let result_tensor: MTensor;
    let result_data: *mut MReal;

    if return_history {
        let out_dims: [MInt; 3] = [(steps + 1) as MInt, rows as MInt, cols as MInt];
        let mut tensor: MTensor = std::ptr::null_mut();
        let err = mt_new(MTYPE_REAL, 3, out_dims.as_ptr(), &mut tensor);
        if err != 0 { return Err(err); }
        result_tensor = tensor;
        result_data = mt_get_real_data(result_tensor);
        std::ptr::copy_nonoverlapping(grid.as_ptr(), result_data, n);
    } else {
        let out_dims: [MInt; 2] = [rows as MInt, cols as MInt];
        let mut tensor: MTensor = std::ptr::null_mut();
        let err = mt_new(MTYPE_REAL, 2, out_dims.as_ptr(), &mut tensor);
        if err != 0 { return Err(err); }
        result_tensor = tensor;
        result_data = mt_get_real_data(result_tensor);
    }

    // Loop
    let mut buf = vec![Complex64::new(0.0, 0.0); n];
    for step in 0..steps {
        lenia_step_internal(
            &mut grid, &kernel_fft, rows, cols, mu, sigma, dt, &mut buf,
            &row_plan, &col_plan, &inv_row_plan, &inv_col_plan
        );
        if return_history {
            std::ptr::copy_nonoverlapping(grid.as_ptr(), result_data.add((step + 1) * n), n);
        }
    }
    if !return_history {
        std::ptr::copy_nonoverlapping(grid.as_ptr(), result_data, n);
    }

    Ok(result_tensor)
}

#[no_mangle]
pub extern "C" fn WolframLibrary_getVersion() -> MInt {
    WOLFRAM_LIBRARY_VERSION
}

#[no_mangle]
pub extern "C" fn WolframLibrary_initialize(_lib_data: WolframLibraryData) -> i32 {
    LIBRARY_NO_ERROR
}

/// lenia_simulate(grid, kernel, steps, mu, sigma, dt) -> result {Real, 2}
#[no_mangle]
pub unsafe extern "C" fn lenia_simulate(
    lib_data: WolframLibraryData,
    _argc: MInt,
    args: *mut MArgument,
    res: MArgument,
) -> i32 {
    let grid = margument_get_tensor(*args.add(0));
    let kernel = margument_get_tensor(*args.add(1));
    let steps = margument_get_integer(*args.add(2)) as usize;
    let mu = margument_get_real(*args.add(3));
    let sigma = margument_get_real(*args.add(4));
    let dt = margument_get_real(*args.add(5));

    match run_lenia_core(lib_data, grid, kernel, steps, mu, sigma, dt, false) {
        Ok(t) => { margument_set_tensor(res, t); 0 },
        Err(e) => e
    }
}

/// lenia_simulate_history(grid, kernel, steps, mu, sigma, dt) -> result {Real, 3}
#[no_mangle]
pub unsafe extern "C" fn lenia_simulate_history(
    lib_data: WolframLibraryData,
    _argc: MInt,
    args: *mut MArgument,
    res: MArgument,
) -> i32 {
    let grid = margument_get_tensor(*args.add(0));
    let kernel = margument_get_tensor(*args.add(1));
    let steps = margument_get_integer(*args.add(2)) as usize;
    let mu = margument_get_real(*args.add(3));
    let sigma = margument_get_real(*args.add(4));
    let dt = margument_get_real(*args.add(5));

    match run_lenia_core(lib_data, grid, kernel, steps, mu, sigma, dt, true) {
        Ok(t) => { margument_set_tensor(res, t); 0 },
        Err(e) => e
    }
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

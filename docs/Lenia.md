# Lenia

`Lenia` is the main simulation function for the Lenia continuous cellular automaton. It evolves an initial grid of densities over a specified number of steps using a circular convolution with a normalized kernel and a growth function.

## Syntax

```wolfram
Lenia[grid, steps]
Lenia[grid, steps, options]
```

## Arguments

* **`grid`**: A 2D matrix (list of lists) or a `NumericArray` of real values in the range `[0.0, 1.0]`. This represents the initial state of the simulation.
* **`steps`**: An integer specifying the number of simulation iterations to perform.

## Options

| Option | Default | Description |
| :--- | :--- | :--- |
| `"Mu"` | `0.15` | The center of the bell-shaped growth function $G$. |
| `"Sigma"` | `0.015` | The width (standard deviation) of the growth function $G$. |
| `"DT"` | `0.1` | The time step for each iteration. Large values may cause instability. |
| `"Radius"` | `13` | The spatial radius of the kernel $K$, measured in grid cells. |
| `"ReturnHistory"` | `False` | If `True`, returns a `NumericArray` containing the grid at every step. If `False`, only returns the final grid. |
| `"LeniaKernel"` | `"Bump"` | Specifies the kernel profile. Supported types include `"Bump"`, `"GaussianRing"`, `"StepRing"`, `"Polynomial"`, `"SmoothLife"`, and `"Sharp"`. |
| `Method` | `"Wolfram"` | The simulation backend. Use `"Wolfram"` for the pure Wolfram Language implementation or `"Rust"` for the accelerated backend. |

## Relationship to Other Functions

`Lenia` acts as a high-level wrapper that orchestrates the following components:

* **[LeniaKernel](file:///Users/arnoudb/github/lenia/Lenia/Kernel/Lenia.wl#L104)**: Used to generate the spatial kernel before the simulation loop begins.
* **[LeniaKernelFFT](file:///Users/arnoudb/github/lenia/Lenia/Kernel/Lenia.wl#L116)**: A convenience function that returns the Fourier transform of the kernel.
* **[LeniaStep](file:///Users/arnoudb/github/lenia/Lenia/Kernel/Lenia.wl#L119)**: The core iteration function used when `Method -> "Wolfram"`. It performs the FFT-based convolution and applies the growth function.
* **[LeniaGrowth](file:///Users/arnoudb/github/lenia/Lenia/Kernel/Lenia.wl#L86)**: The mathematical growth function $G(n) = 2 \cdot \exp(-(n-\mu)^2 / (2\sigma^2)) - 1$.
* **[LeniaSeed](file:///Users/arnoudb/github/lenia/Lenia/Kernel/Lenia.wl#L70)**: Often used to generate the initial `grid` input for `Lenia`. By default, it returns a `NumericArray` for memory efficiency.
* **[LeniaBlob](file:///Users/arnoudb/github/lenia/Lenia/Kernel/Lenia.wl#L80)**: A helper function for generating localized density centers (blobs) on a grid.

### Simulation Flow

1.  **Initialize**: Validate inputs and extract options.
2.  **Kernel Generation**: Create the normalized kernel $K$ using `LeniaKernel` with the specified `"Radius"` and `"LeniaKernel"`.
3.  **Simulation Loop**:
    *   If `Method -> "Rust"`, the entire loop is offloaded to the external library for maximum performance.
    *   If `Method -> "Wolfram"`, the function uses `Nest` (or `NestList` if `"ReturnHistory" -> True`) to repeatedly apply `LeniaStep`.
4.  **Output**: Returns the resulting grid(s) as a `NumericArray` for efficient storage and visualization.

## Examples

### Basic Simulation

Run a 50-step simulation on a 100x100 grid:

```wolfram
grid = LeniaSeed[100, "Ring"];
result = Lenia[grid, 50];
```

### Capturing History

Return all intermediate steps to create an animation:

```wolfram
history = Lenia[grid, 100, "ReturnHistory" -> True];
```

### Using the Rust Backend

For large grids or many steps, the Rust backend is significantly faster:

```wolfram
result = Lenia[grid, 1000, Method -> "Rust"];
```

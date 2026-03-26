# Lenia for Wolfram Language

A professional-grade implementation of the **Lenia** continuous cellular automaton in the Wolfram Language, based on [Bert Chan's original work](https://chakazul.github.io/lenia.html).

## Features

- **Continuous State and Space**: Full support for continuous grid values on [0, 1].
- **FFT-based Convolution**: Fast, correct circular convolution using Fourier transforms with periodic boundary conditions.
- **Gaussian Growth Function**: Standard bell-shaped growth function `G(n) = 2 * Exp[-(n-μ)²/(2σ²)] - 1`.
- **Bump Kernel**: Exponential bump kernel `K(r) = Exp[4 - 1/(r(1-r))]` for `0 < r < 1`.

## Installation

Clone the repository and load the paclet:

```wolfram
PacletDirectoryLoad["path/to/lenia/Lenia"]
Needs["Lenia`"]
```

## Quick Start

### Creating a Viable Seed

Lenia creatures have specific spatial structures. A ring-shaped seed works well:

```wolfram
gridSize = {128, 128};
R = 13;

(* Ring-shaped seed: density peaks at r ~ 0.5R from center *)
grid = N @ Table[
  Module[{r = Sqrt[(i - 64)^2 + (j - 64)^2] / R},
    If[r < 1.0, Exp[-((r - 0.5)^2) / (2 * 0.15^2)], 0.0]
  ],
  {i, 1, 128}, {j, 1, 128}
];
```

### Running the Simulation

```wolfram
result = Lenia[grid, 200, "Mu" -> 0.15, "Sigma" -> 0.015, "Radius" -> R];
ArrayPlot[result, ColorFunction -> "SolarColors"]
```

### Returning History (for Animation)

```wolfram
history = Lenia[grid, 200, "ReturnHistory" -> True];
ListAnimate[ArrayPlot[#, ColorFunction -> "TemperatureMap"] & /@ history]
```

### Using LeniaBlob

For simpler experiments, `LeniaBlob` creates a Gaussian density:

```wolfram
grid = LeniaBlob[{128, 128}, {64, 64}, 10, 1.0];
```

> **Note:** Gaussian blobs may need different `Mu`/`Sigma` values than ring-shaped seeds.

## API Reference

### `Lenia[grid, steps, options]`

Runs the Lenia simulation using FFT-based circular convolution.

| Option | Default | Description |
|--------|---------|-------------|
| `"Mu"` | 0.15 | Growth function center |
| `"Sigma"` | 0.015 | Growth function width |
| `"DT"` | 0.1 | Time step (= 1/T where T=10) |
| `"Radius"` | 13 | Kernel radius R |
| `"ReturnHistory"` | False | Return list of all states |

### `LeniaKernel[gridSize, radius]`

Generates the normalized bump kernel at full grid size (FFT-shifted).

### `LeniaKernelFFT[gridSize, radius]`

Returns `Fourier[LeniaKernel[...], FourierParameters -> {1, -1}]` for use in `LeniaStep`.

### `LeniaStep[grid, kernelFFT, mu, sigma, dt]`

Performs a single Lenia iteration using FFT convolution.

### `LeniaGrowth[n, mu, sigma]`

The bell-shaped growth function: `2 * Exp[-(n-μ)²/(2σ²)] - 1`.

### `LeniaBlob[dims, center, radius, amplitude]`

Creates a Gaussian density blob. `amplitude` defaults to 1.0.

## Troubleshooting

### Why does my simulation decay to zero?

Lenia is sensitive to the match between initial conditions and parameters:

- **Use ring-shaped seeds**: The bump kernel weights a ring at `r ≈ 0.5`. Seeds need density in this ring to trigger growth.
- **Gaussian blobs alone often fail**: Their density is concentrated at the center (`r ≈ 0`), missing the kernel's ring peak.
- **Parameter sensitivity**: Growth only occurs when the convolved potential `n` is within `μ ± ~3σ` of the growth peak.

## License

MIT License.

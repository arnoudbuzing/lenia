# Lenia for Wolfram Language

A high-performance, professional-grade implementation of the **Lenia** continuous cellular automaton for the Wolfram Language, inspired by [Bert Chan's original work](https://chakazul.github.io/lenia.html).

## What is Lenia?

**Lenia** (from the Latin *lenia*, meaning "smooth") is a family of cellular automata that exists in continuous space, time, and states. Unlike traditional discrete CA like Conway's *Game of Life*, Lenia uses:
- **Continuous Grids**: Cell values are real numbers in the range $[0.0, 1.0]$.
- **Circular Convolution**: Neighborhood influences are computed via spatial kernels.
- **Growth Functions**: A bell-shaped function determines how a cell's state changes based on its neighborhood density.

Lenia is capable of producing remarkably lifelike "organisms" with complex behaviors, including stable gliders, rotating colonies, and self-replicating patterns.

## What this Paclet Provides

This paclet brings Lenia to the Wolfram Language with a focus on **performance**, **usability**, and **extensibility**:

- 🚀 **Dual Backends**: Run simulations using the optimized pure Wolfram Language implementation (via `Fourier`) or the ultra-fast **Rust-based backend** for large-scale experiments.
- 🧬 **Organism Seeds**: Built-in library of initial conditions (`LeniaSeed`) to quickly generate "Rings", "Constellations", and other complex lifeforms.
- 💾 **Memory Efficiency**: Transparent support for `NumericArray` objects to handle large grids and long histories without exhausting system RAM.
- 🛠️ **Configurable Kernels**: Support for multiple kernel profiles beyond the standard "Bump", including `"GaussianRing"`, `"SmoothLife"`, and `"StepRing"`.
- 📊 **Rich Visualization**: Seamless integration with Wolfram's visualization tools for creating animations and scientific plots.

---

## Installation

Clone this repository and load the paclet directory:

```wolfram
PacletDirectoryLoad["/path/to/lenia/Lenia"]
Needs["ArnoudBuzing`Lenia`"]
```

## Quick Start

### 1. Generate an Organism
Generate a 128x128 grid with a random complex organism:
```wolfram
grid = LeniaSeed[128];
```

### 2. Run the Simulation
Perform 100 steps of evolution using the default parameters:
```wolfram
result = Lenia[grid, 100];
ArrayPlot[result, ColorFunction -> "SolarColors"]
```

### 3. Animate the History
Capture every step and animate the result:
```wolfram
history = Lenia[grid, 100, "ReturnHistory" -> True];
ListAnimate[ArrayPlot[#, ColorFunction -> "SolarColors", Frame -> False] & /@ history]
```

---

## Advanced Features

### Rust Acceleration
For high-performance simulations, switch to the Rust backend:
```wolfram
result = Lenia[grid, 1000, Method -> "Rust"];
```

### Alternative Kernels
Explore different physics by changing the kernel profile:
```wolfram
result = Lenia[grid, 50, "LeniaKernel" -> "SmoothLife"];
```

### Full Documentation
For a detailed breakdown of all arguments, options, and helper functions, see the [Lenia Documentation](docs/Lenia.md).

---

## API Summary: `Lenia` Options

| Option | Default | Description |
| :--- | :--- | :--- |
| `"Mu"` | `0.15` | Center of the growth function peak. |
| `"Sigma"` | `0.015` | Width of the growth function peak. |
| `"DT"` | `0.1` | Integration time step. |
| `"Radius"` | `13` | Spatial radius of the kernel influence. |
| `"LeniaKernel"` | `"Bump"` | Kernel profile type (e.g., `"GaussianRing"`, `"StepRing"`). |
| `"ReturnHistory"`| `False` | Return all intermediate states as a `NumericArray`. |
| `Method` | `"Wolfram"` | Backend selector: `"Wolfram"` or `"Rust"`. |

---

## License
MIT License.

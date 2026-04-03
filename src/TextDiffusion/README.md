# TextDiffusion - Text-to-CFD Visualization

A C++ application that renders text using TrueType fonts and uses the resulting pattern as initial conditions for a temperature diffusion simulation in OpenFOAM.

## Overview

**TextDiffusion** is an innovative tool that visualizes text through computational fluid dynamics. It:

1. **Renders text** from TrueType font files into a 2D bitmap matrix
2. **Maps the bitmap** to a computational mesh grid
3. **Initializes a temperature field** based on the text pattern
4. **Simulates diffusion** to create an animated text evolution over time

This creates a visually interesting effect where text gradually diffuses across a domain, perfect for visualization, artistic projects, or educational purposes.

## Building

### Prerequisites

- OpenFOAM (v2006 or later recommended)
- C++17 compiler (g++, clang++)
- CMake or OpenFOAM Make system
- TrueType font file (system default: `/usr/share/fonts/truetype/ubuntu/Ubuntu-M.ttf`)

### Build Instructions

Using the OpenFOAM Make system:

```bash
cd src/TextDiffusion/app
wmake
```

## Usage

### Basic Command

```bash
cd src/TextDiffusion/case
textDiffusion
```

### With Custom Parameters

```bash
textDiffusion \
  -case /path/to/case \
  -text "Hello World" \
  -font /path/to/font.ttf \
  -fontsize 48 \
  -scale 1.5
```

### Parameters

| Parameter  | Type   | Default                                         | Description                      |
| ---------- | ------ | ----------------------------------------------- | -------------------------------- |
| `text`     | string | "Hello"                                         | Text to render                   |
| `font`     | path   | `/usr/share/fonts/truetype/ubuntu/Ubuntu-M.ttf` | TrueType font file path          |
| `fontsize` | float  | 30.0                                            | Font size in pixels              |
| `scale`    | scalar | 1.0                                             | Temperature field scaling factor |

### OpenFOAM Configuration

The case directory (`case/`) contains:

**`system/blockMeshDict`**: Defines mesh block structure
- Example: 40×20 cells with horizontal stretch
- Adjustable via `simpleGrading` parameter

**`system/controlDict`**: Simulation time control
- `startTime`: 0
- `endTime`: 0.0004
- `deltaT`: 0.0001

**`system/fvSchemes`**: Discretization schemes
- Convection: Upwind
- Laplacian: Gaussian with 90° correction

**`constant/transportProperties`**: Physical properties
- `DT`: Diffusion coefficient (thermal diffusivity)
- `nu`: Kinematic viscosity

## Physical Model

The application solves the unsteady convection-diffusion equation:

$$\frac{\partial T}{\partial t} + \nabla \cdot (\mathbf{U} T) = \nabla \cdot (D \nabla T)$$

Where:
- **T**: Scalar field (temperature, initialized from text)
- **U**: Velocity field (from PISO solver)
- **D**: Diffusion coefficient (thermal diffusivity)

The PISO algorithm provides temporal accuracy and pressure-velocity coupling.

## License

MIT License - See header in source files

## Author

JackLee (c) 2026

## References

- [OpenFOAM Documentation](https://www.openfoam.com/documentation/)
- [stb_truetype](https://github.com/nothings/stb/blob/master/stb_truetype.h)
- Issa, R.I. (1986). "Solution of the implicitly discretised fluid flow equations by operator-splitting"

# Bouncing Ball Application

## Overview

The Bouncing Ball application is an OpenFOAM-based C++ simulation that models a solid particle (ball) bouncing and interacting with patch boundaries in a fluid domain. The solver tracks particle motion under gravity, handles collisions with walls elastically, and simulates the particle's interaction with the surrounding fluid.

## Configuration

### Particle Properties (`constant/particleProperties`)

Key parameters controlling the particle behavior:

```
rhop          # Particle density (kg/m³)
e             # Elasticity coefficient (0-1, where 1 is perfectly elastic)
mu            # Particle dynamic viscosity (Pa·s)
diameter      # Particle diameter (m)
velocity      # Initial particle velocity (m/s) - vector per particle
position      # Initial particle position (m) - vector per particle
```

## Building

The application uses OpenFOAM's build system (`wmake`):

```bash
cd src/BouncingBall/app
wmake
```

For testing:

```bash
cd src/BouncingBall/app/test
wmake
```

## Running the Simulation

### Individual Particle Simulation

1. Navigate to the case directory:
   ```bash
   cd src/BouncingBall/case
   ```

2. Clean previous run results (optional):
   ```bash
   ./Allclean
   ```

3. Generate the mesh:
   ```bash
   blockMesh
   ```

4. Run the solver:
   ```bash
   patchInteractSolidParticleCloud  # or the executable name from your build
   ```


## Physics

### Motion Equation

The particle motion is governed by:

$$\frac{d\vec{U}_p}{dt} = \vec{g} + \frac{1}{\rho_p V_p}\vec{F}_{drag}$$

where:
- $\vec{U}_p$ = particle velocity
- $\vec{g}$ = gravitational acceleration
- $\rho_p$ = particle density
- $V_p$ = particle volume
- $\vec{F}_{drag}$ = drag force from surrounding fluid

### Collision Model

- **Elasticity Coefficient** ($e$): Controls energy loss during collisions
  - $e = 1$: perfectly elastic (no energy loss)
  - $e = 0$: perfectly inelastic (sticks to surface)
  - $0 < e < 1$: typical bouncing behavior

- Normal component of velocity after collision: $U_n' = -e \cdot U_n$


## License

Copyright (c) 2026 - Released under the MIT License
See LICENSE file for details

## References

- OpenFOAM Documentation: https://www.openfoam.com/documentation

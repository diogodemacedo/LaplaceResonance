# Laplace Resonance and Orbital Dynamics of Jupiter's Moons

Numerical study of the **Io–Europa–Ganymede 1:2:4 Laplace resonance** using GNU Octave. The project models the Jupiter–Io–Europa–Ganymede system and investigates how changing Europa's mass perturbs the orbital eccentricities of the resonant moons.

This repository archives and reconstructs a two-week undergraduate computational-physics mini-project completed in 2019.

## Project overview

The original project had two main goals:

1. Develop a numerical model of the Galilean system containing Jupiter, Io, Europa, and Ganymede.
2. Study how the Laplace-resonant system responds when Europa's mass is artificially increased.

The model uses Newtonian gravity in a simplified coplanar four-body system. Orbital trajectories and osculating eccentricities are propagated numerically and compared between the nominal system and modified-Europa-mass cases.

## Repository structure

```text
.
├── README.md
├── CV_ENTRY.tex
├── GITHUB_SETUP.md
├── report/
│   └── Three-Body_Laplace_Resonance_in_Jupiters_Moons.pdf
├── src/
│   ├── run_laplace_resonance.m
│   ├── simulate_laplace.m
│   ├── gravity_acceleration.m
│   └── orbital_elements_relative.m
└── results/
    └── .gitkeep
```

## Historical provenance and reconstruction note

The original 2019 report was co-authored by **Diogo de Macedo and Stephanie Berntsen Croft**. The report contains the main Octave simulation script in its appendix, but the original source files are no longer available and the appendix calls helper routines (`gravity` and `orbital_params`) that are not included in the PDF.

For that reason, the code in `src/` is a **transparent reconstruction** of the model described in the report rather than a claim to be a byte-for-byte copy of the original 2019 source code. The reconstructed implementation preserves the documented physical setup and methodology while supplying the missing helper functionality needed to make the repository runnable.

## Numerical model

The reconstruction follows the approach documented in the report:

- Newtonian N-body gravitational dynamics
- Jupiter, Io, Europa, and Ganymede
- 2D/coplanar initial conditions
- normalized mass, length, and time scales
- approximate circular initial orbital velocities
- center-of-mass velocity correction
- leapfrog-equivalent / velocity-Verlet time integration
- tracking of osculating semi-major axis and eccentricity relative to Jupiter
- comparison of the nominal Europa mass with an artificially increased Europa mass

The default driver runs two cases:

- `Europa mass factor = 1`
- `Europa mass factor = 100`

for approximately **40 Ganymede orbital periods**.

## Running the simulation

### Requirements

- GNU Octave

No external Octave packages are required.

### Run

From the repository root:

```bash
cd src
octave --quiet run_laplace_resonance.m
```

Generated CSV files and figures are written to `results/`.

## Outputs

For each Europa mass factor, the script generates:

- orbital-trajectory plot
- eccentricity-vs-time plot
- CSV file containing sampled positions, semi-major axes, and eccentricities

Example output names:

```text
results/orbits_mass_factor_1.png
results/eccentricity_mass_factor_1.png
results/laplace_mass_factor_1.csv
results/orbits_mass_factor_100.png
results/eccentricity_mass_factor_100.png
results/laplace_mass_factor_100.csv
```

## Results from the original report

The original study found qualitatively that increasing Europa's mass strongly perturbed the orbital behavior of Io and Ganymede. In the two-orders-of-magnitude mass case, the eccentricity evolution of the system changed substantially, while Europa's own eccentricity decreased by roughly 0.01 in the reported simulation. The report also noted that Ganymede did not settle into the same apparent eccentricity equilibrium within the simulated interval after Europa's mass was increased.

The original analysis was intentionally qualitative because of the short two-week project duration and the computational cost of the simulations.

## Limitations

This is an educational orbital-dynamics model, **not a flight-dynamics or mission-design tool**. Important simplifications include:

- planar dynamics
- simplified initial conditions
- no Jupiter oblateness / J2 terms
- no tidal dissipation
- no relativistic corrections
- no detailed ephemeris initialization
- no uncertainty propagation or orbit determination

The project is useful as a compact demonstration of orbital-mechanics reasoning, numerical integration, scientific programming, and sensitivity analysis.

## Report

The original project report is included in `report/`:

**Three-Body Laplace Resonance in Jupiter's Moons**  
Diogo de Macedo and Stephanie Berntsen Croft, 2019.

Because the report is co-authored, confirm that public redistribution is acceptable to all authors and consistent with any university/course requirements before making the repository public.

## License

No open-source license is attached to this repository by default. The report is co-authored, and the reconstructed code is based on the methodology and partial code preserved in that report. Add a license only after confirming the desired licensing terms with the relevant contributors.

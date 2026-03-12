# URI BEM Analysis Pipeline

![MATLAB](https://img.shields.io/badge/MATLAB-R2021b+-orange?logo=mathworks)
![Python](https://img.shields.io/badge/Python-3.6+-blue?logo=python&logoColor=white)
![Nemoh](https://img.shields.io/badge/BEM-Nemoh-green)
![License](https://img.shields.io/badge/License-MIT-lightgrey)

**Preliminary hydrodynamic analysis of floating structures using boundary element methods.**

An integrated, open-source toolchain that couples the [Nemoh](https://gitlab.com/lheea/Nemoh) BEM solver with custom MATLAB and Python scripts to deliver a complete workflow — from panel-mesh generation through frequency-domain BEM solutions to time-domain Cummins-equation simulations and parametric design sweeps.

Developed at the **University of Rhode Island (URI)**.

---

## Pipeline Overview

| Stage | Tool | Description |
|-------|------|-------------|
| **1. Mesh Generation** | `generate_case.py` | Parametric panel mesh for twin-pontoon catamaran geometry |
| **2. BEM Solve** | Nemoh (pre/Solver/post) | Added mass, radiation damping, and excitation transfer functions |
| **3. Post-Processing** | `remove_irreg_freq.py` | Savitzky–Golay irregular-frequency removal + WAMIT conversion |
| **4. Time-Domain Sim** | `run_nemoh_timeseries.m` | Cummins equation with Prony state-space radiation, irregular waves |
| **5. Parametric Sweep** | `sweep_sea_states.m` | 1 % exceedance force contour maps over *(H_s, T_p)* grid |
| **6. Mooring Analysis** | `compute_mooring_stiffness.m` | Catenary & taut-wire solvers, linearised 6×6 stiffness matrix |

---

## Key Features

- **Full 6-DOF Cummins equation** with state-space radiation convolution via Prony series
- **Wavenumber-space (k-space) wave discretisation** — incommensurate frequencies avoid signal repetition; 128 components suffice for converged statistics
- **JONSWAP / Pierson–Moskowitz** irregular wave synthesis with configurable spectral parameters
- **Catenary and taut-wire mooring** solvers with Newton–Raphson root finding and numerical linearisation
- **Irregular-frequency removal** using Savitzky–Golay filtering with MAD-based outlier detection
- **Exceedance probability analysis** — 1 % and 5 % exceedance values with transient clipping
- **Parametric H_s–T_p sweeps** — batch contour maps for preliminary structural design screening
- **Publication-quality PDF figures** — all plots exported automatically
- **Modular architecture** — easy to swap geometries, add external forces, or couple with other solvers
- **Fully open-source** — Python (stdlib + NumPy/SciPy) and MATLAB with no commercial toolbox dependencies

---

## Example Results

### Frequency-Domain Coefficients

<p align="center">
  <img src="docs/source/figures/irreg_freq_heave_detail.png" alt="Heave added mass — irregular frequency removal" width="70%">
</p>

*Heave added mass before and after irregular-frequency removal, showing clean coefficient curves ready for time-domain use.*

### Added Mass & Damping (All DOFs)

<p align="center">
  <img src="docs/source/figures/irreg_freq_addedmass.png" alt="Added mass matrix" width="70%">
</p>

---

## Quick Start

### Prerequisites

- **Python ≥ 3.6** with NumPy, SciPy, Matplotlib
- **MATLAB R2021b+** (no additional toolboxes required)
- **Nemoh** — pre-compiled Windows binaries included in `bin_windows/`; Linux/macOS users should compile from [source](https://gitlab.com/lheea/Nemoh)

### Steps

1. **Clone the repository**
   ```bash
   git clone https://github.com/<your-org>/uri-bem-pipeline.git
   cd uri-bem-pipeline
   ```

2. **Generate the mesh and Nemoh inputs**
   ```bash
   cd SimplePontoon
   python generate_case.py
   ```

3. **Run the Nemoh BEM solver**
   ```bash
   ../bin_windows/preProcessor
   ../bin_windows/Solver
   ../bin_windows/postProcessor
   ```

4. **Remove irregular frequencies**
   ```bash
   python remove_irreg_freq.py
   ```

5. **Run the time-domain simulation** (MATLAB)
   ```matlab
   run_nemoh_timeseries
   ```

6. **Run the parametric sweep** (MATLAB)
   ```matlab
   sweep_sea_states
   ```

---

## Directory Structure

```
uri-bem-pipeline/
├── README.md
├── LICENSE
├── .gitignore
├── .readthedocs.yaml
│
├── bin_windows/              # Pre-compiled Nemoh executables (Windows)
│   ├── preProcessor.exe
│   ├── Solver.exe
│   └── postProcessor.exe
│
├── SimplePontoon/            # Complete example case
│   ├── generate_case.py          # Mesh generation (Python)
│   ├── remove_irreg_freq.py     # Irregular-frequency removal
│   ├── nemoh_to_wamit.py        # WAMIT format conversion
│   ├── plot_results.m            # Frequency-domain plotting
│   ├── run_nemoh_timeseries.m   # Time-domain solver (MATLAB)
│   ├── sweep_sea_states.m       # Parametric Hs-Tp sweep
│   ├── compute_mooring_stiffness.m  # Mooring stiffness calculator
│   ├── Nemoh.cal                 # Nemoh master input
│   ├── mesh/                     # Panel mesh files
│   ├── Mechanics/                # Mass & stiffness matrices
│   ├── results/                  # Nemoh output (.tec files)
│   ├── Wamit_format/             # Converted WAMIT files
│   └── pdf_outputs/              # Exported figures
│
├── docs/                     # Sphinx documentation (Read the Docs)
│   ├── Makefile
│   ├── make.bat
│   ├── requirements.txt
│   └── source/
│       ├── conf.py
│       ├── index.rst
│       ├── introduction.rst
│       ├── examples.rst
│       ├── references.rst
│       ├── bibliography.bib
│       ├── theory/               # Mathematical foundations
│       ├── pipeline/             # Usage guides
│       ├── api/                  # API reference
│       └── figures/              # Documentation figures
│
└── Nemoh-master/             # Nemoh source (user-provided, gitignored)
```

---

## Documentation

Full documentation is available on **Read the Docs**:

📖 **[URI BEM Analysis Pipeline Documentation](https://uri-bem-pipeline.readthedocs.io/)**

The documentation includes:
- **Theory** — BEM formulation, Cummins equation, Prony method, wave spectra, hydrostatics, mooring
- **Pipeline Guide** — step-by-step usage for each script
- **API Reference** — function signatures and parameter tables
- **Examples** — complete worked example from mesh to sweep

To build locally:
```bash
cd docs
pip install -r requirements.txt
make html
# Open build/html/index.html
```

---

## Citation

If you use this pipeline in academic work, please cite the underlying tools and theory:

```bibtex
@article{Babarit2015,
  author  = {Babarit, A. and Delhommeau, G.},
  title   = {Theoretical and numerical aspects of the open source {BEM}
             solver {NEMOH}},
  journal = {Proc. 11th European Wave and Tidal Energy Conf. (EWTEC)},
  year    = {2015}
}

@article{Cummins1962,
  author  = {Cummins, W. E.},
  title   = {The impulse response function and ship motions},
  journal = {Schiffstechnik},
  volume  = {9},
  pages   = {101--109},
  year    = {1962}
}

@article{PerezFossen2009,
  author  = {Perez, T. and Fossen, T. I.},
  title   = {A {Matlab} toolbox for parametric identification of
             radiation-force models of ships and offshore structures},
  journal = {Modeling, Identification and Control},
  volume  = {30},
  number  = {1},
  pages   = {1--15},
  year    = {2009}
}
```

---

## License

This project is released under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## Acknowledgments

- **[Nemoh](https://gitlab.com/lheea/Nemoh)** — open-source BEM solver developed at École Centrale de Nantes / LHEEA
- **University of Rhode Island** — Ocean Engineering department
- The JONSWAP spectral model is based on the seminal work of Hasselmann et al. (1973)

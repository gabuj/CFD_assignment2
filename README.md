# CFD Assignment 2

> Computational Fluid Dynamics coursework repository focusing on the numerical simulation, discretization, and analysis of fluid flow and transport phenomena.

## 📋 Table of Contents
- [Overview](#overview)
- [Governing Equations & Mathematical Model](#governing-equations--mathematical-model)
- [Numerical Methods & Discretization](#numerical-methods--discretization)
- [Repository Structure](#repository-structure)
- [Prerequisites & Dependencies](#prerequisites--dependencies)
- [Installation & Setup](#installation--setup)
- [Usage](#usage)
- [Results & Discussion](#results--discussion)
- [Author](#author)

---

## 🔍 Overview
This repository contains the source code, numerical solvers, and post-processing scripts for **Assignment 2** of the Computational Fluid Dynamics course. 

The primary objective of this assignment is to simulate [insert specific problem description, e.g., 2D steady-state convection-diffusion / lid-driven cavity flow / compressible nozzle flow] using custom numerical algorithms implemented in [Python / MATLAB / C++].

---

## 📐 Governing Equations & Mathematical Model
The physical system is modeled using the appropriate conservation laws:

- **Continuity Equation (Conservation of Mass):**
  $$\nabla \cdot \mathbf{u} = 0$$

- **Momentum Equations (Navier-Stokes):**
  $$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \mu \nabla^2 \mathbf{u} + \mathbf{f}$$

### Boundary Conditions
- **Inlet:** [e.g., Uniform velocity profile $U_0 = 1.0 \text{ m/s}$]
- **Outlet:** [e.g., Zero-gradient / Convective outlet condition]
- **Walls:** [e.g., No-slip condition ($\mathbf{u} = 0$) on stationary walls]

---

## ⚙️ Numerical Methods & Discretization
- **Spatial Discretization:** [e.g., Second-order central differencing for diffusion terms, Upwind / QUICK scheme for convection terms]
- **Temporal Discretization:** [e.g., First-order Euler explicit / Crank-Nicolson implicit time integration]
- **Pressure-Velocity Coupling:** [e.g., SIMPLE algorithm / Fractional step method]
- **Solver & Convergence:** Iterative linear solvers (e.g., Gauss-Seidel with relaxation) with residual convergence criteria set to $\le 10^{-6}$.

---

## 📂 Repository Structure
```text
CFD_assignment2/
│
├── data/               # Input grid files, meshes, and boundary data
├── src/                # Core source code
│   ├── solver.py       # Main CFD solver implementation
│   ├── mesh.py         # Grid generation and geometry setup
│   └── boundary.py     # Boundary condition handlers
├── scripts/            # Execution and batch runner scripts
├── results/            # Output data files, CSV logs, and flow fields
├── figures/            # Generated plots, residuals, and contour maps
├── report/             # Assignment report (PDF / LaTeX source)
└── README.md           # Project documentation

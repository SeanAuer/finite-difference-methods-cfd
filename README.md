## Finite Difference Methods for Computational Fluid Dynamics
---

### Overview

This set of Jupyter Notebooks is designed to build a foundational understanding of how **Finite Difference Methods (FDM)** are used to computationally solve problems in fluid dynamics.

We focus on a set of commonly encountered governing equations that describe the physics of fluid flows:

- **Potential Flow**  
- **Transport Equations**  
- **The Wave Equation**  
- **The Euler Equations**  
- **The Navier–Stokes Equations**

---

### Objectives

Throughout this course, the goal is to help the reader develop insight into:

- The structure and components of each governing equation  
- Discretization techniques and their application  
- The role of problem dimensionality (**1-D**, **2-D**, **3-D**)  
- Key numerical considerations, including:
  - Truncation error  
  - Stability analysis  
  - Convergence behavior  

Whether you're a student, researcher, or engineer, this series aims to bridge theoretical fluid dynamics with practical numerical implementation.

---

### What Are Finite Difference Methods?

Finite Difference Methods (FDM) are a class of numerical techniques used to approximate the solutions of differential equations by replacing derivatives with difference equations.  

$\frac{\partial \phi}{\partial x} \rightarrow \frac{\Delta \phi}{\Delta x}$  

These methods are especially effective for solving partial differential equations (PDEs) that govern fluid flow, heat transfer, wave propagation, and other physical processes.  

---

### Real-World Applications of FDM

While more advanced solvers today use Finite Volume or Finite Element Methods, FDM remains foundational in many fields, especially for prototyping, reduced-order modeling, and education.

Notable uses of FDM:
- **Early CFD solvers** (e.g., early NASA tools like NASTRAN prototypes)
- **Weather models** from institutions like **NOAA**
- **Educational tools** for aerospace and mechanical engineers
- **Embedded systems** or fast solvers where memory constraints require lightweight numerical schemes

FDM also plays a critical role in validating results from more complex solvers by offering analytically traceable approximations.

---

### How to Use These Notebooks

1. Clone the repository or download the ZIP  
2. Open in JupyterLab or Jupyter Notebook  
3. Follow the recommended notebook order  
4. Ensure you have the required Python libraries (`NumPy`, `Matplotlib`, etc.)

---

### Citation

If you use or build upon this work, please cite:

**Sean Auer**, *Finite Difference Methods for Computational Fluid Dynamics*, 2025.  
GitHub: [https://github.com/Seanauer](https://github.com/Seanauer)  
LinkedIn: *[Insert LinkedIn URL]*

> “An accessible and rigorous breakdown of FDM-based solvers for classic fluid dynamics problems.”

BibTeX:
```bibtex
@misc{auer2025fdm,
  author       = {Sean Auer},
  title        = {Finite Difference Methods for Computational Fluid Dynamics},
  year         = {2025},
  howpublished = {\url{https://github.com/Seanauer/fdm-cfd}},
  note         = {Accessed: YYYY-MM-DD}
}
```
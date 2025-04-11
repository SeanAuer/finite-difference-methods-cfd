# Finite Difference Methods for Computational Fluid Dynamics  

## Table of Contents
- [Overview](#overview)
- [Objectives](#objectives)
- [What Are Finite Difference Methods?](#what-are-finite-difference-methods)
- [Real-World Applications of FDM](#real-world-applications-of-fdm)
- [How to Use These Notebooks](#how-to-use-these-notebooks)
- [Citation](#citation)
---  

## Overview

This set of Jupyter Notebooks is designed to build a foundational understanding of how **Finite Difference Methods (FDM)** are used to computationally solve problems in fluid dynamics.

We focus on a set of commonly encountered governing equations that describe the physics of fluid flows:

- **Potential Flow**  
- **Transport Equations**  
- **The Wave Equation**  
- **The Euler Equations**  
- **The Navier–Stokes Equations**

---

## Objectives

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

## What Are Finite Difference Methods?

Finite Difference Methods (FDM) are a class of numerical techniques used to approximate the solutions of differential equations by replacing derivatives with algebraic difference equations.  

$$
\frac{\partial \phi}{\partial x} \quad \longrightarrow \quad \frac{\Delta \phi}{\Delta x}
$$  

These approximations are grounded in the **physical interpretation of a derivative**—the rate of change of one quantity with respect to another. A classic example is the **forward difference**, derived from the definition of the derivative:

$$
\frac{d f(x)}{d x} = \lim_{\Delta x \rightarrow 0} \frac{f(x + \Delta x) - f(x)}{\Delta x}
$$

To construct a usable numerical method, we **drop the limit** and treat $\Delta x$ as a small but finite step size. This yields an approximation:

$$
\frac{d f(x)}{d x} \approx \frac{f(x + \Delta x) - f(x)}{\Delta x}
$$

This is known as the **forward difference formula**, so called because for $\Delta x$ > 0, each 'step' in x is seen as moving 'forward,' or to the right on an x-axis, and it becomes more accurate as $\Delta x$ becomes smaller.

In numerical applications, the function is assumed to be defined only at a set of **discrete grid points**, such as:

$$
x = [x_0, x_1, x_2, x_3, \dots] = [0, \Delta x, 2\Delta x, 3\Delta x, \dots]
$$

Derivatives are then approximated using values at neighboring grid points. For example, the derivative at a point $x_i$ might be approximated using $f(x_{i+1})$ and $f(x_i)$ (forward difference), or $f(x_{i-1})$ and $f(x_i)$ (backward difference), or a symmetric combination (central difference).

This shift from **continuous calculus** to **discrete algebra** allows us to transform differential equations into systems of equations that can be solved using a computer.


### Finite Difference Approximations

To solve differential equations numerically, we replace calculus-based derivatives with **algebraic approximations** on a discrete grid — this is the foundation of **Finite Difference Methods (FDM)**.

Let $\phi_i = \phi(x_i)$ be the value of a function at point $x_i$, with uniform grid spacing $\Delta x$.

###  **First Derivative: $\frac{\partial \phi}{\partial x}$**

| Method | Approximation | Order |
|--------|---------------|-------|
| **Forward Difference** | $\displaystyle \frac{\phi_{i+1} - \phi_i}{\Delta x}$ | 𝒪($\Delta x$) |
| **Backward Difference** | $\displaystyle \frac{\phi_i - \phi_{i-1}}{\Delta x}$ | 𝒪($\Delta x$) |
| **Central Difference** | $\displaystyle \frac{\phi_{i+1} - \phi_{i-1}}{2\Delta x}$ | 𝒪($\Delta x^2$) |

### **Second Derivative: $\frac{\partial^2 \phi}{\partial x^2}$**

| Method               | Approximation                                                             | Order         |
|----------------------|---------------------------------------------------------------------------|---------------|
| **Forward-biased**   | $\displaystyle \frac{-\phi_{i+2} + 4\phi_{i+1} - 5\phi_i + 2\phi_{i-1}}{\Delta x^2}$ | 𝒪($\Delta x^2$) |
| **Backward-biased**  | $\displaystyle \frac{-\phi_{i-2} + 4\phi_{i-1} - 5\phi_i + 2\phi_{i+1}}{\Delta x^2}$ | 𝒪($\Delta x^2$) |
| **Central Difference** | $\displaystyle \frac{\phi_{i+1} - 2\phi_i + \phi_{i-1}}{\Delta x^2}$              | 𝒪($\Delta x^2$) |
| **Five-Point Stencil** *(optional)* | $\displaystyle \frac{-\phi_{i+2} + 16\phi_{i+1} - 30\phi_i + 16\phi_{i-1} - \phi_{i-2}}{12\Delta x^2}$ | 𝒪($\Delta x^4$) |



**🔍 What does "Order" mean?**  
The *order of a finite difference method* tells us how the error behaves as we shrink the grid spacing $\Delta x$.  
- A method with **𝒪($\Delta x$)** error (first-order) means halving the grid spacing cuts the error in half.
- A method with **𝒪($\Delta x^2$)** error (second-order) means halving the grid spacing cuts the error by a factor of **four**.
- Higher order = faster convergence as the grid gets finer.

> We are not doing calculus on functions — we are using **algebraic approximations** to estimate derivatives using values at nearby grid points.

These finite difference formulas form the building blocks of all the numerical solvers used in this project.

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

## How to Use These Notebooks

1. Clone the repository or download the ZIP  
2. Open in JupyterLab or Jupyter Notebook  
3. Follow the recommended notebook order  
4. Ensure you have the required Python libraries (`NumPy`, `Matplotlib`, etc.)

---

## Citation

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
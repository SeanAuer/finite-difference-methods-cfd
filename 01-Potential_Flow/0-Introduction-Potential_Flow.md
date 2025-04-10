# 1. Introduction to Potential Flow

---

## Table of Contents
1. [Derivation of the Potential Flow Equation](#11-derivation-of-the-potential-flow-equation)  
2. [Assumptions](#12-assumptions)  
3. [Mathematical Characterization](#13-mathematical-characterization)  
4. [Numerical Solution Methods (Finite Difference)](#14-numerical-solution-methods-finite-difference)  
   - [1D Laplace Equation](#141-one-dimensional-laplace-equation)  
   - [2D Laplace Equation](#142-two-dimensional-laplace-equation)  
5. [Solution Methods for Laplace Equation](#15-solution-methods-for-laplace-equation)  
6. [Python Implementation Note](#16-python-implementation-note)

---

## 1.1. Derivation of the Potential Flow Equation

Potential flow is based on the idea that the velocity field $\vec{u}$ can be represented as the gradient of a scalar potential function:

$$
\vec{u} = \nabla \phi
$$

If the flow is **incompressible**, we require:

$$
\nabla \cdot \vec{u} = 0
$$

Substituting $\vec{u} = \nabla \phi$ gives:

$$
\nabla \cdot (\nabla \phi) = \nabla^2 \phi = 0
$$

This is **Laplace’s equation**, the governing equation of potential flow:

$$
\nabla^2 \phi = 0
$$

---

## 1.2. Assumptions

Potential flow theory rests on several simplifying assumptions:

- The fluid is **inviscid** (no viscosity terms)  
- The flow is **incompressible** ($\nabla \cdot \vec{u} = 0$)  
- The flow is **irrotational** ($\nabla \times \vec{u} = 0$)  
- The flow is **steady** ($\frac{\partial \vec{u}}{\partial t} = 0$)  
- The fluid is treated as a **continuous** medium and is **single-phase**

These assumptions reduce the full Navier–Stokes equations to Laplace’s equation.

---

## 1.3. Mathematical Characterization

Let’s characterize Laplace’s equation mathematically:

- **Classification:** Partial Differential Equation (PDE)  
- **Form:** Second-order, linear  
- **Category:** Elliptic  
- **Time-dependence:** Steady-state (no explicit time term)  
- **Independent variables:** Spatial coordinates (e.g., $x$, $y$, $z$)  
- **Dependent variable:** Scalar potential field $\phi(x, y)$

---

## 1.4. Numerical Solution Methods (Finite Difference)

To solve $\nabla^2 \phi = 0$ numerically, we discretize the Laplacian using **finite difference approximations**. We begin with the one-dimensional form and then extend to two dimensions.

---

### 1.4.1. One-Dimensional Laplace Equation

In 1D, Laplace’s equation reduces to:

$$
\frac{d^2 \phi}{dx^2} = 0
$$

#### Central Difference Approximation:

Using the second-order central difference:

$$
\frac{d^2 \phi}{dx^2} \approx \frac{\phi_{i+1} - 2\phi_i + \phi_{i-1}}{\Delta x^2}
$$

Substituting into the PDE:

$$
\frac{\phi_{i+1} - 2\phi_i + \phi_{i-1}}{\Delta x^2} = 0
\quad \Rightarrow \quad
\phi_i = \frac{(\phi_{i+1} + \phi_{i-1})}{2}
$$

This forms the update equation for methods like Jacobi or Gauss–Seidel iteration to be discussed later.

#### Forward Difference (One-Sided Approximation):

We can construct a second-order one-sided stencil using forward points:

$$
\frac{d^2 \phi}{dx^2} \approx \frac{\phi_i - 2\phi_{i+1} + \phi_{i+2}}{\Delta x^2}
$$

This biased stencil is **less stable** and not symmetric, making it unsuitable for interior points in elliptic problems. However, it may be used near domain boundaries when backward points are unavailable.  Similarly, a backwards difference approach may be used near domain boundaries where forward points are unavailable.

---

### 1.4.2. Two-Dimensional Laplace Equation

In 2D Cartesian coordinates, the Laplace equation becomes:

$$
\frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} = 0
$$

#### Central Difference Approximation:

We apply the second-order central difference in both $x$ and $y$ directions:

$$
\nabla^2 \phi \approx 
\frac{\phi_{i+1,j} - 2\phi_{i,j} + \phi_{i-1,j}}{\Delta x^2}
+ \frac{\phi_{i,j+1} - 2\phi_{i,j} + \phi_{i,j-1}}{\Delta y^2}
$$

For a **uniform grid** ($\Delta x = \Delta y$), this simplifies to:

$$
\phi_{i,j} = \frac{1}{4}(\phi_{i+1,j} + \phi_{i-1,j} + \phi_{i,j+1} + \phi_{i,j-1})
$$

This forms the basis for iterative solvers such as:
- Jacobi iteration  
- Gauss–Seidel method  
- Successive Over-Relaxation (SOR)

These algorithms sweep through the grid, updating each interior point using its neighbors.

#### Forward Difference in 2D:

Although forward-biased approximations can be constructed in 2D, they are:
- Asymmetric  
- Less accurate unless more points are used  
- Poorly suited for elliptic PDEs like Laplace’s equation

As a result, **forward differencing is not recommended** for use in 2D potential flow problems.

---

**Note:** These finite difference approximations turn the continuous Laplace equation into a system of algebraic equations. Solving these iteratively on a grid yields an approximate solution to the original PDE.

The actual solvers are implemented in the accompanying Python notebook.

---

## 1.5. Solution Methods for Laplace’s Equation

Once the Laplace equation has been discretized using finite difference methods, we must solve the resulting system of equations. These are typically large, sparse, linear systems of the form:

$$
\mathbf{A} \phi = \mathbf{b}
$$

For Laplace’s equation in 1D or 2D, we generally solve this system using **iterative solvers**. Below are the most common methods.

---

### 1.5.1. Jacobi Iteration

The **Jacobi method** is the simplest iterative solver. It updates all grid points **simultaneously** using only values from the previous iteration.

For a 2D uniform grid:

$$
\phi_{i,j}^{(n+1)} = \frac{1}{4} \left( 
\phi_{i+1,j}^{(n)} + \phi_{i-1,j}^{(n)} + \phi_{i,j+1}^{(n)} + \phi_{i,j-1}^{(n)}
\right)
$$

Key characteristics:
- Easy to implement
- Requires two full arrays (old and new)
- **Slow to converge**, especially for large grids

We will use **Jacobi iteration** in both our 1D and 2D examples to keep the implementation clear and easy to follow. Later notebooks may revisit faster solvers.

---

### 1.5.2. Gauss–Seidel Method

The **Gauss–Seidel method** improves on Jacobi by using the most recently updated values as soon as they become available:

$$
\phi_{i,j}^{(n+1)} = \frac{1}{4} \left( 
\phi_{i+1,j}^{(n)} + \phi_{i-1,j}^{(n+1)} + \phi_{i,j+1}^{(n)} + \phi_{i,j-1}^{(n+1)}
\right)
$$

Key characteristics:
- Faster convergence than Jacobi
- Only requires one array (in-place updates)
- Slightly more complex to implement

Although we won’t implement Gauss–Seidel in this notebook, it’s useful to understand how it improves performance.

---

### 1.5.3. Successive Over-Relaxation (SOR)

**SOR** builds on Gauss–Seidel by introducing a **relaxation factor** $\omega$ to accelerate convergence:

$$
\phi_{i,j}^{(n+1)} = (1 - \omega)\phi_{i,j}^{(n)} + \frac{\omega}{4} \left( 
\phi_{i+1,j}^{(n)} + \phi_{i-1,j}^{(n+1)} + \phi_{i,j+1}^{(n)} + \phi_{i,j-1}^{(n+1)}
\right)
$$

When $\omega = 1$, this reduces to Gauss–Seidel. For $1 < \omega < 2$, convergence can be **significantly faster** if tuned correctly.

Key characteristics:
- Requires tuning of $\omega$ for optimal performance
- Fastest convergence among these three for elliptic problems
- Common in production-grade codes

We will mention SOR later when discussing performance, but it is not implemented here.

---

### 1.5.4. Direct Solvers (e.g., LU Decomposition)

For small systems, we can use **direct solvers** like Gaussian elimination or LU decomposition. These solve $\mathbf{A}\phi = \mathbf{b}$ exactly in a finite number of steps.

However:
- They are **not efficient** for large systems
- Matrix $\mathbf{A}$ is typically sparse and structured, which iterative solvers exploit better

For this reason, we avoid direct methods in favor of iterative approaches.

---

### 1.5.5. Multigrid Methods (Advanced)

**Multigrid solvers** accelerate convergence by solving the system at multiple spatial resolutions. They are extremely efficient for elliptic PDEs but require more setup and advanced implementation.

They are out of scope for this course, but worth exploring in research or production contexts.

---

### Summary: Which Solver Are We Using?

| Method           | Convergence | Memory   | 
|------------------|-------------|----------|
| **_Jacobi_**     | _Slow_      | _2 arrays_ |
| Gauss–Seidel     | Moderate    | 1 array  |
| SOR              | Fast        | 1 array  |
| Direct (LU, etc) | Exact       | High     | 
| Multigrid        | Very Fast   | Varies   | 

In this module, we use **Jacobi iteration** due to its simplicity and transparency, especially for illustrating how finite difference solvers operate. More advanced solvers will be discussed later as performance becomes more important.

---

## 1.6. Python Implementation Note

In the next notebook (`1-Solve-Potential_Flow.ipynb`), we implement a 2D finite difference solver for Laplace’s equation over a simple rectangular domain with Dirichlet boundary conditions.

---

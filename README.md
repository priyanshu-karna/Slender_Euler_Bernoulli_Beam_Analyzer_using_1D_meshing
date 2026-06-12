# Euler-Bernoulli Beam FEM Solver

A Python implementation of the **Euler-Bernoulli Beam Finite Element Method (FEM)** using the **Direct Stiffness Method (DSM)**.

This solver can analyze beam structures with multiple point loads acting at arbitrary locations and supports several common boundary conditions. It computes nodal displacements, rotations, support reactions, and generates smooth deflection, shear force, and bending moment diagrams.

---

## Features

✅ Euler-Bernoulli beam formulation

✅ Direct Stiffness Method (DSM)

✅ Arbitrary number of finite elements

✅ Point loads at any location along the beam

✅ Consistent equivalent nodal load vector formulation

✅ Multiple support conditions:

* Simply Supported
* Cantilever
* Fixed-Fixed

✅ Support reaction calculation

✅ Cubic Hermite interpolation for smooth deflection curves

✅ Shear Force Diagram (SFD)

✅ Bending Moment Diagram (BMD)

✅ Object-oriented Python implementation

---

## Theory

The beam element uses the classical Euler-Bernoulli assumption:

* Plane sections remain plane after deformation.
* Plane sections remain perpendicular to the neutral axis.
* Shear deformation is neglected.
* Small deformation theory is assumed.

Each node possesses two degrees of freedom:

[
{d}
===

[w,\theta]^T
]

where:

* (w) = transverse displacement
* (\theta) = rotation

Each beam element therefore contains four degrees of freedom:

[
{d_e}
=====

[w_1,\theta_1,w_2,\theta_2]^T
]

The element stiffness matrix is:

[
K_e=
\frac{EI}{L_e^3}
\begin{bmatrix}
12 & 6L_e & -12 & 6L_e \
6L_e & 4L_e^2 & -6L_e & 2L_e^2 \
-12 & -6L_e & 12 & -6L_e \
6L_e & 2L_e^2 & -6L_e & 4L_e^2
\end{bmatrix}
]

where:

* (E) = Young's Modulus
* (I) = Second Moment of Area
* (L_e) = Element Length

---

## Point Load Handling

The solver supports point loads applied at arbitrary locations.

### Case 1: Load at Node

If the load coincides with a mesh node, the force is directly applied to the corresponding displacement degree of freedom.

### Case 2: Load Between Nodes

If the load lies inside an element, it is converted into an equivalent nodal force vector using beam shape functions:

[
F_e=
\begin{bmatrix}
F_1 \
M_1 \
F_2 \
M_2
\end{bmatrix}
]

This allows accurate modelling of loads without requiring mesh nodes at load locations.

---

## Boundary Conditions

### Simply Supported

* Vertical displacement fixed at both ends
* Rotations remain free

### Cantilever

* Displacement and rotation fixed at left end

### Fixed-Fixed

* Displacement and rotation fixed at both ends

---

## Outputs

The solver computes:

### Nodal Displacements

Vertical displacement at every node.

### Nodal Rotations

Rotation at every node.

### Support Reactions

Reaction forces and moments are obtained from:

[
R = Kd - F
]

### Deflection Curve

Smooth beam deformation using cubic Hermite interpolation.

### Shear Force Diagram (SFD)

Element shear forces are recovered from nodal displacements.

### Bending Moment Diagram (BMD)

Element end moments are computed directly from beam theory.

---

## Example

```python
from beam_solver import Beam

beam = Beam(
    L=6,
    E=210e9,
    I=1e-6,
    n_elem=250,
    bc_type="simply_supported"
)

beam.add_load(P=1000, pos=2)
beam.add_load(P=1000, pos=4)

d = beam.solve()

beam.def_plotting(d)

beam.get_reactions(d)

beam.plot_sfd_bmd(d)
```

---


Generated Results:

* Deflection Curve
* Support Reactions
* Shear Force Diagram
* Bending Moment Diagram


## Future Improvements

Planned extensions include:

* Distributed loads
* Patch loads
* Variable cross-section beams
* Timoshenko beam elements
* Multiple load cases
* Modal analysis
* Dynamic analysis
* 2D Frame Elements
* General FEM Framework

---

## Author

Priyanshu Karna

Civil Engineering Student

Interests:

* Finite Element Methods
* Computational Mechanics
* Numerical Methods
* Structural Analysis
* Scientific Computing
* Optimization



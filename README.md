# Euler-Bernoulli Beam FEM Solver

A Python implementation of the **Euler-Bernoulli Beam Finite Element Method (FEM)** using the **Direct Stiffness Method (DSM)**.

This solver can analyze beam structures with multiple point loads acting at arbitrary locations and supports several common boundary conditions. It computes nodal displacements, rotations, support reactions, and generates smooth deflection, shear force, and bending moment diagrams.

---

## Features

✅ Euler-Bernoulli beam formulation

✅ Direct Stiffness Method (DSM)

✅ Arbitrary number of finite elements

✅ Point loads at any location along the beam

✅ Distributed loads (Both UDL and UVL) converted to nodal load and moments using 3 point gauss quadrature

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
beam1.add_distributed_load(w1=0, w2=-1000, start_pos=0, end_pos=6)
beam1.add_distributed_load(w1=-1000, w2=-1000, start_pos=4.6, end_pos=5.87)

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



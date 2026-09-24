# Al slab test description, 3x3 cells

Cross sections from Al-6061, subset of the 3 middle groups (2-4) from the
5g_5x test. Small test problem with 3 energy and 3 spatial cells.

Groups are chosen so that group 1 (of this subset, original group 2)
downscatters into groups 2 and 3, while groups 2 and 3 (original groups 3
and 4) have only within-group scattering.

Scattering is considered to be isotropic.

Energy mesh:
- Middle 3 groups of the 5g_5x mesh (groups 24-26 from 36g)
- [7.7426e-6, 4.6415e-6, 2.7825e-6, 1.668e-6] MeV

Cross sections:
- Same source data as 5g_5x (`radiant.jl`), restricted to groups 2-4

Spatial mesh:
- 3 cells, 5.0e-7 cm each

Angular Quadrature:
- s6 double Gauss-Legendre

Boundary conditions:
- Isotropic incoming flux on both sides, same values as 5g_5x groups 2-4,
  truncated to the 3-cell spatial mesh

Independent source:
- 1 in all cells, all groups, all angles

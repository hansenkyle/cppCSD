# Al slab test description, 5x5 cells

Cross sections from Al-6061. Small test problem with 5 energy and 5 spatial cells.

Scattering is considered to be isotropic.

Energy mesh:
- 1 eV - 12.916 eV, log spaced (groups 23-27 from 36g)
- [1.29155e-5, 7.7426e-6, 4.6415e-6, 2.7825e-6, 1.668e-6, 1.00e-6] MeV

Cross sections:
- Generated using `radiant.jl`, using inelastic + elastic collision, Bremsstrahlung, including knock-on production
- Full scattering matrix (relavent for higher-energy groups)

Spatial mesh:
- 5 cells, 5.0e-7 cm each

Angular Quadrature:
- s6 double Gauss-Legendre

Boundary conditions:
- Isotropic incoming flux on both sides
    - Left: psi(E) = 10(1+exp(-3e5 (x-1e-6))) [(cm2sMeVstr)^-1]
    - Right psi(E) = 10                       [(cm2sMeVstr)^-1]

Independent source:
- 1 in all cells, all groups, all angles

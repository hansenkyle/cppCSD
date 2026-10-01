# MCDC Comparison test description

Modified 24-group Al-6061 cross sections:
- first 12 groups only. 0.2236-5.0 MeV, log-spaced energy groups
- diagonal of scattering matrix only (within-group scattering)


CSD implemented in MCDC using excitation reaction with:
- sigma = S/Delta E, where
- Delta E is uniformly (linearly) divided in each energy group, with 2^(10, 11, 12, 13) divisions per group

large delta E: cheap compuation, poorly approximated by CSD
small delta E: expensive, better approximated by CSD.

delta will be parameterized by group width:

All flux units in /(cm^2 s MeV)
Boundary conditions (isotropic):
- Left:
    - psi_in = 1,000, (first energy group)
    -          1000(E-2.9789)/(3.8593-2.9789), second energy group
    -          0    , otherwise
- Right:
    - psi_in = 1, all groups

External source:
- q = 1 /(cm^2 s MeV), all energy groups

MCDC tally: 100-cell grid
deterministic: 40 spatial cells, S32 double Gauss-Legendre quadrature

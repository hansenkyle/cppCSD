# MCDC Comparison test description

Modified 24-group Al-6061 cross sections:
- first 12 groups only. 0.2236-5.0 MeV, log-spaced energy groups
- diagonal of scattering matrix only (within-group scattering)


CSD implemented in MCDC using excitation reaction with Delta_E = delta.
MCDC run using varying values for delta.

large delta: cheap compuation, poorly approximated by CSD
small delta: expensive, better approximated by CSD.

delta will be parameterized by group width:

Delta E_g / C

C will range from 20 to ...


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
deterministic: 10-60 spatial cells, S16 double Gauss-Legendre quadrature

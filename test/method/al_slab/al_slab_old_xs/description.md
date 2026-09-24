# Al slab test description

NOTE: the scattering cross sections used in this test were not generated properly; the test does not represent real Al.

0.465-cm Al-6061 Slab, exposed to isotropic electron flux with spectrum from Upper Van Allen belt.

Scattering is considered to be isotropic.

Spatial mesh:
- 0.456 cm total thickness (from https://www.nasa.gov/wp-content/uploads/2024/02/6.soa-structures-2023.pdf?emrc=65d0f77aa93ee table 6-15)
- 100 cells

Energy mesh:
- 0.01 eV - 1 MeV
- 36 log-spaced groups

Cross sections:
- Generated using RADIANT.jl, with reactions:
    - Elastic collision
    - Inleastic collision, including knock-on production
    - Bremsstrahlung
- Stopping power is negative in some groups; these are set to zero.
- Group-to-group scattering ignored when sigma(g-> g')/sigma_t < 1e-8

Boundary conditions:
- Isotropic incoming flux on one side
- Energy spectrum from from "Van Allen Probes show that the inner radiation zone contains no MeV electrons: ECT/MagEIS data" (Fennel, Claudepierre, et al. 2015), figure 2 is approximately used:
    - psi = 1.5e+5 exp(-0.01 E[keV]) 1/(cm2 s ster keV), E<1000 keV. 0, E>1000 keV.

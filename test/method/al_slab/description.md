# Al shield test description

Al-6061 slab shielding electronics from outer-belt electrons. Cross sections
are in `xs_data/`, for N_g = 6, 12, 24, 36, 48.

Cross sections (`scripts/generate_xs.jl`, Radiant.jl):
- Al-6061: Al 0.9805, Mg 0.0100, Si 0.0060, Fe 0.0035 (weight fractions), 2.70 g/cm^3
- 5 MeV - 10 keV, log-spaced groups
- Inelastic (Moller, primary "S" only -- no knock-on production, knock-on
  energy is deposited locally instead), elastic, bremsstrahlung (as energy loss only;
  photons are not transported)
- Legendre order 1: transport-corrected isotropic scattering
  (sigma_t - sigma_1, sigma_s0 - sigma_1)

Incident spectrum (isotropic):
- Claudepierre et al. (2021), "The Magnetic Electron Ion Spectrometer: A
  Review of On-Orbit Sensor Performance, Data, Operations, and Science",
  Space Sci. Rev. 217, 80, doi:10.1007/s11214-021-00855-2, Fig. 17
  (background-corrected, L = 5.46, 2013-10-18 19:02:30 UT). Values read off
  the plot, so accurate to about a factor of 1.5.

2.5 mm Al in 10 cells. (The CSDA range is about
2 mm at 1 MeV and about 11 mm at 5 MeV, per NIST ESTAR)

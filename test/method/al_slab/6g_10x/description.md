# 0.25 cm Al-6061 test

Spatial mesh:
- 10 uniform spatial cells, totaling 0.25 cm

Angular quadrature:
- S6 double Gauss-Legendre

Energy group structure:
- 6 log-spaced groups, 5 MeV - 10 keV

Cross sections:
- Material:
    - 0.9805 Al
    - 0.0100 Mg
    - 0.0060 Si
    - 0.0035 Fe
- Cross sections generated with RADIANT.jl:
    - Elastic scatter, inelastic scatter, Bremsstrahlung
    - No knock-on electron production
    - Isotropic scattering cross sections, includes downscatter

Boundary conditions:
- Left: Electron energy spectrum in upper van allen belt, from fig. 17 of Claudepierre et al. 2021, The Magnetic Electron Ion Spectrometer: A Review of On-Orbit Sensor Performance, Data, Operations, and Science
- peak ~1e8 at 10 keV
- Right: uniform spectrum
    - 5e+2 /(cm2 s ster MeV)
- Both sides are isotropic

Source:
- Uniform, isotropic source, 



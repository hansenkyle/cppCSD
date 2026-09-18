# Scaled 3g_3x test description

Same structure as 3g_3x (3 groups, 3 cells, group 1 downscatters into
groups 2 and 3, groups 2 and 3 within-group scattering only), but with all
values rescaled to be order-1 so each quantity is easy to eyeball and tell
apart. Not physical cross sections -- a synthetic case for exercising the
method.

Trends preserved from the Al-6061 data:
- Total xs increases with energy (group 1, the highest-energy group, has
  the largest sigma_t).
- Stopping power peaks at the middle group and falls off to either side
  (Bragg-curve shape), for both group-average and group-boundary values.
- Boundary condition magnitude increases with group number.

Spatial mesh:
- 3 cells, width 0.5 each

Energy mesh:
- 3 groups, boundaries [4.00, 3.00, 2.00, 1.00]

Angular Quadrature:
- s6 double Gauss-Legendre (unchanged -- dimensionless)

Boundary conditions:
- Isotropic incoming flux on both sides, magnitude ~1, same per-group shape
  as 3g_3x

Independent source:
- 1 in all cells, all groups, all angles

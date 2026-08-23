# Input Deck Format

An input deck fully describes one problem: spatial mesh, region
materials, energy mesh, per-material cross sections, and convergence
criteria. Read by `Parser::read()` ([../src/parser.h](../src/parser.h),
[../src/parser.cpp](../src/parser.cpp)) into an `InputDeck`.

## Format

YAML. Sizes below are given in terms of:

```
num_cells  = spatial_mesh.size() - 1
num_groups = energy_mesh.size() - 1
```

```yaml
spatial_mesh: [0.0, 1.0, 2.0, 3.0]     # ascending cell boundaries [num_cells + 1 floats]

regions:
  materials: [water, water, lead]      # material name per cell [num_cells strings]

energy_mesh: [5.0, 1.0, 0.5, 0.0]      # MeV, descending group boundaries, ending at 0 [num_groups + 1 floats]

materials:                             # one entry per name used in regions.materials
  water:
    sigma_t: [1.2, 1.0, 0.8]           # total macroscopic xs [num_groups floats]
    sigma_s: [1.1, 0.9, 0.7]           # isotropic scattering xs [num_groups floats]
    stopping_power:
      group_average: [2.0, 1.8, 1.5]           # [num_groups floats]
      group_boundary: [2.2, 1.9, 1.6, 1.3]     # evaluated at group boundaries [num_groups + 1 floats]
  lead:
    sigma_t: [3.2, 3.0, 2.8]
    sigma_s: [2.1, 1.9, 1.7]
    stopping_power:
      group_average: [5.0, 4.8, 4.5]
      group_boundary: [5.2, 4.9, 4.6, 4.3]

convergence:
  max_iters: 200                       # [1 integer]
  epsilon: 1.0e-8                      # [1 float]
```

This is [../test/unit/data/sample_input.yaml](../test/unit/data/sample_input.yaml).
To generate the `materials:` block from cross-section physics rather than
writing it by hand, see [../scripts/generate_xs.jl](../scripts/generate_xs.jl).

## Notes

- `spatial_mesh` must be strictly ascending -- no repeats, no descending
  values. `energy_mesh` must be strictly descending, entirely
  non-negative, and end at exactly `0` -- particles slow down from the
  top group to zero energy.
- `spatial_mesh` and `energy_mesh` are read into a single `Mesh`
  (`deck.mesh`, see [../src/mesh.h](../src/mesh.h)), which also derives
  cell/group widths (`dx`, `dE`) and midpoints (`x_center`, `E_center`).
- `regions.materials` has one entry per cell, not per boundary, and
  every name it uses must appear under `materials:`.
- Group index 0 is the highest energy group. `stopping_power.group_boundary`
  is the only array sized `num_groups + 1`; everything else per material
  is sized `num_groups`.

`Parser::read()` throws `std::runtime_error` naming the offending field
and the size or order it expected.

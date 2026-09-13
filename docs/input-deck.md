# Input Deck Format

An input deck fully describes one problem: spatial mesh, region
materials, energy mesh, per-material cross sections, angular
quadrature, boundary conditions, and convergence criteria. Read by
`InputDeck::read()`
([../src/input_deck.h](../src/input_deck.h),
[../src/input_deck.cpp](../src/input_deck.cpp)).

## Format

YAML. Sizes below are given in terms of:

```
num_cells     = spatial_mesh.size() - 1
num_groups    = energy_mesh.size() - 1
num_ordinates = angular_quadrature.mu.size()
```

```yaml
spatial_mesh: [0.0, 1.0, 2.0, 3.0]     # ascending cell boundaries [num_cells + 1 floats]

regions:
  materials: [water, water, lead]      # material name per cell [num_cells strings]

energy_mesh: [5.0, 1.0, 0.5, 0.0]      # MeV, descending group boundaries, ending at 0 [num_groups + 1 floats]

materials:                             # one entry per name used in regions.materials
  water:
    sigma_t: [1.2, 1.0, 0.8]           # total macroscopic xs [num_groups floats]
    # Sparse group-to-group scattering matrix: only nonzero (from, to)
    # entries need to appear. from/to are 1-indexed, in [1, num_groups].
    scattering:
      - {from: 1, to: 1, value: 1.1}
      - {from: 1, to: 2, value: 0.05}
      - {from: 2, to: 2, value: 0.9}
      - {from: 2, to: 3, value: 0.03}
      - {from: 3, to: 3, value: 0.7}
    stopping_power:
      group_average: [2.0, 1.8, 1.5]           # [num_groups floats]
      group_boundary: [2.2, 1.9, 1.6, 1.3]     # evaluated at group boundaries [num_groups + 1 floats]
  lead:
    sigma_t: [3.2, 3.0, 2.8]
    scattering:
      - {from: 1, to: 1, value: 2.1}
      - {from: 1, to: 2, value: 0.1}
      - {from: 2, to: 2, value: 1.9}
      - {from: 2, to: 3, value: 0.08}
      - {from: 3, to: 3, value: 1.7}
    stopping_power:
      group_average: [5.0, 4.8, 4.5]
      group_boundary: [5.2, 4.9, 4.6, 4.3]

angular_quadrature:
  mu: [-0.9, -0.3, 0.3, 0.9]           # direction cosines, strictly ascending [num_ordinates floats]
  w: [0.5, 0.5, 0.5, 0.5]              # quadrature weights, normalized to sum to 2 [num_ordinates floats]

boundary_conditions:                   # incoming angular flux at each spatial boundary
  left:                                # [num_ordinates][num_groups], each a down/up pair
    down: [[0, 1, 2], [10, 11, 12], [20, 21, 22], [30, 31, 32]]
    up: [[0.5, 1.5, 2.5], [10.5, 11.5, 12.5], [20.5, 21.5, 22.5], [30.5, 31.5, 32.5]]
  right:
    down: [[0, 1, 2], [100, 101, 102], [200, 201, 202], [300, 301, 302]]
    up: [[0.5, 1.5, 2.5], [100.5, 101.5, 102.5], [200.5, 201.5, 202.5], [300.5, 301.5, 302.5]]

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
- A material's `scattering` entries are validated against `num_groups`
  (`from`/`to` in range, non-negative `value`, no duplicate `(from, to)`
  pair) and against each other -- a duplicate entry is almost certainly a
  typo, so it's rejected rather than summed. `CrossSection::scattering` (see
  [../src/cross_section.h](../src/cross_section.h)) is indexed `[cell]`
  rather than `[group][cell]`, since each cell just gets a copy of its
  material's sparse entry list.
- `angular_quadrature.mu` must be strictly ascending. `angular_quadrature.w`
  is rescaled so its entries sum to exactly `2`, regardless of what's
  written in the file; the applied scale factor is logged.
- `boundary_conditions.left`/`right` each hold a down/up pair (`DownUp`,
  see [../src/input_deck.h](../src/input_deck.h)) per ordinate and group --
  down is the lower-energy edge, up the higher-energy edge, matching the
  same convention used elsewhere (e.g. `Corner` in
  [../src/fe_space.h](../src/fe_space.h)).

`InputDeck::read()` throws `std::runtime_error` naming the offending
field and the size or order it expected, catches it internally, and
returns `1`; a `Mesh` construction failure (see
[../src/mesh.h](../src/mesh.h)) is caught the same way.

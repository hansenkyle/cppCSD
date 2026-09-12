# =====================================================================
# generate_xs.jl — Macroscopic multigroup cross sections + stopping
#                   powers for cppCSD, via Radiant.jl
#
# Builds electron cross-section data from one or more Radiant Materials
# and writes them to a CSV intended for later automated conversion to
# cppCSD's YAML input format (that conversion is a separate script, not
# this one). Layout, one block per material after a shared header:
#
#   # metadata comment lines (generator, timestamp, particle,
#   # group structure, material list)
#   energy_mesh_MeV
#   <Ng+1 ascending group boundary values>
#   material,<name>
#   composition
#   <element>,<weight fraction>
#   ...
#   stopping_power_boundary_MeV_cm
#   <Ng+1 values, at group boundaries>
#   group,sigma_t_cm-1,stopping_power_average_MeV_cm
#   <Ng rows, one per group>
#   scattering_matrix_cm-1 (rows=from-group, cols=to-group, l=0 moment)
#   from\to,1,2,...,Ng
#   <Ng rows, one per from-group; zero entries left blank so the
#   sparsity structure of the matrix stays visible>
#
# The scattering matrix is always written in full (Ng x Ng), not just
# its diagonal -- it is not assumed to be sparse or near-diagonal, even
# though in practice only in-group scattering plus a small downscatter
# band tends to be nonzero.
# =====================================================================
using Radiant
using Printf
using Dates

# --------------------------- USER INPUT -----------------------------
# One entry per material. `name` becomes the CSV material key and must
# eventually match the material names used in the deck's
# `regions.materials` list once converted to YAML.
materials_input = [
    (name = "water", density = 1.0, elements = ["H", "O"], wfractions = [0.111894, 0.888106]),
]

Ng         = 20                    # number of energy groups (ignored if custom_energy_bounds is set below)
E_max      = 10.0                  # midpoint energy of highest group [MeV] (ignored if custom_energy_bounds is set)
E_cut      = 0.001                 # cutoff energy [MeV] (ignored if custom_energy_bounds is set)
group_type = "log"                 # "log" or "linear" (ignored if custom_energy_bounds is set)

# --- Custom energy group boundaries (optional) -------------------------
# To use group boundaries that aren't a plain log/linear sweep, list them
# here explicitly instead: Ng+1 boundary energies [MeV], either ascending
# or descending (Radiant accepts either). This OVERRIDES Ng/E_max/E_cut/
# group_type above -- Ng is derived from this vector's length instead.
# Leave this empty ([]) to keep using the log/linear structure above.
custom_energy_bounds = Float64[]   # e.g. Float64[0.001, 0.01, 0.1, 1.0, 10.0]
# -------------------------------------------------------------------------

legendre_order = 1                 # Legendre truncation order used internally by
                                    # Radiant's elastic-scattering decomposition;
                                    # only the l=0 moment is written out.

output_name = "water_20g.csv"  # output filename, written under scripts/xs_data/
# ----------------------------------------------------------------------

output_dir = joinpath(@__DIR__, "xs_data")
mkpath(output_dir)
outfile = joinpath(output_dir, output_name)

particle = Radiant.Electron()

interaction_list = [
    Radiant.Inelastic_Collision(),  # Moller collisional energy loss + catastrophic delta rays
    Radiant.Elastic_Collision(),    # Mott elastic scattering (large-angle part, AFP-decomposed)
    Radiant.Bremsstrahlung(),       # radiative energy loss
]

# --- Build materials ---
material_list = Radiant.Material[]
for m in materials_input
    @assert length(m.elements) == length(m.wfractions) "elements and wfractions must match in length for material '$(m.name)'"
    wfractions_norm = m.wfractions ./ sum(m.wfractions)

    mat = Radiant.Material(m.name)
    mat.set_density(m.density)
    for (el, f) in zip(m.elements, wfractions_norm)
        mat.add_element(el, f)
    end
    push!(material_list, mat)
end

# --- Build cross sections ---
using_custom_bounds = !isempty(custom_energy_bounds)

cs = Radiant.Cross_Sections()
cs.set_source("physics-models")
cs.set_materials(material_list)
cs.set_particles([particle])
if using_custom_bounds
    cs.set_group_structure(custom_energy_bounds)
    Ng = length(custom_energy_bounds) - 1
else
    cs.set_group_structure(group_type, Ng, E_max, E_cut)
end
cs.set_interactions(interaction_list)
cs.set_legendre_order(legendre_order)
cs.build()

# --- Pull data ---
# Radiant numbers groups from the highest energy down to the lowest (group 1
# = E_max). cppCSD's parser requires strictly descending energy boundaries
# (highest → lowest, ending at 0), so we use the data directly without reversal.
Eb = cs.get_energy_boundaries(particle)                 # MeV, size Ng+1, descending (E_max → E_cut)
Σt = cs.get_total(particle)                             # [Ng, Nmat], cm^-1
Σs_moments = cs.get_scattering(particle, particle, legendre_order) # [Nmat, Ng, Ng, legendre_order+1]
S  = cs.get_stopping_powers(particle)                   # [Ng, Nmat], MeV/cm
Sb = cs.get_boundary_stopping_powers(particle)          # [Ng+1, Nmat], MeV/cm

# --- Write CSV ---
format_floatrow(v) = join([@sprintf "%.6e" x for x in v], ",")

open(outfile, "w") do io
    println(io, "# Generated by scripts/generate_xs.jl (Radiant.jl)")
    println(io, "# Generated: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println(io, "# Particle: electron")
    if using_custom_bounds
        println(io, "# Group structure: custom boundaries, Ng=", Ng)
    else
        println(io, "# Group structure: ", group_type, ", Ng=", Ng,
                     ", E_max=", E_max, " MeV, E_cut=", E_cut, " MeV")
    end
    println(io, "# Materials: ", join([m.name for m in materials_input], ", "))
    println(io, "#")
    println(io, "energy_mesh_MeV")
    println(io, format_floatrow(Eb))

    for (imat, m) in enumerate(materials_input)
        sigma_t = Σt[:, imat]
        stopping_power_average = S[:, imat]
        stopping_power_boundary = Sb[:, imat]
        wfractions_norm = m.wfractions ./ sum(m.wfractions)

        # Full l=0 scattering matrix for this material (Radiant already provides
        # data in descending group order: group 1 = highest energy).
        Σs_full = [Σs_moments[imat, f, t, 1] for f in 1:Ng, t in 1:Ng]

        println(io, "#")
        println(io, "material,", m.name)
        println(io, "composition")
        for (el, f) in zip(m.elements, wfractions_norm)
            println(io, el, ",", @sprintf("%.6e", f))
        end
        println(io, "stopping_power_boundary_MeV_cm")
        println(io, format_floatrow(stopping_power_boundary))
        println(io, "group,sigma_t_cm-1,stopping_power_average_MeV_cm")
        for g in 1:Ng
            println(io, g, ",", @sprintf("%.6e", sigma_t[g]), ",", @sprintf("%.6e", stopping_power_average[g]))
        end
        println(io, "scattering_matrix_cm-1 (rows=from-group, cols=to-group, l=0 moment)")
        println(io, "from\\to,", join(1:Ng, ","))
        for f in 1:Ng
            row = [Σs_full[f, t] == 0.0 ? "" : @sprintf("%.6e", Σs_full[f, t]) for t in 1:Ng]
            println(io, f, ",", join(row, ","))
        end
    end
end

println("Wrote cross sections for $(length(materials_input)) material(s) to $outfile")

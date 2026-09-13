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
#   <Ng+1 descending group boundary values>
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

# Legendre truncation order Radiant uses internally to decompose elastic
# scattering; only the l=0 (isotropic/total) moment is ever written out,
# but a higher internal order improves that moment's accuracy. Not a user
# knob -- there's no reason to trade accuracy for speed here.
const RADIANT_LEGENDRE_ORDER = 7

# --------------------------- USER INPUT -----------------------------
# One entry per material. `name` becomes the CSV material key and must
# eventually match the material names used in the deck's
# `regions.materials` list once converted to YAML.
materials_input = [
    (name = "aluminum", density = 2.7, elements = ["Al"], wfractions = [1.0]),
]

# Energy group boundaries [MeV], strictly descending (highest energy
# first, ending at the cutoff) -- cppCSD's parser requires this same
# descending order, so no reversal is ever needed below. Ng is derived
# from this vector's length.
energy_bounds = Float64[100, 10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001]

# Whether to include knock-on/delta-ray PRODUCTION (Moller's "P" interaction
# type, i.e. the ejected secondary electron) in the electron-electron
# scattering matrix, alongside the primary electron's own redirection ("S").
# cppCSD's solver has no production/source-term mechanism yet -- Sigma_t
# only ever counts "S" -- so turning this on without matching solver
# support breaks particle conservation (a group's Sigma_s row sum can
# exceed its Sigma_t). Leave false until delta-ray production is added to
# the solver.
include_knockon_production = false

output_name = "al_27gcc.csv"  # output filename, written under scripts/xs_data/
# ----------------------------------------------------------------------

Ng = length(energy_bounds) - 1

output_dir = joinpath(@__DIR__, "xs_data")
mkpath(output_dir)
outfile = joinpath(output_dir, output_name)

particle = Radiant.Electron()

# Radiant's native soft/catastrophic split (scattering_model = "BFP", the
# default for every interaction below) already keeps particle conservation
# for the "S" (scattering) type on its own: the catastrophic cutoff is
# derived from the energy group structure itself, and only "S" feeds both
# Sigma_t and Sigma_s. "P" (production) is the one type Radiant excludes
# from Sigma_t by design (see include_knockon_production above), since it
# describes a newly created particle rather than redirection of the one
# being tracked.
electron_electron_types = include_knockon_production ? ["S", "P"] : ["S"]
inelastic_collision = Radiant.Inelastic_Collision()
inelastic_collision.set_interaction_types(Dict(
    (Radiant.Positron, Radiant.Positron) => ["S"],
    (Radiant.Positron, Radiant.Electron) => ["P"],
    (Radiant.Electron, Radiant.Electron) => electron_electron_types,
    (Radiant.Proton, Radiant.Proton) => ["S"],
    (Radiant.Proton, Radiant.Electron) => ["P"],
    (Radiant.Alpha, Radiant.Alpha) => ["S"],
    (Radiant.Alpha, Radiant.Electron) => ["P"],
))

interaction_list = [
    inelastic_collision,          # Moller collisional energy loss (+ knock-on production if toggled on)
    Radiant.Elastic_Collision(),  # Mott elastic scattering (large-angle part, AFP-decomposed)
    Radiant.Bremsstrahlung(),     # radiative energy loss
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
cs = Radiant.Cross_Sections()
cs.set_source("physics-models")
cs.set_materials(material_list)
cs.set_particles([particle])
cs.set_group_structure(energy_bounds)
cs.set_interactions(interaction_list)
cs.set_legendre_order(RADIANT_LEGENDRE_ORDER)
cs.build()

# --- Pull data ---
Eb = cs.get_energy_boundaries(particle)                 # MeV, size Ng+1, descending
Σt = cs.get_total(particle)                             # [Ng, Nmat], cm^-1
Σs_moments = cs.get_scattering(particle, particle, RADIANT_LEGENDRE_ORDER) # [Nmat, Ng, Ng, RADIANT_LEGENDRE_ORDER+1]
S  = cs.get_stopping_powers(particle)                   # [Ng, Nmat], MeV/cm
Sb = cs.get_boundary_stopping_powers(particle)          # [Ng+1, Nmat], MeV/cm

# --- Write CSV ---
format_floatrow(v) = join([@sprintf "%.6e" x for x in v], ",")

open(outfile, "w") do io
    println(io, "# Generated by scripts/generate_xs.jl (Radiant.jl)")
    println(io, "# Generated: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println(io, "# Particle: electron")
    println(io, "# Group structure: custom boundaries, Ng=", Ng)
    println(io, "# Materials: ", join([m.name for m in materials_input], ", "))
    println(io, "#")
    println(io, "energy_mesh_MeV")
    println(io, format_floatrow(Eb))

    for (imat, m) in enumerate(materials_input)
        sigma_t = Σt[:, imat]
        stopping_power_average = S[:, imat]
        stopping_power_boundary = Sb[:, imat]
        wfractions_norm = m.wfractions ./ sum(m.wfractions)

        # Full l=0 scattering matrix for this material.
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

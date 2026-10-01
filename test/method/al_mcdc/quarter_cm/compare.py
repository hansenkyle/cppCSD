import h5py
import ldcsd

import numpy as np
import pathlib

from matplotlib import pyplot as plt



LDCSD_DATA = "al.25cm.pyslab.h5"
<<<<<<< Updated upstream
MCDC_DATA = "al_slab2.5cm_12g_512div.h5"
=======
MCDC_DATA = "al_slab2.5cm_12g_8192div.h5"
>>>>>>> Stashed changes



# load cpp data into numpy arrays

run = ldcsd.read(f"{pathlib.Path(__file__).parent.resolve()}/al.25cm.pyslab.h5")
scalar = ldcsd.cell_average(run.scalar_flux)

dx = run.input.dx
xi = run.input.x_center

dE = run.input.dE
Eb = run.input.E_boundary

ld_x = np.zeros(2*len(dx))
for i in range(0, len(dx)):
    ld_x[2*i:2*i+2] = [xi[i]-0.5*dx[i], xi[i]+0.5*dx[i]]

ld_E =np.zeros(2*len(dE))
ld_right_spec = np.zeros(2*len(dE))
for g in range(0, len(dE)):
    ld_E[2*g:2*g+2] = [Eb[g], Eb[g+1]]

    ld_right_spec[2*g:2*g+2] = run.scalar_flux[g, -1, :, 1]


# load mcdc data into numpy arrays

with h5py.File(MCDC_DATA, "r") as f:
    tally = f["tallies/flux"]
    z = tally["grid/z"][:]
<<<<<<< Updated upstream
    mc_phi = tally["flux/mean"][:]
    mc_phi_sdev = tally["flux/mean"][:]
=======
    mc_E_flux = tally["grid/energy"][:]
    mc_phi = tally["flux/mean"][:]
    mc_phi_sdev = tally["flux/sdev"][:]
>>>>>>> Stashed changes

    tally = f["tallies/flux_right"]
    mc_E = tally["grid/energy"][:]
    mc_spectrum = tally["flux/mean"][:]
    mc_spectrum_sdev = tally["flux/sdev"][:]


zi = 0.5*(z[1:] + z[:-1])
dz = (z[1:] - z[:-1])
<<<<<<< Updated upstream
mc_phi *= (0.25/dz)
=======
mc_dE_flux = (mc_E_flux[1:] - mc_E_flux[:-1])
# same normalization as the spectrum: per unit length (1/dz), and E_range/dE in energy
mc_phi *= ((mc_E_flux[-1] - mc_E_flux[0])/mc_dE_flux)[:, None]/dz
>>>>>>> Stashed changes

mc_E_center = 0.5*(mc_E[1:] + mc_E[:-1])
mc_dE = (mc_E[1:] - mc_E[:-1])

mc_spectrum*= ((mc_E[-1]-mc_E[0])/mc_dE)

# plot scalar flux for select groups (to compare spatial distribution)


# plot spectrum at right face of slab (to compare energy distribution)

mc_dx = 1e-7
mc_source = 1e3

<<<<<<< Updated upstream
plt.figure()
plt.plot(xi, 2*np.pi*scalar[-1, :])
plt.plot(zi, mc_source*mc_phi[0]*dE[-1])
plt.title("Flux distribution")
plt.yscale('log')
plt.savefig("spatial.png")
=======
# ldcsd groups run high -> low energy; MC/DC's energy grid runs low -> high
groups = list(range(len(dE)))

plt.figure()
for g in groups:
    line, = plt.plot(xi, 2*np.pi*scalar[g, :], label=f"g{g}")
    plt.plot(zi, mc_source*mc_phi[len(dE)-1-g], '--', color=line.get_color())
plt.title("Flux distribution, (solid: LD, dashed: MC/DC)")
plt.legend(fontsize='small', ncol=2)
plt.yscale('log')
plt.xlabel("x (cm)")
plt.ylabel("phi_g")
plt.savefig("256angle_spatial.png")
>>>>>>> Stashed changes



plt.figure()
plt.plot(mc_E_center, (mc_source/mc_dx)*mc_spectrum[1])
plt.plot(ld_E*1e6, 2*np.pi*ld_right_spec)
<<<<<<< Updated upstream
plt.yscale('log')
plt.title("Spectrum, 0.25 cm through slab")
plt.savefig("spectrum.png")

# compute relative error over each group
=======
plt.xlabel("E (eV)")
plt.ylabel("phi(X)")
plt.yscale('log')
plt.title("Spectrum, right face of slab")
plt.savefig("265angle_spectrum.png")


# compute relative error over each group
plt.show()
>>>>>>> Stashed changes

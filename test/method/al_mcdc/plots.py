import ldcsd
import pathlib
import matplotlib.pyplot as plt
import numpy as np

run = ldcsd.read(f"{pathlib.Path(__file__).parent.resolve()}/runs/latest/al_mcdc_12g_10x.h5")
scalar = ldcsd.cell_average(run.scalar_flux)

dx = run.input.dx
xi = run.input.x_center

dE = run.input.dE
Eb = run.input.E_boundary

x = np.zeros(2*len(dx))
for i in range(0, len(dx)):
    x[2*i:2*i+2] = [xi[i]-0.5*dx[i], xi[i]+0.5*dx[i]]

E =np.zeros(2*len(dE))
right_spec = np.zeros(2*len(dE))
for g in range(0, len(dE)):
    E[2*g:2*g+2] = [Eb[g], Eb[g+1]]

    right_spec[2*g:2*g+2] = run.scalar_flux[g, -1, :, 1]

plt.figure()
plt.title("Cell-centered scalar flux, integrated over energy group")
for g in range(0, 12):
    plt.plot(xi, scalar[g]*run.input.dE[g]*4/1000)
# plt.ylim([1e-6, 3e0])
# plt.yscale('log')

plt.figure()
plt.title("Scalar flux energy spectrum, right face")
plt.plot(E*1e6, right_spec*1e-15/2)
plt.yscale('log')
plt.xscale('log')


plt.show()
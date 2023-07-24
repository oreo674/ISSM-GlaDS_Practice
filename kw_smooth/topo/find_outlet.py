"""
Find the lowest hydraulic potential node for smooth KW geometry
to use as the domain outlet
"""

import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt

# Read in datasets
bed = np.loadtxt('bed_DEM.dat')
sfc = np.loadtxt('sfc_DEM.dat')
contour = np.loadtxt('KW_outline_simplified.txt', delimiter=',')

# Generate phi for different assumed floatation fraction
ffs = np.arange(0.5, 1.1, 0.1)

bed_elev = bed[:, -1].reshape(-1, 1)
sfc_elev = sfc[:, -1].reshape(-1, 1)
bed_xy = bed[:, :2]

rhow = 1000
rhoi = 910
g = 9.81

phi_bed = bed_elev*rhow*g
p_ice = (sfc_elev - bed_elev)*rhoi*g
phi = phi_bed + ffs*p_ice

fig, ax = plt.subplots()
# Interpolate
for i in range(len(ffs)):
    phi_interpolant = interpolate.NearestNDInterpolator(bed_xy, phi[:, i])
    phi_contour = phi_interpolant(contour)
    ax.plot(np.arange(0, len(phi_contour)), phi_contour)

    index_phi_min = np.argmin(phi_contour)
    print('Phi min index:', index_phi_min)

plt.show()
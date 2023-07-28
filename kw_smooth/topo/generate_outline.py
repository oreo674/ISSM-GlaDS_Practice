"""
Generate ISSM outline file from Erik's Kaskawulsh outline
"""

import numpy as np
from matplotlib import pyplot as plt
import cmocean
from scipy import interpolate

# Load DEMs and contour file
bed = np.loadtxt('bed_DEM.dat')
contour = np.loadtxt('contour.dat')
sfc = np.loadtxt('sfc_DEM.dat')
thick = sfc[:, -1] - bed[:, -1]
x = bed[:, 0]
y = bed[:, 1]

# Set ELA
ELA = 2260

# Reshape into 2D arrays
xlen = len(x[x==x[0]])
ylen = len(y[y==y[0]])
thick = thick.reshape((ylen, xlen), order='F')
xx = x.reshape((ylen, xlen), order='F')
yy = y.reshape((ylen, xlen), order='F')
# Ignore where thickness < 75 m, only effects plotting
mask = np.abs(thick)<75
thick[mask] = np.nan

# Plot raw data
fig, ax = plt.subplots(figsize=(6, 3.5))
pc = ax.pcolormesh(xx/1e3, yy/1e3, thick, vmin=0, vmax=1000, cmap=cmocean.cm.deep)
ax.plot(contour[:, 0]/1e3, contour[:, 1]/1e3, 'k.-')
labels = np.arange(0, len(contour[:, 0])).astype(str)
cb = fig.colorbar(pc)
cb.set_label('Thickness (m)')
ax.set_aspect('equal')
ax.set_xlabel('Easting (km)')
ax.set_ylabel('Northing (km)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# FIND BRANCHES TO CUT TO SIMPLIFY THE DOMAIN
# Cut all nodes above the ELA
contour_elev_interp = interpolate.NearestNDInterpolator(sfc[:, :2], sfc[:, 2])
contour_elev = contour_elev_interp(contour)
contour_above_ELA = np.where(contour_elev>=ELA)[0]

# Label all boundary nodes with their index to fit nodes to cut
for i in range(len(contour[:, 0])):
    ax.text(contour[i, 0]/1e3, contour[i, 1]/1e3, labels[i])

# Record ranges of indices to cut
branches_to_cut = [
        np.arange(26, 57),
        np.arange(127, 196),
        np.arange(199, 225),
        np.arange(287, 307),
        contour_above_ELA,
        np.array([333, 334, 435, 68,]) # Chosen manually to smooth out corners, depends on chosen ELA
        ]

# Merge the indices to cut into a single array
indices_to_cut = np.unique(np.concatenate(branches_to_cut))
print('Cutting indices:', indices_to_cut)

# Cut out these indices from the contour array
contour_indices = np.arange(contour.shape[0], dtype=float)
contour_indices[indices_to_cut] = np.nan
contour_simplified = contour[~np.isnan(contour_indices), :]

ax.plot(contour_simplified[:, 0]/1e3, contour_simplified[:, 1]/1e3,
    'r.-')


fig.savefig('KW_simple_thickness.png', dpi=600)

# Save a nicely formatted text file with outline
np.savetxt('KW_outline_simplified.txt', contour_simplified, delimiter=',', fmt='%.3f')

header = """## KWsmooth.outline.exp
## Icon:0
# Points Count  Value
%d 1.000000
# X pos Y pos""" % contour_simplified.shape[0]

print('ISSM header:')
print(header)

np.savetxt('KWsmooth.outline.exp', contour_simplified, header=header, fmt='%.2f', comments='')

plt.show()
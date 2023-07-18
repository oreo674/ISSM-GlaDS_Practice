"""

Synthetic kaskawulsh geometry based on
https://shmip.bitbucket.io/instructions.html

"""

import numpy as np
from matplotlib import pyplot as plt

# CONSTANTS
para = 0.0
# x_max = 30e3
x_max = 100e3
y_scale = 3
z_scale = 3
# z_scale = 2
eps = 1e-16

# DEFINE FUNCTIONS
def surface(x, y):
    return z_scale*(100*(x*6e3/x_max + 200)**(1/4) + 1/60*(x*6e3/x_max) - 100*(200)**(1/4)) + 1e-6

def bed(x, y, para=para):
    return f(x, para) + g(y)*h(x, para)

def f(x, para=para):
    return (surface(x_max, 0) - para*6e3)/x_max**2 * x**2*(6e3/x_max)**0 + para*(x*6e3/x_max)

def g(y):
    return 0.5e-6 * np.abs(y/y_scale)**3

def h(x, para=para):
    num = (-4.5*x/x_max + 5) * (surface(x, 0) - f(x, para))
    denom = surface(x, 0) - f(x, 0.05) + eps
    return num/denom
    
def ginv(x):
    return y_scale*(x/0.5e-6)**(1/3)

def outline(x):
    return ginv( (surface(x, 0) - f(x, 0.05)) / (h(x, 0.05) + eps) )

print(outline(0))

# PLOT FLOWLINE
xx = np.linspace(0, x_max, 101)
fig, ax = plt.subplots()
ax.plot(xx, surface(xx, 0), color='b', label='Surface')
ax.plot(xx, bed(xx, 0), color='g', label='Bed')
ax.plot(xx, surface(xx, 0) - bed(xx, 0), color='k', label='Thickness')
ax.grid()
ax.set_xlabel('x (m)')
ax.set_ylabel('z (m)')
ax.axvline(30e3, color='k', linestyle=':')
ax.legend()

# PLOT CROSS-SECTION
fig2, ax2 = plt.subplots()
xpos = x_max/2
halfwidth = outline(xpos)
yy = np.linspace(-halfwidth, halfwidth, 101)
xx2 = xpos + 0*yy

ax2.plot(yy, surface(xx2, yy), color='b', label='Surface')
ax2.plot(yy, bed(xx2, yy), color='g', label='Bed')
ax2.grid()
ax2.set_xlabel('y (m)')
ax2.set_ylabel('z (m)')
ax2.legend()
ax2.set_title('x = %d m, full width = %.1f m' % (xpos, halfwidth*2))

# PLOT CONTOURS
fig3, ax3 = plt.subplots()
yy = np.linspace(-halfwidth*1.25, halfwidth*1.25, 251)
[xxs, yys] = np.meshgrid(xx, yy)
levels = np.arange(0, 600*z_scale, 200)
upper_outline = outline(xx)
surf = surface(xxs, yys)
b = bed(xxs, yys)
surf[np.abs(yys)>upper_outline] = np.nan
b[np.abs(yys)>upper_outline] = np.nan
ax3.contour(xxs, yys, surf, colors='b', label='Surface', levels=levels)
ax3.contour(xxs, yys, b, colors='g', label='Bed', levels=levels)
ax3.plot(xx, upper_outline, color='k', linewidth=2)
ax3.plot(xx, -upper_outline, color='k', linewidth=2)
ax3.set_xlabel('x (m)')
ax3.set_ylabel('y (m)')

ax3.set_title('Length = %d m, Width = %d m, Depth = %d m, $\Delta$z = %d m' % (
    x_max, 2*np.max(upper_outline), np.nanmax(surf-b), np.nanmax(surf) - np.nanmin(surf)))
fig3.savefig('synth_kask.png', dpi=600)

plt.show()




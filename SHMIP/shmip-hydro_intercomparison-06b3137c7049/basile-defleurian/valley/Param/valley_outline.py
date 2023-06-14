import numpy as np
from expwrite import *

def valley_outline():
	#  y = valley_outline(x)
	# Gives the upper half of the glacier margin.
	# Needs to be kept in sync with valley.m!

	# domain length
	xend = 6.0e3
	step=10.0
	x = np.arange(0,xend+step,step)

	# surf para
	beta = 0.25
	s1 = 100.0
	s2 = 100.0/xend
	sx0 = -200.
	s3 = 0.
	
	# bed para
	g1 = 0.5e-6
	alpha = 3.

	s_xend = s1*(xend-sx0)**beta + s2*xend - s1*(-sx0)**beta + s3
	surf = s1*(x-sx0)**beta + s2*x - s1*(-sx0)**beta + s3	
	# bed:

	f20 = 300./xend
	f10 = (s_xend - f20*xend)/xend**2
	f0 = f10*x**2 + f20*x
	h0 = -4.5*(x/xend) + 5

	y = (((surf-f0)/h0)/g1)**(1/alpha)

	YPos=numpy.hstack((y,-y[::-1]))
	XPos=numpy.hstack((x,x[::-1]))

	PosDict={'x':XPos,
					 'y':YPos}

	expwrite(PosDict,'../Exp/GlacierOutline.exp')

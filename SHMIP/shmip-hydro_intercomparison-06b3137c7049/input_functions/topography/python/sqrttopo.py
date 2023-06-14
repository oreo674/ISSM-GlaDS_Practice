import numpy as np
def sqrttopo(x):
	''' Function defining the ice thickness for the SHMIP exercise.
	This is a square root shape on a flat bed
    
	Input:
	x coordinates (vector with the x coordinate of all points)
	Output:
	bed elevation and ice thickness in that order
	'''
	
	thick_para=6.
	thick_off=5.0e3
	min_thick=1.0
	
	#Geometry
	thickness=min_thick+thick_para*pow(x+thick_off,0.5)-6.0*pow(thick_off,0.5)
	base=np.zeros(np.shape(x))

	return base,thickness

import numpy as np

def SeasonForcing(surf,time,Tparam=None):
	'''Runoff=SeasonForcing(surf,time,Tparam)
	Definition of the seasonal input function for experiment D and F of SHMIP.
	
	Input:
	surf : is the surface topography in meters (one vector)
	time : are the timestamps of the requested timesteps in seconds (one vector)
	Tparam: is the delta T parameter as described in SHMIP exercise

	Output:
	Output is the runoff under the form of a two D array (number of points,number of times)
	'''
	day = 24.*60.*60. #From second to day
	year = 365 *day   # From second to year

	DDF = 0.01/day ## degree day factor m/K/s
	lr = -0.0075 #lapse rate K/m
	
	BackInput = 7.93e-11 #value of the background input

	if Tparam is None:
		Dt = 0.
	else:
		Dt = Tparam
	
	Temp= -16.*np.cos(2.*np.pi*time/year)-5 +Dt

	TotInput=BackInput*np.ones((np.size(surf),np.size(Temp)))
	for t,value in enumerate(Temp):
		PosVal=np.where((surf*lr+value)>0)[0]
		TotInput[PosVal,t]+=(surf[PosVal]*lr+value)*DDF

	return TotInput

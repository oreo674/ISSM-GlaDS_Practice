import numpy as np

def DiurnalForcing(time,RelAmp=None):
	'''Runoff=SeasonForcing(surf,time,Tparam)
	Definition of the diurnal input function for experiment C of the SHMIP exercise
	This is only the Moulin input part, distributed source of run A1 has to be added
	
	Input:
	time : are the timestamps of the requested timesteps (one vector)
	RelAmp: is the relative amplitude of the signal as given in SHMIP 

	Output:
	OUtput is the transformed input which has to be used for every moulin.
	'''
	day = 24.*60.*60. #From second to day
	MI = 0.9 #initial noulin input

	if RelAmp is None:
		ra = 0.
	else:
		ra = RelAmp
	

	Input=MI*(1.-ra*np.sin(2.*np.pi*time/day))
	Input[np.where(Input<0)]=0.0

	return Input

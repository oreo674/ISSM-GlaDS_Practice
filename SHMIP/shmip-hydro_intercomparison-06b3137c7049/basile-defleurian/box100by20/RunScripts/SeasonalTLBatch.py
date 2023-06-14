import numpy as np
from loadmodel import *
from seasonal_runoff import SeasonForcing
from re import split
import organizer as org

# {{{ RunSe(Temp,dt,prefix,RunName):
def RunSe(Temp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	if Org.perform(RunName):

		LoadName='../Models/A5b'+RunName+'.nc'
		md=loadmodel(LoadName)
		md.miscellaneous.name=prefix+RunName

		#Initialization
		lasttime=np.size(md.results.TransientSolution)-1
		md.initialization.sediment_head	 = md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node = md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness	 = md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head			 = md.results.TransientSolution[lasttime].EplHead

		md.results= []
		#Forcing
		daily     = (24.0*3600.0)/31536000.
		times			=	np.arange(0,5+daily,daily)
		secondtime= times*31536000.
		md.basalforcings.groundedice_melting_rate=SeasonForcing(md.geometry.surface,secondtime,Temp)
		md.basalforcings.groundedice_melting_rate=md.basalforcings.groundedice_melting_rate*md.constants.yts
		md.basalforcings.groundedice_melting_rate=np.append(md.basalforcings.groundedice_melting_rate,[times],axis=0)

		#Time
		md.settings.output_frequency=1./(dt*365.)
		md.timestepping.time_step=dt
		md.timestepping.final_time=np.max(times)

		return md,Org
# }}}
# {{{ ReRunSe(Temp,dt,prefix,RunName):
def ReRunSe(Temp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	if Org.perform(RunName):

		loadprefix=''.join(split('(re)',prefix)[:-2])
		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName
		md.hydrology.max_iter=500
		md.hydrology.rel_tol =1.0e-5

		#deal with time and forcing
		lasttime=np.size(md.results.TransientSolution)-1
		while md.results.TransientSolution[lasttime].time%(3600*24)!=0:
			lasttime=lasttime-1#get yourself at the start of a day
		starttime=md.results.TransientSolution[lasttime].time/md.constants.yts
		duration=5.0 #one days
		endtime=starttime+duration
		stepping= (24.0*3600.0)/31536000.
		times=np.arange(starttime,endtime+stepping,stepping)
		secondtime=times*md.constants.yts

		md.basalforcings.groundedice_melting_rate=SeasonForcing(md.geometry.surface,secondtime,Temp)
		md.basalforcings.groundedice_melting_rate=md.basalforcings.groundedice_melting_rate*md.constants.yts
		md.basalforcings.groundedice_melting_rate=np.append(md.basalforcings.groundedice_melting_rate,[times],axis=0)

		#Initialization
		md.initialization.sediment_head	 = md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node = md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness	 = md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head			 = md.results.TransientSolution[lasttime].EplHead
		md.results= []

		#Time
		md.settings.output_frequency=1./(dt*365.)
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=endtime

		return md,Org
# }}}
# {{{ SeasonSub(Temp,dt,prefix,RunName):
def SeasonSub(Temp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	if Org.perform(RunName):

		loadprefix=''.join(split('(sub)',prefix)[:-2])
		md=loadmodel('../Models/'+loadprefix+'rere'+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName

		#deal with time and forcing
		lasttime=np.size(md.results.TransientSolution)-1
		while md.results.TransientSolution[lasttime].time%(3*md.constants.yts)!=0:
			lasttime=lasttime-1#get yourself at the start of a year
		starttime=0.0
		duration=1.0 #one year
		endtime=0.0+duration
		stepping= (24.0*3600.0)/31536000.
		times=np.arange(starttime,endtime+stepping,stepping)
		secondtime=times*md.constants.yts

		md.basalforcings.groundedice_melting_rate=SeasonForcing(md.geometry.surface,secondtime,Temp)
		md.basalforcings.groundedice_melting_rate=md.basalforcings.groundedice_melting_rate*md.constants.yts
		md.basalforcings.groundedice_melting_rate=np.append(md.basalforcings.groundedice_melting_rate,[times],axis=0)

		#Initialization
		md.initialization.sediment_head	 = md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node = md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness	 = md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head			 = md.results.TransientSolution[lasttime].EplHead
		md.results= []

		#Time
		md.settings.output_frequency=1./(dt*365.)
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=endtime

		return md,Org
# }}}

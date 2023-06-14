import numpy as np
from loadmodel import loadmodel
from seasonal_runoff import SeasonForcing
import organizer as org
from re import split
# {{{ def InitRun(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Temp,dt,prefix,RunName):
def InitRun(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Temp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):

		md=loadmodel('../Models/E1Valley.nc')

		md.settings.waitonlock=0
		md.verbose.convergence=True

		md.miscellaneous.name=prefix+RunName

		front=np.nonzero(md.mesh.x<=5.)[0]
		if np.size(front)<1:
			print 'WARNING, no snout given'

		#Hydro Model Parameters
		md.hydrology.isefficientlayer=1
		md.hydrology.sedimentlimit_flag=2
		md.hydrology.rel_tol=1.0e-6
		md.hydrology.penalty_lock=10
		md.hydrology.max_iter=100
		md.hydrology.transfer_flag=1
		md.hydrology.leakage_factor=Leakage
	
		#Forcing
		md.basalforcings.groundedice_melting_rate = 7.93e-11*md.constants.yts*np.ones((md.mesh.numberofvertices))	

		#Sediment Parameters
		md.hydrology.sediment_transmitivity= SedTrans*np.ones((md.mesh.numberofvertices))	
		md.hydrology.sediment_thickness=SedThick	

		#Epl Parameters
		md.hydrology.mask_eplactive_node=np.zeros((md.mesh.numberofvertices))
		md.hydrology.epl_conductivity=EplKond
		md.hydrology.epl_initial_thickness=EplThick
		md.hydrology.epl_colapse_thickness=EplCol
		md.hydrology.eplflip_lock=20

		#Boundaries
		md.hydrology.spcsediment_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcepl_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcsediment_head[front]=md.geometry.base[front]
		md.hydrology.spcepl_head[front]=md.geometry.base[front]

		#Initialisation
		md.initialization.sediment_head=md.geometry.base-SedThick
		md.initialization.epl_thickness=EplThick*np.ones((md.mesh.numberofvertices))
		md.initialization.epl_head=md.geometry.base
		
		#Time
		md.settings.output_frequency=1000
		md.timestepping.time_step=dt
		md.timestepping.final_time=10

		return md,Org
# }}}
# {{{ def SeasonRun(Temp,dt,prefix,RunName):
def SeasonRun(Temp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		md=loadmodel('../Models/Finit'+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName
		
		#Initialization
		lasttime=np.size(md.results.TransientSolution)-1
		md.initialization.sediment_head	 = md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node = md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness	 = md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head			 = md.results.TransientSolution[lasttime].EplHead

		md.results=[]

		#Forcing
		daily      = (24.0*3600.0)/31536000.
		times			 =	np.arange(0,1+daily,daily)
		secondtime = times*md.constants.yts
		md.basalforcings.groundedice_melting_rate = SeasonForcing(md.geometry.surface,secondtime,Temp)
		md.basalforcings.groundedice_melting_rate = md.basalforcings.groundedice_melting_rate*md.constants.yts
		md.basalforcings.groundedice_melting_rate = np.append(md.basalforcings.groundedice_melting_rate,[times],axis=0)
		#Time
		md.settings.output_frequency=1./(365.*dt)
		md.timestepping.time_step=dt
		md.timestepping.final_time=np.max(times)

		return md,Org
# }}}
# {{{ def SeasonReRun(Temp,dt,prefix,RunName):
def SeasonReRun(Temp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		loadprefix=''.join(split('(re)',prefix)[:-2])
		#loadprefix='F5re'
		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName
	
		#Initialization
		lasttime=np.size(md.results.TransientSolution)-1
		while md.results.TransientSolution[lasttime].time%(3600*24)!=0:
			lasttime=lasttime-1#get yourself at the start of a day
		starttime=md.results.TransientSolution[lasttime].time/md.constants.yts
		duration=5.0 #years
		endtime=starttime+duration
		stepping= (24.0*3600.0)/31536000.
		times=np.arange(starttime,endtime+stepping,stepping)
		secondtime=times*md.constants.yts
		print('run from {} to {}'.format(starttime,endtime))

		#Forcing
		md.basalforcings.groundedice_melting_rate =	SeasonForcing(md.geometry.surface,secondtime,Temp)
		md.basalforcings.groundedice_melting_rate=md.basalforcings.groundedice_melting_rate*md.constants.yts
		md.basalforcings.groundedice_melting_rate=np.append(md.basalforcings.groundedice_melting_rate,[times],axis=0)

		md.initialization.sediment_head	 = md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node = md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness	 = md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head			 = md.results.TransientSolution[lasttime].EplHead

		md.results=[]
		#Time
		md.settings.output_frequency=1./(365.*dt)
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=endtime

		return md,Org
# }}}
# {{{ def SeasonSub(Temp,dt,prefix,RunName):
def SeasonSub(Temp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		if prefix=='F1bEnd':
			loadprefix='F1bSub'
			dt			 = 3600.0/(31536000.)
		elif prefix=='F2bEnd':
			loadprefix='F2bSub'
			dt			 = 3600.0/(31536000.)
		elif prefix=='F3bEnd':
			loadprefix='F3bSub'
			dt			 = 3600.0/(6.*31536000.)
		elif prefix=='F4bEnd':
			loadprefix='F4bSub'
			dt			 = 3600.0/(10.*31536000.)
		elif prefix=='F5bEnd':
			loadprefix='F5bSub'
			dt			 = 3600.0/(10.*31536000.)

		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName
	
		#Initialization
		lasttime=np.size(md.results.TransientSolution)-1
		while md.results.TransientSolution[lasttime].time%(md.constants.yts)!=0:
			lasttime=lasttime-1#get yourself at the start of a year
		starttime=0.
		duration=1.0 #years
		endtime=starttime+duration
		stepping= (24.0*3600.0)/31536000.
		times=np.arange(starttime,endtime+stepping,stepping)
		secondtime=times*md.constants.yts
		print('run from {} to {}'.format(starttime,endtime))

		#Forcing
		md.basalforcings.groundedice_melting_rate =	SeasonForcing(md.geometry.surface,secondtime,Temp)
		md.basalforcings.groundedice_melting_rate=md.basalforcings.groundedice_melting_rate*md.constants.yts
		md.basalforcings.groundedice_melting_rate=np.append(md.basalforcings.groundedice_melting_rate,[times],axis=0)

		md.initialization.sediment_head	 = md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node = md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness	 = md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head			 = md.results.TransientSolution[lasttime].EplHead

		md.results=[]
		#Time
		md.settings.output_frequency=1./(365.*dt)
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=endtime

		return md,Org
# }}}

import numpy as np
from loadmodel import loadmodel
from diurnal_runoff import DiurnalForcing
from re import split
import organizer as org

# {{{ RunDi(RelAmp,dt,prefix,RunName):
def RunDi(RelAmp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):

		md=loadmodel('../Models/B5Sub'+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName

		#Forcing
		# stepping=300./md.constants.yts # 5 min forcing
		stepping=3600./md.constants.yts #hourly forcing
		endtime=1.0 #one year
		times=np.arange(0.,endtime+stepping,stepping)
		secondtime=times*md.constants.yts

		inputval=DiurnalForcing(secondtime,RelAmp)
		md.hydrology.basal_moulin_input=np.outer(md.hydrology.basal_moulin_input,inputval)
		md.hydrology.basal_moulin_input=np.append(md.hydrology.basal_moulin_input,[times],axis=0)
		md.basalforcings.groundedice_melting_rate = 7.93e-11*md.constants.yts*np.ones((md.mesh.numberofvertices))
	
		#initialize from results
		lasttime=np.size(md.results.TransientSolution)-1
		md.initialization.sediment_head=md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node=md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness=md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head=md.results.TransientSolution[lasttime].EplHead
		md.results=[]

		#Time
		md.settings.output_frequency=144
		md.timestepping.time_step=dt
		md.timestepping.start_time=0.0
		md.timestepping.final_time=endtime

		return md,Org
# }}}
# {{{ ReRunDi(RelAmp,dt,prefix,RunName):
def ReRunDi(RelAmp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):

		loadprefix=''.join(split('(re)',prefix)[:-2])
		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName

		duration=1.0 #one year

		#deal with time and forcing
		lasttime=np.size(md.results.TransientSolution)-1
		while md.results.TransientSolution[lasttime].time%(3600*24)!=0:
			lasttime=lasttime-1#get yourself at the start of a day
		starttime=md.results.TransientSolution[lasttime].time/md.constants.yts
		endtime=starttime+duration
		stepping=3600./md.constants.yts # 5 min forcing
		times=np.arange(starttime,endtime+stepping,stepping)
		secondtime=times*md.constants.yts

		inputval=DiurnalForcing(secondtime,RelAmp)
		MoulinInput=np.zeros((md.mesh.numberofvertices))
		MoulinInput[np.where(md.hydrology.basal_moulin_input[:-1,0]>0)]=0.9
		md.hydrology.basal_moulin_input=np.outer(MoulinInput,inputval)
		md.hydrology.basal_moulin_input=np.append(md.hydrology.basal_moulin_input,[times],axis=0)
		md.basalforcings.groundedice_melting_rate = 7.93e-11*md.constants.yts*np.ones((md.mesh.numberofvertices))

		#initialize from results
		md.initialization.sediment_head=md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node=md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness=md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head=md.results.TransientSolution[lasttime].EplHead
		md.results=[]

		#Time
		md.settings.output_frequency=144
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=endtime

		return md,Org
# }}}
# {{{ SubDi(RelAmp,dt,prefix,RunName):
def SubDi(RelAmp,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):

		loadprefix=prefix[0:2]+'rerere'
		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName

		duration=1.0/365. #one days

		#deal with time and forcing
		lasttime=np.size(md.results.TransientSolution)-1
		while md.results.TransientSolution[lasttime].time%(3600*24)!=0:
			lasttime=lasttime-1#get yourself at the start of a day
		starttime=0
		endtime=starttime+duration
		stepping=300./md.constants.yts # 5 min forcing
		times=np.arange(starttime,endtime+stepping,stepping)
		secondtime=times*md.constants.yts

		inputval=DiurnalForcing(secondtime,RelAmp)
		MoulinInput=np.zeros((md.mesh.numberofvertices))
		MoulinInput[np.where(md.hydrology.basal_moulin_input[:-1,0]>0)]=0.9
		md.hydrology.basal_moulin_input=np.outer(MoulinInput,inputval)
		md.hydrology.basal_moulin_input=np.append(md.hydrology.basal_moulin_input,[times],axis=0)
		md.basalforcings.groundedice_melting_rate = 7.93e-11*md.constants.yts*np.ones((md.mesh.numberofvertices))

		#initialize from results
		md.initialization.sediment_head=md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node=md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness=md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head=md.results.TransientSolution[lasttime].EplHead
		md.results=[]

		#Time
		md.settings.output_frequency=12
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=endtime

		return md,Org
# }}}

import numpy as np
from model import *
from triangle import *
from setmask import *
from parameterize import *
from TsOptimize import TsFit
from loadmodel import loadmodel
from re import split
import organizer as org

# {{{ RunSt(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName):
def RunSt(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		
		md=triangle(model(),'../Exp/Box100by20.exp',1000.)
		md=setmask(md,'','')
		md=parameterize(md,'../Param/BoxSetUp.py')

		md.settings.waitonlock=0
		md.verbose.convergence=True
	
		md.miscellaneous.name=prefix+RunName

		#Hydro Model Parameters
		md.hydrology.isefficientlayer=1
		md.hydrology.sedimentlimit_flag=2
		md.hydrology.rel_tol=1.0e-5
		md.hydrology.penalty_lock=10
		md.hydrology.max_iter=500
		md.hydrology.transfer_flag=1
		md.hydrology.leakage_factor=Leakage
	
		#Forcing
		md.basalforcings.groundedice_melting_rate = Input*np.ones((md.mesh.numberofvertices))
	
		#Sediment Parameters
		md.hydrology.sediment_transmitivity= SedTrans*np.ones((md.mesh.numberofvertices))	
		md.hydrology.sediment_thickness=SedThick
	
		#Epl Parameters
		md.hydrology.mask_eplactive_node=np.zeros((md.mesh.numberofvertices))
		md.hydrology.epl_conductivity=EplKond
		md.hydrology.epl_initial_thickness=EplThick
		md.hydrology.epl_colapse_thickness=EplCol
		md.hydrology.eplflip_lock=0

		#Initialisation
		md.initialization.sediment_head=md.geometry.base#-SedThick
		md.initialization.epl_head=md.geometry.base
		md.initialization.epl_thickness=EplThick*np.ones((md.mesh.numberofvertices))
		
		#Boundary conditions
		front=np.where(md.mesh.x==0.)[0]
		md.hydrology.spcsediment_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcsediment_head[front]=0.0
		md.hydrology.spcepl_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcepl_head[front]=0.0

		#Time
		md.settings.output_frequency=1/(10*dt)
		md.timestepping.time_step=dt
		md.timestepping.final_time=60

		return md,Org
# }}}
# {{{ ReRunSt(Input,dt,prefix,RunName):
def ReRunSt(Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		
		loadprefix=''.join(split('(re)',prefix)[:-2])
		loadprefix='A5b'
		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName
		last=np.size(md.results.TransientSolution)-1

		#Forcing
		md.basalforcings.groundedice_melting_rate = Input*np.ones((md.mesh.numberofvertices))

		#Epl Parameters
		md.hydrology.mask_eplactive_node=md.results.TransientSolution[last].HydrologydcMaskEplactiveNode

		#Initialisation
		md.initialization.sediment_head=md.results.TransientSolution[last].SedimentHead
		md.initialization.epl_thickness=md.results.TransientSolution[last].HydrologydcEplThickness
		md.initialization.epl_head=md.results.TransientSolution[last].EplHead
		starttime=md.results.TransientSolution[last].time/md.constants.yts
		md.results=[]
	
		#Time
		md.settings.output_frequency=60/(1000*dt)
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=starttime+60.0

		return md,Org
# }}}
# {{{ SubSt(Input,dt,prefix,RunName):
def SubSt(Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		
		loadprefix=prefix[0:2]+'b'
		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName
		last=np.size(md.results.TransientSolution)-1

		#Forcing
		md.basalforcings.groundedice_melting_rate = Input*np.ones((md.mesh.numberofvertices))

		#Epl Parameters
		md.hydrology.mask_eplactive_node=md.results.TransientSolution[last].HydrologydcMaskEplactiveNode

		#Initialisation
		md.initialization.sediment_head=md.results.TransientSolution[last].SedimentHead
		md.initialization.epl_thickness=md.results.TransientSolution[last].HydrologydcEplThickness
		md.initialization.epl_head=md.results.TransientSolution[last].EplHead
		starttime=md.results.TransientSolution[last].time/md.constants.yts
		md.results=[]
	
		#Time
		md.settings.output_frequency=1
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=starttime+5*dt

		return md,Org
# }}}
# {{{ PolyRunSt(SedPolyDeg,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName):
def PolyRunSt(SedPolyDeg,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		
		md=triangle(model(),'../Exp/Box100by20.exp',1000.)
		md=setmask(md,'','')
		md=parameterize(md,'../Param/BoxSetUp.py')

		md.settings.waitonlock=0
		md.verbose.convergence=True
		md.miscellaneous.name=prefix+RunName

		#Hydro Model Parameters
		md.hydrology.isefficientlayer=1
		md.hydrology.sedimentlimit_flag=2
		md.hydrology.rel_tol=1.0e-6
		md.hydrology.penalty_lock=10
		md.hydrology.max_iter=200
		md.hydrology.transfer_flag=1
		md.hydrology.leakage_factor=Leakage
	
		#Forcing
		md.basalforcings.groundedice_melting_rate = Input*np.ones((md.mesh.numberofvertices))
	
		#Sediment Parameters
		PolySed=TsFit(SedPolyDeg)
		md.hydrology.sediment_transmitivity=PolySed(md.mesh.x)
		md.hydrology.sediment_transmitivity[np.where(md.mesh.x<13920.)]=1.025-np.tanh((md.mesh.x[np.where(md.mesh.x<13920.)]+2000.)/4000.)
		md.hydrology.sediment_transmitivity[np.where(md.hydrology.sediment_transmitivity<1.0e-6)]=1.0e-6
		md.hydrology.sediment_thickness=SedThick
	
		#Boundaries
		front=np.nonzero(md.mesh.x==0.)[0]
		md.hydrology.spcsediment_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcsediment_head[front]=0.0
		md.hydrology.spcepl_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcepl_head[front]=0.0

		#Epl Parameters
		md.hydrology.mask_eplactive_node=np.zeros((md.mesh.numberofvertices))
		md.hydrology.epl_conductivity=EplKond
		md.hydrology.epl_colapse_thickness=EplCol
		md.hydrology.eplflip_lock=10
		md.hydrology.epl_initial_thickness=EplThick

		#Initialisation
		md.initialization.sediment_head=md.geometry.base#-SedThick
		md.initialization.epl_head=md.geometry.base
		md.initialization.epl_thickness=EplThick*np.ones((md.mesh.numberofvertices))
	
		#Time
		md.settings.output_frequency=0.1/dt
		md.timestepping.time_step=dt
		md.timestepping.final_time=60.0

		return md,Org
# }}}
# {{{ PolyReRunSt(Input,dt,prefix,RunName):
def PolyReRunSt(Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		
		md=loadmodel('../Models/A5-P2-es10-L1e-07-Ke110.0-ee0.0005-Col5e-07.nc')
		md.miscellaneous.name=prefix+RunName
		last=np.size(md.results.TransientSolution)-1

		#Forcing
		md.basalforcings.groundedice_melting_rate = Input*np.ones((md.mesh.numberofvertices))
	
		#Epl Parameters
		md.hydrology.mask_eplactive_node=md.results.TransientSolution[last].HydrologydcMaskEplactiveNode

		#Initialisation
		md.initialization.sediment_head=md.results.TransientSolution[last].SedimentHead
		md.initialization.epl_thickness=md.results.TransientSolution[last].HydrologydcEplThickness
		md.initialization.epl_head=md.results.TransientSolution[last].EplHead
		md.results=[]

		#Time
		md.settings.output_frequency=0.1/dt
		md.timestepping.time_step=dt
		md.timestepping.final_time=60.0

		return md,Org
# }}}

import numpy as np
from model import *
from loadmodel import *
from setmask import *
from parameterize import *
from bamg import bamg
import organizer as org
from re import split

# {{{ def Mesher(prefix,Slope)

def Mesher(prefix,Slope):
	from valley import valley
	steps=[1]
	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)

	MeshName=prefix+'Valley'
	if Org.perform('Init'):
		md=bamg(model(),'domain','../Exp/GlacierOutline.exp','hmax',200.,'hmin',100.0)
		md=setmask(md,'','')
		md=parameterize(md,'../Param/ValleySetUp.py')

		#Geometry
		md.geometry.base,md.geometry.thickness=valley(md.mesh.x,md.mesh.y,Slope)
		md.geometry.surface=md.geometry.base+md.geometry.thickness
		md.miscellaneous.name='ValleyGeom_'+prefix

		return md,Org

# }}}
# {{{ def Mesher(prefix,Slope)

def Mesher3D(prefix,Slope):
	from valley import valley
	steps=[1]
	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)

	MeshName=prefix+'Valley'
	if Org.perform('Init'):
		md=bamg(model(),'domain','../Exp/GlacierOutline.exp','hmax',20.,'hmin',10.0)
		md=setmask(md,'','')
		md=parameterize(md,'../Param/ValleySetUp.py')

		#Geometry
		md.geometry.base,md.geometry.thickness=valley(md.mesh.x,md.mesh.y,Slope)
		md.geometry.surface=md.geometry.base+md.geometry.thickness
		md.miscellaneous.name='Valley3DGeom_'+prefix
		md.extrude(3,1.)

		return md,Org

# }}}
# {{{ def ValleySlope(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName):
def ValleySlope(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		
		md=loadmodel('../Models/'+prefix[0:2]+'Valley.nc')
		
		md.settings.waitonlock=0
		md.verbose.convergence=True
	
		md.miscellaneous.name=prefix+RunName

		front=np.nonzero(md.mesh.x<=5.)[0]
		if np.size(front)<1:
			print 'WARNING, no snout given'

		rho_ratio=md.materials.rho_ice/md.materials.rho_freshwater

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
		md.hydrology.sediment_transmitivity= SedTrans*np.ones((md.mesh.numberofvertices))	
		md.hydrology.sediment_thickness=SedThick

		#Epl Parameters
		md.hydrology.mask_eplactive_node=np.zeros((md.mesh.numberofvertices))
		md.hydrology.epl_conductivity=EplKond
		md.hydrology.epl_initial_thickness=EplThick
		md.hydrology.epl_colapse_thickness=EplCol
		md.hydrology.eplflip_lock=20
	

		#Initialisation
		md.initialization.sediment_head=0.0*(md.mesh.x*(6000-md.mesh.x)*1.158e-6/0.0319)
#		md.initialization.sediment_head[np.where(md.initialization.sediment_head<md.geometry.base)]=md.geometry.base[np.where(md.initialization.sediment_head<md.geometry.base)]
#md.geometry.base-SedThick
		md.initialization.epl_thickness=EplThick*np.ones((md.mesh.numberofvertices))
		md.initialization.epl_head=md.geometry.base

		#Boundaries
		md.hydrology.spcsediment_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcepl_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcsediment_head[front]=md.geometry.base[front]
		md.hydrology.spcepl_head[front]=md.geometry.base[front]
	
		#Time
		md.settings.output_frequency=1/(1000*dt)
		md.timestepping.time_step=dt
		md.timestepping.final_time=1

		return md,Org
# }}}
# {{{ def ValleyReSlope(Input,dt,prefix,RunName):
def ValleyReSlope(Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):

		loadprefix=''.join(split('(re)',prefix)[:-2])
		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName
		
		md.settings.waitonlock=0
		md.verbose.convergence=True
		rho_ratio=md.materials.rho_ice/md.materials.rho_freshwater

		#Forcing
		md.basalforcings.groundedice_melting_rate = Input*np.ones((md.mesh.numberofvertices))
	
		#Initialisation
		lasttime=np.size(md.results.TransientSolution)-1

		starttime=md.results.TransientSolution[lasttime].time/md.constants.yts
		duration=5.0 #years
		endtime=starttime+duration

		md.initialization.sediment_head=md.results.TransientSolution[lasttime].SedimentHead
		md.hydrology.mask_eplactive_node = md.results.TransientSolution[lasttime].HydrologydcMaskEplactiveNode
		md.initialization.epl_thickness=md.results.TransientSolution[lasttime].HydrologydcEplThickness
		md.initialization.epl_head=md.results.TransientSolution[lasttime].EplHead
		md.results=[]

		#Time
		md.settings.output_frequency=1.0/(100*dt)
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=endtime

		return md,Org
# }}}

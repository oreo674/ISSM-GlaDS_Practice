import numpy as np
import csv
from model import *
from bamg import *
from setmask import *
from parameterize import *
from TsOptimize import TsFit
from scipy import spatial
from re import split
from loadmodel import loadmodel
import organizer as org

# {{{ RunMo(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName)
def RunMo(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		Xs=[]
		Ys=[]
		with open('../Exp/AllMoulins.csv') as csvfile:
			MoulinPos=csv.reader(csvfile,delimiter=',')
			for row in MoulinPos:
				Xs+=[float(row[0])]
				Ys+=[float(row[1])]

		XY=np.zeros((np.size(Xs),2))
		XY[:,0]=Xs
		XY[:,1]=Ys

		md=bamg(model(),'domain','../Exp/Box100by20.exp','hmax',1000.,'RequiredVertices',XY)
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
		md.hydrology.max_iter=100
		md.hydrology.transfer_flag=1
		md.hydrology.leakage_factor=Leakage
	
		#Forcing
		moulinfile='../../../input_functions/source/'+Input
		Xs=[]
		Ys=[]
		Vals=[]
		with open(moulinfile) as csvfile:
			MoulinPos=csv.reader(csvfile,delimiter=',')
			for row in MoulinPos:
				Xs	 +=	[float(row[1])]
				Ys	 +=	[float(row[2])]
				Vals +=	[float(row[3])]

		XY				=	np.zeros((np.size(Xs),2))
		Values		=	np.zeros((np.size(Xs)))
		XY[:,0]		=	Xs[:]
		XY[:,1]		=	Ys[:]
		Values[:]	=	Vals[:]
		gridpoints=(np.vstack((md.mesh.x,md.mesh.y)).T)

		for i,inputval in enumerate(Values):
			distance,index = spatial.KDTree(gridpoints).query(XY[i])
			md.hydrology.basal_moulin_input[index]=inputval
			print 'Moulin is %e meters away from grid point' % distance

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

		#Initialisation
		md.initialization.sediment_head=md.geometry.base-SedThick
		md.initialization.epl_head=md.geometry.base
		md.initialization.epl_thickness=EplThick*np.ones((md.mesh.numberofvertices))

		#Boundary Conditions
		front=np.where(md.mesh.x==0.)[0]
		md.hydrology.spcsediment_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcsediment_head[front]=0.0
		md.hydrology.spcepl_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcepl_head[front]=0.0
	
		#Time
		md.settings.output_frequency=1000
		md.timestepping.time_step=dt
		md.timestepping.final_time=60.0

		return md,Org
# }}}
# {{{ ReRunMo(Input,dt,prefix,RunName)
def ReRunMo(Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):

		loadprefix=''.join(split('(re)',prefix)[:-2])
		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName
		last=np.size(md.results.TransientSolution)-1	
	
		#Epl Parameters
		md.hydrology.mask_eplactive_node=md.results.TransientSolution[last].HydrologydcMaskEplactiveNode

		#Initialisation
		md.initialization.sediment_head=md.results.TransientSolution[last].SedimentHead
		md.initialization.epl_thickness=md.results.TransientSolution[last].HydrologydcEplThickness
		md.initialization.epl_head=md.results.TransientSolution[last].EplHead
		starttime=md.results.TransientSolution[last].time/md.constants.yts
		md.results=[]

		#Time
		md.settings.output_frequency=1/(10*dt)
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=starttime+10.0

		return md,Org
# }}}
# {{{ SubMo(Input,dt,prefix,RunName)
def SubMo(Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):

		if prefix=='B1Sub':
			loadprefix='B1rerererere'
		elif prefix=='B2Sub':
			loadprefix='B2re'
		elif prefix=='B3Sub':
			loadprefix='B3rererere'
		elif prefix=='B4Sub':
			loadprefix='B4rerererere'
		elif prefix=='B5Sub':
			loadprefix='B5rerererere'

		md=loadmodel('../Models/'+loadprefix+RunName+'.nc')
		md.miscellaneous.name=prefix+RunName
		last=np.size(md.results.TransientSolution)-1	
	
		#Epl Parameters
		md.hydrology.mask_eplactive_node=md.results.TransientSolution[last].HydrologydcMaskEplactiveNode

		#Initialisation
		md.initialization.sediment_head=md.results.TransientSolution[last].SedimentHead
		md.initialization.epl_thickness=md.results.TransientSolution[last].HydrologydcEplThickness
		md.initialization.epl_head=md.results.TransientSolution[last].EplHead
		starttime=md.results.TransientSolution[last].time/md.constants.yts
		md.results=[]

		#Time
		md.settings.output_frequency=1/(10*dt)
		md.timestepping.time_step=dt
		md.timestepping.start_time=starttime
		md.timestepping.final_time=starttime+10.0

		return md,Org
# }}}
# {{{ PolyRunMo(SedPolyDeg,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName)
def PolyRunMo(SedPolyDeg,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName):
	steps=[1]

	Org = org.organizer('repository','../Models/','prefix',prefix,'steps',steps)
	
	if Org.perform(RunName):
		Xs=[]
		Ys=[]
		with open('../Exp/AllMoulins.csv') as csvfile:
			MoulinPos=csv.reader(csvfile,delimiter=',')
			for row in MoulinPos:
				Xs+=[float(row[0])]
				Ys+=[float(row[1])]

		XY=np.zeros((np.size(Xs),2))
		XY[:,0]=Xs
		XY[:,1]=Ys

		md=bamg(model(),'domain','../Exp/Box100by20.exp','hmax',1000.,'RequiredVertices',XY)
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
		md.hydrology.max_iter=100
		md.hydrology.transfer_flag=1
		md.hydrology.leakage_factor=Leakage
	
		#Forcing
		moulinfile='../../../input_functions/source/'+Input
		Xs=[]
		Ys=[]
		Vals=[]
		with open(moulinfile) as csvfile:
			MoulinPos=csv.reader(csvfile,delimiter=',')
			for row in MoulinPos:
				Xs	 +=	[float(row[1])]
				Ys	 +=	[float(row[2])]
				Vals +=	[float(row[3])]

		XY				=	np.zeros((np.size(Xs),2))
		Values		=	np.zeros((np.size(Xs)))
		XY[:,0]		=	Xs[:]
		XY[:,1]		=	Ys[:]
		Values[:]	=	Vals[:]
		gridpoints=(np.vstack((md.mesh.x,md.mesh.y)).T)

		for i,inputval in enumerate(Values):
			distance,index = spatial.KDTree(gridpoints).query(XY[i])
			md.hydrology.basal_moulin_input[index]=inputval
			print 'Moulin is %e meters away from grid point' % distance

		md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
	
		#Sediment Parameters
		PolySed=TsFit(SedPolyDeg)
		PolySed[0]=0.95*PolySed[0]
		md.hydrology.sediment_transmitivity=PolySed(md.mesh.x)
		md.hydrology.sediment_transmitivity[np.where(md.mesh.x<13500.)]=1.025-np.tanh((md.mesh.x[np.where(md.mesh.x<13500.)]+2000.)/4000.)
		md.hydrology.sediment_transmitivity[np.where(md.hydrology.sediment_transmitivity<1.0e-6)]=1.0e-6
		md.hydrology.sediment_thickness=SedThick
	
		#Epl Parameters
		md.hydrology.mask_eplactive_node=np.zeros((md.mesh.numberofvertices))
		md.hydrology.epl_conductivity=EplKond
		md.hydrology.epl_initial_thickness=EplThick
		md.hydrology.epl_colapse_thickness=EplCol
		md.hydrology.eplflip_lock=20

		#Initialisation
		md.initialization.sediment_head=md.geometry.base-SedThick
		md.initialization.epl_head=md.geometry.base
		md.initialization.epl_thickness=EplThick*np.ones((md.mesh.numberofvertices))

		#Boundary Conditions
		front=np.where(md.mesh.x==0.)[0]
		md.hydrology.spcsediment_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcsediment_head[front]=0.0
		md.hydrology.spcepl_head=np.nan*np.ones((md.mesh.numberofvertices))
		md.hydrology.spcepl_head[front]=0.0

		#Time
		md.settings.output_frequency=100
		md.timestepping.time_step=dt
		md.timestepping.final_time=60.0

		return md,Org
# }}}

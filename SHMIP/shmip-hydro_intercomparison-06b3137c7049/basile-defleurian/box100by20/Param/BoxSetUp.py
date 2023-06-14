from model import *
import numpy as np
from transient import *
from setflowequation import *
from SetIceSheetBC import SetIceSheetBC

#Start defining model parameters here

#Only dealing with Hydrology
md.transient=transient.setallnullparameters(md.transient)
md.transient.ishydrology=True
md=setflowequation(md,'SSA','all')

md.hydrology=hydrologydc()
md.hydrology=md.hydrology.initialize(md)

#Faking fields asked by transient
md.friction.coefficient				=	np.zeros((md.mesh.numberofvertices))
md.friction.q									=	np.zeros((md.mesh.numberofelements))
md.friction.p									=	np.zeros((md.mesh.numberofelements))
md.initialization.vx					=	np.zeros((md.mesh.numberofvertices))
md.initialization.vy					=	np.zeros((md.mesh.numberofvertices))
md.initialization.vx					=	np.zeros((md.mesh.numberofvertices))
md.initialization.vy					=	np.zeros((md.mesh.numberofvertices))
md.initialization.temperature	=	np.zeros((md.mesh.numberofvertices))
md.initialization.pressure		=	np.zeros((md.mesh.numberofvertices))

#Geometry
md.geometry.thickness=1.0+6.0*(pow(md.mesh.x+5000.0,0.5)-pow(5000.0,0.5))
md.geometry.base=np.zeros((md.mesh.numberofvertices))
md.geometry.surface=md.geometry.base+md.geometry.thickness

#Materials
md.constants.g=9.8
md.materials.rho_ice=910.0
md.materials.rho_freshwater=1000.0
md.materials.latentheat=334000.0
md.materials.rheology_law='None'
md.materials.rheology_n=3.0*np.ones((md.mesh.numberofelements))
md.materials.rheology_B=pow((2./(3**3*2.5e-25)),1./3.)*np.ones((md.mesh.numberofvertices))

#Necessarry fake for the ice model
md=SetIceSheetBC(md)

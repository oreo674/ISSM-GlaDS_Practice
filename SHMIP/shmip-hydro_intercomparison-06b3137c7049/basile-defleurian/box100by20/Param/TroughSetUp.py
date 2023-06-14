from model import *
from transient import *
from SetIceSheetBC import SetIceSheetBC


#Start defining model parameters here

#Geometry
md.geometry.thickness=1.0+(6.0*pow(md.mesh.x+5000.0,0.5)-6.0*pow(5000.0,0.5))
md.geometry.base=numpy.zeros((md.mesh.numberofvertices))
md.geometry.surface=md.geometry.base+md.geometry.thickness

md.initialization.vx=numpy.zeros((md.mesh.numberofvertices,1))
md.initialization.vy=numpy.zeros((md.mesh.numberofvertices,1))
md.initialization.vz=numpy.zeros((md.mesh.numberofvertices,1))
md.initialization.pressure=numpy.zeros((md.mesh.numberofvertices,1))

#Materials
md.materials.rho_ice=910.0
md.materials.rho_freshwater=1000.0
md.initialization.temperature=(273.-20.)*numpy.ones((md.mesh.numberofvertices,1))
md.materials.rheology_B=paterson(md.initialization.temperature)
md.materials.rheology_n=3.*numpy.ones((md.mesh.numberofelements,1))
md.initialization.temperature=md.initialization.temperature

#Friction
md.friction.coefficient=50.*numpy.ones((md.mesh.numberofvertices,1))
md.friction.p=numpy.ones((md.mesh.numberofelements,1))
md.friction.q=numpy.ones((md.mesh.numberofelements,1))

#Numerical parameters
md.stressbalance.viscosity_overshoot=0.3
md.masstransport.stabilization=1.
md.verbose=verbose(0)
md.settings.waitonlock=30
md.timestepping.time_step=1.
md.timestepping.final_time=2.
md.stressbalance.restol=0.05
md.stressbalance.reltol=1.
md.steadystate.reltol=1.
md.stressbalance.abstol=float('nan')



#Necessarry fake for the ice model
md=SetIceSheetBC(md)

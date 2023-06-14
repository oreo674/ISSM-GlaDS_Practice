from netCDF4 import Dataset
import numpy as np
import time
from os import path, remove
from slope import slope
# {{{ Two Layers
def netCDFTwoLay(md,filename): #specific for steady state with no time
	if path.exists(filename):
		print ('File {} allready exist'.format(filename))
		newname=raw_input('Give a new name or "delete" to replace: ')
		if newname=='delete':
			remove(filename)
		else:
			print ('New file name is {}'.format(newname))
			filename=newname

	NCData=Dataset(filename, 'w', format='NETCDF4')
	NCData.title = md.miscellaneous.name
	NCData.meshtype='unstructured'
	NCData.dimension = "2D" ;
	NCData.channels_on_edges = "no" ;
	NCData.institution = "Basile de Fleurian," ;
	NCData.source = "ISSM:svn-dev version on revision 21560" ;
	NCData.references = "http://shmip.bitbucket.io/"
	try:
		last= np.size(md.results.TransientSolution)
	except AttributeError:
		print 'No resuts presents, Exiting'
		NCData.close()
		return
	NumOfFields=len(md.results.TransientSolution[0].__dict__) #getting number of fields from first time

	if len(md.results.TransientSolution[last-1].__dict__)!=NumOfFields and len(md.results.TransientSolution[last-1].__dict__)!=NumOfFields-2:
		print ('Last time is truncated for Run {}, skipping'.format(filename))
		last=last-1

	# 2. Define dimension
	nodenr       = NCData.createDimension('index1',md.mesh.numberofvertices) # number of nodes
	eltnr       = NCData.createDimension('index2',md.mesh.numberofelements) # number of elements
	dim          = NCData.createDimension('dim',md.mesh.dimension()) # number of dimensions
	time         = NCData.createDimension('time',None) # timesteps,

	# 3. define variables
	time       = NCData.createVariable('time','f8',('time'))
	nodes      = NCData.createVariable('coords1','f8',('dim', 'index1'))
	elements   = NCData.createVariable('coords2','f8',('dim', 'index2'))
	H          = NCData.createVariable('H','f8',('index1'))
	B          = NCData.createVariable('B','f8',('index1'))
	N          = NCData.createVariable('N','f8',('time','index1'))
	Ee         = NCData.createVariable('Ee','f8',('time','index1'))
	q          = NCData.createVariable('q','f8',('time','index2'))
	Q          = NCData.createVariable('Q','f8',('time','index1'))

	time.long_name = "time"
	time.units = "s"
	nodes.long_name = "node coordinates"
	nodes.units = "m"
	elements.long_name = "element midpoint"
	elements.units = "m"
	H.long_name = "ice thickness"
	H.units = "m"
	B.long_name = "bed elevation"
	B.units = "m"
	N.long_name = "effective pressure"
	N.units = "Pa"
	Ee.long_name = "EPL thickness"
	Ee.units = "m"
	q.long_name = "water sheet discharge"
	q.units = "m^2/s"
	Q.long_name = "channel discharge"
	Q.units = "m^2/s"


	# 4. Place data
	EltmidX=np.nanmean(md.mesh.x[md.mesh.elements-1],axis=1)
	EltmidY=np.nanmean(md.mesh.y[md.mesh.elements-1],axis=1)

	nodes[:,:]   = np.column_stack((md.mesh.x,md.mesh.y)).T
	elements[:,:]= np.column_stack((EltmidX,EltmidY)).T
	H[:]         = md.geometry.thickness
	B[:]         = md.geometry.base
	for t in range(0,last,1):
		Tres=md.results.TransientSolution[t]
		time[t] = Tres.time
		N[t,:]  = Tres.EffectivePressure.reshape(md.mesh.x.shape)
		Ee[t,:]	= Tres.HydrologydcEplThickness.reshape(md.mesh.x.shape)

		#need some computation for the fluxes
		sedgradx,sedgrady,sedgrad=slope(md,Tres.SedimentHead)

		q[t,:]=sedgradx*np.nanmean(md.hydrology.sediment_transmitivity)

		eltthick=np.mean(Tres.HydrologydcEplThickness[md.mesh.elements-1],axis=1)
		fromelement=np.zeros(np.shape(md.mesh.x))
		fromelement[md.mesh.elements[np.where(Tres.HydrologydcMaskEplactiveElt[:,0]==1),:]-1]=1

		fromelementthick=np.zeros(np.shape(md.mesh.x))
		Connect=np.zeros(np.shape(md.mesh.x))

		for i,element in enumerate(md.mesh.elements):
			for index in element:
				fromelementthick[index-1]+=eltthick[i]*Tres.HydrologydcMaskEplactiveElt[i,0]
				Connect[index-1]+=1.*Tres.HydrologydcMaskEplactiveElt[i,0]
		fromelementthick[np.where(Connect>0)]=fromelementthick[np.where(Connect>0)]/Connect[np.where(Connect>0)]

		Q[t,:]=(Tres.EplHeadSlopeX.reshape(md.mesh.x.shape)*
						md.hydrology.epl_conductivity*fromelementthick*fromelement)



		# Q[t,:]=(Tres.EplHeadSlopeX.reshape(md.mesh.x.shape)
		# 				*md.hydrology.epl_conductivity*
		# 				Tres.HydrologydcMaskEplactiveNode.reshape(md.mesh.x.shape)*
		# 				Tres.HydrologydcEplThickness.reshape(md.mesh.x.shape))

	# 5. Close file
	NCData.close()
# }}}
# {{{ Two Layers last time
def netCDFTwoLayLastTime(md,filename):
	if path.exists(filename):
		print ('File {} allready exist'.format(filename))
		newname=raw_input('Give a new name or "delete" to replace: ')
		if newname=='delete':
			remove(filename)
		else:
			print ('New file name is {}'.format(newname))
			filename=newname

	NCData=Dataset(filename, 'w', format='NETCDF4')
	NCData.title = md.miscellaneous.name
	NCData.meshtype='unstructured'
	NCData.dimension = "2D" ;
	NCData.channels_on_edges = "no" ;
	NCData.institution = "Basile de Fleurian," ;
	NCData.source = "ISSM:svn-dev version on revision 21560" ;
	NCData.references = "http://shmip.bitbucket.io/"
	try:
		last= np.size(md.results.TransientSolution)
	except AttributeError:
		print 'No results presents, Exiting'
		NCData.close()
		return

	# 2. Define dimension
	nodenr       = NCData.createDimension('index1',md.mesh.numberofvertices) # number of nodes
	eltnr       = NCData.createDimension('index2',md.mesh.numberofelements) # number of elements
	dim          = NCData.createDimension('dim',md.mesh.dimension()) # number of dimensions
	time         = NCData.createDimension('time',None) # timesteps,

	# 3. define variables
	time       = NCData.createVariable('time','f8',('time'))
	nodes      = NCData.createVariable('coords1','f8',('dim', 'index1'))
	elements   = NCData.createVariable('coords2','f8',('dim', 'index2'))
	H          = NCData.createVariable('H','f8',('index1'))
	B          = NCData.createVariable('B','f8',('index1'))
	N          = NCData.createVariable('N','f8',('time','index1'))
	Ee         = NCData.createVariable('Ee','f8',('time','index1'))
	q          = NCData.createVariable('q','f8',('time','index2'))
	Q          = NCData.createVariable('Q','f8',('time','index1'))

	time.long_name = "time"
	time.units = "s"
	nodes.long_name = "node coordinates"
	nodes.units = "m"
	elements.long_name = "element midpoint"
	elements.units = "m"
	H.long_name = "ice thickness"
	H.units = "m"
	B.long_name = "bed elevation"
	B.units = "m"
	N.long_name = "effective pressure"
	N.units = "Pa"
	Ee.long_name = "EPL thickness"
	Ee.units = "m"
	q.long_name = "water sheet discharge"
	q.units = "m^2/s"
	Q.long_name = "channel discharge"
	Q.units = "m^2/s"


	# 4. Place data
	EltmidX=np.nanmean(md.mesh.x[md.mesh.elements-1],axis=1)
	EltmidY=np.nanmean(md.mesh.y[md.mesh.elements-1],axis=1)

	nodes[:,:]   = np.column_stack((md.mesh.x,md.mesh.y)).T
	elements[:,:]= np.column_stack((EltmidX,EltmidY)).T
	H[:]         = md.geometry.thickness
	B[:]         = md.geometry.base

	Tres=md.results.TransientSolution[-1]
	N[0,:]  = Tres.EffectivePressure.reshape(md.mesh.x.shape)
	Ee[0,:]	= Tres.HydrologydcEplThickness.reshape(md.mesh.x.shape)

	#need some computation for the fluxes
	sedgradx,sedgrady,sedgrad=slope(md,Tres.SedimentHead)
	q[0,:]=sedgradx*np.nanmean(md.hydrology.sediment_transmitivity)

	#get the element-wise thickness
	eltthick=np.mean(Tres.HydrologydcEplThickness[md.mesh.elements-1],axis=1)
	fromelement=np.zeros(np.shape(md.mesh.x))
	#gett the active elements on their nodes
	fromelement[md.mesh.elements[np.where(Tres.HydrologydcMaskEplactiveElt[:,0]==1),:]-1]=1

	fromelementthick=np.zeros(np.shape(md.mesh.x))
	Connect=np.zeros(np.shape(md.mesh.x))

	for i,element in enumerate(md.mesh.elements):
		for index in element:
			fromelementthick[index-1]+=eltthick[i]*Tres.HydrologydcMaskEplactiveElt[i,0]
			Connect[index-1]+=1.*Tres.HydrologydcMaskEplactiveElt[i,0]

	fromelementthick[np.where(Connect>0)]=fromelementthick[np.where(Connect>0)]/Connect[np.where(Connect>0)]

	Q[0,:]=(Tres.EplHeadSlopeX.reshape(md.mesh.x.shape)*
					md.hydrology.epl_conductivity*fromelementthick*fromelement)
	# Q[0,:]=(Tres.EplHeadSlopeX.reshape(md.mesh.x.shape)*md.hydrology.epl_conductivity*Tres.HydrologydcEplThickness.reshape(md.mesh.x.shape)*Tres.HydrologydcMaskEplactiveNode.reshape(md.mesh.x.shape))

	# 5. Close file
	NCData.close()
# }}}

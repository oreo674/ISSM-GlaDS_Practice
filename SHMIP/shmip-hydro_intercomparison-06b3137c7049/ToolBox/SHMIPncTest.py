#!/usr/bin/env python2
"""Checker for the netCDF files of the GHIP exercise.
You can run the checker by typing : ./GHIPncTest.py yourfile.nc

Your files should not return any ERRORS before sending your results
Some fields are not necessarly in your model so be aware of the WARNINGS
but don't try to make them vanish at all cost
"""

from netCDF4 import Dataset
import numpy as np
import sys
import os

def NCtest(NCFile):
	"""---NCtest---
	Description:
	open a netCDF file and go through the fields to check the compatibility to the experiments

	Needed Parameters:
	-netCDF file name

	Example of use :
	GHIPPlotter.NCtest('yourfile.nc')
	"""
	ErrFlag=0
	WarFlag=0
	# {{{ checking for file existence and format
	try:
		DatFile	 = Dataset(NCFile, mode='r')
	except RuntimeError:
		print("File '{}' does not exist or bad format".format(str(NCFile),))
		ErrFlag+=1
	# }}}
	# {{{ checking the CF related attributes
	CFattrdict={'title':' with the description of the simulation',
							'institution':' with your name and institution',
							'source':' with the model and version that produced the results',
							'references':' with shmip website'}

	for attr in CFattrdict:
		try:
			getattr(DatFile,attr)
		except AttributeError:
			print("attribute '{}', {} is missing from the nc file".format(attr,CFattrdict[attr]))
			ErrFlag+=1

	# }}}
	# {{{ checking the SHMIP related attributes and saving stuff
	attrdict={'meshtype':[' mesh caracteritics', ['structured','unstructured','unstructured Voronoi']],
						'channels_on_edges':[' presence of channels on edges ',['yes','no']]}
	for attr in attrdict:
		try:
			atval=getattr(DatFile,attr)
			if not str(atval) in attrdict[attr][1]:
				print("attribute '{}', stating {} should be one of the following : {}".format(attr,attrdict[attr][0],attrdict[attr][1]))
				ErrFlag+=1
		except AttributeError:
			print("attribute '{}', stating {} is missing from the nc file".format(attr,attrdict[attr][0]))
			ErrFlag+=1

	# }}}
	# {{{ treating the dimensions
	try:
		ischannelonedge=getattr(DatFile,'channels_on_edges')
	except AttributeError:
		ischannelonedge="no"
		print("attribute 'channels_on_edges' is missing, defaulting to 'no'")
	try:
		meshtype=getattr(DatFile,'meshtype')
	except AttributeError:
		meshtype="structured"
		print("attribute 'meshtype', is missing, defaulting to 'structured'")

	meshdim=0
	indexpool={'index1':'number of nodes',
						 'index2':'number of cells',
						 'index_ch':'number of channels'}
	dimdict={'dim':'spatial dimension',
					 'time':'time'}

	for key in DatFile.dimensions.keys():
		if key.startswith('index'):
			try:
				dimdict[str(key)]=indexpool[str(key)]
			except KeyError:
				print("Dimension '{}' in your netCDF is a typo for SHMIP available indexes are [{}]".format(key,indexpool.keys()))
				ErrFlag+=1
	if ischannelonedge=='yes':
		dimdict['index_ch']='number of channels'
		dimdict['n_nodes_ch']='number of nodes per channel'

	for key in dimdict:
		#Checking for presence and naming of the dimensions
		if key=='index_ch' and ischannelonedge=='no':
			print('There is a channel dimension but channels_on_edges is "no"')
			WarFlag+=1
		try:
			node = DatFile.dimensions[key]
			if key=='dim':
				meshdim=len(DatFile.dimensions['dim'])
		except KeyError:
			print("'{}', is an unknown dimension: \n  - check the spelling on the SHMIP website to control that it does not already exist under an other name \n  - contact Mauro and Basile if you need it added to the SHMIP pool".format(key))
			ErrFlag+=1
	# }}}
  # {{{ now dealing with variables
	dimpool={'coords2'   :['cell midpoint coordinates',u'm',(u'dim', u'index2')],
					 'coords_ch' :['channel midpoint coordinates',u'm',(u'dim', u'index_ch')],
					 'connect_ch':['channel connectivity',u'',(u'n_nodes_ch', u'index_ch')]}

	varpool={'h'     :['water sheet thickness',u'm',(u'time', u'index1')],
					 'hstore':['stored water effective layer thickness',u'm',(u'time', u'index1')],
					 'q'     :['water sheet discharge',u'm^2/s',(u'time', u'index2')],
					 'S'     :['channel cross-sectional area',u'm^2',(u'time', u'index_ch')],
					 'Q'     :['channel discharge',u'm^3/s',(u'time', u'index_ch')],
					 'Ee'    :['EPL thickness',u'm',(u'time', u'index1')]}

	vardict={'time'   :['time',u's',(u'time',)],
					 'coords1':['node coordinates',u'm',(u'dim', u'index1')],
					 'B'      :['bed elevation',u'm',(u'index1',)],
					 'H'      :['ice thickness',u'm',(u'index1',)],
					 'N'      :['effective pressure',u'Pa',(u'time', u'index1')]}

	if meshdim==1:
		vardict['W']=['domain width','m',(u'index1',)]

	for key in DatFile.dimensions.keys():
		if key.startswith('index'):
			if key=='index_ch':
				vardict['coords_ch']=dimpool['coords_ch']
				vardict['connect_ch']=dimpool['connect_ch']
			elif 'coords'+str(key[-1]) not in vardict.keys():
				vardict['coords'+str(key[-1])]=dimpool['coords'+str(key[-1])]

	if ischannelonedge=='yes':
		vardict['coords_ch']=dimpool['coords_ch']
		vardict['connect_ch']=dimpool['connect_ch']

	for key in DatFile.variables.keys():
		if key not in vardict.keys():
			try:
				vardict[str(key)]=varpool[str(key)]
			except KeyError:
				print("'{}', is an unknown variable: \n  - check the spelling on the SHMIP website to control that it does not already exist under an other name \n  - contact Mauro and Basile if you need it added to the SHMIP pool".format(key))
				ErrFlag+=1
	
	for key in vardict:
		#Checking for presence and naming of the variables
		try:
			val	= DatFile.variables[key][:]
			if type(val)!=np.ndarray:
				print("Bad format for '{}', type should be 'numpy.ndarray' but it is {}".format(key,str(type(val))))
				ErrFlag+=1
			if not DatFile.variables[key].dimensions==vardict[key][2]:
				print("Bad shape for '{}', the dimensions should be {}, not {}".format(key,vardict[key][2],DatFile.variables[key].dimensions))
				ErrFlag+=1
			try:
				varattr=DatFile.variables[key].getncattr('long_name')
				if not DatFile.variables[key].long_name==vardict[key][0]:
					print("Bad long_name for '{}', it should be '{}', not '{}'".format(key,vardict[key][0],DatFile.variables[key].long_name))
					ErrFlag+=1
			except AttributeError:
				print("Missing long_name for : {}".format(key,))
				ErrFlag+=1
			try:
				varattr=DatFile.variables[key].getncattr('units')
				if not str(DatFile.variables[key].units)==vardict[key][1]:
					print("Bad unit description for'{}', the unit should be '{}', not '{}'".format(key,vardict[key][1],DatFile.variables[key].units))
					ErrFlag+=1
			except AttributeError:
				print("Missing units for : {}".format(key,))
				ErrFlag+=1
		except KeyError:
				print("Bad name or missing variable for '{}', it should be '{}'".format(vardict[key][0],key))
				ErrFlag+=1
	# }}}

	return WarFlag,ErrFlag
	DatFile.close()

if __name__=="__main__":
	if len(sys.argv)==1: # print help if no filename is given:
		print __doc__
		sys.exit(0)
	WarFlag,ErrFlag=NCtest(sys.argv[1])
	if	ErrFlag==0:
		print 'Tests passed, your file should be fine.'
	if WarFlag>0:
		print 'WARNINGS do not necessarily need to be fixed'

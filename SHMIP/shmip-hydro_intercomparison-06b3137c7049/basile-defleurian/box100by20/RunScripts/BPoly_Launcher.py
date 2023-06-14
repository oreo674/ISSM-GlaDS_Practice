from model import *
from solve import *
import organizer as org
from exportVTK import exportVTK
from exportnetcdf_HGIP import netCDFTwoLay
from MoulinTLBatch import *
from export_netCDF import export_netCDF
from os import path

RunType='TwoLayPolyT'
#Uncomment the platform to go on
#clustername = gethostname()
clustername = 'vilje'
if clustername=='vilje':
	cluster=vilje('numnodes',1,'procspernodes',16,'time',30*60);#time in minutes
else:
	cluster=generic('name',clustername,'np',6)

# prefixList = ['B1','B2','B3','B4','B5']
# InputFile	 = ['B1_M.csv','B2_M.csv','B3_M.csv','B4_M.csv','B5_M.csv']

prefixList = ['B1rerere','B5rerere']
InputFile	 = ['B1_M.csv','B5_M.csv']
dt				 = 1.0e-5

# prefixList = ['B3rerere','B4rerere']
# InputFile	 = ['B3_M.csv','B4_M.csv']
# dt				 = 1.0e-6

SedPolyDeg = 2
SedThick	 = 10
EplThick	 = 5.0e-4
EplCol		 = 5.0e-7
EplKond		 = 9.0e1
Leakage		 = 1.0e-7

LaunchOrGet=raw_input('What do you want to do, launch (L), or retrieve (r): ')

for r,run in enumerate(prefixList):
	prefix=run
	MoulinFile=InputFile[r]
		
	RunName='-P'+str(SedPolyDeg)+'-es'+str(SedThick)+'-L'+str(Leakage)+'-Ke'+str(EplKond)+'-ee'+str(EplThick)+'-Col'+str(EplCol)
	Outname='../Models/'+prefix+RunName+'.nc'
	
	if path.exists(Outname):
		print('File {} allready exist, skipping'.format(Outname))
		continue
	
	if LaunchOrGet=='L':
		print 'Launching'
		#md,Org = PolyRunMo(SedPolyDeg,SedThick,EplThick,EplCol,EplKond,Leakage,MoulinFile,dt,prefix,RunName)
		md,Org=ReRunMo(MoulinFile,dt,prefix,RunName)
		md.cluster=cluster
		md = solve(md,'Transient','runtimename',0)
		
	elif LaunchOrGet=='r':
		print 'retreiving'
		#md,Org = PolyRunMo(SedPolyDeg,SedThick,EplThick,EplCol,EplKond,Leakage,MoulinFile,dt,prefix,RunName)
		md,Org=ReRunMo(MoulinFile,dt,prefix,RunName)
		md.cluster=cluster
		md=loadresultsfromcluster(md,prefix+RunName)
		netCDFTwoLay(md,'../Results/netCDF/'+prefix+RunName+'.nc')
		export_netCDF(md, '../Models/'+prefix+RunName+'.nc')
		filename=md.miscellaneous.name
		for extension in ['.bin','.queue','.toolkits']:
			try:
				os.remove(filename+extension)
			except OSError:
				print 'WARNING, no '+extension+'  is present for run '+filename
		
	else:
			print 'Bad Entry, Try again'

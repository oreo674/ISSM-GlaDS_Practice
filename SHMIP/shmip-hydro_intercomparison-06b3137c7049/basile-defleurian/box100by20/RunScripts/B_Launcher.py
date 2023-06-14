from model import *
from solve import *
import organizer as org
from exportVTK import exportVTK
from exportnetcdf_HGIP import netCDFTwoLay
from MoulinTLBatch import *
from os import path
from socket import gethostname
from export_netCDF import export_netCDF

RunType='TwoLay'
#Uncomment the platform to go on
#clustername = gethostname()
clustername = 'vilje'
if clustername=='vilje':
	cluster=vilje('numnodes',1,'procspernodes',8,'time',10*60);#time in minutes
else:
	cluster=generic('name',clustername,'np',8)

# prefixList=['B1','B2','B3','B4','B5']
# InputFile=['B1_M.csv','B2_M.csv','B3_M.csv','B4_M.csv','B5_M.csv']
# dt			 = 1.0e-4

# prefixList=['B1Sub','B2Sub','B3Sub','B4Sub','B5Sub']
# InputFile=['B1_M.csv','B2_M.csv','B3_M.csv','B4_M.csv','B5_M.csv']
# dt			 = 1.0e-4

prefixList=['B1Sub']
InputFile=['B1_M.csv']
dt			 = 1.0e-5

SedTrans = 0.0319
SedThick = 10
EplThick = 5.0e-4
EplCol	 = 5.0e-7
EplKond	 = 9.0e1
Leakage	 = 1.0e-7

LaunchOrGet=raw_input('What do you want to do, launch (L), or retrieve (r): ')

for r,run in enumerate(prefixList):
	prefix=run
	MoulinFile=InputFile[r]
		
	RunName='-Ts'+str(SedTrans)+'-es'+str(SedThick)+'-L'+str(Leakage)+'-Ke'+str(EplKond)+'-ee'+str(EplThick)+'-Col'+str(EplCol)
	Outname='../Models/'+prefix+RunName+'.nc'
	
	if path.exists(Outname):
		print('File {} allready exist, skipping'.format(Outname))
		continue
	
	if LaunchOrGet=='L':
		print 'Launching'
		#md,Org = RunMo(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,MoulinFile,dt,prefix,RunName)
		#md,Org = ReRunMo(MoulinFile,dt,prefix,RunName)
		md,Org = SubMo(MoulinFile,dt,prefix,RunName)
		md.cluster=cluster

		md = solve(md,'Transient','runtimename',0)
		
	elif LaunchOrGet=='r':
		print 'retreiving'
		#md,Org = RunMo(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,MoulinFile,dt,prefix,RunName)
		#md,Org = ReRunMo(MoulinFile,dt,prefix,RunName)
		md,Org = SubMo(MoulinFile,dt,prefix,RunName)
		md.cluster=cluster
		md=loadresultsfromcluster(md,prefix+RunName)
		export_netCDF(md,Outname)
		netCDFTwoLay(md,'../Results/netCDF/'+prefix+RunName+'.nc')
		filename=md.miscellaneous.name
		for extension in ['.bin','.queue','.toolkits']:
			try:
				os.remove(filename+extension)
			except OSError:
				print 'WARNING, no '+extension+'  is present for run '+filename
		
	else:
			print 'Bad Entry, Try again'

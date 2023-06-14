from model import *
from solve import *
import organizer as org
from exportVTK import exportVTK
from exportnetcdf_HGIP import netCDFTwoLay
from DiurnalTLBatch import *
from os import path
from export_netCDF import export_netCDF

RunType='TwoLay'
#Uncomment the platform to go on
#clustername = gethostname()
clustername = 'vilje'
if clustername=='vilje':
	cluster=vilje('numnodes',1,'procspernodes',8,'time',1*60);#time in minutes
else:
	cluster=generic('name',clustername,'np',6)

# prefixList=['C1','C2','C3','C4']
# RelAmpList=[0.25,0.5,1.0,2.0]

prefixList=['C1Sub','C2Sub','C3Sub','C4Sub']
RelAmpList=[0.25,0.5,1.0,2.0]

dt				 = 300.0/31536000. #five minutes timestep with hourly output
SedTrans	 = 0.0319
SedThick	 = 10
EplThick	 = 5.0e-4
EplCol		 = 5.0e-7
EplKond		 = 9.0e1
Leakage		 = 1.0e-7

LaunchOrGet=raw_input('What do you want to do, launch (L), or retrieve (r): ')

for r,run in enumerate(prefixList):
	prefix=run
	RelAmp=RelAmpList[r]

	RunName='-Ts'+str(SedTrans)+'-es'+str(SedThick)+'-L'+str(Leakage)+'-Ke'+str(EplKond)+'-ee'+str(EplThick)+'-Col'+str(EplCol)
	Outname='../Models/'+prefix+RunName+'.nc'
	
	if path.exists(Outname):
		print('File {} allready exist, skipping'.format(Outname))
		continue
	
	if LaunchOrGet=='L':
		print 'Launching'
		if path.exists('./'+prefix+RunName+'.queue'):
			print('File {} allready runing, skipping'.format(Outname))
			continue		
		#md,Org = RunDi(RelAmp,dt,prefix,RunName)
		#md,Org = ReRunDi(RelAmp,dt,prefix,RunName)
		md,Org = SubDi(RelAmp,dt,prefix,RunName)
		md.cluster=cluster
		md = solve(md,'Transient','runtimename',0)
		
	elif LaunchOrGet=='r':
		print 'retreiving'
		#md,Org = RunDi(RelAmp,dt,prefix,RunName)
		#md,Org = ReRunDi(RelAmp,dt,prefix,RunName)
		md,Org = SubDi(RelAmp,dt,prefix,RunName)
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

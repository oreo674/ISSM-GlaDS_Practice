from model import *
from solve import *
import organizer as org
from exportVTK import exportVTK
from exportnetcdf_HGIP import netCDFTwoLay
from SeasonalTLBatch import *
from os import path

RunType='TwoLayPolyT'
#Uncomment the platform to go on
#clustername = gethostname()
clustername = 'vilje'
if clustername=='vilje':
	cluster=vilje('numnodes',1,'procspernodes',16,'time',10*60)#time in minutes
else:
	cluster=generic('name',clustername,'np',6)


prefixList = ['D1','D2','D3','D4','D5']
TempList	 = [-4.0,-2.0,0.0,2.0,4.0]

dt				 = 3600.0/31536000.
SedPolyDeg = 3
SedThick	 = 10
EplThick	 = 5.0e-4
EplCol		 = 5.0e-7
EplKond		 = 12.5e1
Leakage		 = 1.0e-7

LaunchOrGet=raw_input('What do you want to do, launch (L), or retrieve (r): ')

for r,run in enumerate(prefixList):
	prefix=run
	Temp=TempList[r]
		
	RunName='-P'+str(SedPolyDeg)+'-es'+str(SedThick)+'-L'+str(Leakage)+'-Ke'+str(EplKond)+'-ee'+str(EplThick)+'-Col'+str(EplCol)
	Outname='../Models/'+prefix+RunName+'.python'
	
	if path.exists(Outname):
		print('File {} allready exist, skipping'.format(Outname))
		continue
	
	if LaunchOrGet=='L':
		print 'Launching'
		md,Org = RunSe(Temp,dt,prefix,RunName)
#		md,Org = ReRunSe(Temp,dt,prefix,RunName)
		md.cluster=cluster
		md = solve(md,'Transient','runtimename',0)
		
	elif LaunchOrGet=='r':
		print 'retreiving'
		md,Org = RunSe(Temp,dt,prefix,RunName)
#		md,Org = ReRunSe(dt,prefix,RunName)
		md.cluster=cluster
		md=loadresultsfromcluster(md,prefix+RunName)
		export_netCDF(md, '../Models/'+prefix+RunName+'.nc')
		netCDFTwoLay(md,'../Results/netCDF/'+prefix+RunName+'.nc')

		filename=md.miscellaneous.name
		for extension in ['.bin','.queue','.toolkits']:
			try:
				os.remove(filename+extension)
			except OSError:
				print 'WARNING, no '+extension+'  is present for run '+filename
		
	else:
		print 'Bad Entry, Try again'

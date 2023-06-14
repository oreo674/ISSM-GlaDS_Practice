from model import *
from solve import *
import organizer as org
from exportnetcdf_HGIP import netCDFTwoLay
from SteadyTLBatch import *
from os import path
from socket import gethostname
from export_netCDF import export_netCDF

RunType='TwoLay'
#Uncomment the platform to go on
#clustername = gethostname()
clustername = 'vilje'
if clustername=='vilje':
	cluster=vilje('numnodes',1,'procspernodes',8,'time',10.0*60)#time in minutes
else:
	cluster=generic('name',clustername,'np',6)

# prefixList=['A1Sub','A2Sub','A3Sub','A4Sub','A5Sub','A6Sub']
# InputList=[7.93e-11*31536000,
# 					 1.59e-09*31536000,
# 					 5.79e-09*31536000,
# 					 2.50e-08*31536000,
# 					 4.50e-08*31536000,
# 					 5.79e-07*31536000]
# DtVal				=	[1.0e-3]

prefixList=['A6test']
InputList=[5.79e-07*31536000]
DtVal				=	[3600.0/(12*31536000.)]

SedTransVal	=	[0.0319]
SedThickVal	=	[10]
EplThickVal	=	[5.0e-4]
EplKondVal	=	[9.0e1]
LeakageVal	=	[1.0e-7]
ColVal			=	[5.0e-7]

SedTransList = []
SedThickList = []
EplThickList = []
EplColList	 = []
EplKondList	 = []
LeakageList	 = []
DtList			 = []

for trans in SedTransVal:
	for thick in SedThickVal:
		for ethick in EplThickVal:
			for kond in EplKondVal:
				for leak in LeakageVal:
					for Dt in DtVal:
						for ecol in ColVal:
							SedTransList+=[trans]
							SedThickList+=[thick]
							EplThickList+=[ethick]
							EplColList+=[ecol]
							EplKondList+=[kond]
							LeakageList+=[leak]
							DtList+=[Dt]

RunNumber=len(EplThickList)
LaunchOrGet=raw_input('What do you want to do, launch (L), or retrieve (r): ')

for test in range(0,len(prefixList)):
	dt=DtVal[0]
	prefix=prefixList[test]
	Input=InputList[test]

	for run in range(0,RunNumber):
		SedTrans=SedTransList[run]
		SedThick=SedThickList[run]
		EplThick=EplThickList[run]
		EplCol=EplColList[run]
		EplKond=EplKondList[run]
		Leakage=LeakageList[run]
		
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
			md,Org = RunSt(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName)
			#md,Org = ReRunSt(Input,dt,prefix,RunName)
			#md,Org = SubSt(Input,dt,prefix,RunName)
			md.cluster=cluster
			md = solve(md,'Transient','runtimename',0)		
		elif LaunchOrGet=='r':
			print 'retreiving'
			md,Org = RunSt(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName)
			#md,Org = ReRunSt(Input,dt,prefix,RunName)
			#md,Org = SubSt(Input,dt,prefix,RunName)
			md.cluster=cluster
			md=loadresultsfromcluster(md,prefix+RunName)
			export_netCDF(md, Outname)
			netCDFTwoLay(md,'../Results/netCDF/'+prefix+RunName+'.nc')

			filename=md.miscellaneous.name
			for extension in ['.bin','.queue','.toolkits']:
				try:
					os.remove(filename+extension)
				except OSError:
					print 'WARNING, no '+extension+'  is present for run '+filename
		
		else:
			print 'Bad Entry, Try again'

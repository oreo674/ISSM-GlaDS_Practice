from model import *
from solve import *
import organizer as org
from socket import gethostname
from exportnetcdf_HGIP import netCDFTwoLay
from ValleySlope import *
from os import path
from export_netCDF import export_netCDF

RunType='TwoLay'
#Uncomment the platform to go on
#clustername = gethostname()
#clustername = 'vilje'
clustername = 'stallo'
if clustername=='vilje':
	cluster=vilje('numnodes',1,'procspernodes',8,'time',20*60)#time in minutes
elif clustername=='stallo':
	cluster=stallo('numnodes',1,'cpuspernode',8,'time',5*60)#time in minutes
else:
	cluster=generic('name',clustername,'np',6)

# prefixList=['E1newrerere','E2newrerere','E3newrerere','E4newrere','E5newrere']
# SlopeList=[0.05,0,-0.1,-0.5,-0.7]
# dtList=[1.0e-5,1.0e-5,1.0e-5,1.0e-5,1.0e-5]

# prefixList=['E4newrerere','E5newrerere']
# SlopeList=[-0.5,-0.7]
# dtList=[1.0e-5,1.0e-5]


prefixList=['E1stallo']
SlopeList=[0.05]
dtList=[1.0e-5]

SedTrans = 0.0319
SedThick = 10
EplThick = 5.0e-4
EplCol	 = 5.0e-7
EplKond	 = 9.0e1
Leakage	 = 1.0e-7
Input    = 1.158e-6*31536000

LaunchOrGet=raw_input('What do you want to do, mesh (m), launch (L), or retrieve (r): ')

for r,run in enumerate(prefixList):
	prefix=run
	Slope =SlopeList[r]
	dt=dtList[r]

	RunName='-Ts'+str(SedTrans)+'-es'+str(SedThick)+'-L'+str(Leakage)+'-Ke'+str(EplKond)+'-ee'+str(EplThick)+'-Col'+str(EplCol)
	Outname='../Models/'+prefix+RunName+'.nc'
	
	if path.exists(Outname):
		print('File {} allready exist, skipping'.format(Outname))
		continue
	if LaunchOrGet=='m':
		print 'Meshing'
		MeshName='../Models/'+prefix+'Valley'
		md,Org=Mesher(prefix,Slope)
		export_netCDF(md,MeshName+'.nc')

	elif LaunchOrGet=='L':
		print 'Launching'
		if path.exists('./'+prefix+RunName+'.queue'):
			print('File {} allready runing, skipping'.format(Outname))
			continue
		md,Org = ValleySlope(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName)
		#md,Org = ValleyReSlope(Input,dt,prefix,RunName)
		md.cluster=cluster
		md = solve(md,'Transient','runtimename',0)
		
	elif LaunchOrGet=='r':
		print 'retreiving'
		md,Org = ValleySlope(SedTrans,SedThick,EplThick,EplCol,EplKond,Leakage,Input,dt,prefix,RunName)
		#md,Org = ValleyReSlope(Input,dt,prefix,RunName)
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

import numpy
from model import *
from solve import *
from MatlabFuncs import *
from loadmodel import *
import organizer as org
from exportVTK import exportVTK
from exportnetcdf_HGIP import netCDFTwoLay
from SteadyTLBatch import SteadyTLBatch
from os import path

prefix='bdef_4bfit'
RunType='TwoLay'
Input=3.0*0.365

SedTransVal=[2.9e-2]
SedThickVal=[10]
DtVal=[1.0e-3]

EplThickVal=[5.0e-4]
EplKondVal=[8.5e1,8.75e1,9.0e1,9.25e1,9.5e1]
LeakageVal=[1.5e-8,2.0e-8,2.5e-8]

SedTransList=[]
SedThickList=[]
EplThickList=[]
EplKondList=[]
LeakageList=[]
DtList=[]

for trans in SedTransVal:
	for thick in SedThickVal:
		for ethick in EplThickVal:
			for kond in EplKondVal:
				for leak in LeakageVal:
					for Dt in DtVal:
						SedTransList+=[trans]
						SedThickList+=[thick]
						EplThickList+=[ethick]
						EplKondList+=[kond]
						LeakageList+=[leak]
						DtList+=[Dt]

RunNumber=len(SedTransList)

for run in range(0, RunNumber):
	SedTrans=SedTransList[run]
	SedThick=SedThickList[run]
	EplThick=EplThickList[run]
	EplCol=EplThick*1.0e-3
	EplKond=EplKondList[run]
	Leakage=LeakageList[run]
	dt=DtList[run]

	RunName='-Ts'+str(SedTrans)+'-es'+str(SedThick)+'-L'+str(Leakage)+'-Ke'+str(EplKond)+'-ee'+str(EplThick)+'FClong'
	Outname='../Models/'+prefix+RunName

	md = loadmodel(Outname)
		
	outname='../Results/'+prefix+RunName
	netCDFTwoLay(md,'../Results/netCDF/'+prefix+RunName+'.nc')

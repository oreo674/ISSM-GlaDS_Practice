import numpy
from model import *
from loadmodel import *
from exportVTK import exportVTK

RunType='TwoLay'
prefix='B4_fit'
InputList=[1.095]

DtVal=[1.0e-3]
SedTransVal=[2.7e-2]
SedThickVal=[10]

EplThickVal=[1.0e-3,5.0e-4,1.0e-4]
EplKondVal=[1.0e1,5.0e1,10.0e1]
LeakageVal=[1.0e-8,5.0e-8,10.0e-8]

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

	RunName='-Ts'+str(SedTrans)+'-es'+str(SedThick)+'-L'+str(Leakage)+'-Ke'+str(EplKond)+'-ee'+str(EplThick)+'-Col'+str(EplCol)
	Outname='../Models/'+prefix+RunName

	md = loadmodel(Outname)
		
	outname='../Results/VTK/'+prefix+RunName
	exportVTK(outname,md,'mesh','geometry')

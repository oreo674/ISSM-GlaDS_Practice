from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob as gb
import re
import pickle
import collections

def FormatGlobal(NCFile):

	FileList=gb.glob(NCFile)
	ParamDict={}
	#get the different parameter vallue for K and L
	for i,File in enumerate(FileList):
		print File
		ParamList=re.split('(?<!e)[-]|Cont.nc', File)[1:-1]
		for param in ParamList:
			value=float(''.join(re.split('(>?\d)',param)[1:]))
			name=re.split('(>?\d)',param)[0]
			if i==0:
				ParamDict[name]=[value]
			else:
				ParamDict[name]+=[value]

		DatFile	 = Dataset(File, mode='r')
		MeanN    = getattr(DatFile,'meanN')
		MaxN     = getattr(DatFile,'maxN')
		DatFile.close()

		if i==0:
			ParamDict['meanN']=[MeanN]
			ParamDict['maxN']=[MaxN]
		else:			
			ParamDict['meanN']+=[MeanN]
			ParamDict['maxN']+=[MaxN]

	SaveName='../Results/AllRun'+re.split('/',re.split('(?<!e)[-]|.nc', File)[0])[-1]
	with open(SaveName,'w') as handle:
		pickle.dump(ParamDict,handle)

		
def PlotAllSensib(FileName):

	with open(FileName,'r') as handle:
		VarDict=pickle.loads(handle.read())

	fig = plt.figure()
	FigTitle='Sensibility for '+ re.split('/',FileName)[-1]
	fig.suptitle(FigTitle,fontsize=20)
	ParamKeys=['Ke', 'ee', 'Ts', 'L', 'Col', 'es']
	StdValDict={}
	for Key in ParamKeys:
		StdValDict[Key]=collections.Counter(VarDict[Key]).most_common(1)[0][0]
	
	for i,Key in enumerate(ParamKeys):
		index=range(0,len(VarDict['Ts']))
		for Param in ParamKeys:
			if Param != Key:
				index=[j for j in index if VarDict[Param][j]==StdValDict[Param]]

		ax=plt.subplot(2,3,i+1)
		ax.plot([VarDict[Key][i] for i in index],[VarDict['meanN'][i] for i in index] ,color = 'r',marker='.',markersize=10,linestyle=' ')
		Title='Sensibility to '+Key
		ax.set_title(Title, fontsize=18)
	#And draw
	plt.show(block=False)


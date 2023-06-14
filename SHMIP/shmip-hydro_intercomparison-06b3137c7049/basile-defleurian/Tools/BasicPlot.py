from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob as gb
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from TreatMesh import *
from collections import defaultdict
import re

# {{{ PlotMap(NCFile, DataName):

def PlotMap(NCFile, DataName):

	DataNameDict={'N':['Effective Pressure','Pressure [Pa]'],
								'Hsed':['Inneficient sytem head','Water head [m]'],
								'H':['Topo','Ice Thickness [m]'],
								'B':['Topo','Bedrock elevation [m]'],
								'Hepl':['Efficicient system head','Water head [m]'],
								'EplThick':['Efficient system thickness','Thickness [m]']}
	DataString=DataNameDict[DataName][0]
	Units=DataNameDict[DataName][1]

	FileList=gb.glob(NCFile)
	for i,File in enumerate(FileList):
		print File
		DatFile	= Dataset(File, mode='r')
		Xs			= DatFile.variables['coords1'][0,:]
		Ys			= DatFile.variables['coords1'][1,:]
		try:
			data  = DatFile.variables[DataName][-1,:]
		except IndexError:
			data  = DatFile.variables[DataName][0,:]
		DatFile.close()

		#Reshaping data
		lims=[data.min(),data.max()]

		#Now Plotting
		fig = plt.figure(tight_layout=False)
		ax	= fig.add_subplot(111)
		ax2 = fig.add_axes([0.1, 0.1, 0.85, 0.05])

		if DataName=='EplThick':
			norm = mpl.colors.LogNorm(vmin=lims[0], vmax=lims[1])#(vmin=1.0e-7, vmax=1.0)
		else:
			norm = mpl.colors.Normalize(vmin=lims[0], vmax=lims[1])
		cmap	 = plt.cm.coolwarm #plt.cm.RdBu
		cbar	 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='horizontal')

		ax.tricontourf(Xs, Ys, elts.T, data, 128, cmap=cmap,norm=norm)
		ax.set_axes(ax)
		ax.axis('scaled')
		ax.set_xlabel('X axis [km]', fontsize=14)
		ax.set_ylabel('Y axis [km]', fontsize=14)

		Title=DataString+' \n for run \n'+'-'.join(File.split('-')[1:])
		CBTitle=Units
		ax.set_title(Title, fontsize=18)
		ax2.set_title(CBTitle, fontsize=18)

		#And draw
		plt.show(block=False)

# }}}
# {{{ SectionPlot(NCFile, DataName):
def SectionPlot(NCFile, DataName):

	DataNameDict = {'N':['Effective Pressure','Pressure [Pa]'],
								'Hsed':['Inneficient sytem head','Water head [m]'],
								'Hepl':['Efficicient system head','Water head [m]'],
								'EplThick':['Efficient system thickness','Thickness [m]']}
	DataString	 = DataNameDict[DataName][0]
	Units				 = DataNameDict[DataName][1]
	FileList		 = gb.glob(NCFile)

	for File in FileList:
		compare=False
		runtype=re.split('([A-E][1-6])',File)[-2]
		if DataName=='N':
			CompFile = Dataset('/home/bfl022/Dropbox/GHIP_Results/mwer/{}_mwer.nc'.format(runtype), mode='r')
			CompXs	 = CompFile.variables['coords1'][0,:]
			try:
				Compdata = CompFile.variables[DataName][-1,:]
			except IndexError:
				Compdata = CompFile.variables[DataName][0,:]
			CompFile.close()
			compare		 = True
		print File
		DatFile	= Dataset(File, mode='r')
		Xs			= DatFile.variables['coords1'][0,:]
		try:
			data  = DatFile.variables[DataName][-1,:]
		except IndexError:
			data  = DatFile.variables[DataName][0,:]

		DatFile.close()

		#Reshaping data
		lims=[data.min(),data.max()]

		#Now Plotting
		fig = plt.figure(tight_layout=True)
		ax = fig.add_subplot(111)

		if compare:
			ax.plot(CompXs, Compdata,'r.')
		ax.plot(Xs, data,'b.')
		ax.set_axes(ax)
		ax.set_xlabel('Distance from front [m]', fontsize=14)
		ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
		ax.set_ylabel('Effective Pressure [Pa]', fontsize=14)
		Title=DataString+'\n for Run \n'+re.split('/|\.nc',File)[-2]
		ax.set_title(Title, fontsize=18)

		#And draw
		plt.show(block=False)

# }}}
# {{{ SectionStamps(NCFile, DataName):

def SectionStamps(NCFile, DataName):

	DataNameDict = {'N':['Effective Pressure','Pressure [Pa]'],
								'Hsed':['Inneficient sytem head','Water head [m]'],
								'Hepl':['Efficicient system head','Water head [m]'],
								'EplThick':['Efficient system thickness','Thickness [m]']}
	DataString	 = DataNameDict[DataName][0]
	Units				 = DataNameDict[DataName][1]

	FileList	 = gb.glob(NCFile)
	LeakDict	 = {}
	KondDict	 = {}
	ThickDict	 = {}
	coldict		 = {}
	rowdict		 = {}
	leakindex	 = 0
	kondindex	 = 0
	thickindex = 0

	runtype=re.split('([A-E][1-6])',FileList[0])[-2]
	if DataName=='N':
		CompFile = Dataset('/home/bfl022/Dropbox/GHIP_Results/mwer/{}_mwer.nc'.format(runtype), mode='r')
		CompXs	 = CompFile.variables['coords1'][0,:]
		Compdata = CompFile.variables[DataName][:,:]
		try:
			Compdata = CompFile.variables[DataName][-1,:]
		except IndexError:
			Compdata = CompFile.variables[DataName][0,:]
		CompFile.close()

	#get the different parameter vallue for Ke and L
	for File in FileList:
		Leakage=float(re.split('L|(?<!e)[-]|new', File)[4])
		try:
			dummy=LeakDict[Leakage]
		except KeyError:
			LeakDict[Leakage]=leakindex
			leakindex+=1
		EplKond=float(re.split('Ke|(?<!e)[-]|new', File)[5])
		try:
			dummy=KondDict[EplKond]
		except KeyError:
			KondDict[EplKond]=kondindex
			kondindex+=1
		EplThick=float(re.split('ee|(?<!e)[-]|new', File)[6])
		try:
			dummy=ThickDict[EplThick]
		except KeyError:
			ThickDict[EplThick]=thickindex
			thickindex+=1
	DimDict={str(len(ThickDict))+'ThickDict':[len(ThickDict),sorted(ThickDict.keys()),'ee|-Col'],
					 str(len(KondDict))+'KondDict':[len(KondDict),sorted(KondDict.keys()),'Ke|-ee'],
					 str(len(LeakDict))+'LeakDict':[len(LeakDict),sorted(LeakDict.keys()),'L|-Ke']}

	if len(DimDict)>2:
		lenlist=[int(DimDict.keys()[i][0]) for i in [0,1,2]]
		for key in DimDict.keys():
			if key.startswith(str(min(lenlist))):
				del DimDict[key]
				break
	cols=DimDict[DimDict.keys()[0]][0]
	rows=DimDict[DimDict.keys()[1]][0]

	#get the index right to get sorted subplots
	for i,value in enumerate(DimDict[DimDict.keys()[0]][1]):
		coldict[value]=i
	for i,value in enumerate(DimDict[DimDict.keys()[1]][1]):
		rowdict[value]=i

	fig, axar = plt.subplots(cols, rows, sharex=True, sharey=True)
	FigTitle=DataString+' for Runs '+'-'.join(File.split('-')[0:3]).split('/')[-1]
	fig.suptitle(FigTitle,fontsize=20)

	for i,File in enumerate(FileList):
		print File
		DatFile	= Dataset(File, mode='r')
		Xs			= DatFile.variables['coords1'][0,:]
		try:
			data  = DatFile.variables[DataName][-1,:]
		except IndexError:
			data  = DatFile.variables[DataName][0,:]
		DatFile.close()

		colval=float(re.split(DimDict[DimDict.keys()[0]][2], File)[1])
		rowval=float(re.split(DimDict[DimDict.keys()[1]][2], File)[1])

		#Reshaping data
		lims=[data.min(),data.max()]
		axar[coldict[colval],rowdict[rowval]].set_xlim([0,10.0e4])
		axar[coldict[colval],rowdict[rowval]].set_ylim([Compdata.min(),Compdata.max()])
		if DataName=='N' and runtype=='A5':
			axar[coldict[colval],rowdict[rowval]].plot(CompXs, Compdata,'r.')
		axar[coldict[colval],rowdict[rowval]].plot(Xs, data,'b.')
		axar[coldict[max(coldict.keys())],rowdict[rowval]].set_xlabel('X axis [m]', fontsize=14)
		axar[coldict[colval],0].set_ylabel(DataString+' [Pa]', fontsize=14)
		Title='-'.join(File.split('-')[3:]).split('.nc')[0]
		axar[coldict[colval],rowdict[rowval]].set_title(Title, fontsize=18)
		#And draw
		plt.show(block=False)

# }}}
# {{{ GetEnveloppe(NCFile, DataName):

def GetEnveloppe(NCFile, DataName):

	DataNameDict = {'N':['Effective Pressure','Pressure [Pa]'],
								'Hsed':['Inneficient sytem head','Water head [m]'],
								'Hepl':['Efficicient system head','Water head [m]'],
								'EplThick':['Efficient system thickness','Thickness [m]']}
	DataString	 = DataNameDict[DataName][0]
	Units				 = DataNameDict[DataName][1]

	FileList	=	gb.glob(NCFile)
	LeakDict	=	{}
	KondDict	=	{}
	leakindex	=	0
	kondindex	=	0

	compare=False
	runtype=re.split('([A-E][1-6])',FileList[0])[-2]
	if DataName=='N' and runtype=='4B':
		compare		 = True
		RunIndex	 = NCFile.split('_')[1].split('fit')[0]
		CompFile    = Dataset('/home/bfl022/Dropbox/GHIP_Results/mwer/{}_mwer.nc'.format(runtype), mode='r')
		CompXs		 = CompFile.variables['xy'][0,:]
		try:
			Compdata = CompFile.variables[DataName][-1,:]
		except IndexError:
			Compdata = CompFile.variables[DataName][0,:]
		CompFile.close()

	#get the different parameter vallue for K and L
	for File in FileList:
		Leakage=float(re.split('L|(?<!e)[-]|new', File)[4])
		try:
			dummy=LeakDict[Leakage]
		except KeyError:
			LeakDict[Leakage]=leakindex
			leakindex+=1
		EplKond=float(re.split('Ke|(?<!e)[-]|new', File)[5])
		try:
			dummy=KondDict[EplKond]
		except KeyError:
			KondDict[EplKond]=kondindex
			kondindex+=1

		#get the index right to get sorted subplots
		for i,value in enumerate(sorted(LeakDict.keys())):
			LeakDict[value]=i
		for i,value in enumerate(sorted(KondDict.keys())):
			KondDict[value]=i

	fig, axar = plt.subplots(len(KondDict), len(LeakDict), sharex=True, sharey=True,figsize=(len(KondDict)*3, len(LeakDict)*2))
	FigTitle=DataString+' for Runs '+'-'.join(File.split('-')[0:3]).split('/')[-1]
	fig.suptitle(FigTitle,fontsize=20)
	for i,File in enumerate(FileList):
		print File
		DatFile	= Dataset(File, mode='r')
		Xs			= DatFile.variables['xy'][0,:]
		try:
			data  = DatFile.variables[DataName][-1,:]
		except IndexError:
			data  = DatFile.variables[DataName][0,:]
		MeanN   = getattr(DatFile,'meanN')
		MaxN    = getattr(DatFile,'maxN')
		DatFile.close()

		Leakage=float(re.split('L|(?<!e)[-]|new', File)[4])
		EplKond=float(re.split('Ke|(?<!e)[-]|new', File)[5])

		valtab		=	np.vstack((Xs,data)).T
		sortedtab	=	valtab[valtab[:,0].argsort()] #sorting function of Xs
		step			=	1000
		UpEnv			=	np.zeros((60000/step,2))
		LowEnv		=	np.zeros((60000/step,2))
		MeanEnv		=	np.zeros((60000/step,2))

		for j,value in enumerate(range(0,60000,step)):
			Slice=sortedtab[np.where(sortedtab[:,0]<=value),1]
			try:
				UpEnv[j,:]	 = sortedtab[np.argmax(Slice),:]
				LowEnv[j,:]	 = sortedtab[np.argmin(Slice),:]
				MeanEnv[j,:] = np.mean(sortedtab[0:np.shape(Slice)[1]-1,:],axis=0)
				sortedtab		 = sortedtab[np.shape(Slice)[1]:,:]
			except ValueError:
				print('No value skipping step'+str(j)+' for X less than '+str(value))

		#Reshaping data
		lims=[data.min(),data.max()]
		if compare:
			axar[KondDict[EplKond],LeakDict[Leakage]].plot(CompXs, Compdata,color = '0.75',marker='.',linestyle=' ')

		axar[KondDict[EplKond],LeakDict[Leakage]].set_xlim([0,6.0e4])
		axar[KondDict[EplKond],LeakDict[Leakage]].set_ylim([data.min(),max(Compdata.max(),data.max())])
		axar[KondDict[EplKond],LeakDict[Leakage]].plot(UpEnv[:,0], UpEnv[:,1],'b-',linewidth=1.5)
		axar[KondDict[EplKond],LeakDict[Leakage]].plot(LowEnv[:,0], LowEnv[:,1],'b-',linewidth=1.5)
		axar[KondDict[EplKond],LeakDict[Leakage]].plot(MeanEnv[:,0], MeanEnv[:,1],'r-',linewidth=1.5)

#		axar[KondDict[EplKond],LeakDict[Leakage]].set_axes(axar[0,0])
		axar[KondDict[max(KondDict.keys())],LeakDict[Leakage]].set_xlabel('X axis [m]', fontsize=14)
		axar[KondDict[EplKond],0].ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
		axar[KondDict[EplKond],0].set_ylabel(Units, fontsize=14)
		Title='-'.join(File.split('-')[3:]).split('.nc')[0]
		axar[KondDict[EplKond],LeakDict[Leakage]].set_title(Title, fontsize=18)
		# axar[KondDict[EplKond],LeakDict[Leakage]].text(7.0e3,6.0e5, 'MeanN='+str(MeanN))
		# axar[KondDict[EplKond],LeakDict[Leakage]].text(7.0e3,2.0e5, 'MaxN='+str(MaxN))
		#And draw
		plt.show(block=False)

	fig.savefig('../Figures/FitComp'+Title+'.pdf',format='pdf')

# }}}
# {{{ SensibEnveloppe(NCFile, DataName):

def SensibEnveloppe(NCFile, DataName, Param):

	DataNameDict = {'N':['Effective Pressure','Pressure [Pa]'],
								'Hsed':['Inneficient sytem head','Water head [m]'],
								'Hepl':['Efficicient system head','Water head [m]'],
								'EplThick':['Efficient system thickness','Thickness [m]']}
	DataString	 = DataNameDict[DataName][0]
	Units				 = DataNameDict[DataName][1]

	FileList	 = gb.glob(NCFile)
	ParamDict	 = {}
	paramindex = 0
	#get the different parameter vallue for K and L
	for File in FileList:
		ParamVal=float(re.split('(?<!e)[-]|.nc',re.split(Param, File)[1])[0])
		try:
			dummy=ParamDict[ParamVal]
		except KeyError:
			ParamDict[ParamVal]=paramindex
			paramindex+=1
		#get the index right to get sorted subplots
		for i,value in enumerate(sorted(ParamDict.keys())):
			ParamDict[value]=i

	fig				= plt.figure()
	ax				= fig.add_subplot(111)
	jet				= cm = plt.get_cmap('jet')
	cNorm			= mpl.colors.Normalize(vmin=0, vmax=len(FileList))
	scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=jet)

	for i,File in enumerate(FileList):
		ParamVal=float(re.split('(?<!e)[-]|.nc',re.split(Param, File)[1])[0])
		SortedFile=FileList[ParamDict[ParamVal]]
		ParamVal=float(re.split('(?<!e)[-]|.nc',re.split(Param, SortedFile)[1])[0])
		print SortedFile

		DatFile	= Dataset(SortedFile, mode='r')
		Xs			= DatFile.variables['xy'][0,:]
		try:
			data  = DatFile.variables[DataName][-1,:]
		except IndexError:
			data  = DatFile.variables[DataName][0,:]
		MeanN   = getattr(DatFile,'meanN')
		MaxN    = getattr(DatFile,'maxN')
		DatFile.close()

		valtab		=	np.vstack((Xs,data)).T
		sortedtab	=	valtab[valtab[:,0].argsort()] #sorting function of Xs
		step			=	1000
		UpEnv			=	np.zeros((60000/step,2))
		LowEnv		=	np.zeros((60000/step,2))
		MeanEnv		=	np.zeros((60000/step,2))

		for j,value in enumerate(range(0,60000,step)):
			Slice=sortedtab[np.where(sortedtab[:,0]<=value),1]
			try:
				UpEnv[j,:]	 = sortedtab[np.argmax(Slice),:]
				LowEnv[j,:]	 = sortedtab[np.argmin(Slice),:]
				MeanEnv[j,:] = np.mean(sortedtab[0:np.shape(Slice)[1]-1,:],axis=0)
				sortedtab		 = sortedtab[np.shape(Slice)[1]:,:]
			except ValueError:
				print('No value skipping step'+str(j)+' for X less than '+str(value))

		#Reshaping data
		lims=[data.min(),data.max()]
		#Now Plotting
		colorVal = scalarMap.to_rgba(i)
		ax.plot(MeanEnv[:,0], MeanEnv[:,1],color=colorVal,linestyle='-',linewidth=1.5,label=Param+'_'+str(ParamVal))
		ax.set_axes(ax)
		ax.set_xlabel('X axis [m]', fontsize=14)
		ax.set_ylabel('Effective Pressure [Pa]', fontsize=14)
		Title=DataString+'\n for Run \n'+'-'.join(SortedFile.split('-')[1:])
		ax.set_title(Title, fontsize=18)
		ax.legend(loc='center left',bbox_to_anchor=(0.1, 0.2))
		#And draw
		plt.show(block=False)

	fig.savefig('../Figures/Sensib'+Param+'.pdf',format='pdf')

# }}}
# {{{ FitPlot(NCFile, DataName):

def FitPlot(NCFile, DataName):

	DataNameDict = {'N':['Effective Pressure','Pressure [MPa]'],
								'Hsed':['Inneficient sytem head','Water head [m]'],
								'Hepl':['Efficicient system head','Water head [m]'],
								'EplThick':['Efficient system thickness','Thickness [m]']}
	DataString	 = DataNameDict[DataName][0]
	Units				 = DataNameDict[DataName][1]
	FileList		 = gb.glob(NCFile)

	compare	=	False
	runtype	=	re.split('([A-E][1-6])',FileList[0])[-2]
	rundict	=	defaultdict(int)
	xstep=100
	ystep=20
	grid_x, grid_y = np.mgrid[0:xstep:101j, 0:ystep:21j]
	Weights=GetGridWeight(grid_x,grid_y)

	for File in FileList:
		print File
		rundict[File]=re.split('(fit-Ts|-es)',File)[2]

	if DataName=='N':
		CompFile    = Dataset('/home/bfl022/Dropbox/GHIP_Results/mwer/{}_mwer.nc'.format(runtype), mode='r')
		compare     = True
		Compcoords  = CompFile.variables['coords1'][:,:]
		Compdata    = 1.0e-6* CompFile.variables[DataName][-1,:]
		CompFile.close()
		griddeddata	=	griddata(1.0e-3*Compcoords.T, Compdata[:], (grid_x, grid_y), method='linear')
		Compmean    = np.sum(griddeddata*Weights)/np.sum(Weights)
		sorter = Compcoords[0,:].argsort()
		Compdata = Compdata[sorter[::-1]]
		CompX = Compcoords[0,sorter[::-1]]

	fig = plt.figure()
	gs = mpl.gridspec.GridSpec(2,3)
	gs.update(left=0.11, right=0.97, top=0.97, bottom= 0.13, wspace=0.2, hspace=0.25)
	ax = plt.subplot(gs[0,:])
	jet = cm = plt.get_cmap('jet')
	cNorm  = mpl.colors.Normalize(vmin=0, vmax=len(FileList))
	scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=jet)

	MeansN  = np.zeros((2,len(FileList)))
	mean  = np.zeros((2,len(FileList)))
	for i,File in enumerate(sorted(rundict,key=rundict.get)):
		print File
		DatFile	= Dataset(File, mode='r')
		coords  = DatFile.variables['coords1'][:,:]
		data		= 1.0e-6*DatFile.variables[DataName][-1,:]
		DatFile.close()

		griddeddata	=	griddata(1.0e-3*coords.T, data[:], (grid_x, grid_y), method='linear')
		mean[1,i]   = np.sum(griddeddata*Weights)/np.sum(Weights)
		TsVal       = rundict[File]
		mean[0,i]   = float(TsVal)
		sorter = coords[0,:].argsort()
		data = data[sorter[::-1]]
		X = coords[0,sorter[::-1]]


		#Reshaping data
		lims=[data.min(),data.max()]

		#Now Plotting
		colorVal = scalarMap.to_rgba(i)
		ax.plot(X*1.0e-3,data,color=colorVal,label='$Ts=$'+TsVal+' $ms^{-1}$',linewidth=2)
		ax.set_axes(ax)
		ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
		ax.set_xlabel('Distance from front [km]', fontsize=14)
		ax.set_ylabel('Effective Pressure \n [MPa]', fontsize=14)

		#And draw
	if compare:
		ax.plot(CompX*1.0e-3, Compdata,'k--',label='Reference \n Simulation',linewidth=2)
	ax.legend(bbox_to_anchor=(1.0, -0.215),ncol=3,columnspacing=0.5)

	ay = plt.subplot(gs[1,0])
	for i,File in enumerate(FileList):
		colorVal = scalarMap.to_rgba(i)
		ay.plot(mean[0,i],mean[1,i],color=colorVal,marker='.',markersize=15)
		ay.set_xlim([1.0e-2,4.0e-2])
		ay.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
		ay.set_xlabel('Sediment Transmitivity [ms$^{-1}$]', fontsize=14)
		ay.get_xaxis().set_label_coords(0.5,-0.17)
		ay.set_ylabel('Mean effective pressure \n [MPa]', fontsize=14)
		ay.plot([2.0e-2,3.6e-2],[Compmean,Compmean],color='k',linestyle='--')
	fig.savefig('../Figures/RunFit'+runtype+'.pdf',format='pdf')
	plt.show(block=False)


# }}}
# {{{ SensibPlot(NCFile,DataName)

def SensibPlot(NCFile,DataName):

	DataNameDict = {'N':['Effective Pressure','Pressure [MPa]'],
								'Hsed':['Inneficient sytem head','Water head [m]'],
								'Hepl':['Efficicient system head','Water head [m]'],
								'EplThick':['Efficient system thickness','Thickness [m]']}
	DataString	 = DataNameDict[DataName][0]
	Units				 = DataNameDict[DataName][1]

	FileList		 = gb.glob(NCFile)
	LeakDict		 = {}
	KondDict		 = {}
	ThickDict		 = {}
	leakindex		 = 0
	kondindex		 = 0
	thickindex	 = 0

	#get the different parameter vallue for K and L
	for File in FileList:
		Leakage=float(re.split('L|(?<!e)[-]|new', File)[4])
		try:
			dummy=LeakDict[Leakage]
		except KeyError:
			LeakDict[Leakage]=leakindex
			leakindex+=1
		EplKond=float(re.split('Ke|(?<!e)[-]|new', File)[5])
		try:
			dummy=KondDict[EplKond]
		except KeyError:
			KondDict[EplKond]=kondindex
			kondindex+=1
		EplThick=float(re.split('ee|(?<!e)[-]|new', File)[6])
		try:
			dummy=ThickDict[EplThick]
		except KeyError:
			ThickDict[EplThick]=thickindex
			thickindex+=1

		#get the index right to get sorted subplots
		for i,value in enumerate(sorted(LeakDict.keys())):
			LeakDict[value]=i
		for i,value in enumerate(sorted(KondDict.keys())):
			KondDict[value]=i
		for i,value in enumerate(sorted(ThickDict.keys())):
			ThickDict[value]=i

	compare=False
	runtype=re.split('([A-E][1-6])',FileList[0])[-2]
	if DataName=='N':
		CompFile    = Dataset('/home/bfl022/Dropbox/GHIP_Results/mwer/{}_mwer.nc'.format(runtype), mode='r')
		compare       = True
		NumofNodes    = len(CompFile.dimensions['nodenr'])
		Compdata      = np.zeros((2,NumofNodes))
		Compdata[0,:] = CompFile.variables['xy'][0,:]
		Compdata[1,:] = 1.0e-6* CompFile.variables[DataName][:]
		Ys			      = CompFile.variables['xy'][1,:]
		elts		      = CompFile.variables['connect'][:]
		Compmean      = ComputeMean(elts,Compdata[0,:],Ys,Compdata[1,:])
		Compdata      = Compdata[:,Compdata[0,:].argsort()]
		CompFile.close()

	step =100
	window=2000

	NmaxPosPar   = np.NaN*np.ones((len(LeakDict), len(KondDict),len(ThickDict)))
	MaxDrainPar  = np.NaN*np.ones((len(LeakDict), len(KondDict),len(ThickDict)))
	MeanDrainPar = np.NaN*np.ones((len(LeakDict), len(KondDict),len(ThickDict)))
	Param        = np.NaN*np.ones((3,len(FileList)))
	NmaxPar      = np.NaN*np.ones((len(LeakDict), len(KondDict),len(ThickDict)))

	for i,File in enumerate(FileList):
		print File, i
		DatFile	   = Dataset(File, mode='r')
		NumofNodes = len(DatFile.dimensions['nodenr'])
		Xs			   = DatFile.variables['xy'][0,:]
		data       = DatFile.variables[DataName][:]
		DatFile.close()

		Param[0,i]=float(re.split('L|(?<!e)[-]|new', File)[4])
		Param[1,i]=float(re.split('Ke|(?<!e)[-]|new', File)[5])
		Param[2,i]=float(re.split('ee|(?<!e)[-]|new', File)[6])

		if DataName=='N':
			data=1.0e-6*data

		valtab=np.vstack((Xs,data)).T
		sortedtab=valtab[valtab[:,0].argsort()] #sorting function of Xs
		UpEnv=np.zeros((60000/step,2))
		LowEnv=np.zeros((60000/step,2))
		MeanEnv=np.zeros((60000/step,2))

		for k,value in enumerate(range(0,60000,step)):
			maxval=min(value+window/2,60000)
			minval=max(value-window/2,0)
			Slice=sortedtab[np.where(np.logical_and(sortedtab[:,0]>=minval,sortedtab[:,0]<=maxval)),1]
			sortedtab=sortedtab[np.where(sortedtab[:,0]>=minval)]
			UpEnv[k,:]=sortedtab[np.argmax(Slice),:]
			LowEnv[k,:]=sortedtab[np.argmin(Slice),:]
			MeanEnv[k,:]=np.mean(sortedtab[0:np.shape(Slice)[1]-1,:],axis=0)

		x = UpEnv[:,0]
		y = UpEnv[:,1]
		SmoothUp = interp1d(x,y, kind='linear')
		x = LowEnv[:,0]
		y = LowEnv[:,1]
		SmoothLow = interp1d(x,y, kind='linear')
		xx = np.linspace(max(LowEnv[:,0].min(),UpEnv[:,0].min()),min(LowEnv[:,0].max(),UpEnv[:,0].max()), 1000)

		NmaxPar[LeakDict[Param[0,i]],KondDict[Param[1,i]],ThickDict[Param[2,i]]]     = np.amax(data)
		NmaxPosPar[LeakDict[Param[0,i]],KondDict[Param[1,i]],ThickDict[Param[2,i]]]  = Xs[np.argmax(data)]
		MaxDrainPar[LeakDict[Param[0,i]],KondDict[Param[1,i]],ThickDict[Param[2,i]]] = np.amax(SmoothUp(xx)-SmoothLow(xx))
		MeanDrainPar[LeakDict[Param[0,i]],KondDict[Param[1,i]],ThickDict[Param[2,i]]]= np.mean(SmoothUp(xx)-SmoothLow(xx))

		plotshade=mpl.colors.colorConverter.to_rgb('r')
		plotshade=[max(0.6,x) for x in plotshade]

	#And draw
	#Open Figure
	fig = plt.figure(figsize=(9,9))
	gs = mpl.gridspec.GridSpec(4, 3)
	gs.update(left=0.05, right=0.9, top=0.97, bottom= 0.07, wspace=0.1, hspace=0.1)
	jet = cm = plt.get_cmap('jet')
	cNorm  = mpl.colors.Normalize(vmin=0, vmax=max(len(KondDict)*len(ThickDict),len(ThickDict)*len(LeakDict),len(KondDict)*len(LeakDict)))
	scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=jet)

	i=0
	ax1 = plt.subplot(gs[0,0])
	ax2 = plt.subplot(gs[1,0])
	ax3 = plt.subplot(gs[2,0])
	ax4 = plt.subplot(gs[3,0])
	ay1 = plt.subplot(gs[0,1])
	ay2 = plt.subplot(gs[1,1])
	ay3 = plt.subplot(gs[2,1])
	ay4 = plt.subplot(gs[3,1])
	az1 = plt.subplot(gs[0,2])
	az2 = plt.subplot(gs[1,2])
	az3 = plt.subplot(gs[2,2])
	az4 = plt.subplot(gs[3,2])
	for leak in LeakDict:
		for kond in KondDict:
			colorVal = scalarMap.to_rgba(i)
			i+=1
			ax1.plot(sorted(ThickDict.keys()),NmaxPar[LeakDict[leak],KondDict[kond],:],color=colorVal,marker='.',linestyle='-',markersize=15)
			ax1.yaxis.tick_right()
			ax1.yaxis.set_ticks_position('both')
			ax1.yaxis.set_ticklabels([])
			ax1.xaxis.set_ticklabels([])

			ax2.plot(sorted(ThickDict.keys()),NmaxPosPar[LeakDict[leak],KondDict[kond],:],color=colorVal,marker='.',linestyle='-',markersize=15)
			ax2.yaxis.tick_right()
			ax2.yaxis.set_ticks_position('both')
			ax2.yaxis.set_ticklabels([])
			ax2.xaxis.set_ticklabels([])

			ax3.plot(sorted(ThickDict.keys()),MaxDrainPar[LeakDict[leak],KondDict[kond],:],color=colorVal,marker='.',linestyle='-',markersize=15)
			ax3.yaxis.tick_right()
			ax3.yaxis.set_ticks_position('both')
			ax3.yaxis.set_ticklabels([])
			ax3.xaxis.set_ticklabels([])

			ax4.plot(sorted(ThickDict.keys()),MeanDrainPar[LeakDict[leak],KondDict[kond],:],color=colorVal,marker='.',linestyle='-',markersize=15)
			ax4.yaxis.tick_right()
			ax4.yaxis.set_ticks_position('both')
			ax4.yaxis.set_ticklabels([])
			ax4.set_xlabel('Thickness', fontsize=14)

	i=0
	for leak in LeakDict:
		for thick in ThickDict:
			colorVal = scalarMap.to_rgba(i)
			i+=1

			ay1.plot(sorted(KondDict.keys()),NmaxPar[LeakDict[leak],:,ThickDict[thick]],color=colorVal,marker='.',linestyle='-',markersize=15)
			ay1.yaxis.tick_right()
			ay1.yaxis.set_ticks_position('both')
			ay1.yaxis.set_ticklabels([])
			ay1.xaxis.set_ticklabels([])

			ay2.plot(sorted(KondDict.keys()),NmaxPosPar[LeakDict[leak],:,ThickDict[thick]],color=colorVal,marker='.',linestyle='-',markersize=15)
			ay2.yaxis.tick_right()
			ay2.yaxis.set_ticks_position('both')
			ay2.yaxis.set_ticklabels([])
			ay2.xaxis.set_ticklabels([])

			ay3.plot(sorted(KondDict.keys()),MaxDrainPar[LeakDict[leak],:,ThickDict[thick]],color=colorVal,marker='.',linestyle='-',markersize=15)
			ay3.yaxis.tick_right()
			ay3.yaxis.set_ticks_position('both')
			ay3.yaxis.set_ticklabels([])
			ay3.xaxis.set_ticklabels([])

			ay4.plot(sorted(KondDict.keys()),MeanDrainPar[LeakDict[leak],:,ThickDict[thick]],color=colorVal,marker='.',linestyle='-',markersize=15)
			ay4.yaxis.tick_right()
			ay4.yaxis.set_ticks_position('both')
			ay4.yaxis.set_ticklabels([])
			ay4.set_xlabel('Conductivity', fontsize=14)

	i=0
	for thick in ThickDict:
		for kond in KondDict:
			colorVal = scalarMap.to_rgba(i)
			i+=1
			az1.plot(sorted(LeakDict.keys()),NmaxPar[:,KondDict[kond],ThickDict[thick]],color=colorVal,marker='.',linestyle='-',markersize=15)
			az1.yaxis.tick_right()
			az1.yaxis.set_ticks_position('both')
			az1.yaxis.set_label_position('right')
			az1.set_ylabel('Nmax', fontsize=14)
			az1.xaxis.set_ticklabels([])

			az2.plot(sorted(LeakDict.keys()),NmaxPosPar[:,KondDict[kond],ThickDict[thick]],color=colorVal,marker='.',linestyle='-',markersize=15)
			az2.yaxis.tick_right()
			az2.yaxis.set_ticks_position('both')
			az2.yaxis.set_label_position('right')
			az2.set_ylabel('NmaxPos', fontsize=14)
			az2.xaxis.set_ticklabels([])

			az3.plot(sorted(LeakDict.keys()),MaxDrainPar[:,KondDict[kond],ThickDict[thick]],color=colorVal,marker='.',linestyle='-',markersize=15)
			az3.yaxis.tick_right()
			az3.yaxis.set_ticks_position('both')
			az3.yaxis.set_label_position('right')
			az3.set_ylabel('MaxDrain', fontsize=14)
			az3.xaxis.set_ticklabels([])

			az4.plot(sorted(LeakDict.keys()),MeanDrainPar[:,KondDict[kond],ThickDict[thick]],color=colorVal,marker='.',linestyle='-',markersize=15)
			az4.yaxis.tick_right()
			az4.yaxis.set_ticks_position('both')
			az4.yaxis.set_label_position('right')
			az4.set_ylabel('MeanDrain', fontsize=14)
			az4.set_xlabel('Leakage', fontsize=14)

	plt.show(block=False)

# }}}

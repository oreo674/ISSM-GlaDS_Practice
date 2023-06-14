import re
from netCDF4 import Dataset
from model import *
from loadmodel import *
from TreatMesh import *
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import animation
import glob as gb

# {{{ TimePlotting(RunName,PlotType):
def TimePlotting(RunName,PlotType):

	FileList=gb.glob(RunName)
	fig = plt.figure(tight_layout=True)
	ax = fig.add_subplot(111)
	jet = cm = plt.get_cmap('jet')
	cNorm  = mpl.colors.Normalize(vmin=0, vmax=len(FileList))
	scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=jet)

	for i,File in enumerate(FileList):
		print File
		LabelName='_'.join(re.split('(?<!e)[-]',File)[3:6])
		md=loadmodel(File)
		try:
			MaxStep=len(md.results.TransientSolution)
		except AttributeError:
			print ('No transient solution for Run {}, skipping'.format(File))
			continue

		NumOfFields=len(md.results.TransientSolution[0].__dict__)-3 #getting number of fields from first time minus 3 to get read of logs and solution type

		if len(md.results.TransientSolution[MaxStep-1].__dict__)!=NumOfFields:
			print ('Last time is truncated for Run {}, skipping'.format(File))
			MaxStep=MaxStep-1

		TimeData  = np.NaN*np.ones((MaxStep))
		TimeStamp = np.NaN*np.ones((MaxStep))
		Dt        = md.timestepping.time_step
		Elements  = md.mesh.elements
		X         = md.mesh.x
		Y         = md.mesh.y

		for out in range(0,MaxStep,1):

			if PlotType=='Nmean':
				Data           = 1.0e-6*md.results.TransientSolution[out].EffectivePressure
				TimeStamp[out] = md.results.TransientSolution[out].time
				TimeData[out]  = ComputeMean(Elements,X,Y,Data)

		#Now Plotting
		colorVal = scalarMap.to_rgba(i)
		ax.plot(TimeStamp, TimeData,color=colorVal,linestyle='-',linewidth=1.5,label=LabelName)
		ax.set_axes(ax)
		ax.legend(loc=1)
		ax.set_xlabel('Time [year]', fontsize=14)
		ax.set_ylabel('Mean Effective Pressure [MPa]', fontsize=14)

		#And draw
		plt.show(block=False)
	fig.savefig('../Figures/Nevol'+LabelName+'.pdf',format='pdf')
# }}}
# {{{ SaveNevolv(RunName):
def SaveNevolv(RunName):

	FileList=gb.glob(RunName)
	for i,File in enumerate(FileList):
		print File
		LabelName=re.split('/|\.nc',File)[-2]
		SaveName='../Results/Nevol'+LabelName+'.npy'
		DumpName='../Results/Nevol/Dump'+LabelName+'.npy'
		if os.path.exists(SaveName):
			print('File {} allready exist, skipping'.format(SaveName))
			continue

		md=loadmodel(File)
		try:
			MaxStep=len(md.results.TransientSolution)
		except AttributeError:
			print ('No transient solution for Run {}, skipping'.format(File))
			continue

		NumOfFields=len(md.results.TransientSolution[MaxStep-2].__dict__) #getting number of fields from first time minus 2 to get read of logs

		if len(md.results.TransientSolution[MaxStep-1].__dict__)!=NumOfFields:
			print ('Last time is truncated for Run {}, skipping'.format(File))
			MaxStep=MaxStep-1

		TimeData  = np.NaN*np.ones((MaxStep,7))
		Dt        = md.timestepping.time_step
		Elements  = md.mesh.elements
		X         = md.mesh.x
		Y         = md.mesh.y
		Bot=np.where(np.logical_and(md.mesh.x>=600,md.mesh.x<=900))
		Mid=np.where(np.logical_and(md.mesh.x>=3000,md.mesh.x<=3300))
		Top=np.where(np.logical_and(md.mesh.x>=5100,md.mesh.x<=5400))
		# Bot=np.where(np.logical_and(md.mesh.x>=10000,md.mesh.x<=15000))
		# Mid=np.where(np.logical_and(md.mesh.x>=50000,md.mesh.x<=55000))
		# Top=np.where(np.logical_and(md.mesh.x>=90000,md.mesh.x<=95000))

		for out in range(0,MaxStep,1):
			Data           = 1.0e-6*md.results.TransientSolution[out].EffectivePressure
			TimeData[out,0] = md.results.TransientSolution[out].time/md.constants.yts
			TimeData[out,1]  = ComputeMean(Elements,X,Y,Data)
			TimeData[out,2]  = np.nanmean(Data[Bot])
			TimeData[out,3]  = np.nanmean(Data[Mid])
			TimeData[out,4]  = np.nanmean(Data[Top])
			TimeData[out,5]  = np.nanmean(Data)

#		TimeData[:,6]  = np.nanmean(md.basalforcings.groundedice_melting_rate[:-1,1:],axis=0)

		np.save(SaveName,TimeData)
		TimeData=None
		Data=None
		Elements=None
		X=None
		Y=None
		md=None
# }}}
# {{{ PlotNevolv(RunName):
def PlotNevolv(RunName):

	FileList=gb.glob(RunName)
	FileList.sort()
	fig = plt.figure(tight_layout=True)
	ax = fig.add_subplot(111)
	cmap = plt.get_cmap('viridis')
	cNorm  = mpl.colors.Normalize(vmin=0, vmax=len(FileList))
	scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=cmap)

	for i,File in enumerate(FileList):
		print File
		LabelName=re.split('Nevol|.npy',File)[1]
		TimeData=np.load(File)

		#Now Plotting
		colorVal = scalarMap.to_rgba(i)
		ax.plot(TimeData[:,0], TimeData[:,1],color=colorVal,linestyle='-',linewidth=1.5,label=LabelName)
		# ax.plot(TimeData[:,0], TimeData[:,2],color='b',linestyle='-',linewidth=1.5)
		# ax.plot(TimeData[:,0], TimeData[:,3],color='c',linestyle='-',linewidth=1.5)
		# ax.plot(TimeData[:,0], TimeData[:,4],color='r',linestyle='-',linewidth=1.5)
		# ax.plot(TimeData[:,0], TimeData[:,5],color=colorVal,marker='+')
		ax.legend(loc=0)
		ax.set_xlabel('Time [year]', fontsize=14)
		ax.set_ylabel('Mean Effective Pressure [MPa]', fontsize=14)
		ax2=ax.twinx()
#		ax2.plot(TimeData[:,0], TimeData[:,6],color='k',linestyle='-')
		#And draw
		plt.show(block=False)
	fig.savefig('../Figures/Nevol'+LabelName+'.pdf',format='pdf')
# }}}
# {{{ SeasonalPlot(DataName)

def SeasonalPlot(DataName):
	home = os.getenv("HOME")
	names = ['mwer','ogag','bdef']#Add names when they become available
	plotcolor=['r','g','b']#Add a new plot style when you add a new name
	plotmark='.'
	prefix=home+'/Dropbox/GHIP_Results/'

	expe=raw_input('Give the tag of the experiment you want to plot: ')

	if len(expe)==1:
		plotexpe=[expe+'1',expe+'2',expe+'3',expe+'4',expe+'5',expe+'6' ]
	else:
		plotexpe=[expe]

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_axes(ax)
	ax.set_xlabel('Time [s]', fontsize=14)
	ax.set_ylabel('Mean Effective Pressure [MPa]', fontsize=14)
	ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	ax.set_title('Comparison for run '+ expe, fontsize=18)
	#Make space for the legend
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width*0.85 , box.height])

	for i,name in enumerate(names):
		for j,tag in enumerate(plotexpe):
			NCFile=prefix+name+'/box100by20/seasonal/'+tag+'.nc'
			try:
				DatFile	 = Dataset(NCFile, mode='r')
			except RuntimeError:
				print 'File '+ str(NCFile)+' does not exist'
				continue
			Time     = DatFile.variables['time'][:]
			data     = DatFile.variables[DataName][:,:]
			Elts     = DatFile.variables['cellconnect'][:,:]
			x        = DatFile.variables['xy'][0,:]
			y        = DatFile.variables['xy'][1,:]
			if DataName=='N':
				data=1.0e-6*data
			DatFile.close()

			DatMean=np.NaN*np.ones(np.shape(Time))
			for timestep,t in enumerate(Time):
				DatMean[timestep]=ComputeMean(Elts,x,y,data[timestep,:])


			plotshade=mpl.colors.colorConverter.to_rgb(plotcolor[i])
			plotshade=[max(j*0.18,x) for x in plotshade]
			ax.plot(Time, DatMean, linestyle='',color=plotshade,marker=plotmark,label=name+'_'+tag)
			#and place the legend outside on right
			ax.legend(loc='center left',bbox_to_anchor=(1, 0.5))

		plt.show()
#	fig.savefig('../Figures/Seasonal-'+DataName+'-'.join(plotexpe)+'.pdf',format='pdf')

# }}}
# {{{ ElevationMeans(RunName,DataName)
def ElevationMeans(RunName,DataName):
	DatFile	 = Dataset(RunName, mode='r')
	H        = DatFile.varaiable['H'][:]
	Xs			 = DatFile.variables['xy'][0,:]
	Ys			 = DatFile.variables['xy'][1,:]
	elts		 = DatFile.variables['cellconnect'][:].T
	data     = DatFile.variables[DataName][:,:]
	time     = DatFile.variables['time'][:]

	if DataName=='N':
		data=1.0e-6*data
	DatFile.close()

	ElevationLims=[np.min(H),np.max(H)]
	SliceThicness=100.
	SliceNumber=int(floor((ElevationLims[1]-ElevationLims[0])/SliceThickness))

	SliceMeans=np.Nan*np.ones((len(time),SliceNumber))

	for slices in range(0,SliceNumber):
		SliceData=np.zeros(np.shape(data))
		SliceNodesIndex=np.where((H>=slices*SliceThickness) & (H<(slices+1)*SliceThickness))
		SliceX=Xs[SliceNodesIndex]
		SliceY=Ys[SliceNodesIndex]
		SliceData=data[:,SliceNodesIndex]
		# for node in SliceNodesIndex:
		# 	try:
		# 		SliceCells=np.append(SliceCells,where(elts==node)[1])
		# 	except NameError:
		# 		SliceCells=where(elts==node)[1]

# }}}
# {{{ TimeAnimation(RunName,DataName)

def TimeAnimation(RunName,DataName):

	DatFile	 = Dataset(RunName, mode='r')
	Xs			 = DatFile.variables['xy'][0,:]
	Ys			 = DatFile.variables['xy'][1,:]
	elts		 = DatFile.variables['cellconnect'][:].T
	data     = DatFile.variables[DataName][:,:]
	time     = DatFile.variables['time'][:]

	if DataName=='N':
		data=1.0e-6*data
	DatFile.close()

	#Now Plotting
	fig = plt.figure(tight_layout=False)
	ax = fig.add_subplot(111)
	ax2 = fig.add_axes([0.1, 0.1, 0.85, 0.05])
	ax.set_axes(ax)
	ax.set_xlim([0.0,1.0e5])
	ax.set_ylim([0.0,2.0e4])
	ax.set_xlabel('X axis [m]', fontsize=14)
	ax.set_ylabel('Y axis [m]', fontsize=14)
	ax.set_title(DataName+' for the run ', fontsize=18)

	# gathering data limits
	lims=[data.min(),data.max()]

	norm = mpl.colors.Normalize(vmin=lims[0], vmax=lims[1])
	cmap = plt.cm.RdBu
	cbar = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='horizontal')
	# animation function
	def animate(i):
		print i
		z = data[i,:]
		cont = 	ax.tricontourf(Xs, Ys, elts, z, 128, cmap=cmap)
		return cont

	anim = animation.FuncAnimation(fig, animate, frames=len(time))
	anim.save('animation.mp4')
	plt.show()
# }}}

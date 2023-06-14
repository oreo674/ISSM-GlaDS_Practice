#!/usr/bin/env python2
"""Plotter for the GHIP exercise, the plotter relies on the conventions as given on the comparison framework,
Available plot as of now are:
  -Sections (called with ./GHIPPloter.py Sections)
  -Map Plot (called with ./GHIPPloter.py Maps)
  -Eveloppes Plot (called with ./GHIPPloter.py Enveloppes)
  -Evolution of the mean altitudinal vallues (called with ./GHIPPloter.py Evol)
---Sections---
	Description :
	Plot a variable function of the x coordinate for all the participating model and a given experiment.

	The experiment tag is asked later , if only the letter of the experiment is given then all the subexperiments are plotted

---Enveloppes---
	Description :

	Needed parameters :

	The experiment tag is asked later , if only the letter of the experiment is given then all the subexperiments are plotted

	Example of use :
	GHIPPlotter.AllEnvelopes('N')

"""
from netCDF4 import Dataset
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from matplotlib.ticker import SymmetricalLogLocator
from TreatMesh import *
# {{{ AllSections()
def AllSection():
	plotcolor='r'
	plotmark='.'

	#Open Figure
	fig = plt.figure()
	ax = fig.add_subplot(111)
	#	ax.set_axes(ax)
	ax.set_xlabel('Distance from front $[km]$', fontsize=14)
	ax.set_ylabel(TitleList[DataName][0]+' $['+TitleList[DataName][1]+']$', fontsize=14)
	ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	if len(expe)==1:
		ax.set_title('Comparison for run '+ expe, fontsize=18)
	else:
		ax.set_title('Section for run '+ expe, fontsize=18)
	#Make space for the legend
	box = ax.get_position()
	ax.set_position([box.x0-0.02, box.y0, box.width*0.85 , box.height])

	for j,tag in enumerate(plotexpe):
		NCFile=resdir+tag+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			data, coords=GetVar(DatFile,DataName,datacoord,False,tag.endswith('mh1'))
			if datacoord=='coords_ch':
				try:
					coords = DatFile.variables[datacoord][:,:]
				except KeyError:
					connect=DatFile.variables['connect_ch'][:,:]
					coords=ComputeCenter(connect,coords[0,:],coords[1,:])
			DatFile.close()
		except (IOError, RuntimeError):
			print ('File {} does not exist'.format(NCFile))
			continue
		if expe in['A','B','C','D']:
			data=data[np.where(np.logical_and(coords[0,:]>=0,coords[0,:]<=1.0e5))]
			coords=coords[0,np.where(np.logical_and(coords[0,:]>=0,coords[0,:]<=1.0e5))]
		elif expe in['E','F']:
			data=data[np.where(np.logical_and(coords[0,:]>=0,coords[0,:]<=6.0e3))]
			coords=coords[0,np.where(np.logical_and(coords[0,:]>=0,coords[0,:]<=6.0e3))]

		plotshade=mpl.colors.colorConverter.to_rgb(plotcolor)
		plotshade=[max(j*0.18,x) for x in plotshade]

		ax.plot(1.0e-3*coords[0,:], data, linestyle='',color=plotshade,marker=plotmark,label=tag)
		#and place the legend outside on right
		ax.legend(loc='center left',bbox_to_anchor=(1, 0.5))
	plt.show
	return fig

# }}}
# {{{ SynthMap()

def SynthMap():
	axar={}
	if expe in['A','B','C','D']:
		suprows=(max(0,len(plotexpe)-1))/3
		cols=min(3,len(plotexpe))
		ratio=5
	elif expe in['E','F']:
		suprows=(max(0,len(plotexpe)-1))/3
		cols=min(3,len(plotexpe))
		ratio=5.55
	elif len(expe)==2: #plotting only one run
		suprows=0
		cols=1
		if expe.startswith(('A','B','C','D')):
			ratio=5
		if expe.startswith(('E','F')):
			ratio=5.5
	fig, ax = plt.subplots(1+suprows,cols,sharex='col',sharey='row',figsize=(2.5+ratio*cols,1.5+(suprows+1)))

	if len(plotexpe)==1:
		axar['0']=ax
	else:
		for j,val in enumerate(plotexpe):
			axar[str(j)]=ax[j/3,j%3]
	figsizing=fig.get_size_inches()
	vblank=(figsizing[1]-(suprows+1))/((suprows+2)*figsizing[1]) #vert blank space percentage
	hblank=(figsizing[0]-cols*ratio)/((cols/2+2)*figsizing[0]) #horiz blank space percentage
	vtop=vblank
	vbot=vblank
	vspace=vblank
	vplotspace=1.-(vtop+vbot+vspace*suprows)
	height=vplotspace/(1+suprows)
	hleft=hblank*1.0
	hright=hblank*1.3
	if cols>1:
		hspace=hblank/(cols-1)
	else:
		hspace=0
	hplotspace=1.-(hleft+hright+(cols-1.)*hspace)
	width=hplotspace/cols


	# gathering data limits
	tmplims=np.nan*np.ones((len(plotexpe),2))
	for j, tag in enumerate(plotexpe):
		NCFile=resdir+tag+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			data, coords=GetVar(DatFile,DataName,datacoord,False,tag.endswith('mh1'))
			DatFile.close()
			tmplims[j,0]=np.nanmin(data)
			tmplims[j,1]=np.nanmax(data)
		except (RuntimeError,IOError):
			continue

	#defining colormap
	lims=[np.nanmin(tmplims[:,0]),np.nanmax(tmplims[:,1])]
	norm=mpl.colors.SymLogNorm(linthresh=0.2*np.diff(lims), linscale=2,vmin=lims[0], vmax=lims[1])
	#norm=MidpointNormalize(midpoint=2.)
	cmap = plt.cm.viridis

	#now for plotting stuff
	for j, tag in enumerate(plotexpe):
		NCFile=resdir+tag+'.nc'
		#Setting up position title and labels
		axar[str(j)].set_position([hleft+(hspace+width)*(j%3),1-(vtop+height)-(vspace+height)*(j/3),width,height])
		axar[str(j)].set_title('Simulation '+tag, fontsize=18)
		if (j/3)==suprows:
			axar[str(j)].set_xlabel('X axis $[km]$', fontsize=14)
		if (j%3)==0:
			axar[str(j)].set_ylabel('Y axis $[km]$', fontsize=14)

		try:
			#get the data if the file exist
			DatFile	 = Dataset(NCFile, mode='r')
			data, coords=GetVar(DatFile,DataName,datacoord,False,tag.endswith('mh1'))
			if datacoord=='coords_ch':
				try:
					elts = DatFile.variables['connect_ch'][:,:].T.astype(int)
					x    = np.asarray([coords[0,elts[:,0]]*1.0e-3,coords[0,elts[:,1]]*1.0e-3])
					y    = np.asarray([coords[1,elts[:,0]]*1.0e-3,coords[1,elts[:,1]]*1.0e-3])
					scaterplot=False
				except KeyError:
					scaterplot=True
			DatFile.close()
		except (RuntimeError,IOError):
			#Deal with non existant files
			print ('File {} does not exist'.format(NCFile))
			if expe in ['A','B','C','D']:
				axar[str(j)].set_xlim([0,100])
				axar[str(j)].set_ylim([0,20])
			elif expe in ['E','F']:
				axar[str(j)].set_xlim([0,6])
				axar[str(j)].set_ylim([-0.5,0.5])
			continue

		if expe in ['A','B','C','D']:
			Xlims=[1.0e-3*np.min(coords[0,:]),min(100.,1.0e-3*np.max(coords[0,:]))]
		elif expe in ['E','F']:
			Xlims=[1.0e-3*np.min(coords[0,:]),min(6.0,1.0e-3*np.max(coords[0,:]))]
		if np.all(coords[1,:]==10*1.0e3):
			Ylims=[0,20]
		elif np.all(coords[1,:]==0.):
			Ylims=[-0.5,0.5]
		else:
			Ylims=[1.0e-3*np.min(coords[1,:]),1.0e-3*np.max(coords[1,:])]
		axar[str(j)].set_xlim(Xlims)
		axar[str(j)].set_ylim(Ylims)
		#channel plots are different
		if datacoord=='coords_ch':
			dataspan=max(np.diff(lims)[0],100)
			datastep=[lims[0]+0.01*dataspan,lims[0]+0.1*dataspan,lims[0]+0.3*dataspan,lims[0]+0.8*dataspan,lims[1]]
#			datastep=dataspan/5.
			print lims
			for step in range(1,5):
				indices=np.where(np.logical_and(data>datastep[step-1],data<=datastep[step]))[0]
				if len(indices)>0:
					if cols==1:
						legtext='From {:.{prec}f} to {:.{prec}f}'.format(datastep[step-1],datastep[step],prec=0)
					else:
						legtext='From {:.{prec}f} \nto {:.{prec}f}'.format(datastep[step-1],datastep[step],prec=0)
					if scaterplot:
						axar[str(j)].scatter(coords[0,indices[0]]*1.0e3,coords[1,indices[0]]*1.0e3,
																 color='red',s=5.0*step,label=legtext)
					else:
						axar[str(j)].plot(x[:,indices[0]],y[:,indices[0]],
															color='blue',linewidth=step,solid_capstyle='round',
															label=legtext)
					l=axar[str(j)].legend(bbox_to_anchor=(1-hright,0.5*vbot,hright,1-vtop),
																bbox_transform=fig.transFigure,frameon=False,
																mode='expand',fontsize='small',ncol=1,
																title=TitleList[DataName][0]+' $['+TitleList[DataName][1]+']$')
					l.get_title().set_wrap(True)
					l.get_title().set_ha('left')
					l.get_title().set_position((50,0))
				if len(indices)>1:
					if scaterplot:
						axar[str(j)].scatter(coords[0,indices[1:]]*1.0e-3,coords[1,indices[1:]]*1.0e-3,
																 color='red',s=5.0*step)
					else:
						axar[str(j)].plot(x[:,indices[1:]],y[:,indices[1:]],
															color='blue',linewidth=step,solid_capstyle='round')
		#Now dealing with maps
		else:
			if np.all(coords[1,:]==0.) or np.all(coords[1,:]==10*1.0e3):
				axar[str(j)].imshow([data], extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
														origin='lower',norm=norm,cmap=cmap)
			else:
				grid_x, grid_y = np.mgrid[Xlims[0]:Xlims[1]:101j, Ylims[0]:Ylims[1]:101j]
				gridded=griddata(1.0e-3*coords.T, data, (grid_x, grid_y), method='linear')
				axar[str(j)].imshow(gridded.T, extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
														origin='lower',norm=norm,cmap=cmap)

			if cols==1:
				wrapping=True
				horal='center'
				veral='bottom'
				pos=(-0.1,0.5)
			else:
				wrapping=False
				horal='right'
				veral='center'
				pos=(0.0,0.5)
			ax2 = fig.add_axes([0.98-0.02*(4-cols),0.05,0.01*(4-cols), 0.9])
			t=ax2.set_title(TitleList[DataName][0]+' $['+TitleList[DataName][1]+']$',
											fontsize=18,rotation=90,
											multialignment='center',horizontalalignment=horal,verticalalignment=veral,
											wrap=wrapping,position=pos)
			cbar = mpl.colorbar.ColorbarBase(ax2, norm=norm,cmap=cmap, orientation='vertical')
			cbar.solids.set_edgecolor("face")

	#And draw
	plt.show
	return fig

# }}}
# {{{ AllEnvelopes()

def AllEnvelopes():
	plotcolor='r'
	plotmark='.'

	#Open Figure
	fig ,ay= plt.subplots(12,figsize=(16,12))
	gs = mpl.gridspec.GridSpec(4, 3)
	gs.update(left=0.07, right=0.93, top=0.95, bottom= 0.05, wspace=0.07, hspace=0.1)

	ax1 = plt.subplot(gs[0:-1,0:2])
	ax2 = plt.subplot(gs[-1,0:2])
	ay[0] = plt.subplot(gs[0,-1])
	ay[1] = plt.subplot(gs[1,-1])
	ay[2] = plt.subplot(gs[2,-1])
	ay[3] = plt.subplot(gs[3,-1])

	SynthDat  = np.NaN*np.ones((4,len(plotexpe)))

	for j,tag in enumerate(plotexpe):
		NCFile=resdir+tag+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			data, coords=GetVar(DatFile,DataName,datacoord,False,tag.endswith('mh1'))
			if datacoord=='coords_ch':
				try:
					coords = DatFile.variables[datacoord][:,:]
				except KeyError:
					connect=DatFile.variables['connect_ch'][:,:]
					coords=ComputeCenter(connect,coords[0,:],coords[1,:])
				if (np.min(coords[0,:])-np.max(coords[0,:]))==0:
					connect=DatFile.variables['connect_ch'][:,:]
					data, coords=GetVar(DatFile,DataName,datacoord,False,tag.endswith('mh1'))
					coords=ComputeCenter(connect,coords[0,:],coords[1,:])
			DatFile.close()
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))
			continue

		data, coords=GetVar(DatFile,DataName,datacoord,False,tag.endswith('mh1'))
		if datacoord=='coords_ch':
			try:
				coords = DatFile.variables[datacoord][:,:]
			except KeyError:
				connect=DatFile.variables['connect_ch'][:,:]
				coords=ComputeCenter(connect,coords[0,:],coords[1,:])
			if (np.min(coords[0,:])-np.max(coords[0,:]))==0:
				connect=DatFile.variables['connect_ch'][:,:]
				data, coords=GetVar(DatFile,DataName,datacoord,False,tag.endswith('mh1'))
				coords=ComputeCenter(connect,coords[0,:],coords[1,:])
		DatFile.close()

		Xlims=[np.min(coords[0,:]),np.max(coords[0,:])]
		Range=int(Xlims[1]-Xlims[0])
		step =Range/1000
		window=Range/50
		valtab=np.vstack((coords[0,:],data)).T
		sortedtab=valtab[valtab[:,0].argsort()] #sorting function of Xs
		UpEnv=np.zeros((len(np.arange(Xlims[0],Xlims[1],step)),2))
		LowEnv=np.zeros((len(np.arange(Xlims[0],Xlims[1],step)),2))
		MeanEnv=np.zeros((len(np.arange(Xlims[0],Xlims[1],step)),2))
		for k,value in enumerate(np.arange(Xlims[0],Xlims[1],step)):
			maxval=min(value+window/2,Xlims[1])
			minval=max(value-window/2,Xlims[0])
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
		xx = np.linspace(max(np.nanmin(LowEnv[:,0]),np.nanmin(UpEnv[:,0])),min(np.nanmax(LowEnv[:,0]),np.nanmax(UpEnv[:,0])), 1000)
		# compute globals
		SynthDat[0,j] = np.nanmax(data)		   #Nmax
		SynthDat[1,j] = coords[0,np.argmax(data)]*1.0e-3 #NmaxPos
		SynthDat[2,j] = np.nanmax(SmoothUp(xx)-SmoothLow(xx)) #MaxDrain
		SynthDat[3,j] = np.nanmean(SmoothUp(xx)-SmoothLow(xx)) #MeanDrain

		plotshade=mpl.colors.colorConverter.to_rgb(plotcolor)
		plotshade=[max(0.6,x) for x in plotshade]

		if j==3 or len(plotexpe)==1:
			ax1.set_ylabel('{} $[{}]$'.format(TitleList[DataName][0],TitleList[DataName][1]), fontsize=14)
			ax1.set_title('Overview for run '+ tag, fontsize=18)
			ax1.axis([0, 1.0e-3*Xlims[1], np.nanmin(data), SynthDat[0,j]])
			ax1.plot(1.0e-3*coords[0,:], data,color=plotshade,marker=plotmark,linestyle='')
			ax1.plot(1.0e-3*xx, SmoothLow(xx),color=plotcolor,linewidth=1.5,label=name)
			ax1.plot(1.0e-3*xx, SmoothUp(xx),color=plotcolor,linewidth=1.5)
			ax1.xaxis.set_ticklabels([])

			ax2.axis([0, 1.0e-3*Xlims[1], np.nanmin(SmoothUp(xx)-SmoothLow(xx)), SynthDat[2,j]])
			ax2.set_xlabel('Distance from front [km]', fontsize=14)
			ax2.set_ylabel('{} \nSpread $[{}]$'.format(TitleList[DataName][0],TitleList[DataName][1]), fontsize=14)
			ax2.plot(1.0e-3*xx, SmoothUp(xx)-SmoothLow(xx),color=plotcolor,linewidth=1.5,linestyle='-')

	#And draw
	indexes=range(1,len(plotexpe)+1,1)
	for p in range(0,4):
		ay[p].plot(indexes,SynthDat[p,:],color=plotcolor,marker=plotmark,linestyle='',markersize=15)
		ay[p].plot(indexes,SynthDat[p,:],linestyle='-',color=plotshade)
		ay[p].yaxis.tick_right()
		ay[p].yaxis.set_ticks_position('both')
		ay[p].yaxis.set_label_position('right')
		if tag.startswith('A'):
			ay[p].axis([1, 6, 1.1*np.nanmin(SynthDat[p,:])-0.1*np.nanmax(SynthDat[p,:]), 1.1*np.nanmax(SynthDat[p,:])-0.1*np.nanmin(SynthDat[p,:])])
			ay[p].set_xticks([1,2,3,4,5,6])
			ay[3].xaxis.set_ticklabels(['A1','A2','A3','A4','A5','A6'])
		if tag.startswith('B'):
			ay[p].axis([1, 5, 1.1*np.nanmin(SynthDat[p,:])-0.1*np.nanmax(SynthDat[p,:]), 1.1*np.nanmax(SynthDat[p,:])-0.1*np.nanmin(SynthDat[p,:])])
			ay[p].set_xticks([1,2,3,4,5])
			ay[3].xaxis.set_ticklabels(['B1','B2','B3','B4','B5'])
		if tag.startswith('C'):
			ay[p].axis([1, 4, 1.1*np.nanmin(SynthDat[p,:])-0.1*np.nanmax(SynthDat[p,:]), 1.1*np.nanmax(SynthDat[p,:])-0.1*np.nanmin(SynthDat[p,:])])
			ay[p].set_xticks([1,2,3,4])
			ay[3].xaxis.set_ticklabels(['C1','C2','C3','C4'])
		if tag.startswith('D'):
			ay[p].axis([1, 5, 1.1*np.nanmin(SynthDat[p,:])-0.1*np.nanmax(SynthDat[p,:]), 1.1*np.nanmax(SynthDat[p,:])-0.1*np.nanmin(SynthDat[p,:])])
			ay[p].set_xticks([1,2,3,4,5])
			ay[3].xaxis.set_ticklabels(['D1','D2','D3','D4','D5'])
		if tag.startswith('E'):
			ay[p].axis([1, 5, 1.1*np.nanmin(SynthDat[p,:])-0.1*np.nanmax(SynthDat[p,:]), 1.1*np.nanmax(SynthDat[p,:])-0.1*np.nanmin(SynthDat[p,:])])
			ay[p].set_xticks([1,2,3,4,5])
			ay[3].xaxis.set_ticklabels(['E1','E2','E3','E4','E5'])
		if tag.startswith('F'):
			ay[p].axis([1, 5, 1.1*np.nanmin(SynthDat[p,:])-0.1*np.nanmax(SynthDat[p,:]), 1.1*np.nanmax(SynthDat[p,:])-0.1*np.nanmin(SynthDat[p,:])])
			ay[p].set_xticks([1,2,3,4,5])
			ay[3].xaxis.set_ticklabels(['F1','F2','F3','F4','F5'])
	ay[0].set_ylabel('Maximum \n{} $[{}]$'.format(TitleList[DataName][0],TitleList[DataName][1]),
									 multialignment='center',position=(-0.1,0.5),
									 fontsize=14,wrap=True)
	ay[1].set_ylabel('Position of the maximum\n {} $[km]$'.format(TitleList[DataName][0],),
									 multialignment='center',position=(-0.1,0.5),
									 fontsize=14,wrap=True)
	ay[2].set_ylabel('Maximum spread for \n {}'.format(TitleList[DataName][0],),
									 multialignment='center',position=(-0.1,0.5),
									 fontsize=14,wrap=True)
	ay[3].set_ylabel('Mean spread for \n {}'.format(TitleList[DataName][0],),
									 multialignment='center',position=(-0.1,0.5),
									 fontsize=14,wrap=True)
	plt.show
	return fig

# }}}
# {{{ def MeanEvol():

def MeansEvol():
	axar={}
	cols=len(plotexpe)
	ratio=2
	fig, ax = plt.subplots(2,cols,sharex='col',sharey='row',figsize=(5+ratio*cols,6))
	if len(plotexpe)==1:
		axar['0']=[ax[0],ax[1]]
	else:
		for j,val in enumerate(plotexpe):
			axar[str(j)]=[ax[0,j],ax[1,j]]
	figsizing=fig.get_size_inches()
	vblank=(figsizing[1]-4.5)/(3*figsizing[1]) #vert blank space percentage
	vtop=vblank*1.
	vbot=vblank*1.8
	vspace=vblank*0.2
	vplotspace=1.-(vtop+vbot+vspace)
	height=vplotspace/2.
	hblank=(figsizing[0]-cols*ratio)/((cols/2+2)*figsizing[0]) #horiz blank space percentage
	hleft=hblank*0.9
	hright=hblank*1.1
	if cols>1:
		hspace=hblank/(cols-1)
	else:
		hspace=0
	hplotspace=1.-(hleft+hright+(cols-1.)*hspace)
	width=hplotspace/cols

	# gathering data limits
	tmplims=np.nan*np.ones((len(plotexpe),2))
	for j, tag in enumerate(plotexpe):
		NCFile=resdir+tag+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			data, coords=GetVar(DatFile,DataName,datacoord,True,tag.endswith('mh1'))
			DatFile.close()
			tmplims[j,0]=np.nanmin(data)
			tmplims[j,1]=np.nanmax(data)
		except (RuntimeError,IOError):
			continue

	lims=[np.nanmin(tmplims[:,0]),np.nanmax(tmplims[:,1])]
	MaxTime=0.0
	MinTime=1.0e3

	for j,tag in enumerate(plotexpe):
		NCFile=resdir+tag+'.nc'
		#define the plots
		axar[str(j)][0].set_position([hleft+(hspace+width)*j,1-(vtop+0.5*height),width,0.5*height])
		axar[str(j)][1].set_position([hleft+(hspace+width)*j,vbot,width,1.5*height])
		axar[str(j)][0].set_title('Simu. '+tag, fontsize=18)
		axar[str(j)][1].set_xlabel('Surface Elevation $[m]$', fontsize=14)
		plt.setp(axar[str(j)][1].xaxis.get_majorticklabels(), rotation=45)
		if j==0:
			axar[str(j)][0].set_ylabel(TitleList[DataName][0]+' $['+TitleList[DataName][1]+']$', fontsize=14,wrap='True')
			if tag.startswith('C'):
				axar[str(j)][1].set_ylabel('Time $[hours]$', fontsize=14)
			else:
				axar[str(j)][1].set_ylabel('Time $[years]$', fontsize=14)

		try:
			#get the data if the file exist
			DatFile	 = Dataset(NCFile, mode='r')
			time     = DatFile.variables['time'][:]
			bed      = DatFile.variables['B'][:]
			data, coords=GetVar(DatFile,DataName,datacoord,True,tag.endswith('mh1'))
			thick, nodes=GetVar(DatFile,'H','coords1',False,tag.endswith('mh1'))
			surf=thick+bed
			if datacoord=='coords_ch':
				try:
					if len(DatFile.dimensions['dim'])==1:
						Xs		 = DatFile.variables[datacoord][:]
						if np.max(Xs>=100*1.0e3):
							Ys		 = 10*1.0e3*np.ones(np.shape(Xs))
						else:
							Ys		 = np.zeros(np.shape(Xs))
						coords=np.vstack((Xs,Ys))
					else:
						coords = DatFile.variables[datacoord][:,:]
				except KeyError:
					connect=DatFile.variables['connect_ch'][:,:]
					coords=ComputeCenter(connect,coords[0,:],coords[1,:])
				if (np.min(coords[0,:])-np.max(coords[0,:]))==0:
					connect=DatFile.variables['connect_ch'][:,:]
					data, coords=GetVar(DatFile,DataName,datacoord,True,tag.endswith('mh1'))
					coords=ComputeCenter(connect,coords[0,:],coords[1,:])
			DatFile.close()
		except (RuntimeError,IOError):
			#Deal with non existant files
			print ('File {} does not exist'.format(NCFile))
			continue

		#define time unit function of run
		if tag.startswith('C'):
			data=data[-25:,:]
			time=time[-25:]
			time=(time-min(time))/(3600.)
		else:
			# data=data[-365:,:]
			# time=time[-365:]
			time=(time-min(time))/(365.*24.*3600.)

		#define altitudinal stepping on the first file only
		if j==0:
			Xlims=[1.0e-3*np.min(coords[0,:]),1.0e-3*np.max(coords[0,:])]
			Ylims=[0,20]#1.0e-3*np.min(coords[1,:]),1.0e-3*np.max(coords[1,:])]
			grid_x, grid_y = np.mgrid[Xlims[0]:Xlims[1]:101j, Ylims[0]:Ylims[1]:101j]
			griddedsurf=griddata(1.0e-3*nodes.T, surf, (grid_x, grid_y), method='nearest')
			AltiBounds=np.asarray([np.nanmin(griddedsurf[:,0]),np.nanmax(griddedsurf[:,0])])
			StepNum=20.
			StepSize=np.diff(AltiBounds)/StepNum
			HalfStep= StepSize/2.
			Steps=np.arange(AltiBounds[0]+HalfStep,AltiBounds[1],StepSize)
			StepLims=np.zeros(np.size(Steps)+1)
			StepMeans=np.zeros((np.size(time),np.size(Steps)))
			StepLims[0]=AltiBounds[0]
			plot_tab=np.nan*np.ones((len(time),len(Steps),len(plotexpe)))

		#grid weights to take the sides into account, just for rectangle and given steps for now
		if tag.startswith(('E','F')):
			Weights=np.ones(np.shape(grid_x))
		else:
			Weights_x=np.zeros(np.shape(grid_x))
			Weights_x[np.where(grid_x+0.5<(np.max(grid_x)))]=0.5
			Weights_x[np.where(grid_x-0.5>(np.min(grid_x)))]=Weights_x[np.where(grid_x-0.5>(np.min(grid_x)))]+0.5
			Weights_y=np.zeros(np.shape(grid_y))
			Weights_y[np.where(grid_y+0.5<(np.max(grid_y)))]=0.5
			Weights_y[np.where(grid_y-0.5>(np.min(grid_y)))]=Weights_y[np.where(grid_y-0.5>(np.min(grid_y)))]+0.5
			Weights=Weights_x*Weights_y

		for t,date in enumerate(time):
			griddeddata=griddata(1.0e-3*coords.T, data[t,:], (grid_x, grid_y), method='nearest')
			for s,step in enumerate(Steps):
				StepLims[s+1]=step+HalfStep
				Onslice=np.where(np.logical_and(griddedsurf[:,0]>step-HalfStep,griddedsurf[:,0]<step+HalfStep))
				Onbound=np.where(np.logical_or(griddedsurf[:,0]==step-HalfStep,griddedsurf[:,0]==step+HalfStep))
				StepSum=np.sum(griddeddata[Onslice,:]*Weights[Onslice,:])+np.sum(griddeddata[Onbound,:]*0.5*Weights[Onbound,:])
				StepMeans[t,s]=StepSum/(np.sum(Weights[Onslice,:])+np.sum(0.5*Weights[Onbound,:]))

		MaxTime=max(MaxTime,np.max(time))
		MinTime=min(MinTime,np.min(time))
		plot_tab[:,:,j]=StepMeans[:,:]-np.mean(StepMeans,axis=0)

		#plotting data function of surf elev at last time
		axar[str(j)][0].plot(griddedsurf,griddeddata,linestyle='',color='r',marker='.')
		axar[str(j)][0].axis([StepLims[0], StepLims[-1],lims[0],lims[1]])


	#defining colormap
	lims=[np.nanmin(plot_tab),np.nanmax(plot_tab)]
	norm=mpl.colors.SymLogNorm(linthresh=0.2*np.diff(lims), linscale=2,vmin=lims[0], vmax=lims[1])
	#	tickpos=[lims[0],-0.2*np.diff(lims),-0.1*np.diff(lims),0.0,0.1*np.diff(lims),0.2*np.diff(lims),lims[1]]
	#norm=MidpointNormalize(midpoint=2.)
	cmap = plt.cm.viridis

	#and plot evol
	for j,tag in enumerate(plotexpe):
		axar[str(j)][1].imshow(plot_tab[:,:,j], cmap=cmap, norm=norm,aspect='auto',
													 extent=(StepLims[0], StepLims[-1],MinTime,MaxTime))
		axar[str(j)][1].axis([StepLims[0], StepLims[-1],MinTime,MaxTime])

	#and colorbar last
	ax2 = fig.add_axes([0.98-0.02*(7-cols),0.05,0.01*(7-cols), 0.9])
	ax2.set_title('EBM '+TitleList[DataName][0]+' $['+TitleList[DataName][1]+']$',
								fontsize=18,rotation=90,
								multialignment='center',horizontalalignment='right',verticalalignment='center',
								position=(0.0,0.5))
	cbar = mpl.colorbar.ColorbarBase(ax2, norm=norm,cmap=cmap,orientation='vertical')
	cbar.solids.set_edgecolor("face")

	plt.show
	return fig

# }}}
# {{{ def GetVar(DatFile,AllTime,cism_case)

def GetVar(DatFile,name,coords,AllTime,cism_case):
	#getting data array
	if np.size(DatFile.variables[name].dimensions)==1:
		data = DatFile.variables[name][:]
	else:
		if AllTime:
			data = DatFile.variables[name][:,:]
		else:
			lasttime = np.shape(DatFile.variables[name])[0]-1
			data     = DatFile.variables[name][lasttime,:]
	if name=='N':
		data=1.0e-6*data
	#Hack for some remnant mask that get to inf
	data[np.where(data>1.0e30)]=0.0
	#and getting the corresponding coordinates
	if len(DatFile.dimensions['dim'])==1:
		if coords=='coords_ch':
			try:
				connect= DatFile.variables['coords_ch']
				Xs		 = DatFile.variables['coords1'][:]
			except KeyError:
				Xs		 = DatFile.variables[coords][:]
		else:
			Xs		 = DatFile.variables[coords][:]
		if np.max(Xs>=100*1.0e3):
			Ys		 = 10*1.0e3*np.ones(np.shape(Xs))
		else:
			Ys		 = np.zeros(np.shape(Xs))
	elif coords=='coords_ch':
		try:
			connect= DatFile.variables['coords_ch']
			if cism_case:
				Xs		 = DatFile.variables['coords3'][0,:]
				Ys		 = DatFile.variables['coords3'][1,:]
			else:
				Xs		 = DatFile.variables['coords1'][0,:]
				Ys		 = DatFile.variables['coords1'][1,:]
		except KeyError:
			print('No channel connectivity exists (or badly named) we will work with channel midpoints')
			Xs		 = DatFile.variables[coords][0,:]
			Ys		 = DatFile.variables[coords][1,:]
	else:
		Xs		 = DatFile.variables[coords][0,:]
		Ys		 = DatFile.variables[coords][1,:]

	coords = np.vstack((Xs,Ys))
	return data, coords

# }}}

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

if __name__=="__main__":
	PlotList={'Sections': AllSection,
						'Maps':SynthMap,
						'Enveloppes':AllEnvelopes,
						'Evol':MeansEvol}

	TitleList={'h':['distributed water layer thickness',u'm'],
						 'hstore':['stored water: effective layer thickness',u'm'],
						 'S':['channel cross-sectional area',u'm^2'],
						 'Q':['channel discharge',u'm^3/s'],
						 'q':['distributed water discharge (abs)',u'm^2/s'],
						 'H':['ice thickness',u'm'],
						 'B':['bedrock elevation',u'm'],
						 'N':['effective pressure',u'MPa'],
						 'Ee':['EPL thickness',u'm']}

	if len(sys.argv)==1: # print help if no filename plot name is given:
		print __doc__
		sys.exit(0)

	resdir=os.getcwd()+'/'

	figdir=resdir+'../figures'

	name = raw_input('Give the name of the experimenter as it appears in the file name: ')
	DataName=raw_input('Give the variable you want to plot: ')
	expe=raw_input('Give the tag of the experiment you want to plot: ')


	if len(expe)==1:
		if expe in ['A']:
			plotexpe=[expe+'1_'+name,expe+'2_'+name,expe+'3_'+name,expe+'4_'+name,expe+'5_'+name,expe+'6_'+name ]
			nameroot='-'+DataName+'-'+name+'_'+expe+'1'+expe+'2'+expe+'3'+expe+'4'+expe+'5'+expe+'6'
		elif expe in ['B','D','E','F']:
			plotexpe=[expe+'1_'+name,expe+'2_'+name,expe+'3_'+name,expe+'4_'+name,expe+'5_'+name ]
			nameroot='-'+DataName+'-'+name+'_'+expe+'1'+expe+'2'+expe+'3'+expe+'4'+expe+'5'
		elif expe in ['C']:
			plotexpe=[expe+'1_'+name,expe+'2_'+name,expe+'3_'+name,expe+'4_'+name]
			nameroot='-'+DataName+'-'+name+'_'+expe+'1'+expe+'2'+expe+'3'+expe+'4'
	else:
		plotexpe=[expe+'_'+name]
		nameroot='-'+DataName+'-'+name+'_'+expe
	#gathering necessarry info on what to plot
	for j,tag in enumerate(plotexpe):
		NCFile=resdir+tag+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			for key in DatFile.variables[DataName].dimensions:
				if key.startswith('index'):
					if key=='index_ch':
						datacoord='coords_ch'
					else:
						datacoord='coords'+str(key[-1])
					break
			DatFile.close()
		except (RuntimeError,IOError):
			continue

	os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
	if sys.argv[1] == 'All':
			for i,plot in enumerate(PlotList.values()):
				fig=plot()

				try:
					#fig.savefig(figdir+'/PDFs/'+name+'/'+PlotList.keys()[i]+nameroot+'.pdf',format='pdf')
					fig.savefig(figdir+'/'+name+'/'+PlotList.keys()[i]+nameroot+'.png',format='png')
				except IOError:
					print('No Figure directory existing where needed, creating it here "{}"'.format(os.path.abspath(figdir)))
					#os.makedirs(figdir+'/PDFs/'+name)
					os.makedirs(figdir+'/'+name)
					#fig.savefig(figdir+'/PDFs/'+name+'/'+PlotList.keys()[i]+nameroot+'.pdf',format='pdf')
					fig.savefig(figdir+'/'+name+'/'+PlotList.keys()[i]+nameroot+'.png',format='png')
	else:
		try:
			fig=PlotList[sys.argv[1]]()
			try:
				#fig.savefig(figdir+'/PDFs/'+name+'/'+sys.argv[1]+nameroot+'.pdf',format='pdf')
				fig.savefig(figdir+'/'+name+'/'+sys.argv[1]+nameroot+'.png',format='png')
			except IOError:
				print('No Figure directory existing where needed, creating it here "{}"'.format(os.path.abspath(figdir)))
				#os.makedirs(figdir+'/PDFs/'+name)
				os.makedirs(figdir+'/'+name)
				#fig.savefig(figdir+'/PDFs/'+name+'/'+sys.argv[1]+nameroot+'.pdf',format='pdf')
				fig.savefig(figdir+'/'+name+'/'+sys.argv[1]+nameroot+'.png',format='png')

		except KeyError:

			print 'Bad keyword for the figure production, run without keyword to get the help and check the spelling.'

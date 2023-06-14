#!/usr/bin/env python2
"""Plotter for the GHIP exercise, the plotter relies on the conventions as given on the comparison framework,
Available plot as of now are:
  -Sections (called with $:GHIPPloter AllSection)
  -Map Plot (called with $:GHIPPloter SynthMap)
"""

from netCDF4 import Dataset
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
from matplotlib.ticker import SymmetricalLogLocator
from matplotlib.ticker import MaxNLocator
from TreatMesh import *

# {{{ AllSections()

def AllSection():
	"""---Allsection---
	Description :
	Plot a variable function of the x coordinate for all the participating model and a given experiment.

	The experiment tag is asked later , if only the letter of the experiment is given then all the subexperiments are plotted
	"""
	#Open Figure

	fig ,ax= plt.subplots(1,len(plotexpe),sharex=True,sharey=True,figsize=(16,12))
	fig.patch.set_alpha(0.0)
	# ax.set_xlabel('Distance from front [km]', fontsize=14)
	# ax.set_ylabel('Effective Pressure [MPa]', fontsize=14)
	# ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	# ax.set_title('Comparison for run '+ expe, fontsize=18)
	# #Make space for the legend
	# box = ax.get_position()
	# ax.set_position([box.x0, box.y0, box.width*0.85 , box.height])

	maxval = np.nan
	minval = np.nan
	norm = mpl.colors.Normalize(vmin=0, vmax=len(names)-1)
	cmap = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.Accent)

	for i,name in enumerate(names):
		for j,tag in enumerate(plotexpe):
			NCFile=resdir+name+'/'+tag+name+'.nc'
			try:
				DatFile	 = Dataset(NCFile, mode='r')
				for key in DatFile.variables[DataName].dimensions:
					if key.startswith('index'):
						if key=='index_ch':
							datacoord='coords_ch'
						else:
							datacoord='coords'+str(key[-1])
				data, coords=GetVar(DatFile,DataName,datacoord,False,tag.endswith('mh1'))
				DatFile.close()
			except (RuntimeError,IOError):
				print 'File '+ str(NCFile)+' does not exist'
				continue 
		
			plotshade = cmap.to_rgba(i)
			# plotshade=mpl.colors.colorConverter.to_rgb(plotcolor[i])
			# plotshade=[max(j*0.18,x) for x in plotshade]
			minval=np.nanmin([np.nanmin(data),minval])
			maxval=np.nanmax([np.nanmax(data),maxval])
			ax.plot(coords[0,:], data, linestyle='',color=plotshade,marker=plotmark,label=name+'_'+tag)
	#and place the legend outside on right
	ax.legend(loc='center left',bbox_to_anchor=(1, 0.5))
	ax.axis([0,100e3,minval,maxval])
	print 'done'
	plt.show()
	return fig

# }}}
# {{{ SynthMap()
def SynthMap():
	totplots = (len(names))*(len(plotexpe))
	cols=len(plotexpe)
	suprows=len(names)-1	
	if expe in['A','B','C','D']:
		ratio=5
	elif expe in['E','F']:
		ratio=5.55
	fig, axar = plt.subplots(1+suprows, cols, sharex='col', sharey='row', figsize=(2.5+ratio*cols,1.5+(suprows+1)))
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
	tmplims=np.nan*np.ones((len(plotexpe)*len(names)+1,2))
	for i,name in enumerate(names):
		for j, tag in enumerate(plotexpe):
			NCFile=prefix+name+'/'+tag+name+'.nc'
			try:
				DatFile	 = Dataset(NCFile, mode='r')
				data, coords=GetVar(DatFile,DataName,datacoord,False)
				DatFile.close()
				tmplims[(j+1)*(i+1),0]=np.nanmin(data)
				tmplims[(j+1)*(i+1),1]=np.nanmax(data)
			except (RuntimeError,IOError):
				continue

	#defining colormap
	lims=[np.nanmin(tmplims[:,0]),np.nanmax(tmplims[:,1])]
	norm=mpl.colors.SymLogNorm(linthresh=0.2*np.diff(lims), linscale=2,vmin=lims[0], vmax=lims[1])
	#norm=MidpointNormalize(midpoint=2.)
	cmap = plt.cm.viridis
	
	for i,name in enumerate(names):
		for j, tag in enumerate(plotexpe):
			NCFile=prefix+name+'/'+tag+name+'.nc'
			print('Treating {}'.format(NCFile))
			axar[i,j].set_position([hleft+(hspace+width)*j,1-(vtop+height+(vspace+height)*i),width,height])
			if i==len(names)-1:
				axar[i,j].set_xlabel('X axis [m]', fontsize=14)
			if j==0:
				axar[i,j].set_ylabel('Y axis [m]', fontsize=14)
			axar[i,j].set_title('Simulation '+tag+name, fontsize=18)

			try:
				#get the data if the file exist
				DatFile	 = Dataset(NCFile, mode='r')
				data, coords=GetVar(DatFile,DataName,datacoord,False)
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
				print 'File '+ str(NCFile)+' does not exist'
				if expe in ['A','B','C','D']:
					axar[i,j].set_xlim([0,100])
					axar[i,j].set_ylim([0,20])
				elif expe in ['E','F']:
					axar[i,j].set_xlim([0,6])
					axar[i,j].set_ylim([-0.5,0.5])
				continue

			Xlims=[1.0e-3*np.min(coords[0,:]),1.0e-3*np.max(coords[0,:])]
			Ylims=[1.0e-3*np.min(coords[1,:]),1.0e-3*np.max(coords[1,:])]
			axar[i,j].set_xlim(Xlims)
			axar[i,j].set_ylim(Ylims)
			
			#Now Plotting
			grid_x, grid_y = np.mgrid[Xlims[0]:Xlims[1]:101j, Ylims[0]:Ylims[1]:101j]
			gridded=griddata(1.0e-3*coords.T, data, (grid_x, grid_y), method='linear')
			axar[i,j].imshow(gridded.T, extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]), 
													origin='lower',norm=norm,cmap=cmap)
			if i==len(names)-1 and j==len(plotexpe)-1:
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

				tick_locator = mpl.ticker.MaxNLocator(nbins=7)
				cbar.locator = tick_locator
				cbar.update_ticks()
	#And draw
	plt.show
	return fig

# }}}
# {{{ AllEnvelopes()

def AllEnvelopes():
	"""---AllEnvelopes---
	Description :

	Needed parameters :

	The experiment tag is asked later , if only the letter of the experiment is given then all the subexperiments are plotted
	
	Example of use :
	GHIPPlotter.AllEnvelopes('N')
	"""
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


	for name in names:
		SynthDat  = np.NaN*np.ones((4,len(plotexpe)))
		for j,tag in enumerate(plotexpe):
			NCFile=resdir+name+'/'+tag+name+'.nc'
			try:
				DatFile	 = Dataset(NCFile, mode='r')
				for key in DatFile.variables[DataName].dimensions:
					if key.startswith('index'):
						if key=='index_ch':
							datacoord='coords_ch'
						else:
							datacoord='coords'+str(key[-1])
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
			
			
			if expe in ['A','B','C','D']:
				Xlims=[max(0,np.min(coords[0,:])),min(1.0e5,np.max(coords[0,:]))]
			elif expe in ['E','F']:
				Xlims=[max(0,np.min(coords[0,:])),min(6.0e3,np.max(coords[0,:]))]
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
		
			plotshade=mpl.colors.colorConverter.to_rgb(namestyle[name])
			plotshade=[max(0.6,x) for x in plotshade]

			if j==3 or len(plotexpe)==1:
				ax1.set_ylabel('{} $[{}]$'.format(TitleList[DataName][0],TitleList[DataName][1]), fontsize=14)
				ax1.set_title('Overview for run '+ tag, fontsize=18)
				ax1.axis([0, 1.0e-3*Xlims[1], 0.0, 5])
				ax1.plot(1.0e-3*coords[0,:], data,color=plotshade,marker=plotmark,linestyle='')		
				ax1.plot(1.0e-3*xx, SmoothLow(xx),color=namestyle[name],linewidth=1.5,label=name)
				ax1.plot(1.0e-3*xx, SmoothUp(xx),color=namestyle[name],linewidth=1.5)
				ax1.xaxis.set_ticklabels([])

				ax2.axis([0, 1.0e-3*Xlims[1], -1.0, 1.0])
				ax2.set_xlabel('Distance from front [km]', fontsize=14)
				ax2.set_ylabel('{} \nSpread $[{}]$'.format(TitleList[DataName][0],TitleList[DataName][1]), fontsize=14)
				ax2.plot(1.0e-3*xx, SmoothUp(xx)-SmoothLow(xx),color=namestyle[name],linewidth=1.5,linestyle='-')
				
		#And draw
		indexes=range(1,len(plotexpe)+1,1)
		limits=[[0,14],[0,110],[0,7],[0,0.04]]

		for p in range(0,4):
			ay[p].plot(indexes,SynthDat[p,:],color=namestyle[name],marker=plotmark,linestyle='',markersize=15)
			ay[p].plot(indexes,SynthDat[p,:],linestyle='-',color=plotshade)
			ay[p].yaxis.tick_right()
			ay[p].yaxis.set_ticks_position('both')
			ay[p].yaxis.set_label_position('right')
		if tag.startswith('A'):
			ay[p].axis([1, 6, limits[p][0], limits[p][1]])
			ay[p].set_xticks([1,2,3,4,5,6])
			if p==3:
				ay[p].xaxis.set_ticklabels(['A1','A2','A3','A4','A5','A6'])
			else:
				ay[p].xaxis.set_ticklabels([])
		if tag.startswith('B'):
			ay[p].axis([1, 5, limits[p][0], limits[p][1]])
			ay[p].set_xticks([1,2,3,4,5])
			if p==3:
				ay[p].xaxis.set_ticklabels(['B1','B2','B3','B4','B5'])
			else:
				ay[p].xaxis.set_ticklabels([])
		if tag.startswith('C'):
			ay[p].axis([1, 4, limits[p][0], limits[p][1]])
			ay[p].set_xticks([1,2,3,4])
			if p==3:
				ay[p].xaxis.set_ticklabels(['C1','C2','C3','C4'])
			else:
				ay[p].xaxis.set_ticklabels([])
		if tag.startswith('D'):
			ay[p].axis([1, 5, limits[p][0], limits[p][1]])
			ay[p].set_xticks([1,2,3,4,5])
			if p==3:
				ay[p].xaxis.set_ticklabels(['D1','D2','D3','D4','D5'])
			else:
				ay[p].xaxis.set_ticklabels([])
		if tag.startswith('E'):
			ay[p].axis([1, 5, limits[p][0], limits[p][1]])
			ay[p].set_xticks([1,2,3,4,5])
			if p==3:
				ay[p].xaxis.set_ticklabels(['E1','E2','E3','E4','E5'])
			else:
				ay[p].xaxis.set_ticklabels([])
		if tag.startswith('F'):
			ay[p].axis([1, 5, limits[p][0], limits[p][1]])
			ay[p].set_xticks([1,2,3,4,5])
			if p==3:
				ay[p].xaxis.set_ticklabels(['F1','F2','F3','F4','F5'])
			else:
				ay[p].xaxis.set_ticklabels([])

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
# {{{ def GetVar(DatFile,AllTime)
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

if __name__=="__main__":
	PlotList={'Sections': AllSection,
						'Maps':SynthMap,
						'Enveloppes':AllEnvelopes}

	TitleList={'h':['distributed water layer thickness',u'm'],
						 'hstore':['stored water: effective layer thickness',u'm'],
						 'S':['channel cross-sectional area',u'm^2'],
						 'Q':['channel discharge',u'm^3/s'],
						 'q':['distributed water discharge (abs)',u'm^2/s'],
						 'H':['ice thickness',u'm'],
						 'B':['bedrock elevation',u'm'],
						 'N':['effective pressure',u'MPa']}
	if len(sys.argv)==1: # print help if no filename plot name is given:
		print __doc__
		sys.exit(0)

	resdir=os.getcwd()+'/'
	figdir=resdir+'/figures'

	DataName=raw_input('Give the variable you want to plot: ')
	expe=raw_input('Give the tag of the experiment you want to plot: ')

	# for trying in [1,2]:
	# 	if trying==1:
	# 		names = ['jsb','mh2','mh1','mw','og','sb']#,'og2'
	# 		suffix='all1'
	# 	elif trying==2:
	# 		names = ['as','bf','cdf','db','id','jd']
	# 		suffix='all2'

	names = ['db','cdf','id', 'bf','sb','jd','jsb','as','mh1','mh2','mw','og','og2'	]

	namestyle={'as':'r',
						 'bf':'b',
						 'cdf':'g',
						 'db':'k',
						 'id':'c',
						 'jd':'m',
						 'jsb':'r',
						 'mh2':'b',
						 'mh1':'g',
						 'mw':'k',
						 'og':'c',
						 'og2':'y',
						 'sb':'m'}#Add a new plot style when you add a new name
	plotmark='.'
	if len(expe)==1:
		if expe in ['A']:
			plotexpe=[expe+'1_',expe+'2_',expe+'3_',expe+'4_',expe+'5_',expe+'6_']
		elif expe in ['B','D','E','F']:
			plotexpe=[expe+'1_',expe+'2_',expe+'3_',expe+'4_',expe+'5_']
		elif expe in ['C']:
			plotexpe=[expe+'1_',expe+'2_',expe+'3_',expe+'4_']
		else:
			plotexpe=[expe+'_']

	os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
	if sys.argv[1] == 'All':

		for i,plot in enumerate(PlotList.values()):
			fig=plot()
			try:
				fig.savefig(figdir+'all'+PlotList.keys()[i]+'.png',format='png')
			except IOError:
				print('No Figure directory existing where needed, creating it here "{}"'.format(os.path.abspath(figdir)))
				os.makedirs(figdir+'/all')
				fig.savefig(figdir+'/all/'+PlotList.keys()[i]+'.png',format='png')
	else:
		try:
			print 'here'
			fig=PlotList[sys.argv[1]]()
			try:
				fig.savefig(figdir+'/all/'+sys.argv[1]+'.png',format='png')
			except IOError:
				print('No Figure directory existing where needed, creating it here "{}"'.format(os.path.abspath(figdir)))
				os.makedirs(figdir+'/all')
				fig.savefig(figdir+'/all/'+sys.argv[1]+'.png',format='png')

		except KeyError:
			print 'Bad keyword for the figure production, run without keyword to get the help and check the spelling.'

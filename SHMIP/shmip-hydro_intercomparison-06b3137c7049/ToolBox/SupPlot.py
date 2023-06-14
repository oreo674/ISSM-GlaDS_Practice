#!/usr/bin/env python2
"""
Max
"""
import os
import sys
import csv
import string
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from glob import glob
from matplotlib.ticker import MaxNLocator
from matplotlib import cbook
from matplotlib import ticker
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
from matplotlib.patches import ConnectionPatch
from matplotlib.collections import PatchCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d.art3d import PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from netCDF4 import Dataset

from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'size':9})
rc('legend',scatterpoints=1)

# {{{ Transient(run)

def Transient(runnum):
	tag=expe+str(runnum)

	width=17.8/2.54 #two columns JoG
	plotrow=3
	plotcol=2

	if expe in ['C','D']:
		ratio=1.8
		#0.3 is the width/heigh ratio of the maps (without borders)
		height=0.3*width*ratio
		Xlims=[0,100]
		Xticks=[10,30,50,70,90]
		Xlabels=['$10$','$30$','$50$','$70$','$90$']
		Ylims=[0,20]
		Yticks=[5,10,15]
		Ylabels=['$5$','$10$','$15$']
		midcoord=10.0e3

		#horizontal buffer and related spaces
		hbuff=0.12
		hleft=0.45*hbuff
		hright=0.45*hbuff
		hspace=0.1*hbuff
		#vertical buffer and related spaces
		vbuff=1.-((1.-hbuff)/ratio)
		vbot=0.39*vbuff
		vtop=0.17*vbuff
		vspace=(0.44*vbuff)/2.


	if expe=='F':
		ratio=1.95
		#1/6 is the width/heigh ratio of the maps (without borders)
		height=0.25*width*ratio
		Xlims=[0,6]
		Xticks=[1,2,3,4,5]
		Xlabels=['$1$','$2$','$3$','$4$','$5$']
		Ylims=[-0.5,0.5]
		Yticks=[-0.3,0,0.3]
		Ylabels=['$-0.3$','$0$','$0.3$']
		midcoord=0.0

		#horizontal buffer and related spaces
		hbuff=0.16
		hleft=0.45*hbuff
		hright=0.45*hbuff
		hspace=0.1*hbuff
		#vertical buffer and related spaces
		vbuff=1.-((1.-hbuff)/ratio)
		vbot=0.42*vbuff
		vtop=0.16*vbuff
		vspace=(0.42*vbuff)/2.
		#subplot sizes

	plotwidth=(1.-hleft-hright-hspace)/plotcol
	plotheight=(1.0-vbot-vtop-2*vspace)/plotrow

	AmpList={'C1':0.25,'C2':0.5,'C3':1.0,'C4':2.0,
					 'D1':-4.0,'D2':-2.0,'D3':0.0,'D4':2.0,'D5':4.0,
					 'F1':-6.0,'F2':-3.0,'F3':0.0,'F4':3.0,'F5':6.0}
	BandProp={'Low':[10000.,600.,'b','$Low$','lower band'],
						'Mid':[50000.,3000.,'c','$Mid$','middle band'],
						'Top':[90000.,5100.,'r','$Top$','highest band']}

	if expe=='C':
		moulinfile='/home/bfl022/Model/GHIP/hydro_intercomparison/input_functions/source/B5_M.csv'
		MoulinInput=np.zeros((101))
		with open(moulinfile) as csvfile:
			MoulinPos=csv.reader(csvfile,delimiter=',')
			for row in MoulinPos:
				Xs=int(row[1])
				MoulinInput[Xs/1000]+=float(row[3])

	for name in names:
		PlotFlux=False
		NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
		print NCFile
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			time     = DatFile.variables['time'][:]
			data,coords,datacoordtype=GetVar(DatFile,DataName,name,True)
			thick,coords,thickcoordtype=GetVar(DatFile,'H',name,False)
			base,coords,basecoordtype=GetVar(DatFile,'B',name,False)
			DatFile.close()
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))
			continue
		#===layout and some axis stuff
		layout,ax=plt.subplots(plotrow,plotcol,figsize=(width,height))
		for i in[0,2]:
			for j in [0,1]:
				#===X axis
				ax[i,j].set_xlim(Xlims)
				ax[i,j].set_xticks(Xticks)
				ax[i,j].xaxis.set_ticklabels(Xlabels)
				ax[i,j].set_xlabel('$x$ (km)',labelpad=2)
				if i==2:
					ax[i,j].xaxis.tick_top()
					ax[i,j].xaxis.set_ticks_position('both')
					ax[i,j].xaxis.set_label_position('top')
					ax[i,j].tick_params(axis='x', which='major', pad=3)
					ax[i,j].set_xlabel('$x$ (km)',labelpad=4)
				#===Y axis
				ax[i,j].set_ylim(Ylims)
				ax[i,j].set_yticks(Yticks)
				ax[i,j].yaxis.set_ticklabels(Ylabels)
				ax[i,j].set_ylabel('$y$ (km)',labelpad=2)
				if j==1:
					ax[i,j].yaxis.tick_right()
					ax[i,j].yaxis.set_ticks_position('both')
					ax[i,j].yaxis.set_label_position('right')

		#===defining colormap
		lims=[np.nanmin(data),np.nanmax(data)]
		# norm,cmap=CenteredSymLogNorm(plt.cm.bwr,TreshRatio=0.2, linlenght=1,minval=lims[0], maxval=lims[1])
		# norm=mpl.colors.SymLogNorm(linthresh=0.2*np.diff(lims), linscale=2,vmin=lims[0], vmax=lims[1])
		# cmap = plt.cm.viridis
		norm=mpl.colors.Normalize(vmin=lims[0], vmax=lims[1])
		cmap=shiftedColorMap(plt.cm.bwr, start=0, midpoint= 1-lims[1]/(lims[1]-lims[0]), stop=1.0, name='shiftedcmap')
		#===geometry related
		if expe=='F':
			BenchLength=6.0e3
			BenchStep=60.
			BenchLines=np.arange(0,BenchLength+BenchStep,BenchStep)
			s2 = 100.0/BenchLength
			s_xend = 100.0*(BenchLength+200)**0.25 + s2*BenchLength - 100.0*200**0.25
			BenchSurf = 100.0*(BenchLines+200)**0.25 + s2*BenchLines - 100.0*200**0.25
			f10 = (s_xend - 0.05*BenchLength)/BenchLength**2
			f0 = f10*BenchLines**2 + 0.05*BenchLines
			h0 = -4.5*(BenchLines/BenchLength) + 5
			BenchWidth = 2.0*(((BenchSurf-f0)/h0)/0.5e-6)**(1/3.)
			boundaries=np.vstack((np.hstack((BenchLines,np.flipud(BenchLines)))*1.0e-3,
														np.hstack((BenchWidth,np.flipud(-BenchWidth)))*1.0e-3*0.5)).T
			#same but with data shape
			tmp0= 100.0*(coords[0,:]+200)**0.25 + s2*coords[0,:] - 100.0*200**0.25
			tmp1 = f10*coords[0,:]**2 + 0.05*coords[0,:]
			tmp2 = -4.5*(coords[0,:]/BenchLength) + 5
			dataWidth = 2.0*(((tmp0-tmp1)/tmp2)/0.5e-6)**(1/3.)


		#===for the center time line
		ax[1,0].set_position([hleft+(0.25*plotwidth),
													vbot+plotheight+0.9*vspace,
													1.2*plotwidth,0.8*plotheight])
		ax[1,1].set_position([hleft+(0.25*plotwidth),
													vbot+1.9*plotheight+0.9*vspace,
													1.2*plotwidth,0.4*plotheight])
		ax[1,1].set_zorder(2)
		ax[1,0].set_zorder(2)

		if expe=='C':
			time=24.+(time-max(time))/3600.
			if name in ['id','jd','mh2','mw','mw_prime']:
				time=np.insert(time,0,0)
			Cycle=DiurnalForcing(time*3600.,AmpList[tag])
			Forcing=np.outer(MoulinInput,Cycle)
			for mid in [0,1]:
				ax[1,mid].set_xlim([0 ,24])
				ax[1,mid].set_xticks([3,6,9,12,15,18,21,24])
				ax[1,mid].xaxis.set_ticklabels([])
				ax[1,mid].yaxis.set_major_locator(MaxNLocator(2))
				if mid==0:
					ax[1,mid].yaxis.set_major_locator(MaxNLocator(4,prune='upper'))
					ax[1,mid].xaxis.set_ticklabels(['$3$','$6$','$9$','$12$','$15$','$18$','$21$','$24$'])
					ax[1,mid].set_xlabel('Time (h)',labelpad=-1,x=0.6)

		else:
			time=1+(time-max(time))/(365.*24.*3600.)
			time=time[np.where(time>=0)][:]
			surf=thick+base
			Forcing=SeasonForcing(surf,time*3600.*24.*365.,AmpList[tag])
			for mid in [0,1]:
				ax[1,mid].set_xlim([0.0 ,1])
				ax[1,mid].set_xticks([0.0,0.16,0.33,0.49,0.66,0.83,1.0])
				ax[1,mid].xaxis.set_ticklabels([])
				ax[1,mid].yaxis.set_major_locator(MaxNLocator(2))
				if mid==0:
					ax[1,mid].yaxis.set_major_locator(MaxNLocator(4,prune='upper'))
					ax[1,mid].xaxis.set_ticklabels(['$0$','$2$','$4$','$6$','$8$','$10$','$12$'])
					ax[1,mid].set_xlabel('Time (month)',labelpad=-1,x=0.6)

		Band=np.NaN*np.ones((len(time),3,3))
		Mean=np.NaN*np.ones((len(time)))
		Meanwidth=np.NaN*np.ones((3))
		Force=np.NaN*np.ones((len(time),3))
		MeanForce=np.NaN*np.ones((len(time)))

		nodeselect=[]
		if name=='db':
			for t,step in enumerate(time):
				Mean[t]=np.nanmean(data[t,:])
				if expe=='C':
					MeanForce[t]=np.nansum(Forcing[:,t])
				elif expe=='D':
					MeanForce[t]=np.nanmean(Forcing[:,t])*100.0e3*20.0e3
				elif expe=='F':
					MeanForce[t]=np.nanmean(Forcing[:,t]*dataWidth)*6.0e3
		else:
			for b,bandname in enumerate(['Low','Mid','Top']):
				if expe=='F':
					nodeselect.append(np.where(np.logical_and(coords[0,:]>=BandProp[bandname][1],
																										coords[0,:]<=BandProp[bandname][1]+300))[0])
				else:
					nodeselect.append(np.where(np.logical_and(coords[0,:]>=BandProp[bandname][0],
																										coords[0,:]<=BandProp[bandname][0]+5000))[0])
			if name in ['id','jd','mh2','mw','mw_prime'] and expe=='C':
				data=np.insert(data,0,data[-1,:],axis=0)

			for t,step in enumerate(time):
				if expe=='C':
					MeanForce[t]=np.nansum(Forcing[:,t])
				elif expe=='D':
					MeanForce[t]=np.nanmean(Forcing[:,t])*100.0e3*20.0e3
				elif expe=='F':
					MeanForce[t]=np.nanmean(Forcing[:,t]*dataWidth)*6.0e3
				for b,bandname in enumerate(['Low','Mid','Top']):
					Band[t,b,0]=np.nanmin(data[t,nodeselect[b]])
					Band[t,b,1]=np.nanmean(data[t,nodeselect[b]])
					Band[t,b,2]=np.nanmax(data[t,nodeselect[b]])
					Mean[t]=np.nanmean(data[t,:])

		if name=='db':
			pressure=ax[1,0].plot(time, Mean[:], linestyle='-',linewidth=1.5,color='g',label='Mean $N$')[0]
			recharge=ax[1,1].plot(time, MeanForce, linestyle='-',linewidth=1.5,color='k',label='Total recharge')[0]
			tenPspread=0.1*(np.nanmax(Mean)-np.nanmin(Mean))
			ax[1,0].set_ylim([np.nanmin(Mean)-tenPspread ,np.nanmax(Mean)+tenPspread])
			ax[1,0].set_ylabel('$N$ (MPa)',rotation=0,ha='right',va='center')
			ax[1,1].set_ylabel('Recharge \n (m$^3$s$^{-1}$)',rotation=0,ha='right',va='center')
			plt.legend(handles=[recharge,pressure],bbox_to_anchor=(hleft+1.45*plotwidth+hspace,vbot+plotheight+0.8*vspace,0.55*plotwidth,1.3*plotheight),
								 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.2,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,fontsize=9,
								 borderaxespad=0,mode='expand', frameon=False)
		else:
			pressure=[]
			for b,bandname in enumerate(['Low','Mid','Top']):
				ax[1,0].plot(time, Band[:,b,1], linestyle='-',linewidth=1.5,color=BandProp[bandname][2])
				recharge=ax[1,1].plot(time, MeanForce, linestyle='-',linewidth=1.5,color='k',label='Total recharge')[0]
				pressure.append(ax[1,0].fill_between(time,Band[:,b,0],Band[:,b,2],
																						 facecolor=BandProp[bandname][2],
																						 edgecolor=BandProp[bandname][2],alpha=0.3,
																						 label='$N$ for '+BandProp[bandname][4]))
			tenPspread=0.1*(np.nanmax(data)-np.nanmin(data))
			ax[1,0].set_ylim([np.nanmin(data)-tenPspread ,np.nanmax(data)+tenPspread])
			ax[1,0].set_ylabel('$N$ (MPa)',rotation=0,ha='right',va='center')
			ax[1,1].set_ylabel('Recharge \n (m$^3$s$^{-1}$)',rotation=0,ha='right',va='center')
			plt.legend(handles=[recharge,pressure[0],pressure[1],pressure[2]],
								 bbox_to_anchor=(hleft+1.45*plotwidth+hspace,vbot+plotheight+0.8*vspace,0.55*plotwidth,1.3*plotheight),
								 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.6,columnspacing=0.2,
								 borderpad=0.3,handletextpad=0.2,fontsize=9,borderaxespad=0,
								 handler_map={pressure[0]:FillandLineHandler('b'),
															pressure[1]:FillandLineHandler('c'),
															pressure[2]:FillandLineHandler('r')},mode='expand', frameon=False)

		#===N Map
		if expe in['D','F']:
			#get a set of dates
			dates=[min(np.where(MeanForce>1.6e-1)[0]),#first melt day
						 max(np.where(MeanForce>1.6e-1)[0]),#last melt day
						 min(np.where(Mean==np.nanmin(Mean))[0]),#minimum effective pressure
						 min(np.where(MeanForce==np.nanmax(MeanForce))[0])]#maximum water input

			datedict={str(dates[0]):'First melt, day {}'.format(str(dates[0])),
								str(dates[1]):'Last melt, day {}'.format(str(dates[1])),
								str(dates[2]):'Minimum $N$, day {}'.format(str(dates[2])),
								str(dates[3]):'Maximum recharge, day {}'.format(str(dates[3]))}
		elif expe=='C':
			#get a set of dates
			minInput=max(np.where(MeanForce==np.nanmin(MeanForce))[0])
			maxInput=min(np.where(MeanForce==np.nanmax(MeanForce))[0])
			dates=[minInput,#first melt day
						 maxInput,#maximum water input
						 min(np.where(Mean==np.nanmin(Mean))[0]),#minimum effective pressure
						 int(minInput+0.5*(maxInput-minInput))]#mid input period

			datedict={str(dates[0]):'Minimum recharge, {} hours'.format(str(dates[0])),
								str(dates[1]):'Maximum recharge, {} hours'.format(str(dates[1])),
								str(dates[2]):'Minimum $N$, {} hours'.format(str(dates[2])),
								str(dates[3]):'Median recharge, {} hours'.format(str(dates[3]))}
		dates=np.sort(dates)


		#loop on the dates
		for m in np.arange(0,4):
			row=2-2*(m/2)
			col=(np.mod(m,2))
			ax[row,col].set_position([hleft+col*(hspace+plotwidth),
																vbot+row*(plotheight+vspace),
																plotwidth,plotheight])

			# get the limit and plot vertical date line
			centerbotYlims=ax[1,0].get_ylim()
			ax[1,0].plot((time[dates[m]],time[dates[m]]), centerbotYlims, linestyle='-',linewidth=1.5,color=(0.7,0.7,0.7))
			centertopYlims=ax[1,1].get_ylim()
			ax[1,1].plot((time[dates[m]],time[dates[m]]), centertopYlims, linestyle='-',linewidth=1.5,color=(0.7,0.7,0.7))

			if np.all(coords[1,:]==midcoord):
				if name in ['mh1','cdf']:
					timedata=np.flipud(data[dates[m],:])
				else:
					timedata=(data[dates[m],:])
				mapim=ax[row,col].imshow(np.asarray([timedata]), extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
														 origin='lower',norm=norm,cmap=cmap,aspect='auto')
			else:
				griddeddata,gridpoint_x,gridpoint_y=Griding(data[dates[m],:],coords)
				mapim=ax[row,col].imshow(griddeddata.T, extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
													 origin='lower',norm=norm,cmap=cmap)

			#clip and outline Bench like
			if expe=='F':
				patch = Polygon(boundaries,transform=ax[row,col].transData)
				mapim.set_clip_path(patch)
				ax[row,col].plot(BenchLines*1.0e-3,0.5*BenchWidth*1.0e-3,color=(0.7,0.7,0.7))
				ax[row,col].plot(BenchLines*1.0e-3,-0.5*BenchWidth*1.0e-3,color=(0.7,0.7,0.7))
				if name!='db':
					lowband=mpl.patches.Rectangle([0.6, -0.55], 0.3, 1.1, facecolor='b',alpha=1,
																				edgecolor='b', transform=ax[row,col].transData)
					midband=mpl.patches.Rectangle([3, -0.55], 0.3, 1.1, facecolor='c',alpha=1,
																				edgecolor='c', transform=ax[row,col].transData)
					topband=mpl.patches.Rectangle([5.1, -0.55], 0.3, 1.1, facecolor='r',alpha=1,
																				edgecolor='r', transform=ax[row,col].transData)
			else:
				if name!='db':
					lowband=mpl.patches.Rectangle([10, -1], 5, 22, facecolor='b',alpha=1,
																				edgecolor='b', transform=ax[row,col].transData)
					midband=mpl.patches.Rectangle([50, -1], 5, 22, facecolor='c',alpha=1,
																				edgecolor='c', transform=ax[row,col].transData)
					topband=mpl.patches.Rectangle([85, -1], 5, 22, facecolor='r',alpha=1,
																				edgecolor='r', transform=ax[row,col].transData)
			if name!='db':
				lowband.set_zorder(0)
				midband.set_zorder(0)
				topband.set_zorder(0)
				layout.patches.append(lowband)
				layout.patches.append(midband)
				layout.patches.append(topband)

			#===Create the arrow
			# 1. Get transformation operators for axis and figure
			axtimetr = ax[1,0].transData # transformation from time axis
			axuptimetr = ax[1,1].transData # transformation from top time axis
			axmaptr = ax[row,col].transData # transformation from map axis
			figtr = layout.transFigure.inverted() # Display -> Figure
			# 2. Transform arrow point from axis to figure coordinates
			if row==2:
				startpoint = figtr.transform(axuptimetr.transform((time[dates[m]], centertopYlims[1])))
				if col==1:
					labtext=ax[1,1].text(1.,-0.2,datedict[str(dates[m])],fontsize=9,va='center',ha='right',transform=ax[row,col].transAxes,zorder=5)
					endpoint = figtr.transform(axmaptr.transform((0.3*Xlims[1], Ylims[0])))
				else:
					labtext=ax[1,1].text(0.,-0.2,datedict[str(dates[m])],fontsize=9,va='center',ha='left',transform=ax[row,col].transAxes,zorder=5)
					endpoint = figtr.transform(axmaptr.transform((0.7*Xlims[1], Ylims[0])))
			else:
				startpoint = figtr.transform(axtimetr.transform((time[dates[m]], centerbotYlims[0])))
				if col==1:
					labtext=ax[1,1].text(1.,1.2,datedict[str(dates[m])],fontsize=9,va='center',ha='right',transform=ax[row,col].transAxes,zorder=5)
					endpoint = figtr.transform(axmaptr.transform((0.7*Xlims[1], Ylims[1])))
				else:
					labtext=ax[1,1].text(0.,1.2,datedict[str(dates[m])],fontsize=9,va='center',ha='left',transform=ax[row,col].transAxes,zorder=5)
					endpoint = figtr.transform(axmaptr.transform((0.5*Xlims[1], Ylims[1])))
			# 4. Create the patch
			arrow = mpl.patches.FancyArrowPatch(startpoint, endpoint,
																					transform=layout.transFigure,
																					fc = (0.7,0.7,0.7), ec=(0.7,0.7,0.7),connectionstyle="arc3,rad=0.0",#"arc3,rad=0.2",
																					arrowstyle='simple,tail_width=0.1,head_length=0.7', alpha = 1,shrinkA=0,
																					mutation_scale = 10.)

			# 5. Add patch to list of objects to draw onto the figure
			layout.patches.append(arrow)
			arrow.set_zorder(1)
			ax[row,col].set_zorder(0)


		#===Colorbar
		barNax = layout.add_axes([hleft,0.25*vbot,2.*plotwidth+hspace,0.2*plotheight])
		Nbar = mpl.colorbar.ColorbarBase(barNax, norm=norm,cmap=cmap, orientation='horizontal')
		Nbar.locator = ticker.MaxNLocator(nbins=9)
		Nbar.update_ticks()
		Nbar.set_label('Effective pressure (MPa)', labelpad=-25,ha='center',va='bottom')
		Nbar.solids.set_edgecolor("face")

		layout.savefig(figdir+'/'+name+'/Transient_'+DataName+'_'+name+'_'+tag+'.pdf',format='pdf')
	plt.show
	return 0

# }}}

# {{{ TransientFlux(run)

def TransientFlux(runnum):
	tag=expe+str(runnum)
	width=17.8/2.54 #two columns JoG
	plotrow=3
	plotcol=2

	if expe in ['C','D']:
		ratio=1.8
		#0.3 is the width/heigh ratio of the maps (without borders)
		height=0.3*width*ratio
		Xlims=[0,100]
		Xticks=[10,30,50,70,90]
		Xlabels=['$10$','$30$','$50$','$70$','$90$']
		Ylims=[0,20]
		Yticks=[5,10,15]
		Ylabels=['$5$','$10$','$15$']
		midcoord=10.0e3

		#horizontal buffer and related spaces
		hbuff=0.12
		hleft=0.45*hbuff
		hright=0.45*hbuff
		hspace=0.1*hbuff
		#vertical buffer and related spaces
		vbuff=1.-((1.-hbuff)/ratio)
		vbot=0.39*vbuff
		vtop=0.17*vbuff
		vspace=(0.44*vbuff)/2.


	if expe=='F':
		ratio=1.95
		#1/6 is the width/heigh ratio of the maps (without borders)
		height=0.25*width*ratio
		Xlims=[0,6]
		Xticks=[1,2,3,4,5]
		Xlabels=['$1$','$2$','$3$','$4$','$5$']
		Ylims=[-0.5,0.5]
		Yticks=[-0.3,0,0.3]
		Ylabels=['$-0.3$','$0$','$0.3$']
		midcoord=0.0

		#horizontal buffer and related spaces
		hbuff=0.16
		hleft=0.45*hbuff
		hright=0.45*hbuff
		hspace=0.1*hbuff
		#vertical buffer and related spaces
		vbuff=1.-((1.-hbuff)/ratio)
		vbot=0.42*vbuff
		vtop=0.16*vbuff
		vspace=(0.42*vbuff)/2.
		#subplot sizes

	plotwidth=(1.-hleft-hright-hspace)/plotcol
	plotheight=(1.0-vbot-vtop-2*vspace)/plotrow

	AmpList={'C1':0.25,'C2':0.5,'C3':1.0,'C4':2.0,
					 'D1':-4.0,'D2':-2.0,'D3':0.0,'D4':2.0,'D5':4.0,
					 'F1':-6.0,'F2':-3.0,'F3':0.0,'F4':3.0,'F5':6.0}
	BandProp={'Low':[10000.,600.,'b','$Low$','lower band'],
						'Mid':[50000.,3000.,'c','$Mid$','middle band'],
						'Top':[90000.,5100.,'r','$Top$','highest band']}

	if expe=='C':
		moulinfile='/home/bfl022/Model/GHIP/hydro_intercomparison/input_functions/source/B5_M.csv'
		MoulinInput=np.zeros((101))
		with open(moulinfile) as csvfile:
			MoulinPos=csv.reader(csvfile,delimiter=',')
			for row in MoulinPos:
				Xs=int(row[1])
				MoulinInput[Xs/1000]+=float(row[3])

	for name in names:
		NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
		print NCFile
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			time     = DatFile.variables['time'][:]
			data,coords,datacoordtype=GetVar(DatFile,DataName,name,True)
			thick,coords,thickcoordtype=GetVar(DatFile,'H',name,False)
			base,coords,basecoordtype=GetVar(DatFile,'B',name,False)

			if name in['mh2','mw','mw_prime','og','og_prime','mh1','cdf','rh','bf']:
				EffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'Q',name,True)
				EffFlux=np.abs(EffFlux)
				if effcoordtype=='coords_ch' or name=='rh':
					try:
						elts = DatFile.variables['connect_ch'][:,:].T.astype(int)
						if name=='mh1':
							fluxcoords=DatFile.variables['coords3'][:,:]
						else:
							fluxcoords=DatFile.variables['coords1'][:,:]
						x    = np.asarray([fluxcoords[0,elts[:,0]]*1.0e-3,fluxcoords[0,elts[:,1]]*1.0e-3])
						y    = np.asarray([fluxcoords[1,elts[:,0]]*1.0e-3,fluxcoords[1,elts[:,1]]*1.0e-3])
						scaterplot=False
					except KeyError:
						scaterplot=True
			elif name=='id':
				Channelisation,Ineffcoords,Inefcoordtype=GetVar(DatFile,'R',name,True)
				EffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'Q',name,True)
				EffFlux=EffFlux*Channelisation
				elts = DatFile.variables['connect_ch'][:,:].T.astype(int)
				fluxcoords=DatFile.variables['coords1'][:][0]
				x   = np.asarray([fluxcoords[elts[:,0]]*1.0e-3,fluxcoords[elts[:,1]]*1.0e-3])
				if expe in ['E','F']:
					y = np.zeros(np.shape(x))
				else:
					y = 10*np.ones(np.shape(x))
				scaterplot=False
			elif name=='as':
				IneffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'q',name,True)
				Channelisation,Effcoords,Effcoordtype=GetVar(DatFile,'Dc',name,True)
				EffFlux=IneffFlux*Channelisation
			elif name=='db':
				IneffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'q',name,True)
				Channelisation,Effcoords,Effcoordtype=GetVar(DatFile,'R',name,True)
				EffFlux=IneffFlux*Channelisation
			elif name=='sb':
				IneffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'q',name,True)
				Conductivity,Ineffcoords,Condcoordtype=GetVar(DatFile,'T',name,True)
				ischannel=np.ones(np.shape(Conductivity))
				ischannel[np.where(Conductivity==np.nanmin(Conductivity))]=0
				EffFlux=IneffFlux*ischannel
			elif name=='sb_old':
				IneffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'q',name,True)
				Conductivity,Ineffcoords,Condcoordtype=GetVar(DatFile,'K',name,True)
				ischannel=np.ones(np.shape(Conductivity))
				ischannel[np.where(Conductivity==np.nanmin(Conductivity))]=0
				EffFlux=IneffFlux*ischannel
			else:
				DatFile.close()
				continue
			DatFile.close()
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))
			continue

		#===layout and some axis stuff
		layout,ax=plt.subplots(plotrow,plotcol,figsize=(width,height))
		for i in[0,2]:
			for j in [0,1]:
				#===X axis
				ax[i,j].set_xlim(Xlims)
				ax[i,j].set_xticks(Xticks)
				ax[i,j].xaxis.set_ticklabels(Xlabels)
				ax[i,j].set_xlabel('$x$ (km)',labelpad=2)
				if i==2:
					ax[i,j].xaxis.tick_top()
					ax[i,j].xaxis.set_ticks_position('both')
					ax[i,j].xaxis.set_label_position('top')
					ax[i,j].tick_params(axis='x', which='major', pad=3)
					ax[i,j].set_xlabel('$x$ (km)',labelpad=4)
				#===Y axis
				ax[i,j].set_ylim(Ylims)
				ax[i,j].set_yticks(Yticks)
				ax[i,j].yaxis.set_ticklabels(Ylabels)
				ax[i,j].set_ylabel('$y$ (km)',labelpad=2)
				if j==1:
					ax[i,j].yaxis.tick_right()
					ax[i,j].yaxis.set_ticks_position('both')
					ax[i,j].yaxis.set_label_position('right')

		#===geometry related
		if expe=='F':
			BenchLength=6.0e3
			BenchStep=60.
			BenchLines=np.arange(0,BenchLength+BenchStep,BenchStep)
			s2 = 100.0/BenchLength
			s_xend = 100.0*(BenchLength+200)**0.25 + s2*BenchLength - 100.0*200**0.25
			BenchSurf = 100.0*(BenchLines+200)**0.25 + s2*BenchLines - 100.0*200**0.25
			f10 = (s_xend - 0.05*BenchLength)/BenchLength**2
			f0 = f10*BenchLines**2 + 0.05*BenchLines
			h0 = -4.5*(BenchLines/BenchLength) + 5
			BenchWidth = 2.0*(((BenchSurf-f0)/h0)/0.5e-6)**(1/3.)
			boundaries=np.vstack((np.hstack((BenchLines,np.flipud(BenchLines)))*1.0e-3,
														np.hstack((BenchWidth,np.flipud(-BenchWidth)))*1.0e-3*0.5)).T
			#same but with data shape
			tmp0= 100.0*(coords[0,:]+200)**0.25 + s2*coords[0,:] - 100.0*200**0.25
			tmp1 = f10*coords[0,:]**2 + 0.05*coords[0,:]
			tmp2 = -4.5*(coords[0,:]/BenchLength) + 5
			dataWidth = 2.0*(((tmp0-tmp1)/tmp2)/0.5e-6)**(1/3.)


		#===for the center time line
		ax[1,0].set_position([hleft+0.25*plotwidth,
													vbot+plotheight+0.9*vspace,
													1.2*plotwidth,0.8*plotheight])
		ax[1,1].set_position([hleft+0.25*plotwidth,
													vbot+1.9*plotheight+0.9*vspace,
													1.2*plotwidth,0.4*plotheight])
		ax[1,1].set_zorder(2)
		ax[1,0].set_zorder(2)

		if expe=='C':
			time=24.+(time-max(time))/3600.
			if name in ['id','jd','mh2','mw','mw_prime']:
				time=np.insert(time,0,0)
			Cycle=DiurnalForcing(time*3600.,AmpList[tag])
			Forcing=np.outer(MoulinInput,Cycle)
			for mid in [0,1]:
				ax[1,mid].set_xlim([0 ,24])
				ax[1,mid].set_xticks([3,6,9,12,15,18,21,24])
				ax[1,mid].xaxis.set_ticklabels([])
				ax[1,mid].yaxis.set_major_locator(MaxNLocator(2))
				if mid==0:
					ax[1,mid].yaxis.set_major_locator(MaxNLocator(4,prune='upper'))
					ax[1,mid].xaxis.set_ticklabels(['$3$','$6$','$9$','$12$','$15$','$18$','$21$','$24$'])
					ax[1,mid].set_xlabel('Time (h)',labelpad=-1,x=0.6)
		else:
			time=1+(time-max(time))/(365.*24.*3600.)
			time=time[np.where(time>=0)][:]
			#time=time[np.where(time>0)][:]
			surf=thick+base
			Forcing=SeasonForcing(surf,time*3600.*24.*365.,AmpList[tag])
			for mid in [0,1]:
				ax[1,mid].set_xlim([0.0 ,1])
				ax[1,mid].set_xticks([0.0,0.16,0.33,0.49,0.66,0.83,1.0])
				ax[1,mid].xaxis.set_ticklabels([])
				ax[1,mid].yaxis.set_major_locator(MaxNLocator(2))
				if mid==0:
					ax[1,mid].yaxis.set_major_locator(MaxNLocator(4,prune='upper'))
					ax[1,mid].xaxis.set_ticklabels(['$0$','$2$','$4$','$6$','$8$','$10$','$12$'])
					ax[1,mid].set_xlabel('Time (month)',labelpad=-1,x=0.6)

		Band=np.NaN*np.ones((len(time),3,3))
		Mean=np.NaN*np.ones((len(time)))
		FluxBand=np.NaN*np.ones((len(time),3,3))
		FluxMean=np.NaN*np.ones((len(time)))
		Meanwidth=np.NaN*np.ones((3))
		Force=np.NaN*np.ones((len(time),3))
		MeanForce=np.NaN*np.ones((len(time)))

		channelflux,sheetflux,lines=GetFluxRatio(name,tag,time)
		for t,step in enumerate(time):
			if expe=='C':
				MeanForce[t]=np.nansum(Forcing[:,t])
			elif expe=='D':
				MeanForce[t]=np.nanmean(Forcing[:,t])*100.0e3*20.0e3
			elif expe=='F':
				MeanForce[t]=np.nanmean(Forcing[:,t]*dataWidth)*6.0e3

		nodeselect=[]
		bandloc=[]
		if name=='db':
			for t,step in enumerate(time):
				Mean[t]=np.nanmean(data[t,:])
				FluxMean[t]=np.nanmean(EffFlux[t,:])
		else:
			for b,bandname in enumerate(['Low','Mid','Top']):
				if expe=='F':
					bandloc.append(np.where(lines==BandProp[bandname][1])[0])
					nodeselect.append(np.where(np.logical_and(coords[0,:]>=BandProp[bandname][1],
																										coords[0,:]<=BandProp[bandname][1]+300))[0])
				else:
					bandloc.append(np.where(lines==BandProp[bandname][0])[0])
					nodeselect.append(np.where(np.logical_and(coords[0,:]>=BandProp[bandname][0],
																										coords[0,:]<=BandProp[bandname][0]+5000))[0])
			if name in ['id','jd','mh2','mw','mw_prime'] and expe=='C':
				data=np.insert(data,0,data[-1,:],axis=0)
				EffFlux=np.insert(EffFlux,0,EffFlux[-1,:],axis=0)

			for t,step in enumerate(time):
				Mean[t]=np.nanmean(data[t,:])
				FluxMean[t]=np.nanmean(EffFlux[t,:])
				for b,bandname in enumerate(['Low','Mid','Top']):
					Band[t,b,0]=np.nanmin(data[t,nodeselect[b]])
					Band[t,b,1]=np.nanmean(data[t,nodeselect[b]])
					Band[t,b,2]=np.nanmax(data[t,nodeselect[b]])

		if name=='db':
			recharge=ax[1,1].plot(time, MeanForce, linestyle='-',linewidth=1.5,color='k',label='Total recharge')[0]
			ax[1,1].set_ylabel('Recharge \n (m$^3$s$^{-1}$)',rotation=0,ha='right',va='center')
			pressure=ax[1,0].plot(time, np.nanmean(channelflux,axis=1), linestyle='-',linewidth=1.5,color='g',label='Mean efficient discharge')[0]
			tenPspread=0.1*(np.nanmax(channelflux)-np.nanmin(channelflux))
			ax[1,0].set_ylim([np.nanmin(channelflux)-tenPspread ,np.nanmax(channelflux)+tenPspread])
			ax[1,0].set_ylabel('Discharge \n (m$^3$s$^{-1}$)',rotation=0,ha='right',va='center')
			MidLegend=plt.legend(handles=[recharge,pressure],bbox_to_anchor=(hleft+1.45*plotwidth+hspace,vbot+plotheight,0.8*plotwidth,1.3*plotheight),
													 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.2,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,fontsize=9,
													 borderaxespad=0,mode='expand', frameon=False)
		else:
			pressure=[]
			for b,bandname in enumerate(['Low','Mid','Top']):
				recharge=ax[1,1].plot(time, MeanForce, linestyle='-',linewidth=1.5,color='k',label='Total recharge')[0]
				ax[1,1].set_ylabel('Recharge \n (m$^3$s$^{-1}$)',rotation=0,ha='right',va='center')
				pressure.append(ax[1,0].plot(time, channelflux[:,bandloc[b]], linestyle='-',linewidth=1.5,color=BandProp[bandname][2],label='Eff. discharge for '+BandProp[bandname][4])[0])
			tenPspread=0.1*(np.nanmax(channelflux[:,bandloc[0]])-np.nanmin(channelflux[:,bandloc[2]]))
			ax[1,0].set_ylim([np.nanmin(channelflux[:,bandloc[2]])-tenPspread ,np.nanmax(channelflux[:,bandloc[0]])+tenPspread])
			ax[1,0].set_ylabel('Discharge \n (m$^3$s$^{-1}$)',rotation=0,ha='right',va='center')
			MidLegend=plt.legend(handles=[recharge,pressure[0],pressure[1],pressure[2]],
													 bbox_to_anchor=(hleft+1.45*plotwidth+hspace,vbot+plotheight+0.8*vspace,0.55*plotwidth,1.3*plotheight),
													 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.6,columnspacing=0.2,
													 borderpad=0.3,handletextpad=0.2,fontsize=9,borderaxespad=0,mode='expand', frameon=False)

		#===N Map
		if expe in['D','F']:
			#get a set of dates
			dates=[min(np.where(MeanForce>1.6e-1)[0]),#first melt day
						 max(np.where(MeanForce>1.6e-1)[0]),#last melt day
						 min(np.where(Mean==np.nanmin(Mean))[0]),#minimum effective pressure
						 min(np.where(MeanForce==np.nanmax(MeanForce))[0])]#maximum water input
			datedict={str(dates[0]):'First melt, day {}'.format(str(dates[0])),
								str(dates[1]):'Last melt, day {}'.format(str(dates[1])),
								str(dates[2]):'Minimum $N$, day {}'.format(str(dates[2])),
								str(dates[3]):'Maximum recharge, day {}'.format(str(dates[3]))}
		elif expe=='C':
			#get a set of dates
			minInput=max(np.where(MeanForce==np.nanmin(MeanForce))[0])
			maxInput=min(np.where(MeanForce==np.nanmax(MeanForce))[0])
			dates=[minInput,#first melt day
						 maxInput,#maximum water input
						 min(np.where(Mean==np.nanmin(Mean))[0]),#minimum effective pressure
						 int(minInput+0.5*(maxInput-minInput))]#mid input period
			datedict={str(dates[0]):'Minimum recharge, {} hours'.format(str(dates[0])),
								str(dates[1]):'Maximum recharge, {} hours'.format(str(dates[1])),
								str(dates[2]):'Minimum $N$, {} hours'.format(str(dates[2])),
								str(dates[3]):'Median recharge, {} hours'.format(str(dates[3]))}
		dates=np.sort(dates)
		Fluxlims=[np.nanmin(EffFlux[dates[:],:]),np.nanmax(EffFlux[dates[:],:])]
		#Fluxlims=[np.nanmin(EffFlux),np.nanmax(EffFlux)]
		fluxnorm=mpl.colors.SymLogNorm(linthresh=0.1*np.diff(Fluxlims), linscale=1,vmin=Fluxlims[0], vmax=Fluxlims[1])
		#fluxnorm=mpl.colors.Normalize(vmin=Fluxlims[0], vmax=Fluxlims[1])
		fluxcmap = plt.cm.viridis
		lims=[np.nanmin(data),np.nanmax(data)]
		#norm,cmap=CenteredSymLogNorm(plt.cm.bwr,TreshRatio=0.2, linlenght=1,minval=lims[0], maxval=lims[1])
		norm=mpl.colors.Normalize(vmin=lims[0], vmax=lims[1])
		cmap=shiftedColorMap(plt.cm.bwr, start=0, midpoint= 1-lims[1]/(lims[1]-lims[0]), stop=1.0, name='shiftedcmap')

		#loop on the dates
		for m in np.arange(0,4):
			row=2-2*(m/2)
			col=(np.mod(m,2))
			ax[row,col].set_position([hleft+col*(hspace+plotwidth),
																vbot+row*(plotheight+vspace),
																plotwidth,plotheight])

			# get the limit and plot vertical date line
			centerbotYlims=ax[1,0].get_ylim()
			ax[1,0].plot((time[dates[m]],time[dates[m]]), centerbotYlims, linestyle='-',linewidth=1.5,color=(0.7,0.7,0.7))
			centertopYlims=ax[1,1].get_ylim()
			ax[1,1].plot((time[dates[m]],time[dates[m]]), centertopYlims, linestyle='-',linewidth=1.5,color=(0.7,0.7,0.7))
			#===Channel discharge
			if effcoordtype=='coords_ch' or name=='rh':
				#=== We have a channel overlay, first plot effective pressure
				#===defining colormap
				if np.all(coords[1,:]==midcoord):
					if name in ['mh1','cdf']:
						timedata=np.flipud(data[dates[m],:])
					else:
						timedata=(data[dates[m],:])
						mapim=ax[row,col].imshow(np.asarray([timedata]), extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
																		 origin='lower',norm=norm,cmap=cmap,aspect='auto')
				else:
					griddeddata,gridpoint_x,gridpoint_y=Griding(data[dates[m],:],coords)
					mapim=ax[row,col].imshow(griddeddata.T, extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
																	 origin='lower',norm=norm,cmap=cmap)
				#===and then go on with channels
				indices=[]
				legtext=[]
				indices.append(np.where(np.logical_and(EffFlux[dates[m],:]>=1.0e-3,EffFlux[dates[m],:]<=1.0))[0])
				legtext.append(r'From  $1.0 \times 10^{-3}$ to $1.0$')
				indices.append(np.where(EffFlux[dates[m],:]>1.0)[0])
				legtext.append(r'More than $1.0$')
				#put first the legend and first point
				for index in [0,1]:
					if len(indices[index])>0:
						if scaterplot:
							ax[row,col].scatter(fluxcoords[0,indices[index][0]]*1.0e3,fluxcoords[1,indices[index][0]]*1.0e3,
																	color='black',marker='s',linewidths=0,s=0.5+index,label=legtext[index])
						else:
							ax[row,col].plot(x[:,indices[index][0]],y[:,indices[index][0]],
															 color='black',linewidth=0.5+index,solid_capstyle='round',
															 label=legtext[index])

					#and then go on without legend
					if len(indices[index])>1:
						if scaterplot:
							ax[row,col].scatter(fluxcoords[0,indices[index][1:]]*1.0e-3,fluxcoords[1,indices[1:]]*1.0e-3,
																	color='black',marker='s',linewidth=0,s=0.5+index)
						else:
							ax[row,col].plot(x[:,indices[index][1:]],y[:,indices[index][1:]],
															 color='black',linewidth=0.5+index,solid_capstyle='round')
					if dates[m]==dates[2]:
						l=ax[row,col].legend(bbox_to_anchor=(hleft,0,plotwidth,0.5*vbot),
																 bbox_transform=layout.transFigure,frameon=False,
																 mode='expand',fontsize=9,ncol=4,
																 title='Efficient Discharge (m$^3$s$^{-1}$)',
																 borderaxespad=0,labelspacing=0.2,columnspacing=0.2,
																 borderpad=0.3,handletextpad=0.2)

			#===Now dealing with maps
			else:
				if np.all(fluxcoords[1,:]==0.) or np.all(fluxcoords[1,:]==10*1.0e3):
					mapim=ax[row,col].imshow([EffFlux[dates[m],:]], extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
																		 origin='lower',norm=fluxnorm,cmap=fluxcmap)
				else:
					griddedflux,fluxpoint_x,fluxpoint_y=Griding(EffFlux[dates[m],:],fluxcoords)
					mapim=ax[row,col].imshow(abs(griddedflux).T, extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
																		 origin='lower',norm=fluxnorm,cmap=fluxcmap)

			#clip and outline Bench like
			if expe=='F':
				patch = Polygon(boundaries,transform=ax[row,col].transData)
				mapim.set_clip_path(patch)
				ax[row,col].plot(BenchLines*1.0e-3,0.5*BenchWidth*1.0e-3,color=(0.7,0.7,0.7))
				ax[row,col].plot(BenchLines*1.0e-3,-0.5*BenchWidth*1.0e-3,color=(0.7,0.7,0.7))
				if name!='db':
					lowband=mpl.patches.Rectangle([0.6, -0.55], 0.3, 1.1, facecolor='b',alpha=1,
																				edgecolor='b', transform=ax[row,col].transData)
					midband=mpl.patches.Rectangle([3, -0.55], 0.3, 1.1, facecolor='c',alpha=1,
																				edgecolor='c', transform=ax[row,col].transData)
					topband=mpl.patches.Rectangle([5.1, -0.55], 0.3, 1.1, facecolor='r',alpha=1,
																				edgecolor='r', transform=ax[row,col].transData)
			else:
				if name!='db':
					lowband=mpl.patches.Rectangle([10, -1], 5, 22, facecolor='b',alpha=1,
																				edgecolor='b', transform=ax[row,col].transData)
					midband=mpl.patches.Rectangle([50, -1], 5, 22, facecolor='c',alpha=1,
																				edgecolor='c', transform=ax[row,col].transData)
					topband=mpl.patches.Rectangle([85, -1], 5, 22, facecolor='r',alpha=1,
																				edgecolor='r', transform=ax[row,col].transData)
			if name!='db':
				lowband.set_zorder(0)
				midband.set_zorder(0)
				topband.set_zorder(0)
				layout.patches.append(lowband)
				layout.patches.append(midband)
				layout.patches.append(topband)

			#===Create the arrow
			# 1. Get transformation operators for axis and figure
			axtimetr = ax[1,0].transData # transformation from time axis
			axuptimetr = ax[1,1].transData # transformation from top time axis
			axmaptr = ax[row,col].transData # transformation from map axis
			fluxtr = layout.transFigure.inverted() # Display -> Figure
			# 2. Transform arrow point from axis to figure coordinates
			if row==2:
				startpoint = fluxtr.transform(axuptimetr.transform((time[dates[m]], centertopYlims[1])))
				if col==1:
					labtext=ax[1,1].text(1.,-0.2,datedict[str(dates[m])],fontsize=9,va='center',ha='right',transform=ax[row,col].transAxes,zorder=5)
					endpoint = fluxtr.transform(axmaptr.transform((0.3*Xlims[1], Ylims[0])))
				else:
					labtext=ax[1,1].text(0.,-0.2,datedict[str(dates[m])],fontsize=9,va='center',ha='left',transform=ax[row,col].transAxes,zorder=5)
					endpoint = fluxtr.transform(axmaptr.transform((0.7*Xlims[1], Ylims[0])))
			else:
				startpoint = fluxtr.transform(axtimetr.transform((time[dates[m]], centerbotYlims[0])))
				if col==1:
					labtext=ax[1,1].text(1.,1.2,datedict[str(dates[m])],fontsize=9,va='center',ha='right',transform=ax[row,col].transAxes,zorder=5)
					endpoint = fluxtr.transform(axmaptr.transform((0.3*Xlims[1], Ylims[1])))
				else:
					labtext=ax[1,1].text(0.,1.2,datedict[str(dates[m])],fontsize=9,va='center',ha='left',transform=ax[row,col].transAxes,zorder=5)
					endpoint = fluxtr.transform(axmaptr.transform((0.7*Xlims[1], Ylims[1])))
			# 4. Create the patch
			fluxarrow = mpl.patches.FancyArrowPatch(startpoint, endpoint,
																							transform=layout.transFigure,
																							fc = (0.7,0.7,0.7), ec=(0.7,0.7,0.7),connectionstyle="arc3,rad=0.0",#"arc3,rad=0.2",
																							arrowstyle='simple,tail_width=0.1,head_length=0.7', alpha = 1,shrinkA=0,
																							mutation_scale = 10.)

			# 5. Add patch to list of objects to draw onto the figure
			layout.patches.append(fluxarrow)
			fluxarrow.set_zorder(1)
			ax[row,col].set_zorder(0)

		#===colorbar
		if effcoordtype=='coords_ch' or name=='rh':
			plt.gca().add_artist(MidLegend)
			barPax = layout.add_axes([hleft+hspace+plotwidth,0.25*vbot,plotwidth,0.2*plotheight])
			Pbartext='Effective pressure'
			Pbarunit='MPa'
			offseter=0
			binnumber=5
		else:
			barPax = layout.add_axes([hleft,0.25*vbot,2.*plotwidth+hspace,0.2*plotheight])
			Pbartext='Efficient discharge'
			Pbarunit='m$^3$s$^{-1}$'
			binnumber=9
			if np.nanmean(EffFlux)!=0.:
				offseter=np.int(np.log10(np.nanmax(EffFlux)))-1
			else:
				offseter=0

		fmt=OOMFormatter(offseter, "%1.1f")
		if effcoordtype=='coords_ch' or name=='rh':
			Pbar = mpl.colorbar.ColorbarBase(barPax, norm=norm,cmap=cmap, orientation='horizontal',format=fmt)
		else:
			Pbar = mpl.colorbar.ColorbarBase(barPax, norm=fluxnorm,cmap=fluxcmap, orientation='horizontal',format=fmt)
		offset = Pbar.ax.get_xaxis().get_offset_text()
		Pbar.locator = ticker.MaxNLocator(nbins=binnumber)
		Pbar.update_ticks()
		if offseter!=0:
			offtext=r'$\times 10^{{{}}}$'.format(str(offseter))
			barPtitle=Pbar.set_label(Pbartext+' ('+offtext+' '+Pbarunit+')', labelpad=-22,ha='center',va='bottom')
		else:
			barPtitle=Pbar.set_label(Pbartext+' ('+Pbarunit+')', labelpad=-22,ha='center',va='bottom')
		offset.set_visible(False)
		Pbar.solids.set_edgecolor("face")

		layout.savefig(figdir+'/'+name+'/TransientFlux_'+name+'_'+tag+'.pdf',format='pdf')
	plt.show
	return 0

# }}}

# {{{ Steady(runnum)

def Steady(runnum):
	tag=expe+str(runnum)

	width=17.8/2.54 #two columns JoG
	plotrow=2
	plotcol=2

	if expe in ['A','B']:
		ratio=1.6
		#0.2 is the width/heigh ratio of the maps (without borders)
		height=0.2*width*ratio
		Xlims=[0,100]
		Xticks=[10,30,50,70,90]
		Xlabels=['$10$','$30$','$50$','$70$','$90$']
		Ylims=[0,20]
		Yticks=[5,10,15]
		Ylabels=['$5$','$10$','$15$']
		midcoord=10.0e3

	if expe=='E':
		ratio=1.7
		#1/6 is the width/heigh ratio of the maps (without borders)
		height=(1./6.)*width*ratio
		Xlims=[0,6]
		Xticks=[1,2,3,4,5]
		Xlabels=['$1$','$2$','$3$','$4$','$5$']
		Ylims=[-0.5,0.5]
		Yticks=[-0.3,0,0.3]
		Ylabels=['$-0.3$','$0$','$0.3$']
		midcoord=0.0

		BenchLength=6.0e3
		BenchStep=60.
		BenchLines=np.arange(0,BenchLength+BenchStep,BenchStep)
		s2 = 100.0/BenchLength
		s_xend = 100.0*(BenchLength+200)**0.25 + s2*BenchLength - 100.0*200**0.25
		BenchSurf = 100.0*(BenchLines+200)**0.25 + s2*BenchLines - 100.0*200**0.25
		f10 = (s_xend - 0.05*BenchLength)/BenchLength**2
		f0 = f10*BenchLines**2 + 0.05*BenchLines
		h0 = -4.5*(BenchLines/BenchLength) + 5
		BenchWidth = 2.0*(((BenchSurf-f0)/h0)/0.5e-6)**(1/3.)
		upstream=np.zeros(np.shape(BenchWidth))
		for w,val in enumerate(BenchWidth):
			upstream[w]=np.sum(BenchWidth[:-w])
		upstream[0]=np.sum(BenchWidth)

	#horizontal buffer and related spaces
	hbuff=0.2
	hleft=0.4*hbuff
	hright=0.5*hbuff
	hspace=0.1*hbuff
	#vertical buffer and related spaces
	vbuff=1.-((1.-hbuff)/ratio)
	vbot=0.45*vbuff
	vtop=0.45*vbuff
	vspace=0.1*vbuff
	#subplot sizes
	plotwidth=(1.-hleft-hright-hspace)/plotcol
	plotheight=(1.0-vbot-vtop-vspace)/plotrow

	Influx={'A1':7.93e-11,
					'A2':1.59e-09,
					'A3':5.79e-09,
					'A4':2.50e-08,
					'A5':4.50e-08,
					'A6':5.79e-07,
					'E1':1.158e-6,
					'E2':1.158e-6,
					'E3':1.158e-6,
					'E4':1.158e-6,
					'E5':1.158e-6,}
	if expe=='B':
		#gathering mouins input (once and for all)
		inputval=np.zeros((101))
		moulinfile='/home/bfl022/Model/GHIP/hydro_intercomparison/input_functions/source/'+tag+'_M.csv'
		with open(moulinfile) as csvfile:
			MoulinPos=csv.reader(csvfile,delimiter=',')
			for row in MoulinPos:
				Xs =int(row[1])
				inputval[:Xs/1000]+=float(row[3])

	for name in names:
		NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
		print NCFile
		#===Loading data
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			data,coords,datacoordtype=GetVar(DatFile,DataName,name,False)
			if name in['mh2','mw','mw_prime','og','og_prime','mh1','cdf','rh','bf','id']:
				EffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'Q',name,False)
				EffFlux=np.abs(EffFlux)
				if effcoordtype=='coords_ch' or name=='rh':
					try:
						elts = DatFile.variables['connect_ch'][:,:].T.astype(int)
						if name=='id':
							fluxcoords=DatFile.variables['coords1'][:][0]
							x    = np.asarray([fluxcoords[elts[:,0]]*1.0e-3,fluxcoords[elts[:,1]]*1.0e-3])
							if expe in ['E','F']:
								y    = np.zeros(np.shape(x))
							else:
								y    = 10*np.ones(np.shape(x))
						else:
							if name=='mh1':
								fluxcoords=DatFile.variables['coords3'][:,:]
							else:
								fluxcoords=DatFile.variables['coords1'][:,:]
							x    = np.asarray([fluxcoords[0,elts[:,0]]*1.0e-3,fluxcoords[0,elts[:,1]]*1.0e-3])
							y    = np.asarray([fluxcoords[1,elts[:,0]]*1.0e-3,fluxcoords[1,elts[:,1]]*1.0e-3])
						scaterplot=False
					except KeyError:
						scaterplot=True
			elif name=='as':
				EffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'Dc',name,False)
			elif name=='db':
				EffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'R',name,False)
			elif name=='sb':
				EffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'T',name,False)
			elif name=='sb_old':
				EffFlux,fluxcoords,effcoordtype=GetVar(DatFile,'K',name,False)
			DatFile.close()
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))
			continue
		if name=='jsb' and expe=='E':
			tmp0= 100.0*(coords[0,:]+200)**0.25 + s2*coords[0,:] - 100.0*200**0.25
			tmp1 = f10*coords[0,:]**2 + 0.05*coords[0,:]
			tmp2 = -4.5*(coords[0,:]/BenchLength) + 5
			dataWidth = 2.0*(((tmp0-tmp1)/tmp2)/0.5e-6)**(1/3.)
			data[np.where(abs(coords[1,:])>0.55*dataWidth)]=np.nan
		if name=='rh' and tag=='E4':
			EffFlux[np.where(abs(fluxcoords[0,:])>740)]=np.nan
			data[np.where(abs(coords[0,:])>740)]=np.nan

		if expe=='E':
			boundaries=np.vstack((np.hstack((BenchLines,np.flipud(BenchLines)))*1.0e-3,
														np.hstack((BenchWidth,np.flipud(-BenchWidth)))*1.0e-3*0.5)).T

		#===layout and some axis stuff
		layout,ax=plt.subplots(plotrow,plotcol,figsize=(width,height))
		for i in [0,1]:
			for j in [0,1]:
				ax[i,j].set_position([hleft+i*(hspace+plotwidth),
															vbot+j*(plotheight+vspace),
															plotwidth,plotheight])
				ax[i,j].set_xlim(Xlims)
				ax[i,j].set_xticks(Xticks)
				ax[i,0].xaxis.set_ticklabels(Xlabels)
				ax[i,1].xaxis.set_ticklabels([])

				ax[i,1].set_ylim(Ylims)
				ax[i,1].set_yticks(Yticks)
				ax[i,1].yaxis.set_ticklabels(Ylabels)
				ax[1,j].yaxis.tick_right()
				ax[1,j].yaxis.set_ticks_position('both')
				ax[1,j].yaxis.set_label_position('right')

				ax[i,0].set_xlabel('$x$ (km)',labelpad=2)
				ax[i,1].set_ylabel('$y$ (km)',labelpad=2)
		ax[0,0].set_ylabel('$y$ (km)',labelpad=2)
		ax[0,0].set_ylabel('$N$ (MPa)',labelpad=2)
		ax[0,0].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))

		ax[1,0].set_ylabel('Discharge \n (m$^3$s$^{-1}$)',labelpad=2)
		ax[1,0].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))

		# #===defining colormap
		# cmap = plt.cm.viridis
		# if name=='db':
		# 	cmap = plt.cm.get_cmap('viridis', 1)

		# #===N Map
		# lims=[np.nanmin(data),np.nanmax(data)]
		# if name=='db':
		# 	lims=[0.9*np.nanmean(data),1.1*np.nanmean(data)]
		# norm=mpl.colors.Normalize(lims[0], lims[1])
		if name=='db':
			cmap = plt.cm.get_cmap('bwr', 1)
			lims=[0.9*np.nanmean(data),1.1*np.nanmean(data)]
			norm=mpl.colors.Normalize(lims[0], lims[1])
		else:
			lims=[np.nanmin(data),np.nanmax(data)]
			norm=mpl.colors.Normalize(vmin=lims[0], vmax=lims[1])
			cmap=shiftedColorMap(plt.cm.bwr, start=0, midpoint= 1-lims[1]/(lims[1]-lims[0]), stop=1.0, name='shiftedcmap')
			#norm,cmap=CenteredSymLogNorm(plt.cm.bwr,TreshRatio=0.2, linlenght=1,minval=lims[0], maxval=lims[1])

		if np.all(coords[1,:]==midcoord):
			if name in ['mh1','cdf']:
				data=np.flipud(data)
				coords=np.fliplr(coords)
			if len(np.shape(data))>1:
				data=np.nanmean(data,axis=0)
			mapim=ax[0,1].imshow(np.asarray([data]), extent=(Xlims[0],np.nanmax(coords[0])*1.0e-3,Ylims[0],Ylims[1]),
													 origin='lower',norm=norm,cmap=cmap,aspect='auto')
		else:
			griddeddata,gridpoint_x,gridpoint_y=Griding(data,coords)
			mapim=ax[0,1].imshow(griddeddata.T, extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
													 origin='lower',norm=norm,cmap=cmap)
			meandata=np.nanmean(griddeddata,axis=1)
			maxdata=np.nanmax(griddeddata,axis=1)
			mindata=np.nanmin(griddeddata,axis=1)

		barNax = layout.add_axes([hleft,1-0.6*vtop,plotwidth,0.2*plotheight])
		barNtitle=barNax.set_title('Effective pressure (MPa)',fontsize=9)
		if name=='db':
			Nbar = mpl.colorbar.ColorbarBase(barNax, norm=norm,cmap=cmap, orientation='horizontal',ticks=[np.nanmean(data)],format=FormatScalarFormatter('%.2f'))
			Nbar.ax.xaxis.set_ticklabels(str(np.nanmean(data)))
		else:
			Nbar = mpl.colorbar.ColorbarBase(barNax, norm=norm,cmap=cmap, orientation='horizontal')
			Nbar.locator = ticker.MaxNLocator(nbins=9)
		Nbar.update_ticks()
		Nbar.solids.set_edgecolor("face")

		#===N Section
		if np.all(coords[1,:]==midcoord):#for 1D
			Npress=ax[0,0].plot(coords[0,:]*1.0e-3, data, linestyle='-',linewidth=1.5,color='b',label='Effective pressure')[0]
			ax[0,0].legend(handles=[Npress],bbox_to_anchor=(hleft,0.,0.2,0.07),bbox_transform=layout.transFigure,ncol=1,
										 labelspacing=0.,columnspacing=0.,borderpad=0.,handletextpad=0.2,fontsize=9,handler_map={},
										 borderaxespad=0,mode='expand', frameon=False)
			if name=='db':
				ax[0,0].set_ylim(lims)
				ax[0,0].set_yticks([np.nanmean(data)])
				ax[0,0].yaxis.set_ticklabels([str(np.nanmean(data))])
				ax[0,0].yaxis.set_major_formatter(FormatScalarFormatter('%.2f'))

		else:
			Npress=ax[0,0].fill_between(gridpoint_x[:,0], maxdata,mindata, facecolor='b',edgecolor='b',alpha=0.3,label='Effective pressure')
			ax[0,0].plot(gridpoint_x[:,0], meandata, linestyle='-',linewidth=1.5,label=namestyles[name][0],color='b')
			ax[0,0].legend(handles=[Npress],
										 bbox_to_anchor=(hleft,0.,0.2,0.07),
										 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.,columnspacing=0.,borderpad=0.,handletextpad=0.2,fontsize=9,
										 borderaxespad=0,handler_map={Npress:FillandLineHandler('b')},mode='expand', frameon=False)

		#===Flux ratio Map
		if name in['mh2','mw','mw_prime','og','og_prime','mh1','cdf','rh','bf','id']:
			Fluxlims=[np.nanmin(EffFlux),np.nanmax(EffFlux)]
			if name=='cdf':
				EffFlux=np.flipud(EffFlux)
				coords=np.fliplr(fluxcoords)
			if effcoordtype=='coords_ch' or name=='rh':
				#===and then go on with channels
				indices=[]
				legtext=[]
				indices.append(np.where(np.logical_and(EffFlux>=1.0e-3,EffFlux<=1.0))[0])
				legtext.append(r'From  $1.0 \times 10^{-3}$ to $1.0$')
				indices.append(np.where(EffFlux>1.0)[0])
				legtext.append(r'More than $1.0$')
				#put first the legend and first point
				for index in [0,1]:
					if len(indices[index])>0:
						if scaterplot:
							ax[1,1].scatter(fluxcoords[0,indices[index][0]]*1.0e3,fluxcoords[1,indices[index][0]]*1.0e3,
															color='black',marker='s',linewidths=0,s=0.5+index,label=legtext[index])
						else:
							ax[1,1].plot(x[:,indices[index][0]],y[:,indices[index][0]],
													 color='black',linewidth=0.5+index,solid_capstyle='round',
													 label=legtext[index])

					#and then go on without legend
					if len(indices[index])>1:
						if scaterplot:
							ax[1,1].scatter(fluxcoords[0,indices[index][1:]]*1.0e-3,fluxcoords[1,indices[index][1:]]*1.0e-3,
															color='blue',marker='s',linewidth=0,s=0.5+index)
						else:
							ax[1,1].plot(x[:,indices[index][1:]],y[:,indices[index][1:]],
													 color='blue',linewidth=0.5+index,solid_capstyle='round')
				l=ax[1,1].legend(bbox_to_anchor=(hleft+hspace+plotwidth,1-vtop,plotwidth,vtop),
												 bbox_transform=layout.transFigure,frameon=False,
												 mode='expand',fontsize=9,ncol=4,
												 title='Efficient Discharge (m$^3$s$^{-1}$)',
												 borderaxespad=0,labelspacing=0.2,columnspacing=0.2,
												 borderpad=0.3,handletextpad=0.2)

			#Now dealing with maps
			else:
				fluxnorm=mpl.colors.SymLogNorm(linthresh=0.1*np.diff(Fluxlims), linscale=1,vmin=Fluxlims[0], vmax=Fluxlims[1])
				if np.all(fluxcoords[1,:]==0.) or np.all(fluxcoords[1,:]==10*1.0e3):
					fluxmap=ax[1,1].imshow([EffFlux], extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
																 origin='lower',norm=fluxnorm,cmap=mpl.cm.viridis)
				else:
					griddedflux,fluxpoint_x,fluxpoint_y=Griding(EffFlux,fluxcoords)
					fluxmap=ax[1,1].imshow(abs(griddedflux).T, extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
												origin='lower',norm=fluxnorm,cmap=mpl.cm.viridis)

				barPax = layout.add_axes([hleft+hspace+plotwidth,1-0.6*vtop,plotwidth,0.2*plotheight])
				if np.nanmean(EffFlux)!=0.:
					offseter=np.int(np.log10(np.nanmax(EffFlux)))-1
				else:
					offseter=0
				if name=='rh' and tag in ['A4','A5','E3','E2']:
					fmt=OOMFormatter(offseter, "%1.f")
				elif name=='cdf' and tag in ['A5']:
					fmt=OOMFormatter(offseter, "%1.f")
				else:
					fmt=OOMFormatter(offseter, "%1.1f")
				Pbar = mpl.colorbar.ColorbarBase(barPax, norm=fluxnorm,cmap=mpl.cm.viridis, orientation='horizontal',format=fmt)
				offset = Pbar.ax.get_xaxis().get_offset_text()
				if offseter!=0:
					offtext=r'$\times 10^{{{}}}$'.format(str(offseter))
					barPtitle=barPax.set_title('Efficient Discharge ('+offtext+' m$^3$s$^{-1}$)',fontsize=9)
				else:
					barPtitle=barPax.set_title('Efficient Discharge (m$^3$s$^{-1}$)',fontsize=9)
				offset.set_visible(False)
				Pbar.locator = ticker.MaxNLocator(nbins=5)
				Pbar.update_ticks()
				Pbar.solids.set_edgecolor("face")

		elif name in ['as','db','sb','sb_old']:
			if name=='db':
				Fluxlims=[0.9*np.nanmean(EffFlux),1.1*np.nanmean(EffFlux)]
			else:
				Fluxlims=[np.nanmin(EffFlux),np.nanmax(EffFlux)]
			fluxnorm=mpl.colors.Normalize(Fluxlims[0],Fluxlims[1])
			if np.all(fluxcoords[1,:]==0.) or np.all(fluxcoords[1,:]==10*1.0e3):
				fluxmap=ax[1,1].imshow([EffFlux], extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
											 origin='lower',norm=fluxnorm,cmap=mpl.cm.viridis)
			else:
				griddedflux,fluxpoint_x,fluxpoint_y=Griding(EffFlux,fluxcoords)
				fluxmap=ax[1,1].imshow(abs(griddedflux).T, extent=(Xlims[0],Xlims[1],Ylims[0],Ylims[1]),
											 origin='lower',norm=fluxnorm,cmap=mpl.cm.viridis)

			barPax = layout.add_axes([hleft+hspace+plotwidth,1-0.6*vtop,plotwidth,0.2*plotheight])
			if np.nanmean(EffFlux)!=0.:
				offseter=np.int(np.log10(np.nanmax(EffFlux)))-1
			else:
				offseter=0
			fmt=OOMFormatter(offseter, "%1.1f")
			if name=='db':
				Pbar = mpl.colorbar.ColorbarBase(barPax, norm=fluxnorm,cmap=mpl.cm.viridis, orientation='horizontal',ticks=[np.nanmean(EffFlux)])
				if np.nanmean(EffFlux)!=1:
					tickstring=r'${}\times 10^{{{}{}}}$'.format('{:1.2e}'.format(np.nanmean(EffFlux))[0:4],'{:1.2e}'.format(np.nanmean(EffFlux))[5],'{:1.2e}'.format(np.nanmean(EffFlux))[-1:])
				else:
					tickstring='$1.0$'
				Pbar.ax.xaxis.set_ticklabels([tickstring])
				barPtitle=barPax.set_title('Chanelisation ratio',fontsize=9)
			else:
				Pbar = mpl.colorbar.ColorbarBase(barPax, norm=fluxnorm,cmap=mpl.cm.viridis, orientation='horizontal',format=fmt)
				offset = Pbar.ax.get_xaxis().get_offset_text()
				if offseter!=0:
					offtext=r'$\times 10^{{{}}}$'.format(str(offseter))
					barPtitle=barPax.set_title('Chanelisation ratio ('+offtext+')',fontsize=9)
				else:
					barPtitle=barPax.set_title('Chanelisation ratio',fontsize=9)
				offset.set_visible(False)
				Pbar.locator = ticker.MaxNLocator(nbins=8)
				Pbar.update_ticks()
			Pbar.solids.set_edgecolor("face")

		#===Flux Section
		channelflux,sheetflux,lines=GetFluxRatio(name,tag,[])
		if name=='rh' and tag=='E4':
			channelflux[np.where(lines>740)]=np.nan

		if name=='db':
			ax[1,0].plot(lines*1.0e-3,sheetflux[0]*np.ones(np.shape(lines)),'c',label='Inef. discharge')
			ax[1,0].plot(lines*1.0e-3,channelflux[0]*np.ones(np.shape(lines)),'r',label='Eff. discharge')
			ax[1,0].plot(lines*1.0e-3,(sheetflux+channelflux)[0]*np.ones(np.shape(lines)),color='k',label='Tot. discharge')
			if expe=='A':
				ax[1,0].plot(lines*1.0e-3,100.0e3*np.ones(np.shape(lines))*Influx[tag]*20.0e3,color='y',linestyle='--',label='Total recharge')
				ax[1,0].set_ylim([0.0 ,1.1*np.nanmax(100.0e3*Influx[tag]*20.0e3)])
			elif expe=='B':
				ax[1,0].plot(lines*1.0e-3,inputval,color='y',linestyle='--',label='Total recharge')
				ax[1,0].set_ylim([0.0 ,1.1*np.nanmax(inputval)])
			elif expe=='E':
				ax[1,0].plot(lines*1.0e-3,Influx[tag]*upstream[0]*np.ones(np.shape(lines))*BenchStep,color='y',linestyle='--',label='Total recharge')
				ax[1,0].set_ylim([0.0 ,1.1*np.nanmax(Influx[tag]*upstream*BenchStep)])

		else:
			if np.nansum(sheetflux)>0.0:
				ax[1,0].plot(lines*1.0e-3,sheetflux,'c',label='Inef. discharge')
			if np.nansum(channelflux)>0.0:
				ax[1,0].plot(lines*1.0e-3,channelflux,'r',label='Eff. discharge')
			if np.nansum(sheetflux)>0.0 and np.nansum(channelflux)>0.0:
				ax[1,0].plot(lines*1.0e-3,sheetflux+channelflux,color='k',label='Total discharge')

			if expe=='A':
				ax[1,0].plot(lines*1.0e-3,(100.0e3-lines)*Influx[tag]*20.0e3,color='y',linestyle='--',label='Total recharge')
				ax[1,0].set_ylim([0.0 ,1.1*np.nanmax(100.0e3*Influx[tag]*20.0e3)])
			elif expe=='B':
				ax[1,0].plot(lines*1.0e-3,inputval+(100.0e3-lines)*Influx['A1']*20.0e3,color='y',linestyle='--',label='Total recharge')
				ax[1,0].set_ylim([0.0 ,1.3*np.nanmax(inputval)])
			elif expe=='E':
				ax[1,0].plot(lines*1.0e-3,Influx[tag]*upstream*BenchStep,color='y',linestyle='--',label='Total recharge')
				ax[1,0].set_ylim([0.0 ,1.1*np.nanmax(Influx[tag]*upstream*BenchStep)])

		(lines, labels)= ax[1,0].get_legend_handles_labels()
		leg1=ax[1,0].legend(handles=lines,labels=labels,bbox_to_anchor=(hleft+0.2,0.,len(lines)*(2*plotwidth+hspace-0.2)/4.,0.07),
												bbox_transform=layout.transFigure,ncol=4,labelspacing=0.,columnspacing=0.,borderpad=0.,handletextpad=0.2,fontsize=9,
												borderaxespad=0,mode='expand', frameon=False)
		if expe=='E':
			patch = Polygon(boundaries,transform=ax[0,1].transData)
			mapim.set_clip_path(patch)
			patch = Polygon(boundaries,transform=ax[1,1].transData)
			if effcoordtype!='coords_ch' :
				fluxmap.set_clip_path(patch)
			for i in [0,1]:
				ax[i,1].plot(BenchLines*1.0e-3,0.5*BenchWidth*1.0e-3,color=(0.7,0.7,0.7))
				ax[i,1].plot(BenchLines*1.0e-3,-0.5*BenchWidth*1.0e-3,color=(0.7,0.7,0.7))

		layout.savefig(figdir+'/'+name+'/Steady_'+DataName+'_'+name+'_'+tag+'.pdf',format='pdf')
		layout.savefig(figdir+'/'+name+'/Steady_'+DataName+'_'+name+'_'+tag+'.png',format='png')
		plt.show
	return 0

# }}}

# {{{ SectionsDiff()

def SectionsDiff(frame):

	width=8.6/2.54
	layout,ax=plt.subplots(9,1,sharex=True,figsize=(width,width*2))
	if frame<3:
		hleft=0.1
	else:
		hleft=0.15
	hright=0.02
	width=(1-hleft-hright)
	vside=0.045
	vspace=0.01
	height=0.85/9

	k=-1

	#gathering mouins input (once and for all)
	inputval=np.zeros(101)
	moulinfile='/home/bfl022/Model/GHIP/hydro_intercomparison/input_functions/source/B'+str(frame)+'_M.csv'
	with open(moulinfile) as csvfile:
		MoulinPos=csv.reader(csvfile,delimiter=',')
		for row in MoulinPos:
			Xs =int(row[1])
			inputval[Xs/1000]+=float(row[3])

	for i,name in enumerate(names):
		#checking existence of the two needed runs
		if len(glob(resdir+name+'/B'+str(frame)+'_'+name+'.nc'))>0 and len(glob(resdir+name+'/A5_'+name+'.nc'))>0:
			NCDiff=resdir+name+'/A5_'+name+'.nc'
			DiffFile = Dataset(NCDiff, mode='r')
			diffdata,diffcoords,diffcoordtype=GetVar(DiffFile,DataName,name,False)
			DiffFile.close()
			k+=1
			#Now work on it
			NCFile=resdir+name+'/B'+str(frame)+'_'+name+'.nc'
			try:
				DatFile	 = Dataset(NCFile, mode='r')
				data,coords,datacoordtype=GetVar(DatFile,DataName,name,False)
				DatFile.close()
			except (RuntimeError,IOError):
				print ('File {} does not exist'.format(NCFile))
				continue

			griddeddata,gridpoint_x,gridpoint_y=Griding(data,coords)
			griddeddiff,gridpoint_x,gridpoint_y=Griding(diffdata,diffcoords)
			floatation=1.0e-6*910*9.81*(1.0+6.*pow(gridpoint_x*1.0e3+5.0e3,0.5)-6.0*pow(5.0e3,0.5))

			#			normdiff=100.*((griddeddata-griddeddiff)/floatation)
			normdiff=(griddeddata-griddeddiff)

			meannormdiff=np.nanmean(normdiff,axis=1)
			maxnormdiff=np.nanmax(normdiff,axis=1)
			minnormdiff=np.nanmin(normdiff,axis=1)

			ax[k].set_position([hleft,vside+(8-k)*(height+vspace),width,height])
			#axis stuff
			ax[k].set_xlim([0 ,100])
			ax[k].set_xticks([10,30,50,70,90])
			ax[k].xaxis.set_ticklabels([])

			if k==0:
				ax[k].set_title('B'+str(frame),fontsize=9,position=(0.5,1.))
			if k==8:
				ax[k].xaxis.set_ticklabels(['$10$','$30$','$50$','$70$','$90$'],fontsize=9)
				ax[k].tick_params(axis='x', which='major', pad=2)
				ax[k].set_xlabel('$x$ (km)',labelpad=0)


			ax2=ax[k].twinx()
			ax2.set_position(ax[k].get_position())
			ax2.bar(gridpoint_x[:,0],inputval[:],align='center',label='input',color='r',log=True,edgecolor='none')

			ax2.set_xlim([0 ,100])
			ax2.set_ylim([0.5 ,120])
			ax2.yaxis.set_ticklabels([])
			ax2.yaxis.set_ticks_position('none')
			if frame<4:
				ax2.set_title(namestyles[name][0],fontsize=9,
											va='top',ha='right',
											position=(0.95,0.2))
			else:
				ax2.set_title(namestyles[name][0],fontsize=9,
											va='top',ha='left',
											position=(0.9,0.85))

			ax[k].fill_between(gridpoint_x[:,0], maxnormdiff,minnormdiff, facecolor='b',edgecolor='b',alpha=0.3)
			ax[k].plot(gridpoint_x[:,0], meannormdiff, linestyle='-',linewidth=1.5,label=namestyles[name][0],color='b')

			ax[k].set_ylim([np.nanmin(minnormdiff[5:]) ,np.nanmax(maxnormdiff[5:])])
			ax[k].set_zorder(ax[k].get_zorder()+1)
			ax[k].patch.set_visible(False)
			ax[k].yaxis.set_ticks_position('both')

			if frame>2:
				ax[k].set_ylim([-0.5,0.6])
				ax[k].set_yticks([-0.4,-0.2,0,0.2,0.4])
				ax[k].yaxis.set_ticklabels(['$-0.4$','$-0.2$','$0$','$0.2$','$0.4$'],fontsize=9)
			else:
				ax[k].set_ylim([-1,4])
				ax[k].set_yticks([-1,0,1,2,3])
				ax[k].yaxis.set_ticklabels(['$-1$','$0$','$1$','$2$','$3$'],fontsize=9)


			ax[k].tick_params(axis='y', which='major', pad=2)

			ax[k].text(0.05,0.7,string.lowercase[k],fontsize=9,va='bottom',ha='left',transform=ax[k].transAxes)
	ax[0].text(0.1*hleft,vside+height*((k+1)/2.)+((k)/2.)*vspace,'Effective pressure difference (MPa)',fontsize=9,va='center',ha='left',rotation=90,transform=layout.transFigure)

	plt.show
	return layout

# }}}

# {{{ DiurnalComp(dummy)

def DiurnalComp(dummy):

	for i,name in enumerate(names):
		print name
		try:
			if name=='mw_prime':
				NCComp=resdir+'mw/B5_mw.nc'
			else:
				NCComp=resdir+name+'/B5_'+name+'.nc'
			CompFile = Dataset(NCComp, mode='r')
			compdata,coords,compcoordstype=GetVar(CompFile,DataName,name,False)
			CompFile.close()
			griddedcomp,gridpoint_x,gridpoint_y=Griding(compdata,coords)
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCComp))
			continue
		layout,ax=ExpeLayout()
		for j,tag in enumerate(plotexpe):
			NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
			try:
				DatFile	 = Dataset(NCFile, mode='r')
				try:
					DatFile.variables[DataName].dimensions
					data,coords,datacoordstype=GetVar(DatFile,DataName,name,True)
					DatFile.close()
				except (KeyError):
					print ('File {} does not have {}'.format(NCFile,DataName))
					continue
			except (RuntimeError,IOError):
				print ('File {} does not exist'.format(NCFile))
				continue

			temporalmean=np.nanmean(data,axis=0)
			temporalmax=np.nanmax(data,axis=0)
			temporalmin=np.nanmin(data,axis=0)
			griddeddata,gridpoint_x,gridpoint_y=Griding(temporalmean,coords)
			griddedmax,gridpoint_x,gridpoint_y=Griding(temporalmax,coords)
			griddedmin,gridpoint_x,gridpoint_y=Griding(temporalmin,coords)
			meandata=np.nanmean(griddeddata,axis=1)
			maxdata=np.nanmean(griddedmax,axis=1)
			mindata=np.nanmean(griddedmin,axis=1)
			meancomp=np.nanmean(griddedcomp,axis=1)

			if j==0:
				ax[j].fill_between(gridpoint_x[:,0], maxdata,mindata, facecolor='b',edgecolor='b',alpha=0.3)
				ax[j].plot(gridpoint_x[:,0], meandata, linestyle='-',linewidth=1.5,label=namestyles[name][0],color='b')
				ax[j].plot(gridpoint_x[:,0], meancomp, linestyle='-',linewidth=1.5,label='B5',color='k')
				ax[j].set_xlim([0 ,100])
			else:
				# for t in np.arange(0,np.shape(data)[0]):
				# 	ax[j].scatter(coords[0]*1.0e-3,data[t,:],s=1,marker='.',color='r',alpha=0.3)
				ax[j].fill_between(gridpoint_x[:,0], maxdata,mindata, facecolor='b',edgecolor='b',alpha=0.3)
				ax[j].plot(gridpoint_x[:,0], meandata, linestyle='-',linewidth=1.5,color='b')
				ax[j].plot(gridpoint_x[:,0], meancomp, linestyle='-',linewidth=1.5,color='k')
				ax[j].set_xlim([0 ,100])

		ax[0].legend(bbox_to_anchor=(0.55,1), loc="upper center",  bbox_transform=layout.transFigure,ncol=4,
								 labelspacing=0.2,columnspacing=0.2,borderpad=0.3,handletextpad=0.2)
		layout.savefig(figdir+'/'+name+'/DiurnalComp_'+DataName+'_'+name+'.pdf',format='pdf')
	plt.show()
	return layout

# }}}

# {{{ BandEvol()

def BandEvol(sub):
	if plotexpe[0].startswith('C'):
		width=17.8/2.54
		height=width
		plotrow=11
		plotcol=3
		layout,ax=plt.subplots(plotrow,plotcol,figsize=(width,height))
		Exptag='C'
		baselist=sub
		hleft=0.45/width
		hright=0.45/width
		hspace=0.4/width#0.15
		vside=0.3/height
		vbot=0.9*vside
		vtop=2*vside-vbot
		vspace=0.05/height
		plotwidth=(1-hleft-hright-(hspace*(plotcol-1.)))/(plotcol-0.5)
		plotheight=(1.0-(2*vside)-(vspace*(plotrow-1)))/plotrow

	elif plotexpe[0].startswith('D'):
		width=17.8/2.54
		height=1.*width
		plotrow=11
		plotcol=3
		layout,ax=plt.subplots(plotrow,plotcol,figsize=(width,height))
		Exptag='D'
		baselist=sub
		hleft=0.45/width
		hright=0.3/width
		hspace=0.4/width#0.15
		vside=0.45/height
		vbot=0.55*vside
		vtop=2*vside-vbot
		vspace=0.05/height
		plotwidth=(1-hleft-hright-(hspace*(plotcol-1.)))/(plotcol-0.5)
		plotheight=(1.0-(2*vside)-(vspace*(plotrow-1)))/plotrow

	elif plotexpe[0].startswith('F'):
		width=17.8/2.54
		height=1.1*width
		plotrow=9
		plotcol=3
		layout,ax=plt.subplots(plotrow,plotcol,figsize=(width,height))
		Exptag='F'
		baselist=sub
		hleft=0.42/width
		hright=0.43/width
		hspace=0.4/width#0.15
		vside=0.47/height
		vbot=0.55*vside
		vtop=2*vside-vbot
		vspace=0.05/height
		plotwidth=(1-hleft-hright-(hspace*(plotcol-1.)))/(plotcol-0.5)
		plotheight=(1.0-(2*vside)-(vspace*(plotrow-1)))/plotrow

	k=-1
	IDS=[]
	EDS=[]
	pressure=[]

	AmpList={'C1':0.25,'C2':0.5,'C3':1.0,'C4':2.0,
					 'D1':-4.0,'D2':-2.0,'D3':0.0,'D4':2.0,'D5':4.0,
					 'F1':-6.0,'F2':-3.0,'F3':0.0,'F4':3.0,'F5':6.0}
	BandProp={'Low':[10000.,600.,'b','$Low$','lower band'],
						'Mid':[50000.,3000.,'c','$Mid$','middle band'],
						'Top':[90000.,5100.,'r','$Top$','highest band']}

	if Exptag=='C':
		moulinfile='/home/bfl022/Model/GHIP/hydro_intercomparison/input_functions/source/B5_M.csv'
		MoulinInput=np.zeros((101))
		with open(moulinfile) as csvfile:
			MoulinPos=csv.reader(csvfile,delimiter=',')
			for row in MoulinPos:
				Xs=int(row[1])
				MoulinInput[Xs/1000]+=float(row[3])

	for i,name in enumerate(names):
		print name
		if len(glob(resdir+name+'/'+Exptag+'*_'+name+'.nc'))>0:
			framelist=baselist-1
			k+=1
			#=====Axes for right column left axis===========
			ax[k+1,plotcol-1].set_position([hleft+(plotwidth+hspace)*(plotcol-1),
																		vbot+(plotrow-2-k)*(plotheight+vspace),
																		0.5*plotwidth,plotheight])

			ax[k+1,plotcol-1].yaxis.tick_right()
			ax[k+1,plotcol-1].yaxis.set_ticks_position('both')
			ax[k+1,plotcol-1].yaxis.set_label_position('right')
			ax[k+1,plotcol-1].ticklabel_format(axis='both')
			ax[k+1,plotcol-1].tick_params(axis='y', which='major', pad=2,top=False)
			ax[k+1,plotcol-1].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))

			ax[k+1,plotcol-1].set_xlim([0.8,len(plotexpe)+0.2])
			ax[k+1,plotcol-1].xaxis.set_ticklabels([])
			ax[plotrow-1,plotcol-1].tick_params(axis='x', which='major', pad=2)
			ax[plotrow-1,plotcol-1].set_xlabel('sub-experiment',labelpad=0)

			ax[k+1,plotcol-1].text(0.15,0.8,string.uppercase[k],fontsize=9,va='center',ha='center',transform=ax[k+1,plotcol-1].transAxes)

			#=====Axes for right column right axis===========
			ax2=ax[k+1,plotcol-1].twinx()
			ax2.set_position(ax[k+1,plotcol-1].get_position())
			ax2.yaxis.set_major_locator(MaxNLocator(5,prune='both'))
			ax2.spines['right'].set_color('b')
			ax2.tick_params(axis='y', colors='b')
			ax2.set_xlim([0.8,len(plotexpe)+0.2])

			if Exptag=='C':
				ax2.set_ylim([-0.01 ,0.04])
				if name=='id':
					ax2.set_ylim([4 ,14])
				elif name in ['jd','og','as','mw']:
					ax2.set_ylim([0 ,4])
				elif name in ['mw_prime','mh2','sb']:
					ax2.set_ylim([0 ,0.2])
			elif Exptag=='D':
				ax2.set_ylim([0.5 ,8])
			elif Exptag=='F':
				ax2.set_ylim([-0.1 ,5])

			#==== Core computation==============
			for j,tag in enumerate(plotexpe):
				print tag
				NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
				try:
					DatFile	 = Dataset(NCFile, mode='r')
					time     = DatFile.variables['time'][:]
					data,coords,datacoordtype=GetVar(DatFile,DataName,name,True)
					thick,coords,thickcoordtype=GetVar(DatFile,'H',name,False)
					base,coords,basecoordtype=GetVar(DatFile,'B',name,False)
					DatFile.close()
					if name=='db':
						RefFile						= Dataset(resdir+'mw/'+tag+'_mw.nc', mode='r')
						reftime						= RefFile.variables['time'][:]
						refdata,refcoords,refcoordtype	=	GetVar(RefFile,DataName,'mw',True)
						RefFile.close()
				except (RuntimeError,IOError):
					print ('File {} does not exist'.format(NCFile))
					if j==framelist:
						layout.delaxes(ax[k+1,0])
						layout.delaxes(ax[k+1,1])
					continue
				if Exptag=='C':
					time=24.+(time-max(time))/3600.
					Cycle=DiurnalForcing(time*3600.,AmpList[tag])
					Forcing=np.outer(MoulinInput,Cycle)
					if name=='db':
						reftime=24.+(reftime-max(reftime))/3600.
				else:
					time=1+(time-max(time))/(365.*24.*3600.)
					time=time[np.where(time>=0)][:]
					surf=thick+base
					Forcing=SeasonForcing(surf,time*3600.*24.*365.,AmpList[tag])
					if name=='db':
						reftime=1.+(reftime-max(reftime))/(365.*24.*3600.)
				if Exptag=='F':
					#===============computing glacier width============
					geomx=coords[0,:]
					benchsurf=100.0*(geomx+200.)**0.25+ 100.0*geomx/6.0e3-100.0*200.0**0.25
					surfend=100.0*(6.0e3+200.0)**0.25+100.0-100.0*200.0**0.25
					f10 = (surfend - 0.05*6.0e3)/6.0e3**2
					f0 = f10*geomx**2 + 300.*geomx/6.0e3
					h0 = -4.5*(geomx/6.0e3) + 5
					benchwidth=2.0*(((benchsurf-f0)/h0)/0.5e-6)**(1/3.)
				Band=np.NaN*np.ones((len(time),3,3))
				Mean=np.NaN*np.ones((len(time)))
				Force=np.NaN*np.ones((len(time),3))
				MeanForce=np.NaN*np.ones((len(time)))
				nodeselect=[]
				if name=='db':
					RefMean=np.NaN*np.ones((len(reftime)))
					for t,step in enumerate(time):
						Mean[t]=np.nanmean(data[t,:])
						if Exptag=='C':
							MeanForce[t]=np.nansum(Forcing[:,t])
						elif Exptag=='D':
							MeanForce[t]=np.nanmean(Forcing[:,t])*100.0e3*20.0e3
						elif Exptag=='F':
							MeanForce[t]=np.nanmean(Forcing[:,t]*benchwidth)*6.0e3
					for t,step in enumerate(reftime):
						RefMean[t]=np.nanmean(refdata[t,:])
				else:
					for b,bandname in enumerate(['Low','Mid','Top']):
						if Exptag=='F':
							nodeselect.append(np.where(np.logical_and(coords[0,:]>=BandProp[bandname][1],
																												coords[0,:]<=BandProp[bandname][1]+300))[0])
						else:
							nodeselect.append(np.where(np.logical_and(coords[0,:]>=BandProp[bandname][0],
																												coords[0,:]<=BandProp[bandname][0]+5000))[0])
					for t,step in enumerate(time):
						if Exptag=='C':
							MeanForce[t]=np.nansum(Forcing[:,t])
						elif Exptag=='D':
							MeanForce[t]=np.nanmean(Forcing[:,t])*100.0e3*20.0e3
						elif Exptag=='F':
							MeanForce[t]=np.nanmean(Forcing[:,t]*benchwidth)*6.0e3

						for b,bandname in enumerate(['Low','Mid','Top']):
							Band[t,b,0]=np.nanmin(data[t,nodeselect[b]])
							Band[t,b,1]=np.nanmean(data[t,nodeselect[b]])
							Band[t,b,2]=np.nanmax(data[t,nodeselect[b]])
							Mean[t]=np.nanmean(data[t,:])
							if Exptag=='C':
								Force[t,b]=np.nansum(Forcing[int(BandProp[bandname][0]*1.0e-3+1):,t])
							elif Exptag=='D':
								upstream=np.where(coords[0,:]>=BandProp[bandname][0])[0]
								Force[t,b]=np.nanmean(Forcing[upstream,t])*20.0e3*(100.0e3-BandProp[bandname][0])
							elif Exptag=='F':
								upstream=np.where(coords[0,:]>=BandProp[bandname][1])[0]
								Force[t,b]=np.nanmean(Forcing[upstream,t])*np.nanmean(benchwidth[upstream])*(6.0e3-BandProp[bandname][1])

				#=======plotting right column==============
				print('treating {} of {}'.format(tag,name))
				timeshift=time[np.where(np.nanmin(Mean)==Mean)[0][0]]-time[np.where(np.nanmax(MeanForce)==MeanForce)[0][0]]
				if Exptag=='C':
					if timeshift<-12:
						timeshift=timeshift+24
				elif Exptag in ['D','F']:
					timeshift=12.*timeshift
				if j==0 and k==0:
					col2a=ax[1,plotcol-1].scatter(int(tag[-1]), timeshift, marker='*',color='k',label='Effective pressure lag')
					col2b=ax2.scatter(int(tag[-1]),np.nanmax(Mean)-np.nanmin(Mean), marker='+',color='b',label='Effective pressure amplitude')
				else:
					ax[k+1,plotcol-1].scatter(int(tag[-1]), timeshift, marker='*',color='k')
					ax2.scatter(int(tag[-1]),np.nanmax(Mean)-np.nanmin(Mean), marker='+',color='b')

				#==============setting limits================
				if Exptag=='C':
					ax[k+1,plotcol-1].set_ylim([-12 ,12])
					ax[k+1,plotcol-1].set_xticks([1,2,3,4])
					ax[plotrow-1,plotcol-1].xaxis.set_ticklabels(['C1','C2','C3','C4'],fontsize=9)
				elif Exptag=='D':
					ax[k+1,plotcol-1].set_ylim([-1.8,1.44])
					ax[k+1,plotcol-1].set_xticks([1,2,3,4,5])
					ax[plotrow-1,plotcol-1].xaxis.set_ticklabels(['D1','D2','D3','D4','D5'],fontsize=9)
				elif Exptag=='F':
					ax[k+1,plotcol-1].set_ylim([-3 ,2])
					ax[k+1,plotcol-1].set_xticks([1,2,3,4,5])
					ax[plotrow-1,plotcol-1].xaxis.set_ticklabels(['F1','F2','F3','F4','F5'],fontsize=9)

				if j==framelist:
					#========patch on right to show what we have on left======
					if Exptag=='C':
						rect = [Rectangle((framelist+0.5,-12), 1.0, 24)]
					elif Exptag=='D':
						rect = [Rectangle((framelist+0.5,-1.8), 1.0, 3.24)]
					elif Exptag=='F':
						rect = [Rectangle((framelist+0.5,-3), 1.0, 5)]

					pc = PatchCollection(rect, facecolor=(0.8,0.8,0.8), alpha=1,
															 edgecolor=(0.8,0.8,0.8),zorder=0)
					ax[k+1,plotcol-1].add_collection(pc)

					#================setting up left plots====================
					for count in[0,1]:
						if Exptag=='C':
							ax[k+1,count].text(0.05,0.85,string.lowercase[(k+1)+count*(plotrow-1)],fontsize=9,va='center',ha='center',transform=ax[k+1,count].transAxes)
							ax[k+1,count].yaxis.set_major_locator(MaxNLocator(5,prune='both'))
						else:
							ax[k+1,count].text(0.05,0.85,string.lowercase[(k+2)+count*(plotrow-1)],fontsize=9,va='center',ha='center',transform=ax[k+1,count].transAxes)
							ax[k+1,count].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
						ax[k+1,1].set_title(namestyles[name][0],fontsize=9,va='top',ha='right',position=(0.97,0.8))

						ax[k+1,count].set_position([hleft+count*(plotwidth+hspace),vbot+(plotrow-2-k)*(plotheight+vspace),plotwidth,plotheight])
						ax[k+1,count].tick_params(axis='y', which='major', pad=2,top=False)

						ax[k,count].xaxis.set_ticklabels([])
						if Exptag=='C':
							ax[k+1,count].set_xlim([1 ,24])
							ax[k+1,count].set_xticks([3,6,9,12,15,18,21,24])
							ax[k+1,0].set_ylim([0.4,2.6])
							ax[k+1,1].set_ylim([0 ,120])
							if name=='id':
								ax[k+1,0].set_ylim([-18 ,12])
							elif name in ['jd','og','as','mw']:
								ax[k+1,0].set_ylim([-2 ,3.5])
						elif Exptag=='D':
							ax[k+1,count].set_xlim([0.25 ,1])
							ax[k+1,count].set_xticks([0.33,0.49,0.66,0.83,1.0])
							if k+1==plotrow-1 or (k+1==plotrow-2 and sub>3):
								ax[k+1,count].xaxis.set_ticklabels(['$4$','$6$','$8$','$10$','$12$'])
								ax[k+1,count].set_xlabel('Time (month)',labelpad=0)
							else:
								ax[k+1,count].xaxis.set_ticklabels([])
							if framelist==0:
								ax[k+1,0].set_ylim([-5,12])
								maxdischarge=400
							else:
								ax[k+1,0].set_ylim([-5,10])
								if framelist==1:
									maxdischarge=600
								elif framelist==2:
									maxdischarge=1000
								elif framelist==3:
									maxdischarge=1400
								elif framelist==4:
									maxdischarge=1800
							ax[k+1,1].set_ylim([0,maxdischarge])
							ax[0,0].set_ylim([0,maxdischarge])
							ax[0,1].set_ylim([0,maxdischarge])
						elif Exptag=='F':
							ax[k+1,count].set_xlim([0.25 ,1])
							ax[k+1,count].set_xticks([0.33,0.49,0.66,0.83,1.0])
							if k+1==plotrow-1:
								ax[k+1,count].xaxis.set_ticklabels(['$4$','$6$','$8$','$10$','$12$'])
								ax[k+1,count].set_xlabel('Time (month)',labelpad=0)
							else:
								ax[k+1,count].xaxis.set_ticklabels([])
							if name in ['db','jsb']:
								ax[k+1,0].set_ylim([0,2])
							else:
								ax[k+1,0].set_ylim([0,6])
							ax[k+1,1].set_ylim([0 ,2+2*framelist])
							ax[0,0].set_ylim([0 ,2+2*framelist])
							ax[0,1].set_ylim([0 ,2+2*framelist])

						if name=='bf':
							if Exptag=='C':
								ax[0,count].text(0.05,0.85,string.lowercase[0],fontsize=9,va='center',ha='center',transform=ax[0,count].transAxes)
							else:
								ax[0,count].text(0.05,0.85,string.lowercase[count],fontsize=9,va='center',ha='center',transform=ax[0,count].transAxes)
							ax[0,count].set_position([hleft+count*(plotwidth+hspace),vbot+(plotrow-1)*(plotheight+vspace),plotwidth,plotheight])
							ax[0,count].yaxis.set_major_locator(MaxNLocator(5,prune='both'))
							if Exptag=='C':
								ax[0,count].set_xlim([1 ,24])
								ax[0,count].set_xticks([3,6,9,12,15,18,21,24])
								ax[plotrow-1,count].xaxis.set_ticklabels(['$3$','$6$','$9$','$12$','$15$','$18$','$21$','$24$'])
								ax[plotrow-1,count].set_xlabel('Time (h)',labelpad=0)
							elif Exptag=='D':
								ax[0,count].set_xlim([0.25 ,1])
								ax[0,count].set_xticks([0.33,0.49,0.66,0.83,1.0])
							elif Exptag=='F':
								ax[0,count].set_xlim([0.25 ,1])
								ax[0,count].set_xticks([0.33,0.49,0.66,0.83,1.0])
							ax[plotrow-1,count].tick_params(axis='x', which='major', pad=2)
							if count==0:
								if Exptag=='C':
									topleft=ax[0,0].plot(time, MeanForce, linestyle='-',linewidth=1.5,color='k',label='Applied forcing')
								elif Exptag in['D','F']:
									topleft=ax[0,0].plot(time, MeanForce, linestyle='-',linewidth=1.5,color='k',label='Total recharge')
							else:
								if Exptag=='C':
									ax[0,1].plot(time, MeanForce, linestyle='-',linewidth=1.5,color='k')
								elif Exptag in ['D','F']:
									topright=[]
									for b,bandname in enumerate(['Low','Mid','Top']):
										topright.append(ax[0,1].plot(time, Force[:,b], linestyle='-',linewidth=1.5,color=BandProp[bandname][2],label='Recharge upstream of '+BandProp[bandname][4])[0])

					channelflux,sheetflux,lines=GetFluxRatio(name,tag,time)
					print('for {} in {} max channelflux is {}'.format(name,tag,np.nanmax(channelflux)))
					print('for {} in {} max sheetflux is {}'.format(name,tag,np.nanmax(sheetflux)))
					totflux=channelflux+sheetflux
					if name=='db':
						dbN=ax[k+1,0].plot(time, Mean[:], linestyle='-',linewidth=1.5,color='g',label='Mean effective pressure')[0]
						# IDS.append(ax[k+1,1].plot(time, np.nanmean(sheetflux,axis=1), linestyle='--',linewidth=1.5,color='g',label='Mean inefficient syst. flux')[0])
						# EDS.append(ax[k+1,1].plot(time, np.nanmean(channelflux,axis=1), linestyle=':',linewidth=1.5,color='g',label='Mean efficient syst. flux')[0])
						IDS.append(ax[k+1,1].plot(time, np.nanmean(sheetflux,axis=1), dashes=[6, 2],linewidth=1.5,color='g',label='Mean inefficient syst. flux')[0])
						EDS.append(ax[k+1,1].plot(time, np.nanmean(channelflux,axis=1), dashes=[2, 2],linewidth=1.5,color='g',label='Mean efficient syst. flux')[0])
					else:
						for b,bandname in enumerate(['Low','Mid','Top']):
							ax[k+1,0].plot(time, Band[:,b,1], linestyle='-',linewidth=1.5,color=BandProp[bandname][2])
							pressure.append(ax[k+1,0].fill_between(time,Band[:,b,0],Band[:,b,2],
																										 label='Effective pressure for '+BandProp[bandname][4],
																										 facecolor=BandProp[bandname][2],edgecolor=BandProp[bandname][2],alpha=0.3))
							if Exptag=='C':
								bandloc=np.where(lines==10000)[0]
								if k==0:
									# IDS=ax[1,1].plot(time, sheetflux[:,bandloc], linestyle='--',linewidth=1.5,color='b',label='Inefficient system discharge')[0]
									# EDS=ax[1,1].plot(time, channelflux[:,bandloc], linestyle=':',linewidth=1.5,color='b',label='Efficient system discharge')[0]
									IDS=ax[1,1].plot(time, sheetflux[:,bandloc], dashes=[6, 2],linewidth=1.5,color='b',label='Inefficient system discharge')[0]
									EDS=ax[1,1].plot(time, channelflux[:,bandloc], dashes=[2, 2],linewidth=1.5,color='b',label='Efficient system discharge')[0]
								else:
									# ax[k+1,1].plot(time, sheetflux[:,bandloc], linestyle='--',linewidth=1.5,color='b')
									# ax[k+1,1].plot(time, channelflux[:,bandloc], linestyle=':',linewidth=1.5,color='b')
									ax[k+1,1].plot(time, sheetflux[:,bandloc], dashes=[6, 2],linewidth=1.5,color='b')
									ax[k+1,1].plot(time, channelflux[:,bandloc], dashes=[2, 2],linewidth=1.5,color='b')
							elif Exptag in['D','F']:
								if Exptag=='D':
									bandloc=np.where(lines==BandProp[bandname][0])[0]
								else:
									bandloc=np.where(lines==BandProp[bandname][1])[0]
								if k==1:
									# IDS.append(ax[2,1].plot(time, sheetflux[:,bandloc], linestyle='--',linewidth=1.5,color=BandProp[bandname][2],label='Inef. syst. flux '+BandProp[bandname][4])[0])
									# EDS.append(ax[2,1].plot(time, channelflux[:,bandloc], linestyle=':',linewidth=1.5,color=BandProp[bandname][2],label='Eff. syst. flux '+BandProp[bandname][4])[0])
									IDS.append(ax[2,1].plot(time, sheetflux[:,bandloc], dashes=[6, 2],linewidth=1.5,color=BandProp[bandname][2],label='Inef. syst. flux '+BandProp[bandname][4])[0])
									EDS.append(ax[2,1].plot(time, channelflux[:,bandloc], dashes=[2, 2],linewidth=1.5,color=BandProp[bandname][2],label='Eff. syst. flux '+BandProp[bandname][4])[0])
								else:
									# ax[k+1,1].plot(time, sheetflux[:,bandloc], linestyle='--',linewidth=1.5,color=BandProp[bandname][2])
									# ax[k+1,1].plot(time, channelflux[:,bandloc], linestyle=':',linewidth=1.5,color=BandProp[bandname][2])
									ax[k+1,1].plot(time, sheetflux[:,bandloc],dashes=[6, 2],linewidth=1.5,color=BandProp[bandname][2])
									ax[k+1,1].plot(time, channelflux[:,bandloc],dashes=[2, 2],linewidth=1.5,color=BandProp[bandname][2])
					ax[k+1,0].plot(time, np.zeros(np.shape(time)), linestyle='--',linewidth=0.5,color='k')

	layout.delaxes(ax[0,-1])

	#================adding axis labels====================
	if Exptag=='F':
		spacer =0.3*hspace
	elif Exptag=='D':
		spacer =0.25*hspace
	elif Exptag=='C':
		spacer =0.25*hspace

	# adding a hline spacer on top
	hsep = ConnectionPatch(xyA=(0.5*hleft,1.-(vtop+0.5*vspace+plotheight)), xyB=(1-(0.5*hright),1.-(vtop+0.5*vspace+plotheight)),
												 coordsA="figure fraction", coordsB="figure fraction",
												 linewidth=0.5,color="black")
	plt.gca().add_artist(hsep)

	ax[-1,-1].text(spacer,1-vtop-0.5*plotheight,'Recharge (m$^3$s$^{-1}$)',fontsize=9,va='center',ha='center',rotation=90,transform=layout.transFigure)
	ax[-1,-1].text(spacer,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Effective pressure (MPa)',fontsize=9,va='center',ha='center',rotation=90,transform=layout.transFigure)
	ax[-1,-1].text(hleft+plotwidth+spacer,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Discharge (m$^3$s$^{-1}$)',fontsize=9,va='center',ha='center',rotation=90,transform=layout.transFigure)
	ax[-1,-1].text(1-0.2*hright,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Amplitude (MPa)',color='b',fontsize=9,va='center',ha='center',rotation=90,transform=layout.transFigure)
	if Exptag=='C':
		ax[-1,-1].text(hleft+2*plotwidth+hspace+spacer,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Time lag (h)',fontsize=9,va='center',ha='center',rotation=90,transform=layout.transFigure)
	elif Exptag in['D','F']:
		ax[-1,-1].text(hleft+2*plotwidth+hspace+spacer,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Time lag (month)',fontsize=9,va='center',ha='center',rotation=90,transform=layout.transFigure)

	print np.size(pressure)
	if Exptag=='C':
		L=plt.legend(handles=[pressure[0],pressure[1],pressure[2],topleft[0]],
								 bbox_to_anchor=(hleft,1-0.85*vtop,hspace+2.0*plotwidth,0.8*vtop),
								 bbox_transform=layout.transFigure,ncol=2,labelspacing=0.2,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,fontsize=9,
								 borderaxespad=0,handler_map={pressure[0]:FillandLineHandler('b'),pressure[1]:FillandLineHandler('c'),pressure[2]:FillandLineHandler('r')},mode='expand', frameon=False)
 		M=plt.legend(handles=[IDS,EDS,col2a,col2b],
								 bbox_to_anchor=(hleft+2.*plotwidth+0.85*hspace,1-0.9*vtop,0.5*plotwidth+hright+hspace,0.8*vtop),
								 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.2,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,fontsize=9,
								 borderaxespad=0,handler_map={pressure[2]:FillandLineHandler('r')},mode='expand', frameon=False)
	elif Exptag=='D':
		L=plt.legend(handles=[pressure[0],pressure[1],pressure[2],dbN,topleft[0],topright[0],topright[1],topright[2]],
								 bbox_to_anchor=(hleft,1-0.79*vtop,2*plotwidth+hspace,0.8*vtop),
								 bbox_transform=layout.transFigure,ncol=2,labelspacing=0.2,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,fontsize=9,
								 borderaxespad=0,handler_map={pressure[0]:FillandLineHandler('b'),pressure[1]:FillandLineHandler('c'),pressure[2]:FillandLineHandler('r')},mode='expand',frameon=False)
		M=plt.legend(handles=[EDS[0],EDS[1],EDS[2],IDS[0],IDS[1],IDS[2],col2a,col2b],
								 bbox_to_anchor=(hleft+2*plotwidth+0.9*hspace,1-0.79*vtop,0.5*plotwidth+hright+hspace,0.8*vtop),
								 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.2,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,fontsize=9,
								 borderaxespad=0,mode='expand',frameon=False)#
	elif Exptag=='F':
		L=plt.legend(handles=[pressure[0],pressure[1],pressure[2],dbN,topleft[0],topright[0],topright[1],topright[2]],
								 bbox_to_anchor=(hleft,1-0.79*vtop,2*plotwidth+hspace,0.8*vtop),
								 bbox_transform=layout.transFigure,ncol=2,labelspacing=0.2,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,fontsize=9,
								 borderaxespad=0,handler_map={pressure[0]:FillandLineHandler('b'),pressure[1]:FillandLineHandler('c'),pressure[2]:FillandLineHandler('r')},mode='expand',frameon=False)
		M=plt.legend(handles=[EDS[0],EDS[1],EDS[2],EDS[3],IDS[0],IDS[1],IDS[2],IDS[3],col2a,col2b],
								 bbox_to_anchor=(hleft+2*plotwidth+0.9*hspace,1-0.79*vtop,0.5*plotwidth+hright+hspace,0.8*vtop),
								 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.2,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,fontsize=9,
								 borderaxespad=0,mode='expand',frameon=False)#

	plt.gca().add_artist(L)
	plt.gca().add_artist(M)

	plt.show
	return layout

# }}}

# {{{ GetVar(DatFile,var,name,AllTime):

def GetVar(DatFile,var,name,AllTime):
	#getting data array
	if np.size(DatFile.variables[var].dimensions)==1:
		data = DatFile.variables[var][:]
	else:
		if AllTime:
			data = DatFile.variables[var][:,:]
			time = DatFile.variables['time'][:]
			step=(time[-1]-time[-2])/3600.
			if step>7:
				time=1+(time-max(time))/(365.*24.*3600.)
				data=data[np.where(time>=0)][:]
			else:
				time=24.+(time-max(time))/(3600.)
				data=data[np.where(time>-0.5)][:]
			if name=='id' and var=='N':
				data[:,-1]=data[:,-2]
		else:
			if DatFile.variables[var].dimensions[-1]=='time':
				data = DatFile.variables[var][:,-1]
			else:
				data = DatFile.variables[var][-1,:]
			if name=='id' and var=='N':
				data[-1]=data[-2]

	if var=='N':
		data=1.0e-6*data

	for key in DatFile.variables[var].dimensions:
		if key.startswith('index'):
			if key=='index_ch':
				datacoord='coords_ch'
			else:
				datacoord='coords'+str(key[-1])
			break
	if len(DatFile.dimensions['dim'])==1:
		if DatFile.variables[datacoord].dimensions[-1]=='dim':
			Xs		 = DatFile.variables[datacoord][:].T-94000
		else:
			Xs		 = DatFile.variables[datacoord][:]
		if np.max(Xs>=10*1.0e3):
			Ys		 = 10*1.0e3*np.ones(np.shape(Xs))
		else:
			Ys		 = np.zeros(np.shape(Xs))
	else:
		Xs		 = DatFile.variables[datacoord][0,:]
		Ys		 = DatFile.variables[datacoord][1,:]

	coords = np.vstack((Xs,Ys))
	return data, coords, datacoord

# }}}

# {{{ Griding(data,coords):

def Griding(data,coords):

	Xlims=[1.0e-3*np.min(coords[0,:]),min(100.,1.0e-3*np.max(coords[0,:]))]
	if np.all(coords[1,:]==10*1.0e3):
		Ylims=[0,20]
	else:
		Ylims=[1.0e-3*np.min(coords[1,:]),1.0e-3*np.max(coords[1,:])]

	grid_x, grid_y = np.mgrid[Xlims[0]:Xlims[1]:101j, Ylims[0]:Ylims[1]:101j]
	#taking the mean on all times before gridding

	# if expe in ['E','F']:
	# 	if len(np.shape(data))>1:
	# 		griddedsurf=griddata(1.0e-3*coords.T,np.nanmean(data,axis=0), (grid_x, grid_y), method='linear')
	# 	else:
	# 		griddedsurf=griddata(1.0e-3*coords.T,data, (grid_x, grid_y), method='linear')
	# else:
	if len(np.shape(data))>1:
		griddedsurf=griddata(1.0e-3*coords.T,np.nanmean(data,axis=0), (grid_x, grid_y), method='nearest')
	else:
		griddedsurf=griddata(1.0e-3*coords.T,data, (grid_x, grid_y), method='nearest')

	return griddedsurf,grid_x,grid_y

# }}}

# {{{ GetFluxRatio(name,tag,time):

def GetFluxRatio(name,tag,time):
	if tag.startswith('F') or tag.startswith('E'):
		Length=6.0e3
		Step=60.
		Lines=np.arange(0,Length+Step,Step)
		s2 = 100.0/Length
		s_xend = 100.0*(Length+200)**0.25 + s2*Length - 100.0*200**0.25
		surf = 100.0*(Lines+200)**0.25 + s2*Lines - 100.0*200**0.25
		f10 = (s_xend - 0.05*Length)/Length**2
		f0 = f10*Lines**2 + 0.05*Lines
		h0 = -4.5*(Lines/Length) + 5
		Width = 2.0*(((surf-f0)/h0)/0.5e-6)**(1/3.)
	else:
		Length=100.0e3
		Step=1000.
		Lines=np.arange(0,Length+Step,Step)
		Width=20.0e3*np.ones(np.shape(Lines))

	if len(time)>0:
		TimeVar=True
		SheetFlux=np.zeros((len(time),len(Lines)))
		ChannelFlux=np.zeros((len(time),len(Lines)))
	else:
		TimeVar=False
		SheetFlux=np.nan*np.ones(len(Lines))
		ChannelFlux=np.nan*np.ones(len(Lines))

	if name in['mh1','mh2','mw','mw_prime','og','og_prime','cdf','db','sb','jd','jsb','as','as_2','id','rh','bf']:
		NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			if name in['mh2','mw','mw_prime','og','og_prime','mh1']:
				IneffFlux,Ineffcoords,Innefcoordtype=GetVar(DatFile,'q',name,TimeVar)
				EffFlux,Effcoords,Effcoordtype=GetVar(DatFile,'Q',name,TimeVar)
				if name in ['mh2','mw','mw_prime'] and expe=='C' and sys.argv[1]=='TransientFlux':
					IneffFlux=np.insert(IneffFlux,0,IneffFlux[-1,:],axis=0)
					EffFlux=np.insert(EffFlux,0,EffFlux[-1,:],axis=0)
				if name=='mh1':
					Xs=DatFile.variables['coords3'][0,:]
				else:
					Xs=DatFile.variables['coords1'][0,:]
				connect=DatFile.variables['connect_ch'][:,:]
				intconnect=np.asarray(connect,dtype=int)
				DatFile.close()
			elif name in ['cdf','rh']:
				if name=='cdf':
					IneffFlux,Ineffcoords,Innefcoordtype=GetVar(DatFile,'q',name,TimeVar)
				EffFlux,Effcoords,Effcoordtype=GetVar(DatFile,'Q',name,TimeVar)
				DatFile.close()
			elif name in['jd','jsb']:
				IneffFlux,Ineffcoords,Innefcoordtype=GetVar(DatFile,'q',name,TimeVar)
				if expe=='C' and sys.argv[1]=='TransientFlux':
					IneffFlux=np.insert(IneffFlux,0,IneffFlux[-1,:],axis=0)
				DatFile.close()
			elif name in ['as','as_2']:
				IneffFlux,Ineffcoords,Innefcoordtype=GetVar(DatFile,'q',name,TimeVar)
				Channelisation,Effcoords,Effcoordtype=GetVar(DatFile,'Dc',name,TimeVar)
				EffFlux=IneffFlux*Channelisation
				IneffFlux=IneffFlux*(1-Channelisation)
				DatFile.close()
			elif name=='sb':
				IneffFlux,Ineffcoords,Ineffcoordtype=GetVar(DatFile,'q',name,TimeVar)
				Conductivity,Ineffcoords,Ineffcoordtype=GetVar(DatFile,'T',name,TimeVar)
				ischannel=np.ones(np.shape(Conductivity))
#				ischannel[np.where(Conductivity==np.nanmin(Conductivity))]=0
				ischannel[np.where(Conductivity<0.1)]=0
				EffFlux=IneffFlux*ischannel
				IneffFlux=IneffFlux*(1-ischannel)
				DatFile.close()
			elif name=='bf':
				IneffFlux,Ineffcoords,Ineffcoordtype=GetVar(DatFile,'q',name,TimeVar)
				EffFlux,Effcoords,Effcoordtype=GetVar(DatFile,'Q',name,TimeVar)
				Thickness,Effcoords,Thickcoordtype=GetVar(DatFile,'Ee',name,TimeVar)
				isactive=np.zeros(np.shape(EffFlux))
				isactive[np.where(EffFlux>0)]=1.
				DatFile.close()
			elif name=='db':
				IneffFlux,Ineffcoords,Inefcoordtype=GetVar(DatFile,'q',name,TimeVar)
				Channelisation,Effcoords,Effcoordtype=GetVar(DatFile,'R',name,TimeVar)
				EffFlux=IneffFlux*Channelisation
				IneffFlux=IneffFlux*(1-Channelisation)
				DatFile.close()
			elif name=='id':
				Channelisation,Ineffcoords,Inefcoordtype=GetVar(DatFile,'R',name,TimeVar)
				EffFlux,Effcoords,Effcoordtype=GetVar(DatFile,'Q',name,TimeVar)
				if expe=='C' and sys.argv[1]=='TransientFlux':
					Channelisation=np.insert(Channelisation,0,Channelisation[-1,:],axis=0)
					EffFlux=np.insert(EffFlux,0,EffFlux[-1,:],axis=0)
				IneffFlux=EffFlux*(1-Channelisation)
				EffFlux=EffFlux*Channelisation
				Xs=DatFile.variables['coords1'][0,:]
				connect=DatFile.variables['connect_ch'][:,:]
				intconnect=np.asarray(connect,dtype=int)
				DatFile.close()

		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))

		for i,loc in enumerate(Lines):
			#===mh2=mw=og===
			if name in['mh2','mw','mw_prime','og','og_prime','mh1']:
				sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																				 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				if name=='mh1':
					channelloc=np.where((Xs[intconnect[0,:]]-loc+1)*(Xs[intconnect[1,:]]-loc+1)<0)[0]
				else:
					channelloc=np.where((Xs[intconnect[0,:]]-loc)*(Xs[intconnect[1,:]]-loc)<0)[0]
				if loc in [min(Lines),max(Lines)] and len(channelloc)==0:
					channelloc=np.where((Xs[intconnect[0,:]]-loc)*(Xs[intconnect[1,:]]-loc)==0)[0]
				if TimeVar:
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*Width[i]
					if len(channelloc)>0:
						ChannelFlux[:,i]=np.nansum(abs(EffFlux[:,channelloc]),axis=1)
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*Width[i]
					if len(channelloc)>0:
						ChannelFlux[i]=np.nansum(abs(EffFlux[channelloc]))
			#===cdf===
			elif name=='cdf':
				sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																				 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				channelloc=np.where(np.logical_and(Effcoords[0,:]<(loc+0.5*Step),
																					 Effcoords[0,:]>=(loc-0.5*Step)))[0]
				if TimeVar:
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*Width[i]
					if len(channelloc)>0:
						ChannelFlux[:,i]=np.nanmean(abs(EffFlux[:,channelloc]),axis=1)*20.0
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*Width[i]
					if len(channelloc)>0:
						ChannelFlux[i]=np.nanmean(abs(EffFlux[channelloc]))*20.0
			#===rh===
			elif name=='rh':
				channelloc=np.where(np.logical_and(Effcoords[0,:]<(loc+0.5*Step),
																					 Effcoords[0,:]>=(loc-0.5*Step)))[0]
				if len(channelloc)>0:
					if TimeVar:
						ChannelFlux[:,i]=np.nanmean(abs(EffFlux[:,channelloc]),axis=1)
					else:
						ChannelFlux[i]=np.nanmean(abs(EffFlux[channelloc]))
			#===id===
			elif name=='id':
				channelloc=np.where((Xs[intconnect[0,:]]-loc)*(Xs[intconnect[1,:]]-loc)<0)[0]
				if TimeVar:
					if len(channelloc)>0:
						ChannelFlux[:,i]=np.nanmean(abs(EffFlux[:,channelloc]),axis=1).T
						SheetFlux[:,i]=np.nanmean(abs(IneffFlux[:,channelloc]),axis=1).T
				else:
					if len(channelloc)>0:
						ChannelFlux[i]=np.sum(abs(EffFlux[channelloc])).T
						SheetFlux[i]=np.sum(abs(IneffFlux[channelloc])).T
			#===as==
			elif name in ['as','as_2']:
				sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																				 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				if TimeVar:
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*Width[i]
					ChannelFlux[:,i]=np.nanmean(EffFlux[:,sheetloc],axis=1)*Width[i]
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*Width[i]
					ChannelFlux[i]=np.nanmean(EffFlux[sheetloc])*Width[i]
			#===jd==
			elif name=='jd':
				sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																				 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				if TimeVar:
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*Width[i]
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*Width[i]
			#===db===
			elif name=='db':
				sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																				 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				if TimeVar:
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*np.nanmean(Width)
					ChannelFlux[:,i]=np.nanmean(EffFlux[:,sheetloc],axis=1)*np.nanmean(Width)
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*np.nanmean(Width)
					ChannelFlux[i]=np.nanmean(EffFlux[sheetloc])*np.nanmean(Width)
			#===jsb=
			elif name=='jsb':
				sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																				 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				if TimeVar:
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*np.nanmean(Width)
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*np.nanmean(Width)
			#===sb===
			elif name in ['sb','sb_old']:
				#for sheet flux
				sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																				 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				if TimeVar:
					channelratio=np.nansum(ischannel[:,sheetloc],axis=1)/len(sheetloc)
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*Width[i]*(1-channelratio[:])
					ChannelFlux[:,i]=np.nanmean(EffFlux[:,sheetloc],axis=1)*Width[i]*channelratio[:]
				else:
					channelratio=np.nansum(ischannel[sheetloc])/len(sheetloc)
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*Width[i]*(1-channelratio)
					ChannelFlux[i]=np.nanmean(EffFlux[sheetloc])*Width[i]*channelratio
			#===bf===
			elif name=='bf':
				#for sheet flux
				sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																				 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				channelloc=np.where(np.logical_and(Effcoords[0,:]<(loc+0.5*Step),
																					 Effcoords[0,:]>=(loc-0.5*Step)))[0]
				if TimeVar:
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*Width[i]
					ChannelFlux[:,i]=Width[i]*np.nanmean(abs(EffFlux[:,channelloc]),axis=1)

					ChannelFlux[np.where(np.isnan(ChannelFlux[:,i])),i]=0.0
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*Width[i]
					ChannelFlux[i]=Width[i]*np.nanmean(abs(EffFlux[channelloc]))
					if np.isnan(ChannelFlux[i]):
						ChannelFlux[i]=0.0

	return ChannelFlux,SheetFlux,Lines

# }}}

# {{{ ExpeLayout()

def ExpeLayout():

	fig,ax=plt.subplots(len(plotexpe),1,sharex='col',figsize=(6,12))
	hleft=0.1
	width=0.8
	vside=0.08
	vspace=0.0
	height=0.9/len(plotexpe)

	for i,expe in enumerate(plotexpe):
		ax[i].set_position([hleft,1-(vside+(i+1)*height),width,height])
		ax[i].set_ylabel(expe, fontsize=15)
	return fig,ax

# }}}

# {{{ DiurnalForcing(time,RelAmp=None):

def DiurnalForcing(time,RelAmp=None):
	'''Runoff=SeasonForcing(surf,time,Tparam)
	Definition of the diurnal input function for experiment C of the SHMIP exercise
	This is only the Moulin input part, distributed source of run A1 has to be added
	Input:
	time : are the timestamps of the requested timesteps (one vector)
	RelAmp: is the relative amplitude of the signal as given in SHMIP
	Output:
	OUtput is the transformed input which has to be used for every moulin.
	'''

	day = 24.*60.*60. #From second to day
	MI = 0.9 #initial noulin input
	if RelAmp is None:
		ra = 0.
	else:
		ra = RelAmp

		Input=MI*(1.-ra*np.sin(2.*np.pi*time/day))
		Input[np.where(Input<0)]=0.0

	return Input

# }}}

# {{{ SeasonForcing(surf,time,Tparam=None):

def SeasonForcing(surf,time,Tparam=None):
	'''Runoff=SeasonForcing(surf,time,Tparam)
  '''
	day = 24.*60.*60. #From second to day
	year = 365 *day   # From second to year
	DDF = 0.01/day ## degree day factor m/K/s
	lr = -0.0075 #lapse rate K/m
	BackInput = 7.93e-11 #value of the background input

	if Tparam is None:
		Dt = 0.
	else:
		Dt = Tparam

	Temp= -16.*np.cos(2.*np.pi*time/year)-5 +Dt

	TotInput=BackInput*np.ones((np.size(surf),np.size(Temp)))
	for t,value in enumerate(Temp):
		PosVal=np.where((surf*lr+value)>0)[0]
		TotInput[PosVal,t]+=(surf[PosVal]*lr+value)*DDF

	return TotInput

# }}}

# {{{ColorMapping():

def ColorMapping():
	cdict1 = {'alpha': [(0.0, 1.0, 1.0), (0.05, 1.0, 1.0), (0.1, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
					 'blue': [(0.0, 0.96, 0.96), (0.05, 1.0, 1.0), (0.1, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
					 'green': [(0.0, 0.74, 0.74), (0.05, 1.0, 1.0), (0.1, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
					 'red': [(0.0, 0.16, 0.16), (0.05, 1.0, 1.0), (0.1, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)]}
	blue_red1 = colors.LinearSegmentedColormap('BlueRed1', cdict1)
	return blue_red1

# }}}

# {{{ def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
	'''
	Function to offset the "center" of a colormap. Useful for
	data with a negative min and positive max and you want the
	middle of the colormap's dynamic range to be at zero

	Input
	-----
	cmap : The matplotlib colormap to be altered
	start : Offset from lowest point in the colormap's range.
	Defaults to 0.0 (no lower ofset). Should be between
	0.0 and `midpoint`.
	midpoint : The new center of the colormap. Defaults to
	0.5 (no shift). Should be between 0.0 and 1.0. In
	general, this should be  1 - vmax/(vmax + abs(vmin))
	For example if your data range from -15.0 to +5.0 and
	you want the center of the colormap at 0.0, `midpoint`
	should be set to  1 - 5/(5 + 15)) or 0.75
	stop : Offset from highets point in the colormap's range.
	Defaults to 1.0 (no upper ofset). Should be between
	`midpoint` and 1.0.
	'''
	cdict = {'red': [],
					 'green': [],
					 'blue': [],
					 'alpha': []}
	print midpoint
	if midpoint<=0:
		# regular index to compute the colors
		reg_index = np.linspace(0.5, stop, 257)
		# shifted index to match the data
		shift_index = np.linspace(0.0, 1.0, 257, endpoint=True)
	else:
		# regular index to compute the colors
		reg_index = np.linspace(start, stop, 257)
		# shifted index to match the data
		shift_index = np.hstack([np.linspace(0.0, midpoint, 128, endpoint=False),
														 np.linspace(midpoint, 1.0, 129, endpoint=True)])

	for ri, si in zip(reg_index, shift_index):
		r, g, b, a = cmap(ri)

		cdict['red'].append((si, r, r))
		cdict['green'].append((si, g, g))
		cdict['blue'].append((si, b, b))
		cdict['alpha'].append((si, a, a))

	newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
	plt.register_cmap(cmap=newcmap)

	return newcmap

# }}}

# {{{def CenteredSymLogNorm(cmap,linthresh=0.2*np.diff(lims), linscale=2,vmin=lims[0], vmax=lims[1]):

def CenteredSymLogNorm(cmap, TreshRatio=0.1, linlenght=1 , minval=-1., maxval=1.):
	'''
	Function to define a zero centered colorbar (for diverging colorbars)
	with a SymLog normalisation
	linscalelenght is the half lenght of the linear portion (in log decades)
	'''
	extent=maxval-minval
	linearlenght=TreshRatio*extent
	#		first define the norm
	newnorm=mpl.colors.SymLogNorm(linthresh=linearlenght, linscale=linlenght, vmin=minval, vmax=maxval)
	# and compute midpoint position:

	if maxval>linearlenght:
		posextent=linlenght+np.log10(maxval-linearlenght)
		colsections=2
	else:
		posextent=maxval/linearlenght
		colsections=1
	if minval<-linearlenght:
		negextent=linlenght+np.log10(-minval-linearlenght)
		colsections+=2
	else:
		negextent=-minval/linearlenght
		colsections+=1
	midpoint=negextent/(negextent+posextent)


	cdict = {'red': [],
					 'green': [],
					 'blue': [],
					 'alpha': []}

	# regular index to compute the colors
	colorvalues=int(257/colsections)*colsections

	if minval<-linearlenght:
		neglinpoint=(negextent-linlenght)/(negextent+posextent)
		reg_index= np.hstack([0.0,np.geomspace(1.0e-3, neglinpoint, (colorvalues/colsections)-1, endpoint=False),
															np.linspace(neglinpoint,0.5, colorvalues/colsections, endpoint=False)])
		shift_index = np.hstack([0.0,np.geomspace(1.0e-3, neglinpoint, (colorvalues/colsections)-1, endpoint=False),
																 np.linspace(neglinpoint,midpoint, colorvalues/colsections, endpoint=False)])
	elif minval<0:
		reg_index = np.hstack([np.linspace(0.0, 0.5, colorvalues/colsections, endpoint=False)])
		shift_index = np.hstack([np.linspace(0.0, midpoint, colorvalues/colsections, endpoint=False)])
	else:
		reg_index=[]
		shift_index=[]

	if maxval>linearlenght:
		if minval>0:
			poslinpoint=(negextent+linlenght)/(negextent+posextent)
			reg_index = np.hstack([reg_index,
														 np.linspace(0.5, poslinpoint, colorvalues/colsections, endpoint=False),
														 np.geomspace(poslinpoint, 1.0, colorvalues/colsections, endpoint=True)])
			shift_index = np.hstack([shift_index,
															 np.linspace(midpoint, poslinpoint, colorvalues/colsections, endpoint=False),
															 np.geomspace(poslinpoint, 1.0, colorvalues/colsections, endpoint=True)])
		else:
			poslinpoint=(negextent+linlenght)/(negextent+posextent)
			reg_index = np.hstack([reg_index,
														 np.linspace(0.5, poslinpoint, colorvalues/colsections, endpoint=False),
														 np.geomspace(poslinpoint, 1.0, colorvalues/colsections, endpoint=True)])
			shift_index = np.hstack([shift_index,
															 np.linspace(midpoint, poslinpoint, colorvalues/colsections, endpoint=False),
															 np.geomspace(poslinpoint, 1.0, colorvalues/colsections, endpoint=True)])
	else:
		reg_index = np.hstack([reg_index_index,
													 np.linspace(0.5, 1.0, colorvalues/colsections, endpoint=True)])
		shift_index = np.hstack([shift_index,
														 np.linspace(midpoint, 1.0, colorvalues/colsections, endpoint=True)])

	print minval
	print maxval
	print

	for ri, si in zip(reg_index, shift_index):
		r, g, b, a = cmap(ri)

		cdict['red'].append((si, r, r))
		cdict['green'].append((si, g, g))
		cdict['blue'].append((si, b, b))
		cdict['alpha'].append((si, a, a))

	newmap = mpl.colors.LinearSegmentedColormap('LogNormCentered', cdict)
	return newnorm,newmap

# }}}

class FillandLineHandler(object):
	def __init__(self, color):
		self.color=color
	def legend_artist(self, legend, orig_handle, fontsize, handlebox):
		x0, y0 = handlebox.xdescent, handlebox.ydescent
		width, height = handlebox.width, handlebox.height
		patch = plt.Rectangle([x0, y0], width, height, facecolor=self.color,alpha=0.3,
													edgecolor=self.color, transform=handlebox.get_transform())
		patch2 = plt.Line2D((x0,x0+width), (y0+0.5*height,y0+0.5*height),
												color=self.color, transform=handlebox.get_transform())
		handlebox.add_artist(patch)
		handlebox.add_artist(patch2)
		return patch

class MidpointNormalize(colors.Normalize):
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y))

class OOMFormatter(ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % self.format

class FormatScalarFormatter(ticker.ScalarFormatter):
	def __init__(self, fformat="%1.1f", offset=True, mathText=True):
		self.fformat = fformat
		ticker.ScalarFormatter.__init__(self,useOffset=offset,
																		useMathText=mathText)
	def _set_format(self, vmin, vmax):
		self.format = self.fformat
		if self._useMathText:
			self.format = '$%s$' % self.format

if __name__=="__main__":
	PlotList={'SecDiff':SectionsDiff,
						'BandEvol':BandEvol,
						'SteadyState':Steady,
						'Transient':Transient,
						'TransientFlux':TransientFlux,
						'DiurnalComp':DiurnalComp}


	TitleList={'N':['effective pressure',u'MPa']}
	if len(sys.argv)==1: # print help if no filename plot name is given:
		print __doc__
		sys.exit(0)

	resdir=os.getcwd()+'/'
	figdir=resdir+'/figures'

	DataName='N'
	expe=raw_input('Give the tag of the experiment you want to plot: ')

#	names=['bf']
	names =['db','id','rh','cdf','jd','jsb','as','sb','bf','mh1','mh2','og','og_prime','mw','mw_prime']#,,'as_2'

	namestyles = {'db':['\\boldmath$db$',(0.,0.,0.,1.)],
								'rh':['\\boldmath$rh$',(0.7,0.7,0.,1.)],
								'cdf':['\\boldmath$cdf$',(0.3,0.3,0.3,1.)],
								'id':['\\boldmath$id$',(0.6,0.6,0.6,1.)],
								'jd':['\\boldmath$jd$',(0.,0.,0.5,1.)],
								'jsb':['\\boldmath$jsb$',(0.,0.,0.7,1.)],
								'as':['\\boldmath$as$',(0.,0.,1.,1.)],
								'as_2':['\\boldmath$as\'$',(0.,0.,1.,1.)],
								'sb':['\\boldmath$sb$',(0.7,0.7,0.,1.)],
								'bf':['\\boldmath$bf$',(0.5,0.5,0.,1.)],
								'mh1':['\\boldmath$mh_1$',(0.,1.,1.,1.)],
								'mh2':['\\boldmath$mh_2$',(0.,0.5,0.5,1.)],
								'og':['\\boldmath$og$',(0.,1.,1.,0.3)],
								'og_prime':['\\boldmath$og\'$',(0.,0.,1.,1.)],
								'mw':['\\boldmath$mw$',(0.,.7,0.7,0.6)],
								'mw_prime':['\\boldmath$mw\'$',(0.,.7,0.7,0.6)],
								'mw_A':['\\boldmath$mwA$',(0.,.7,0.7,0.6)],
								'mw_A_ev':['\\boldmath$mwA\'$',(0.,.7,0.7,0.6)],
								'sb_old':['\\boldmath$sbO$',(0.,.7,0.7,0.6)]}


	if expe in ['A']:
		plotexpe=[expe+'1',expe+'2',expe+'3',expe+'4',expe+'5',expe+'6']
		subs=[1,2,3,4,5,6]
	elif expe in ['D','E','F']:
		plotexpe=[expe+'1',expe+'2',expe+'3',expe+'4',expe+'5']
		subs=[1,2,3,4,5]
	elif expe in ['C']:
		plotexpe=[expe+'1',expe+'2',expe+'3',expe+'4']
		subs=[1,2,3,4]
	elif expe in['B']:
		plotexpe=[expe+'1',expe+'4']
		subs=[1,2,3,4,5]
	os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
	try:
		#This is for suplementaries (need to put sub in select frame above)
		if sys.argv[1]=='DiurnalComp':
			subs=[0]
			plotexpe=plotexpe[:-1]
		for subsuite in subs:
	 		fig=PlotList[sys.argv[1]](subsuite)
			if fig!=0:
				try:
					fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+str(subsuite)+'.png',format='png')
					fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+str(subsuite)+'.pdf',format='pdf')
				except IOError:
					print('No Figure directory existing where needed, creating it here "{}"'.format(os.path.abspath(figdir)))
					os.makedirs(figdir+'/all')
					fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+str(subsuite)+'.png',format='png')
					fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+str(subsuite)+'.pdf',format='pdf')
	except KeyError:

		print 'Bad keyword for the figure production, run without keyword to get the help and check the spelling.'

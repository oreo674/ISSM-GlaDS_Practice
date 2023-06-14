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
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d.art3d import PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from netCDF4 import Dataset

from matplotlib import rc
sizeoffont=18
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'size':sizeoffont})
rc('legend',scatterpoints=1)

# {{{ MultiSec()

def MultiSec():
	layout,az=SliceLayout(13)
	midpoint=10.
	k=-1
	
	for i,name in enumerate(names):
		if len(glob(resdir+name+'/A*_'+name+'.nc'))>0:
			k+=1
			verts=[]
			limitedverts=[]
			for j,tag in enumerate(plotexpe):
				NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
				try:
					DatFile	 = Dataset(NCFile, mode='r')
					data,coords=GetVar(DatFile,DataName,name,False)
					DatFile.close()
				except (RuntimeError,IOError):
					print ('File {} does not exist'.format(NCFile))
					continue 
			
				griddeddata,gridpoint_x,gridpoint_y=Griding(data,coords)
				gridx=gridpoint_x[:,0]
				meandata=np.nanmean(griddeddata,axis=1)
				verts.append(list(zip(gridx,j*np.ones(np.shape(meandata)),meandata)))
				limitedverts=verts

				x=np.tile(gridpoint_x[:,0],(2,1))
				y=np.tile(j*np.ones(np.shape(gridpoint_x[:,0])),(2,1))
				z=np.vstack((meandata,np.zeros(np.shape(meandata))))

				channelflux,sheetflux,lines=GetFluxRatio(name,tag,[])
				if all(np.isnan(channelflux)):
					ratio=np.zeros(np.shape(sheetflux))
				elif all(np.isnan(sheetflux)):
					ratio=100.*np.ones(np.shape(channelflux))
				else:
					ratio=100.*channelflux/(sheetflux+channelflux)
				print ('ratio min is {} and max is {} for {} {}'.format(np.nanmin(ratio),np.nanmax(ratio),name,tag))
				
				color_dimension = np.tile(ratio,(2,1)) # change to desired fourth dimension
				norm=colors.Normalize(0,100)
				cmap= ColorMapping()#mpl.cm.bwr
				m = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
				m.set_array([])
				fcolors = m.to_rgba(color_dimension)

				az[k].plot_surface(x,y,z,rstride=1,cstride=1,facecolors=fcolors,vmin=0.0,vmax=1.0,shade=False)

			if name=='sb':
				az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','w','k','w','k'], linewidth=1, linestyles='-'))
			elif name=='bf':
				az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','w','k','w'], linewidth=1, linestyles='-'))
			elif name=='cdf':
				az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','k','w','k'], linewidth=1, linestyles='-'))
			elif name=='as':
				az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','k','k','w','k'], linewidth=1, linestyles='-'))
			elif name in ['mw']:
				az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','w','k','w','k'], linewidth=1, linestyles=['-','-',':','-',':','-']))
			else:
				az[k].add_collection3d(Line3DCollection(verts, colors='k', linewidth=1, linestyles='-'))

			if name in ['sb','bf','cdf','as']:
				az[k].set_title(namestyles[name][0],rotation=10,
												va='bottom',ha='left',
												position=(0.16,0.7),color='r')
				az[k].text(90,5,13,string.lowercase[k],va='center',ha='center')
			else:
				az[k].set_title(namestyles[name][0],rotation=10,
												va='bottom',ha='left',
												position=(0.16,0.7))
				az[k].text(90,5,13,string.lowercase[k],va='center',ha='center')

			az[k].invert_yaxis()
			az[k].view_init(elev=25., azim=245)
			
	barax = layout.add_axes([0.62,0.1,0.35,0.05])
	bartitle=barax.set_title('Percentage of Flux \n in the efficient system')
	cbar = mpl.colorbar.ColorbarBase(barax, norm=norm,cmap=cmap, orientation='horizontal')
	cbar.solids.set_edgecolor("face")
	plt.show
	return layout

# }}}
# {{{ MwerSec()

def MwerSec():
	rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'size':15})
	layout = plt.figure(figsize=(6,5))
	az=layout.add_subplot(1, 1, 1, projection='3d')
	layout.subplots_adjust(left=0., right=0.8, top=1.03, bottom=0.03, wspace=0.0, hspace=0.0)
	az.patch.set_visible(False) 
	az.set_xlabel('$x$ (km)')
	az.tick_params(axis='x', which='major')#, pad=0)
	az.set_xlim3d(0, 100)
	az.set_xticks([10,50,90])
	az.set_xticklabels(['$10$','$50$','$90$'],ha='center',va='bottom',rotation=10)

	az.set_ylim3d(0,len(plotexpe))
	az.set_ylim3d(-0.1,5.1)
	az.set_yticks([0,1,2,3,4,5])
	az.set_yticklabels(['$1$','$2$','$3$','$4$','$5$','$6$'],ha='right',va='center')#
	az.set_ylabel('A')#,labelpad=-10)
	az.tick_params(axis='y', which='major')# pad=-5)

	az.zaxis.set_rotate_label(False)
	az.set_zlabel('N (MPa)',rotation=90)#,labelpad=-5
	az.set_zlim3d(0, 14)
	az.set_zticks([0,2,4,6,8,10,12,14])
	az.set_zticklabels(['$0$','$2$','$4$','$6$','$8$','$10$','$12$','$14$'],ha='right',va='center')
	az.tick_params(axis='z', which='major')#pad=-3
	
	az.set_frame_on(False)

	verts=[]
	for j,tag in enumerate(plotexpe):
		NCFile=resdir+'mw/'+tag+'_mw.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			data,coords=GetVar(DatFile,DataName,'mw',False)
			DatFile.close()
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))
			continue 
			
		griddeddata,gridpoint_x,gridpoint_y=Griding(data,coords)
		gridx=gridpoint_x[:,0]
		meandata=np.nanmean(griddeddata,axis=1)
		verts.append(list(zip(gridx,j*np.ones(np.shape(meandata)),meandata)))
		limitedverts=verts

		x=np.tile(gridpoint_x[:,0],(2,1))
		y=np.tile(j*np.ones(np.shape(gridpoint_x[:,0])),(2,1))
		z=np.vstack((meandata,np.zeros(np.shape(meandata))))
		print('shape of verts {}'.format(np.shape(verts)))
		channelflux,sheetflux,lines=GetFluxRatio('mw',tag,[])
		ratio=100.*channelflux/(sheetflux+channelflux)
				
		color_dimension = np.tile(ratio,(2,1)) # change to desired fourth dimension
		norm=colors.Normalize(0,100)
		cmap= ColorMapping()#mpl.cm.bwr
		m = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
		m.set_array([])
		fcolors = m.to_rgba(color_dimension)

		az.plot_surface(x,y,z,rstride=1,cstride=1,facecolors=fcolors,vmin=0.0,vmax=1.0,shade=False)
	
		az.add_collection3d(Line3DCollection(verts, colors=['k','k','w','k','w','k'], linewidth=2, linestyles=['-','-',':','-',':','-']))
		az.set_title(namestyles['mw'][0],rotation=10,
									va='bottom',ha='left',
									position=(0.16,0.7))

	az.invert_yaxis()
	az.view_init(elev=25., azim=245)
			
	barax = layout.add_axes([0.76,0.05,0.1,0.8])
	bartitle=barax.set_title('Percentage of Flux \n in the efficient system', y=1.02)
	cbar = mpl.colorbar.ColorbarBase(barax, norm=norm,cmap=cmap, orientation='vertical')
	cbar.solids.set_edgecolor("face")
	plt.show
	return layout

# }}}
# {{{ MwerMeanSec()

def MwerMeanSec():
	width=17.8/(5*2.54)
	layout = plt.figure(figsize=(width,width))
	az=layout.add_subplot(1, 1, 1, projection='3d')
	layout.subplots_adjust(left=0., right=1.1, top=1.14, bottom=0.03, wspace=0.0, hspace=0.0)
	az.patch.set_visible(False) 
	az.tick_params(axis='x', which='major', pad=-3)
	az.set_xlim3d(0, 100)
	az.set_xticks([10,50,90])
	az.set_xticklabels(['$10$','$50$','$90$'],ha='center',va='center',rotation=10)

	az.set_ylim3d(0,len(plotexpe))
	az.set_ylim3d(-0.1,5.1)
	az.set_yticks([0,1,2,3,4,5])
	az.set_yticklabels(['','$2$','','$4$','','$6$'],ha='right',va='center')#
	az.set_ylabel('A',labelpad=-10)
	az.tick_params(axis='y', which='major', pad=-5)

	az.zaxis.set_rotate_label(False)
	az.set_zlim3d(0, 14)
	az.set_zticks([0,2,4,6,8,10,12,14])
	az.set_zticklabels(['$0$','','$4$','','$8$','','$12$',''],ha='right',va='center')
	az.tick_params(axis='z', which='major',pad=-5)
	
	az.set_frame_on(False)

	verts=[]
	for j,tag in enumerate(plotexpe):
		NCFile=resdir+'mw/'+tag+'_mw.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			data,coords=GetVar(DatFile,DataName,'mw',False)
			DatFile.close()
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))
			continue 
			
		meandata=np.nanmean(data)
		verts.append(list(zip([0,100],[j,j],[meandata,meandata])))
		limitedverts=verts

		x=np.tile([0,100],(2,1))
		y=np.tile([j,j],(2,1))
		z=np.vstack(([meandata,meandata],[0,0]))
		channelflux,sheetflux,lines=GetFluxRatio('mw',tag,[])
		ratio=np.nanmean(100.*channelflux/(sheetflux+channelflux))
				
		color_dimension = np.tile([ratio,ratio],(2,1)) # change to desired fourth dimension
		norm=colors.Normalize(0,100)
		cmap= ColorMapping()#mpl.cm.bwr
		m = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
		m.set_array([])
		fcolors = m.to_rgba(color_dimension)

		az.plot_surface(x,y,z,rstride=1,cstride=1,facecolors=fcolors,vmin=0.0,vmax=1.0,shade=False)
	
		az.add_collection3d(Line3DCollection(verts, colors=['k','k','w','k','w','k'], linewidth=1, linestyles=['-','-',':','-',':','-']))
		az.set_title(namestyles['mw'][0],rotation=10,
									va='bottom',ha='left',
									position=(0.16,0.7))

	az.invert_yaxis()
	az.view_init(elev=25., azim=245)

	plt.show
	return layout

# }}}
# {{{ BandEvol()

def BandEvol():
	width=18.
	height=8.
	plotrow=8
	plotcol=3
	layout,ax=plt.subplots(plotrow,plotcol,figsize=(width,height))
	Exptag='D'
	baselist=[3]#[1,3]
	hleft=1.1/width
	hright=6./width
	hspace=1./width#0.15
	vside=0.5/height
	vbot=vside
	vtop=0.1*vside
	vspace=0.05/height
	plotwidth=(1-hleft-hright-(2*hspace))/(plotcol-0.5)
	plotheight=(1.0-(1.1*vside)-(vspace*(plotrow-1)))/plotrow

	k=-1

	AmpList={'D1':-4.0,'D2':-2.0,'D3':0.0,'D4':2.0,'D5':4.0,}
	BandProp={'Low':[10000.,600.,'b','$Low$','lower band'],
						'Mid':[50000.,3000.,'c','$Mid$','middle band'],
						'Top':[90000.,5100.,'r','$Top$','highest band']}

	for i,name in enumerate(names):
		if len(glob(resdir+name+'/'+Exptag+'*_'+name+'.nc'))>0:
			if name=='mh1':
				framelist=[1]
			else:
				framelist=baselist
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
			ax[plotrow-1,plotcol-1].tick_params(axis='x', which='major', pad=5)
			ax[plotrow-1,plotcol-1].set_xlabel('sub-experiment',labelpad=0)

			#=====Axes for right column right axis===========
			ax2=ax[k+1,plotcol-1].twinx()
			ax2.set_position(ax[k+1,plotcol-1].get_position())
			ax2.yaxis.set_major_locator(MaxNLocator(5,prune='both'))
			ax2.spines['right'].set_color('b')
			ax2.tick_params(axis='y', colors='b')
			ax2.set_xlim([0.8,len(plotexpe)+0.2])

			ax2.set_ylim([0.5 ,8])

			#==== Core computation==============
			for j,tag in enumerate(plotexpe):
				NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
				try:
					DatFile	 = Dataset(NCFile, mode='r')
					time     = DatFile.variables['time'][:]
					data,coords=GetVar(DatFile,DataName,name,True)
					thick,coords=GetVar(DatFile,'H',name,False)
					base,coords=GetVar(DatFile,'B',name,False)
					DatFile.close()
				except (RuntimeError,IOError):
					print ('File {} does not exist'.format(NCFile))
					if j in framelist:
						layout.delaxes(ax[k+1,0])
						layout.delaxes(ax[k+1,1])
					continue 
				time=1+(time-max(time))/(365.*24.*3600.)
				time=time[np.where(time>=0)][:]
				surf=thick+base
				Forcing=SeasonForcing(surf,time*3600.*24.*365.,AmpList[tag])

				Band=np.NaN*np.ones((len(time),3,3))
				Mean=np.NaN*np.ones((len(time)))
				Force=np.NaN*np.ones((len(time),3))
				MeanForce=np.NaN*np.ones((len(time)))
				Meanwidth=np.NaN*np.ones((3))
				nodeselect=[]
				if name=='db':
					for t,step in enumerate(time):
						Mean[t]=np.nanmean(data[t,:])
						MeanForce[t]=np.nanmean(Forcing[:,t])*100.0e3*20.0e3
				else:
					for b,bandname in enumerate(['Low','Mid','Top']):
						nodeselect.append(np.where(np.logical_and(coords[0,:]>=BandProp[bandname][0],
																											coords[0,:]<=BandProp[bandname][0]+5000))[0])
					for t,step in enumerate(time):
						MeanForce[t]=np.nanmean(Forcing[:,t])*100.0e3*20.0e3
							
						for b,bandname in enumerate(['Low','Mid','Top']):
							Band[t,b,0]=np.nanmin(data[t,nodeselect[b]])
							Band[t,b,1]=np.nanmean(data[t,nodeselect[b]])
							Band[t,b,2]=np.nanmax(data[t,nodeselect[b]])
							Mean[t]=np.nanmean(data[t,:])
							upstream=np.where(coords[0,:]>=BandProp[bandname][0])[0]
							Force[t,b]=np.nanmean(Forcing[upstream,t])*20.0e3*(100.0e3-BandProp[bandname][0])

				#=======plotting right column==============
				print('treating {} of {}'.format(tag,name))
				timeshift=time[np.where(np.nanmin(Mean)==Mean)[0][0]]-time[np.where(np.nanmax(MeanForce)==MeanForce)[0][0]]
				timeshift=12.*timeshift
				if j==0 and k==0:
					col2a=ax[1,plotcol-1].plot(int(tag[-1]), timeshift, marker='*',color='k',label='Effective pressure lag',ms=10,mew=2,linewidth=0)[0]
					col2b=ax2.plot(int(tag[-1]),np.nanmax(Mean)-np.nanmin(Mean), marker='+',color='b',label='Effective pressure amplitude',ms=10,mew=2,linewidth=0)[0]
				else:
					ax[k+1,plotcol-1].plot(int(tag[-1]), timeshift, marker='*',color='k',ms=10,mew=2)
					ax2.plot(int(tag[-1]),np.nanmax(Mean)-np.nanmin(Mean), marker='+',color='b',ms=10,mew=2,linewidth=0)

				#==============setting limits================
				ax[k+1,plotcol-1].set_ylim([-1.8,1.44])
				ax[k+1,plotcol-1].set_xticks([1,2,3,4,5])
				ax[plotrow-1,plotcol-1].xaxis.set_ticklabels(['D1','D2','D3','D4','D5'])
						
				for frame in framelist:
					if j==frame:
						#========patch on right to show what we have on left======
						rect = [Rectangle((frame+0.5,-1.8), 1.0, 3.24)]
						pc = PatchCollection(rect, facecolor=(0.8,0.8,0.8), alpha=1,
																 edgecolor=(0.8,0.8,0.8),zorder=0)
						ax[k+1,plotcol-1].add_collection(pc)

						#================setting up left plots====================
						for count in[0,1]:
							ax[k+1,1].set_title(namestyles[name][0],va='top',ha='right',position=(0.97,0.8))
							ax[k+1,0].set_title(namestyles[name][0],va='center',ha='right',position=(0.97,0.1))

							ax[k+1,count].set_position([hleft+count*(plotwidth+hspace),vbot+(plotrow-2-k)*(plotheight+vspace),plotwidth,plotheight])
							ax[k+1,count].tick_params(axis='y', which='major', pad=2,top=False)
							ax[k+1,count].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
							ax[k,count].xaxis.set_ticklabels([])
							ax[k+1,count].set_xlim([0.25 ,1])
							ax[k+1,count].set_xticks([0.25,0.5,0.75,1.0])
							ax[k+1,0].set_ylim([-5,10])
							ax[k+1,1].set_ylim([0,1200])
							ax[0,0].set_ylim([0,1200])
							ax[0,1].set_ylim([0,1200])
							
							if k==2:
								ax[0,count].set_position([hleft+count*(plotwidth+hspace),vbot+(plotrow-1)*(plotheight+vspace),plotwidth,plotheight])
								ax[0,count].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
								ax[0,count].set_xlim([0.25 ,1])
								ax[0,count].set_xticks([0.25,0.5,0.75,1.0])
								ax[plotrow-1,count].xaxis.set_ticklabels(['$0.25$','$0.5$','$0.75$','$1.0$'])
								ax[plotrow-1,count].set_xlabel('Time in years',labelpad=0)
								ax[plotrow-1,count].tick_params(axis='x', which='major', pad=5)
								if count==0:
									topleft=ax[0,0].plot(time, MeanForce, linestyle='-',linewidth=2.,color='k',label='Total recharge')
								else:
									topright=[]
									for b,bandname in enumerate(['Low','Mid','Top']):
										topright.append(ax[0,1].plot(time, Force[:,b], linestyle='-',linewidth=2.,color=BandProp[bandname][2],label='Recharge upstream of '+BandProp[bandname][4])[0])

						channelflux,sheetflux,lines=GetFluxRatio(name,tag,time)
						totflux=channelflux+sheetflux
						print('for {} in {} max channelflux is {}'.format(name,tag,np.nanmax(channelflux)))
						print('for {} in {} max sheetflux is {}'.format(name,tag,np.nanmax(sheetflux)))
						if k==0:
							pressure=[]
							IDS=[]
							EDS=[]
						if name=='db':
							pressure.append(ax[k+1,0].plot(time, Mean[:], linestyle='-',linewidth=2.,color='g',label='Mean effective Pressure')[0])
							ax[k+1,0].plot(time, np.zeros(np.shape(time)), linestyle=':',linewidth=2.,color='k',)
							IDS.append(ax[k+1,1].plot(time, np.nanmean(sheetflux,axis=1), linestyle='--',linewidth=2.,color='g',label='Mean inefficient syst. flux')[0])
							EDS.append(ax[k+1,1].plot(time, np.nanmean(channelflux,axis=1), linestyle=':',linewidth=2.,color='g',label='Mean efficient syst. flux')[0])
							print np.shape(pressure)
						else:
							if name=='mw':
								pressure.append(ax[1,0].plot(time, Mean[:], linestyle='-.',linewidth=2.,color='g',label='Mean effective Pressure for $mw$')[0])
							for b,bandname in enumerate(['Low','Mid','Top']):
								ax[k+1,0].plot(time, Band[:,b,1], linestyle='-',linewidth=2.,color=BandProp[bandname][2])
								if k==1:
									pressure.append(ax[k+1,0].fill_between(time,Band[:,b,0],Band[:,b,2],
																												 label='Effective pressure for '+BandProp[bandname][4],
																												 facecolor=BandProp[bandname][2],edgecolor=BandProp[bandname][2],alpha=0.3))
								else:
									ax[k+1,0].fill_between(time,Band[:,b,0],Band[:,b,2],
																				 facecolor=BandProp[bandname][2],edgecolor=BandProp[bandname][2],alpha=0.3)
								bandloc=np.where(lines==BandProp[bandname][0])[0]
								if k==1:
									IDS.append(ax[2,1].plot(time, sheetflux[:,bandloc], linestyle='--',linewidth=2.,color=BandProp[bandname][2],label='Inef. syst. flux in '+BandProp[bandname][4])[0])
									EDS.append(ax[2,1].plot(time, channelflux[:,bandloc], linestyle=':',linewidth=2.,color=BandProp[bandname][2],label='Eff. syst. flux in '+BandProp[bandname][4])[0])
								else:
									ax[k+1,1].plot(time, sheetflux[:,bandloc], linestyle='--',linewidth=2.,color=BandProp[bandname][2])
									ax[k+1,1].plot(time, channelflux[:,bandloc], linestyle=':',linewidth=2.,color=BandProp[bandname][2])
							ax[k+1,0].plot(time, np.zeros(np.shape(time)), linestyle=':',linewidth=2.,color='k',)
	layout.delaxes(ax[0,-1])

	#================adding axis labels====================
	ax[-1,-1].text(0.015,1-vtop-0.5*plotheight,'Recharge\n (m$^3$s$^{-1}$)',va='center',ha='center',rotation=90,transform=layout.transFigure,color='r')
	ax[-1,-1].text(0.015,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Effective Pressure (MPa)',va='center',ha='center',rotation=90,transform=layout.transFigure)
	ax[-1,-1].text(hleft+plotwidth+0.25*hspace,1-vtop-0.5*plotheight,'Recharge',va='center',ha='center',rotation=90,transform=layout.transFigure,color='r')
	ax[-1,-1].text(hleft+plotwidth+0.25*hspace,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Discharge (m$^3$s$^{-1}$)',va='center',ha='center',rotation=90,transform=layout.transFigure)
	ax[-1,-1].text(hleft+2.5*plotwidth+2.5*hspace,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Amplitude (MPa)',color='b',va='center',ha='center',rotation=90,transform=layout.transFigure)
	ax[-1,-1].text(hleft+2*plotwidth+1.25*hspace,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Time lag (month)',va='center',ha='center',rotation=90,transform=layout.transFigure)

	L=plt.legend(handles=[pressure[1],pressure[2],pressure[3],pressure[0],pressure[4],topleft[0],topright[0],topright[1],topright[2],IDS[0],IDS[1],IDS[2],IDS[3],EDS[0],EDS[1],EDS[2],EDS[3],col2a,col2b],
							 bbox_to_anchor=(1+hspace-hright,vbot,0.88*hright,1.-vbot),numpoints=1,loc=3,
							 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.3,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,
							 borderaxespad=0,handler_map={pressure[1]:FillandLineHandler('b'),pressure[2]:FillandLineHandler('c'),pressure[3]:FillandLineHandler('r')})
	plt.setp(L.texts, family='Arial')
	plt.gca().add_artist(L)

	plt.show
	return layout

# }}}
# {{{ BandEvolOneCol()

def BandEvolOneCol():
	width=12.
	height=9.
	plotrow=8
	plotcol=1
	layout,ax=plt.subplots(plotrow,plotcol,figsize=(width,height))
	Exptag='D'
	hleft=0.1
	hright=0.5
	vside=0.05
	vbot=1.8*vside
	vtop=0.2*vside
	vspace=0.05/height
	plotwidth=(1-hleft-hright)
	plotheight=(1.0-(2.*vside)-(vspace*(plotrow-1)))/plotrow

	k=-1

	AmpList={'D1':-4.0,'D2':-2.0,'D3':0.0,'D4':2.0,'D5':4.0,}
	BandProp={'Low':[10000.,600.,'b','$Low$','lower band'],
						'Mid':[50000.,3000.,'c','$Mid$','middle band'],
						'Top':[90000.,5100.,'r','$Top$','highest band']}
					
	pressure=[]
	IDS=[]
	EDS=[]
	for i,name in enumerate(names):
		if len(glob(resdir+name+'/'+expe+'_'+name+'.nc'))>0:
			k+=1
			#==== Core computation==============
			NCFile=resdir+name+'/'+expe+'_'+name+'.nc'
			try:
				DatFile	 = Dataset(NCFile, mode='r')
				time     = DatFile.variables['time'][:]
				data,coords=GetVar(DatFile,DataName,name,True)
				if name=='mw':
					thick,coords=GetVar(DatFile,'H',name,False)
					base,coords=GetVar(DatFile,'B',name,False)
				DatFile.close()
			except (RuntimeError,IOError):
				print ('File {} does not exist'.format(NCFile))
				continue 
			time=1+(time-max(time))/(365.*24.*3600.)
			time=time[np.where(time>=0)][:]
			if name=='mw':
				surf=thick+base
				Forcing=SeasonForcing(surf,time*3600.*24.*365.,AmpList[expe])

			Band=np.NaN*np.ones((len(time),3,3))
			Mean=np.NaN*np.ones((len(time)))
			Force=np.NaN*np.ones((len(time),3))
			Meanwidth=np.NaN*np.ones((3))
			nodeselect=[]
			if name=='db':
				for t,step in enumerate(time):
					Mean[t]=np.nanmean(data[t,:])
			else:
				for b,bandname in enumerate(['Low','Mid','Top']):
					nodeselect.append(np.where(np.logical_and(coords[0,:]>=BandProp[bandname][0],
																										coords[0,:]<=BandProp[bandname][0]+5000))[0])
				for t,step in enumerate(time):
					for b,bandname in enumerate(['Low','Mid','Top']):
						Band[t,b,0]=np.nanmin(data[t,nodeselect[b]])
						Band[t,b,1]=np.nanmean(data[t,nodeselect[b]])
						Band[t,b,2]=np.nanmax(data[t,nodeselect[b]])
						Mean[t]=np.nanmean(data[t,:])
						if name=='mw':
							upstream=np.where(coords[0,:]>=BandProp[bandname][0])[0]
							Force[t,b]=np.nanmean(Forcing[upstream,t])*20.0e3*(100.0e3-BandProp[bandname][0])

			#================setting up left plots====================
			ax[k+1].set_title(namestyles[name][0],va='center',ha='right',position=(0.97,0.1))
			
			ax[k+1].set_position([hleft,vbot+(plotrow-2-k)*(plotheight+vspace),plotwidth,plotheight])
			ax[k+1].tick_params(axis='y', which='major', pad=2,top=False)
			ax[k+1].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
			ax[k].xaxis.set_ticklabels([])
			ax[k+1].set_xlim([0.25 ,1])
			ax[k+1].set_xticks([0.25,0.5,0.75,1.0])
#			ax[k+1].set_ylim([-5,10])
			ax[0].set_ylim([0,1200])

			if k==0:
				ax[0].set_position([hleft,vbot+(plotrow-1)*(plotheight+vspace),plotwidth,plotheight])
				ax[0].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
				ax[0].set_xlim([0.25 ,1])
				ax[0].set_xticks([0.25,0.5,0.75,1.0])
				ax[plotrow-1].xaxis.set_ticklabels(['$0.25$','$0.5$','$0.75$','$1.0$'])
				ax[plotrow-1].set_xlabel('Time in years',labelpad=0)
				ax[plotrow-1].tick_params(axis='x', which='major', pad=5)
			if name=='db':
				pressure.append(ax[k+1].plot(time, Mean[:], linestyle='-',linewidth=2.,color='g',label='Mean effective Pressure')[0])
			else:
				if name=='mw':
					pressure.append(ax[1].plot(time, Mean[:], linestyle='-.',linewidth=2.,color='g',label='Mean effective Pressure for $mw$')[0])
					topright=[]
					for b,bandname in enumerate(['Low','Mid','Top']):
						topright.append(ax[0].plot(time, Force[:,b], linestyle='-',linewidth=2.,color=BandProp[bandname][2],label='Recharge upstream of '+BandProp[bandname][4])[0])
				for b,bandname in enumerate(['Low','Mid','Top']):
					ax[k+1].plot(time, Band[:,b,1], linestyle='-',linewidth=2.,color=BandProp[bandname][2])
					if k==1:
						pressure.append(ax[k+1].fill_between(time,Band[:,b,0],Band[:,b,2],
																								 label='Effective pressure for '+BandProp[bandname][4],
																								 facecolor=BandProp[bandname][2],edgecolor=BandProp[bandname][2],alpha=0.3))
					else:
						ax[k+1].fill_between(time,Band[:,b,0],Band[:,b,2],
																 facecolor=BandProp[bandname][2],edgecolor=BandProp[bandname][2],alpha=0.3)
			print ('plotting 0 line')
			ax[k+1].plot(time, np.zeros(np.shape(time)), linestyle=':',linewidth=2.,color='k',)
	#================adding axis labels====================
	ax[-1].text(0.03,1-vtop-0.5*plotheight,'Recharge\n (m$^3$s$^{-1}$)',va='center',ha='center',rotation=90,transform=layout.transFigure,color='r')
	ax[-1].text(0.03,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Effective Pressure (MPa)',va='center',ha='center',rotation=90,transform=layout.transFigure)

	L=plt.legend(handles=[pressure[1],pressure[2],pressure[3],pressure[0],pressure[4],topright[0],topright[1],topright[2]],
							 bbox_to_anchor=(hleft+plotwidth+0.05,0.66,0.88*hright,1.-vbot),numpoints=1,loc=3,
							 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.3,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,
							 borderaxespad=0,handler_map={pressure[1]:FillandLineHandler('b'),pressure[2]:FillandLineHandler('c'),pressure[3]:FillandLineHandler('r')})
	plt.setp(L.texts, family='Arial')
	plt.gca().add_artist(L)

	plt.show
	return layout

# }}}
# {{{ BandEvolChannel()

def BandEvolChannel():
	width=12.
	height=9.
	plotrow=8
	plotcol=1
	layout,ax=plt.subplots(plotrow,plotcol,figsize=(width,height))
	Exptag='D'
	hleft=0.1
	hright=0.5
	vside=0.05
	vbot=1.8*vside
	vtop=0.2*vside
	vspace=0.05/height
	plotwidth=(1-hleft-hright)
	plotheight=(1.0-(2.*vside)-(vspace*(plotrow-1)))/plotrow

	k=-1

	AmpList={'D1':-4.0,'D2':-2.0,'D3':0.0,'D4':2.0,'D5':4.0,}
	BandProp={'Low':[10000.,600.,'b','$Low$','lower band'],
						'Mid':[50000.,3000.,'c','$Mid$','middle band'],
						'Top':[90000.,5100.,'r','$Top$','highest band']}
					
	pressure=[]
	IDS=[]
	EDS=[]
	for i,name in enumerate(names):
		if len(glob(resdir+name+'/'+expe+'_'+name+'.nc'))>0:
			k+=1
			#==== Core computation==============
			NCFile=resdir+name+'/'+expe+'_'+name+'.nc'
			try:
				DatFile	 = Dataset(NCFile, mode='r')
				time     = DatFile.variables['time'][:]
				data,coords=GetVar(DatFile,'S',name,True)
				if name=='mw':
					thick,coords=GetVar(DatFile,'H',name,False)
					base,coords=GetVar(DatFile,'B',name,False)
				DatFile.close()
			except (RuntimeError,IOError):
				print ('File {} does not exist'.format(NCFile))
				continue 
			time=1+(time-max(time))/(365.*24.*3600.)
			time=time[np.where(time>=0)][:]
			if name=='mw':
				surf=thick+base
				Forcing=SeasonForcing(surf,time*3600.*24.*365.,AmpList[expe])

			Band=np.NaN*np.ones((len(time),3,3))
			Mean=np.NaN*np.ones((len(time)))
			Force=np.NaN*np.ones((len(time),3))
			Meanwidth=np.NaN*np.ones((3))
			nodeselect=[]
			for b,bandname in enumerate(['Low','Mid','Top']):
				nodeselect.append(np.where(np.logical_and(coords[0,:]>=BandProp[bandname][0],
																									coords[0,:]<=BandProp[bandname][0]+5000))[0])
			for t,step in enumerate(time):
				for b,bandname in enumerate(['Low','Mid','Top']):
					Band[t,b,0]=np.nanmin(data[t,nodeselect[b]])
					Band[t,b,1]=np.nanmean(data[t,nodeselect[b]])
					Band[t,b,2]=np.nanmax(data[t,nodeselect[b]])
					Mean[t]=np.nanmean(data[t,:])
					if name=='mw':
						upstream=np.where(coords[0,:]>=BandProp[bandname][0])[0]
						Force[t,b]=np.nanmean(Forcing[upstream,t])*20.0e3*(100.0e3-BandProp[bandname][0])

			#================setting up left plots====================
			ax[k+1].set_title(namestyles[name][0],va='center',ha='right',position=(0.97,0.1))
			
			ax[k+1].set_position([hleft,vbot+(plotrow-2-k)*(plotheight+vspace),plotwidth,plotheight])
			ax[k+1].tick_params(axis='y', which='major', pad=2,top=False)
			ax[k+1].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
			ax[k].xaxis.set_ticklabels([])
			ax[k+1].set_xlim([0.25 ,1])
			ax[k+1].set_xticks([0.25,0.5,0.75,1.0])
			ax[k+1].set_ylim([-5,105])
			ax[0].set_ylim([0,1200])

			if k==0:
				ax[0].set_position([hleft,vbot+(plotrow-1)*(plotheight+vspace),plotwidth,plotheight])
				ax[0].yaxis.set_major_locator(MaxNLocator(5,prune='upper'))
				ax[0].set_xlim([0.25 ,1])
				ax[0].set_xticks([0.25,0.5,0.75,1.0])
				ax[plotrow-1].xaxis.set_ticklabels(['$0.25$','$0.5$','$0.75$','$1.0$'])
				ax[plotrow-1].set_xlabel('Time in years',labelpad=0)
				ax[plotrow-1].tick_params(axis='x', which='major', pad=5)
				IDS=[]
				EDS=[]
				Rat=[]

			channelflux,sheetflux,lines=GetFluxRatio(name,expe,time)
			totflux=channelflux+sheetflux
			ratio=100.*channelflux/totflux
			print('for {} in {} max channelflux is {}'.format(name,expe,np.nanmax(channelflux)))
			print('for {} in {} max sheetflux is {}'.format(name,expe,np.nanmax(sheetflux)))

			if name=='db':
				# IDS.append(ax[k+1].plot(time, np.nanmean(sheetflux,axis=1), linestyle='--',linewidth=2.,color='g',label='Mean inefficient syst. flux')[0])
				# EDS.append(ax[k+1].plot(time, np.nanmean(channelflux,axis=1), linestyle=':',linewidth=2.,color='g',label='Mean efficient syst. flux')[0])
				Rat.append(ax[k+1].plot(time, np.nanmean(ratio,axis=1), linestyle='-',linewidth=2.,color='g',label='Channelisation ratioin '+BandProp[bandname][4])[0])
			else:
				if name=='mw':
					topright=[]
					for b,bandname in enumerate(['Low','Mid','Top']):
						topright.append(ax[0].plot(time, Force[:,b], linestyle=':',linewidth=2.,color=BandProp[bandname][2],label='Recharge upstream of '+BandProp[bandname][4])[0])
				for b,bandname in enumerate(['Low','Mid','Top']):
					bandloc=np.where(lines==BandProp[bandname][0])[0]
					if k==1:
						# IDS.append(ax[2].plot(time, sheetflux[:,bandloc], linestyle='--',linewidth=2.,color=BandProp[bandname][2],label='Inef. syst. flux in '+BandProp[bandname][4])[0])
						# EDS.append(ax[2].plot(time, channelflux[:,bandloc], linestyle=':',linewidth=2.,color=BandProp[bandname][2],label='Eff. syst. flux in '+BandProp[bandname][4])[0])
						Rat.append(ax[2].plot(time, ratio[:,bandloc], linestyle='-',linewidth=2.,color=BandProp[bandname][2],label='Channelisation ratio in '+BandProp[bandname][4])[0])
					else:
						# ax[k+1].plot(time, sheetflux[:,bandloc], linestyle='--',linewidth=2.,color=BandProp[bandname][2])
						# ax[k+1].plot(time, channelflux[:,bandloc], linestyle=':',linewidth=2.,color=BandProp[bandname][2])
						ax[k+1].plot(time, ratio[:,bandloc], linestyle='-',linewidth=2.,color=BandProp[bandname][2])

	#================adding axis labels====================
	ax[-1].text(0.03,1-vtop-0.5*plotheight,'Recharge\n (m$^3$s$^{-1}$)',va='center',ha='center',rotation=90,transform=layout.transFigure,color='r')
	ax[-1].text(0.03,vbot+plotheight*((k+1)/2.)+((k)/2.)*vspace,'Channelisation ratio \%',va='center',ha='center',rotation=90,transform=layout.transFigure)

	L=plt.legend(handles=[topright[0],topright[1],topright[2],Rat[0],Rat[1],Rat[2]],
							 bbox_to_anchor=(hleft+plotwidth+0.05,0.66,0.88*hright,1.-vbot),numpoints=1,loc=3,
							 bbox_transform=layout.transFigure,ncol=1,labelspacing=0.3,columnspacing=0.2,borderpad=0.3,handletextpad=0.2,
							 borderaxespad=0,)
	plt.setp(L.texts, family='Arial')
	plt.gca().add_artist(L)

	plt.show
	return layout

# }}}
# {{{ BandVal()

def BandVal():

	layout,ax=NameLayout()
	indexes=np.arange(len(plotexpe))+1
	k=-1
	
	for i,name in enumerate(names):
		if len(glob(resdir+name+'/A*_'+name+'.nc'))>0:
			print ('treating for {}'.format(name))
			LowBand=np.NaN*np.ones((len(plotexpe),3))
			MidBand=np.NaN*np.ones((len(plotexpe),3))
			TopBand=np.NaN*np.ones((len(plotexpe),3))
			DomainMean=np.NaN*np.ones((len(plotexpe)))
			k+=1
			for j,tag in enumerate(plotexpe):
				NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
				try:
					DatFile	 = Dataset(NCFile, mode='r')
					data,coords=GetVar(DatFile,DataName,name,False)
					DatFile.close()
				except (RuntimeError,IOError):
					print ('File {} does not exist'.format(NCFile))
					continue 

				griddeddata,gridpoint_x,gridpoint_y=Griding(data,coords)
				DomainMean[j]=np.nanmean(data)
				if name!='db':
					lowloc=np.where(np.logical_and(gridpoint_x>=10,gridpoint_x<=15))
					midloc=np.where(np.logical_and(gridpoint_x>=50,gridpoint_x<=55))
					toploc=np.where(np.logical_and(gridpoint_x>=90,gridpoint_x<=95))
					LowBand[j,0]=np.nanmin(griddeddata[lowloc])
					MidBand[j,0]=np.nanmin(griddeddata[midloc])
					TopBand[j,0]=np.nanmin(griddeddata[toploc])
					LowBand[j,1]=np.nanmean(griddeddata[lowloc])
					MidBand[j,1]=np.nanmean(griddeddata[midloc])
					TopBand[j,1]=np.nanmean(griddeddata[toploc])
					LowBand[j,2]=np.nanmax(griddeddata[lowloc])
					MidBand[j,2]=np.nanmax(griddeddata[midloc])
					TopBand[j,2]=np.nanmax(griddeddata[toploc])

			if k==1:
				bands=[]
				ax[k].plot(indexes,LowBand[:,1],color='b',marker='.',linestyle='-',markersize=7,mew=2,linewidth=2)
				bands.append(ax[k].fill_between(indexes,LowBand[:,0],LowBand[:,2],facecolor='b',alpha=0.3,edgecolor='b',label='Lower band'))
				ax[k].plot(indexes,MidBand[:,1],color='c',marker='.',linestyle='-',markersize=7,mew=2,linewidth=2)
				bands.append(ax[k].fill_between(indexes,MidBand[:,0],MidBand[:,2],facecolor='c',alpha=0.3,edgecolor='c',label='Middle band'))
				ax[k].plot(indexes,TopBand[:,1],color='r',marker='.',linestyle='-',markersize=7,mew=2,linewidth=2)
				bands.append(ax[k].fill_between(indexes,TopBand[:,0],TopBand[:,2],facecolor='r',alpha=0.3,edgecolor='r',label='Highest band'))
				bands.append(ax[k].plot(indexes,DomainMean,color='g',marker='+',linestyle='-',markersize=5,mew=2,label='Domain mean',linewidth=2)[0])
			else:
				ax[k].plot(indexes,DomainMean,color='g',marker='+',linestyle='-',markersize=7,mew=2,linewidth=2)
				if name=='mw':
					ax[0].plot(indexes,DomainMean,color='g',linestyle=':',markersize=7,mew=2,linewidth=2)
				if name!='db':
					ax[k].plot(indexes,LowBand[:,1],color='b',marker='.',linestyle='-',markersize=7,mew=2,linewidth=2)
					ax[k].fill_between(indexes,LowBand[:,0],LowBand[:,2],facecolor='b',alpha=0.3,edgecolor='b')
					ax[k].plot(indexes,MidBand[:,1],color='c',marker='.',linestyle='-',markersize=7,mew=2,linewidth=2)
					ax[k].fill_between(indexes,MidBand[:,0],MidBand[:,2],facecolor='c',alpha=0.3,edgecolor='c')
					ax[k].plot(indexes,TopBand[:,1],color='r',marker='.',linestyle='-',markersize=7,mew=2,linewidth=2)
					ax[k].fill_between(indexes,TopBand[:,0],TopBand[:,2],facecolor='r',alpha=0.3,edgecolor='r')
			
			
	L=plt.legend(handles=[bands[2],bands[1],bands[0],bands[3]],bbox_to_anchor=(10.5/14,0.6,3./14,0.3),bbox_transform=layout.transFigure,
							 mode='expand',loc=2,borderaxespad=0,fontsize=sizeoffont,
							 handler_map={bands[0]:FillandLineHandler('b'),bands[1]:FillandLineHandler('c'),bands[2]:FillandLineHandler('r')})
	plt.setp(L.texts, family='Arial')
	plt.gca().add_artist(L)

	plt.show
	return layout

# }}}
# {{{ Mapping()

def Mapping():
	layout,ax=MapLayout()

	#defining colormap
	if DataName=='N':
		norm=mpl.colors.Normalize(vmin=0, vmax=2)
	elif DataName=='Q':
		norm=mpl.colors.LogNorm(vmin=1.0e-6, vmax=30)
	else:
		norm=mpl.colors.Normalize(vmin=0, vmax=10)

	cmap = plt.cm.viridis
	k=-1
	for i,name in enumerate(names):
		NCFile=resdir+name+'/'+plotexpe[0]+'_'+name+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			data,coords=GetVar(DatFile,DataName,name,False)
			DatFile.close()
			k+=1
		except (RuntimeError,IOError):
			continue

		print('for {} min {} is {} and max is {}'.format(name,DataName,np.nanmin(abs(data)),np.nanmax(abs(data))))
		if np.all(coords[1,:]==10*1.0e3):
			if name in ['mh1','cdf']:
				data=np.flipud(data)
			if len(np.shape(data))>1:
				data=np.nanmean(data,axis=0)
			ax[k].imshow(np.asarray([data]).T, extent=(0,20,0,100),
									 origin='lower',norm=norm,cmap=cmap,aspect='auto')
			ax[k].imshow(np.zeros(np.shape((np.asarray([data]).T))), extent=(0,20,0,100),
									 origin='lower',norm=norm,cmap=cmap,aspect='auto')
			ax[k].imshow(np.asarray([data]).T, extent=(9.9,10.1,0,100),
									 origin='lower',norm=norm,cmap=cmap,aspect='auto')
		else:
			griddeddata,gridpoint_x,gridpoint_y=Griding(data,coords)
			ax[k].imshow(abs(griddeddata), extent=(0,20,0,100),
									 origin='lower',norm=norm,cmap=cmap,aspect='auto')

		if k==0:
			if DataName=='N':
				ax[0].annotate('Effective pressure maps for simulation '+plotexpe[0]+' (in MPa)',xy=(0.5,0.96),va='center',ha='center',
											 textcoords='figure fraction',xycoords='figure fraction',fontsize=18)
			elif DataName=='S':
				ax[0].annotate('Efficient system cross sectional area for simulation '+plotexpe[0]+' (in m$^2$)',xy=(0.5,0.96),va='center',ha='center',
											 textcoords='figure fraction',xycoords='figure fraction',fontsize=18)

		ax2 = layout.add_axes([(9.4/10.),0.12,0.2/10, 0.75])
		cbar = mpl.colorbar.ColorbarBase(ax2, norm=norm,cmap=cmap, orientation='vertical')
		cbar.solids.set_edgecolor("face")
		plt.show
	return layout

# }}}
# {{{ InputCurve()
def InputCurve():
	fig = plt.figure(figsize=(6,3))
	ax = fig.add_subplot(111)

	AmpList={'D1':[-4.0,(0,0,0)],
					 'D2':[-2.0,(0.2,0.2,0.2)],
					 'D3':[0.0,(0.4,0.4,0.4)],
					 'D4':[2.0,'r'],
					 'D5':[4.0,(0.6,0.6,0.6)]}
	for tag in plotexpe:
		NCFile=resdir+'mw/'+tag+'_mw.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			time     = DatFile.variables['time'][:]
			thick,coords=GetVar(DatFile,'H','mw',False)
			base,coords=GetVar(DatFile,'B','mw',False)
			DatFile.close()
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))
			continue 
		time=1+(time-max(time))/(365.*24.*3600.)
		time=time[np.where(time>=0)][:]
		surf=thick+base
		
		Forcing=SeasonForcing(surf,time*3600.*24.*365.,AmpList[tag][0])
		MeanForce=np.nanmean(Forcing,axis=0)*100.0e3*20.0e3
		ax.set_position([0.15,0.22,0.65,0.7])
		ax.plot(time, MeanForce, linestyle='-',linewidth=2.,color=AmpList[tag][1],label=tag)
		ax.set_xlim([0.0 ,1])
		ax.set_xticks([0.0,0.25,0.5,0.75,1.0])
		ax.set_ylabel('Recharge (m$^3$s$^{-1}$)')
		ax.set_xlabel('Time in years')
		ax.legend(loc=(1.02,0),labelspacing=0.12,columnspacing=0.2,borderpad=0.2,handletextpad=0.1,
							 borderaxespad=0)

	return fig
# }}}
# {{{ OneSection()
def OneSection():
	layout,ay=SectionLayout()
	k=-1
	for i,name in enumerate(names):
		NCFile=resdir+name+'/A1_'+name+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			DatFile.variables[DataName].dimensions
			data,coords=GetVar(DatFile,DataName,name,False)
			DatFile.close()
			k+=1
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))
			continue 

		griddeddata,gridpoint_x,gridpoint_y=Griding(data,coords)
		meandata=np.nanmean(griddeddata,axis=1)
		ay[k/3,np.mod(k,3)].plot(gridpoint_x[:,0], meandata, linestyle='-',linewidth=2,color='r')
		
	return layout
# }}}
# {{{ GetVar(DatFile,var,name,AllTime):

def GetVar(DatFile,var,name,AllTime):
	#getting data array
	try:
		dimension=np.size(DatFile.variables[var].dimensions)
		locvar=var
	except KeyError:
		if name=='bf' and var=='S':
			dimension=np.size(DatFile.variables['Ee'].dimensions)
			locvar='Ee'
		else:
			dimension=np.size(DatFile.variables['N'].dimensions)
			locvar='N'
	if dimension==1:
		data = DatFile.variables[locvar][:]
	else:
		if AllTime:
			data = DatFile.variables[locvar][:,:]
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
			if DatFile.variables[locvar].dimensions[-1]=='time':
				data = DatFile.variables[locvar][:,-1]
			else:
				data = DatFile.variables[locvar][-1,:]
			if name=='id' and var=='N':
				data[-1]=data[-2]
	if var=='N':
		data=1.0e-6*data
	elif locvar!=var and locvar!='Ee':
		data=np.nan*data
	elif locvar=='Ee':
		data=1000.*data

	for key in DatFile.variables[locvar].dimensions:
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
	return data, coords

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

	if name in['mh1','mh2','mw','mw_prime','og','og_2','cdf','db','sb','jd','jsb','as','as_2','id','rh','bf']:
		NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			if name in['mh2','mw','mw_prime','og','og_2','mh1','id']:
				if name!='id':
					IneffFlux,Ineffcoords=GetVar(DatFile,'q',name,TimeVar)
				EffFlux,Effcoords=GetVar(DatFile,'Q',name,TimeVar)
				if name=='mh1':
					Xs=DatFile.variables['coords3'][0,:]
				else:
					Xs=DatFile.variables['coords1'][0,:]
				connect=DatFile.variables['connect_ch'][:,:]
				intconnect=np.asarray(connect,dtype=int)
				DatFile.close()
			elif name in ['cdf','rh']:
				if name=='cdf':
					IneffFlux,Ineffcoords=GetVar(DatFile,'q',name,TimeVar)
				EffFlux,Effcoords=GetVar(DatFile,'Q',name,TimeVar)
				DatFile.close()
			elif name in['jd','jsb']:
				IneffFlux,Ineffcoords=GetVar(DatFile,'q',name,TimeVar)
				DatFile.close()
			elif name in ['as','as_2']:
				IneffFlux,Ineffcoords=GetVar(DatFile,'q',name,TimeVar)
				Channelisation,Effcoords=GetVar(DatFile,'Dc',name,TimeVar)
				EffFlux=IneffFlux*Channelisation
				IneffFlux=IneffFlux*(1-Channelisation)
				DatFile.close()
			elif name=='sb':
				IneffFlux,Ineffcoords=GetVar(DatFile,'q',name,TimeVar)
				Conductivity,Ineffcoords=GetVar(DatFile,'K',name,TimeVar)
				ischannel=np.ones(np.shape(Conductivity))
				ischannel[np.where(Conductivity==np.nanmin(Conductivity))]=0
				EffFlux=IneffFlux*ischannel
				IneffFlux=IneffFlux*(1-ischannel)
				DatFile.close()
			elif name=='bf':
				IneffFlux,Ineffcoords=GetVar(DatFile,'q',name,TimeVar)
				EffFlux,Effcoords=GetVar(DatFile,'Q',name,TimeVar)
				Thickness,Effcoords=GetVar(DatFile,'Ee',name,TimeVar)
				isactive=np.zeros(np.shape(EffFlux))
				isactive[np.where(EffFlux>0)]=1.
				DatFile.close()
			elif name=='db':
				IneffFlux,Ineffcoords=GetVar(DatFile,'q',name,TimeVar)
				Channelisation,Effcoords=GetVar(DatFile,'R',name,TimeVar)
				EffFlux=IneffFlux*Channelisation
				IneffFlux=IneffFlux*(1-Channelisation)
				DatFile.close()
		except (RuntimeError,IOError):
			print ('File {} does not exist'.format(NCFile))

		for i,loc in enumerate(Lines):
			#===mh2=mw=og===
			if name in['mh2','mw','mw_prime','og','og_2','mh1']:
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
						ChannelFlux[:,i]=np.nanmean(abs(EffFlux[:,channelloc]),axis=1)*20.
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*Width[i]
					if len(channelloc)>0:
						ChannelFlux[i]=np.nanmean(abs(EffFlux[channelloc]))*20.
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
				else:
					if len(channelloc)>0:
						ChannelFlux[i]=np.sum(abs(EffFlux[channelloc])).T
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
			#===jsb=db===
			elif name in ['jsb','db']:
				sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																				 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				if TimeVar:
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*np.nanmean(Width)
					if name=='db':
						ChannelFlux[:,i]=np.nanmean(EffFlux[:,sheetloc],axis=1)*np.nanmean(Width)
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*np.nanmean(Width)
					if name=='db':
						ChannelFlux[i]=np.nanmean(EffFlux[sheetloc])*np.nanmean(Width)
			#===sb===
			elif name=='sb':
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
		ax[i].set_ylabel(expe)
	return fig,ax

# }}}
# {{{ SliceLayout(frames,Exp) 

def SliceLayout(frames):

	width=17.8/2.54
	fig = plt.figure(figsize=(width,np.ceil(frames/5.)/5.*width))
	ax=[]
	fig.subplots_adjust(left=0.04, right=1, top=1.03, bottom=0.03, wspace=0.05, hspace=-0.10)
	for i in np.arange(0,frames):
		ax.append(fig.add_subplot(int(np.ceil(frames/5.)), 5, i+1, projection='3d'))
		ax[i].patch.set_visible(False) 
		if i in [8,9,10,11,12]:
			ax[i].set_xlabel('$x$ (km)',labelpad=-11)
		ax[i].tick_params(axis='x', which='major', pad=0)
		ax[i].set_xlim3d(0, 100)
		ax[i].set_xticks([10,50,90])
		ax[i].set_xticklabels(['$10$','$50$','$90$'],ha='center',va='bottom',rotation=10)

		ax[i].set_ylim3d(0,len(plotexpe))
		ax[i].set_ylim3d(-0.1,5.1)
		ax[i].set_yticks([0,1,2,3,4,5])
		ax[i].set_yticklabels(['','$2$','','$4$','','$6$'],ha='right',va='center')#
		if i in [8,9,10,11,12]:
			ax[i].set_ylabel('A',labelpad=-10)
		ax[i].tick_params(axis='y', which='major', pad=-5)

		ax[i].zaxis.set_rotate_label(False)
		if i in [0,5,10]:
			ax[i].set_zlabel('N (MPa)',labelpad=0,rotation=90)
		ax[i].set_zlim3d(0, 14)
		ax[i].set_zticks([0,2,4,6,8,10,12,14])
		ax[i].set_zticklabels(['$0$','','$4$','','$8$','','$12$',''],ha='right',va='center')
		ax[i].tick_params(axis='z', which='major', pad=-3)

		ax[i].set_frame_on(False)

	plt.show
	return fig,ax

# }}}
# {{{ NameLayout() 

def NameLayout():

	fig,ax=plt.subplots(1,13,sharey='row',figsize=(13,4))
	hleft=0.5/14.
	width=1./14.
	vside=0.1
	height=0.8
	k=-1

	for i,name in enumerate(names):
		if len(glob(resdir+name+'/A*_'+name+'.nc'))>0:
			k+=1
			ax[k].set_position([hleft+k*width,vside,width,height])
			ax[k].set_ylim([-0.5 ,14])
			ax[k].set_xlim([0.5 ,len(plotexpe)+0.5])
			ax[k].set_xticks(np.arange(len(plotexpe))+1)
			ax[k].xaxis.set_ticklabels(['1','2','3','4','5','6'])
			ax[k].set_title(namestyles[name][0])
			if k==0:
				ax[k].set_ylabel('Mean Effective Pressure in Band',labelpad=0)
			if k==12:
				ax[k].yaxis.tick_right()
				ax[k].yaxis.set_ticks_position('both')
				ax[k].yaxis.set_label_position('right')
				ax[k].set_ylabel('Mean Effective Pressure in Band',labelpad=0)

	return fig,ax

# }}}
# {{{ MapLayout() 
def MapLayout():
	plotnum=9
	fig,ax=plt.subplots(1,plotnum,sharey='row',figsize=(plotnum+1,5))
	hleft=0.7/(plotnum+1.)
	width=0.95/(plotnum+1.)
	vside=0.12
	height=0.75
	k=-1

	for i,name in enumerate(names):
		if len(glob(resdir+name+'/A5_'+name+'.nc'))>0:
			k+=1
			ax[k].set_position([hleft+k*width,vside,width,height])
			ax[k].set_xlim([0,20])
			ax[k].set_ylim([0,100])
			ax[k].set_xticks([5,10,15])
			ax[k].set_yticks([10,30,50,70,90])
			ax[k].set_title(namestyles[name][0])
			ax[k].set_xlabel('$y$ (km)')
			if k==0:
				ax[k].set_ylabel('$x$ (km)')
			if k==11:
				ax[k].yaxis.tick_right()
				ax[k].yaxis.set_ticks_position('both')
				ax[k].yaxis.set_label_position('right')

	return fig,ax

# }}}
# {{{ SectionLayout() 
def SectionLayout():
	fig,ax=plt.subplots(3,3,sharey='row',figsize=(16,4))
	hleft=0.05
	hspace=0.01
	width=0.3
	vbot=0.1
	vtop=0.02
	vspace=0.04
	height=0.25
	k=-1

	for i,name in enumerate(names):
		if len(glob(resdir+name+'/A1_'+name+'.nc'))>0:
			k+=1
			ax[k/3,np.mod(k,3)].set_position([hleft+(np.mod(k,3))*(width+hspace),1-(vtop+height+(k/3)*(height+vspace)),width,height])
			ax[k/3,np.mod(k,3)].set_ylim([0,15])
			ax[k/3,np.mod(k,3)].set_xlim([0,100])
			ax[k/3,np.mod(k,3)].set_yticks([0,5,10])
			ax[k/3,np.mod(k,3)].set_yticklabels([])
			ax[k/3,np.mod(k,3)].set_xticks([10,30,50,70,90])
			ax[k/3,np.mod(k,3)].set_xticklabels([])
			ax[k/3,np.mod(k,3)].set_title(namestyles[name][0],
																		 va='top',ha='left',
																		 position=(0.05,0.8))
			if np.mod(k,3)==0:
				ax[k/3,np.mod(k,3)].set_yticklabels(['$0$','$5$','$10$'])
				if (k/3)==1:
					ax[k/3,np.mod(k,3)].set_ylabel('Effective Pressure (MPa)')
			if np.mod(k,3)==2:
				ax[k/3,np.mod(k,3)].yaxis.tick_right()
				ax[k/3,np.mod(k,3)].yaxis.set_ticks_position('both')
				ax[k/3,np.mod(k,3)].yaxis.set_label_position('right')
				ax[k/3,np.mod(k,3)].set_yticklabels(['$0$','$5$','$10$'])
			if k/3==2:
				ax[k/3,np.mod(k,3)].set_xticklabels(['$10$','$30$','$50$','$70$','$90$'])
				ax[k/3,np.mod(k,3)].set_xlabel('$x$ coordinate (km)',labelpad=0)

	return fig,ax

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
# {{{ ColorMapping():
def ColorMapping():
	cdict1 = {'alpha': [(0.0, 1.0, 1.0), (0.05, 1.0, 1.0), (0.1, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
					 'blue': [(0.0, 0.96, 0.96), (0.05, 1.0, 1.0), (0.1, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
					 'green': [(0.0, 0.74, 0.74), (0.05, 1.0, 1.0), (0.1, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
					 'red': [(0.0, 0.16, 0.16), (0.05, 1.0, 1.0), (0.1, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)]}
	blue_red1 = colors.LinearSegmentedColormap('BlueRed1', cdict1)
	return blue_red1
# }}}
# {{{ class FillandLineHandler(object):

class FillandLineHandler(object):
	def __init__(self, color):
		self.color=color
	def legend_artist(self, legend, orig_handle, fontsize, handlebox):
		x0, y0 = handlebox.xdescent, handlebox.ydescent
		width, height = handlebox.width, handlebox.height
		patch = plt.Rectangle([x0, y0], width, height, facecolor=self.color,alpha=0.3,
													edgecolor=self.color, transform=handlebox.get_transform())
		patch2 = plt.Line2D((x0,x0+width), (y0+0.5*height,y0+0.5*height),
												color=self.color, transform=handlebox.get_transform(),linewidth=2)
		handlebox.add_artist(patch)
		handlebox.add_artist(patch2)
		return patch

# }}}
if __name__=="__main__":
	PlotList={'MultiSec':MultiSec,
						'MwerSec':MwerSec,
						'MwerMeanSec':MwerMeanSec,
						'BandVal':BandVal,
						'Map':Mapping,
						'BandEvol':BandEvol,
						'BandEvolOneCol':BandEvolOneCol,
						'BandEvolChannel':BandEvolChannel,
						'Input':InputCurve,
						'Section':OneSection}

	TitleList={'N':['effective pressure',u'MPa']}
	if len(sys.argv)==1: # print help if no filename plot name is given:
		print __doc__
		sys.exit(0)

	resdir=os.getcwd()+'/'
	figdir='/home/bfl022/Documents/ecrits/FiguresPool/ExpeSHMIP/'

	DataName='S'

	names =['db','rh','cdf','id','jd','jsb','as','as_2','sb','bf','mh1','mh2','og','mw','mw_prime','og_2']#
	namestyles = {'db':['$db$',(0.,0.,0.,1.)],
								'cdf':['$cdf$',(0.3,0.3,0.3,1.)],
								'id':['$id$',(0.6,0.6,0.6,1.)],
								'bf':['$bf$',(0.5,0.5,0.,1.)],
								'sb':['$sb$',(0.7,0.7,0.,1.)],
								'rh':['$rh$',(0.7,0.7,0.,1.)],
								'jd':['$jd$',(0.,0.,0.5,1.)],
								'jsb':['$jsb$',(0.,0.,0.7,1.)],
								'as':['$as$',(0.,0.,1.,1.)],
								'as_2':['$as_2$',(0.,0.,1.,1.)],
								'mh1':['$mh2$',(0.,1.,1.,1.)],
								'mh2':['$mh2_2$',(0.,0.5,0.5,1.)],
								'mw':['$mw$',(0.,.7,0.7,0.6)],
								'mw_prime':['$mw_2$',(0.,.7,0.7,0.6)],
								'og':['$og$',(0.,1.,1.,0.3)],
								'og_2':['$og_2$',(0.,0.,1.,1.)]}

	if sys.argv[1] in ['MultiSec','MwerSec','MwerMeanSec','BandVal']:
		expe='A'
		plotexpe=['A1','A2','A3','A4','A5','A6']
	elif sys.argv[1]=='Map':
		expe='A5'
		plotexpe=['A5']
		names=names[4:]
	elif sys.argv[1] in['BandEvol','Input','Section','BandEvolOneCol','BandEvolChannel']:
		plotexpe=['D1','D2','D3','D4','D5']
		if sys.argv[1]=='Section':
			expe='A1'
			names=['db','id','jd','jsb','sb','bf','mh2','og','mw']#
		elif sys.argv[1]=='BandEvolChannel':
			expe='D4'
			names =['db','id','sb','bf','mh1','mh2','og','mw','mw_prime','og_2']
		else:
			expe='D'
			names =['db','id','jd','jsb','sb','bf','mw']#
	os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
	try:
		fig=PlotList[sys.argv[1]]()
		fig.savefig(figdir+'PlotPres'+sys.argv[1]+'_'+DataName+'_'+expe+'.png',format='png')
		fig.savefig(figdir+'PlotPres'+sys.argv[1]+'_'+DataName+'_'+expe+'.pdf',format='pdf')
	except KeyError:

		print 'Bad keyword for the figure production, run without keyword to get the help and check the spelling.'

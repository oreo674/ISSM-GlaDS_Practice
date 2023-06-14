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

# {{{ MultiSec()

def MultiSec(dummy):
	if plotexpe[0].startswith('A'):
		layout,az=SliceLayout(13,'A')
		midpoint=10.
		Exptag='A'
	elif plotexpe[0].startswith('E'):
		layout,az=SliceLayout(11,'E')
		midpoint=0.
		Exptag='E'
	k=-1

	for i,name in enumerate(names):
		if len(glob(resdir+name+'/'+Exptag+'*_'+name+'.nc'))>0:
			k+=1
			verts=[]
			limitedverts=[]
			for j,tag in enumerate(plotexpe):
				NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
				try:
					DatFile	 = Dataset(NCFile, mode='r')
					data,coords=GetVar(DatFile,DataName,name,False)
					thick,coords=GetVar(DatFile,'H',name,False)
					DatFile.close()
				except (RuntimeError,IOError):
					print ('File {} does not exist'.format(NCFile))
					continue

				griddeddata,gridpoint_x,gridpoint_y=Griding(data,coords)
				griddedthick,gridpoint_x,gridpoint_y=Griding(thick,coords)

				gridx=gridpoint_x[:,0]

				if Exptag=='E':
					databand=np.where(abs(gridpoint_y[0,:])<=bandwith*1.0e-3)[0]
					meandata=np.nanmean(griddeddata[:,databand],axis=1)
				else:
					meandata=np.nanmean(griddeddata,axis=1)

				# if name=='rh' and tag=='E4':
				# 	meandata=meandata[:13]
				# 	gridx=gridx[:13]

				verts.append(list(zip(gridx,j*np.ones(np.shape(meandata)),meandata)))
				if name=='rh':
					if tag=='E2':
						limitedverts.append(list(zip(gridx[:63],j*np.ones(np.shape(meandata[:64])),meandata[:64])))
					else:
						limitedverts.append(list(zip(gridx,j*np.ones(np.shape(meandata)),meandata)))
				elif name=='sb':
					if tag in['E1','E2','E3','E5']:
						limitedverts.append(list(zip(gridx,j*np.ones(np.shape(meandata)),meandata)))
					elif tag=='E4':
						limitedverts.append(list(zip(gridx[:17],j*np.ones(np.shape(meandata[:17])),meandata[:17])))
						limitedverts.append(list(zip(gridx[76:],j*np.ones(np.shape(meandata[76:])),meandata[76:])))
				elif name=='jsb':
					if tag in ['E1','E2','E3','E4']:
						if tag=='E1':
							toplim=84
						elif tag=='E2':
							toplim=90
						elif tag=='E3':
							toplim=93
						elif tag=='E4':
							toplim=97
						limitedverts.append(list(zip(gridx[:toplim],j*np.ones(np.shape(meandata[:toplim])),meandata[:toplim])))
					elif tag=='E5':
						limitedverts.append(list(zip(gridx,j*np.ones(np.shape(meandata)),meandata)))
				else:
					limitedverts=verts
				x=np.tile(gridx,(2,1))
				y=np.tile(j*np.ones(np.shape(gridx)),(2,1))
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

			if Exptag=='A':
				if name=='sb':
					az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','w','k','w','k'], linewidth=1, linestyles='-'))
				elif name=='bf':
					az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','w','k','w'], linewidth=1, linestyles='-'))
				elif name=='cdf':
#					az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','k','w','k'], linewidth=1, linestyles='-'))
					az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','w','k'], linewidth=1, linestyles='-')) #waiting for A2
				elif name=='as':
					az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','w','k','k','k'], linewidth=1, linestyles='-'))
				elif name=='jd':
					az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','k','k','w','k'], linewidth=1, linestyles='-'))
				elif name in ['mw']:
					az[k].add_collection3d(Line3DCollection(verts, colors='k', linewidth=1, linestyles='-'))
					az[k].add_collection3d(Line3DCollection(verts, colors=['k','k','w','k','w','k'], linewidth=1, linestyles=['-','-','--','-','--','-']))
				else:
					az[k].add_collection3d(Line3DCollection(verts, colors='k', linewidth=1, linestyles='-'))
			elif Exptag=='E':
				az[k].add_collection3d(Line3DCollection(limitedverts, colors='k', linewidth=1, linestyles='-'))
				az[k].add_collection3d(Line3DCollection(verts, colors='k', linewidth=1, linestyles=(0,(1,1))))

			if Exptag=='A':
				if name in ['sb','bf','cdf','as','jd']:
					az[k].set_title(namestyles[name][0],fontsize=9,rotation=10,
													va='bottom',ha='left',
													position=(0.16,0.7),color='r')
					az[k].text(90,5,13,string.lowercase[k],fontsize=9,va='center',ha='center')
				else:
					az[k].set_title(namestyles[name][0],fontsize=9,rotation=10,
													va='bottom',ha='left',
													position=(0.16,0.7))
					az[k].text(90,5,13,string.lowercase[k],fontsize=9,va='center',ha='center')

			elif Exptag=='E':
				az[k].set_title(namestyles[name][0],fontsize=9,rotation=10,
												va='bottom',ha='left',
												position=(0.16,0.75))
				az[k].text2D(0.65,0.85,string.lowercase[k],fontsize=9,va='center',ha='right',transform=az[k].transAxes)

			if k==1 and tag.startswith('E'):
				geoverts=[]
				for p,para in enumerate([0.05,0.,-0.1,-0.5,-0.7]):
					geomx=np.arange(0,6.1,0.1)
					surf   = 100.*(geomx*1.0e3+200.)**0.25+ 100.*geomx*1.0e3/6.0e3 - 2.0e10**0.25 + 1.0
					s_xend = surf[-1]
					f_func  = para*(geomx*1.0e3) + (geomx*1.0e3)**2. * (s_xend-para*6.0e3)/6.0e3**2.
					f_Bench = 0.05*geomx*1.0e3 + (geomx*1.0e3)**2. * (s_xend-0.05*6.0e3)/6.0e3**2.
					h_func  = 0.5e-6 * (5 - 4.5*(geomx*1.0e3)/6.0e3) * (surf-f_func)/(surf-f_Bench)
					bed= f_func + h_func
					gridx=np.append(geomx,np.flipud(geomx[:-1]))
					geom=np.append(bed,np.flipud(surf[:-1]))
					geoverts.append(list(zip(gridx,p*np.ones(np.shape(gridx)),geom)))
				az[10].add_collection3d(Poly3DCollection(geoverts, facecolor=(0.8,0.8,0.8),edgecolor='k', linewidth=1, linestyles='-'))
				az[10].add_collection3d(Line3DCollection(geoverts, colors='k', linewidth=1, linestyles=':'))
				az[10].text2D(0.25,0.8,string.lowercase[10],fontsize=9,va='center',ha='right',transform=az[10].transAxes)
				az[10].invert_yaxis()
				az[10].view_init(elev=25., azim=245)
				az[10].set_title('geometry',fontsize=9,rotation=-53,
												va='bottom',ha='left',
												position=(0.17,0.04))
				az[10].text(90,5,13,string.lowercase[11],fontsize=9,va='center',ha='center')


			az[k].invert_yaxis()
			az[k].view_init(elev=25., azim=245)

	if Exptag=='A':
		barax = layout.add_axes([0.55,0.3,0.4,0.05])
	elif Exptag=='E':
		barax = layout.add_axes([0.1,0.15,0.55,0.05])
	bartitle=barax.set_title('Percentage of Flux in the efficient system',fontsize=9)
	cbar = mpl.colorbar.ColorbarBase(barax, norm=norm,cmap=cmap, orientation='horizontal')
	cbar.solids.set_edgecolor("face")
	plt.show
	return layout

# }}}

# {{{ SectionsDiff()

def SectionsDiff(dummy):

	width=8.6/2.54
	layout,ax=plt.subplots(9,2,sharex=True,figsize=(width,width*2))
	hleft=0.1
	hright=0.15
	hspace=0.02
	width=(1-hleft-hright-hspace)/2.
	vside=0.045
	vspace=0.01
	height=0.85/9

	k=-1

	#gathering mouins input (once and for all)
	inputval=np.zeros((101,len(plotexpe)))
	for t,tag in enumerate(plotexpe):
		moulinfile='/home/bfl022/Model/GHIP/hydro_intercomparison/input_functions/source/'+tag+'_M.csv'
		with open(moulinfile) as csvfile:
			MoulinPos=csv.reader(csvfile,delimiter=',')
			for row in MoulinPos:
				Xs =int(row[1])
				inputval[Xs/1000,t]+=float(row[3])

	for i,name in enumerate(names):
		#checking existence of the two needed runs
		if len(glob(resdir+name+'/B*_'+name+'.nc'))>0 and len(glob(resdir+name+'/A5_'+name+'.nc'))>0:
			NCDiff=resdir+name+'/A5_'+name+'.nc'
			DiffFile = Dataset(NCDiff, mode='r')
			diffdata,diffcoords=GetVar(DiffFile,DataName,name,False)
			DiffFile.close()
			k+=1
			#Now work on it
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
				griddeddiff,gridpoint_x,gridpoint_y=Griding(diffdata,diffcoords)
				floatation=1.0e-6*910*9.81*(1.0+6.*pow(gridpoint_x*1.0e3+5.0e3,0.5)-6.0*pow(5.0e3,0.5))#this is actually ice pressure

				#		normdiff=100.*((griddeddata-griddeddiff)/floatation)
				normdiff=(griddeddata-griddeddiff)

				meannormdiff=np.nanmean(normdiff,axis=1)
				maxnormdiff=np.nanmax(normdiff,axis=1)
				minnormdiff=np.nanmin(normdiff,axis=1)

				ax[k,j].set_position([hleft+j*(width+hspace),vside+(8-k)*(height+vspace),width,height])
				#axis stuff
				ax[k,j].set_xlim([0 ,100])
				ax[k,j].set_xticks([10,30,50,70,90])
				ax[k,j].xaxis.set_ticklabels([])

				if k==0:
					ax[k,j].set_title(tag,fontsize=9,position=(0.5,1.))
				if k==8:
					ax[k,j].xaxis.set_ticklabels(['$10$','$30$','$50$','$70$','$90$'],fontsize=9)
					ax[k,j].tick_params(axis='x', which='major', pad=2)
					ax[k,j].set_xlabel('$x$ (km)',labelpad=0)

				ax[k,j].yaxis.tick_right()
				ax[k,j].yaxis.set_ticks_position('both')
				ax[k,j].yaxis.set_label_position('right')

				ax2=ax[k,j].twinx()
				ax2.set_position(ax[k,j].get_position())
				ax2.bar(gridpoint_x[:,0],inputval[:,j],align='center',label='input',color='r',log=True,edgecolor='none')

				ax2.set_xlim([0 ,100])
				ax2.set_ylim([0.5 ,120])
				ax2.yaxis.set_ticklabels([])
				ax2.yaxis.set_ticks_position('none')
				if j==0:
					ax2.set_title(namestyles[name][0],fontsize=9,
												va='top',ha='right',
												position=(0.95,0.2))
				else:
					ax2.set_title(namestyles[name][0],fontsize=9,
												va='bottom',ha='right',
												position=(0.95,0.65))
					# ax2.set_title(namestyles[name][0],fontsize=9,
					# 							va='bottom',ha='right',
					# 							position=(0.95,0.6))

				ax[k,j].fill_between(gridpoint_x[:,0], maxnormdiff,minnormdiff, facecolor='b',edgecolor='b',alpha=0.3)
				ax[k,j].plot(gridpoint_x[:,0], meannormdiff, linestyle='-',linewidth=1.5,label=namestyles[name][0],color='b')

				ax[k,j].set_ylim([np.nanmin(minnormdiff[5:]) ,np.nanmax(maxnormdiff[5:])])
				ax[k,j].set_zorder(ax[k,0].get_zorder()+1)
				ax[k,j].patch.set_visible(False)
				ax[k,j].yaxis.set_ticks_position('both')

				# if j==0:
				# 	if name in ['id','jd','jsb','as']:
				# 		ax[k,j].set_ylim([-7,17])
				# 		ax[k,j].set_yticks([-5,0,5,10,15])
				# 		ax[k,j].yaxis.set_ticklabels(['$-5$','$0$','$5$','$10$','$15$'],fontsize=9)
				# 	else:
				# 		ax[k,j].set_ylim([-7 ,30])
				# 		ax[k,j].set_yticks([-7,0,7,14,21,28])
				# 		ax[k,j].yaxis.set_ticklabels(['$-7$','$0$','$7$','$14$','$21$','$28$'],fontsize=9)
				# else:
				# 	if name=='id':
				# 		ax[k,j].set_ylim([-1.5,1.5])
				# 		ax[k,j].set_yticks([-1.4,-0.7,0,0.7,1.4])
				# 		ax[k,j].yaxis.set_ticklabels(['$-1.4$','$-0.7$','$0$','$0.7$','$1.4$'],fontsize=9)
				# 	else:
				# 		ax[k,j].set_ylim([-12 ,12])
				# 		ax[k,j].set_yticks([-10,-5,0,5,10])
				# 		ax[k,j].yaxis.set_ticklabels(['$-10$','$-5$','$0$','$5$','$10$'],fontsize=9)
				# 	ax[k,j].tick_params(axis='y', which='major', pad=2)
				# 	ax[k,j].yaxis.tick_right()
				# 	ax[k,j].yaxis.set_ticks_position('both')
				# 	ax[k,j].yaxis.set_label_position('right')

				if j==0:
					ax[k,j].set_ylim([-1,4])
					ax[k,j].set_yticks([-1,0,1,2,3])
					ax[k,j].yaxis.set_ticklabels(['$-1$','$0$','$1$','$2$','$3$'],fontsize=9)
				else:
					# ax[k,j].set_ylim([-1,4])
					# ax[k,j].set_yticks([-1,0,1,2,3])
					# ax[k,j].yaxis.set_ticklabels(['$-1$','$0$','$1$','$2$','$3$'],fontsize=9)
					ax[k,j].set_ylim([-0.5,0.6])
					ax[k,j].set_yticks([-0.4,-0.2,0,0.2,0.4])
					ax[k,j].yaxis.set_ticklabels(['$-0.4$','$-0.2$','$0$','$0.2$','$0.4$'],fontsize=9)
					ax[k,j].tick_params(axis='y', which='major', pad=2)
					ax[k,j].yaxis.tick_right()
					ax[k,j].yaxis.set_ticks_position('both')
					ax[k,j].yaxis.set_label_position('right')

				ax[k,j].text(0.05,0.7,string.lowercase[(k)+j*9],fontsize=9,va='bottom',ha='left',transform=ax[k,j].transAxes)
	ax[0,0].text(0.1*hleft,vside+height*((k+1)/2.)+((k)/2.)*vspace,'Effective pressure difference (MPa)',fontsize=9,va='center',ha='left',rotation=90,transform=layout.transFigure)
	ax[0,1].text(1-(0.1*hright),vside+height*((k+1)/2.)+((k)/2.)*vspace,'Effective pressure difference (MPa)',fontsize=9,va='center',ha='right',rotation=90,transform=layout.transFigure)
	plt.show
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
		baselist=[2]
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
		baselist=[2]
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
		baselist=[3]
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
		if len(glob(resdir+name+'/'+Exptag+'*_'+name+'.nc'))>0:
			if name=='mh1' and Exptag=='D':
				framelist=[1]
			# elif name=='sb' and Exptag=='D':
			# 	framelist=[4]
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
			ax[plotrow-1,plotcol-1].tick_params(axis='x', which='major', pad=2)
			ax[plotrow-1,plotcol-1].set_xlabel('run',labelpad=0)

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
				elif name in ['jd','og','as','mw','mw_A']:
					ax2.set_ylim([0 ,4])
				elif name in ['mw_prime','mh2','mw_A_ev','sb']:
					ax2.set_ylim([0 ,0.2])
			elif Exptag=='D':
				ax2.set_ylim([0.5 ,8])
			elif Exptag=='F':
				ax2.set_ylim([-0.1 ,5])

			#==== Core computation==============
			for j,tag in enumerate(plotexpe):
				NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
				print NCFile
				try:
					DatFile	 = Dataset(NCFile, mode='r')
					time     = DatFile.variables['time'][:]
					data,coords=GetVar(DatFile,DataName,name,True)
					thick,coords=GetVar(DatFile,'H',name,False)
					base,coords=GetVar(DatFile,'B',name,False)
					DatFile.close()
					if name=='db':
						RefFile						= Dataset(resdir+'mw/'+tag+'_mw.nc', mode='r')
						reftime						= RefFile.variables['time'][:]
						refdata,refcoords	=	GetVar(RefFile,DataName,'mw',True)
						RefFile.close()
				except (RuntimeError,IOError):
					print ('File {} does not exist'.format(NCFile))
					if j in framelist:
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
				Meanwidth=np.NaN*np.ones((3))
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

				for frame in framelist:
					if j==frame:
						#========patch on right to show what we have on left======
						if Exptag=='C':
							rect = [Rectangle((frame+0.5,-12), 1.0, 24)]
						elif Exptag=='D':
							rect = [Rectangle((frame+0.5,-1.8), 1.0, 3.24)]
						elif Exptag=='F':
							rect = [Rectangle((frame+0.5,-3), 1.0, 5)]

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
								elif name in ['jd','og','as','mw','mw_A']:
									ax[k+1,0].set_ylim([-2 ,3.5])
							elif Exptag=='D':
								ax[k+1,count].set_xlim([0.25 ,1])
								ax[k+1,count].set_xticks([0.33,0.49,0.66,0.83,1.0])
								ax[k+1,0].set_ylim([-5,10])
								ax[k+1,1].set_ylim([0,1000])
								ax[0,0].set_ylim([0,1000])
								ax[0,1].set_ylim([0,1000])
							elif Exptag=='F':
								ax[k+1,count].set_xlim([0.25 ,1])
								ax[k+1,count].set_xticks([0.33,0.49,0.66,0.83,1.0])
								if name in ['db','jsb']:
									ax[k+1,0].set_ylim([0,2])
								else:
									ax[k+1,0].set_ylim([0,6])
								ax[k+1,1].set_ylim([0 ,8])
								ax[0,0].set_ylim([0 ,8])
								ax[0,1].set_ylim([0 ,8])

							if k==2:
								if Exptag=='C':
									ax[0,count].text(0.05,0.85,string.lowercase[0],fontsize=9,va='center',ha='center',transform=ax[0,count].transAxes)
								else:
									ax[0,count].text(0.05,0.85,string.lowercase[count],fontsize=9,va='center',ha='center',transform=ax[0,count].transAxes)
								ax[0,count].set_position([hleft+count*(plotwidth+hspace),vbot+(plotrow-1)*(plotheight+vspace),plotwidth,plotheight])
								ax[0,count].yaxis.set_major_locator(MaxNLocator(5,prune='both'))
								if Exptag=='C':
									#ax[0,count].set_title(tag,position=(0.5,0.95),fontsize=9)
									# ax[0,count].spines['left'].set_color('r')
									# ax[0,count].tick_params(axis='y', colors='r')
									ax[0,count].set_xlim([1 ,24])
									ax[0,count].set_xticks([3,6,9,12,15,18,21,24])
									ax[plotrow-1,count].xaxis.set_ticklabels(['$3$','$6$','$9$','$12$','$15$','$18$','$21$','$24$'])
									ax[plotrow-1,count].set_xlabel('Time (h)',labelpad=0)
								elif Exptag=='D':
									ax[0,count].set_xlim([0.25 ,1])
									ax[0,count].set_xticks([0.33,0.49,0.66,0.83,1.0])
									ax[plotrow-1,count].xaxis.set_ticklabels(['$4$','$6$','$8$','$10$','$12$'])
									ax[plotrow-1,count].set_xlabel('Time (month)',labelpad=0)
								elif Exptag=='F':
									ax[0,count].set_xlim([0.25 ,1])
									ax[0,count].set_xticks([0.33,0.49,0.66,0.83,1.0])
									ax[plotrow-1,count].xaxis.set_ticklabels(['$4$','$6$','$8$','$10$','$12$'])
									ax[plotrow-1,count].set_xlabel('Time (month)',labelpad=0)
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
						pressure=[]
						if name=='db':
							dbN=ax[k+1,0].plot(time, Mean[:], linestyle='-',linewidth=1.5,color='g',label='Mean effective pressure')[0]
							# IDS.append(ax[k+1,1].plot(time, np.nanmean(sheetflux,axis=1), linestyle='--',linewidth=1.5,color='g',label='Mean inefficient syst. flux')[0])
							# EDS.append(ax[k+1,1].plot(time, np.nanmean(channelflux,axis=1), linestyle=':',linewidth=1.5,color='g',label='Mean efficient syst. flux')[0])
							IDS.append(ax[k+1,1].plot(time, np.nanmean(sheetflux,axis=1),dashes=[6, 2],linewidth=1.5,color='g',label='Mean inefficient syst. flux')[0])
							EDS.append(ax[k+1,1].plot(time, np.nanmean(channelflux,axis=1),dashes=[2, 2],linewidth=1.5,color='g',label='Mean efficient syst. flux')[0])
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
		if tag.startswith('E'):
			if name not in ['db','rh','cdf','id']:
				Width[np.where(Width>2*bandwith)]=2*bandwith
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

	if name in['mh1','mh2','mw','mw_prime','og','og_prime','cdf','db','sb','jd','jsb','as','as_2','id','rh','bf','mw_A','mw_A_ev']:
		NCFile=resdir+name+'/'+tag+'_'+name+'.nc'
		try:
			DatFile	 = Dataset(NCFile, mode='r')
			if name in['mh2','mw','mw_prime','og','og_prime','mh1','mw_A','mw_A_ev']:
				IneffFlux,Ineffcoords=GetVar(DatFile,'q',name,TimeVar)
				EffFlux,Effcoords=GetVar(DatFile,'Q',name,TimeVar)
				if name=='mh1':
					Xs=DatFile.variables['coords3'][0,:]
					Ys=DatFile.variables['coords3'][1,:]
				else:
					Xs=DatFile.variables['coords1'][0,:]
					Ys=DatFile.variables['coords1'][1,:]
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
				Conductivity,Ineffcoords=GetVar(DatFile,'T',name,TimeVar)
				ischannel=np.ones(np.shape(Conductivity))
#				ischannel[np.where(Conductivity==np.nanmin(Conductivity))]=0
				ischannel[np.where(Conductivity<0.1)]=0
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
			elif name=='id':
				Channelisation,Ineffcoords=GetVar(DatFile,'R',name,TimeVar)
				EffFlux,Effcoords=GetVar(DatFile,'Q',name,TimeVar)
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
			if name in['mh2','mw','mw_prime','og','og_prime','mh1','mw_A','mw_A_ev']:
				if tag.startswith('E'):
					sheetloc=np.where(np.logical_and(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																													Ineffcoords[0,:]>=(loc-0.5*Step)),
																					 abs(Ineffcoords[1,:])<bandwith))[0]
					if name=='mh1':
						channelloc=np.where(np.logical_and((Xs[intconnect[0,:]]-loc+1)*(Xs[intconnect[1,:]]-loc+1)<0,
																							 np.logical_and(abs(Ys[intconnect[0,:]])<bandwith,abs(Ys[intconnect[1,:]])<bandwith)))[0]
					else:
						channelloc=np.where(np.logical_and((Xs[intconnect[0,:]]-loc)*(Xs[intconnect[1,:]]-loc)<0,
																							 np.logical_and(abs(Ys[intconnect[0,:]])<bandwith,abs(Ys[intconnect[1,:]])<bandwith)))[0]
					if loc in [min(Lines),max(Lines)] and len(channelloc)==0:
						channelloc=np.where(np.logical_and((Xs[intconnect[0,:]]-loc)*(Xs[intconnect[1,:]]-loc)==0,
																							 np.logical_and(abs(Ys[intconnect[0,:]])<bandwith,abs(Ys[intconnect[1,:]])<bandwith)))[0]
				else:
					sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
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
				if tag.startswith('E'):
					sheetloc=np.where(np.logical_and(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																													Ineffcoords[0,:]>=(loc-0.5*Step)),
																					 abs(Ineffcoords[1,:])<bandwith))[0]
					channelloc=np.where(np.logical_and(np.logical_and(Effcoords[0,:]<(loc+0.5*Step),
																														Effcoords[0,:]>=(loc-0.5*Step)),
																						 abs(Effcoords[1,:])<bandwith))[0]
				else:
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
				if tag.startswith('E'):
					channelloc=np.where(np.logical_and(np.logical_and(Effcoords[0,:]<(loc+0.5*Step),
																														Effcoords[0,:]>=(loc-0.5*Step)),
																						 abs(Effcoords[1,:])<bandwith))[0]
				else:
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
				if tag.startswith('E'):
					sheetloc=np.where(np.logical_and(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																													Ineffcoords[0,:]>=(loc-0.5*Step)),
																					 abs(Ineffcoords[1,:])<bandwith))[0]
				else:
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
				if tag.startswith('E') or tag.startswith('F'):
					sheetloc=np.where(np.logical_and(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																													Ineffcoords[0,:]>=(loc-0.5*Step)),
																					 abs(Ineffcoords[1,:])<bandwith))[0]
				else:
					sheetloc=np.where(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																					 Ineffcoords[0,:]>=(loc-0.5*Step)))[0]
				if TimeVar:
					SheetFlux[:,i]=np.nanmean(IneffFlux[:,sheetloc],axis=1)*Width[i]
				else:
					SheetFlux[i]=np.nanmean(IneffFlux[sheetloc])*Width[i]
			#===jsb=db===
			elif name in ['jsb','db']:
				if tag.startswith('E'):
					sheetloc=np.where(np.logical_and(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																													Ineffcoords[0,:]>=(loc-0.5*Step)),
																					 abs(Ineffcoords[1,:])<bandwith))[0]
				else:
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
				if tag.startswith('E'):
					sheetloc=np.where(np.logical_and(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																													Ineffcoords[0,:]>=(loc-0.5*Step)),
																					 abs(Ineffcoords[1,:])<bandwith))[0]
				else:
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
				if tag.startswith('E'):
					sheetloc=np.where(np.logical_and(np.logical_and(Ineffcoords[0,:]<(loc+0.5*Step),
																													Ineffcoords[0,:]>=(loc-0.5*Step)),
																					 abs(Ineffcoords[1,:])<bandwith))[0]
					channelloc=np.where(np.logical_and(np.logical_and(Effcoords[0,:]<(loc+0.5*Step),
																														Effcoords[0,:]>=(loc-0.5*Step)),
																						 abs(Effcoords[1,:])<bandwith))[0]
				else:
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

# {{{ SliceLayout(frames,Exp)

def SliceLayout(frames,Exp):

	width=17.8/2.54
	if Exp=='A':
		fig = plt.figure(figsize=(width,width))
	elif Exp=='E':
		fig = plt.figure(figsize=(width,np.ceil(frames/5.)/5.*width))
	ax=[]
	fig.subplots_adjust(left=0.03, right=0.99, top=1.02, bottom=0.02, wspace=0.05, hspace=-0.10)
	for i in np.arange(0,frames):
		if  Exp=='A':
			if i>8:
				location=i+4
			elif i>6:
				location=i+2
			else:
				location=i+1

			ax.append(fig.add_subplot(4, 4, location, projection='3d'))
		elif i==10 and Exp=='E':
			ax.append(fig.add_subplot(int(np.ceil(frames/5.)), 5, 15, projection='3d'))
		else:
			ax.append(fig.add_subplot(int(np.ceil(frames/5.)), 5, i+1, projection='3d'))
		ax[i].patch.set_visible(False)
		if i in [3,6,9,10,11,12] and Exp=='A':
			ax[i].set_xlabel('$x$ (km)',labelpad=-11,fontsize=9)
		if i in [5,6,7,8,10] and Exp=='E':
			ax[i].set_xlabel('$x$ (km)',labelpad=-11,fontsize=9)
		ax[i].tick_params(axis='x', which='major', pad=0)
		if Exp=='A':
			ax[i].set_xlim3d(0, 100)
			ax[i].set_xticks([10,50,90])
			ax[i].set_xticklabels(['$10$','$50$','$90$'],ha='center',va='bottom',fontsize=9,rotation=10)
		elif Exp=='E':
			ax[i].set_xlim3d(0, 6)
			ax[i].set_xticks([1,2,3,4,5])
			ax[i].set_xticklabels(['$1$','$2$','$3$','$4$','$5$'],ha='center',va='bottom',fontsize=9,rotation=10)

		ax[i].set_ylim3d(0,len(plotexpe))
		if Exp=='A':
			ax[i].set_ylim3d(-0.1,5.1)
			ax[i].set_yticks([0,1,2,3,4,5])
			ax[i].set_yticklabels(['$1$','$2$','$3$','$4$','$5$','$6$'],fontsize=9,ha='right',va='center')#
			#if i in [9,10,11,12]:
			ax[i].set_ylabel('A',labelpad=-10,fontsize=9)
		elif Exp=='E':
			ax[i].set_ylim3d(-0.1, 4.1)
			ax[i].set_yticks([0,1,2,3,4])
			ax[i].set_yticklabels(['$1$','$2$','$3$','$4$','$5$'],fontsize=9,ha='right',va='center')#
			if i in [6,7,8,10,11]:
				ax[i].set_ylabel('E',labelpad=-10,fontsize=9)
		ax[i].tick_params(axis='y', which='major', pad=-5)

		ax[i].zaxis.set_rotate_label(False)
		if i in [0,5,10] and Exp=='E':
			if i==10:
				ax[i].set_zlabel('Elevation (m)',fontsize=9,rotation=90)
			else:
				ax[i].set_zlabel('N (MPa)',labelpad=-5,fontsize=9,rotation=90)
		if i in [0,4,7,9] and Exp=='A':
			ax[i].set_zlabel('N (MPa)',labelpad=-5,fontsize=9,rotation=90)
		if Exp=='A':
			ax[i].set_zlim3d(0, 14)
			ax[i].set_zticks([0,2,4,6,8,10,12,14])
			ax[i].set_zticklabels(['$0$','$2$','$4$','$6$','$8$','$10$','$12$','$14$'],ha='right',va='center',fontsize=9)
		elif Exp=='E':
			if i==10:
				ax[i].set_zlim3d(-1000,600)
				ax[i].set_zticks([-800,-400,0,400])
				ax[i].set_zticklabels(['$-800$','$-400$','$0$','$400$'],ha='right',va='center',fontsize=9)
			else:
				ax[i].set_zlim3d(0,4)
				ax[i].set_zticks([0,1,2,3,4])
				ax[i].set_zticklabels(['$0$','$1$','$2$','$3$','$4$'],ha='right',va='center',fontsize=9)
				#ax[i].set_zlim3d(-3,6)
		ax[i].tick_params(axis='z', which='major', pad=-3)

		ax[i].set_frame_on(False)

	plt.show
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

# {{{ColorMapping():reduce white band and move allong

def ColorMapping():
	cdict1 = {'alpha': [(0.0, 1.0, 1.0), (0.10, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
					 'blue': [(0.0, 0.96, 0.96), (0.10, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
					 'green': [(0.0, 0.74, 0.74), (0.10, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
					 'red': [(0.0, 0.16, 0.16), (0.10, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)]}#, (0.12, 1.0, 1.0)
	blue_red1 = colors.LinearSegmentedColormap('BlueRed1', cdict1)
	return blue_red1

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


if __name__=="__main__":
	PlotList={'MultiSec':MultiSec,
						'SecDiff':SectionsDiff,
						'BandEvol':BandEvol}

	TitleList={'N':['effective pressure',u'MPa']}
	if len(sys.argv)==1: # print help if no filename plot name is given:
		print __doc__
		sys.exit(0)

	resdir=os.getcwd()+'/'
	figdir=resdir+'/figures'

	DataName='N'
	expe=raw_input('Give the tag of the experiment you want to plot: ')

	names =['db','id','rh','cdf','jd','jsb','as','sb','bf','mh1','mh2','og','og_prime','mw','mw_prime']#,,'as_2'
#	names =['db','rh','cdf','id','jd','jsb','as','sb','bf','mh1','mh2','mw_prime','mw_A_ev','og_prime','mw','mw_A']#,,'as_2'
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
								'mw_A_ev':['\\boldmath$mwA\'$',(0.,.7,0.7,0.6)]}

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
	os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))
	bandwith=100.
	try:
	#This is for suplementaries (need to put sub in select frame above)
	 # for subsuite in subs:
	 # 		fig=PlotList[sys.argv[1]](subsuite)
	 # 		try:
	 # 			fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+str(subsuite)+'.png',format='png')
	 # 			fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+str(subsuite)+'.pdf',format='pdf')
	 # 		except IOError:
	 # 			print('No Figure directory existing where needed, creating it here "{}"'.format(os.path.abspath(figdir)))
	 # 			os.makedirs(figdir+'/all')
	 # 			fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+str(subsuite)+'.png',format='png')
	 # 			fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+str(subsuite)+'.pdf',format='pdf')
	#This is for publication figure with the focus suite defined above
		fig=PlotList[sys.argv[1]](2)
		try:
			fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+'.png',format='png')
			fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+'.pdf',format='pdf')
		except IOError:
			print('No Figure directory existing where needed, creating it here "{}"'.format(os.path.abspath(figdir)))
			os.makedirs(figdir+'/all')
			fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+'.png',format='png')
			fig.savefig(figdir+'/all/'+sys.argv[1]+'_'+DataName+'_'+expe+'.pdf',format='pdf')

	except KeyError:

		print 'Bad keyword for the figure production, run without keyword to get the help and check the spelling.'

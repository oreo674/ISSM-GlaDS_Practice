#!/usr/bin/env python2

"""
Reads files.
"""
from netCDF4 import Dataset
import os
import sys
import string
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from SHMIPPlotter import *

# times in seconds
hour = 3600
day = 24*hour
year = 365*day

# scenarios, their number of runs, the indices which should be kept, and the unit of time
scens = {"A" :( 6, [0], 1),
         "B": (5, [0], 1),
         "C": (4, [-1,-12], hour),
         "D": (5, range(-200,-50), day),
         "E": (5, [0], 1),
         "F": (5, [-180], day)}

def participants(dir_or_out):
        """
        Retrun abbrev of all participants.
        """
        try:
                names = os.listdir(dir_or_out)
                return [name for name in names if len(name)==4]
        except:
                return dir_or_out.keys()

def readone(fl, var, coords, inds=[0]):
        """
        Read one variable and its coordiantes (well they need to be specified, correctly by hand)
        """
        ds = Dataset(fl, mode='r')
        data, (x,y) = GetVar(ds, var, coords, True, False)
        #if len(inds)>3: raise Exception()
        if len(data.shape)==1:
                return data, x, y, 0
        elif len(data.shape)==2:
                t = ds.variables['time'][inds]
                return data[inds,:], x, y, t
        else:
                error("ha")


def readall(dir):
        """
        Read all netcdf files for coordinates, effective pressure.  And maybe more in the future.
        """
        out = {}

        for name in participants(dir):
                # if name=="id": continue
                outout = {}
                out[name] = outout
                for fl in os.listdir(os.path.join(dir,name)):
                        if os.path.splitext(fl)[1]!=".nc" or fl[0]=='.': continue
                        run = string.replace(os.path.splitext(fl)[0], "_"+name, "")
                        if not 0<int(run[1])<scens[run[0]][1]: continue
                        ds = Dataset(os.path.join(dir,name,fl), mode='r')
                        N,x,y,t = readone(os.path.join(dir,name,fl), "N", "coords1", scens[run[0]][1])
                        H,x,y,t = readone(os.path.join(dir,name,fl), "H", "coords1")
                        B,x,y,t = readone(os.path.join(dir,name,fl), "B", "coords1")
                        outout[run] = {"N":N, "x":x, "y":y, 't':t, "H":H, "B":B}


        return out

def select_in_band(x, xband):
        """
        Returns indices of all points in x-band.
        """
        inds = []
        for i,xx in enumerate(x):
                if xband[0]<=xx<=xband[1]:
                        inds.append(i)
        return inds

def band_all(out, xbands_sqrt=([1e4,1.5e4],[5e4,5.5e4],[8.5e4,9e4]),
             xbands_bench=([1e3,1.5e3],[3e3,3.5e3],[5e3,5.5e3]) ):
        """
        Takes min,mean,max of N in the three bands.
        """
        for name in out:
                for run in out[name]:
                        outout = out[name][run]
                        if run[0]=="E" or run[0]=="F":
                                inds = map(lambda xband: select_in_band(outout["x"],xband), xbands_bench)
                        else:
                                inds = map(lambda xband: select_in_band(outout["x"],xband), xbands_sqrt)
                        outout["band-inds"] = inds
                        outout["band-means"] = map(lambda ind: np.mean(outout["N"][:,ind],1), inds)
                        outout["band-max"] = map(lambda ind: np.max(outout["N"][:,ind],1), inds)
                        outout["band-min"] = map(lambda ind: np.min(outout["N"][:,ind],1), inds)
        return


def plotall_xy_point_cloud(out):
        """
        Plot the mesh points.
        """
        for name in out:
                for run in out[name]:
                        outout = out[name][run]
                        x,y = outout["x"],outout["y"]
                        plt.plot(x,y,".")
                        plt.title(name+" "+run)
                        plt.show()
                        raw_input()
                        break

        return

def plot_band_Ns(out):
        """
        Plots the banded N for all participants.  One plot per scenario-run.
        """
        labs = ["bottom", "middle", "top"]
        npart = len(participants(out))
        for scen in scens:
                nruns,inds,timescale = scens[scen]
                ninds = len(inds)

                for ii in range(ninds): # loop over all plotting times.
                        fig, axs = plt.subplots(1,npart, sharey=True, sharex=True, figsize=(12, 4))
                        fig.subplots_adjust(wspace=0.01, left=0.05, top=0.93)
                        axs[0].set_ylabel("N (MPa)")

                        time = out['bf']["%s%i" % (scen,1)]['t']
                        if inds[0]!=0 or ninds>1:
                                axs[npart//2].set_xlabel("Scenario %s# at time %i" % (scen, time[inds[ii]]/timescale) )
                        else:
                                axs[npart//2].set_xlabel("Scenario %s#" % scen)
                        for i,name in enumerate(participants(out)): #enumerate(['mw']): #
                                x = np.zeros(nruns)
                                ys = [np.nan*np.zeros(nruns), np.nan*np.zeros(nruns), np.nan*np.zeros(nruns)]
                                yerrs = [np.nan*np.zeros([2,nruns]), np.nan*np.zeros([2,nruns]), np.nan*np.zeros([2,nruns])]
                                for run in xrange(0,nruns):
                                        key = "%s%i" % (scen,run+1)
                                        if not out[name].has_key(key): continue
                                        x[run] = run+1
                                        for jj in range(0,len(labs)):
                                                ys[jj][run] = out[name][key]["band-means"][jj][ii]
                                                yerrs[jj][:,run] = np.array([out[name][key]["band-min"][jj][ii], out[name][key]["band-max"][jj][ii]])


                                for jj in range(0,len(labs)):
                                        #axs[i].errorbar(x, ys[jj], yerr=yerrs[jj], label=labs[jj])
                                        axs[i].fill_between(x, yerrs[jj][0,:], yerrs[jj][1,:], alpha=0.3)
                                        axs[i].plot(x, ys[jj], '-o', label=labs[jj])

                                #axs[i].legend()
                                axs[i].set_xticks(range(1,nruns+1,2))
                                axs[i].set_title(name)

                        axs[-1].legend(loc=[1,0.8])
                        if ninds==1:
                                plt.savefig(os.path.join(dir,"figures/all", "%s.pdf" % scen))
                                plt.savefig(os.path.join(dir,"figures/all", "%s.png" % scen))
                        else:
                                plt.savefig(os.path.join(dir,"figures/all", "%s_%i.pdf" % (scen,ii)))
                                plt.savefig(os.path.join(dir,"figures/all", "%s_%i.png" % (scen,ii)))
        #plt.show()
        plt.close('all')


def plot_tuning(out):
        """
        Plot the tuning runs A3 and A5 for all
        """

def plot_one_band_Ns(out, band=1):
        """
        Plot all model N into one plot for one band.
        """
        labs = ["bottom", "middle", "top"]
        npart = len(participants(out))
        for scen in scens:
                nruns,inds,timescale = scens[scen]
                ninds = len(inds)
                for ii in range(ninds): # loop over all plotting times.
                        fig, ax = plt.subplots(1, 1, sharey=True, sharex=True, figsize=(12, 4))
                        ax.set_ylabel("N (MPa)")
                        time = out['bf']["%s%i" % (scen,1)]['t']
                        if inds[0]!=0 or ninds>1:
                                ax.set_xlabel("Scenario %s# at time %i" % (scen, time[inds[ii]]/timescale) )
                        else:
                                ax.set_xlabel("Scenario %s#" % scen)
                        for i,name in enumerate(participants(out)): #enumerate(['mw']): #
                                x = np.zeros(nruns)
                                ys = [np.nan*np.zeros(nruns), np.nan*np.zeros(nruns), np.nan*np.zeros(nruns)]
                                yerrs = [np.nan*np.zeros([2,nruns]), np.nan*np.zeros([2,nruns]), np.nan*np.zeros([2,nruns])]
                                for run in xrange(0,nruns):
                                        key = "%s%i" % (scen,run+1)
                                        if not out[name].has_key(key): continue
                                        x[run] = run+1
                                        jj = band
                                        ys[jj][run] = out[name][key]["band-means"][jj][ii]
                                        yerrs[jj][:,run] = np.array([out[name][key]["band-min"][jj][ii], out[name][key]["band-max"][jj][ii]])



                                #ax.errorbar(x, ys[jj], yerr=yerrs[jj], label=labs[jj])
                                jj = band
                                ax.fill_between(x, yerrs[jj][0,:], yerrs[jj][1,:], alpha=0.3)
                                ax.plot(x, ys[jj], '-o', label=name)

                                #ax.legend()
                                ax.set_xticks(range(1,nruns+1))
                                ax.set_title("At %s band" % labs[band])

                        ax.legend(loc=[1,0.1])
                        if ninds==1:
                                plt.savefig(os.path.join(dir,"figures/all", "%s_band%i.pdf" % (scen,band)))
                        else:
                                plt.savefig(os.path.join(dir,"figures/all", "%s_%i_band%i.pdf" % (scen,ii,band)))
        #plt.show()
        plt.close('all')



if __name__=="__main__":
        dir = "/home/mauro/projects/hydro_intercomparsion/shmip-files"
        out = readall(dir)
        # band_all(out)
        # plot_band_Ns(out)
        # plot_one_band_Ns(out, band=0)
        # plot_one_band_Ns(out, band=1)

#!/usr/bin/env python

###############################################################################
# Copyright (C) 2017-2019 Potsdam-Institute for Climate Impact Reasearch (PIK), 
# Author: Torsten Albrecht (albrecht@pik-potsdam.de)
# License: GNU AFFERO GENERAL PUBLIC LICENSE version 3
# 
# This script aggregates the total score and plots the results.
###############################################################################


import os
import netCDF4 as nc

import matplotlib.pyplot as plt
from matplotlib import cm, colors, rcParams, colorbar
import numpy as np
import csv

#on pik cluster
#source activate python_for_pism_calib

## this hack is needed to import config.py from the project root
#project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)


class compare_score(object):

  def __init__(self,enums,paleo,snp):

    self.ensnums = enums
    self.paleo = paleo
    self.showplots = snp[0]
    self.printtopdf = snp[1]
    self.printout = snp[2]

    self.res = cf.resolution #km resolution
    #runtime=cf.runtime #yrs
    self.outpath = cf.output_data_path

    self.param_dict_file = cf.param_dict_file
    #self.slvol0 = 64.4475 #m old (wrong) PISM calculation
    #self.slvol0=56.801 #m 15km boot
    #self.slvol0=58.528 #m 15km
    #self.slvol0 = 56.8412 #m 16km
    self.slvol0 = 56.85 #m 16km

    self.normscore = cf.normscore
    self.grvol = cf.groundedvol

    self.msp = cf.maxscoreplot

    #get all individual scores
    self.best_scores = []
    self.best_score_ind = cf.best_score_ind
    self.score = self.get_ensemble_score(self.ensnums)
    self.score_names = cf.score_names

    #get the individual combination of parameter values for each ensemble member
    self.param_ens,self.param_dict = self.get_parameter_values_for_ensemble(self.param_dict_file)

    ### plot paramter - score analysis (for 4 parameter dimensions defined)
    if self.showplots or self.printtopdf:

      #clean figure data
      for fignum in plt.get_fignums():
        plt.close(fignum)

      min_enum=min(sorted(self.score.keys()))
      max_enum=max(sorted(self.score.keys()))
      self.printname = self.outpath+"plots/placeholder_"+str(min_enum)+"-"+str(max_enum)+"_"+str(int(self.res))+"km.pdf"
      if paleo:
        self.printname = self.printname.replace(".nc","_paleo.nc")

      #plot slvol timeseries weighted by score
      self.plot_score_slvol(self.score)

      #plot scores as scatter for each metric
      self.plot_score_scatter_metric(self.score,self.param_dict,self.param_ens)

      #plot scores in paramter matrix 
      #(first two paramters as inner panel matrices, third and fourth paramter as outer matrix dimension)
      self.plot_score_over_parameter_space(self.score,self.param_dict,self.param_ens)

      #plot mean score for each of the four parameters 
      self.plot_mean_score_per_parameter(self.score,self.param_dict,self.param_ens)

      #plot combined score for each pair combination of the four parameters
      self.plot_score_over_parameter_pairs(self.score,self.param_dict,self.param_ens)

      #if self.showplots:
      plt.show()



  #############################################################
  def get_ensemble_score(self,ensnums):


    stats_dict={}
    ### read misfits for individual dta types from txt-file
    for i,enum in enumerate(ensnums):
        statfile = self.outpath+"stats/le_ens"+str(enum)+"_"+str(int(self.res))+"km.txt"

        try:
            savestat = open(statfile, 'r')
            savestat.readline() #skip measure titles
         
            #savestat.readline().rstrip("\n").split("         ")
            #  if scname in self.score_names:

            scores = [np.float(i) for i in savestat.readline().rstrip("\n").split(" ")[1:]]

            if self.paleo: #add paleo scores
              savestat.readline() #skip measure titles
              for paleo_score in savestat.readline().rstrip("\n").split(" ")[1:]:
                scores.append(np.float(paleo_score))
            savestat.close()

            stats_dict[str(enum)] = scores

            #if str(enum) in ['6165']: #,'6164','6152','6173','6125','6200','6217']:
            #  print enum,scores

        except:
            print 'No statsfile for ensemble member '+str(enum)


    ###  get median
    median_misfits = np.nanmedian(stats_dict.values(),axis=0)
    #median_misfits = np.nanmean(stats_dict.values(),axis=0)

    ### get individual scores weighted by median
    single_scores = {k: np.exp(-np.array(v)/median_misfits) for k, v in stats_dict.iteritems()}

    for k, v in stats_dict.iteritems():
      print k,v[0:2],median_misfits[0:2],single_scores[k][0:2]


    ### inter data weighting ######################

    ### get total score as the product of individual scores (no inter-data-type weighting here)
    #self.best_score_ind="TOTAL":
    if self.normscore:
      total_scores = {k: [np.prod(v)] for k, v in single_scores.iteritems()}
      self.best_scores = [ k for k, v in sorted(total_scores.items(), key=lambda l: l[1], reverse=True)]
      #print self.best_scores
      total_scores = {k: v / total_scores[self.best_scores[0]][0] for k, v in total_scores.iteritems()} 
    else:
      total_scores = {k: [np.prod(v)**(1.0/len(scores))] for k, v in single_scores.iteritems()}
      self.best_scores = [ k for k, v in sorted(total_scores.items(), key=lambda l: l[1], reverse=True)]
    #self.score_names.index(self.best_score_ind)
    #print self.best_scores

    #for enum in ['6165','6164','6152','6173','6125','6200','6217']:
    #  print single_scores[enum]
    #  print enum,total_scores[enum]

    ### get the total score weighted by inter-data-weights according to Briggs et al., 2014
    #idws = np.ones(len(scores))
    #idws = np.array([0.0375,0.0375,0.225,0.075,0.075,0.1,0.1,0.22,0.1,0.03]) # np.sum(idws)=1
    #single_scores_w = {k: np.exp(-np.array(v*idws)/median_misfits) for k, v in stats_dict.iteritems()}
    #total_scores  = {k: [np.prod(v)] for k, v in single_scores_w.iteritems()}
    #total_scores  = {k: [np.prod(v*idws)] for k, v in single_scores_w.iteritems()}

    #self.best_scores = [ k for k, v in sorted(total_scores.items(), key=lambda l: l[1], reverse=True)]

    ### normalize ######################

    ### normalize to probability sum(P)=1
    #total_scores = {k: [v / sum(total_scores.values())] for k, v in total_scores.iteritems()}

    ### normalize to median
    #total_scores = {k: [v / (np.exp(-1.0))**(len(scores))] for k, v in total_scores.iteritems()}

    ### normalize to best total score
    #best_score = np.nanmax(total_scores.values(),axis=0) #[0]
    #total_scores = {k: v / best_score for k, v in total_scores.iteritems()}

    stattotfile = self.outpath+"stats/le_all06"+"_"+str(int(self.res))+"km.txt"
    if os.path.isfile(stattotfile):
      os.system("rm "+stattotfile)
    savestat = open(stattotfile, 'a')
    for ind in self.best_scores:
        #wrtxt = ind+","+str(total_scores[ind][0])+"\n"
        #print ind,single_scores[ind],len(single_scores[ind])

        wrtxt = ind+","+str(total_scores[ind][0])
        for indscore in single_scores[ind]:
          wrtxt += ","+str(indscore)
        wrtxt+="\n"

        savestat.write(wrtxt)
    savestat.close()


    #with open(stattotfile,'wb') as f:
      #w = csv.writer(f)
      #w.writerow(total_scores.keys())
      #w.writerow(total_scores.values())


    ### add total score to list of individual scores
    total_scores = {k: list(np.concatenate((v,single_scores[k]))) for k, v in total_scores.iteritems()}

    return total_scores

  ###################################################################
  def get_parameter_values_for_ensemble(self,pfile):

    param_dict={}
    param_ens={}
    param_names=[]
    param_diff={}

    with open(pfile) as f:
      for l,line in enumerate(f.readlines()):

        param_vals={}
        params=line.rstrip('\n').split('"')

        #header with parameter names
        if l==0:
          for pn in params[0].split(' '):
            #if pn!='ens_member':
              param_names.append(pn)

        #parameter values in lines below
        else:
          #FIXME: if numbers instead of hashes
          ensmem=str(params[1].split(' ')[1])
          for pi,pn in enumerate(params[2].split(' ')):
            #if pn!='':
              try:
                param_vals[param_names[pi]]=np.float(pn)
              except:
                "no paramter val"
                #  print "no paramter val for "+param_names[pi]
              #corrections:
              if param_names[pi]=='visc':
                param_vals[param_names[pi]]=np.float(pn)*1e21

          #span param_dict
          if l==1:
            param_diff=param_vals
          else:
            for pn,pval in param_vals.items():
              if pval!=param_diff[pn]:

                if pn not in param_dict.keys():
                  param_dict.setdefault(pn, [])

                if pval not in param_dict[pn]:
                  param_dict[pn].append(pval)

          param_ens[ensmem]=param_vals

    #add paramter vals from first run
    for pn,pvals in param_dict.items():
      #param_dict[pn].append(param_diff[pn])
      param_dict[pn].insert(0,param_diff[pn])

    print param_dict
    #print param_ens

    return (param_ens,param_dict)



  def plot_score_slvol(self,ens_scr):

      #import analyse_timeseries as at; reload(at)
      #import glob
      #ts_path_names = os.path.join(cf.resultpath,cf.pism_file_name.replace(str(cf.fillnum),str("*")))
      #ensemble_members = glob.glob(ts_path_names)
      #ts_data = at.Ensemble_TimeSeries(ensemble_members)

      #print ens_scr # dict of ensemble members and individual score(s)

      #fig24, [ax24a,ax24b] = plt.subplots(2, 1,figsize=(7, 6),num=24) #sharex='col', sharey='row'
      fig24, ax24b = plt.subplots(1, 1,figsize=(9, 4),num=24) #sharex='col', sharey='row'

      ax24a = ax24b.twinx()
      #ax24a.grid('on')

      plotgradient=False
      plottopresent=False

      if self.grvol:
        plt.title("grounded ice volume anomaly")
        ax24a.set_ylabel("(PD obs - mod) in mio km3")
        volunit="mio. km3"
        volscale=1.0e15
        varname="ice_volume_glacierized_grounded"
        printname = self.printname.replace("placeholder","grvolall")
        slvol0 = 26.2914 #16km
        varmax = 10
      else:
        #plt.title("equivalent global-mean sea-level contribution (ESL)")
        if plottopresent:
          ax24a.set_ylabel("(PD - ESL) in m sea-level equivalent")
        else:
          #ax24a.set_ylabel("(PD obs - ESL) in m sea-level equivalent")
          ax24a.set_ylabel("sea-level anomaly in m")
          ax24b.set_ylabel("sea-level relevant volume in m")
        volunit="m SLE"
        volscale=1.0
        varname="slvol"
        printname = self.printname.replace("placeholder","slvolall")
        slvol0 = self.slvol0
        varmax = 20

      ax24a.plot([-125,0],[0,0],color="k",linewidth=1,linestyle="dashed")
      #ax24a.axvline(-14.35,color="k",linewidth=1,linestyle="dotted")
      #ax24a.set_xlabel("time [kyr]")
      ax24a.set_xlabel("time in kyr")
      #x24a.axis([-125,0,-varmax,varmax])
      ax24a.axis([-125,0,-varmax/2,varmax])
      ax24a.set_xticks(np.arange(-120,varmax,varmax))
      
      if True: #include penultimate glacial cycle
        ax24a.plot([-210,0],[0,0],color="k",linewidth=1,linestyle="dashed")
        #ax24a.axis([-210,0,-varmax,varmax])
        ax24a.axis([-210,0,-varmax/2,varmax])
        ax24a.set_xticks(np.arange(-200,25,25))
        #print cf.start_from_file
      
      
      if True: #show deglaciation period
        #ax24a.plot([-15,0],[0,0],color="k",linewidth=1,linestyle="dashed")
        #ax24b.plot([-15,0],[slvol0,slvol0],color="k",linewidth=1,linestyle="dashed")
        #ax24a.axis([-15,0,-varmax/2,varmax])
        ax24b.axis([-15,0,-varmax/2+slvol0,varmax+slvol0]) 
        #ax24a.set_xticks(np.arange(-14,2,2))
        ax24b.set_xticks(np.arange(-14,2,2))
        if plotgradient:
          ax24a.axis([-15,0,-0.02,0.02])
      

      #ax24a.yaxis.grid(True)
      #ax24a.grid(True, which='both')
      ax24b.set_ylim([-varmax/2+slvol0,varmax+slvol0])

      for i in [2529,2533,2534,2535]:
          #initfile='/p/tmp/albrecht/pism18/pismOut/forcing/'
          #initfile='/p/projects/pism/albrecht/pism_paleo_ensemble/ensemble_2400-2540/'
          #TODO
          initfile=cf.resultpath2+'forcing'+str(i)+'_TPSO/results/ts_forcing_'+str(cf.resolution)+'km_'+str(cf.runtime)+'yrs.nc'
          if os.path.exists(initfile):
            tsdata = nc.Dataset(initfile, 'r')
            slvol = np.squeeze(tsdata.variables[varname][:])/volscale
            sltime = np.squeeze(tsdata.variables["time"][:])/cf.seconds_per_year*1e-3
            tsdata.close()
            if plottopresent:
              slvol0=slvol[-1]
            if plotgradient:
              ax24a.plot(sltime,np.gradient(slvol),color=cf.colorscheme[1],linewidth=1,zorder=2,alpha=0.8)
            else:
              #ax24a.plot(sltime,slvol0-slvol,color=cf.colorscheme[1],linewidth=1,zorder=0,alpha=0.8)
              ax24a.plot(sltime[0:85000],slvol[0:85000]-slvol0,color=cf.colorscheme[1],linewidth=1,zorder=2,alpha=0.8)
            #print i,sltime[-5000],slvol[-5000]

      ### aggreagated score normalized to probability sum(P)=1
      agg_score = {k: v[0] for k, v in ens_scr.iteritems()}
      agg_score_prop = {k: v / sum(agg_score.values()) for k, v in agg_score.iteritems()}

      #best_score = np.nanmax(agg_score.values(),axis=0)
      best_score_mem = self.best_scores[0]
      best_score = agg_score[best_score_mem]
      worst_score = agg_score[self.best_scores[-1]]

      sl_data = {}

      #loop over ensemble members/scores
      for l,(enum,scores) in enumerate(sorted(ens_scr.items())):
        
        ### read each netcdf timeseries file ###############
        ts_file_name = cf.pism_file_name.replace(str(cf.fillnum),str(enum))
        tsfile = os.path.join(cf.resultpath,ts_file_name)
        tsfile = tsfile.replace('paleo.nc','timeseries.nc') #'timeseries-slvol.nc'
        
        if os.path.exists(tsfile):
          tsdata = nc.Dataset(tsfile, 'r')
          slvol = np.squeeze(tsdata.variables[varname][:])/volscale
          sltime = np.squeeze(tsdata.variables["time"][:])/cf.seconds_per_year*1e-3
          tsdata.close()
          sl_data[enum]=slvol
          if plottopresent:
              slvol0=slvol[-1]

          ### plot sl-curve with alpha according to score
          alphamin = 0.03
          alphamax = 0.9
          alphaval = alphamax+((alphamax-alphamin)*(scores[0]-best_score))/(best_score-worst_score)
          if plotgradient:
            ax24a.plot(sltime,np.gradient(slvol),color="k",linewidth=0.5,alpha=alphaval,zorder=1)
          else:
            #ax24a.plot(sltime,slvol0-slvol,color="k",linewidth=0.5,alpha=alphaval,zorder=1)
            ax24a.plot(sltime[::10],slvol[::10]-slvol0,color="k",linewidth=0.5,alpha=alphaval,zorder=1)

          ### best/reference run #############################
          #if enum=='1170':#reference
          if enum == best_score_mem:
          #if enum in self.best_scores[0:3]:
            if plotgradient:
              ax24a.plot(sltime,np.gradient(slvol),color=cf.colorscheme[3],linewidth=2,zorder=3,label="best")
            else:
              ax24a.plot(sltime,slvol-slvol0,color=cf.colorscheme[3],linewidth=1.5,zorder=3,label="best")
              #ax24a.plot(sltime,slvol0-slvol,color=cf.colorscheme[3],linewidth=1.5,zorder=3,label="best")

          ### ensemble mean ###################################
          if l==0:
            ensemble_mean_slvol = slvol * agg_score_prop[enum]
            tslen = len(slvol)
            sltime0 = sltime
            ensemble_stddev_slvol = np.zeros_like(slvol)
          else:
            ensemble_mean_slvol += slvol[-tslen:] * agg_score_prop[enum]
      if not plotgradient:
        if plottopresent:
          slvol0=ensemble_mean_slvol[-1]
        ax24a.plot(sltime0,ensemble_mean_slvol-slvol0,color=cf.colorscheme[2],linewidth=2.0,zorder=2,label="mean")
        #ax24a.plot(sltime0,slvol0-ensemble_mean_slvol,color=cf.colorscheme[2],linewidth=2.0,zorder=2,label="mean")

      ### ensemble standard deviation ########################
      for l,(enum,sl) in enumerate(sorted(sl_data.items())):
          ensemble_stddev_slvol += (ensemble_mean_slvol - sl[-tslen:])**2 * agg_score_prop[enum]
      ensemble_stddev_slvol = np.sqrt(ensemble_stddev_slvol)
      #stnd_upper = slvol0-ensemble_mean_slvol+ensemble_stddev_slvol
      #stnd_lower = slvol0-ensemble_mean_slvol-ensemble_stddev_slvol
      stnd_upper = ensemble_mean_slvol-slvol0+ensemble_stddev_slvol
      stnd_lower = ensemble_mean_slvol-slvol0-ensemble_stddev_slvol
      if not plotgradient:
        ax24a.fill_between(sltime0[::10],stnd_upper[::10],stnd_lower[::10],color=cf.colorscheme[2],zorder=0,alpha=0.4,label="stdev")
        #ax24a.plot(sltime0,stnd_upper,color=cf.colorscheme[2],linewidth=2,zorder=2,linestyle="dotted",label="stdev")
        #ax24a.plot(sltime0,stnd_lower,color=cf.colorscheme[2],linewidth=2,zorder=2,linestyle="dotted")

      legend24 = ax24a.legend(loc="lower left", shadow=False, fontsize=11)
      rcParams['legend.frameon'] = 'False'

      ### min and max values ##################################
      pm = str(r"$\pm$")
      pd_mean_sl   = ensemble_mean_slvol[-1]
      #pd_span      = str(np.around(stnd_lower[-1],decimals=1))
      #pd_span     += " - "+str(np.around(stnd_upper[-1],decimals=1))+" \n "+volunit
      pd_span      = str(np.around(pd_mean_sl-slvol0,decimals=1)) + pm 
      pd_span      += str(np.around(ensemble_stddev_slvol[-1],decimals=1))#+" \n "+volunit

      lgm_mean_sl  = np.max(ensemble_mean_slvol)
      lgm_mean_k   = np.argmax(ensemble_mean_slvol)
      lgm_mean_k2  = np.argmin((sltime0+15)**2)
      lgm_mean_t   = sltime0[lgm_mean_k]
      #lgm_span     = str(np.around(stnd_lower[lgm_mean_k2],decimals=1))
      #lgm_span    += " - "+str(np.around(stnd_upper[lgm_mean_k2],decimals=1))+" \n "+volunit
      lgm_span      = str(np.around(ensemble_mean_slvol[lgm_mean_k2]-slvol0,decimals=1)) + pm 
      lgm_span     += str(np.around(ensemble_stddev_slvol[lgm_mean_k2],decimals=1)) #+" \n "+volunit

      eeam_mean_sl = np.min(ensemble_mean_slvol[-tslen:-tslen/2])
      eeam_mean_k  = np.argmin(ensemble_mean_slvol[-tslen:-tslen/2])
      eeam_mean_k2 = np.argmin((sltime0+120)**2)
      eeam_mean_t  = sltime0[eeam_mean_k]
      #eeam_span    = str(np.around(stnd_lower[eeam_mean_k2],decimals=1))
      #eeam_span   += " - "+str(np.around(stnd_upper[eeam_mean_k2],decimals=1))+" \n "+volunit
      eeam_span     = str(np.around(ensemble_mean_slvol[eeam_mean_k2]-slvol0,decimals=1)) + pm 
      eeam_span     += str(np.around(ensemble_stddev_slvol[eeam_mean_k2],decimals=1)) #+" \n "+volunit

      degl_mean_k10 = np.argmin((sltime0+10)**2)
      degl_mean_k5 = np.argmin((sltime0+5)**2)


      #print pd_mean_sl,lgm_mean_sl,lgm_mean_t,eeam_mean_sl,eeam_span,eeam_mean_t
      #print slvol0,ensemble_mean_slvol[-1],ensemble_stddev_slvol[-1]

      print "Reconstructions of sea-level contributions:"
      print "period","mean sl","std sl","time","mean anomaly"
      print "\nPD:",pd_mean_sl,ensemble_stddev_slvol[-1],sltime0[-1],pd_mean_sl-slvol0
      print "LGM:",lgm_mean_sl,ensemble_stddev_slvol[lgm_mean_k],lgm_mean_t,lgm_mean_sl-slvol0
      print "LIG:",eeam_mean_sl,ensemble_stddev_slvol[eeam_mean_k],eeam_mean_t,eeam_mean_sl-slvol0,"\n"
      print "LGM2:",ensemble_mean_slvol[lgm_mean_k2],ensemble_stddev_slvol[lgm_mean_k2],sltime0[lgm_mean_k2],ensemble_mean_slvol[lgm_mean_k2]-slvol0
      print "LIG2:",ensemble_mean_slvol[eeam_mean_k2],ensemble_stddev_slvol[eeam_mean_k2],sltime0[eeam_mean_k2],ensemble_mean_slvol[eeam_mean_k2]-slvol0,"\n"
      print "DEG10:",ensemble_mean_slvol[degl_mean_k10],ensemble_stddev_slvol[degl_mean_k10],sltime0[degl_mean_k10],ensemble_mean_slvol[degl_mean_k10]-slvol0
      print "DEG5:",ensemble_mean_slvol[degl_mean_k5],ensemble_stddev_slvol[degl_mean_k5],sltime0[degl_mean_k5],ensemble_mean_slvol[degl_mean_k5]-slvol0,"\n"



      #ax24a.text(eeam_mean_t,6.0*varmax/8.0,eeam_span,color=cf.colorscheme[2],zorder=3,fontsize=10)
      #ax24a.text(-35,varmax/2.0,lgm_span,color=cf.colorscheme[2],zorder=3,fontsize=10)
      #ax24a.text(-20,3.0*varmax/4.0,pd_span,color=cf.colorscheme[2],zorder=3,fontsize=10)

      ax24a.text(eeam_mean_t,-3.0*varmax/8.0,eeam_span,color=cf.colorscheme[2],zorder=3,fontsize=11)
      ax24a.text(-30,7.0*varmax/8.0,lgm_span,color=cf.colorscheme[2],zorder=3,fontsize=11)
      ax24a.text(-21,-3.0*varmax/8.0,pd_span,color=cf.colorscheme[2],zorder=3,fontsize=11)


      if self.printtopdf:
        #printname = self.printname.replace("placeholder","slvolall")
        plt.savefig(printname, format='pdf')
        plt.savefig(printname.replace(".pdf",".png"), format='png',dpi=300)



  def plot_score_over_parameter_space(self,ens_scr,par_space,ens_par):

      #print ens_scr # dict of ensemble members and individual score(s)
      #print ens_par #dict of ensemble members and individual paramter values
      #print par_space #dict of paramter diminsions and possible values

      #loop over ensemble members/scores
      pos_per_param={}
      for enum,escr in sorted(ens_scr.items()):

        #match parameter combination of ensemble member with location in parameter space
        pos_par_array=[]
        for pnum,pname in enumerate(par_space):
          try:
            pos_par_array.append(par_space[pname].index(np.float(ens_par[enum][pname])))
          except:
            print "Position not found for "+str(enum)
        pos_per_param[enum]=[pos_par_array,escr]
      #print pos_per_param 
 
      #find parameter matrix position of ensemble member
      alen=len(par_space[list(par_space)[0]])
      blen=len(par_space[list(par_space)[1]])
      clen=len(par_space[list(par_space)[2]])
      dlen=len(par_space[list(par_space)[3]])
      #print alen,blen,clen,dlen
      
      par_field=np.zeros([clen*dlen,alen,blen])
      par_names=list(np.zeros_like(par_field))

      for enum,point in sorted(pos_per_param.items()):

        apos=point[0][0]
        bpos=point[0][1]
        cpos=point[0][2]
        dpos=point[0][3]
        epos=cpos*clen+dpos

        par_field[epos,apos,bpos]=point[1][0] #total score
        par_names[epos][apos][bpos]=np.int(enum)

        #print top3 ensemble members
        if enum in self.best_scores[0:5]:
          print "best score members: \n"
          print enum, point[1], ens_par[enum], "\n"

      #plot matrix
      fig20 = plt.figure(20,figsize=(8, 9))

      for sub in xrange(clen*dlen):
        ax20=plt.subplot(clen,dlen,sub+1)
        if self.normscore:
          cs20 = ax20.pcolor(par_field[sub],cmap=cm.YlOrRd,norm=colors.LogNorm(vmin=0.01, vmax=self.msp))
          tcks = np.array([0,0.02,0.05,0.1,0.2,0.4,0.8])
        else:
          cs20 = ax20.pcolor(par_field[sub],cmap=cm.hot_r,vmin=0,vmax=self.msp)
          tcks = np.arange(0,self.msp,0.2)

        #write ensemblenumber to matrix location
        for i in xrange(alen):
         for j in xrange(blen):
           enumstr = str(np.int(par_names[sub][j][i]))
           if enumstr in self.best_scores[0:3]:
             ax20.text(i+0.15,j+0.4,enumstr,color="w",fontsize=8)
           else: #bright backround
             ax20.text(i+0.15,j+0.4,enumstr,color="k",fontsize=8)

        #inner ticks, shifted by 0.5
        plt.xticks(range(blen), par_space[list(par_space)[1]], fontsize=8)
        plt.yticks(range(alen), par_space[list(par_space)[0]], fontsize=8)

        ax20.xaxis.set(ticks=np.arange(0.5, blen), ticklabels=par_space[list(par_space)[1]])
        ax20.yaxis.set(ticks=np.arange(0.5, alen), ticklabels=par_space[list(par_space)[0]])
        
        #outer labels
        if (sub)/clen>=dlen-1:
          ax20.set_xlabel(par_space[list(par_space)[3]][sub%dlen])
        if sub%dlen==0: 
        #if sub%3==0:
          ax20.set_ylabel(par_space[list(par_space)[2]][(sub+1)/dlen])
        else: 
          ax20.set_yticks([])
        
      #define colorbar
      #cbaxes = fig20.add_axes([0.97, 0.1, -0.04, 0.4])
      cbaxes = fig20.add_axes([0.97, 0.1, -0.05, 0.4])
      cb = plt.colorbar(cs20,cax=cbaxes,ticks=tcks,format='%.2f') #,extend='max')
      cb.set_label('total score',multialignment="left")
      cb.outline.set_linewidth(0)

      ### print plot to pdf file
      if self.printtopdf:
        printname = self.printname.replace("placeholder","paramscore")
        plt.savefig(printname, format='pdf')
        plt.savefig(printname.replace(".pdf",".png"), format='png',dpi=300)


  def plot_mean_score_per_parameter(self,ens_scr,par_space,ens_par):

      num_of_params=len(par_space) #number of parameter dimensions: 4
      alen=len(par_space[list(par_space)[0]]) #assuming the same number of vals for each paramter dimension
      num_of_mems=len(ens_scr) #number of ensemble members
     
 
      mean_score=np.zeros([num_of_params,alen])
      stnd_score=np.zeros([num_of_params,alen])

      for i,(num,scr) in enumerate(ens_scr.iteritems()):
        for j,(pn,par) in enumerate(par_space.iteritems()):
          for k,pval in enumerate(par):
            if np.float(ens_par[num][pn])==np.float(pval):
              mean_score[j,k]+=scr[0] #total score
      mean_score /= ((num_of_mems+1)/alen)
      mean_score /= mean_score.sum(axis=1) #normalize to %


      for i,(num,scr) in enumerate(ens_scr.iteritems()):
        for j,(pn,par) in enumerate(par_space.iteritems()):
          for k,pval in enumerate(par):
            if np.float(ens_par[num][pn])==np.float(pval):
              stnd_score[j,k]+=(scr[0]-mean_score[j,k])**2.0
      stnd_score = np.sqrt(stnd_score/(num_of_mems/alen))
      #print stnd_score


      mean_per_par=np.zeros([num_of_params])
      stnd_per_par=np.zeros([num_of_params])

      for j,(pn,par) in enumerate(par_space.iteritems()):
        for k,pval in enumerate(par):
          mean_per_par[j]+=pval*mean_score[j,k]
      mean_per_par /= np.sum(mean_score,axis=1)

      for j,(pn,par) in enumerate(par_space.iteritems()):
        for k,pval in enumerate(par):
          stnd_per_par[j]+=mean_score[j,k]*(pval-mean_per_par[j])**2.0
      stnd_per_par = np.sqrt(stnd_per_par/np.sum(mean_score,axis=1))

      fig21 = plt.figure(21,figsize=(7, 8))
      bin_width=np.zeros(num_of_params)

      for sub,(pn,par) in enumerate(par_space.iteritems()):

        diffax=(max(par)-min(par))*0.2
        minax=min(par)-diffax
        maxax=max(par)+diffax
        bin_width[sub]=(maxax-minax)/10.0
        if sub==0:
          print "mean per parameter value:"
        print "values:",par_space[pn],"\nmeans:",mean_score[sub],"\nsigmas:",stnd_score[sub],"\n"
        ax21=plt.subplot(num_of_params,1,sub+1)
        plt.bar(par_space[pn],mean_score[sub],bin_width[sub],color='grey',edgecolor="k",alpha=0.8)
        #plt.errorbar(par_space[pn],mean_score[sub],stnd_score[sub],fmt='.',markersize='10',color='b',capsize=4, elinewidth=1,alpha=0.4) #linestyle='None', marker='^'
        plt.errorbar(mean_per_par[sub],np.mean(mean_score[sub,:]),xerr=stnd_per_par[sub],yerr=0.0,fmt='.',markersize='10',color='r',capsize=5, elinewidth=1)    
        ax21.set_xticks(par_space[pn])

        ax21.set_ylabel(pn)
        #ax21.set_yticks([0,0.25,0.5,0.75,1.0])
        if self.normscore:
          ax21.set_yticks(np.arange(0,0.6,0.1))
          ax21.axis([minax,maxax,0.0,0.5])
          #ax21.axis([minax,maxax,0.0,0.25])
          #ax21.set_yticks(np.arange(0,1.2,0.2))
          #ax21.axis([minax,maxax,0.0,1.0])
          #ax21.grid("on")
        else:
          ax21.set_yticks(np.arange(0,self.msp,0.2))
          ax21.axis([minax,maxax,0.0,0.6])

        #ax21.axhline(np.exp(-1.0),color='g',linestyle='dashed',alpha=0.4,linewidth=0.5) #median

      ### print plot to pdf file
      if self.printtopdf:
        printname = self.printname.replace("placeholder","parammean")
        plt.savefig(printname, format='pdf')
        plt.savefig(printname.replace(".pdf",".png"), format='png',dpi=300)


  def plot_score_over_parameter_pairs(self,ens_scr,par_space,ens_par):

      num_of_params=np.int(len(par_space)) #number of parameter dimensions: 4
      alen=np.int(len(par_space[list(par_space)[0]])) #assuming the same number of vals for each paramter dimension
      num_of_mems=np.int(len(ens_scr)) #number of ensemble members

      num_of_pairs=np.int(num_of_params*(num_of_params-1)/2.0)
      pairs=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]

      #loop over ensemble members/scores
      pos_per_param={}
      for enum,escr in sorted(ens_scr.items()):

        #match parameter combination of ensemble member with location in parameter space
        pos_par_array=[]
        for pnum,pname in enumerate(par_space):
         
          pos_par_array.append(par_space[pname].index(np.float(ens_par[enum][pname])))
        pos_per_param[enum]=[pos_par_array,escr]

      #mean score for each pair of parameter combinations
      par_field=np.zeros([num_of_pairs,alen,alen])
      for i,pair in enumerate(pairs):
        par1=pair[0]
        par2=pair[1]

        for enum,point in sorted(pos_per_param.items()):
          k=point[0][par1]
          l=point[0][par2]
          par_field[i,k,l]+=ens_scr[enum][0] #total score
      par_field/=(alen)**2

      #print np.max(par_field)
      #par_field/=np.max(par_field) #normalize to best pair

      #for i,pf in enumerate(par_field):
      #  par_field[i]/=np.sum(pf)
      par_field/=(np.sum(par_field)/num_of_pairs) #normalize to sum 1 for each pair

      fig22 = plt.figure(22,figsize=(8, 8))
      for sub,pair in enumerate(pairs):

        #locations in plat matrix
        subpos=(-pair[1]+(num_of_params-1))+(num_of_params-1)*(pair[0])+1
        ax22=plt.subplot(num_of_params-1,num_of_params-1,subpos)
        if self.normscore:
          cs22 = ax22.pcolor(par_field[sub],cmap=cm.YlOrRd,norm=colors.LogNorm(vmin=0.01, vmax=self.msp)) #cm.gist_heat_r
          tcks = np.array([0,0.02,0.05,0.1,0.2,0.4,0.8])
        else:
          cs22 = ax22.pcolor(par_field[sub],cmap=cm.hot_r,vmin=0,vmax=self.msp)
          tcks = np.arange(0,self.msp,0.2)

        pname0=list(par_space)[pair[0]]
        pname1=list(par_space)[pair[1]]

        #ticks
        plt.xticks(range(alen), par_space[pname1], fontsize=8)
        plt.yticks(range(alen), par_space[pname0], fontsize=8)
        ax22.xaxis.set(ticks=np.arange(0.5, alen), ticklabels=par_space[pname1])
        ax22.yaxis.set(ticks=np.arange(0.5, alen), ticklabels=par_space[pname0])
        if subpos in [2,3,5]:
          ax22.yaxis.set(ticks=[])

        #labels
        if subpos in [3,5,7]:
          ax22.set_xlabel(pname1)
        if subpos in [1,4,7]:
          ax22.set_ylabel(pname0)
   
      #colorbar
      #cbaxes = fig22.add_axes([0.77, 0.1, -0.04, 0.4])
      cbaxes = fig22.add_axes([0.77, 0.1, -0.05, 0.4])
      cb = plt.colorbar(cs22,cax=cbaxes,ticks=tcks,format='%.2f') #,extend='max')
      cb.set_label('mean score',multialignment="left")
      cb.outline.set_linewidth(0)

      ### print plot to pdf file
      if self.printtopdf:
        printname = self.printname.replace("placeholder","parampairs")
        plt.savefig(printname, format='pdf')
        plt.savefig(printname.replace(".pdf",".png"), format='png',dpi=300)


  def plot_score_scatter_metric(self,ens_scr,par_space,ens_par):

      #print par_space
      #print ens_scr
      #print ens_par
      
      num_of_params=np.int(len(par_space)) #number of parameter dimensions: 4
      num_of_metrics=np.int(len(ens_scr[list(ens_scr)[0]])) #total and individual scores: 5(PD)+5(Paleo)=11
      #print num_of_params,num_of_metrics

      best_score_mem = self.best_scores[0]
      best_score_mem2 = self.best_scores[1]
      best_score_mem3 = self.best_scores[2]
      print best_score_mem,best_score_mem2,best_score_mem3,self.best_scores[3],self.best_scores[4]

      score_names = self.score_names

      fig23 = plt.figure(23,figsize=(int(1.6*num_of_metrics), int(1.6*num_of_params)))
      plt.clf()
      #plt.title("Aggregated Score        Present Day Scores        Paleo Scores")

      mean_scr_param=np.zeros([num_of_metrics,num_of_params,4])


      for i,(num,scr) in enumerate(ens_scr.iteritems()): #ensemble members
        for k,metric in enumerate(scr):
          for l,(pn,par) in enumerate(par_space.iteritems()):

            subpos=k+(l)*num_of_metrics+1
            #print k,l,subpos,metric
            ax23=plt.subplot(num_of_params,num_of_metrics,subpos)
            if i==0: #need to be set only once
              if k==0:
                ax23.set_ylabel(pn)
                ax23.set_yticks([0,0.5,1])
                ax23.tick_params(axis='y', labelsize=7)
              else:
                ax23.set_yticks([])
              if not self.normscore:
                ax23.axhline(np.exp(-1.0),color='g',linestyle='dashed',alpha=0.2,linewidth=0.5,zorder=0) #median
              if l==num_of_params-1:
                ax23.set_xlabel(score_names[k])

            for j,pval in enumerate(par):
              if pval==ens_par[num][pn]:
                #print i,num,k,metric,l,j,pval,subpos
                ax23.plot(pval,metric,"k.",markersize=2,zorder=1,alpha=0.5) #score points
                if i==0: # set axis already for the first esemble member
                  diffax=(max(par)-min(par))*0.2
                  minax=min(par)-diffax
                  maxax=max(par)+diffax
                  #ax23.axis([minax,maxax,0,1.0]) #1.25
                  ax23.axis([minax,maxax,0,1.05]) #1.25
                  ax23.set_xticks(par_space[pn])
                  ax23.tick_params(axis='x', labelsize=6)

                if best_score_mem == num: #show best run scores
                  ax23.plot(pval,metric,"r.",markersize=9,zorder=4)
                elif best_score_mem2 == num: #show second best run scores
                  ax23.plot(pval,metric,"g.",markersize=9,zorder=3,alpha=0.5)
                elif best_score_mem3 == num: #show third best run scores
                  ax23.plot(pval,metric,"b.",markersize=9,zorder=2,alpha=0.3)
                #elif ens_par[num]['sia_e']==7.0: 
                #  ax23.plot(pval,metric,"y.",markersize=9,zorder=2,alpha=0.3)

                #calculate score means per parameter values
                #m = np.where( par==pval )
                m = par.index( pval )
                mean_scr_param[k,l,m]+=metric

      mean_scr_param/=(np.float(len(ens_scr))/4.0)
      for k,metric in enumerate(scr):
        for l,(pn,par) in enumerate(par_space.iteritems()):
          subpos=k+(l)*num_of_metrics+1
          ax23=plt.subplot(num_of_params,num_of_metrics,subpos)
          ax23.plot(par,mean_scr_param[k,l,:],linestyle="dashed",color="grey",alpha=0.5,zorder=0) #,marker="x")


      if self.printtopdf:
        printname = self.printname.replace("placeholder","scorescatter")
        plt.savefig(printname, format='pdf')
        plt.savefig(printname.replace(".pdf",".png"), format='png',dpi=300)

#!/usr/bin/env python

###############################################################################
# Copyright (C) 2017-2019 Potsdam-Institute for Climate Impact Reasearch (PIK), 
# Author: Torsten Albrecht (albrecht@pik-potsdam.de)
# License: GNU AFFERO GENERAL PUBLIC LICENSE version 3
# 
# This script calls the ensemble analysis. You can set flags below.
###############################################################################

do_paleo    = True
cal_individual_misfits = True

show_plots  = False
save_to_pdf = True
print_info  = True


###############################################################################
# encoding=utf8


import numpy as np
import matplotlib.pyplot as plt
import os

import config as cf; reload(cf)
import tools as tl; reload(tl)
import get_score as gs; reload(gs)
import get_paleo as ps; reload(ps)
import compare_score as cs; reload(cs)

### all ensemblenumbers over which ensemble shall be analyzed

#ens_exp=2
ens_min = cf.ens_min
ens_max = cf.ens_max

####################################################################################
ensnums=np.arange(ens_min,ens_max+1)
shownprint = [show_plots,save_to_pdf,print_info]

if save_to_pdf:
  if not os.path.exists(cf.output_data_path+"plots"):
    os.makedirs(cf.
output_data_path+"plots")

### load observations (Bedmap2, Rignot velocities)
obsfile      = cf.obsfile 
velobsfile   = cf.velobsfile 
#observations = [xobs,yobs,Mx,My,mobs,Bobs,Hobs,hobs,lonobs,latobs,cellarea,velobs]
observations = tl.load_observations(obsfile,velobsfile)

### calculate the calibration measures
if cal_individual_misfits:
  
  for ensnum in ensnums:

    print "ensemble number: "+str(ensnum)

    #get the PD scores and write to stats/txt-file
    pds = gs.pd_score(ensnum,observations,shownprint)

    #get the paleo scores and write to stats/txt-file
    if do_paleo:
      pls = ps.paleo_score(ensnum,observations,shownprint)

    if show_plots:
      plt.show()

### get the total score(from stats/txt-files) and plot analysis 
ens = cs.compare_score(ensnums,do_paleo,shownprint)
score = ens.score

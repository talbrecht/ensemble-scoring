
###############################################################################
# Copyright (C) 2017-2019 Potsdam-Institute for Climate Impact Reasearch (PIK), 
# Author: Torsten Albrecht (albrecht@pik-potsdam.de)
# License: GNU AFFERO GENERAL PUBLIC LICENSE version 3
# 
# This file contains user defined paths and settings.
###############################################################################

import os
import pwd


#settings
resolution=16 #km resolution
runtime=210000 #yrs

ens_min=6000
ens_max=6255

ensemblenumber='06'

approach="a"
#approach="b"

normscore=True

groundedvol=False

maxscoreplot=1.0

best_score_ind = "TOTAL"

score_names = [best_score_ind,"TOTE","TOTI","TOTDH","TOTVEL","TOTGL","TOTUPL","TROUGH","ELEV","EXT"] #"RSL","DSLV"
# THROUGH only works with extra_paleo.nc files, which did not enter the data publication for memory limitation

###############################################################################

#constants
seconds_per_year = 365.0*24.0*3600.0
#seconds_per_year = 3.1556926e7
rhoi  = 910.0 #PISM default
rhosw = 1028.0
rhofw = 1.0e3
km3_to_msle = 361.0e3*rhosw/rhoi

#mask values
moc=4 #mask value open ocean
mfl=3 #mask value floating
mgr=2 #mask value grounded
mif=0 # mask value ice-free bedrock


colorscheme=["#9C7A18","#532C80","#166B28","#630718","#496588","#8ee6c1"]
#brown, violet, green, red, bluegrey, cyan


###############################################################################

username = pwd.getpwuid(os.getuid()).pw_name
author= username+"@pik-potsdam.de"

basedir="/p/tmp/albrecht/paleo_ensemble/"
output_data_path = os.path.expanduser(basedir+'out/')

resultpath = basedir+'datapub/model_data2/' #TODO
#resultpath = '/p/projects/pism/albrecht/paleo_ensemble/'


fillnum=9999
pism_file_name = "pism1.0_paleo"+ensemblenumber+"_"+str(fillnum)+"/paleo.nc"

param_dict_file = basedir+"datapub/model_data/aggregated_data/pism1.0_paleo"+ensemblenumber+"_"+str(ens_min)+".csv"

inputpath      = basedir+'init_data/'

obsfile       = inputpath+"pism_boot/result_boot_"+str(resolution)+"km_nowostok.nc"

velobsfile    = inputpath+"rignotvel_initmip"+str(resolution)+"km_filled.nc"

gps_data_file = basedir+"datapub/model_data/aggregated_data/gps/gps_whitehouse_"+str(resolution)+"km.p"

trans_data_file = basedir+"datapub/model_data/aggregated_data/gps/transect_lonlat_"+str(resolution)+"km.p"

raised_file   = inputpath+"Raised/raised_glmask_"+str(resolution)+"km_pollard.nc"

briggs_data   = basedir+"datapub/model_data/aggregated_data/briggs/paleodata_"+str(resolution)+"_briggs.p"

briggs_input  = inputpath+"Briggs"



#### No edits needed below that line. ####
#output_file_name = os.path.join(resultpath,output_file_name)
project_root = os.path.dirname(os.path.abspath(__file__))


#!/usr/bin/env python

###############################################################################
# Copyright (C) 2017-2019 Potsdam-Institute for Climate Impact Reasearch (PIK), 
# Author: Torsten Albrecht (albrecht@pik-potsdam.de)
# License: GNU AFFERO GENERAL PUBLIC LICENSE version 3
# 
# This script executes the ensemble analysis for paleo data.
###############################################################################


import os, pickle
import matplotlib.pyplot as plt
from matplotlib import cm, colors, rcParams
import numpy as np

#on pik cluster
#source activate python_for_pism_calib

## this hack is needed to import config.py from the project root
#project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import tools as tl; reload(tl)


class paleo_score(object):

  def __init__(self,enum,obs,snp):

    self.ensnum = str(enum)
    self.showplots = snp[0]
    self.printtopdf = snp[1]
    self.printout = snp[2]

    self.res = cf.resolution #km resolution
    #runtime=cf.runtime #yrs

    # do gps calculations or load from pickle file
    self.calc_gps=True
    if os.path.exists(cf.gps_data_file):
      self.calc_gps=False

    self.calc_trans=True
    if os.path.exists(cf.trans_data_file):
      self.calc_trans=False

    self.calc_briggs=True
    if os.path.exists(cf.briggs_data):
      self.calc_briggs=False

    self.score_names = cf.score_names
    self.pd_score_names =["TOTUPL","TROUGH","ELEV","EXT","RSL"] #available paleo measures
    self.score_text=""

    self.appr = cf.approach

    # pathnames ################################################################
    #workpath     = cf.workpath
    resultpath   = cf.resultpath
    self.outpath = cf.output_data_path

    self.gps_data_file = cf.gps_data_file
    self.trans_data_file = cf.trans_data_file
    self.raised_file = cf.raised_file

    self.briggs_outfile = cf.briggs_data
    self.briggs_compilation = cf.briggs_input

    pism_file_name  = cf.pism_file_name.replace(str(cf.fillnum),str(enum))
    pismfile        = os.path.join(resultpath,pism_file_name)
    pism_extra_file = pismfile.replace("paleo.nc","extra_paleo.nc")

    global x,y,Mx,My,mobs,Bobs,lonobs,latobs #only the one used in the functions below
    [x,y,Mx,My,mobs,Bobs,Hobs,hobs,lonobs,latobs,cellarea,velobs,velstnd] = obs


    if os.path.exists(pismfile):

      ### reconstructgion of PD uplift rates in different sites of WAIS
      uplcalc = self.get_uplcalc(tl.get_data(pismfile,['dbdt']))

      ### paleo data (grounding line reconstructions along transects, surface elevation)
      troughcalc = self.get_troughcalc(pismfile)

      ### AntICEdat
      briggsdat = self.import_briggs()
      extime = tl.get_data(pism_extra_file,['time'])
      thk,topg,usurf = tl.get_data(pism_extra_file,['thk','topg','usurf'])

      elevcalc = self.get_elevcalc(briggsdat['elev'],extime,[thk,topg])
      extcalc = self.get_extcalc(briggsdat['ext'],extime,[thk,topg,usurf])
      rslcalc = self.get_rslcalc(briggsdat['rsl'],extime,[topg,usurf])

      score_choice={} #indicated which measures to consider
      score_choice["TOTUPL"] = uplcalc
      score_choice["TROUGH"] = troughcalc
      score_choice["ELEV"]   = elevcalc
      score_choice["EXT"]    = extcalc
      score_choice["RSL"]    = rslcalc

      ### write scores to txt file 
      #scores = [uplcalc,throughcalc,elevcalc,extcalc,rslcalc]
      self.score_text=""
      scores=[]
      for pdsc in self.pd_score_names:
        if pdsc in self.score_names:
          self.score_text += "         "+pdsc
          scores.append(score_choice[pdsc])

      self.print_to_txtfile(scores)

    else:
      print "\nThere is no PISM resultfile "+pismfile+" !\nSkip..."
      print "\n###############################################\n"


  ##############################################################
  ### after Pollard et al., 2016
  ### Large ensemble modeling of the last deglacial retreat of the West 
  ### Antarctic Ice Sheet: comparison of simple and advanced statistical techniques, 
  ### Geosci. Model Dev., 9, 1697-1723, 
  ### https://doi.org/10.5194/gmd-9-1697-2016
  ################################################################
  def get_uplcalc(self,var):


    #     Calcs rms difference of modern uplift rates from each obs site (zrms), 
    #     using closest model grid pt (iupl,jupl).

    #     1. Whitehouse et al,GJI, Table S2, whitehouse_table_s2.dat,tabupl.
    #        https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fj.1365-246X.2012.05557.x&file=GJI_5557_sm_SupportingInformation.zip 

    #FIXME: use intra-data-type weighting

    dbdt=var[0] #uplift rate variable
    scalemm=1.0 #1000.0 im mm

   
    gps_data={}
    gps_data[1]=['FTP1',-78.93,162.57,2.14,2.77]
    gps_data[2]=['ROB1',-77.03,163.19,7.54,2.59]
    gps_data[3]=['TNB1',-74.7,164.1,-0.23,0.81]
    gps_data[4]=['MCM4_AV',-77.85,166.76,0.77,0.31]
    gps_data[5]=['MBL1_AV',-78.03,204.98,3.28,1.09]
    gps_data[6]=['W01_AV',-87.42,210.57,-2.8,1.17]
    gps_data[7]=['MBL2',-76.32,215.7,0.23,4.13]
    gps_data[8]=['MBL3',-77.34,218.13,0.09,1.97]
    gps_data[9]=['W09',-82.68,255.61,4.54,2.59]
    gps_data[10]=['W06A',-79.63,268.72,-2.2,2.42]
    gps_data[11]=['W07_AV',-80.32,278.57,3.61,1.58]
    gps_data[12]=['W05_AV',-80.04,279.44,4.86,1.01]
    gps_data[13]=['HAAG',-77.04,281.71,3.47,0.71]
    gps_data[14]=['W08A/B',-75.28,287.82,1.31,1.28]
    gps_data[15]=['W02_AV',-85.61,291.45,2.17,1]
    gps_data[16]=['OHIG',-63.32,302.1,3.8,1]
    gps_data[17]=['PALM',-64.78,295.95,0.08,1.87]
    gps_data[18]=['ROTB',-67.57,291.87,1.5,1.9]
    gps_data[19]=['SMRT',-68.12,292.9,-0.22,1.93]
    gps_data[20]=['FOS1',-71.31,291.68,2.14,0.4]
    gps_data[21]=['BREN',-72.67,296.97,3.85,1.6]
    gps_data[22]=['W04_AV',-82.86,306.8,3.42,0.84]
    gps_data[23]=['BELG',-77.86,325.38,2.97,1.47]
    gps_data[24]=['W03_AV',-81.58,331.6,-2.47,1.28] 
    gps_data[25]=['SVEA',-74.58,348.78,2.07,1.95]

    gps_data[26]=['ABOA',-73.04,346.59,1.4,0.84]
    gps_data[27]=['VESL',-71.67,357.16,1.06,0.45]
    gps_data[28]=['MAIT_AV',-70.77,11.84,0.12,0.4]
    gps_data[29]=['SYOG',-69.01,39.58,2.26,0.36]
    gps_data[30]=['MAW1',-67.61,62.87,0.06,0.39]
    gps_data[31]=['A368',-74.29,66.79,0.39,1.0]
    gps_data[32]=['A351',-72.91,74.91,0.81,1.28]
    gps_data[33]=['DAV1',-68.58,77.97,-0.94,0.49]
    gps_data[34]=['CAS1',-66.28,110.52,1.18,0.43]
    gps_data[35]=['DUM1',-66.67,140.0,-0.79,0.46]

    #after Larsen-B break-off 2002
    #gps_data[16]=['OHIG_EL1',-63.32,302.1,6.78,0.83]
    #gps_data[16]=['OHIG_EL2',-63.32,302.1,3.22,1.68]
    #gps_data[17]=['PALM_EL',-64.78,295.95,8.75,0.64]
    #gps_data[18]=['ROTB_EL',-67.57,291.87,6.96,0.66]
    #gps_data[19]=['SMRT_EL',-68.12,292.9,3.15,1.94]

    #intra data weighting
    lonlat_dat={}
    for p,val in gps_data.items():
    #coordinates of gps locations
        lonlat_coords={}
        lonlat_coords['lon']=val[2]
        lonlat_coords['lat']=val[1]
        lonlat_dat[p]=lonlat_coords

    dlon=10.0 #degree
    dlat=5.0 #degree
    self.add_intra_weights(lonlat_dat,dlon,dlat)

    #radius = 6.371220e6
    #distcrit = radius * (5./180.)
    iupl=np.zeros(len(gps_data)+1,dtype=int)
    jupl=np.zeros_like(iupl,dtype=int)
    idw=np.zeros(len(gps_data)+1)

    # find gps location on the grid
    if self.calc_gps:

      if not os.path.exists(self.outpath+"gps"):
          os.makedirs(self.outpath+"gps")

      #intra data weighting
      lonlat_dat={}
      for p,val in gps_data.items():
      #coordinates of gps locations
        lonlat_coords={}
        lonlat_coords['lon']=val[2]
        lonlat_coords['lat']=val[1]
        lonlat_dat[p]=lonlat_coords

      dlon=10.0 #degree
      dlat=5.0 #degree
      self.add_intra_weights(lonlat_dat,dlon,dlat)

      gps_save=gps_data
      for p,val in gps_data.items():
        #coordinates of gps locations
        iupl[p],jupl[p] = tl.get_closest_point(val[2],val[1],lonobs,latobs)

        gps_save[p].append(np.int(iupl[p]))
        gps_save[p].append(np.int(jupl[p]))
        #print lonlat_dat[p]['weight']
        idw[p] = lonlat_dat[p]['weight']
        gps_save[p].append(idw[p])

      #print gps_data
      print "dump to gps pickle file"
      pickle.dump( gps_save, open( self.gps_data_file, "wb" ) )

    else:

      #print "read pickle file"
      gps_save = pickle.load( open( self.gps_data_file, "rb" ) )

      for p,val in gps_save.items():
        iupl[p] = np.int(val[5])
        jupl[p] = np.int(val[6])
        idw[p] = val[7]

    zdif = 0.
    nrms = 0
    zrms = np.zeros_like(iupl,dtype=float)
    weidat = np.ones_like(iupl,dtype=float)
    totupl = 0.


    for p,val in gps_data.items():
        if self.appr =="a":
          weidat[p]=idw[p]*(val[4])**(-2.0)
        elif self.appr =="b":
          weidat[p]=1.0  #approach B : no uncertainty, no intra-data-weighting

    for p,val in gps_data.items():
        zdif = dbdt[iupl[p],jupl[p]]*scalemm-np.float(val[3])
        zrms[p] = (zdif**2)*weidat[p]


    if self.appr =="a":
      totupl = np.nanmean(zrms)
      printline= '\nUPL (misfit uplift GPS stations)'
      #print np.sum(idw) #sum of weights should equal unity

    elif self.appr =="b":
      totupl = np.sqrt(np.nanmean(zrms))
      printline= '\nUPL (rmse uplift GPS stations in mm/yr)'

    if self.printout:
      print printline
      print totupl
    
    #return (totupl)

    ###############################################################
    if self.showplots or self.printtopdf:  #FIXME: use basemap

      ###map of PD upflift rates compared to PGS locations #######
      fig12 = plt.figure(12,figsize=(9, 9))
      plt.clf()
      ax12=plt.subplot(111)
      ax12.axis("equal")
      ax12.axis([0,Mx,0,My])
      ax12.axis("off")

      #map of uflift rate
      dbdtm=np.ma.array(dbdt,mask = mobs == 4)
      #dbdt.mask = mobs  == 4
      relbd=20
      ticks12=np.arange(-relbd,relbd+5,5)
      cs1  = ax12.contourf(dbdtm*scalemm,ticks12,alpha=1.0,cmap=cm.PuOr_r,extend='both')#cm.RdBu_r
      ## colorbar dbdt
      #cbaxes = fig12.add_axes([0.15, 0.25, -0.04, 0.5])
      cbaxes12 = fig12.add_axes([0.15, 0.23, 0.3, 0.02])
      #ticks = np.arange(int(mn),int(mx),4)
      #cb = plt.colorbar(cs1,cax=cbaxes,ticks=np.arange(-30,35,5))
      ticks12b=np.arange(-relbd,relbd+10,10)
      cb = plt.colorbar(cs1,cax=cbaxes12, orientation="horizontal",ticks=ticks12b)
      cb.set_label('rel. bed elevation change [mm/yr]',multialignment="left")
      cb.outline.set_linewidth(0)


      #bedmap2 contour
      cs7 = ax12.contour(mobs,[2.5],colors="w",linewidths=1.0,alpha=1.,)

      #GPS locations
      #cmrb = cm.get_cmap('RdBu_r')
      cmrb = cm.get_cmap('PuOr_r')
      colors=[cmrb(1.*i/len(ticks12)) for i in range(len(ticks12)+1)]
      for p,val in gps_data.items():
        colorind=0
        for ti,t in enumerate(ticks12):
          if val[3]>t:
            colorind=min(ti+1,len(ticks12))
        ax12.plot(jupl[p],iupl[p],'ro',markerfacecolor=colors[colorind],markersize=20)
        ax12.text(jupl[p]-1,iupl[p]+2,str(val[3]))
        #if p>1:
        # ax12.plot([jupl[p-1],jupl[p]],[iupl[p-1],iupl[p]],"r-")

      ### print plot to pdf file
      printname = self.outpath+"plots/gpsmap_ensemble"+self.ensnum+"_"+str(int(self.res))+"km.pdf"
      if self.printtopdf:
        plt.savefig(printname, format='pdf')
        #plt.savefig(printname.replace(".pdf",".png"), format='png',dpi=300)


      ###GPS upflift rates lined-up and compared to PISM results ###############
      fig13 = plt.figure(13,figsize=(13, 10))
      plt.clf()
      ax13=plt.subplot(111)
      ax13.axis([0,len(gps_data)+1,-15,15])
      ax13.set_xticks(np.arange(1,len(gps_data)+1,1))
      ax13.set_title("GPS data at locations around WAIS")

      #rates with standard deviation
      for p,val in gps_data.items():
        ax13.plot(p,val[3],"bo",label="gps obs")
        ax13.plot([p,p],[val[3]-val[4],val[3]+val[4]],"k-")
        ax13.plot(p,dbdt[iupl[p],jupl[p]]*scalemm,"ro",label="uplift model")

        #draw a line and legend
        if p==1:
          legend13 = ax13.legend(loc='lower left', shadow=True, fontsize=11)
        else: 
          ax13.plot([p-1,p],[dbdt[iupl[p-1],jupl[p-1]]*scalemm,dbdt[iupl[p],jupl[p]]*scalemm],"r")

      ax13.text(1,12,"Ross west \n(1-4)")
      ax13.text(4.5,12,"Ross east \n(5-8)")
      ax13.text(8,12,"WAIS divide \n(9)")
      ax13.text(12,12,"Ronne \n(10-14)")
      ax13.text(17,12,"Peninsula \n(15-20)")
      ax13.text(22,12,"Filchner \n(21-24)")
      ax13.text(30,12,"EAIS \n(25-37)")

      ax13.axhline(0,linestyle="dashed")

      ### print plot to pdf file
      printname = self.outpath+"plots/gpsanomaly_ensemble"+self.ensnum+"_"+str(int(self.res))+"km.pdf"
      if self.printtopdf:
        plt.savefig(printname, format='pdf')
        #plt.savefig(printname.replace(".pdf",".png"), format='png',dpi=300)

      #plt.show()

    return (totupl)



  ##############################################
  def get_troughcalc(self,fn):

    rhor = cf.rhoi/cf.rhosw

    dkm=10 #km steps
    #raised_data = pickle.load( open( raised_data_file, "rb" ) )
    point_names=["Ross","Pine Island","Ronne","Amery"]

    ### define transects as lon,lat coordinates
    points_ross=np.array([[-112.5,-79.8],[-119,-80.8],[-152.8,-82.35],[-166.2,-82.3],[-175.1,-81.4],[-178.2,-79.9],[-175.9,-77.8],[-173.6,-73.7]])
    points_pine=np.array([[-89.1,-75.4],[-96.7,-75.4],[-99.8,-75.2],[-103.7,-74.7],[-106.3,-74.2],[-107.4,-73.3],[-107.8,-72.0],[-107.0,-70.2]])
    points_ronne=np.array([[-91.5,-82.8],[-65.7,-82.0],[-62.1,-81.4],[-58.84,-80.77],[-56.5,-79.4],[-55.5,-77.9],[-53.9,-76.3],[-49.4,-75.1],[-42.0,-72.7]])
    #points_ronne=np.array([[-91.5,-82.8],[-65.7,-82.0],[-62.1,-81.4],[-59.3,-81.45],[-56.5,-79.4],[-55.5,-77.9],[-53.9,-76.3],[-49.4,-75.1],[-42.0,-72.7]])
    points_amery=np.array([[63.5,-76.8],[66.4,-75.4],[65.8,-74.4],[66.1,-73.8],[67.6,-72.85],[68.6,-72.1],[70.2,-71.1],[71.5,-69.9],[73.7,-68.7],[70.6,-66.4]])

    #points on 20km grid
    #points_ross=[[89,119],[96,116],[121,103],[130,99],[136,93],[138,85],[135,74],[130,51]] #ross, gl [121,103]
    #points_pine=[[60,141],[61,131],[60,126],[59,120],[57,116],[53,113],[46,110],[36,108]] #pine island, gl [60,126]
    #points_ronne=[[101,139],[101,158],[101,161],[100,164],[92,172],[86,177],[80,184],[78,193],[77,210]] #ronne, gl [101,158]
    #points_amery=[[204,172],[213,172],[218,175],[221,176],[227,176],[231,176],[237,175],[245,175],[252,173],[263,183]] #amery, gl [227,176]

    #find closest points on grid
    if self.calc_trans:
      if not os.path.exists(self.outpath+"gps"):
          os.makedirs(self.outpath+"gps")

      transects=[]
      for xi,tr in enumerate([points_ross,points_pine,points_ronne,points_amery]):
        transx=[]
        for yi,lonlat in enumerate(tr):
          ip,jp = tl.get_closest_point(lonlat[0],lonlat[1],lonobs,latobs)
          transx.append([jp,ip])
        transects.append(transx)

      print "dump to trans coordinates pickle file"
      pickle.dump( transects, open( self.trans_data_file, "wb" ) )

    else:

      #print "read pickle file"
      transects = pickle.load( open( self.trans_data_file, "rb" ) )

    #raised_snaps=[20,15,10,5,0]
    raised_snaps=[0,5,10,15,20]
    raised_snaps=[0,20]
    raised_option="A"
    raised_mask={}

    ### scale to 15km grid
    for snap in raised_snaps:
      if snap>0:
        raised_mask[snap]=tl.get_data(self.raised_file,["raised_"+raised_option+"_"+str(snap)+"ka"])[0]
      else:
        raised_mask[snap]=tl.get_data(self.raised_file,["bedmap2_"+str(snap)+"ka"])[0] #mobs
    

    ### get points along transects first ########################################
    trans_points = tl.get_points_along_transect(x,y,dkm,transects)

    glp=np.zeros([len(raised_snaps),len(trans_points)])
    glp0=np.zeros_like(glp)
    gldiff=np.zeros([len(raised_snaps),len(trans_points)])
    sea_level=np.zeros(len(raised_snaps))

    ######################################################################
    if self.showplots or self.printtopdf:  #FIXME: use basemap

      #clean figure window
      if 16 in plt.get_fignums():
        #fig16.clf()
        plt.close(16)
        plt.close(17)

      fig16, ax16all = plt.subplots(2, 2,figsize=(13, 7),num=16) #sharex='col', sharey='row'
      #fig16.canvas.set_window_title('Figure16')

      fig17, ax17 = plt.subplots(figsize=(9, 7),num=17)

      #not fill in legend
      rcParams['legend.frameon'] = 'False'

      Bobsm=np.ma.array(Bobs)
      Bobsm.mask = (Bobs<-2500)
      cs1  = ax17.contourf(Bobsm,np.arange(-2000,0,100),alpha=0.8,cmap=cm.Blues_r,extend='both')
      cs7 = ax17.contour(mobs,[2.5,3.5],colors="w",linewidths=.9,alpha=1.0)

      cbaxes = fig17.add_axes([0.15, 0.25, -0.04, 0.5])
      #ticks = np.arange(int(mn),int(mx),4)
      cb = plt.colorbar(cs1,cax=cbaxes,ticks=np.arange(-2000,0,200))
      cb.set_label('bed elevation [m]',multialignment="left")
      cb.outline.set_linewidth(0)

      ax17.axis("equal")
      ax17.axis([0,Mx,0,My])
      ax17.axis("off")

      for j,points in enumerate(transects):
        for i,p in enumerate(points):
          if i>0:
            ax17.plot([points[i-1][0],p[0]],[points[i-1][1],p[1]],color="#777777")


    ###########################################################################
    #### find snapshot index for raised data comparison
    for count,snap in enumerate(raised_snaps):
      time_raised=snap*(-1000) #years

      ### load snapshots
      #snap_add=str(time_raised).zfill(8)
      snap_add=str(time_raised)
      if time_raised==0:
        sfn=fn
      else:
        sfn=fn.replace("paleo.nc","snapshots_")+snap_add+".000.nc"
      #sfn=fn.replace("result_","snap_")+"_"+snap_add+".000.nc"

      m,B,H,h = tl.get_data(sfn,['mask','topg','thk','usurf'])

      #print np.shape(H),np.shape(h),np.shape(B)

      rm=raised_mask[snap]
      ##################################################

      if self.showplots or self.printtopdf:

        cs8 = ax17.contour(m,[2.5],colors=cf.colorscheme[count],linewidths=1.2,alpha=1.,labels=str(snap))
        if snap>0:
          cs9 = ax17.contour(rm,[0.5],colors=cf.colorscheme[count],linewidths=0.8,alpha=0.8,linestyles="dotted")

      ###############################################################
      ### loop over four troughs
      for l,po in enumerate(trans_points):

        Bav=np.zeros(len(po))
        hav=np.zeros(len(po))
        Hav=np.zeros(len(po))
        mav=np.zeros(len(po))
        hab=np.zeros(len(po))

        rav=np.zeros(len(po))

        sea_level[count]=h[0,0] #sea_level
        glfound=False
        gl0found=False

        ### loop over equidistant points along trough
        for k,p in enumerate(po):

          i=np.int(np.floor(p[1]))
          j=np.int(np.floor(p[0]))
          di=p[0]-np.floor(p[0])
          dj=p[1]-np.floor(p[1])

          #interpolation
          Bav[k] = tl.get_interpolated_value(i,j,di,dj,B)
          hav[k] = tl.get_interpolated_value(i,j,di,dj,h)
          Hav[k] = tl.get_interpolated_value(i,j,di,dj,H)
          rav[k] = tl.get_interpolated_value(i,j,di,dj,rm)
          #mav[k] = get_interpolated_value(i,j,di,dj,m)
          mav[k]=m[np.int(i),np.int(j)]
          hab[k]=rhor*Hav[k]-(sea_level[count]-Bav[k])

          #find point on transects closest to raised data
          if k>0 and gl0found==False:
            if rav[k]>0.0 and rav[k-1]==0.0:
              glp0[count,l]=k-0.5
              #print snap,point_names[l],i,j,k
              ax17.plot(int(j),int(i),'bx',alpha=0.7,markersize=5)
              gl0found=True

          #check for grounding line
          if k>0 and glfound==False:
            if hab[k]<0.0 and hab[k-1]>=0.0:
              lam=hab[k]/(hab[k]-hab[k-1])
              lam=max(0.0,min(lam,1.0))
              glp[count,l]=k-lam
              #print snap,point_names[l],i,j,k-lam
              ax17.plot(int(j),int(i),'rx',alpha=0.7,markersize=5)
              glfound=True


        #no reconstructions available
        if point_names[l]=='Amery' and snap in [5,10,15]:
          glp0[count,l]=np.nan

        #reconstructions highly uncertain (Check Pollard)
        if point_names[l]=='Ross' and snap in [10,15]:
          glp0[count,l]=np.nan

        gldiff[count,l]=(dkm*(glp[count,l]-glp0[count,l]))**2
        #print snap,point_names[l],dkm*glp[count,l],dkm*glp0[count,l],np.sqrt(gldiff[count,l])

        ##################################################################

        if self.showplots or self.printtopdf:

          ax16=ax16all[l/2,l%2]

          #Normalizing
          kmdist=np.arange((-glp0[0,l])*dkm,(len(po)-glp0[0,l])*dkm,dkm)
          kmdist=kmdist[0:len(po)]
          #ice shelves
          ax16.plot(kmdist,hav,cf.colorscheme[count],label=str(snap)+" kyr BP")
          ax16.plot(kmdist,hav-Hav,cf.colorscheme[count])
          #bed
          ax16.plot(kmdist,Bav,"k",linewidth=1.0,alpha=0.6)


          #grounding line position along transect
          ax16.plot(kmdist[np.int(glp[count,l])],-1900,color=cf.colorscheme[count],marker="^")
          try: #if available in data
            ax16.plot(kmdist[np.int(glp0[count,l])],1900,color=cf.colorscheme[count],marker="v")
          except:
            donothing=1
            #print "No glp0 value "+point_names[l]+" for -"+str(snap)+" kyr"

          #ax16.fill_between(kmdist,hav-Hav,Bav,color="#000099")
          #ax16.fill_between(kmdist,Bav,-3000,color="#996633")
          #ax16.fill_between(kmdist,3000,hav,color="#e9e8e5")

          if count==0:
            ax16.axis([kmdist[0],kmdist[-1],-2000,2000])
            ax16.set_yticks(np.arange(-2000,2000,1000))
            ax16.set_xticks([-250,0,250,500,750])
            ax16.set_title(point_names[l])
            ax16.axhline(sea_level[count],color="k",linestyle="dashed",linewidth=0.5)
            if l>0:
              ax16.text(-250,-200,"mean: "+str(np.int(np.sqrt(np.nanmean(gldiff,axis=0)[l])))+" km")
            if l>1:
              ax16.set_xlabel("distance to PD GL in km", fontsize=11)

    #print glp0,glp

    #labeling
    ax16=ax16all[0,0]
    legend16 = ax16.legend(loc='upper right', fontsize=9)
    ax16.text(-600,-200,"mean: "+str(np.int(np.sqrt(np.nanmean(gldiff,axis=0)[0])))+" km")
    ax16.text(-600,-1900,"PISM GL position", fontsize=9)
    ax16.text(-600,1800,"RAISED GL position", fontsize=9)
   
    if self.printtopdf:
        printname = self.outpath+"plots/le_trough_sideview_ensemble"+self.ensnum+"_"+str(int(self.res))+"km.pdf"
        fig16.savefig(printname, format='pdf')


    #print gldiff
    #print np.sqrt(np.nanmean(gldiff,axis=0)),"mean over snaps per basin (Ross, PI, Ronne, Amery)"
    #print np.sqrt(np.nanmean(gldiff,axis=1)),"mean over basins per snap (PD left, LGM right)"

    #approach A
    if self.appr =="a":
      sigma_trough=30.0
      rmsgl=np.nanmean(gldiff/sigma_trough**2)
      printline= '\nTROUGH (misfit GL in 4 sections at LGM and PD snapshots)'
    #approach B
    elif self.appr =="b":
      rmsgl=np.sqrt(np.nanmean(gldiff))
      printline= '\nTROUGH (rmse GL in 4 sections at LGM and PD snapshots in km)'

    if self.printout:
      print printline
      print rmsgl

    return (rmsgl)

  ##############################################################
  ### AntICEdat
  ##############################################################
  def add_intra_weights(self,dat,dlon,dlat):

    relative_site_weight={}
    count_per_site={}
    dat_weight={}

    for ids in dat.keys():
      count_region=0
      count_site=0
      lon0=dat[ids]["lon"]
      lat0=dat[ids]["lat"]
      site_num=str(ids)[0:4]

      for ids2 in dat.keys():
        lon2=dat[ids2]["lon"]
        lat2=dat[ids2]["lat"]
        site_num2=str(ids2)[0:4]
        #fixme: average over four shifted lon-lat rectangles?!
        if np.abs(lon0-lon2) <= dlon and np.abs(lat0-lat2) <= dlat:
          count_region+=1
        if site_num==site_num2:
          count_site+=1

      relative_site_weight[site_num] = np.sqrt(count_site)/np.sqrt(count_region)
      count_per_site[site_num] = count_site 
    #print count_per_site,len(dat),relative_site_weight


    ### weighted per point
    total_site_weight = sum(relative_site_weight.values())
    for ids in dat.keys():
      site_num=str(ids)[0:4]
      dat_weight[ids] = relative_site_weight[site_num]/total_site_weight/count_per_site[site_num]
      dat[ids]["weight"] = np.around(dat_weight[ids],decimals=7)

    ### check sum of weights = 1
    #print sum(dat_weight.values())

    #for site,weight in relative_site_weight.items():
        #print site,count_per_site[site],weight,weight/total_site_weight

    return dat_weight


  ################################################################
  def read_xlsxfile(self,data_file,sheet_name,sheet_range,col_num,col_types):

    import openpyxl as px

    wb = px.load_workbook(data_file)
    ws = wb.get_sheet_by_name(name = sheet_name)

    col_id = col_num[0]
    col_ids = col_num[1]

    col_lat = col_num[2]
    col_lon = col_num[3]

    col_val = col_num[4]
    col_valu = col_num[5]

    col_age = col_num[6]
    col_ageu = col_num[7]

    col_type=col_num[8]
    col_proc=col_num[9]

    col_reg=col_num[10]
    col_ref=col_num[11]


    data_dict={}
    ##### get the data
    for rownum,rowval in enumerate(ws.iter_rows(sheet_range)):
      data = {}
      site_id = rowval[col_id].value
      site_type = rowval[col_type].value
      #print site_id,site_type,col_types,site_type in col_types

      #if site_id!=None:
      if site_type in col_types and site_id!=None: #isinstance(site_id,int)

        if rowval[col_ids].value != None:
          site_id=str(site_id)+str(rowval[col_ids].value)

        lon0 = np.float(rowval[col_lon].value)
        lat0 = np.float(rowval[col_lat].value)
        data['lon'] = lon0
        data['lat'] = lat0

        try:
          data['val'] = np.float(rowval[col_val].value)
          data['valu'] = np.float(rowval[col_valu].value)
        except:
          data['val'] = None
          data['valu'] = None

        data['age'] = np.float(rowval[col_age].value)
        data['ageu'] = np.float(rowval[col_ageu].value)

        data['type'] = str(rowval[col_proc].value)

        print sheet_name,rownum,site_id

        ilonlat,jlonlat = tl.get_closest_point(lon0,lat0,lonobs,latobs)
        data['j_grid'] = jlonlat
        data['i_grid'] = ilonlat

        di,dj = tl.find_interpolation(ilonlat,jlonlat,lon0,lat0,lonobs,latobs)

        data['di_grid'] = np.around(di,decimals=6)
        data['dj_grid'] = np.around(dj,decimals=6)

        data_dict[int(site_id)] = data

    return data_dict


  ################################################################
  def import_briggs(self):

      # get briggs data compilation from xlsx file
    if self.calc_briggs:

      if not os.path.exists(self.outpath+"briggs"):
          os.makedirs(self.outpath+"briggs")


      EXT_data_file = self.briggs_compilation+"/mmc3.xlsx"
      ext_data = self.read_xlsxfile(EXT_data_file,'marineCores',"A4:U36",[0,20,3,4,20,20,14,15,19,19,2,17],['M','G']) # 20 is none
      #print ext_data #[2101]

      RSL_data_file = self.briggs_compilation+"/mmc2.xlsx"
      rsl_data = self.read_xlsxfile(RSL_data_file,'rsl',"A3:AF111",[0,8,6,7,17,18,26,27,31,11,2,29],['D'])
      #print rsl_data #[91011]

      ELEV_data_file = self.briggs_compilation+"/mmc4.xlsx"
      elev_data = self.read_xlsxfile(ELEV_data_file,'paleoThickness',"A3:R147",[0,8,6,7,12,13,10,11,17,17,2,15],['D'])
      #print elev_data[11011]

      dlon=10.0 #degree
      dlat=5.0 #degree
      self.add_intra_weights(elev_data,dlon,dlat)
      self.add_intra_weights(rsl_data,dlon,dlat)
      
      briggs_data={}
      briggs_data['ext']=ext_data
      briggs_data['rsl']=rsl_data
      briggs_data['elev']=elev_data

      print "dump to gps pickle file"
      pickle.dump( briggs_data, open( self.briggs_outfile, "wb" ) )

    else:
      #print "read pickle file"
      briggs_data = pickle.load( open( self.briggs_outfile, "rb" ) )

    """
    if self.showplots:

      plt.rcParams['contour.negative_linestyle'] = 'solid'
      fig15 = plt.figure(15,figsize=(13, 10))
      plt.clf()
      ax=plt.subplot(111)
      im1 = ax.imshow(mobs,origin="lower",cmap=cm.Greys_r,vmin=0,vmax=4,interpolation="nearest")
      ax.contour(lonobs,np.linspace(-180,180,37),colors="k",linewidths=0.5)
      ax.contour(latobs,np.linspace(-90,90,37),colors="k",linewidths=0.5)

      for data_type,data_vals in briggs_data.items():

        if data_type=="ext":
          colorplot="ro"
          plt.text(10,10, '%s' % (data_type),color="r",fontsize=11)
        elif data_type=="rsl":
          colorplot="go"
          plt.text(10,30, '%s' % (data_type),color="g",fontsize=11)
        elif data_type=="elev":
          colorplot="bo"
          plt.text(10,50, '%s' % (data_type),color="b",fontsize=11)

        for ids,each in data_vals.items():
          j=each['i_grid']
          i=each['j_grid']
          ax.plot(i,j,colorplot,markersize=4)
          #plt.text(j+1,i+2, '%d' % (ids),color="k",fontsize=4)

      ax.axis("equal")
      ax.axis([0,Mx,0,My])
      ax.axis("off")
      plt.show()
    """

    return briggs_data


  ################################################################
  def get_elevcalc(self,dat,tex,vars):

      Hex=vars[0]
      Bex=vars[1]
      tex=tex[0]/cf.seconds_per_year
      snaps = np.shape(Hex)[0]

      s_struct=100.0
      mean_misfit=0

      for ids in dat.keys():
        #print ids,dat[ids]

        i       =  dat[ids]['i_grid']
        j       =  dat[ids]['j_grid']
        di      =  dat[ids]['di_grid']
        dj      =  dat[ids]['dj_grid']
        sigma_h =  (dat[ids]['valu'])**2 + s_struct**2
        dat_age =  -dat[ids]['age']
        dh      =  np.sqrt(2*sigma_h)
        dat_wgt =  dat[ids]['weight']

        H_int = tl.get_interpolated_value(i,j,di,dj,Hex)
        b_int = tl.get_interpolated_value(i,j,di,dj,Bex[-1]) #PD: sl and bad adjustment removed
        h_int = H_int + b_int
        h_thin = np.copy(h_int)
        #print b_int,np.shape(h_int),H_int[-1],h_int[-1]

        #misfit of surface elevation
        misfit_int = (h_int-dat[ids]["val"])**2/sigma_h

        misfit_all=np.zeros(snaps)

        misfitmin=10000.0
        h_at_age=0
        h_lgm = h_int[-1]


        # exclude periods of growing surface elevation
        for l,t in enumerate(reversed(tex)):
          lmax=len(tex)-1

          if h_int[lmax-l] + dh < h_lgm:
              h_thin[0:lmax-l] = 0.0
              break
          else:
              h_lgm = max(h_lgm,h_int[lmax-l])


        for l,t in enumerate(tex):

          #misfit over h and t
          if h_thin[l] > 0.0:
            misfit_all[l]=misfit_int[l]+((t-dat_age)/dat[ids]['ageu'])**2
          else:
            misfit_all[l]=np.NaN


          if t < 0.0:
            #time interpolation to datum
            if t < dat_age and tex[l+1] >= dat_age:
              t_before = dat_age - t
              t_after  = tex[l+1] - dat_age 
              h_at_age=(h_int[l]*t_after + h_int[l+1]*t_before) / (tex[l+1]-t)
              #print dat_age,t_before,t_after,t,tex[l+1],tex[l+1]-t,h_int[l],h_int[l+1],h_at_age


        misfitmin = np.nanmin(misfit_all)
        mint = tex[np.nanargmin(misfit_all)]
        mean_misfit += misfitmin*dat_wgt

        #print ids,mint,misfitmin

        """
        if self.showplots:
          fig16 = plt.figure(16,figsize=(8, 4));plt.clf()
          ax3=plt.subplot(111)
          ax3b=ax3.twinx()

          #vals and uncertainty
          ax3.axhline(dat[ids]['val'],color="k")
          ax3.axhline(h_lgm,color="b",linestyle="dotted")
          ax3.axhline(h_lgm+dh,color="b",linestyle="dashed",alpha=0.3)
          ax3.axhline(h_lgm-dh,color="b",linestyle="dashed",alpha=0.3)

          ax3.axvline(-dat[ids]['age'],color="k")
          ax3.axvline(-dat[ids]['age']-dat[ids]['ageu'],color="k",linestyle="dashed")
          ax3.axvline(-dat[ids]['age']+dat[ids]['ageu'],color="k",linestyle="dashed")

          #ax3.axvline(mint,color="b")
          ax3b.plot(tex[:],misfit_int[:],"m-")
          ax3b.plot(tex[:],misfit_all[:],"c-")
          ax3b.plot(mint,misfitmin,"ko",markersize=4)

          ax3.plot(tex[:],h_int[:],"g-")
          ax3.plot(tex[:],h_thin[:],"r.",markersize=2)
          #ax3.plot(tex[:],h_int[:]+dh,"y-",linestyle="dashed")
          #ax3.plot(tex[:],h_int[:]-dh,"y-",linestyle="dashed")
          ax3.plot(dat_age,h_at_age,"bo",markersize=4)

          ax3.plot(tex[:],Hex[:,i,j]+Bex[-1,i,j],"r-",alpha=0.2)
          ax3.plot(tex[:],Hex[:,i+np.int(np.sign(di)),j]+Bex[-1,i+np.int(np.sign(di)),j],"r-",alpha=0.2)
          ax3.plot(tex[:],Hex[:,i,j+np.int(np.sign(dj))]+Bex[-1,i,j+np.int(np.sign(dj))],"r-",alpha=0.2)
          ax3.plot(tex[:],Hex[:,i+np.int(np.sign(di)),j+np.int(np.sign(dj))]+Bex[-1,i+np.int(np.sign(di)),j+np.int(np.sign(dj))],"r-",alpha=0.2)

          ax3b.axis([-31000,0,0,100.0])
          ax3.axis([-31000,0,0,3500.0])
          plt.text(-30000, 50, '%d: score: %.3f, %d,%d' % (ids,np.sqrt(misfitmin),j,i))
          plt.show()
        """

      if self.printout:
          print '\nELEV (mean of minimum deviation in surface elevation and date)'
          print mean_misfit

      return (mean_misfit)


  ################################################################
  def get_extcalc(self,dat,tex,vars):
      
      Hex=vars[0]
      Bex=vars[1]
      sl=vars[2] # in 0,0 corner
      tex=tex[0]/cf.seconds_per_year
      snaps = np.shape(Hex)[0]

      retreat0=-35000.0
      s_struct_glr=250.0
      s_struct_omc=500.0
      mean_misfit=0.0
      count_ooc=0


      for ids in dat.keys():
        #print ids,dat[ids]

        i       =  dat[ids]['i_grid']
        j       =  dat[ids]['j_grid']
        di      =  dat[ids]['di_grid']
        dj      =  dat[ids]['dj_grid']
        sigma   =  (dat[ids]['ageu'])**2
        dat_age =  -dat[ids]['age']

        retreat=retreat0
        retreat_fine=retreat0

        mask_int=np.zeros(snaps)
        H_int = tl.get_interpolated_value(i,j,di,dj,Hex)
        B_int = tl.get_interpolated_value(i,j,di,dj,Bex)
      
        for k,t in enumerate(tex[0:-1]):  
        #for k in xrange(snapshots-1):
        
          mask_int[k] = tl.get_mask_value(H_int[k],B_int[k],sl[k,0,0])
          mask_int[k+1] = tl.get_mask_value(H_int[k+1],B_int[k+1],sl[k,0,0])
          if mask_int[k] == cf.mgr and mask_int[k+1] > mask_int[k]:
            retreat = 0.5 * (tex[k+1] + t)

          ### interpolate to 100yr resolution!
          tsteps = np.arange(t,tex[k+1],100.0)
          mi = np.zeros_like(tsteps)
          mi[0] = mask_int[k]
          for l,tint in enumerate(tsteps):
              t_before = tint - t
              t_after  = tex[k+1] - tint 
              Hi = (H_int[k]*t_after + H_int[k+1]*t_before) / (tex[k+1]-t)
              Bi = (B_int[k]*t_after + B_int[k+1]*t_before) / (tex[k+1]-t)
              sli = (sl[k,0,0]*t_after + sl[k+1,0,0]*t_before) / (tex[k+1]-t)
              mi[l] = tl.get_mask_value(Hi,Bi,sli)
              #print tint,mi[l],t_before,t_after,Hi,Bi
              if l > 0:
                if mi[l-1] == cf.mgr and mi[l] > mi[l-1]:
                  retreat_fine = 0.5 * (tint + tsteps[l-1])
                  #print tint,mi[l]

        retreat = max(retreat,retreat0) #>-35ka
        retreat_fine = max(retreat_fine,retreat0) #>-35ka
        retreat = retreat_fine

        # misfit = model - obs
        if retreat == retreat0:
          misfit = (retreat0)**2
          count_ooc+=1
        else:
          misfit = (retreat - dat_age)**2.0

        if dat[ids]['type']=="G": #two-way GLR
          sigma += (s_struct_glr)**2
        elif dat[ids]['type']=="M": #one-way OMC
          sigma += (s_struct_omc)**2

          modelled_event_before_obs = (dat_age > retreat and retreat>retreat0)
          only_open_conditions = (retreat == retreat0 and mask_int[-1] > cf.mgr)
          if (modelled_event_before_obs or only_open_conditions): #one-sided gaussian!
            sigma += (s_struct_omc)**2
            #print "oneway",dat_age,modelled_event_before_obs,retreat,only_open_conditions,count_ooc
        
        misfit/=sigma
        mean_misfit+=misfit

        """
        ### check in plots
        if self.showplots:
          fig17 = plt.figure(17,figsize=(8, 4))
          ax2=plt.subplot(111)
          #ax2.plot(tex[:],mex[:,i,j],"g.",markersize=4) # i,j mask value
          ax2.plot(tex[:],mask_int[:],"m.",markersize=2) # interpolated mask value
          ax2.axvline(dat_age,color="b") # event in paleo records
          ax2.axvline(dat_age + dat[ids]['ageu'],color="b",linestyle="dashed")
          ax2.axvline(dat_age - dat[ids]['ageu'],color="b",linestyle="dashed")
          ax2.axvline(retreat,color="r") # event in model
          #ax2.axvline(retreat_fine,color="g") # fine event in model

          if dat[ids]['type']=="M": # model uncertainty
            if modelled_event_before_obs:
              ax2.axvline(retreat+np.sqrt(sigma),color="r",linestyle="dashed")
            else:
              ax2.axvline(retreat-np.sqrt(sigma),color="r",linestyle="dashed")
          else:
            ax2.axvline(retreat-np.sqrt(sigma),color="r",linestyle="dashed")
            ax2.axvline(retreat+np.sqrt(sigma),color="r",linestyle="dashed")
          
          ax2.axis([-35000,0,1.5,4.5])
          plt.text(-30000, 3.5, '%d: score: %.3f, %d,%d' % (ids,np.sqrt(misfit),j,i))
          plt.show()
      """

      mean_misfit = mean_misfit / len(dat)

      if self.printout:
          print '\nEXT (mean of gaussian age error in grounding line transition datum)'
          print mean_misfit

      return (mean_misfit)


  ################################################################
  def get_rslcalc(self,dat,tex,vars):
      
      Bex=vars[0]
      sl=vars[1] # in 0,0 corner
      tex=tex[0]/cf.seconds_per_year
      snaps = np.shape(Bex)[0]
      #print dat

      mean_misfit=0.0
      sigma_low=2.0
      sigma_up=50.0

      # see Table 2 in Briggs et al., 2013
      #rsl_type={}
      #rsl_type['1']="Two-way"
      #rsl_type['2a']="One-way" #RSL minimum (RSL above)
      #rsl_type['2b']="One-way" #RSL maximum (RSL was below)
      #rsl_type['3']="Two-way adjusted" 
      #rsl_type['4a']="Two-way adjusted" #RSL minimum (RSL above)
      #rsl_type['4b']="Two-way adjusted" #RSL maximum (RSL was below)

      for ids in dat.keys():
        #print ids,dat[ids]

        i            =  dat[ids]['i_grid']
        j            =  dat[ids]['j_grid']
        di           =  dat[ids]['di_grid']
        dj           =  dat[ids]['dj_grid']
        sl_dat       = dat[ids]['val']
        sl_sigma     = dat[ids]['valu']
        sigma_age    = 2.0*dat[ids]['ageu']
        dat_age      = -dat[ids]['age']
        lower_tbound = dat_age - sigma_age
        upper_tbound = dat_age + sigma_age
        dat_weight   = dat[ids]['weight']
        dtype        = dat[ids]['type']

        ### interpolation to lon lat location
        B_int = tl.get_interpolated_value(i,j,di,dj,Bex)
        rsl_mod = (sl[:,0,0]-B_int[:])-(sl[-1,0,0]-B_int[-1])

        ### find minimum model-obs difference within age uncertainty (confidence interval)
        ltb = np.argmin((tex-lower_tbound)**2.0) #index of closest snapshot
        utb = np.argmin((tex-upper_tbound)**2.0)+1
        drsl = np.min(np.abs(rsl_mod[ltb:utb]-sl_dat))
        min_snap = ltb+np.argmin(np.abs(rsl_mod[ltb:utb]-sl_dat))
        mint = tex[min_snap]
        sign_dsl = np.sign(rsl_mod[min_snap]-sl_dat)
        #print rsl_mod[ltb:utb],sl_dat,drsl,np.min(np.abs(rsl_mod[ltb:utb]-sl_dat)),drsl,mint

        ### one and two sided error model
        sigma_dat = max(sl_sigma+1.0,sigma_low)

        sigma_sided ={}
        sigma_sided[1] = sigma_sided[-1] = sigma_dat
        if dtype in ["2a","4a"]: #"minimum"
          sigma_sided[1] = max(sigma_up,sigma_dat)
        elif dtype in ["2b","4b"]: #"maximum"
          sigma_sided[-1] = max(sigma_up,sigma_dat)

        #print dat[ids]['type'],sl_sigma,sl_dat,dat_age,sign_dsl,sigma_sided[sign_dsl]

        #approach a
        #if self.appr =="a":
        misfit = (drsl/sigma_sided[sign_dsl])**2.0
        mean_misfit += misfit*dat_weight

        """
        if self.showplots:
          fig18 = plt.figure(18,figsize=(8, 4));plt.clf()
          ax4=plt.subplot(111)
          ax4.axhline(sl_dat,color="k")
          ax4.axvline(dat_age,color="k")
          ax4.axvline(mint,color="b")
          ax4.plot(tex[:],rsl_mod[:],"r-")

          ax4.axvline(lower_tbound,color="k",linestyle="dashed")
          ax4.axvline(upper_tbound,color="k",linestyle="dashed")
          #ax4.plot(tex[:],rsl_mod[:]+sigma_sided[1],"r",linestyle="dotted")
          #ax4.plot(tex[:],rsl_mod[:]-sigma_sided[-1],"r",linestyle="dotted")
          ax4.axhline(sl_dat+sigma_sided[1],color="k",linestyle="dotted")
          ax4.axhline(sl_dat-sigma_sided[-1],color="k",linestyle="dotted")

          ax4.axis([-35000,0,-100,200])
          plt.text(-30000, 100, '%d: score: %.3f, %d,%d' % (ids,np.sqrt(misfit),j,i))
          plt.show()
          """

      if self.printout:
          print '\nRSL (mean square error of regional sea-level change at datum)'
          print mean_misfit

      return (mean_misfit)


  #####################################################################
  def print_to_txtfile(self,sc):

        if not os.path.exists(self.outpath+"stats"): 
          os.makedirs(self.outpath+"stats")

        statfile = self.outpath+"stats/le_ens"+self.ensnum+"_"+str(int(self.res))+"km.txt"
        if self.printout:
          #print "\n...Write paleo score to "+statfile
          print "\n###############################################"
        #if os.path.isfile(statfile):
        #  os.system("rm "+statfile)
        savestat = open(statfile, 'a')
        savetext = self.score_text
        savestat.write(savetext)
        savetext = "\n"+self.ensnum+":"
        for sco in sc:
          savetext+=" "+str(sco)
        savetext+="\n"
        savestat.write(savetext)
        savestat.close()


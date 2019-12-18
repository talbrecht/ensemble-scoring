#!/usr/bin/env python

###############################################################################
# Copyright (C) 2017-2019 Potsdam-Institute for Climate Impact Reasearch (PIK), 
# Author: Torsten Albrecht (albrecht@pik-potsdam.de)
# License: GNU AFFERO GENERAL PUBLIC LICENSE version 3
# 
# This script executes the ensemble analysis for present day data.
###############################################################################



import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors

#on pik cluster
#source activate python_for_pism_calib

## this hack is needed to import config.py from the project root
#project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#if project_root not in sys.path: sys.path.append(project_root)
import config as cf; reload(cf)
import tools as tl; reload(tl)


class pd_score(object):

  def __init__(self,enum,obs,snp):

    self.ensnum = str(enum)
    self.showplots = snp[0]
    self.printtopdf = snp[1]
    self.printout = snp[2]

    ### settings

    self.res=cf.resolution #km resolution

    self.appr = cf.approach

    self.score_names = cf.score_names
    self.pd_score_names =["DSLV","TOTE","TOTI","TOTDH","TOTVEL","TOTGL"] #available measures
    self.score_text=""

    # do iterative ice rise exclude method for mean GL distance
    ex_ir=True

    #pathnames
    #workpath     = cf.workpath
    resultpath   = cf.resultpath
    self.outpath = cf.output_data_path

    pism_file_name = cf.pism_file_name.replace(str(cf.fillnum),str(enum))
    pismfile = os.path.join(resultpath,pism_file_name)

    global Mx,My,mobs,Bobs,Hobs,cellarea,velobs,velstnd #,hobs #only the one used in the functions below
    [x,y,Mx,My,mobs,Bobs,Hobs,hobs,lonobs,latobs,cellarea,velobs,velstnd] = obs

    if os.path.exists(pismfile):

      thk,topg,usurf,mask,velsurf,okmask = tl.get_data(pismfile,['thk','topg','usurf','mask','velsurf_mag','ocean_kill_mask'])

      ### total ice volume and sea-level relevant volume
      totvolcalc = self.get_totvolcalc([thk])
      slvolcalc = self.get_slvolcalc([thk,topg,usurf])

      ### mismatch in grounded and floating areas
      toteicalc = self.get_toteicalc([mask,okmask])

      ### rms error in ice thickness
      dhcalc = self.get_dhcalc([mask,thk,okmask])
      #dhcalc = self.get_dhcalc([mask,usurf,okmask])
  
      ### grounding line distance (around whole Antarctica)
      mgldcalc = self.get_gl_dist([mask],ex_ir)
      #mgldcalc=1.0 #if skfmm not available

      ### grounding line distance (in Ross along transect), Maris et al. 2014
      #TODO?!

      ### rms error in surface velocity
      velcalc = self.get_velcalc([mask,velsurf])


      score_choice={} #indicated which measures to consider
      score_choice["DSLV"]=slvolcalc[1]
      score_choice["TOTE"]=toteicalc[0]
      score_choice["TOTI"]=toteicalc[3]
      score_choice["TOTDH"]=dhcalc[1]
      score_choice["TOTVEL"]=velcalc[1]
      score_choice["TOTGL"]=mgldcalc


      ### write scores to txt file 
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


  ##############################################
  # calculate total volume anomaly
  def get_totvolcalc(self,var):

    H=var[0] #ice thickness variable

    volobs = np.sum(np.sum(cellarea * Hobs)) * 1e-9 #mio km3
    volens = np.sum(np.sum(cellarea * H)) * 1e-9 #mio km3
    dvol   = np.sum(np.sum(cellarea * (H-Hobs))) * 1e-9 #mio km3

    if self.printout:
      print '\nTOTV (anomaly in total ice volume in mio. km3)'
      #print volens,dvol
      print dvol

    return (volens,dvol)

  ##############################################
  # calculate sea-level equivalent volume anomaly
  def get_slvolcalc(self,var):

    H=var[0] #ice thickness variable
    B=var[1] #bed topography variable
    sl=var[2][0,0] #sea-level

    volobs = 0.
    volens = 0.

    dvol = 0.
    rhori = cf.rhosw/cf.rhoi

    for i in range(Mx):
      for j in range(My):

        #grounded in model
        if (H[i,j] > (sl-B[i,j]) * rhori and H[i,j] > 10.0): #as in iMreport in PISM
          volgr = H[i,j] * cellarea[i,j]
          volbf = (sl-B[i,j]) * cellarea[i,j]
          if sl < B[i,j]:
            volens += volgr
          else:
            volens += (volgr - rhori * volbf) 
      
        #grounded in observations
        if (Hobs[i,j] > (0.0-Bobs[i,j]) * rhori and Hobs[i,j] > 10.0): #as in iMreport in PISM
          volgrobs = Hobs[i,j] * cellarea[i,j]
          volbfobs = (sl-Bobs[i,j]) * cellarea[i,j]
          if B[i,j] > 0.0:
            volobs += volgrobs
          else:
            volobs += (volgrobs - rhori * volbfobs)
 
    volens *= 1e-3 / cf.km3_to_msle #in m sle
    volobs *= 1e-3 / cf.km3_to_msle #in m sle

    dvol = volens - volobs

    if self.printout:
      print '\nSLVOL (anomaly in sea-level equivalent ice volume in m SLE)'
      #print volens,dvol
      print dvol

    return (volens,np.abs(dvol))


  ##############################################
  # calculate grounded/floating area misfit
  def get_toteicalc(self,var):

    #Calculates area mismatch in (i) for grounded ice sheet
    #and (ii) for floating ice shelves

    m=var[0] #mask variable
    okm=var[1] #ocean_kill_mask variable

    mf=np.ones_like(m)*np.NaN
    mg=np.ones_like(m)*np.NaN

    mgr = cf.mgr
    mfl = cf.mfl

    tota = np.sum(np.sum(cellarea)) #FIXME: in Pollard nested grid of WAIS
    #tota = np.sum(np.sum(cellarea[ okm==0.0 ]))
    #tota_gr = np.sum(np.sum(cellarea[ mobs == mgr ]))
    #tota_fl = np.sum(np.sum(cellarea[ mobs == mfl ]))

    #FIXME: add sectors
        #if (alond(i,j) .gt.-120..and. alond(i,j).lt. -90.) then
          #ASE (PIG/THW):
        #else if (alond(i,j) .gt. 150..or.  alond(i,j).lt.-120.) then
          #Ross:
        #else if (alond(i,j) .gt.-90. .and. alond(i,j).lt.  0.) then
          #Weddell:

    #grounded
    cond_gr_p = (mobs != mgr) & (m == mgr) #positive misfit (grounded where obs is not grounded)
    cond_gr_m = (mobs == mgr) & (m != mgr) #negative misfit (not grounded where obs is grounded)
    cond_gr_0 = (mobs == mgr) & (m == mgr) #for the plot
    zdareag = np.sum(np.sum(cellarea[ cond_gr_p | cond_gr_m ])) #total misfit
    zdareagp = np.sum(np.sum(cellarea[ cond_gr_p ]))
    zdareagm = np.sum(np.sum(cellarea[ cond_gr_m ]))

    #floating
    cond_fl_p = (mobs != mfl) & (m == mfl)
    cond_fl_m = (mobs == mfl) & (m != mfl)
    cond_fl_0 = (mobs == mfl) & (m == mfl)
    zdareaf = np.sum(np.sum(cellarea[ cond_fl_p | cond_fl_m ]))
    zdareafp = np.sum(np.sum(cellarea[ cond_fl_p ]))
    zdareafm = np.sum(np.sum(cellarea[ cond_fl_m ]))

    mg[ cond_gr_p ] = 1.0
    mg[ cond_gr_m ] = -1.0
    mg[ cond_gr_0 ] = 0.0

    mf[ cond_fl_p ] = 1.0
    mf[ cond_fl_m ] = -1.0
    mf[ cond_fl_0 ] = 0.0


    #approach A
    if self.appr =="a":
      sigwid=30.0
      totb=np.sqrt(tota)*sigwid
      #totbgr=np.sqrt(tota_gr)*sigwid
      #totbfl=np.sqrt(tota_fl)*sigwid

      zdareag=(zdareag/totb)**2
      zdareagp=(zdareagp/totb)**2
      zdareagm=(zdareagm/totb)**2

      zdareaf=(zdareaf/totb)**2
      zdareafp=(zdareafp/totb)**2
      zdareafm=(zdareafm/totb)**2

      #zdareaf=(zdareaf/totbfl)**2
      #zdareafp=(zdareafp/totbfl)**2
      #dareafm=(zdareafm/totbfl)**2


      print_measure='mismatch'

    #approach B
    elif self.appr =="b":

      zdareag=zdareag/tota
      zdareagp=zdareagp/tota
      zdareagm=zdareagm/tota

      zdareaf=zdareaf/tota
      zdareafp=zdareafp/tota
      zdareafm=zdareafm/tota

      print_measure='rmse'

    if self.showplots or self.printtopdf:  #FIXME: use basemap

      ###map of PD surface elevation compared to Bedmap2 #######
      fig10, [ax10,ax10b] = plt.subplots(1, 2,figsize=(11, 7))      
      #plt.clf()
      plt.title("mismatch of modeled and observed mask")
      #diffm = (m-mobs)
      #diffmm=np.ma.array(diffm,mask = ( (mobs == moc) & (m == moc) ))
      #'+1' mean grounded became floating or floating became ice free ocean and -'1' vice versa

      #ax10=plt.subplot(121, aspect='equal')
      ax10.axis("equal")
      ax10.axis([0,Mx,0,My])
      ax10.axis("off")

      relm=1
      tcks=np.arange(-relm-1,relm+1,1)
      cs1  = ax10.contourf(mg,tcks,alpha=0.9,cmap=cm.RdYlBu_r) #,extend="both")
      ## colorbar diffm
      #cbaxes = fig10.add_axes([0.15, 0.25, -0.04, 0.5])
      #cb = plt.colorbar(cs1,ticks=tcks,orientation='horizontal')
      cb = plt.colorbar(cs1,orientation='horizontal',shrink=0.6, aspect=20)


      cb.set_ticks(tcks+.5)
      cb.set_ticklabels(tcks+1)
      cb.set_label('mismatch grounded mask',multialignment="left")
      cb.outline.set_linewidth(0)
      #bedmap2 contour
      cs7 = ax10.contour(mobs,[2.5,3.5],colors="w",linewidths=.5,alpha=1.,)
      #cs7b = ax10.contour(mobs,[1.5,2.5],colors="k",linewidths=.5,alpha=1.,)

      #ax10b=plt.subplot(122, aspect='equal')
      ax10b.axis("equal")
      ax10b.axis([0,Mx,0,My])
      ax10b.axis("off")
      cs1b  = ax10b.contourf(mf,tcks,alpha=0.9,cmap=cm.RdYlBu_r) #,extend="both")
      cb2 = plt.colorbar(cs1b,orientation='horizontal',shrink=0.6, aspect=20)
      cb2.set_ticks(tcks+.5)
      cb2.set_ticklabels(tcks+1)
      cb2.set_label('mismatch floating mask',multialignment="left")
      cb2.outline.set_linewidth(0)
      #bedmap2 contour
      cs7b = ax10b.contour(mobs,[2.5,3.5],colors="w",linewidths=.5,alpha=1.,)

      #ax10b.text(Mx/2,10,"all: "+str(np.around(zdareaf,decimals=1))+"\npos: "+str(np.around(zdareafp,decimals=1))+"\nneg: "+str(np.around(zdareafm,decimals=1)))

      plt.subplots_adjust(wspace = -0.2, hspace = 0.0)

      ### print plot to pdf file
      printname = self.outpath+"plots/mismatch_ensemble"+self.ensnum+"_"+str(int(self.res))+"km.pdf"
      if self.printtopdf:
        plt.savefig(printname, format='pdf')

    if self.printout:
      print '\nTOTE ('+print_measure+' grounded area)'
      #print zdareag,zdareagp,zdareagm,totb,np.sqrt(zdareag)*totb,np.sqrt(zdareag)
      print zdareag
      print '\nTOTI ('+print_measure+' floating area)'
      #print zdareaf,zdareafp,zdareafm,totb,np.sqrt(zdareaf)*totb,np.sqrt(zdareaf)
      #print self.appr,np.sum(np.sum(cellarea[ cond_fl_p | cond_fl_m ])),np.sum(np.sum(cellarea[ cond_fl_p ])),np.sum(np.sum(cellarea[ cond_fl_m ])),tota,totb
      print zdareaf

    return (zdareag,zdareagp,zdareagm,zdareaf,zdareafp,zdareafm)


  ################################################################
  def get_dhcalc(self,var):

    #Calculates rms error in (i) modern grounded ice thicknesses
    #and (ii) modern bedrock elevations over all domain.

    #Pollrad: WAIS only (as in BT13 Fig. 1):
    #if ( (alond(i,j) .gt. 170. .or. alond(i,j).lt. -30.) .and.(alatd(i,j) .gt. -86.) ) then

    m = var[0] #mask variable
    h = var[1] #ice thickness
    okm=var[2] #ocean_kill_mask variable
    ho = Hobs
    #ho = hobs

    diffh = (h-ho)


    if self.showplots or self.printtopdf:  #FIXME: use basemap

      ###map of PD surface elevation compared to Bedmap2 #######
      fig11 = plt.figure(11,figsize=(9, 9))
      plt.clf()
      ax11=plt.subplot(111)
      ax11.axis("equal")
      ax11.axis([0,Mx,0,My])
      ax11.axis("off")
      
      diffhm=np.ma.array(diffh,mask = ( (mobs == cf.moc) & (m == cf.moc) ))
      topg = np.ma.array(Bobs, mask = (mobs < cf.moc))
      topg.mask[topg < -3000] = True
      relh=1000
      ticks11=np.arange(-relh,relh+250,250)
      cs2  = ax11.contourf(topg,[-4750,-4250,-3750,-3250,-2750,-2250,-1750,-1250,-750,-250],alpha=0.3,cmap=cm.Greys_r)
      cs1  = ax11.contourf(diffhm,ticks11,alpha=1.0,cmap=cm.RdBu_r,extend='both')
      cs7a = ax11.contour(m,[2.5,3.5],colors="k",linewidths=.5,alpha=1.,)
      ## colorbar diffh
      cbaxes11 = fig11.add_axes([0.15, 0.23, 0.3, 0.02])
      #cbaxes = fig11.add_axes([0.15, 0.25, -0.04, 0.5])
      #cb = plt.colorbar(cs1,cax=cbaxes,ticks=np.arange(-relh,relh+100,100))
      ticks11b=np.arange(-relh,relh+500,500)
      cb = plt.colorbar(cs1,cax=cbaxes11, orientation="horizontal",ticks=ticks11b)
      #cb.set_label('surface elevation anomaly [m]',multialignment="left")
      cb.set_label('ice thickness anomaly [m]',multialignment="left")
      cb.outline.set_linewidth(0)

      #bedmap2 contour
      cs7 = ax11.contour(mobs,[2.5],colors="w",linewidths=1.0,alpha=1.,)

      ### print plot to pdf file
      #printname = self.outpath+"plots/diffsurf_ensemble"+self.ensnum+"_"+str(int(self.res))+"km.pdf"
      printname = self.outpath+"plots/diffthk_ensemble"+self.ensnum+"_"+str(int(self.res))+"km.pdf"
      if self.printtopdf:
        plt.savefig(printname, format='pdf')
        #plt.savefig(printname.replace(".pdf",".png"), format='png',dpi=300)

    #all grid points (incl. open ocean)
    zdhb_a = np.sum(np.sum((cellarea*(h-ho)**2)))
    #zahb = np.sum(np.sum(cellarea))
    zahb = np.sum(np.sum(cellarea[ okm==0.0 ]))

    #obs grounded only
    dhcm = (cellarea*(h-ho)**2)[mobs==cf.mgr]
    zdh_a = np.sum(np.sum(dhcm))
    zah = np.sum(np.sum(cellarea[mobs==cf.mgr]))

    #FIXME: wais only

    #approach A
    if self.appr =="a":
      sigh = 10.0
      totdh  = zdh_a  /zah  /(sigh**2)
      totdhb = zdhb_a /zahb /(sigh**2)
      printline='\nTOTDH (misfit ice thickness)'

    #approach B
    elif self.appr =="b":
      totdh  = np.sqrt(zdh_a  / zah )
      totdhb = np.sqrt(zdhb_a / zahb)
      printline='\nTOTDH (rmse ice thickness in m)'

    if self.printout:
      #print zdh_a,zah,zdhb_a,zahb
      print printline
      #print totdh,totdhb
      print totdhb

    return (totdh,totdhb)


  ################################################################
  def get_gl_dist(self,m,exir):

    import scipy.ndimage
    import skfmm

    if exir: #excluding ice rises

      mni = 1

      glmask = np.zeros_like(m[0]) 
      glmask[m[0] == cf.mif] = 2 #remove ice free bedrock, FIXME: floating points in inner ice shield
      glmwir = tl.mask_without_icerises(m[0],(Mx/2,My/2),mni)
      glmask[glmwir == mni] = 1
      glmask[glmwir != mni] = -1

      glmaskobs = np.zeros_like(mobs)
      mobs_copy = np.copy(mobs)
      #glmaskobs[mobs == mif] = 2 #remove ice free bedrock, FIXME: floating points in inner ice shield
      mobs_copy[Mx/2,My/2] = cf.mgr
      glmobs = tl.mask_without_icerises(mobs_copy,(Mx/2,My/2),mni)
      glmaskobs[glmobs == mni] = 1
      glmaskobs[glmobs != mni] = -1


    else: #not excluding ice rises

      glmask=m[0]    
      glmask[glmask <= cf.mgr] = -1
      glmask[glmask > cf.mgr] = 1

      glmaskobs=mobs
      glmaskobs[glmaskobs <= cf.mgr] = -1
      glmaskobs[glmaskobs > cf.mgr] = 1

    distanceobs = skfmm.distance(glmaskobs)*self.res #km

    #if showplots: 
    if True: #needs to be always true in order to calculate GL countour line
      fig14 = plt.figure(14,figsize=(6, 7))
      plt.clf()

      ax14a=plt.subplot(111)
      isa = ax14a.imshow(distanceobs,cmap=cm.RdBu_r,vmin=-1500,vmax=1500)
      cba = plt.colorbar(isa,ticks=np.arange(-1500,2000,500),orientation='horizontal')
      cba.set_label('sign distance from observed GL',multialignment="left")
      cba.outline.set_linewidth(0)
      ax14a.contour(distanceobs,0,colors='k',linewidth=1)
      cs = ax14a.contour(glmask,[0.0],colors='r',linewidth=1)
      plt.ylim(plt.ylim()[::-1])
      ax14a.axis("equal")
      ax14a.axis([0,Mx,0,My])
      ax14a.axis("off")

      ### print plot to pdf file
      #printname = self.outpath+"plots/signdistgl_ensemble"+self.ensnum+"_"+str(int(self.res))+"km.pdf"
      #if self.printtopdf:
      #  plt.savefig(printname, format='pdf')

    #calculate mean along the grounding line(s)
    mean_dist_gl = 0
    cnt_p = 0
    for p in cs.collections[0].get_paths()[:]:
      v = p.vertices
      cx = v[:,0]
      cy = v[:,1]
      lenc=len(cx)
      cnt_p += lenc
      diffdistint = scipy.ndimage.map_coordinates(distanceobs, [cy, cx], order=1)
      for i in xrange(lenc):
        #print i,cx[i],cy[i],diffdistint[i]
        #if i>0:
        #  print np.sqrt((cx[i]-cx[i-1])**2+(cy[i]-cy[i-1])**2)
        mean_dist_gl += (diffdistint[i])**2
        #mean_dist_gl += np.abs(diffdistint[i])
    #mean_dist_gl=np.sqrt(mean_dist_gl/cnt_p)
    mean_dist_gl=mean_dist_gl/cnt_p


    ax14a.text(Mx*2.5/5.0,My/20.0,"RS mean GL dist: "+str(np.around(np.sqrt(mean_dist_gl),decimals=1))+" km",fontsize=10)
    ### print plot to pdf file
    printname = self.outpath+"plots/signdistgl_ensemble"+self.ensnum+"_"+str(int(self.res))+"km.pdf"
    if self.printtopdf:
        plt.savefig(printname, format='pdf')

    if self.appr =="a":
      sigwid=30.0
      mean_dist_gl/=(sigwid**2)
      printline='\nTOTGL (misfit GL distance)'

    #approach B
    elif self.appr =="b":
      mean_dist_gl=np.sqrt(mean_dist_gl)
      printline='\nTOTGL (rmse GL distance in km)'

    if self.printout:
      print printline
      #print mean_dist_gl,cnt_p
      print mean_dist_gl
    
    return (mean_dist_gl)


  ################################################################
  def get_velcalc(self,var):

      m=var[0] #mask variable
      v=var[1] #surface speed variable

      mfl=cf.mfl
      mgr=cf.mgr

      #all grounded and floating points in obs and mod)
      cond_all = ((mobs <= mfl) & (m <= mfl) & (velobs > 0.0))
      #cond_all = ((mobs <= mfl) & (m <= mfl) & (velobs >= 100.0))
      if self.appr =="a":
        dvcm = (cellarea*((v-velobs)/velstnd)**2)[ cond_all ]
        #print np.mean(velstnd),np.max(velstnd),np.min(velstnd)

        #for i in xrange(Mx):
        #  for j in xrange(My):
        #    if (mobs[i,j]<=mfl and m[i,j]<=mfl and velobs[i,j]>0.0):
        #      xt=i
        #      yt=j
        #      mse=cellarea[xt,yt]*((v[xt,yt]-velobs[xt,yt])/velstnd[xt,yt])**2
        #      if mse>100000:
        #        print i,j,cellarea[xt,yt],v[xt,yt],velobs[xt,yt],velstnd[xt,yt],mse
        #print cellarea[xt,yt],v[xt,yt],velobs[xt,yt],velstnd[xt,yt],cellarea[xt,yt]*((v[xt,yt]-velobs[xt,yt])/velstnd[xt,yt])**2
        #print len(dvcm)
      elif self.appr =="b":
        dvcm = (cellarea*(v-velobs)**2)[ cond_all ]

      zdvb_a = np.sum(np.sum(dvcm))
      zavb = np.sum(np.sum(cellarea[ cond_all ]))
      #print zdvb_a,zavb
      #print np.shape(dvcm),np.shape(cellarea[ cond_all ])

      #all grounded or floating points in mod and grounded in obs)
      cond_gr = ((mobs == mgr) & (m <= mfl))
      #cond_gr = ((mobs == mgr) & (m <= mfl) & (velobs >=100.0))
      if self.appr =="a":
        dvcmgr = (cellarea*((v-velobs)/velstnd)**2)[ cond_gr ]
      elif self.appr =="b":
        dvcmgr = (cellarea*(v-velobs)**2)[ cond_gr ]
      zdv_a = np.sum(np.sum(dvcmgr))
      zav = np.sum(np.sum(cellarea[ cond_gr ]))

      #approach a
      if self.appr =="a":
        #sigwid = 10.0
        totdv  = (zdv_a  / zav ) #/ sigwid**2
        totdvb = (zdvb_a / zavb) #/ sigwid**2
        printline='\nTOTVEL (misfit velocity magnitude)'

      #approach B
      elif self.appr =="b":
        totdv  = np.sqrt(zdv_a  / zav )
        totdvb = np.sqrt(zdvb_a / zavb)
        printline='\nTOTVEL (rmse velocity magnitude in m/yr)'

      if self.showplots or self.printtopdf: 

        csurf_anom = np.ma.array(v - velobs)
        csurf_anom.mask = ((m>mfl) | (mobs>mfl))

        vel_fl = np.ma.array(np.copy(v))
        velobs_fl = np.ma.array(np.copy(velobs))
        vel_fl.mask = ((m!=mfl) | (mobs!=mfl))
        velobs_fl.mask = ((m!=mfl) | (mobs!=mfl))

        fig15 = plt.figure(15,figsize=(8, 8))
        plt.clf()
        ax15=plt.subplot(111)

        maxvel=3000.0
        s = 3 #skip in scatter plot
        #ax15.plot(v,velobs,"k.",alpha=0.2)
        ax15.loglog(velobs[::s,::s],v[::s,::s],"k.", basex=10,alpha=0.2)
        ax15.loglog(velobs_fl[::s,::s],vel_fl[::s,::s],"g.", basex=10,alpha=0.2)
        ax15.plot([0.1,maxvel],[0.1,maxvel],"r-")
        ax15.set_xlabel("Rignot surface velocity magnitude (m/yr)",fontsize=14)
        ax15.set_ylabel("PISM surface velocity magnitude (m/yr)",fontsize=14)
        ax15.text(1e2,5e-1,"RMSE: "+str(np.around(totdvb,decimals=1))+" m/yr",fontsize=14)
        #ax15.axis("equal")
        #ax15.axis("off")
        ax15.axis([0.1,maxvel,0.1,maxvel])

        ### print plot to pdf file
        printname = self.outpath+"plots/velscatter_ensemble"+self.ensnum+"_"+str(int(self.res))+"km.pdf"
        if self.printtopdf:
          #plt.savefig(printname, format='pdf')
          plt.savefig(printname.replace(".pdf",".png"), format='png',dpi=300)

      if self.printout:
        print printline
        #print totdv,totdvb
        print totdvb

      return (totdv,totdvb)


  #####################################################################
  def print_to_txtfile(self,sc):

        if not os.path.exists(self.outpath+"stats"): 
          os.makedirs(self.outpath+"stats")

        statfile = self.outpath+"stats/le_ens"+self.ensnum+"_"+str(int(self.res))+"km.txt"
        if self.printout:
          print "\n...Write PD score to "+statfile
          #print "###############################################"
        if os.path.isfile(statfile):
          os.system("rm "+statfile)
        savestat = open(statfile, 'a')
        savetext = self.score_text
        savestat.write(savetext)
        savetext = "\n"+self.ensnum+":"
        for sco in sc:
          savetext+=" "+str(sco)
        savetext+="\n"
        savestat.write(savetext)
        savestat.close()

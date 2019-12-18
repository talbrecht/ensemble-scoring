#!/usr/bin/env python

###############################################################################
# Copyright (C) 2017-2019 Potsdam-Institute for Climate Impact Reasearch (PIK), 
# Author: Torsten Albrecht (albrecht@pik-potsdam.de)
# License: GNU AFFERO GENERAL PUBLIC LICENSE version 3
# 
# This script provides tools for the ensemble analysis
###############################################################################



from netCDF4 import Dataset as NC
import numpy as np

import config as cf; reload(cf)

#import matplotlib.pyplot as plt
#from matplotlib import cm, colors
#import matplotlib as mpl

################################################################

def ncload(filename):
  
  try:
    data = NC(filename, 'r')
  except Exception:
    print """Specify NetCDF input file."""
    print filename
    exit(-1)
  return data


def get_data(fn,variables):

  vardata = []
  fndata = ncload(fn)
  for varname in variables:
    var = np.squeeze(fndata.variables[varname][:])
    vardata.append(var)

  return vardata

######### load observational data  ###############################
def load_observations(of,ovf):

  ### Bedmap2
  bedmap = ncload(of)

  #global xobs,yobs,Mx,My,Hobs,hobs,Bobs,mobs,cellarea,latobs,lonobs

  xobs = np.squeeze(bedmap.variables["x"][:])
  yobs = np.squeeze(bedmap.variables["y"][:])
  latobs = np.squeeze(bedmap.variables["lat"][:])
  lonobs = np.squeeze(bedmap.variables["lon"][:])
  Hobs = np.squeeze(bedmap.variables["thk"][:])
  Bobs = np.squeeze(bedmap.variables["topg"][:])
  hobs = np.ma.array(np.squeeze(bedmap.variables["usurf"][:]))
  mobs = np.squeeze(bedmap.variables["mask"][:])
  cellarea = np.squeeze(bedmap.variables["cell_area"][:])
  
  #obsres = np.diff(xobs)[0]/1000.0 #(x[1]-x[0])/1000.0
  Mx=np.shape(mobs)[0]
  My=np.shape(mobs)[1]
  bedmap.close()


  ### Rignot velocities
  veldata = ncload(ovf)
  #global velobs
  velobs = np.ma.array(np.squeeze(veldata.variables["v_magnitude"][:]))
  velstnd = np.ma.array(np.squeeze(veldata.variables["v_std"][:]))
  velstnd[velstnd<1.5]=1.5 #lower bound
  veldata.close()


  return [xobs,yobs,Mx,My,mobs,Bobs,Hobs,hobs,lonobs,latobs,cellarea,velobs,velstnd]

################################################################
def get_xyz_coordinates(lonval,latval):

  xc = np.cos (latval*np.pi/180.)*np.cos (lonval*np.pi/180.)
  yc = np.cos (latval*np.pi/180.)*np.sin (lonval*np.pi/180.)
  zc = np.sin (latval*np.pi/180.)

  return (xc,yc,zc)



################################################################
def get_closest_point(lonval,latval,lonobs,latobs):

      zdistmin = 1.e20
      iupl=0
      jupl=0
      zxs,zys,zzs = get_xyz_coordinates(lonval,latval)
      Mx,My = np.shape(lonobs)

      #find the closest point on grid (not interpolated!)
      for i in range(Mx):
        for j in range(My):
          zx,zy,zz = get_xyz_coordinates(lonobs[i,j],latobs[i,j])
          zdist = np.sqrt ((zxs-zx)**2 + (zys-zy)**2 + (zzs-zz)**2)
          if (zdist<zdistmin):
            iupl = i
            jupl = j
            zdistmin = zdist

      return (iupl,jupl)

#### calculate interpolated value
def get_interpolated_value(i,j,di,dj,val):

  i2 = np.int(i+np.sign(di))
  j2 = np.int(j+np.sign(dj))
  di=np.abs(di)
  dj=np.abs(dj)
  #print np.shape(val),np.ndim(val)
  
  if np.ndim(val)==2:
    return (1.0-dj)*((1.0-di)*val[i,j] + di*val[i,j2]) + dj*((1.0-di)*val[i2,j] + di*val[i2,j2])
  elif np.ndim(val)==3:
    return (1.0-dj)*((1.0-di)*val[:,i,j] + di*val[:,i,j2]) + dj*((1.0-di)*val[:,i2,j] + di*val[:,i2,j2])
  else:
    print "No interpolation possible for dimensions of order "+str(np.dim(val))
    return 0


################################################################
def find_interpolation(icl,jcl,lon0,lat0,lonobs,latobs):


  lonij = lonobs[icl,jcl]
  latij = latobs[icl,jcl]

  xi,yi,zi = get_xyz_coordinates(lonij,latij)
  x0,y0,z0 = get_xyz_coordinates(lon0,lat0)

  #direct grid neighbors
  xn,yn,zn = get_xyz_coordinates(lonobs[icl+1,jcl],latobs[icl+1,jcl])
  xs,ys,zs = get_xyz_coordinates(lonobs[icl-1,jcl],latobs[icl-1,jcl])
  xe,ye,ze = get_xyz_coordinates(lonobs[icl,jcl+1],latobs[icl,jcl+1])
  xw,yw,zw = get_xyz_coordinates(lonobs[icl,jcl-1],latobs[icl,jcl-1])

  dist0=np.zeros([2,2])
  dist0[0,0] = np.sqrt ((xn-x0)**2 + (yn-y0)**2 + (zn-z0)**2)
  dist0[1,0] = np.sqrt ((xs-x0)**2 + (ys-y0)**2 + (zs-z0)**2)
  dist0[0,1] = np.sqrt ((xe-x0)**2 + (ye-y0)**2 + (ze-z0)**2)
  dist0[1,1] = np.sqrt ((xw-x0)**2 + (yw-y0)**2 + (zw-z0)**2)

  distb=np.zeros([2,2])
  distb[0,0] = np.sqrt ((xn-xi)**2 + (yn-yi)**2 + (zn-zi)**2)
  distb[1,0] = np.sqrt ((xs-xi)**2 + (ys-yi)**2 + (zs-zi)**2)
  distb[0,1] = np.sqrt ((xe-xi)**2 + (ye-yi)**2 + (ze-zi)**2)
  distb[1,1] = np.sqrt ((xw-xi)**2 + (yw-yi)**2 + (zw-zi)**2)

  nni = np.int(np.sign(dist0[1,0]-dist0[0,0]))
  nnj = np.int(np.sign(dist0[1,1]-dist0[0,1]))

  #mapping of 1-> 0 and -1-> 1
  indx = np.int(0.5*(1-nni))
  indy = np.int(0.5*(1-nnj))
  #print nni,nnj,indx,indy

  dista = np.sqrt ((x0-xi)**2 + (y0-yi)**2 + (z0-zi)**2)
  distb1 = distb[indx,0]
  distb2 = distb[indy,1]
  distc1 = dist0[indx,0]
  distc2 = dist0[indy,1]

  #formula of heron
  s1 = 0.5* (dista + distb1 + distc1)
  s2 = 0.5* (dista + distb2 + distc2)
  di = 2.0 * np.sqrt(s1*(s1-dista)*(s1-distb1)*(s1-distc1)) / distb1 / distb2
  dj = 2.0 * np.sqrt(s2*(s2-dista)*(s2-distb2)*(s2-distc2)) / distb2 / distb1

  """
  fig14 = plt.figure(14,figsize=(8, 4)),plt.clf()
  ax3=plt.subplot(111)
  ax3.plot(lon0,lat0,"ro",markersize=15)
  ax3.plot(lonij,latij,"go")
  ax3.plot(lonobs[icl+nni,jcl+nnj],latobs[icl+nni,jcl+nnj],"bo")
  ax3.plot(lonobs[icl+nni,jcl],latobs[icl+nni,jcl],"bo")
  ax3.plot(lonobs[icl,jcl+nnj],latobs[icl,jcl+nnj],"bo")
  lonint=get_interpolated_value(icl,jcl,nni*di,nnj*dj,lonobs)
  latint=get_interpolated_value(icl,jcl,nni*di,nnj*dj,latobs)
  ax3.plot(lonint,latint,"yo")
  plt.show()
  """

  return (nni*di,nnj*dj)



#############################################################
def get_mask_value(H,B,sl):

    intmask=0.0
    rhor = cf.rhoi/cf.rhosw
    if (H > 0.0 and H*rhor > (sl-B)):
      intmask = cf.mgr #2.0 # grounded
    elif H > 0.0 and H*rhor <= (sl-B):
      intmask = cf.mfl #3.0 # floating
    elif (H == 0.0 and 0.0 <= (sl-B)):
      intmask = cf.moc #4.0 #icefree ocean

    return intmask

############################################################
def mask_without_icerises(data, start_coords, fill_value):
    """
    Flood fill algorithm from https://gist.github.com/JDWarner/1158a9515c7f1b1c21f1
    
    Parameters
    ----------
    data : (M, N) ndarray of uint8 type
        Image with flood to be filled. Modified inplace.
    start_coords : tuple
        Length-2 tuple of ints defining (row, col) start coordinates.
    fill_value : int
        Value the flooded area will take after the fill.
        
    """
    xsize, ysize = data.shape
    orig_value = data[start_coords[0], start_coords[1]]
    
    stack = set(((start_coords[0], start_coords[1]),))
    if fill_value == orig_value:
        raise ValueError("Filling region with same value "
                         "already present is unsupported. "
                         "Did you already fill this region?")

    while stack:
        x, y = stack.pop()

        if data[x, y] == orig_value:
            data[x, y] = fill_value
            if x > 0:
                stack.add((x - 1, y))
            if x < (xsize - 1):
                stack.add((x + 1, y))
            if y > 0:
                stack.add((x, y - 1))
            if y < (ysize - 1):
                stack.add((x, y + 1))

    return data

#####################################################################
def get_points_along_transect(xobs,yobs,dkm,transects):

  trans_points=[]
  for l,points in enumerate(transects):

    pkm=[]
    for pp in points:
      pkm.append([xobs[np.int(pp[0])],yobs[np.int(pp[1])]])

    dp=np.float(dkm)/np.float(cf.resolution)
    dpr=dp

    po=[]
    po.append(points[0])
    pox=po[0][0]
    poy=po[0][1]

    count=0
    d=0

    while count<len(points)-1:

      xdir=np.sign(points[count+1][0]-points[count][0])
      ydir=np.sign(points[count+1][1]-points[count][1])
      my=np.float(points[count][1]-points[count+1][1])
      mx=np.float(points[count][0]-points[count+1][0])

      if mx==0:
        dpx=0.0
        dpy=ydir*dpr
        mnew=np.nan
      elif my==0:
        dpy=0.0
        dpx=xdir*dpr
        mnew=0.0
      else:
        mnew=my/mx
        dpy=ydir*np.sqrt((mnew*dpr)**2/(1.0+mnew**2))
        dpx=dpy/mnew

      dx=dpx*(xobs[np.int(points[count][0]+xdir)]-xobs[np.int(points[count][0])])
      dy=dpy*(yobs[np.int(points[count][1]+ydir)]-yobs[np.int(points[count][1])])

      px=dpx+pox
      py=dpy+poy

      if (xdir*(px-points[count+1][0])<=0.0 and ydir*(py-points[count+1][1])<=0.0):
        po.append([px,py])
        dpr=dp
        pox=po[d+1][0]
        poy=po[d+1][1]
      else: #overshoot, next point, take rest distance
        dpr=np.sqrt((px-points[count+1][0])**2+(py-points[count+1][1])**2)
        pox=points[count+1][0]
        poy=points[count+1][1]
        count+=1
        d-=1
      d+=1

    trans_points.append(po)

  return trans_points

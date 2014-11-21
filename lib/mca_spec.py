# -*- coding: utf-8 -*-
from __future__ import division
import scipy
from scipy import *
from scipy import signal
from numpy import linalg
import numpy as np
#if scipy.__version__ !='0.3.2':
#   import numpy
#from numpy.oldnumeric import *
import copy
import math
import time

def rebin3(a, (m, n, o)):
   if (m,n,o)==(1,1,1): return a
   M, N ,O= a.shape
   print a.shape
   print (M/m,m,N/n,n,O/o,o)
   print a[:M-M%(M/m),:N-N%(N/n),:O-O%(O/o)].shape
   a = reshape(a[:M-M%(M/m),:N-N%(N/n),:O-O%(O/o)], (M/m,m,N/n,n,O/o,o))
   return sum(sum(sum(a,5),3),1) / float(m*n*o)

def rebin2(a, (m, n)):
   M, N= a.shape
   #print a.shape
   #print (M/m,m,N/n,n)
   #print a[:M-M%(M/m),:N-N%(N/n)].shape
   a = reshape(a[:M-M%(M/m),:N-N%(N/n)], (M/m,m,N/n,n))
   return sum(sum(a,3),1) / float(m*n)

def list2float(l):
   for i in xrange(len(l)):
      l[i]=float(l[i])
   return l

def vcross(v,w):
   return reshape(cross(v.flat,w.flat),v.shape)

def rotmat(angle,axis):
   """ Create a rotation matrix from angle (in radian !) around an arbitrary vector.
   http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
   """
   x,y,z=axis.flat[0],axis.flat[1],axis.flat[2]
   s=sin(angle/2.)/sqrt(x*x+y*y+z*z)
   a,b,c,d=cos(angle/2.),s*x,s*y,s*z     #quaternion
   
   return array([[a**2+b**2-c**2-d**2 , 2*b*c-2*a*d         , 2*a*c+2*b*d        ],
                 [2*a*d+2*b*c         , a**2-b**2+c**2-d**2 , 2*c*d-2*a*b        ],
                 [2*b*d-2*a*c         , 2*a*b+2*c*d         , a**2-b**2-c**2+d**2]])

class SpecSetup:
   """ Class to hold basic spec parameters to be used as default parameters to functions in the spec module.
   """
   specFile=None
   mcaBaseName=None
   beamline=None
   filterColumn=None
   filterCoefficient=None
   monitorColumn=None
   nRebinMCA=None
   nRebinX=None
   psd0=None
   psdPixelSizeOverDistance=None
   psdMin=None
   psdMax=None
   mapType=None
   def Reset(self):
      self.specFile=None
      self.mcaBaseName=None
      self.beamline=None
      self.filterColumn=None
      self.filterCoefficient=None
      self.monitorColumn=None
      self.nRebinMCA=None
      self.nRebinX=None
      self.psd0=None
      self.psdPixelSizeOverDistance=None
      self.psdMin=None
      self.psdMax=None
      self.mapType=None
      
   
gSpecSetup=SpecSetup()
   

class SpecGeometry:
   """ Spec parameters (lattice parameters, orientation reflections and matrix as deduced from
   a Spec header, and current 'P' values for motors.
   
   """
   G={}# raw Geometrical parameters parameters (G0= mode, sector, azH,azK,azL, zoneH0,zoneK0,zoneL0,
       #                                        G1= lattice parameters (direct and reciprocal space), 1st reflection angles, 2nd reflection, reflection HKL's
       #                                        G3=  UB*2pi matrix , G4= H K L and ??????
   P=[]# raw 'P' parameters - motor values, depends on beamlines !
   q0=zeros((3,1),float32) # current Q vector, in HKL coordinates
   psd0=257.6            # position of the direct beam on the psd
   psdPixelSizeOverDistance=0.000115    # PSD pixel size, divided by the distance between sample and the PSD
   psdor="out" #orientation of psd: "out" horizontal (+k), "in" horizontal (-k), "up" vertical (+i), "down" vertical (-i)
   warning=True
   def __init__(self,headers):
      self.headers=headers
      for k in headers.keys():
         if k[0]=='G':
            self.G[k]=list2float(headers[k])
         if k[0]=='P':
            self.P=list2float(headers[k])
         if k[0]=='Q':
            self.q0[0][0],self.q0[1][0],self.q0[2][0]=float(headers[k][0]),float(headers[k][1]),float(headers[k][2])
      self.lattice  =[self.G["G1"][0],self.G["G1"][1],self.G["G1"][2],self.G["G1"][3]*pi/180,self.G["G1"][4 ]*pi/180,self.G["G1"][5 ]*pi/180]
      self.lattice_r=[self.G["G1"][6],self.G["G1"][7],self.G["G1"][8],self.G["G1"][9]*pi/180,self.G["G1"][10]*pi/180,self.G["G1"][11]*pi/180]
      # Spec orientation matrix
      self.ub=array([[self.G["G3"][0],self.G["G3"][1],self.G["G3"][2]],[self.G["G3"][3],self.G["G3"][4],self.G["G3"][5]],[self.G["G3"][6],self.G["G3"][7],self.G["G3"][8]]])/(2*pi)
      self.wavelength=self.G["G4"][3]
      #Current motor positions, as a dictionnary
      self.angles0={}
      for i in xrange(len(self.P)):# :TODO: Some motor names may include spaces, and break conversion
         name=headers["O"][i]
         if name=="phi":name="Phi"
         if name=="chi":name="Chi"
         if name=="mu":name="Mu"
         if not self.angles0.has_key(name):# Kludge - e.g. if one motor is "Theta" and another is "Theta def."...
            self.angles0[name]=self.P[i]*pi/180
   def Z(self,angles):
      """ Get the orientation matrix due to the sample motor angles
      """
      if self.headers["geom"]=="kappapsic":
         m_phi=rotmat(angles["Phi"],array([[0],[0],[-1]]))
         m_chi=rotmat(angles["Chi"],array([[0],[1],[0]]))
         m_eta=rotmat(angles["Eta"],array([[0],[0],[-1]]))
         m_mu =rotmat(angles["Mu"] ,array([[1],[0],[0]]))
         return dot(m_mu,dot(m_eta,dot(m_chi,m_phi)))
      elif self.headers["geom"]=="dafs" or self.headers["geom"]=="dafsd2ps" or self.headers["geom"]=="dafsd2psold":
         m_phi=rotmat(angles["Phi"],array([[0],[0],[-1]]))
         m_chi=rotmat(angles["Chi"],array([[0],[1],[0]]))
         m_eta=rotmat(angles["Eta"],array([[0],[0],[-1]]))
         m_mu =rotmat(angles["Mu"] ,array([[1],[0],[0]]))
         return dot(m_eta,dot(m_chi,dot(m_mu,m_phi)))
      elif self.headers["geom"]=="gonio" or self.headers["geom"]=="gmci":
         #print "Z gonio: ome=%6.2f, the=%6.2f"%(angles["Omega"]*180/pi,angles["Theta"]*180/pi)
         m_ome=rotmat(angles["Omega"],array([[0],[0],[-1]]))
         m_the=rotmat(angles["Theta"] ,array([[1],[0],[0]]))
         return dot(m_the,m_ome)
      #elif self.headers["geom"]=="dafs":
      #else:
      #   except("!! Unknown spec geometry:",self.headers.geom)
   def QL(self,angles,psd=None):
      """ Get the QL coordinates of the point on Ewalds sphere seen by the detector, in Angstroem-1
      """
      if self.headers["geom"]=="kappapsic":
         ql=array([[sin(angles["Del"])],[cos(angles["Del"])*cos(angles["Nu"])],[cos(angles["Del"])*sin(angles["Nu"])]])
      elif self.headers["geom"]=="dafs" or self.headers["geom"]=="dafsd2ps" or self.headers["geom"]=="dafsd2psold":
         ql=array([[sin(angles["Del"])],[cos(angles["Del"])*cos(angles["Nu"])],[cos(angles["Del"])*sin(-angles["Nu"])]])
      elif self.headers["geom"]=="gonio" or self.headers["geom"]=="gmci":
         #print "QL gonio: psi=%6.2f, phi=%6.2f"%(angles["Psi"]*180/pi,angles["Phi"]*180/pi)
         ql=array([[sin(angles["Psi"])]  ,[cos(angles["Psi"])  *cos(angles["Phi"])],[cos(angles["Psi"])  *sin(angles["Phi"])]])
      #elif self.headers["geom"]=="dafs":
      #else:
      #   except("!! Unknown spec geometry:",self.headers.geom)
      if psd !=None:# this is for an *horizontal* PSD
         if isscalar(psd):
            if self.psdor=="out":
               ql+=array([[0.0],[0.0],[self.psdPixelSizeOverDistance*(psd-self.psd0)]])
            elif self.psdor=="in":
               ql-=array([[0.0],[0.0],[self.psdPixelSizeOverDistance*(psd-self.psd0)]])
            elif self.psdor=="up":
               ql+=array([[self.psdPixelSizeOverDistance*(psd-self.psd0)],[0.0],[0.0]])
            elif self.psdor=="down":
               ql-=array([[self.psdPixelSizeOverDistance*(psd-self.psd0)],[0.0],[0.0]])
            else:raise RunTimeError
            ql/=linalg.norm(ql)
         else:
            da=zeros((3,len(psd.flat)),float32)
            if self.psdor=="out":
               da[2,:]= self.psdPixelSizeOverDistance*(psd.flatten()-self.psd0)
            elif self.psdor=="in":
               da[2,:]=-self.psdPixelSizeOverDistance*(psd.flatten()-self.psd0)
            elif self.psdor=="up":
               da[0,:]= self.psdPixelSizeOverDistance*(psd.flatten()-self.psd0)
            elif self.psdor=="down":
               da[0,:]=-self.psdPixelSizeOverDistance*(psd.flatten()-self.psd0)
            else:raise RunTimeError
            ql=ql+da
            qln=sqrt((ql**2).sum(axis=0))
            ql[0,:]/=qln
            ql[1,:]/=qln
            ql[2,:]/=qln
      ql[1,:]-=1.0 #reciprocal space reference frame has origin on Ewald's sphere, not at its center
      return ql/self.wavelength
   def Angles2HKL(self,angles,psd=None,verbose=False):
      """
      Convert angles and psd pixel to hkl values
      psd can be an array of pixel values
      """
      m=dot(self.Z(angles),self.ub)
      im=linalg.inv(m)
      ql=self.QL(angles,psd=psd)
      h=dot(im,ql)
      if self.headers["geom"]=="gonio" or self.headers["geom"]=="gmci":
         if self.warning:
            print "WARNING gn,onio/gmci: rotating h/k"
            self.warning=False
         # WHY ??? Not using Busing & Levy axis convention ?
         horig=copy.deepcopy(h)
         h[0,:]=-horig[1,:]
         h[1,:]=horig[0,:]
      if verbose:
         if self.headers["geom"]=="gonio" or self.headers["geom"]=="gmci":
            print "Phi=%6.2f Psi=%6.2f Omega=%6.2f Theta=%6.2f H=%5.2f K=%5.2f L=%5.2f"%(angles["Phi"]*180/pi,angles["Psi"]*180/pi,angles["Omega"]*180/pi,angles["Theta"]*180/pi,h[0],h[1],h[2])
            cu = cos(angles["Psi"])
            su = sin(angles["Psi"])
            co = cos(angles["Omega"])
            so = sin(angles["Omega"])
            ca = cos(angles["Theta"])
            sa = sin(angles["Theta"])
            cd = cos(angles["Phi"] - angles["Theta"])
            sd = sin(angles["Phi"] - angles["Theta"])
            hphi=array([[ -co * su - so * ca + so * cu * cd],[so * su - co * ca + co * cu * cd],[ sa + cu * sd]])/self.wavelength
            print hphi
            h=dot(linalg.inv(self.ub),hphi)
            print "Phi=%6.2f Psi=%6.2f Omega=%6.2f Theta=%6.2f H=%5.2f K=%5.2f L=%5.2f"%(angles["Phi"]*180/pi,angles["Psi"]*180/pi,angles["Omega"]*180/pi,angles["Theta"]*180/pi,h[0],h[1],h[2])
      return h


def ReadSpec(fname,scan):
   f=open(fname,'r')
   s="#S %i"%scan
   print s
   title=0
   headers={}
   headers["geom"]=""
   while 1:
      title=f.readline()
      if s == title[0:len(s)]:
         break;
      if "#O0"==title[0:3]:#Motor names, in the order reported in 
         #:TODO: some motors have spaces in their names !!
         headers["O"]=title.split()[1:]
      else:
         if "#O"==title[0:2]:#Motor names, continued
            headers["O"]+=title.split()[1:]
      if headers["geom"]=="" and "#C"==title[0:2]:# first comment gives used geometry / spec session name (gnio, psic, dafs)
         headers["geom"]=title.split()[1]
      if len(title)==0:
         break;
   headers["S"]=title[2:]
   s="#L"
   coltit=0
   while 1:
      coltit=f.readline()
      if s == coltit[0:len(s)]:
         break;
      if len(coltit)==0:
         break;
      if coltit[0:2]=="#D" or coltit[0:2]=="#G" or coltit[0:2]=="#Q":
         coltit=coltit.split()
         headers[coltit[0][1:]]=coltit[1:]
      if coltit[0:3]=="#P0":
         headers["P"]=coltit.split()[1:]
      else :
         if coltit[0:2]=="#P":
            headers["P"]+=coltit.split()[1:]
   d={}
   coltit=coltit.split()
   for i in xrange(1,len(coltit)):
   	d[coltit[i]]=[]
   while 1:
      l=f.readline()
      if len(l)<2:
         break;
      if l[0]=="#":continue
      if l[0]=="@":
         while l[-2]=='\\':l=f.readline()
         continue
      l=l.split()
      for i in xrange(1,len(coltit)):
         d[coltit[i]].append(float(l[i-1]))
   nb=len(d[coltit[1]])
   for i in xrange(1,len(coltit)):
      a=zeros(nb,float32)
      for j in xrange(nb):
        a[j]=d[coltit[i]][j]
      d[coltit[i]]=copy.deepcopy(a)
   f.close()
   # Transform #O and #P headers into a dictionnary of motor positions, *before* the beginning of the scan
   if headers.has_key("P") and headers.has_key("O"):
      if len(headers["P"])==len(headers["O"]):
         headers["motors"]={}
         for i in xrange(len(headers["P"])):
            headers["motors"][headers["O"][i]]=float(headers["P"][i])
   return headers,d

def ReadSpecHKLMAD(fname,sn,fout="HKLMAD.dat",ycol="vct1_4",normcol="vct1_3",filtercol=None,filtercoeff=10.):
   headers,d=ReadSpec(fname,sn[0])
   t=headers["S"]
   n=len(d["H"])
   ne=len(sn)
   a=zeros([n,ne+3],float32)
   a[:,0]=d["H"]
   a[:,1]=d["K"]
   a[:,2]=d["L"]
   for i in xrange(ne):
      headers,d=ReadSpec(fname,sn[i])
      t=headers["S"]
      a[:,3+i]=d[ycol]
      if (normcol!="") and (normcol!=None):
         a[:,3+i] /= d[normcol]
      if filtercol!=None:
         a[:,3+i]*=filtercoeff**d[filtercol]
      #a[:,3+i]=d["Mcacorr"]/d["Mon2"]
   f=open(fout,'w')
   f.write("%i\n"%n)
   for i in xrange(n):
      for j in xrange(ne+3):
         f.write("%20.12f "%a[i,j])
      f.write("\n")
   print "Total:%i scan imported"%ne

def ReadSpec2D(fname,sn,ix="H",iy="L",ycol="vct1_4",filtercol=None,filtercoeff=10.,normcol=None,resizex=None,mcabasename=None,psdmin=None,psdmax=None):
   """ Read a 2D map from series of 1D scans
   
   If resize!=None, resize to this length of array. Non-observed values are masked.
   This can be used for a series of scans with an increasing numer of points
   
   If mcabasename!=None, the ycol will be ignored, and the intensity will be integrated
   from the original MCA data, using the optionnal [psdmin,psdmax] range
   
   """
   headers,d=ReadSpec(fname,sn[0])
   t=headers["S"]
   nx=len(d[ix])
   ny=len(sn)
   a=zeros([nx,ny],float32)
   x=zeros([nx,ny],float32)
   y=zeros([nx,ny],float32)
   if psdmin==None:psdmin=0
   if psdmax==None:psdmin=-1
   if resizex != None:
      # use a masked array, all values are true by default
      nx=resizex
      a=zeros([nx,ny],float32)
      x=zeros([nx,ny],float32)
      y=zeros([nx,ny],float32)
      a=ma.array(a,mask=a>1.0)
      x=ma.array(x,mask=a>1.0)
      y=ma.array(y,mask=a>1.0)
   for i in xrange(ny):
      headers,d=ReadSpec(fname,sn[i])
      t=headers["S"]
      nxtmp=len(d[ix])
      if mcabasename==None:
         a[:nxtmp,i]=d[ycol]
      else:
         for n in xrange(nxtmp):
            f=mcabasename%(sn[i],n)
            a[n,i]=ReadMCA(f)[psdmin:psdmax].sum()
      if (normcol != "") and (normcol != None):
         a[:nxtmp,i]/=d[normcol]
      if filtercol!=None:
         if d.has_key(filtercol):
            a[:nxtmp,i] *= filtercoeff**(d[filtercol])
         else:
            if headers["motors"].has_key(filtercol):
                  a[:nxtmp,i]*=filtercoeff**headers["motors"][filtercol]
            else:
               print "filter column not found !!!"
               print headers["motors"].keys()," ->NOT FOUND:",iy
      if d.has_key(iy):
         y[:nxtmp,i]=d[iy]
      else:
         if headers["motors"].has_key(iy):
            y[:nxtmp,i]=headers["motors"][iy]
         else:
            print headers["motors"].keys()," ->NOT FOUND:",iy
      if d.has_key(ix):
         x[:nxtmp,i]=d[ix]
      else:
         if headers["motors"].has_key(ix):
            x[:nxtmp,i]=headers["motors"][ix]
         else:
            print headers["motors"].keys()," ->NOT FOUND:",ix
      if resizex != None:
        x.mask[nxtmp:,i]=True
        y.mask[nxtmp:,i]=True
        a.mask[nxtmp:,i]=True
   print "Total:%i scan imported"%ny
   return a,x,y

def ReadSpecMesh(fname,sn,ycol="vct1_4",filtercol=None,filtercoeff=10.,normcol=None):
   """ Read a 2D map from a spec mesh
   
   If mcabasename!=None, the ycol will be ignored, and the intensity will be integrated
   from the original MCA data, using the optionnal [psdmin,psdmax] range
   """
   headers,d=ReadSpec(fname,sn)
   t=headers["S"]
   print t
   ix=t.split()[-9]
   iy=t.split()[-5]
   a=d[ycol]
   x=d[ix]
   y=d[iy]
   nx=int(t.split()[-6])+1
   #ny=int(t.split()[-2])+1
   ny=len(a)//nx
   if filtercol!=None:
     if d.has_key(filtercol):
       a *= filtercoeff**(d[filtercol])
     else:
       if headers["motors"].has_key(filtercol):
         a[:nxtmp,i]*=filtercoeff**headers["motors"][filtercol]
       else:
         print "filter column not found !!!"
         print headers["motors"].keys()," ->NOT FOUND:",iy
   if (normcol != "") and (normcol != None):
     a/=d[normcol]
   # Crop - if scan was stopped before the end
   a=a[:nx*ny]
   x=x[:nx*ny]
   y=y[:nx*ny]
   return a.reshape((ny,nx)),x.reshape((ny,nx)),y.reshape((ny,nx)),ix,iy

def ReadMCA(filename):
   #print "ReadMCA: reading %s"%filename
   f=open(filename,'r')
   nb=0
   pointperline=16
   scan=zeros(1,'f')
   i=0
   while True:
      l=f.readline()
      s=l.split()
      if len(s)>0:
         if s[0]=="#@CHANN":
            nb=float(s[1])
            scan=zeros(nb,'f')
         if s[0]=="@A":
            l=l[2:len(l)-2]
            s=l.split()
            for c in s:
               scan[i]=int(c)
               i=i+1
            break
   while i<nb:
      l=f.readline()
      if l[len(l)-2]=="\\":
         l=l[:len(l)-2]
      s=l.split()
      for c in s:
         scan[i]=int(c)
         i=i+1
         if(i>=nb):
            break
   return scan

def ReadMCA2D (mcabasename,scanid,specfilename,filtercol=None,filtercoeff=10.,normcol="Mon2",nrebin_mca=None,psdmin=None,psdmax=None):
   """ Read 2D MCA scan, rebin it and use Delaunay interpolation to obtain a representation along orthonormal coordinates.
   
       psdmin: first pixel included in the data (to avoid null regions on the borders)
       psdmax: first pixel *not* included in the data => e.g. keep a[:,psdmin:psdmax]
       filtercoeff: value of the factor when changing filters
       normcol: monitor column in the Spec scan
       filtercol: filter column in the spec file
       
       Rebinning (numbers should be exact dividers of the recorded number of pixels):
         nrebin_mca: final number of desired pixel in the PSD direction
   """
   scan2d=[]
   headers,scan=ReadSpec(specfilename,scanid)
   title=headers["S"]
   norm=scan[normcol]
   if filtercol!=None:
      norm/=filtercoeff**scan[filtercol]
   nb=len(norm)
   for n in xrange(nb): 			     	
      f=mcabasename
      s=ReadMCA(f)/norm[n]
      if n==0:
         scan2d=zeros([nb,len(s)],float32)
      scan2d[n,:]=s
   if psdmin==None:
      psdmin=0
   if psdmax==None:
      psdmax=scan2d.shape[1]
   if nrebin_mca!=None:
      s2dr=zeros((scan2d.shape[0],nrebin_mca),float32)
      step=(psdmax-psdmin)/nrebin_mca
      psdmax=psdmin+nrebin_mca*step #Make sure of the exact number of integration intervals
      psdidx=linspace(psdmin+step/2.0,psdmax-step/2.0,nrebin_mca)
      for i in xrange(nrebin_mca):
         s2dr[:,i]=scan2d[:,i*step:(i+1)*step].sum(axis=1)
      return headers,scan,s2dr
   else:
      return headers,scan,scan2d

def ReadMCA2D_complete(filename):
	""" Read all carto in a single MCA file
	return: scan3d --> numpy array of many 2D array. Each 2D array is a complete carto"""
	f = open(filename, 'r')
	scan3d   = []
	scan2d   = []
	scan1d   = []
	compteur = 0
	point_nb = 0
	for line in f.readlines():
		if not line.strip():
			continue
		else:
			if line.startswith("#S"):
				point_nb_actual = int(line.split()[6])
				if point_nb_actual != point_nb or compteur>point_nb:
					scan2d.append(scan1d)
					scan3d.append(scan2d)
					scan1d = []
					scan2d = []
					compteur = 0
					point_nb = point_nb_actual
			elif line.startswith("#"):
				continue
			else:
				l = line.split()
				if l[0]=="@A":
					compteur = compteur+1
					scan2d.append(scan1d)
					scan1d = []
					line = line[2:len(line)-2]
					#print line
					l = line.split()
					for i in l:
						scan1d.append(int(i))
				else:
					if line[len(line)-2]=="\\":
						#print line
						line = line[:len(line)-2]
					#print "###"
					#print line
					l = line.split()
					for i in l:
						scan1d.append(int(i))
	scan2d.append(scan1d) #dernier scan
	scan3d.append(scan2d) #derniere carto
	scan3d = scan3d[1:]
	for i in range(len(scan3d)):
		scan3d[i] = np.asarray(scan3d[i][1:])
	scan3d = np.asarray(scan3d)
	f.close()
	return scan3d

def Interp2d(x0,y0,z0,x1=None,y1=None,gridsize=None,nbknots=20,s=None,smargin=0.1,ks=3):
   """ Interpolate 2D data z0(x0,y0) to a regular grid, assuming x0 and z0 present smooth variations other
   the 2D array, which can be approximated using splines.
   
   x0,z0,y0: 2D arrays giving the values of z0(x0,y0) - x0, y0 must be slowly varying functions of the indices
             of the arrays.
   x1,y1: 1D arrays defining the x,y grid over which the interpolation will be computed
   gridsize: if x1 or y1 is not given, x1=linspace(s0.min(),x0.max(),gridsize)
   nbknots: the 2D spline used to calculate the indices at which a given (x,y) can be foud will use nbknots*nbknots
   s: smooth factor
   """
   x0,y0=y0,x0 # Someday I will understand the order of indices...
   if x1==None:
      if gridsize==None:gridsize=400
      x1=linspace(x0.min(),x0.max(),gridsize)
   if y1==None:
      if gridsize==None:gridsize=400
      y1=linspace(y0.min(),y0.max(),gridsize)
   
   splx0=zeros((nbknots,nbknots),float32)
   sply0=zeros((nbknots,nbknots),float32)
   splix =zeros((nbknots,nbknots),float32)
   spliy =zeros((nbknots,nbknots),float32)
   for i in xrange(nbknots):
      ix=int(i*(z0.shape[0]-1)/float(nbknots))
      splix[i,:]=ix
      for j in xrange(nbknots):
         iy=int(j*(z0.shape[1]-1)/float(nbknots))
         splx0[i,j]=x0[ix,iy]
         sply0[i,j]=y0[ix,iy]
         spliy[i,j]=iy

   x0min,x0max,y0min,y0max=x0.min(),x0.max(),y0.min(),y0.max()
   #Spline which can give the point indices (ix,iy) in the array, given x and y coordinates
   dx0,dy0=x0max-x0min,y0max-y0min
   # Supply by hand xe and ye to avoid flat splines immediately beyond the known values
   xb,xe=x0min-dx0*smargin,x0max+dx0*smargin
   yb,ye=y0min-dy0*smargin,y0max+dy0*smargin
   tx,fpx,ier,msg=interpolate.bisplrep(splx0.flatten(),sply0.flatten(),splix.flatten(),xb=xb,xe=xe,yb=yb,ye=ye,full_output=1,kx=ks,ky=ks,s=s)
   print "#1, residual=",fpx," error=",ier,splix.min(),splix.max()
   ty,fpy,ier,msg=interpolate.bisplrep(splx0.flatten(),sply0.flatten(),spliy.flatten(),xb=xb,xe=xe,yb=yb,ye=ye,full_output=1,kx=ks,ky=ks,s=s)
   print "#2, residual=",fpy," error=",ier,spliy.min(),spliy.max()
   #ty,fp,ier,msg=interpolate.bisplrep(splx0.flatten(),sply0.flatten(),spliy.flatten(),full_output=1,s=200,task=-1,tx=splx0.flatten(),ty=sply0.flatten())
   #print "#3, residual=",fp," error=",ier,spliy.min(),spliy.max()
   if True: #fpx+fpy>0.1:
      # Try to get back x and y to check !
      x0test,y0test,ixtest,iytest=splx0,sply0,splix,spliy
      dx0bis=x0test*0
      dy0bis=y0test*0
      for i in xrange(x0test.shape[0]):
         for j in xrange(x0test.shape[1]):
            dx0bis[i,j]=interpolate.bisplev(x0test[i,j],y0test[i,j],tx)-splix[i,j]
            dy0bis[i,j]=interpolate.bisplev(x0test[i,j],y0test[i,j],ty)-spliy[i,j]
      maxdx=abs(dx0bis).max()
      maxdy=abs(dy0bis).max()
      errx=sqrt(((dx0bis)**2).mean())
      erry=sqrt(((dy0bis)**2).mean())
      print "<dix^2>=%6.4f%%, <diy^2>=%6.3f%% , maxdx=%6.3f%%, maxdy=%6.3f%%"%(errx*100,erry*100,maxdx,maxdy)
      fpx,fpy=maxdx,maxdy
   
   #Grid coordinates where we want the values
   ix1=interpolate.bisplev(x1,y1,tx)
   iy1=interpolate.bisplev(x1,y1,ty)
   ix1i=floor(ix1).astype(Int)                # :CHECK: Is casting floor to integer safe ??
   iy1i=floor(iy1).astype(Int)
   imask=(ix1i>=0)*(ix1i<(z0.shape[0]-1))*(iy1i>=0)*(iy1i<(z0.shape[1]-1))
   # Bilinear interpolation
   dix=ix1-ix1i
   diy=iy1-iy1i
   i00=(ix1i  )*z0.shape[1]+(iy1i  )
   i01=(ix1i  )*z0.shape[1]+(iy1i+1)
   i10=(ix1i+1)*z0.shape[1]+(iy1i  )
   i11=(ix1i+1)*z0.shape[1]+(iy1i+1)
   i00=0+i00*imask    # does take() accept masked array input ?
   i01=0+i01*imask
   i10=0+i10*imask
   i11=0+i11*imask
   z00=z0.flatten().take(i00)
   z01=z0.flatten().take(i01)
   z10=z0.flatten().take(i10)
   z11=z0.flatten().take(i11)
   z1=ma.masked_where(imask==False,z00*(1-dix)*(1-diy)+z01*(1-dix)*diy+z10*dix*(1-diy)+z11*dix*diy)
   return y1,x1,z1,fpx+fpy


def ReadMCA2DRebin(mcabasename,scanid,specfilename,filtercol=None,filtercoeff=10.,normcol="Mon2",nrebin_mca=None,psdmin=None,psdmax=None,psd0=None,psdPixelSizeOverDistance=None,psdor="out"):
   """ Read 2D MCA scan, rebin it and use Delaunay interpolation to obtain a representation along orthonormal coordinates.
   
       psdmin: first pixel included in the data (to avoid null regions on the borders)
       psdmax: first pixel *not* included in the data => e.g. keep a[:,psdmin:psdmax]
       filtercoeff: value of the factor when changing filters
       normcol: monitor column in the Spec scan
       filtercol: filter column in the spec file
       
       Rebinning (numbers should be exact dividers of the recorded number of pixels):
         nrebin_mca: final number of desired pixel in the PSD direction
         nrebin_x: final number of desired pixel in the x (non-PSD) direction. If 0 (default), no rebin is done
   """
   headers,s,s2d=ReadMCA2D(mcabasename,scanid,specfilename,filtercol=filtercol,filtercoeff=filtercoeff,normcol=normcol)
   if psdmin==None:
      psdmin=0
   if psdmax==None:
      psdmax=s2d.shape[1]
   s2d=s2d[:,psdmin:psdmax]
   s2dr=s2d
   psdidx=arange(s2d.shape[1])
   if nrebin_mca!=None:
      s2dr=zeros((s2d.shape[0],nrebin_mca),float32)
      step=s2d.shape[1]/nrebin_mca
      psdmax=psdmin+nrebin_mca*step #Make sure of the exact number of integration intervals
      psdidx=linspace(psdmin+step/2.0,psdmax-step/2.0,nrebin_mca)
      for i in xrange(nrebin_mca):
         s2dr[:,i]=s2d[:,i*step:(i+1)*step].sum(axis=1)
         
   
   g=SpecGeometry(headers)
   print specfilename, g.headers["S"]
   g.psd0=psd0
   g.psdPixelSizeOverDistance=psdPixelSizeOverDistance
   g.psdor=psdor
   h,k,l=s2dr*0,s2dr*0,s2dr*0
   nx,ny=h.shape
   angles0=g.angles0
   #This should work for both normal and hklscan, as long as the latter include all motor positions
   for ix in xrange(nx):
      # Find columns corresponding to motor names - convert all assuming they are angles in degrees
      for tit in s.keys():
         if angles0.has_key(tit):
            angles0[tit]=s[tit][ix]*pi/180
      hkl=g.Angles2HKL(angles0,psd=psdidx)
      h[ix,:],k[ix,:],l[ix,:]=hkl[0,:],hkl[1,:],hkl[2,:]
   
   return h,k,l,s2dr,g,s

def ReadMCA3DRebin(mcabasename,scanid,specfilename,filtercol=None,filtercoeff=10.,normcol="Mon2",nrebin_mca=None,psdmin=None,psdmax=None,psd0=None,psdPixelSizeOverDistance=None,psdor="out"):
   """ Read 3D MCA scan, rebin it and use Delaunay interpolation to obtain a representation along orthonormal coordinates.
   
       psdmin: first pixel included in the data (to avoid null regions on the borders)
       psdmax: first pixel *not* included in the data => e.g. keep a[:,psdmin:psdmax]
       filtercoeff: value of the factor when changing filters
       normcol: monitor column in the Spec scan
       filtercol: filter column in the spec file
       
       Rebinning (numbers should be exact dividers of the recorded number of pixels):
         nrebin_mca: final number of desired pixel in the PSD direction
         nrebin_x: final number of desired pixel in the x (non-PSD) direction. If 0 (default), no rebin is done
   """
   h0,k0,l0,s2dr0,g0,s0=ReadMCA2DRebin(mcabasename,scanid[0],specfilename,filtercol=filtercol,filtercoeff=filtercoeff,normcol=normcol,nrebin_mca=nrebin_mca,psdmin=psdmin,psdmax=psdmax,psd0=psd0,psdPixelSizeOverDistance=psdPixelSizeOverDistance,psdor=psdor)
   n=len(scanid)
   if n==1:return h0,k0,l0,s2dr0,g0,s0
   h=zeros((n,s2dr0.shape[0],s2dr0.shape[1]),dtype=float32)
   k=zeros((n,s2dr0.shape[0],s2dr0.shape[1]),dtype=float32)
   l=zeros((n,s2dr0.shape[0],s2dr0.shape[1]),dtype=float32)
   sr=zeros((n,s2dr0.shape[0],s2dr0.shape[1]),dtype=float32)
   h[0],k[0],l[0],sr[0]=h0,k0,l0,s2dr0
   for i in xrange(1,n):
      h0,k0,l0,s2dr0,g0,s0=ReadMCA2DRebin(mcabasename,scanid[i],specfilename,filtercol=filtercol,filtercoeff=filtercoeff,normcol=normcol,nrebin_mca=nrebin_mca,psdmin=psdmin,psdmax=psdmax,psd0=psd0,psdPixelSizeOverDistance=psdPixelSizeOverDistance,psdor=psdor)
      h[i],k[i],l[i],sr[i]=h0,k0,l0,s2dr0
   
   return h,k,l,sr,g0,s0

def ReadSpecHKLMAD_MCA2D(mcabasename,sn,specfilename,fout="HKLMAD.dat",psdmin=None,psdmax=None,nrebin_mca=None,filtercol=None,filtercoeff=None,normcol=None,psd0=None,psdPixelSizeOverDistance=None,psdor="out"):
   """ Create HKLMAD file from HKL-scans recorded with a PSD in MCA scans (ID1 data)
   
       psdmin: first pixel included in the data (to avoid null regions on the borders)
       psdmax: last pixel  included in the data
       filtercoeff: value of the factor when changing filters
       normcol: monitor column in the Spec scan
       filtercol: filter column in the spec file
       
       Rebinning (numbers should be exact dividers of the recorded number of pixels):
         nrebin_mca: final number of desired pixel in the PSD direction
   """
   ne=len(sn)
   h,k,l,d,g=[],[],[],[],[]
   for i in xrange(ne):
      tmph,tmpk,tmpl,tmpd,tmpg,tmps=ReadMCA2DRebin(mcabasename,sn[i],specfilename,filtercol=filtercol,filtercoeff=filtercoeff,normcol=normcol,nrebin_mca=nrebin_mca,psdmin=psdmin,psdmax=psdmax,psd0=psd0,psdPixelSizeOverDistance=psdPixelSizeOverDistance,psdor=psdor)
      h.append(tmph)
      k.append(tmpk)
      l.append(tmpl)
      d.append(tmpd)
      g.append(tmpg)
      print tmph.shape,tmpk.shape,tmpl.shape,tmpd.shape
   n=len(h[0].flatten())
   f=open(fout,'w')
   f.write("%i\n"%(n))
   nec=int(ne/2)
   for i in xrange(n):
      f.write("%8.5f %8.5f%8.5f "%(h[nec].flat[i],k[nec].flat[i],l[nec].flat[i]))
      for j in xrange(ne):
         f.write("%20.12f "%(d[j].flat[i]))
      f.write("\n")
   print "Total:%i scan (%i x %i) imported"%(ne,h[0].shape[0],h[0].shape[1])

def ReadMAD2D(fin="FTFAFNDPHI.dat.calc",npsd=None):
   f=open(fin,"r")
   nbRefl=int(f.readline())
   nbx=npsd
   nby=nbRefl/nbx
   d = array([float(x) for x in f.read().split()], float32)
   f.close()
   d.shape = (len(d)/13,13)  # 13 colonnes H K L FT sigFT FA sigFA FN sigFN PhiT-PhiA sig(PhiT-PhiA) PhiN-PhiA sig(PhiN-PhiA)
   h          =reshape(d[:,0 ],(nby,nbx))
   k          =reshape(d[:,1 ],(nby,nbx))
   l          =reshape(d[:,2 ],(nby,nbx))
   ft         =reshape(d[:,3 ],(nby,nbx))
   sigft      =reshape(d[:,4 ],(nby,nbx))
   fa         =reshape(d[:,5 ],(nby,nbx))
   sigfa      =reshape(d[:,6 ],(nby,nbx))
   fn         =reshape(d[:,7 ],(nby,nbx))
   sigfn      =reshape(d[:,8 ],(nby,nbx))
   PhiTPhiA   =reshape(d[:,9 ],(nby,nbx))%360
   sigPhiTPhiA=reshape(d[:,10],(nby,nbx))
   PhiNPhiA   =reshape(d[:,11],(nby,nbx))%360
   sigPhiNPhiA=reshape(d[:,12],(nby,nbx))
   PhiTPhiA -= (abs(PhiTPhiA-360)<90)*360
   PhiNPhiA -= (abs(PhiNPhiA-360)<90)*360
   return h,k,l,ft,sigft,fa,sigfa,fn,sigfn,PhiTPhiA,sigPhiTPhiA

################################# XPAD ###############################
def importXPADI0(specfiles=[],verbose =False):
  xpad_i0={}
  for specfile in specfiles:
    ll=open(specfile,'r').readlines()
    for i in xrange(len(ll)):
      l=ll[i]
      if l[:7]=="#XPADCT":
        if len(ll)>(i+7):
          l1=ll[i+6]
          l2=ll[i+7]
          if len(l2)>6:
            if l1[:2]=="#L" and l2[:6]!="#ABORT":
              l1=l1.split()[1:]
              l2=l2.split()
              d={}
              for j in xrange(len(l1)):
                k=l1[j]
                v=l2[j]
                if v.find('(')>=0: d[k]=float(l2[j].split('(')[0])
                else: d[k]=float(l2[j])
                imgname=l.split()[1].split(':')[1]+".txt"
                if verbose:print imgname,d
                i+=7
                xpad_i0[imgname]=d
  return xpad_i0

def XPADcalcFlatFieldDarkMask(dark,flat,medfilt_width=21):
  print "XPADcalcFlatFieldDarkMask - median filter will take a few seconds..."
  # flat field calculated using a large median filter
  flatmed=signal.medfilt2d(flat,medfilt_width)
  #Include in the dark mask pixels counting less than half or more than 50%
  #than the median filtered value
  mask=(flat<(.5*flatmed))+(flat>(1.5*flatmed))
  # Also include pixels counting in the dark...
  mask+=(dark>5)
  return flatmed/(flat+1e-8)*(mask==0),mask>0

def XPADCorrFlatMask(d,flat,mask,fastmask=False):
  dcorr=d*flat
  if mask!=None:
    if fastmask: dcorr*=(mask==0)
    else:
      h,w=dcorr.shape
      idx=nonzero(mask)
      # Use a weighted average from neighbours, the weight multiplied by a factor
      # 0.3 with each unit departure from the original pixel
      for i in xrange(len(idx[0])):
        x0,y0=idx[0][i],idx[1][i]
        #print "XPADCorrFlatMask : %d/%d  -> %d,%d"%(i,len(idx[0]),x0,y0)
        for dxy in xrange(1,10+1):#try to find at least 4 pixels to average from
          x=arange(-dxy,dxy+1,dtype=int)
          y=arange(-dxy,dxy+1,dtype=int)[:,newaxis]
          weight=0.3**sqrt(x**2+y**2)
          x,y=(x+x0+0*y)%h,(y+y0+0*x)%w
          idx1=(x*w+y).flatten()
          #print idx1
          weight2=weight.flatten()*(take(mask,idx1)==False)
          nb=len(flatnonzero(weight2))
          if nb>4:
            dcorr[x,y]=(take(dcorr,idx1)*weight2).sum()/weight2.sum()
            break
  return dcorr
    

def readXPAD2D(fname,shape=(960,560),corrBigPixels=True,dark=None,corrOverlap=5,flat=None,mask=None,fastmask=None):
   f=open(fname)
   a=fromfile(f,sep=" ",dtype=float32)
   a=a.reshape(shape)
   if dark!=None: a-=dark
   
   if corrBigPixels: # correct for pixels 2.42 larger than others
    for i in xrange(1,6+1):
      a[:,80*i-1:80*i+1]/=2.42
   
   if corrOverlap>0:# remove overlapping pixels ?
     for i in xrange(1,7+1):
       a[120+(120-corrOverlap)*(i-1):120+(120-corrOverlap)*i]=a[120*i+corrOverlap:120*(i+1)]
     a[:120+7*(120-corrOverlap)]
   
   if flat!=None:a=XPADCorrFlatMask(a,flat,mask,fastmask)
   return a

def calcMaskMedian(d,nsigma=3,savemask=None):
   dm=d*0
   n=d.shape[0]
   for i in xrange(n):
      print "Median filtering to generate mask... %d/%d"%(i,n)
      dm[i]=signal.medfilt(d[i],3)
   wdiff=(d-dm)/sqrt(abs(dm)+1)
   themask=wdiff>nsigma
   themask=themask.sum(axis=0)>(n/5.+1)
   nbmask=(themask==0).sum()
   print "found %d pixels in mask (%6.4f%%)"%(nbmask,float(nbmask*100)/len(themask.flat))
   if savemask!=None:
      f=open(savemask,'w')
      cPickle.dump(themask,f,protocol=2)
      f.close()
   return themask

def calcMaskDark(savemask=None):
   # Import I0
   xpad_i0=importXPADI0()
   # find all dark images
   ll=os.listdir('./')
   dark_l=[]
   for k,v in xpad_i0.iteritems():
      i=ll.count(k)
      if i>0 and v['Seconds']>=1 and (v['vct1_3']/v['Seconds'])<10:
         i=ll.index(k)
         print "dark:",k
         dark_l.append(readXPAD2D(k))
   mask=dark_l[0]*0
   for d in dark_l:
      #imshow(d>2)
      mask+=(d>=2)
   mask=(mask<3) # exclude pixels which have had at least 3 times a count during a dark exposure
   nbmask=(mask==0).sum()
   print "found %d pixels in mask (%6.4f%%)"%(nbmask,float(nbmask*100)/len(mask.flat))
   if savemask!=None:
      f=open(savemask,'w')
      cPickle.dump(mask,f,protocol=2)
      f.close()
   return mask


def readXPAD3D(l,prefix=None,rebin=None,xpad_i0=None,corrBigPixels=True,dark=None,corrOverlap=5,flat=None,mask=None,fastmask=False):
   d=array([])
   n=len(l)
   ddark=0
   for i in xrange(len(l)):
      ix=l[i]
      i0=1
      if xpad_i0!=None:i0*=xpad_i0["%s%d.txt"%(prefix,ix)]['vct1_3']
      fname=prefix+"%d.dat"%ix
      tmp=readXPAD2D(fname,corrBigPixels=True,dark=dark,corrOverlap=corrOverlap,flat=flat,mask=mask,fastmask=fastmask)/i0
      if len(d)==0:
         nx,ny,nz=n,tmp.shape[0],tmp.shape[1]
         if rebin!=None:
            ny/=rebin[1]
            nz/=rebin[2]
         d=zeros((nx,ny,nz),dtype=float32)
      print ix,d.shape,tmp.min(),tmp.max()
      if rebin!=None: d[i,:,:]=rebin2(tmp,(rebin[1],rebin[2]))
      else :          d[i,:,:]=tmp
   if rebin!=None:
      if rebin[0]!=1:
         return rebin3(d,(rebin[0],1,1))
   return d

#! /usr/bin/python3
#
from __future__ import print_function
from numpy import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import *

hc=197.3269 # MeV fm
amu=931.494 #MeV
mp=938.27
mn=939.5653
mHe6=4*mn+2*mp-4.878*6
mHe4=2*mn+2*mp-7.073*4
mH3=2*mn+mp-2.8272*3

#
# Replace part of file strings
#
def replace_line_infile(ff, line_num, strings,ffo=None):
   """ replace entire lines with new one in file 
       line number is counted from zero !
   """
   f=open(ff,'r')
   lines=f.readlines()
   if (strings[-1]=='\n'):
     lines[line_num]=strings
   else:
     lines[line_num]=strings+'\n' # new line
   f.close()
   if not ffo:  # ffo is empty
     ffo=ff
   f=open(ffo,'w')
   for i in lines:
     f.write(i)
   f.close()

def replace_line_column_infile(ff,line_num,col_num,strings,ffo=None):
  """ replace column in line with new one in file
      line number and column number counts from zero"""
  f=open(ff,'r')
  lines=f.readlines()
  words=lines[line_num].split()
  words[col_num]=strings
  # new replaced line
  ss=''
  for i in words:
    ss=ss+i+'  '  
  lines[line_num]=ss+'\n'
  f.close()
  if not ffo: # ffo is empty
    ffo=ff
  f=open(ffo,'w')
  for i in lines:
    f.write(i)
  f.close()

def tabulate(inx,iny,xnew,fill_value=None):
  """
    convert array (x,y) into a regular table of
    [xnew,ynew] in equal steps, 
    fill_value= (y_min, y_max) : set y=y_min if x< inx_min
                                     y=y_max if x> inx_max
                'extrapolate'  : use extrapolation if x is out of range of inx         
    return [xnew,ynew] array                    
  """
  if fill_value:
    f=interp1d(inx,iny,bounds_error=False,fill_value=fill_value)
  else :
    f=interp1d(inx,iny)

  out=np.array([ [i,f(i)] for i in xnew])
  return out

def convert_table(fname,x_range,fill_value=None):
  """
   load data from file and convert it into a table in equal step size
   f_name : file name string
   x_range: array [xini,xfinal,xstep ]
   fill_value: option for extrapolation (y< value, y> value) or 'extrapolate'
  """
  xnew=np.arange(x_range[0],x_range[1],x_range[2])
  dat=np.loadtxt(fname)
  out=tabulate(dat[:,0],dat[:,1],x_range,fill_value=fill_value)
  return out 

def fresco_input_formfactor(dat,x_range,fname,fill_value=None,comment=''):
  """
    Generate input potential form factor format 
    converted from 'fname' file
  """
  xnew=np.arange(x_range[0],x_range[1],x_range[2])
  npoints=len(xnew) 
  rstep=x_range[2]
  rfirst=xnew[0]
  out=tabulate(dat[:,0],dat[:,1],x_range,fill_value=fill_value)  
  ff=open(fname,'w')
  ff.write(comment+'\n')
  ff.write('%i  %f  %f\n'%(npoints,rstep,rfirst))
  for i in out[:,1]:
    ff.write('%f\n'%(i))

  return  

def to_lab_energy(mp,mt,ecm):
   """ convert ecm to elab
      for projectile mass mp
          target mass mt
      all are in MeV
   """
   mu=mp*mt/(mp+mt*1.0)
   return mp/mu*ecm

def to_cm_energy(mp,mt,elab):
   mu=mp*mt/(mp+mt*1.0)
   return mu/mp*elab


def BW_form(x,er,gamma):
   """ Breit-Wigner form for resonance
   """
   return 0.25*gamma**2/((x-er)**2+0.25*gamma**2)


def find_res(ecm,phase,guess=1.0,guess2=2.0):
  """ find the position of resonance and width
    from energy and phase array (degree)
    Assume the case with Breit-Wigner form 
    background phase shift=0.
    units are degree.
    center of mass energy in MeV

    resonance energy is defined by passing 90 degree phase shift
    width is guessed from 45 degree phase shift
    But, then cross section is fitted with BW-form  

  """
  # first make the phase shift positive, continuous 
  for i in range(len(phase)) :
    if  phase[i]< 0 : 
       phase[i] = phase[i]+180. 
  # now find the position of 90 degree
  f=interp1d(ecm,phase-90.,'cubic')
  # brentq method for solution
  #er=brentq(f,guess,guess2)  
  #  root finding
  sol=root(f,guess) 
  er=sol.x  # resonance energy
  # now search width 
  #  first guess width from result
  f=interp1d(ecm,phase-45.,'cubic')
  sol=root(f,er)
  gam=abs(sol.x-er)/2.
  # now curve fit with BW_form
  f=interp1d(ecm,phase,'cubic')
  xn=arange(er-gam,er+gam,2*gam/100.)
  ydata=sin(f(xn)*pi/180.)**2
  popt,pcov=curve_fit(BW_form,xn,ydata,p0=(er,gam))
  
  xn=arange(ecm[0],ecm[-1],(ecm[-1]-ecm[0])/200.)
  ydata=sin(f(xn)*pi/180.)**2
  gdata=BW_form(xn,popt[0],popt[1])
  plt.plot(xn,ydata,xn,gdata)
  plt.title('E=%f +i %f'%(popt[0],popt[1]))
  plt.savefig('phase.png')
  return popt[0],popt[1] 

def clean_comm(fname):
  # remove all comments @ from the file
  f=open(fname,'r')
  lines=f.readlines()
  f.close() 
  f=open(fname+'x','w')
  for i in lines:
     if '@' in i:
        i='#'+i
     if ('END'in i) or ( '&' in i ):
        i='#'+i+'\n\n'
     f.write(i)
  f.close()
  return 

def chck_fresco_out(fname=''):
  """ test whether the FRESCO ended normally
      by checking 
      'Total CPU '
      and 
      'Recommended RNL: non-local width' as 'OK'     
  """
  ff=open(fname,'r')
  ll=ff.readlines()
  ff.close()
  chck=0
  for i in ll[:-10:-1]: #last 9 lines
    if 'Total CPU ' in i:
       chck=1
    if (chck==1) and ('Recommended RNL: non-local width' in i):
       if ': OK' in i:
         chck=2
  if chck==0 :
    print('ERROR in %s file'%fname)
    return 1
  elif chck==1 :
    print('ERROR in %s file'%fname)
    return 1
  elif chck==2 :  
    print('ok')
    return 0

def read_fresco_res(fname):
  """ read fresco results *.res files
      assume different data are separated by a blank line
      
      return dictionary 
  """
  ff=open(fname,'r')
  lines=ff.readlines()
  ff.close()
  out={}
  j=0
  ll=[]
  for i in lines:
    w=i.split()
    if len(w)!=0 and (w[0][0] in ['#','!']):
      continue
    if len(w)!=0 :
      ll.append([ float(k) for k in w])
    if len(w)==0 and len(ll)==0 : # continuos blank
      continue
    if (len(w)==0 and len(ll)!=0) or (i ==lines[-1]): #if met blank, it means end of data
      out[j]=ll[:]
      j=j+1
      ll=[]
  return out 

#=====================MAIN===================================== 

""" This routine only works for a specific cutoff value """
if __name__ == '__main__':
   import matplotlib.pyplot as plt
   import numpy as np
   import sys 
   from subprocess import call
  # remove all comments @ from the file
   if (sys.argv[1]=='res'):
     clean_comm(str(sys.argv[2]))
     fname=sys.argv[2]+'x'
     out=np.loadtxt(fname)
     en=out[:,0];phase=out[:,3]
     ecm=to_cm_energy(2,4,en)
     er,gam=find_res(ecm,phase,guess=0.8,guess2=2.0)
   elif (sys.argv[1]=='clean'):
     clean_comm(str(sys.argv[2]))
   else :
   #---default case
   # usage: test.py [infile] [inf]
   #        convert fort.16 into .res file
     infile=str(sys.argv[1])
     print(infile)
     inf = infile.split('.')[0]
     if (len(sys.argv)==3):
        inf = str(sys.argv[2])
     call("fresco < "+infile+" > "+inf+".out",shell=True)
     clean_comm('fort.16')
     call('mv fort.16x '+inf+'.res',shell=True)
     call('rm fort.*',shell=True)
     # call fresco
     # cleanup file
     # rename file 

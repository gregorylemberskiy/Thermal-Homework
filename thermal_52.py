#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as mlab
np.set_printoptions(precision=3,threshold='nan')
from scipy.optimize import fsolve
from scipy.integrate import quad

"""
Generates a plot of vant der waals isotherm, with appropriate maxwell construction.  
"""

def vanderwaals(v):
     t = .95
     return (8.*t)/(3.*v-1)-3./(v**2)

def Eq2(pars):
     v1, v3 = pars 
     print pars
     t = .95
     a = vanderwaals(v1) - vanderwaals(v3)
     b = quad(vanderwaals,v1,v3)[0] - vanderwaals(v1)*(v3-v1)
     return [a,b]

def main():
     ## Plotting Figures
     fig1   = plt.figure(figsize=(11,8))
     ax1    = fig1.add_subplot(111)
#     ax2    = fig1.add_subplot(111)

     v = np.arange(.333, 10, .001)
     p_iso = vanderwaals(v)          
     
     pars = [.5, 5.]
     vopt = fsolve(Eq2, pars)

     pv1 = vanderwaals(vopt[0])
     pv2 = vanderwaals(vopt[1])
     print "pv1 and pv2", pv1,pv2
     mxwllcnstrg = np.arange(vopt[0],vopt[1],.001)

     ax1.set_title('Van Der Waals',fontsize=24)
     ax1.plot(v,p_iso)
     ax1.plot(mxwllcnstrg,len(mxwllcnstrg)*[pv1])
     ax1.axvline(vopt[0],linestyle = 'dashed')
     ax1.axvline(vopt[1],linestyle = 'dashed')
     ax1.set_xlabel('$v$',fontsize=24)
     ax1.set_ylabel('$p$',fontsize=24)
     ax1.grid(True)
     ax1.set_xlim(.33333,5)
     ax1.set_ylim(0,2)
     plt.show()

if __name__ == "__main__":
     main()

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

def vanderwaals(v,t):
     return (8.*t)/(3.*v-1)-3./(v**2)

def Equations(pars,t):
     v1, v2 = pars 
     print pars
     a = vanderwaals(v1,t) - vanderwaals(v2,t)
     b = quad(vanderwaals,v1,v2,args=(t))[0] - vanderwaals(v1,t)*(v2-v1)
     return [a,b]

def G(v,t): 
     return 8.*t/(3*(3.*v-1.))-6./v-(8.*t/3.)*np.log(3.*v-1) 
#     return t/(3.*v-1.)-9./(4.*v)-t*np.log(3.*v-1)  # Kelvin's Answer
 
def main():
     ## Plotting Figures
     fig1   = plt.figure(figsize=(11,8))
     ax1    = fig1.add_subplot(121)
     ax2    = fig1.add_subplot(122)

     v = np.arange(.333, 10, .001)
     t = .95
     p_iso = vanderwaals(v,t)          
     
     pars = [.5, 5.]
     vopt = fsolve(Equations, pars,args=(t))

     pv1 = vanderwaals(vopt[0],t)
     pv2 = vanderwaals(vopt[1],t)
     print "pv1 and pv2", pv1, pv2

     mxwllcnstrg = np.arange(vopt[0], vopt[1], .001)
     gibbs = G(v,.95)

     ax1.set_title('Van Der Waals',fontsize=24)
     ax1.plot(v,p_iso)
     ax1.plot(mxwllcnstrg,len(mxwllcnstrg)*[pv1])
     ax1.axvline(vopt[0],linestyle = 'dashed')
     ax1.text(vopt[0],.5," v1 = " + str('%.3f' %vopt[0]))
     ax1.text(vopt[1],.5," v2 = " + str('%.3f' %vopt[1]))
     ax1.text(1.,1., "Vapor Pressure = " + str('%.3f' % pv1))

     ax1.axvline(vopt[1],linestyle = 'dashed')
     ax1.set_xlabel('$v$',fontsize=24)
     ax1.set_ylabel('$p$',fontsize=24)
     ax1.grid(True)
     ax1.set_xlim(.33333,3)
     ax1.set_ylim(0,1.5)

     ax2.set_title('Gibbs Free Energy ',fontsize=24)
     ax2.plot(p_iso,gibbs)
     ax2.axvline(pv1,linestyle = 'dashed')

     ax2.set_xlabel('$P/Pc$',fontsize=24)
     ax2.set_ylabel('$G$',fontsize=24)
     ax2.grid(True)
     ax2.set_xlim(.7,.9)
     ax2.set_ylim(-6.8,-6.2)
#     ax2.set_ylim(-2.50,-2.4)
#     plt.savefig("Thermal_52.png")
     plt.show()

if __name__ == "__main__":
     main()

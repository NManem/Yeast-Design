# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:19:53 2020

@author: Nitesh Manem
"""
#Constants
a1 = 6.67
a3 = 2.22
a2 = 2.20
k1 = 5
k2 = 0.38
KE = 0.02
KS = 0.5
K = 0.51

def ra(Aw, S):
    return k1*Aw*S/(KS+S)

def rb(Aw, E):
    return k2*Aw*E/(KE+E)

def rc(Bw):
    return K*Bw

def dBwdt(t, Aw, Bw, S, E):
    return 2*ra(Aw, S) + 2*rb(Aw,E)-rc(Bw)

def dAwdt(t, Aw, Bw, S, E):
    return -1*ra(Aw,S)-rb(Aw,E)+rc(Bw)

def dEdt(t, Aw, Bw, S, E):
    return a2*ra(Aw, S)-a3*rb(Aw, E)

def dSdt(t, Aw, Bw, S, E):
    return -1*a1*ra(Aw,S)
    

# Finds value of y for a given x using step size h 
# and initial value y0 at x0. 


def rungeKutta(Aw0, Bw0, E0, S0, t0, t, h): 
    AwList = []
    BwList = []
    EList = []
    SList = []
    tList = []
    
    # Count number of iterations using step size or 
    # step height h 
    n = (int)((t - t0)/h)  
    
    AwList.append(Aw0)
    BwList.append(Bw0)
    EList.append(E0)
    SList.append(S0)
    tList.append(t0)
    
    # Iterate for number of iterations 
    Aw = Aw0
    Bw = Bw0
    E = E0
    S = S0
    for i in range(1, n + 1): 
        "Apply Runge Kutta Formulas to find next value"
        k1A = h * dAwdt(t0, Aw, Bw, S, E)
        k1B = h * dBwdt(t0, Aw, Bw, S, E)
        k1E = h * dEdt(t0, Aw, Bw, S, E)
        k1S = h * dSdt(t0, Aw, Bw, S, E)
        
        k2A = h * dAwdt(t0 + 0.5*h, Aw+k1A/2, Bw+k1B/2, S+k1S/2, E+k1E/2)
        k2B = h * dBwdt(t0 + 0.5*h, Aw+k1A/2, Bw+k1B/2, S+k1S/2, E+k1E/2)
        k2E = h * dEdt(t0 + 0.5*h, Aw+k1A/2, Bw+k1B/2, S+k1S/2, E+k1E/2)
        k2S = h * dSdt(t0 + 0.5*h, Aw+k1A/2, Bw+k1B/2, S+k1S/2, E+k1E/2)
        
        k3A = h * dAwdt(t0 + 0.5*h, Aw+k2A/2, Bw+k2B/2, S+k2S/2, E+k2E/2)
        k3B = h * dBwdt(t0 + 0.5*h, Aw+k2A/2, Bw+k2B/2, S+k2S/2, E+k2E/2)
        k3E = h * dEdt(t0 + 0.5*h, Aw+k2A/2, Bw+k2B/2, S+k2S/2, E+k2E/2)
        k3S = h * dSdt(t0 + 0.5*h, Aw+k2A/2, Bw+k2B/2, S+k2S/2, E+k2E/2) 
        
        k4A = h * dAwdt(t0 + h, Aw+k3A, Bw+k3B, S+k3S, E+k3E)
        k4B = h * dBwdt(t0 + h, Aw+k3A, Bw+k3B, S+k3S, E+k3E)
        k4E = h * dEdt(t0 + h, Aw+k3A, Bw+k3B, S+k3S, E+k3E)
        k4S = h * dSdt(t0 + h, Aw+k3A, Bw+k3B, S+k3S, E+k3E) 
         
        # Update next values 
        Aw = Aw + (1.0 / 6.0)*(k1A + 2 * k2A + 2 * k3A + k4A)
        Bw = Bw + (1.0 / 6.0)*(k1B + 2 * k2B + 2 * k3B + k4B)
        E = E + (1.0 / 6.0)*(k1E + 2 * k2E + 2 * k3E + k4E)
        S = S + (1.0 / 6.0)*(k1S + 2 * k2S + 2 * k3S + k4S)
        
        
        # Update next value of t 
        t0 = t0 + h 
        
        AwList.append(Aw)
        BwList.append(Bw)
        EList.append(E)
        SList.append(S)
        tList.append(t0)
        
    return [tList, AwList, BwList, SList, EList] 
  
# Driver method 
    
import matplotlib.pyplot as plt
import numpy as np

Aw0 = 0.005
Bw0 = 0.005
E0 = 0.01
S0 = 10000


t0 = 0.0001
t = 40
h = 0.00001



tList, AwList, BwList, SList, EList = rungeKutta(Aw0, Bw0, E0, S0, t0, t, h)

XList = np.add(AwList,BwList)

print(len(XList))
print(XList[-1])

fig, ax1 = plt.subplots()
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Ethanol and Glucose Concentration (g/L)')
ax1.scatter(tList, EList, color = 'r', marker = '.', label = "Ethanol")
ax1.scatter(tList, SList, color = 'g', marker = '.', label = 'Glucose')

ax2 = ax1.twinx()
ax2.set_ylabel('Biomass Concentration (g/L)')
ax2.scatter(tList, XList, color = 'b', marker = '.', label = 'Biomass')

fig.tight_layout()
ax1.legend(loc=6)
ax2.legend(loc=7)
plt.show()

print('Max value of ethanol', max(EList))

print('Done!')
  
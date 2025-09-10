'''
Tyrosine pKa1/carboxyl=2.20,pKa2/ammonium=9.11,pKa3/phenol=10.07,pI=5.66
use conservation of charge and mass to calculate concetration
set:"concentration of H+"=H, "concentration of Tyrosine zwitterion"=x
Newtonian (fsolve) is applied to solve equations
parameters:0.1M*10 mL Tyr-HCl, 0.1M NaOH as titrant
aim:export a theoretical titration curve of 0.1M Tyr-HCl solution
    mark equivalence points 
    mark pH at equivalence points and 'best buffer points'
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#data
K_CARBOXYL=10**(-2.20)
K_AMMONIUM=10**(-9.11)
K_PHENOL=10**(-10.07)
K_W=10**(-14)
C_INIT_TYROSINE=0.1 #unit:mol/l
C_BASE=0.1 
V_INIT_TYROSINE=10 #unit:ml
V_BASE_LIST=[]
for i in range (0,801):
    V_BASE_LIST.append(i/20)
PH_LIST=[]

#define function for solving equation
def chargebalance(ph,v_b):#used in fsolve,v_b=volume of base added
    H=10**(-ph)
    V_TOTAL=V_INIT_TYROSINE+v_b #ignore volume change in mixing
    C_TYR_TOTAL=1/V_TOTAL #1 means 1mmol Tyr-HCl
    C_NA=C_BASE*v_b/V_TOTAL
    C_CL=C_INIT_TYROSINE*V_INIT_TYROSINE/V_TOTAL
    #chloride is a counter ion in initial TyrH+ solution. n(Cl)=1mmol
    x=C_TYR_TOTAL/(H/K_CARBOXYL+1+K_AMMONIUM/H+K_AMMONIUM*K_PHENOL/(H**2))
    #x=CONCENTRATION OF ZWITTERION,BASED ON DEFINITION OF Ka and mass balance
    positive=(H+(x*H/K_CARBOXYL)+C_NA)*V_TOTAL
    negative=(K_W/H+x*K_AMMONIUM/H+(2*x*K_AMMONIUM*K_PHENOL/H**2))*V_TOTAL+1
    #amounts of substances, unit:mmol
    return positive-negative

#iteration at every drop  
PH_0=-np.log10(np.sqrt(K_CARBOXYL*C_INIT_TYROSINE))
#use approx formula to obtain the first Newtonian starting point
for v_base in V_BASE_LIST:
    sol=fsolve(chargebalance,PH_0,args=(v_base,))
    PH=sol[0]
    PH_0=PH
    #use the pH value as a Newtonian starting point for next equation
    #because pH is continuous and monotonic
    PH_LIST.append(PH)


#plotting part
plt.figure()
plt.plot(V_BASE_LIST,PH_LIST,label='Titration curve',color='b')
plt.title('Theoretical titration curve of Tyrosine')
plt.xlim(0,40)
plt.ylim(0,14)
plt.xlabel('v_base(mL)')
plt.ylabel('pH')
#mark equivalence points where there is only one kind of Tyrosine
#these pH are calculated previously 
plt.plot(10,5.66,marker='o',label='1st equivalence pH=5.66')
plt.plot(20,9.59,marker='o',label='2nd equivalence pH=9.59')
plt.plot(30,11.20,marker='o',label='3rd equivalence pH=11.20')
#mark the pkas and pi
plt.axhline(y=5.66,linestyle='dotted',color='r')#pI
plt.text(32,5.70,'pI=5.66')
plt.axhline(y=2.20,linestyle='dotted',color='r')#pK1
plt.text(32,2.25,'pKa1=2.20')
plt.axhline(y=9.11,linestyle='dotted',color='r')#pK2
plt.text(32,9.15,'pKa2=9.11')
plt.axhline(y=10.07,linestyle='dotted',color='r')#pK3
plt.text(32,10.15,'pKa3=10.07')

plt.grid()
plt.legend()
plt.show()












            
            

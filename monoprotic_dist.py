#data'base' part
import pandas as pd
import math
import matplotlib.pyplot as plt
#monoprotic acid part
this_acid=input('Acid name:')
#slice a dataframe and generate a list [name,pka1]
pka=float(input('pKa:'))

k1=math.pow(10,-pka)
pH=[]#x-axis
dis0=[]
dis1=[]
#y-axis,dis0=HA,dis1=A anion
#calculation
for i in range(0,1401):
    ph=i/100
    pH.append(ph)
    ch=math.pow(10,-ph) #[H+]
    d0=ch/(ch+k1)
    dis0.append(d0)
    d1=k1/(ch+k1)
    dis1.append(d1)
#plot
plt.figure()
#tagged pH=pKa
xpoint=[pka]
ypoint=[0.5]
plt.plot(pH,dis0,label='HA')
plt.plot(pH,dis1,label='A$^-$')
plt.plot(xpoint,ypoint,marker='o')
plt.axvline(pka,linestyle='dotted',color='r')
plt.text(pka,0.5,'  ('+str(pka)+',0.5)')
plt.legend()
plt.grid()
plt.title('Distribution-pH diagram for '+this_acid)
plt.xlabel('pH')
plt.ylabel('Distribution coefficient delta')
plt.show()







































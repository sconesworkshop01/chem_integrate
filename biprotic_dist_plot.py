#data'base' part
import pandas as pd
import math
import matplotlib.pyplot as plt
#all conjugated acids have 'protonated' in the name
#H2S pKa2 from Pubchem
#diprotic acid part
acid_name=input('Acid name:')
#slice a dataframe and generate a list [name,pka1,pka2]
pka1=float(input('pKa1='))
pka2=float(input('pKa2='))

k1=math.pow(10,-pka1)
k2=math.pow(10,-pka2)
pH=[]#x-axis
h2a=[]
ha=[]
a=[]
#y-axis,di0=[H2A],di1=[HA-] anion,di2=[A2-]
#calculation
for i in range(0,1401):
    ph=i/100
    pH.append(ph)
    ch=math.pow(10,-ph) #[H+]
    
    di0=math.pow(ch,2)/(math.pow(ch,2)+k1*ch+k1*k2)
    h2a.append(di0)
    di1=(ch*k1)/(math.pow(ch,2)+k1*ch+k1*k2)
    ha.append(di1)
    di2=(k1*k2)/(math.pow(ch,2)+k1*ch+k1*k2)
    a.append(di2)
#plot
plt.figure()
#tagged pH=pKa
xpoint=[pka1,pka2]
ypoint=[0.5,0.5]
plt.plot(pH,h2a,label='H$_2$A')
plt.plot(pH,ha,label='HA$^{-}$')
plt.plot(pH,a,label='A$^{2-}$')
#pH=pKa labels
plt.plot(pka1,0.5,marker='o')
plt.plot(pka2,0.5,marker='o')
plt.axvline(pka1,color='r',linestyle='dotted')
plt.axvline(pka1,color='r',linestyle='dotted')
plt.text(pka1-2.5,0.5,'  ('+str(pka1)+',0.5)')
plt.text(pka2,0.5,'  ('+str(pka2)+',0.5)')

#isoelectric points
plt.axvline((pka1+pka2)/2,color='g',linestyle='--')
plt.text((pka1+pka2)/2,0.8,'  pI='+str((pka1+pka2)/2))

plt.legend()
plt.grid()
plt.title('Distribution-pH diagram for '+acid_name)
plt.xlabel('pH')
plt.ylabel('Distribution coefficient delta')
plt.show()







































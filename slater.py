'''slater rule calculation
group(1s)(2s,2p)(3s,3p)(3d)(4s,4p)(4d)(4f)(5s,5p)(5d)(5f)(6s,6p)
num_e:2   8      8      10  8      10  14  8      10  14  8
def shielding constant:
group=input
An amount of 0.35 /e within the same group
the groups outside group does not shield 
if group=='1s':
    0.30/e
elif group in [ns, np](n>=2):
    an amount of 0.85/e with principal quantum number (n–1),
    an amount of 1.00/e with principal quantum number (n–2) or less.
if group in [d] or group in [f]:
    an amount of 1.00/e for other groups smaller
    #This includes both electrons with a smaller principal quantum number than n
    and "n and a smaller l"
    
aim:input atomic number and orbital name,output shielding constant and Z*
'''
atomic_num=int(input('Atomic number:'))#this is an int var. e.g.12,34,56
orbital=input('orbital:')#This is a str var with len=2. e.g. '1s','2p'
slater_list=[['1s'],['2s','2p'],['3s','3p'],['3d'],
             ['4s','4p'],['4d'],['4f'],['5s','5p'],['5d'],['5f'],
             ['6s','6p'],['6d'],['7s','7p']]#only 118 elements!

def get_config(Z):#Z is int
    config=[]
    remaining=Z
    n=1
    l_str='spdfghij'#index=l
    #aufbau principle:smallest"n+l"fills first,
    #if n+l is the same,small n first(e.g. 3p-4s-"3d-4p-5s")
    for A in range(1,8):
        for n in range(1,A+1):
            l=A-n
            if l>=n:#angular momentum QM<=(n-1)
                continue
            else:
                if remaining<=0:
                    break
                max_e=4*l+2
                orbital_name=str(n)+l_str[l]
                real_e=min(remaining,max_e)
                remaining-=real_e
                config.append([orbital_name,real_e])
            if remaining<=0:
                break
            n+=1
    cfg = dict(config)
    #consider Hund principle
    for n in (3, 4, 5, 6):
        d, s = f'{n}d', f'{n+1}s'
        if cfg.get(d, 0) == 4 and cfg.get(s, 0) >= 1:
            cfg[d] += 1
            cfg[s] -= 1
        elif cfg.get(d, 0) == 9 and cfg.get(s, 0) >= 1:
            cfg[d] += 1
            cfg[s] -= 1

    if Z == 41:          # Nb
        cfg['4d'], cfg['5s'] = 4, 1
    elif Z == 46:        # Pd
        cfg['4d'], cfg['5s'] = 10, 0
    elif Z == 78:        # Pt
        cfg['5d'], cfg['6s'] = 9, 1

    def sort_key(item):
        name = item[0]#name of an orbital
        n_value = int(name[0])
        l_value={'s':0,'p':1,'d':2,'f':3}.get(name[1],4)
        return (n_value, l_value)
    config = [[orb, e] for orb, e in cfg.items() if e > 0]
    #turn the dict into list for output
    config.sort(key=lambda x: (int(x[0][0]), 'spdf'.index(x[0][1])))
    return config
c=get_config(atomic_num)
print("Electron configuration is:",c)

def slater(Z,orb_id):#Z:int,orb_id:str,output shielding constant
    target_n = int(orb_id[0])          
    target_l = orb_id[1]  
    this_cfg=dict(c)
    for j in range(0,len(slater_list)):
        if orb_id in slater_list[j]:
            group_no=j #find the index of the highest group
    layer_pop = [0] * len(slater_list)
    for orb, e in c:
        idx = next(i for i, g in enumerate(slater_list) if orb in g)
        layer_pop[idx] += e

    sigma = 0.0
    for idx in range(len(slater_list)):#layer no. in slater_list
        this_n = int(slater_list[idx][0][0])   #principle QM
        if idx == group_no:     #same layer               
            cnt = layer_pop[idx] - 1
            sigma += (0.30 if orb_id == '1s' else 0.35) * cnt
        elif idx > group_no: #no_shield
            continue
        else: #inner layers
            if target_l in 'sp':               # s/p
                if this_n == target_n - 1:     # (n-1) 
                    sigma += 0.85 * layer_pop[idx]
                elif this_n <= target_n - 2:   # <=(n-2) 
                    sigma += 1.00 * layer_pop[idx]
            else:                              # d/f
                sigma += 1.00 * layer_pop[idx]
    return sigma
slater=slater(atomic_num,orbital)
print('The shielding coefficient is:',slater)
print("The effective atomic number is:",round(atomic_num-slater,2))
    
    

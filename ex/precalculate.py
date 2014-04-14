import statePlanar as sp
import itertools
import sys
import numpy as np
import copy
import os

d=5


syndrome_file = "%s/precalculatedDecoder/%dsyndromes_xyz.txt"%(os.environ['HOME'],d)
results_file='%s/precalculatedDecoder/size%d_xyz.txt'%(os.environ['HOME'],d)

L=sp.PlanarLattice(d)
#L.add_errors(0.5)
#L.show()

stabilisers= L.x_stabiliser_indices + L.z_stabiliser_indices

qubits = L.qubit_indices
'''
syndrome_list=[]
# loop over all possible Y-only error configurations
for x in itertools.product([0,1,2,3],repeat =len(qubits)):
    L.array=np.zeros((d,
                      d),dtype='uint8')
    L.array[zip(*L.qubit_indices)]=np.array(x)
    syndrome = L.get_syndrome()
    if syndrome not in syndrome_list:
        syndrome_list+=[syndrome]

'''
syndrome_list=[]
for n in range(len(stabilisers)+1):
    for syndrome in itertools.combinations(stabilisers,n):
        syndrome_list+=[syndrome]




test=0
results=[]
p=0.1
plists=[]

f=open(syndrome_file,"w")

for syndrome in syndrome_list:
    print syndrome
    m_init = L.matching_from_syndrome(syndrome)
    
    m_0=copy.copy(m_init)
    L.array=m_init
    L.add_z()
    m_1= copy.copy(L.array)
    L.add_x()
    m_3=copy.copy(L.array)
    L.add_z()
    m_2=copy.copy(L.array)
    
    plist=[]
    e_list=[]
    for m in [m_0,m_1,m_2,m_3]:
        
        error_freq=[0 for i in range(len(L.qubit_indices)+1)]

        
        for n in range(len(stabilisers)+1):
            for stabs in itertools.combinations(stabilisers,n):
    
                L.array=copy.copy(m)
                for (x,y) in stabs:
                    L.apply_stabiliser(x,y)
              
             
                
                error_count=0
                er_test=False
                
                for er in [x for y in L.array.tolist() for x in y]:
                    if er!=0: error_count+=1
                

                    
                  #  if er==3 or er==2:
                   #     er_test=True
                    #    break
                    #elif er==1:
                     #   error_count+=1
                      
                
                #if er_test==True: continue       
                error_freq[error_count]+=1
                
                
        sumv=0
        for i in range(len(L.qubit_indices)+1):
            sumv+=((1-p)**5)*error_freq[i]*((p/(1-p))**i)
        
        plist+=[sumv]
        e_list+=[ error_freq]

    results+=[e_list]
    for logi in e_list:
       
        for el in logi:
            f.write("%d "%(el,))
        f.write(",")
    f.write("\n")
f.close()
print results


'''
p=0.1
pre=(1-p)**(len(L.qubit_indices))
r=p/((1-p))
print len(L.qubit_indices)
maxvals=[]
for syn in results: #    
    val=[]
    for log_e in syn:
        sumv=0

        for i in range(len(L.qubit_indices)+1):
            sumv+=log_e[i]*(r**i)

        #val+=[sum([log_e[i]*(r**i) for i in range(len(L.qubit_indices)+1)])]
        val+=[sumv]

    maxvals+=[pre*max(val)]
    
print maxvals,sum(maxvals)



'''


data = []
for p in [0.01*i for i in range(100)]:
     
    pre=(1-p)**(len(L.qubit_indices))
    r=p/(3*(1-p))

    maxvals=[]
    for syn in results: #
    
        val=[]
        for log in syn:
       
            val+=[sum([log[i]*(r**i) for i in range(len(L.qubit_indices)+1)])]
        maxvals+=[pre*max(val)]

    
    data+=[[p, max(p/3,(1-p)),sum(maxvals)]]


f=open(results_file,'w')
for a,b,c in data:
    f.write("%f %f\n"%(a,c))
f.close()
    
#import matplotlib.pyplot as plt

#plt.plot(zip(*data)[0],zip(*data)[1],zip(*data)[0],zip(*data)[2])
#plt.show()
            



    
'''
#loops over possible error configurations
for n in range(len(qubits)+1):
    for error_config in itertools.combinations(qubits,n):
        break
        

# loops over all possible syndromes
for n in range(1,len(stabilisers)+1):
    for syndrome in itertools.combinations(stabilisers,n):

        break;

 '''       
        
#L.match_from_syndrome(syndrome)

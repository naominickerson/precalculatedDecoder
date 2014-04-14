## Created by Naomi Nickerson
##
## This is code to
##
## Qubits are encoded as a number in {0,1,2,3}
## 0 = no error
## 1 = Z error
## 2 = X error
## 3 = Y error

import statePlanar as sp
import itertools
import sys
import numpy as np
import copy

## LATTICE SIZE
d=3
size = 3

## ERROR TYPE
error_type = "y"

## LIES
lies = True

######################

if error_type not in {"y","xyz"}:
    print "the error type ",error_type," isn't supported"
    sys.exit(0)    


lie_type = "_lies" if lies==True else ""

file1='%dsyndromes_%s%s.txt'%(size,error_type,lie_type)
file2='%dlies_%s%s.txt'%(size,error_type,lie_type)
file3='%dresults_%s%s.txt'%(size,error_type,lie_type)



L=sp.PlanarLattice(d)
stabilisers= L.x_stabiliser_indices + L.z_stabiliser_indices
qubits = L.qubit_indices


## THIS CREATES A RESTRICTED LIST OF POSSIBLE SYNDROMES, ONLY THOSE
## POSSIBLE GIVEN THE PARTICULAR ERROR MODEL


syndrome_list=[]

;2if(error_type == "y" and lies == False):
    # loop over all possible Y-only error configurations (w/o lies)
    for x in itertools.product([0,1,2,3],repeat =len(qubits)):
        L.array=np.zeros((d,d),dtype='uint8')
        L.array[zip(*L.qubit_indices)]=np.array(x)
        syndrome = L.get_syndrome()
        if syndrome not in syndrome_list:
            syndrome_list+=[syndrome]
else:
## THIS CREATES A LIST OF ALL POSSIBLE SYNDROMES (for when lies are included)

    for n in range(len(stabilisers)+1):
        for syndrome in itertools.combinations(stabilisers,n):
            syndrome_list+=[syndrome]
    syndrome_sets=[set(s) for s in syndrome_list]



test=0
results=[]
p=0.1
plists=[]

f=open(file1,"w")

# FOR EVERY POSSIBLE SYNDROME
for syndrome in syndrome_list:

    ## Create initial possible matching
    m_init = L.matching_from_syndrome(syndrome)
    L.array=m_init

    ## CHECK LOGICAL STATE HERE AND SET IT TO ZERO
    if L.measure_z(m_init): L.add_z()
    if L.measure_x(m_init): L.add_x()
    m_init=copy.copy(L.array)

    
    ## CREATE THE FOUR DIFFERENT LOGICAL STATES m_0,1,2,3    
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

        ## initialise a list to count errors     
        error_freq=[0 for i in range(len(L.qubit_indices)+1)]

        ## loop over all possible stabiliser combinations (of all lengths n)    
        for n in range(len(stabilisers)+1):
            for stabs in itertools.combinations(stabilisers,n):
    
                L.array=copy.copy(m)
                for (x,y) in stabs:
                    L.apply_stabiliser(x,y)
                          
                error_count=0
                er_test=False

                ## for every qubit error state in the updated array 
                for er in [x for y in L.array.tolist() for x in y]:

                    ## VERSION 1: Full channel noise                    
                    if error_type == "xyz":
                        if er!=0: error_count+=1
                    

                    ## VERSION 2: Single channel noise
                    ## To restrict to Y errors only if the error type is
                    ## 1 or 2 (X or Z) for any
                    ## qubit we break out of the loop and skip to the
                    ## next stabilizer configuration. 
                    if error_type =="y":
                        if er==1 or er==2:
                            er_test=True
                            break                    
                        elif er==3:
                            error_count+=1
                
                if er_test==True: continue       
                error_freq[error_count]+=1
                

        ## we now have a list 'error_freq' which has counted all the different
        ## possible errors for this particular logical error.
        e_list+=[error_freq]

        #=================================================#
        ## This additional section calculates the success rate
        ## for specific error rates specified above
        #sumv=0
        #for i in range(len(L.qubit_indices)+1):
        #    sumv+=((1-p)**5)*error_freq[i]*((p/(1-p))**i)      
        #plist+=[sumv]


        #=================================================#

    # accumulate all the error counting results (one for each possible syndrome)
    results+=[e_list]


    #===========================================#
    ## write the syndrome results to file1
    for logi in e_list:
        for el in logi:
            f.write("%d "%(el,))
        f.write(",")
    f.write("\n")
f.close()




## GENERATE 'slies_list' which for each observed syndrome
## contains the number of lies needed if each other
## syndrome was the 'real' one
if lies == True:
    slies_list=[]
    for s_observed in syndrome_sets:
        slies=[]
        for s_real in syndrome_list:
             slies+=[len(s_observed.symmetric_difference(s_real))]
        slies_list+=[slies]

    f=open(file2,"w")

    for x in slies_list:
        for y in x:
            f.write("%d "%(y,))
        f.write("\n")
    f.close()



#===================================================================#
#===================End of calculations=============================#
#===================================================================#



### results: list of error configs compatible with each syndrome
### list of lies between pairs of syndromes



if lies == True:

    output=[]
    prange=[0.001*i for i in range(70)]
    qrange=[0.005*i for i in range(40)]


    for p in prange:

        if error_type == "xyz":
            r=(p/(3*(1-p)))
        elif error_type == "y":
            r=(p/((1-p)))
       

        pre_p=(1-p)**len(L.qubit_indices)

        syndrome_probs=[]
        for j in range(len(results)):
            val=[]  
            for logical in results[j]:
           
                val+=[pre_p*sum([logical[k]*(r**k) for k in range(len(L.qubit_indices)+1)])]
               
            syndrome_probs+=[np.array(val)]

 
        for q in qrange:
            

            qr=q/(1-q)
            pre_q=(1-q)**len(stabilisers)

            #for observed syndrome
            results_array2=[]
            for i in range(len(syndrome_probs)):
                #for possible real syndrome
                results_array=np.zeros(4)
                for j in range(len(syndrome_probs)):
                    results_array+=(qr**(slies_list[i][j]))*syndrome_probs[j]
                results_array2+=[pre_q*max(results_array)]
            #print results_array2
            
            #print pre_q*results_array
        
            output+=[[q,p,sum(results_array2)-max(p,1-p)]]
             
    f=open(file3,'w')
    for a,b,c in output:
        f.write("%f %f %f\n"%(a,b,c))
    f.close()


    ### Save results




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

## Calculate performance data for PERFECT MEASUREMENTS 
if lies == False: 

    data = []
    for p in [0.01*i for i in range(100)]:
         
        pre=(1-p)**(len(L.qubit_indices))

        if error_type ==  "y":
            r=p/((1-p))
        else:
            r = p/(3*(1-p)) #for xyz errors 

        maxvals=[]
        for syn in results: #
        
            val=[]
            for log in syn:
           
                val+=[sum([log[i]*(r**i) for i in range(len(L.qubit_indices)+1)])]
            maxvals+=[pre*max(val)]

        data+=[[p, max(p,(1-p)),sum(maxvals)]]



    f=open(file3,'w')
    for a,b,c in data:
        f.write("%f %f\n"%(a,c))
    f.close()


    
#import matplotlib.pyplot as plt
#plt.plot(zip(*data)[0],zip(*data)[1],zip(*data)[0],zip(*data)[2])
#plt.show()

            

         
        

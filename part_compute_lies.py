
import statePlanar as sp
import itertools
import sys
import numpy as np
import copy
import os

#print sys.argv[1:]

try: 
    d,error_type,p,qend,nsteps,results_file = sys.argv[1:]
    d=int(d)
    p = float(p)
    nsteps = int(nsteps)
    qend = float(qend)

except: 
    print "Usage: [d] [error_type] [p] [qend] [nsteps] [results_file]"
    sys.exit(-1)


L=sp.PlanarLattice(d)
stabilisers= L.x_stabiliser_indices + L.z_stabiliser_indices
qubits = L.qubit_indices

home = os.environ['HOME']

f = open("%s/precalculatedDecoder/%dsyndromes_%s_lies.txt"%(home,d,error_type),"r");
results = [[[int(z) for z in y.split()]
        for y in x.split(",") if y!="\n"]
       for x in f.readlines()];
f.close();

f= open("%s/precalculatedDecoder/size%d_%s_lies.txt"%(home,d,error_type),"r")
slies_list = [[int(y) for y in x.split()] for x in f.readlines()]
f.close();



qmultiplier = qend/nsteps;

output=[]
prange=[p]
qrange=[qmultiplier*i for i in range(nsteps)]



#results_file = "%s/precalculatedDecoder/%doutput_%s_lies_test.txt"%(home,d,error_type)
f=open(results_file,'w')
f.write("results for %d-code with %s errors and lies\n"%(d,error_type))
f.close()

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
        
        res = sum(results_array2)
        newoutput=[q,p,res,res-max(p,1-p)]
        f=open(results_file,"a");
        f.write("%f %f %f\n"%(q,p,res))
        f.close()
        #output+=[newoutput]
        


#f.close()


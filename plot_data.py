import importResults
import matplotlib.pyplot as plt

'''
## PERFECT MEASUREMENTS, Y ERRORS
ax = plt.gca()
data1 = [0.01*max([x,100-x]) for x in range(100)]
data3 = importResults.results("results/3results_y.txt")
data5 = importResults.results("results/5results_y.txt")

plt.plot([0.01*x for x in range(100)],data1,label="1-code")
plt.plot(zip(*data3)[0],zip(*data3)[1],label = "3-code")
plt.plot(zip(*data5)[0],zip(*data5)[1],label = "5-code")

plt.legend(loc=1)
plt.ylabel("p_d")
plt.xlabel("p")
plt.ylim([0.4,1])

plt.savefig("planar_y.png")
#plt.show()
plt.close()


## PERFECT MEASUREMENTS, XYZ ERRORS

ax = plt.gca()

data1 = [0.01*max([x/3.,100-x]) for x in range(100)]
data3 = importResults.results("results/3results_xyz.txt")
data5 = importResults.results("results/5results_xyz.txt")

plt.plot([0.01*x for x in range(100)],data1,label = "1-code")
plt.plot(zip(*data3)[0],zip(*data3)[1],label ="3-code")
plt.plot(zip(*data5)[0],zip(*data5)[1],label="5-code")

plt.legend(loc=1)
plt.ylabel("p_d")
plt.xlabel("p")
plt.ylim([0.,1])

plt.savefig("planar_xyz.png")
#plt.show()
plt.close()
'''

## LIES 

from scipy.interpolate import Rbf
import numpy as np



filename = "5results_xyz_lies_2.txt"
data = importResults.results("results/%s"%(filename,))

### Can combine  multiple data files here: 
data = data+importResults.results("results/5results_xyz_lies.txt")

## pull out parameters from name 
size = int(filename[0])
error_type = filename.split('_')[1]

#filter data
filtereddata = []
for q,p,z in data:
    if q>0.01: continue
    if p>0.25: continue
    filtereddata+=[[q,p,z]]

data = filtereddata

qlist = list(zip(*data)[0])
qlist.sort()
qrange = []
for q in qlist:
    if q not in qrange: qrange.append(q)

plist = list(zip(*data)[1])
plist.sort()
prange=[]
for p in plist:
    if p not in prange:prange.append(p)


zarray = [[0 for _ in prange] for _ in qrange]
for q,p,z in data:
    zarray[qrange.index(q)][prange.index(p)]=z




zrelative = []
for q,p,z in data:
    zrelative+=[z-max(p,1-p)]


#print z
Q,P=np.meshgrid(qrange,prange)

q,p,z=zip(*data)
rbf=Rbf(q,p,zrelative,function='linear')
Z=rbf(Q,P)






#matplotlib.rcParams['contour.negative_linestyle']=(6,0)
plt.figure()
CS=plt.contour(Q,P,Z,[0.005*(i+1) for i in range(8)],colors = ['k']*8)
#CS=plt.contour(Q,P,Z,20,colors=['w']*0+['k']*20)
#zc=CS.collections[11]

cset3 = plt.contour(Q, P, Z, (0,),
                colors = 'k',
                linewidths = 4,
                hold='on')
plt.clabel(cset3,inline=1,fontsize=10,)

#plt.setp(zc,linewidth=4)
plt.clabel(CS,inline=1,fontsize=10,)
plt.ylim([0,0.12])
plt.xlim([0,0.007])

plt.xlabel("q")
plt.ylabel("p")
plt.savefig("planar_%dcode_%s_lies.png"%(size,error_type))
plt.show()


## LIES, XYZ ERRORS, 5-code

## LIES Y ERRORS, 5-code 


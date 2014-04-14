
import statePlanar as sp
import itertools
import sys
import numpy as np
import copy

## LATTICE SIZE
d=5

L=sp.PlanarLattice(d)
stabilisers= L.x_stabiliser_indices + L.z_stabiliser_indices
qubits = L.qubit_indices




f = open("5syndromes_y_lies.txt","r");
results = [[[int(z) for z in y.split()]
        for y in x.split(",") if y!="\n"]
       for x in f.readlines()[:10]];
results = np.array(results);
f.close();

f= open("size5_y_lies.txt","r")
slies_list = [[int(y) for y in x.split()] for x in f.readlines()[:10]]
f.close();

print results

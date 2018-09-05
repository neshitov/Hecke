from sage.all import *
import itertools
import json
import pickle
from progressbar import progressbar
import Integral_motivic_decomposition as md

L=RootSystem(['G',2]).root_lattice()
I_P=[2]
Cartan, W, WP, reduced,w0,e,s,sP,max_length,length_list, block_dim, pivots, raising,S,a,Hasse=md.init(L,I_P)
basis={}


ind=md.endomorphism_basis_indices(0)
diags=[md.augmented_diagonal(x[0],x[1],0) for x in ind]
diags=md.remove_duplicates(diags)
print('basis diagonals')
for x in diags:
    print(x)
print('all')

#for x in diags:
#    print(x)

p=diags[-2]
p=md.square(p)
q=md.id_complement(p)
pro1=[md.project(x,0,p) for x in diags]
pro2=[md.project(x,0,q) for x in diags]
g=md.add(pro1[0],pro1[3])

#g=md.basis_idempotents(diags+[p,md.zero_diagonal()],md.identity_diagonal())
#print(g)
#print('orthogonal idempotents',len(g))
#print(g)
m=matrix(ZZ,[[0,-1,0,2,3,2],
[0, 0, -3, -6, -6, -3],
[0, 0, 0, -6, -6, 0],
[0, 0, 0, -2, -4, -2],
[0, 0, 0, 0, 0, -6],
[0, 0, 0, 0, 0, 18],
[0, 0, 0, 0, 0, 32],
[0, 0, 0, 0, 18, 18],
[0, 0, 0, 4, 6, 2],
[0, 0, 0, 6, 0, -6],
[0, 1, 1, 0, -1, -1],
[1, 1, 1, 1, 1, 1]])
print(m.echelon_form())

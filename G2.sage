from sage.all import *
import itertools
L=RootSystem(['G',2]).root_lattice()
#charactersitic of the coefficient field
p=2
# The set of simple roots corresponding to the Parabolic subgroup P
I_P=[2]
Cartan=CartanMatrix(L)
W=WeylGroup(L)
# list of reduced words decompositions for w in W
reduced=dict([(w,(w.reduced_word())) for w in W])
w0=W.long_element()
e=W.unit()
s=W.simple_reflections()
#set of simmple reflections that generate P
sP=[s[i] for i in I_P]
#function that checks if w is a minimal representative for its right coset wW_P
def is_minimal(w):
    for i in I_P:
        if (w*s[i]).length()<w.length():
            ans=False
            break
    else:
        ans=True
    return ans
# The set of minimal right coset representatives wW_P:
WP=[w for w in W if is_minimal(w)]
max_length=max([w.length() for w in WP])
# length_list[i] = list of all elements of WP of length i
length_list=[[w for w in WP if w.length()==i] for i in range(0, max_length+1)]
# block_dimension[i]=dimension of the block corresponding to the length i
block_dimension=[len(length_list[i]) for i in range(0, max_length+1)]
#pivots returns the list of pivot elements in WP
def rules(w):
    return [(s[i]*w).length()<w.length() or (s[i]*w not in WP) for i in I_P]
pivots=[w for w in WP if all(rules(w))]

print('number of pivots')
print(len(pivots))


# for every pivot w lovering[w] gives a list of reflections in I_P such that l(s[i]w)<l(w)
raising=dict([(w,[i for i in I_P if (s[i]*w).length()>w.length()]) for w in pivots])
# length_list[i] = list of all elements of WP of length i
length_list=[[w for w in WP if w.length()==i] for i in range(0, max_length+1)]
a=['a%s'%r for r in range(1, (L.rank())+1)]
S=PolynomialRing(GF(p),a)
print(S)
S.inject_variables()
print(Cartan)
#n-th simple root as the element of the polynomial ring
def alpha(n):
    return S.gens_dict()['a%s'%n]
zero=alpha(1)-alpha(1)
unit=alpha(1)//alpha(1)
#action of the n-th simple reflection on the polynomial ring in simple roots
def reflection(S,f,n):
    return f(tuple([alpha(i)-Cartan[n-1,i-1]*alpha(n) for i in range(1, (L.rank()+1))]))
#action of the differential operator Delta_i
def Delta(S,f,i):
    return (f-reflection(S,f,i))//alpha(i)
#extends value from pivots to all entries
def compute(f):
    ans=f
    def descend(w):
        for i in I_P:
            if (s[i]*w).length()<w.length():
                ans[s[i]*w]=-Delta(S,ans[w],i)
                descend(s[i]*w)
    for w in pivots:
        descend(w)
    return ans


#returns list of monomials of degree deg
def monomial(deg):
    part=list(WeightedIntegerVectors(deg, [1 for i in range(0,L.rank())]))
    return [S({tuple(a):1}) for a in part]

#computes the matrix of the linear operator given by n-th simple reflection on the space
#of homogeneous polynomials of degree deg matrix is understood in terms of right multiplication
def reflection_matrix(deg,n):
    mons=monomial(deg)
    size=len(mons)
    m=matrix(GF(p),size,size)
    for i in range(0,size):
        for j in range(0,size):
            m[i,j]=reflection(S,mons[i],n).monomial_coefficient(mons[j])
    return m

monomials=dict([(w,monomial(w.length())) for w in WP])
#returns the list of vectors that represent basis of invariant polynomials in standard monomial basis
def invariant_basis_vector(w):
    size=len(monomials[w])
    if raising[w]==[]:
        big=matrix(size, lambda i,j: 0)
    else:
        big=block_matrix(1,len(raising[w]),[reflection_matrix(w.length(),i)-matrix.identity(size) for i in raising[w]])
    return kernel(big).basis()

#returns the vector space
def invariant_vector_space(w):
    mons=monomial(w.length())
    size=len(mons)
    if raising[w]==[]:
        big=matrix(size, lambda i,j: 0)
    else:
        big=block_matrix(1,len(raising[w]),[reflection_matrix(w.length(),i)-matrix.identity(size) for i in raising[w]])
    return kernel(big)

inv_basis_vector=dict([(w,0) for w in pivots])
for w in pivots:
    inv_basis_vector=dict([(w,invariant_basis_vector(w)) for w in pivots])
    print(w.length(),len(inv_basis_vector))

inv_basis_space=dict([(w,invariant_vector_space(w)) for w in pivots])


#from the tuple of coefficients returns the linear combination of monomials
def vector_to_polynomial(w,x):
    return sum(monomials[w][i]*x[i] for i in range(0,len(monomials[w])))

#gives expansion of the invariant vector in the basis
def polynomial_to_vector(w,f):
    x=tuple(f.monomial_coefficient(monomials[w][i]) for i in range (0,len(monomials[w])))
    return (inv_basis_space[w]).coordinate_vector(x)

#gives the ith coefficient in the basis expansion
def polynomial_coefficient(w,f,i):
    x=tuple(f.monomial_coefficient(monomials[w][i]) for i in range (0,len(monomials[w])))
    return (inv_basis_space[w]).coordinate_vector(x)[i]

def endomorphism_matrix(first):
    M=dict([((v,w),alpha(1)-alpha(1)) for v in WP for w in WP])
    def iterativema(v,w):
        if w==e:
             return first[v]
        else:
            j=reduced[w][0]
            if (s[j]*v).length()==v.length()-1:
                return reflection(S,iterativema(s[j]*v,s[j]*w),j)+Delta(S,iterativema(v,s[j]*w),j)
            elif (s[j]*v).length()==v.length()+1:
                return Delta(S,iterativema(v,s[j]*w),j)
    for w in WP:
        for v in WP:
            M[(v,w)]=iterativema(v,w)
    return M

#return the basis element in the endomorphism space
def basis(w,i):
    x=dict([(v,zero) for v in WP])
    x[w]=vector_to_polynomial(w,inv_basis_vector[w][i])
    return endomorphism_matrix(compute(x))

#for each element of the basis returns the list of diagonal blocks of the basis endomorphism
def basis_diag_blocks(w,i):
    return [dict([((x,y),basis(w,i)[(x,y)]) for x in length_list[j] for y in length_list[j]]) for j in range (0, max_length+1)]

#def block_print(w,i):
#    M=basis_diag_blocks(w,i)
#    for j in range (0, max_length+1):
#        print('length',j)
#        for x in length_list[j]:
#            for y in length_list[j]:
#                print(M[j][(x,y)]),
#            print('\n')


def block_print(w,i):
    M=basis_diag_blocks(w,i)
    for j in range (0, max_length+1):
        print(j,M[j])

# returns the composition of two endomorphisms in terms of matrices
def mult(f,g):
    M=dict([((v,w),zero) for v in WP for w in WP])
    for v in WP:
        for w in WP:
            M[(v,w)]=sum(f[(v,r)]*g[(r,w)] for r in WP)
    return M

inv_dim=dict([(w,inv_basis_space[w].dimension()) for w in pivots])
indices=[(w,i) for w in pivots for i in range(0,inv_dim[w])]

#dimension of the space of endomorphisms

#returns the 3d array c[i,j,k] in F_p such that f_if_j= sum c[i,j,k]f_k where
#f_k are endomorphisms from the basis

def array((x1,x2),(y1,y2),(w,i)):
    return polynomial_coefficient(w,(mult(basis(x1,x2),basis(y1,y2))[(w,e)]),i)
def structure_matrix():
    ans=dict([((a,b,c),0) for a in indices for b in indices for c in indices ])
    for a in indices:
        for b in indices:
            for c in indices:
                ans[(a,b,c)]=array(a,b,c)
                print(ans[(a,b,c)])

#print(t)
for i in range (0,max_length+1):
    print(i,[reduced[w] for w in length_list[i]])
print(block_dimension)

print('number of basis endomorphisms:',len(indices))

#out=open('Dropbox/Hecke project/blocksF4(123).txt','w')

for w in pivots:
    for i in range(0,inv_dim[w]):
        #block_print(w,i)
        print [[[basis(w,i)[(x,y)] for x in length_list[j]] for y in length_list[j]] for j in range (0, max_length+1)]
#g=structure_matrix()

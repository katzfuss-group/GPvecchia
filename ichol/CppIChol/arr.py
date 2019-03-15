import pdb
import numpy as np
from scipy.sparse import csr_matrix
import sys
sys.path.append("/home/marcin/pyMRA/pyMRA")
from MRATools import dispMat


def dot_prod(u1, iter1, u2, iter2, inds, vals):

    result = 0.0

    while(u1<=iter1 and u2<=iter2):
        if(inds[u1]==inds[u2]):
            result += vals[u1]*vals[u2]
            u1 += 1
            u2 += 1
        elif(inds[u1]<inds[u2]):
            u1 += 1
        else:
            u2 += 1
    return result




def ichol(ptrs, inds, vals):

    N = len(ptrs)-1
    for i in range(N):
        for j in range(ptrs[i], ptrs[i+1]):

            u1 = ptrs[i]
            u2 = ptrs[inds[j]]

            if(j==2):
                print(u2)

            dp = dot_prod( u1, ptrs[i+1]-2,
                           u2, ptrs[inds[j] + 1] - 2,
                           inds,
                           vals)

            if inds[j] < i:
                
                vals[j] = (vals[j] - dp)/ vals[ ptrs[inds[j] + 1] -1 ]
                
            elif inds[j]==i:                
                vals[j] = np.sqrt(vals[j] - dp)
            else:
                print("ERROR")                
    
    return vals






#np.set_printoptions(precision=3, suppress=True)

np.random.seed(10)
N = 3
Nval = 10

if N==7:
    S = np.matrix(np.zeros((N,N)))
    S = S + np.matrix(np.eye(N))
    S[0,:] = S[:,0] = 1
    S[1,(3, 4)] = S[(4,3),1] = 1
    S[2,(5, 6)] = S[(5,6),2] = 1
    Srev = S[::-1,::-1]
else:
    S = np.ones((N,N))

vals = np.random.normal(size=N**2)
M = np.matrix(vals.reshape((N,N)))

M = np.linalg.inv(np.multiply(S,M*M.T)) + 10*np.eye(N)
cM = np.linalg.cholesky(M)


inds = np.triu_indices(N,1)
M[inds] = 0
Mcsr = csr_matrix(M)

ptrs = Mcsr.indptr
inds = Mcsr.indices
vals = Mcsr.data

vals = ichol(ptrs, inds, vals)

pdb.set_trace()
print(vals)

cMcsr = csr_matrix((vals, inds, ptrs), shape=(N,N))
diff = np.max(np.abs(cMcsr.todense() - cM))
#print(diff)

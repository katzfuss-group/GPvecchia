import logging
import pdb
from scipy.spatial.distance import cdist
import numpy as np
import sys
import numpy as np
sys.path.append("/home/marcin/pyMRA/pyMRA/")
from MRATools import dispMat, filterNNZ, Matern

fill=1
N = 4*fill + 3
lambd = 0.1

np.set_printoptions(precision=3)


def ichol_ip(A, S):

    n = A.shape[0]
    A = np.array(A)
    for i in range(n):
        for j in range(i):
            if S[i,j]!=0:
                A[i,j] = (A[i,j] - sum(A[i,:j]*A[j,:j]))/A[j,j]
            else:
                logging.debug("Skipping entry (%d,%d)" % (i,j))
        A[i,i] = np.sqrt(A[i,i] - sum(A[i,:i]**2))

    for i in range(n):
        for j in range(i+1,n):
            A[i,j] = 0
    return np.matrix(A)



def ichol(L, S):

    n = L.shape[0]
    for i in range(n):
        for j in range(i):
            if S[i,j]!=0:
                L[i,j] = (L[i,j] - sum(L[i,:j]*L[j,:j]))/L[j,j]
            else:
                logging.debug("Skipping entry (%d,%d)" % (i,j))
        L[i,i] = np.sqrt(L[i,i] - sum(L[i,:i]**2))
    return np.matrix(L)





def icholCompare(Sigma, S):

    n = Sigma.shape[0]
    SigmaTil = np.multiply(Sigma,S)
    Li = np.zeros((n,n))
    L = np.zeros((n,n))
    for i in range(n):
        for j in range(i):

            L[i,j] = (SigmaTil[i,j] - sum(L[i,:j]*L[j,:j]))/L[j,j]
            
            if S[i,j]!=0:
                Li[i,j] = (Sigma[i,j] - sum(Li[i,:j]*Li[j,:j]))/Li[j,j]

            if Li[i,j] != L[i,j]:
                pdb.set_trace()
                
        L[i,i] = np.sqrt(SigmaTil[i,i] - sum(L[i,:i]**2))
        Li[i,i] = np.sqrt(Sigma[i,i] - sum(Li[i,:i]**2))

    
    return np.matrix(Li), np.matrix(L)



		    



if __name__=='__main__':


    locs = np.linspace(0, 1, N).reshape((N,1))
    mra_ord = np.concatenate([np.arange(3,3+fill), np.array([1]), np.arange(3+fill, 3+2*fill), np.array([0]), np.arange(3+2*fill, 3+3*fill), np.array([2]), np.arange(3+3*fill,3+4*fill)])
    ord_vec = np.argsort(mra_ord)
    locsord = locs[ord_vec,:]#[::-1]


    locsord = np.array([[0.3676648], [0.5539822], [0.1852040], [0.4492080]])
    D = np.matrix(cdist(locsord, locsord))
    Theta_mra = np.matrix(np.exp(-D/lambd))
    Theta = Matern(locsord, nu=0.5, l=0.1)


    Ac = ichol_ip(Theta, np.ones((N,N)))

    pdb.set_trace() 

    
    Theta_mrized = Theta[:,0] * Theta[0,:]
    m1 = np.zeros((N,2))


    
    S = np.matrix(np.eye(N))
    S[abs(A_mra)>1e-6] = 1
    S[0,:] = 1; S[:,0] = 1
    #S[N-1,:] = 1; S[:,N-1] = 1;
    locsord = locsord.ravel()
    S[locsord<locsord[0],1] = 1; S[locsord>locsord[0],2] = 1;
    S[1,locsord<locsord[0]] = 1; S[2,locsord>locsord[0]] = 1; 
    #S[locsord<locsord[N-1],N-2] = 1; S[locsord>locsord[N-1],N-3] = 1;
    #S[N-2,locsord<locsord[N-1]] = 1; S[N-3,locsord>locsord[N-1]] = 1; 
    
    Annz = set(map(tuple, np.vstack(np.where(np.abs(A_mra)>1e-6)).T))
    Snnz = set(map(tuple, np.vstack(np.where(S)).T))
   
    if len(Annz - Snnz):
        print("error")

    #create \hat{Sigma}
    Theta_mrized = np.linalg.inv(np.multiply(A,S))
    dispMat(np.linalg.cholesky(Theta_mrized))
    





            
        
    P = np.matrix(np.zeros((N,N)))
    #for i in range(N):
    #    P[i,N-1-i] = 1


    U = np.linalg.cholesky(A)
    L = P*(np.linalg.inv(U).T)*P
    
    pdb.set_trace()



        

    Li = ichol(Theta, S)      
    L = np.linalg.cholesky(Theta)
    
    #dispMat(filterNNZ(Li, tol=1e-10))
    #dispMat(filterNNZ(L, tol=1e-10))
    print(np.max(np.abs(Li - L)))


    
    #LinvT = np.linalg.inv(L).T
    #PLitP = P * LinvT * P
    #print(np.max(np.abs(P*A*P - PLitP * PLitP.T)))

import numpy as np 
import calfem.core as cfc

def uppgift1():
    NELEM, NDOF = 3, 4

    ex = np.array([2,
                   4,
                   6
                ])
    ey = np.array([4,
                   6,
                   8
                ])
    
    edof = np.array([
        [1, 2],
        [2, 3],
        [3, 4]
    ])


    A = 10
    k = 5
    Q = 100
    L = 2

    K = np.zeros((NDOF,NDOF))

    kei = np.zeros((2,2))
    for i in np.arange(0,NELEM):
        kei = cfc.spring1e(A*k/L)
        cfc.assem(edof[i,:], K, kei)
    
    print(K)

    F = np.zeros((NDOF,1))
    F[3] = -15*A
    bc_dof = np.array([1])
    bc_val = np.array([0])
    F_load = [[100], [200], [200], [100]]

    F=F+F_load

    a,r = cfc.solveq(K, F, bc_dof, bc_val)
    print("temperature")
    print(a)
    print("heat flow")
    print(r)
uppgift1()
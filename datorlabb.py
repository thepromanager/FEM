import numpy as np 
import calfem.core as cfc
import matplotlib.pyplot as plt
import geom2

def uppgift1():
    NELEM, NDOF = 20, 21
    L = 6/NELEM

    ex = np.array([2+ x*L for x in range(0,NELEM)
                ])
    ex = np.array([2+ x*L for x in range(1,NELEM+1)
                ])
    
    edof = np.array([[i,i+1] for i in range(1,NELEM+1)
    ])


    A = 10
    k = 5
    Q = 100

    K = np.zeros((NDOF,NDOF))

    kei = np.zeros((2,2))
    for i in np.arange(0,NELEM):
        kei = cfc.spring1e(A*k/L)
        cfc.assem(edof[i,:], K, kei)
    
    print(K)

    F = np.zeros((NDOF,1))
    F[NELEM] = -15*A
    bc_dof = np.array([1])
    bc_val = np.array([0])
    F_load = Q*L*np.ones((NDOF,1))
    F_load[0], F_load[-1] = Q*L*0.5, Q*L*0.5

    F=F+F_load

    a,r = cfc.solveq(K, F, bc_dof, bc_val)
    print("temperature")
    print(a)
    print("heat flow")
    print(r)

    x = np.arange(2,8.1,L)
    plt.plot(x,a)
    plt.show()

def uppgift2():
    coord, edof, dofs, bdofs, bc, bc_value, element_markers = geom2.generate_mesh(False)
    k = 1
    
    nDofs = np.size(dofs)
    ex, ey = cfc.coordxtr(edof, coord, dofs)

    K = np.zeros([nDofs, nDofs])

    for eltopo, elx, ely in zip(edof, ex, ey):
        Ke = cfc.flw2te(elx,ely,[1],k*np.eye(2))
        K = cfc.assem(eltopo, K, Ke)

    a,r = cfc.solveq(K,np.zeros((nDofs,1)),bc, bc_value)

    print(a)
    

uppgift2()
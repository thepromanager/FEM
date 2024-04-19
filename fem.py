import numpy as np 
import calfem.core as cfc

#define truss structure:

NELEM, NDOF = 3,8
L = 1

# en rad per element, beskriver nodernas läge?
ex = np.array([
            [-L, 0],
            [0, 0],
            [L, 0]
            ])

ey = np.array ([
            [0, -L],
            [0, -L],
            [0, -L]
 ])

#vilka degrees of freedom varje element är connectat med
edof = np.array([
    [1, 2, 7, 8],
    [3, 4, 7, 8], 
    [5, 6, 7, 8]
])

# Bar parameters
E = 210e9 #elasticitetsmodul
A = np.array ([1e-4, 1e-4, 1e-4]) #tvärsnittsareor

#define stiffness matrix:
K = np.zeros((NDOF,NDOF))

#compute stiffness matrix. räkna ut element stiffness matrix för varje element (kei) med cfc.bar2e. 
#assembla med cfc.assem. 
kei = np.zeros((4,4))
for i in np.arange(0,NELEM):
    kei = cfc.bar2e(ex[i, :], ey[i, :], [E, A[i]])
    cfc.assem(edof[i, :], K, kei)

print(K)
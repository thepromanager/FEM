import numpy as np 
import calfem.core as cfc
import matplotlib.pyplot as plt
import geom2
import matplotlib as mpl

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu

import plantml
g = cfg.geometry()

# define parameters
scale=1
Height=200*scale
bigWidth=400*scale
smallWidth=300*scale
circleRadius=25*scale
circleSpacing = 87.5*scale
DEPTH = 1600*scale


# Mesh data
el_sizef, el_type, dofs_pn = 10, 2, 1
mesh_dir = "./"

# boundary markers
MARKER_T_293=0 #293 K
MARKER_T_277=1
MARKER_T_285=2
MARKER_QN_0=3


# add points
g.point([0, 0], 0)
g.point([smallWidth, 0], 1)
g.point([bigWidth, Height], 2)
g.point([0, Height], 3)
g.point([0, Height/2+circleRadius], 4)
g.point([0, Height/2], 5)
g.point([circleRadius, Height/2], 6)
g.point([0, Height/2-circleRadius], 7)


# define lines / circle segments

g.spline([0, 1], 0, marker=MARKER_QN_0)
g.spline([1, 2], 1, marker=MARKER_QN_0)
g.spline([2, 3], 2, marker=MARKER_T_293)
g.spline([3, 4], 3, marker=MARKER_QN_0)
g.circle([4, 5, 6], 4, marker=MARKER_T_285)
g.circle([6, 5, 7], 5, marker=MARKER_T_285)
g.spline([7, 0], 6, marker=MARKER_QN_0)

pointIndex=8
splineIndex=7

def circle(x,y,Marker,pointIndex,splineIndex):
	g.point([x, y], pointIndex)
	g.point([x-circleRadius, y], pointIndex+1)
	g.point([x, y+circleRadius], pointIndex+2)
	g.point([x+circleRadius, y], pointIndex+3)
	g.point([x, y-circleRadius], pointIndex+4)

	g.circle([pointIndex+1, pointIndex, pointIndex+2], splineIndex, marker=Marker)
	g.circle([pointIndex+2, pointIndex, pointIndex+3], splineIndex+1, marker=Marker)
	g.circle([pointIndex+3, pointIndex, pointIndex+4], splineIndex+2, marker=Marker)
	g.circle([pointIndex+4, pointIndex, pointIndex+1], splineIndex+3, marker=Marker)
	pointIndex=pointIndex+5
	splineIndex=splineIndex+4
	return (pointIndex,splineIndex)
(pointIndex,splineIndex)= circle(circleSpacing,Height/2,MARKER_T_277,pointIndex,splineIndex)
(pointIndex,splineIndex)= circle(circleSpacing*2,Height/2,MARKER_T_285,pointIndex,splineIndex)
(pointIndex,splineIndex)= circle(circleSpacing*3,Height/2,MARKER_T_277,pointIndex,splineIndex)

# define surface
g.surface([0, 1, 2, 3,4,5,6],holes=[[i,i+1,i+2,i+3] for i in (7,11,15)])

# Create Mesh
mesh = cfm.GmshMeshGenerator(g, mesh_dir=mesh_dir)
mesh.el_size_factor = el_sizef
mesh.el_type = el_type
mesh.dofs_per_node = dofs_pn
coords, edof, dofs, bdofs, element_markers = mesh.create()

#Material Parameters
k = 80/1000
alpha_c = 120/1000_000
alpha_n = 40/1000_000

#Apply boundry condition
bc, bc_value = np.array([], 'i'), np.array([], 'f')
#bc, bc_value = cfu.applybc(bdofs, bc, bc_value, MARKER_T_285, 285, 1)
#bc, bc_value = cfu.applybc(bdofs, bc, bc_value, MARKER_T_277, 277, 1)
#bc, bc_value = cfu.applybc(bdofs, bc, bc_value, MARKER_T_293, 293, 1)

#Making K Matrix
nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

K = np.zeros([nDofs, nDofs])


def NtN(ex: np.array, ey: np.array, s: float): # First node is not on boundry
    if not ex.shape == (3,) or not ey.shape == (3,):
        raise Exception("Incorrect shape of ex or ey: {0}, {1} but should be (3,)".format(ex.shape, ey.shape))
    Cmat = np.vstack((np.ones((3, )), ex, ey))

F=np.zeros((nDofs,1))
for  i in range(0,np.shape(ex)[0]):
    Ke = cfc.flw2te(ex[i,:],ey[i,:],[DEPTH], k*np.eye(2))
    #print(cfc.flw2te(np.array([850,1000,1000]),np.array([750,700,800]),[1],np.eye(2)))
    boundryTemp=0
    nodes=[]
    for marker in [0,1,2]:
        temp=0
        nodesOnBoundry=[]
        for j in range(3):
            node = edof[i,:][j]
            if(node in bdofs[marker]):
                temp = [293,277,285][marker]
                nodesOnBoundry.append(j)
        if(len(nodesOnBoundry)>1):
            boundryTemp=temp
            nodes=nodesOnBoundry
    if(boundryTemp):

        alpha = (alpha_n)*(boundryTemp==293)+alpha_c*(boundryTemp==277 or boundryTemp==285)

        #temp=plantml.plantml(ex[i,:],ey[i,:],alpha*DEPTH)
        #temp=plantml.plantml(np.array([850,1000,1000]),np.array([750,700,800]),1)
        firstPoint=np.array([ex[i,:][nodes[0]],ey[i,:][nodes[0]]])

        secondPoint=np.array([ex[i,:][nodes[1]],ey[i,:][nodes[1]]])
        L=np.linalg.norm(firstPoint-secondPoint)
        
        temp = np.array([[0, 0, 0],[0,0,0],[0,0,0]])
        temp[nodes[0]][nodes[0]]=L/3
        temp[nodes[1]][nodes[1]]=L/3
        temp[nodes[0]][nodes[1]]=L/6
        temp[nodes[1]][nodes[0]]=L/6
        temp=alpha*DEPTH*temp
        Ke += temp

        node1 = edof[i,:][nodes[0]]
        node2 = edof[i,:][nodes[1]]
        F[node1]=F[node1]+L/2*alpha*DEPTH*boundryTemp
        F[node2]=F[node2]+L/2*alpha*DEPTH*boundryTemp

    K = cfc.assem(edof[i,:], K, Ke)
#making F

new_arr=F[np.where(F!=0)]
print(np.size(new_arr))
print(nDofs)


#solve
#a,r = cfc.solveq(K,F,bc, bc_value)
a=np.linalg.solve(K,F)
a=a.T.tolist()[0]
mi, ma = 10000000,-100000
for t in a:
    if t>ma:
        ma=t
    if t<mi:
        mi=t
print(mi,ma)

fig, ax = plt.subplots()

cfv.draw_geometry(
    g,
    label_curves=True,
    title="Geometry: Computer Lab Exercise 2"
)
cfv.drawMesh(
            coords=coords,
            edof=edof,
            dofs_per_node=mesh.dofsPerNode,
            el_type=mesh.elType,
            filled=False,
            title="Example 01",
            
        )

cfv.draw_nodal_values_shaded(a,coords,edof,title=None,dofs_per_node=mesh.dofs_per_node,el_type=mesh.el_type, draw_elements=False)
cfv.colorbar()
plt.show()
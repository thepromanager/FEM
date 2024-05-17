import numpy as np 
import calfem.core as cfc
import matplotlib.pyplot as plt
import matplotlib as mpl

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv
import calfem.utils as cfu
from plantml import plantml

g = cfg.geometry()

# Geometry Parameters
scale=1/1000
Height=200*scale
bigWidth=400*scale
smallWidth=300*scale
circleRadius=25*scale
circleSpacing = 87.5*scale
DEPTH = 1600*scale

#Simulation Parameters
type_of_charging = 1 #1 eller 2, motsvarar Q = f1(t) resp Q = f2(t)
timestep = 10 #definiera i sekunder
totTime = 920*(type_of_charging==1)+600*(type_of_charging==2) # 920 (charge 1) 600 ()


#Thermic Parameters
k = 80 #/1000
alpha_c = 120 #/1000_000
alpha_n = 40 #/1000_000
rho = 540 #/ 1000_000_000 #kg/mm^(3)
c_p = 3600 #J/(kgK)

#Mechanical Parameters
young_E = 5 * 10**9
poission_v = 0.36
expansion_a = 60 * 10**(-6)


# Mesh data
el_sizef, el_type, dofs_pn = 10*scale, 2, 1
mesh_dir = "./"
dofs_pn2=2


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

g.spline([0, 1], 0, marker=MARKER_QN_0) # also u = 0
g.spline([1, 2], 1, marker=MARKER_QN_0) # also u = 0
g.spline([2, 3], 2, marker=MARKER_T_293) # also t = 0
g.spline([3, 4], 3, marker=MARKER_QN_0) # also u = 0
g.circle([4, 5, 6], 4, marker=MARKER_T_285)
g.circle([6, 5, 7], 5, marker=MARKER_T_285)
g.spline([7, 0], 6, marker=MARKER_QN_0) # also u = 0

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

# Create Mesh2
mesh2 = cfm.GmshMeshGenerator(g, mesh_dir=mesh_dir)
mesh2.el_size_factor = el_sizef
mesh2.el_type = el_type
mesh2.dofs_per_node = dofs_pn2
coords2, edof2, dofs2, bdofs2, element_markers2 = mesh2.create()

#Making K Matrix
nDofs = np.size(dofs)
ex, ey = cfc.coordxtr(edof, coords, dofs)

K = np.zeros([nDofs, nDofs])
KC = np.zeros([nDofs, nDofs])
C = np.zeros([nDofs, nDofs]) #C-matrisen i Cda/dt + K*a = F
F=np.zeros((nDofs,1))
first_F=np.zeros((nDofs,1))


for  i in range(0,np.shape(ex)[0]):
    Ke = cfc.flw2te(ex[i,:],ey[i,:],[DEPTH], k*np.eye(2)) #gör normala Ke
    boundryTemp=0
    temperature=0 #initialisera temperaturer
    nodes=[]
    for marker in [0,1,2]:
        tempList = np.array([293,277,285]) #lista över temperaturer
        temperature = tempList[marker] #sätt temperaturen vid detta ställe
        nodesOnBoundry=[] 
        
        for j in range(3):
            #för varje frihetsgrad för detta element
            node = edof[i,:][j] #frihetsgraden
            
            if(node in bdofs[marker]):
                nodesOnBoundry.append(j) #ifall den frihetsgraden är vid en kant
                
        if(len(nodesOnBoundry)>1):
            boundryTemp=temperature
            nodes=nodesOnBoundry
            #sparar temperaturen och index för frihetsgraderna vid randen i edof[i,:]
            
    if(boundryTemp != 0):
        
        alpha = (alpha_n)*(boundryTemp==293)+alpha_c*(boundryTemp==277 or boundryTemp==285)

        firstPoint=np.array([ex[i,:][nodes[0]],ey[i,:][nodes[0]]])

        secondPoint=np.array([ex[i,:][nodes[1]],ey[i,:][nodes[1]]])
        L=np.linalg.norm(firstPoint-secondPoint)
        
        node1 = edof[i,:][nodes[0]]
        node2 = edof[i,:][nodes[1]]
        
        #gör integralen över randen för K_c
        KC[node1-1][node1-1] += L/3*alpha*DEPTH
        KC[node2-1][node2-1] += L/3*alpha*DEPTH
        KC[node1-1][node2-1] += L/6*alpha*DEPTH
        KC[node2-1][node1-1] += L/6*alpha*DEPTH

        #gör integralen över randen för F
        F[node1-1]=F[node1-1]+L/2*alpha*DEPTH*boundryTemp
        F[node2-1]=F[node2-1]+L/2*alpha*DEPTH*boundryTemp
        
        
        
    
    
    #assemblar styvhetsmatrisen
    K = cfc.assem(edof[i,:], K, Ke)
    
    #assemblar C
    C = cfc.assem(edof[i,:], C, plantml(ex[i,:],ey[i,:], rho*DEPTH*c_p))
    
print("klar med konstanta grejer")
K = K + KC

np.copyto(first_F,F) #spara den konstanta delen av F


#själva timesteppingen:
maximum_a_values = np.array([])
minimum_a_values = np.array([])
max_T_T0 = np.array([])
inverse_matrix = np.linalg.inv(C+timestep*K)
time_list = np.array([x*timestep for x in range(0,(totTime+1)//timestep)])
a = 293*np.ones((nDofs,1)) #initialisera med konstant temperatur i batteriet, T=293K
m = 0 #printar tiden lite då och då
Cs=[]
for i in range(0,np.shape(ex)[0]):     
    el_node1,el_node2,el_node3 = edof[i,:][0],edof[i,:][1],edof[i,:][2]
    node1_coord = np.array([ex[i,:][0],ey[i,:][0]])
    node2_coord = np.array([ex[i,:][1],ey[i,:][1]])
    node3_coord = np.array([ex[i,:][2],ey[i,:][2]])
    L1,L2,L3 = node2_coord-node1_coord,node1_coord-node3_coord, node3_coord-node2_coord

    cross1,cross2,cross3 = np.cross(L1,L2),np.cross(L1,L3),np.cross(L2,L3)
    C1,C2,C3 = np.linalg.norm(cross1),np.linalg.norm(cross2),np.linalg.norm(cross3)
    
    Cs.append((C1,C2,C3))
for current_time in time_list:
    
    
    np.copyto(F,first_F) 
    
    if type_of_charging == 1:
        
        future_Q = 100*np.exp(-144*((600-(current_time+timestep))/3600)**2)
             
    elif type_of_charging == 2:
   
        future_Q = 88.42 if (600-(current_time+timestep)>0) else 0

        
    else:
        
        print("fel typ av uppladdning")
    future_Q = future_Q * 1000 #/ 1_000_000_000#fixar enheter till mm^(-3)
    

    for i in range(0,np.shape(ex)[0]):
        
        el_node1,el_node2,el_node3 = edof[i,:][0],edof[i,:][1],edof[i,:][2]
        #node1_coord = np.array([ex[i,:][0],ey[i,:][0]])
        #node2_coord = np.array([ex[i,:][1],ey[i,:][1]])
        #node3_coord = np.array([ex[i,:][2],ey[i,:][2]])
        #L1,L2,L3 = node2_coord-node1_coord,node1_coord-node3_coord, node3_coord-node2_coord

        #cross1,cross2,cross3 = np.cross(L1,L2),np.cross(L1,L3),np.cross(L2,L3)
        #D1,D2,D3 = np.linalg.norm(L1),np.linalg.norm(L2),np.linalg.norm(L3)
        #C1,C2,C3 = np.linalg.norm(cross1),np.linalg.norm(cross2),np.linalg.norm(cross3)

        #utgår här ifrån konstant fördelning av Q, hoppas d e okej lol
        #F[el_node1-1] = F[el_node1-1]+(D1*D2/2)*DEPTH*future_Q
        #F[el_node2-1] = F[el_node2-1]+(D1*D3/2)*DEPTH*future_Q
        #F[el_node3-1] = F[el_node3-1]+(D3*D2/2)*DEPTH*future_Q #OBS osäker om d bör vara + eller -
        (C1,C2,C3)=Cs[i]

        F[el_node1-1] = F[el_node1-1]+(C1/6)*DEPTH*future_Q
        F[el_node2-1] = F[el_node2-1]+(C2/6)*DEPTH*future_Q
        F[el_node3-1] = F[el_node3-1]+(C3/6)*DEPTH*future_Q
        
        
    
    
    
    maximum_a_values = np.append(maximum_a_values, np.max(a))
    minimum_a_values = np.append(minimum_a_values, np.min(a))
    max_T_T0 = np.append(max_T_T0,np.abs(max(a-293*np.ones((nDofs,1)), key=abs)))
    a = np.matmul(inverse_matrix, np.matmul(C,a)+timestep*F) #ett steg framåt
    
    if m%50 == 0:
        
        print(f"tid: {current_time}s")
    m+=1
    

print(np.max(a),np.min(a))


# REMOVE LATER
nDofs2 = np.size(dofs2)
ex2, ey2 = cfc.coordxtr(edof2, coords2, dofs2)
Kmech = np.zeros([nDofs2, nDofs2])
ptype=2
D = np.eye(3)*young_E
badD = cfc.hooke(ptype,young_E,poission_v)

F=np.zeros((nDofs2,1))

D[0,0]=badD[0,0]
D[1,0]=badD[1,0]
D[0,1]=badD[0,1]
D[1,1]=badD[1,1]
D[2,2]=badD[3,3]
print(D)
for  i in range(0,np.shape(ex2)[0]):
    
    Ke = cfc.plante(ex2[i,:],ey2[i,:],[ptype,DEPTH], D) #
    #Kmech = cfc.assem(edof2[i,:], Kmech, Ke)

    el_node1,el_node2,el_node3 = edof[i,:][0],edof[i,:][1],edof[i,:][2]
    temp1,temp2,temp3 = a[el_node1-1],a[el_node2-1],a[el_node3-1]

    #sigx = young_E/(2*poission_v-1)*expansion_a*T
    #sigy = young_E/(2*poission_v-1)*expansion_a*T
    #sigz = young_E/(2*poission_v-1)*expansion_a*T
    #tauxy = 0
    deltaT=np.abs(293-(temp1+temp2+temp3)/3)
    sig = np.matmul(D,expansion_a*deltaT*np.array([1,1,0]))
    fe=cfc.plantf(ex2[i,:],ey2[i,:],[ptype,DEPTH],np.array([sig]))
    (Kmech,F)=cfc.assem(edof2[i,:],Kmech,Ke,F,fe)


bc, bc_value = np.array([], 'i'), np.array([], 'f')
bc, bc_value = cfu.applybc(bdofs2, bc, bc_value, MARKER_QN_0, 0, 0)

#f=np.zeros((nDofs2,1))
u,rowan = cfc.solveq(Kmech,F,bc, bc_value)
ed=cfc.extract_eldisp(edof2,u)
vonMises=[]
maxVonMises=0
for  i in range(0,np.shape(ex2)[0]):
    [es,et]= cfc.plants(ex2[i,:],ey2[i,:],[ptype,DEPTH], D, ed[i,:])

    el_node1,el_node2,el_node3 = edof[i,:][0],edof[i,:][1],edof[i,:][2]
    temp1,temp2,temp3 = a[el_node1-1],a[el_node2-1],a[el_node3-1]
    deltaT=(temp1+temp2+temp3)/3-293
    Depsilon = (expansion_a*young_E*deltaT/(1-2*poission_v))[0]

    sigzz = ((young_E/(1-2*poission_v))*((poission_v/(1+poission_v))*(et[0][0]+et[0][1])-expansion_a*deltaT))[0]
    sigzz-=-Depsilon
    sigxx = es[0][0] - Depsilon
    sigyy = es[0][1] - Depsilon
    vM = np.sqrt(sigxx*sigxx+sigyy*sigyy+sigzz*sigzz-sigxx*sigzz-sigxx*sigyy-sigyy*sigzz)
    vonMises.append(vM)
    if(vM>maxVonMises):
        maxVonMises=vM

print(maxVonMises)


#u=u.T.tolist()[0]





    
plt.plot(time_list, maximum_a_values,label = "maximum temperatures")
plt.plot(time_list, minimum_a_values, label = "minimum temperatures")

plt.title(f"temperature as function of time using $f_{type_of_charging}(t)$")
plt.legend()
plt.show()

plt.plot(time_list, max_T_T0, label = "$max|T-T_0|$") 
#print("MAXT0T1 index 50:",max_T_T0[50])
#f1: 2.153903685658406, f2: 
plt.show()

fig, ax = plt.subplots()

cfv.drawMesh(
            coords=coords,
            edof=edof,
            dofs_per_node=mesh.dofsPerNode,
            el_type=mesh.elType,
            filled=False,
            title="konvektion + transienter",
            
        )

#cfv.draw_nodal_values_shaded(a,coords,edof,title=None,dofs_per_node=mesh.dofs_per_node,el_type=mesh.el_type, draw_elements=False)
cfv.draw_element_values(vonMises,coords2,edof2,title=None,dofs_per_node=mesh2.dofs_per_node,el_type=mesh2.el_type,draw_elements=False)
cfv.colorbar()
plt.show()
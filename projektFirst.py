import numpy as np 
import calfem.core as cfc
import matplotlib.pyplot as plt
import geom2
import matplotlib as mpl

import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis_mpl as cfv

g = cfg.geometry()

# define parameters
R = 2
L = 10

Height=200
bigWidth=400
smallWidth=300
circleRadius=25
circleSpacing = 87.5


# Mesh data
el_sizef, el_type, dofs_pn = 0.5, 2, 1
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
	print(x,y,pointIndex,splineIndex)
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
print(pointIndex,splineIndex)

# define surface
g.surface([0, 1, 2, 3,4,5,6],holes=[[7,8,9,10]])
#g.surface([7,8,9,10])



fig, ax = plt.subplots()
cfv.draw_geometry(
    g,
    label_curves=True,
    title="Geometry: Computer Lab Exercise 2"
)

plt.show()
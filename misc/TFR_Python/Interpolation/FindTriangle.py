#####################################################################################
#                                                                                   #
#                             Copyright (c) 2014-2016                               #
#				      Manuel Zündel (zuendel@imes.mavt.ethz.ch)                     #
#                       Alexander E. Ehret and Edoardo Mazza                        #
#				       Experimental Continuum Mechanics Group                       #
#				    Institute of Mechanical Systems, ETH Zürich                     #
#				                All rights reserved.                                #
#                                                                                   #
#####################################################################################

#####################################################################################
# For a node with coordinates XY, return the triangle in which the point lies.      #
# Also compute the coodrdinates  xi and nu of the point in the local coordinate		# 
# system of the triangle.													        #
#####################################################################################

import numpy as np

def FindTriangle( triangles,points,XY ):
	
	# Initialize variables and paramenters
	tol=1e-3

	index=-1
	xe_t=False;
	xi=-1
	nu=-1
	
	# Go through all triangles
	for t in range(len(triangles)):
		triangle=triangles[t];
		
		P1=points[triangle.item(0)]
		P2=points[triangle.item(1)]
		P3=points[triangle.item(2)]
		
		
		# is point inside the bounding box of the triangle vertices?
		if XY[0,0]>min(P1[0,0],P2[0,0],P3[0,0])-tol and XY[0,0]<max(P1[0,0],P2[0,0],P3[0,0])+tol and XY[0,1]>min(P1[0,1],P2[0,1],P3[0,1])-tol and XY[0,1]<max(P1[0,1],P2[0,1],P3[0,1])+tol:
			
			# Compute Matrices needed for the computation of local coordinates
			J=np.matrix([[-P1.item(0)+P2.item(0),-P1.item(1)+P2.item(1)],[-P1.item(0)+P3.item(0),-P1.item(1)+P3.item(1)]])
			if np.linalg.det(J)!=0:
				
				M=np.matrix([[P2.item(0)-P1.item(0),P3.item(0)-P1.item(0)],[P2.item(1)-P1.item(1),P3.item(1)-P1.item(1)]])
				b=np.matrix([XY.item(0)-P1.item(0),XY.item(1)-P1.item(1)]).T
				
				# Compute local coordinates
				xe=np.linalg.solve(M,b)
				xi=xe.item(0);
				nu=xe.item(1);
				
				# Is point inside triangle?
				if xi>=0-tol and xi<=1+tol and nu>=0-tol and nu<=1+tol and nu+xi <=1+tol:
					index=t;
					xe_t=xe;

					break
	
	return [index,[xi,nu]]
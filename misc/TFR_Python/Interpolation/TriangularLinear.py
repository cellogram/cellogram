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

# Linear Triangular Interpolation Method 

import numpy as np
import Interpolation

def TriangularLinear(triangles,points,U,XY):
	
	# Get Triangle
	FQ=Interpolation.FindTriangle( triangles,points,XY )
	q=FQ[0]
	
	if q<0:
		return False
	else:
	
		# Get Local Coordinades
		P1=points[triangles[q].item(0)]
		P2=points[triangles[q].item(1)]
		P3=points[triangles[q].item(2)]
		J=np.matrix([[-P1.item(0)+P2.item(0),-P1.item(1)+P2.item(1)],[-P1.item(0)+P3.item(0),-P1.item(1)+P3.item(1)]])
		
		P1def=points[triangles[q].item(0)]+np.matrix([U[triangles[q].item(0)].item(0),U[triangles[q].item(0)].item(1)])
		P2def=points[triangles[q].item(1)]+np.matrix([U[triangles[q].item(1)].item(0),U[triangles[q].item(1)].item(1)])
		P3def=points[triangles[q].item(2)]+np.matrix([U[triangles[q].item(2)].item(0),U[triangles[q].item(2)].item(1)])
		Jdef=np.matrix([[-P1def.item(0)+P2def.item(0),-P1def.item(1)+P2def.item(1)],[-P1def.item(0)+P3def.item(0),-P1def.item(1)+P3def.item(1)]])
		
		# Check that the triangle did not flip (if flipped the linear triangular  interpolation is impossible)
		if np.sign(np.linalg.det(J))*np.sign(np.linalg.det(Jdef))>0:
			xi=FQ[1][0]
			nu=FQ[1][1]
			
			U_int=[]
			# Compute interpolated value for the displacement field
			for i in range( U[0].size):
				U1=U[triangles[q].item(0)].item(i)
				U2=U[triangles[q].item(1)].item(i)
				U3=U[triangles[q].item(2)].item(i)

				u=U1*(1-xi-nu)+U2*xi+U3*nu
				U_int.append(u)

			return U_int
			
		else:
			return False
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

# RBF interpolation methods for 3D displacement fields

import numpy as np
import Interpolation
import sys
import scipy.interpolate

def setRBF(points,U,eps,f='thin-plate'):
	# Compute RBF weigths
	RBFX = scipy.interpolate.Rbf(points[:,0], points[:,1], U[:,0],function=f, epsilon=eps)
	RBFY = scipy.interpolate.Rbf(points[:,0], points[:,1], U[:,1],function=f, epsilon=eps)
	RBFZ = scipy.interpolate.Rbf(points[:,0], points[:,1], U[:,2],function=f, epsilon=eps)
	return[RBFX, RBFY,RBFZ]
	
	
def RBFInterp(XY,RBF):
	# Compute interpolated value
	ux= float(RBF[0](XY[0], XY[1]))
	uy= float(RBF[1](XY[0], XY[1]))
	uz= float(RBF[2](XY[0], XY[1]))
	return [ux,uy,uz]		
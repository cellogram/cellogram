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

# Check if the CAE model is large enough compared to the given set of QD
	
def CheckModelSize(points,CAEBox):
	x_max=0
	y_max=0
	y_min=0
	x_min=0
	
	for xy in points:
		if xy[0,0]>x_max:
			x_max=xy[0,0]				
		if xy[0,0]<x_min:
			x_min=xy[0,0]				
		if xy[0,1]>y_max:
			y_max=xy[0,1]				
		if xy[0,1]<y_min:
			y_min=xy[0,1]
			
		
	if x_min<CAEBox[0]:	
		return [False,[x_min,y_min,x_max,y_max]]	
	if y_min<CAEBox[1]:
		return [False,[x_min,y_min,x_max,y_max]]
	if x_max>CAEBox[2]:
		return [False,[x_min,y_min,x_max,y_max]]
	if y_max>CAEBox[3]:
		return [False,[x_min,y_min,x_max,y_max]]
		
	return [True,[x_min,y_min,x_max,y_max]]
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

# Get all the points within the radius R to point C

import numpy as np
	
def GetPCircle( points,C,R ):
    points_in_C=[]
    indexes=[]
    i=1
    for p in range(len(points)):
        point=points[p]
        if np.linalg.norm(point-C)<R:
            points_in_C.append(point)
            indexes.append(p)
    return indexes

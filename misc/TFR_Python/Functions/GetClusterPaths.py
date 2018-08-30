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
# The recursive function GetClusterPaths iteratively finds all possible paths 		#
# starting from the point with the given index. This function is needed for the 	#
# Adaptive Meshing Scheme															#
#####################################################################################

import numpy as np
import copy

def GetClusterPaths(neighbours,cluster_indexes,index,paths,Closed):
	# Make a hard Copy of the cluster indexes
	cluster_indexes=copy.copy(cluster_indexes)
	
	# Get Neighbours
	neigh=neighbours[index]
	cluster_indexes.append(index)
	found=0

	if len(cluster_indexes)==1:
		# print neigh
		if len(neigh)>0:
			neigh=[neigh[0]]
		
	
	# Go through all neighbours of the indexed QD
	for n in neigh:
		
		# If the point is the same as the first point in the path, close the path and save it
		if n==cluster_indexes[0] and len(cluster_indexes)>3:
			paths.append(cluster_indexes)
			Closed.append(1)

		# Apply this function iteratively to the neighbour (if we did not passed already through the point earlier)
		elif n>0 and n not in cluster_indexes and sum(Closed)<1:	
			GetClusterPaths(neighbours,cluster_indexes,n,paths,Closed)
			found=1
	# If no suitable neighbours where found, save the path
	if found ==0 and cluster_indexes not in paths:
		
		paths.append(cluster_indexes)
		Closed.append(0)
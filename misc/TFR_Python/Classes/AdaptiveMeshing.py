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
# AdaptiveMeshing Class																#
# This Class contains the algorith for the selection of the refined meshing regions.#
# The key-step is the definition of the non-convex perimeters around the single     #
# regions. 																	        #
#####################################################################################

import numpy as np
import Functions
import copy
import csv

class AdaptiveMeshing():

	def __init__(self,UI,CAE,points,displacements,triangles):
		# Class Initialization Method
		print 'Initializing Adaptive Meshing...'

		# Store internal variables
		self.UI=UI
		self.CAE=CAE
		self.points=points
		self.displacements=displacements
		self.triangles=triangles

		# Define the absolute value of the threshold from the relative one
		if UI['relative_threshold']==True:
			max_disp=0
			for u in displacements:
				max_disp=max(max_disp,np.linalg.norm(u))

			self.threshold=UI['threshold']*max_disp
		else:
			self.threshold=UI['threshold']


	def FindBoundaries(self):
		# This method find the boundary perimenter of the regions above the threshold

		map_t=[]
		points_t=[]
		U_t=[]
		conn_t=[]

		# select points over the threshold
		for i in range(len(self.points)):
			point=self.points[i]
			u=self.displacements[i]
			if np.linalg.norm(u)>self.threshold:
				map_t.append(i)
				points_t.append(point)
				U_t.append(u)
				conn_t.append([])

		# select triangles over the threshold
		triangles_t=[]
		for triangle in self.triangles:
			p1=triangle.item(0)
			p2=triangle.item(1)
			p3=triangle.item(2)
			if p1 in map_t and p2 in map_t and p3 in map_t:
				triangles_t.append(triangle)

		# build connection map for the points
		for triangle in triangles_t:
			p1=triangle.item(0)
			p1_t=map_t.index(p1)
			p2=triangle.item(1)
			p2_t=map_t.index(p2)
			p3=triangle.item(2)
			p3_t=map_t.index(p3)

			if p2 not in conn_t[p1_t]:
				conn_t[p1_t].append(p2)
			if p3 not in conn_t[p1_t]:
				conn_t[p1_t].append(p3)

			if p1 not in conn_t[p2_t]:
				conn_t[p2_t].append(p1)
			if p3 not in conn_t[p2_t]:
				conn_t[p2_t].append(p3)

			if p1 not in conn_t[p3_t]:
				conn_t[p3_t].append(p1)
			if p2 not in conn_t[p3_t]:
				conn_t[p3_t].append(p2)


		# select points at the border of the point clusters
		map_tt=[]
		points_tt=[]
		conn_tt=[]

		for i in range(len(points_t)):
			point=points_t[i]
			if len(conn_t[i]) in [2,3,4,5]:
				map_tt.append(map_t[i])
				points_tt.append(point)

		boundaries=[]
		counter=0 # safety counter

		#find boudary paths (until no points over the treshold are left)
		while len(points_tt)>2 and counter<1000:
			neigh=[]

			# Get neighbours of all points
			for i in range(len(points_tt)):
				point=points_tt[i]
				neighbours=Functions.GetPCircle( points_tt,point,self.UI['dot_distance']*1.1)
				neighbours.remove(i)
				neigh.append(neighbours)



			# Initialize Arrays for closed paths (perimenters)
			Closed=[]
			cl=[]
			paths=[]

			# Call recursive path search function
			Functions.GetClusterPaths(neigh,cl,0,paths,Closed)
			closed_paths=[]
			closed_paths_len=[]

			# Search for close paths
			for path in paths:
				neigh_start=copy.copy(neigh[path[0]])

				if path[-1] in neigh_start and len(path)>2:
					closed_paths.append(path)
					closed_paths_len.append(len(path))

			# If closed paths have been found, select the one with the longest perimenter and eliminate all connected QDs from list
			if len(closed_paths)>0:
				path=closed_paths[closed_paths_len.index(max(closed_paths_len))]


				boundary=[]
				for c in path:
					boundary.append(points_tt[c])

				points_ttt=[]
				for j in range(len(points_tt)):
					if j not in path:
						points_ttt.append(points_tt[j])

				points_tt=points_ttt

				boundaries.append(boundary)

			else:
				points_tt=points_tt[1:len(points_tt)]

			counter=counter+1
			# Continue, until all paths hae been found

		print('Found '+str(len(boundaries))+' Regions for refined meshing')
		self.boundaries=boundaries;

	def Partition(self):
		# Partition substrate  with all the found perimeters
		prt=1
		for boundary in self.boundaries:
			b=copy.copy(boundary)
			b.append(b[0])
			self.CAE.PartitionSketch(b,'Partition face-'+str(prt),self.UI['refined_mesh_size'])
			prt=prt+1


	def Mesh(self):
		# Call the CAE mesher
		self.CAE.Mesh(self.UI['mesh_growt'],self.UI['quadratic_formulation'])

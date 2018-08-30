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
# BCApplication Class																#
# This Class contains the algorith for the selection of the refined meshing regions.#
# The key-step is the definition of the non-convex perimeters around the single     #
# regions. 																	        #
#####################################################################################


import numpy as np
np.seterr(all='ignore')
import Interpolation
import Functions
import os

class BC_Application():
	def __init__(self,UI,CAE,points,displacements,triangles):
		# Class Initialization
		self.UI=UI
		self.CAE=CAE
		self.points=points
		self.displacements=displacements
		self.triangles=triangles


	def Apply(self):
		# Boundary Condition Application Method
		print('Applying boundary condions...')

		# Get displacement field size
		x_max=0
		y_max=0
		y_min=0
		x_min=0
		for xy in self.points:
			if xy[0,0]>x_max:
				x_max=xy[0,0]
			if xy[0,0]<x_min:
				x_min=xy[0,0]
			if xy[0,1]>y_max:
				y_max=xy[0,1]
			if xy[0,1]<y_min:
				y_min=xy[0,1]

		i=1
		interp=[]
		tol=1e-4

		# If reversed triangles should be ignored, search them and add them to the IgnoredRegions List
		IgnoredRegions=[]
		if self.UI['IgnoreSwitchedTriangles']==True:
			for T in self.triangles:

				# A triangle is reversed if the cross product changes sign after deformation
				V1=self.points[T[0,1]]-self.points[T[0,0]]
				V2=self.points[T[0,2]]-self.points[T[0,0]]
				A0=V1[0,0]*V2[0,1]-V1[0,1]*V2[0,0]

				d1= self.displacements[T[0,0]]
				d2= self.displacements[T[0,1]]
				d3= self.displacements[T[0,2]]
				v1=self.points[T[0,1]]+d2[0,0:2]-self.points[T[0,0]]-d1[0,0:2]
				v2=self.points[T[0,2]]+d3[0,0:2]-self.points[T[0,0]]-d1[0,0:2]
				A1=v1[0,0]*v2[0,1]-v1[0,1]*v2[0,0]

				if A1/A0<0:

					P1= self.points[T[0,0]]
					P2= self.points[T[0,1]]
					P3= self.points[T[0,2]]
					IgnoredRegions.append([np.array([P1[0,0],P1[0,1]]),np.array([P2[0,0],P2[0,1]]),np.array([P3[0,0],P3[0,1]])])



		# Prepare RBF Interpolation if needed
		if 'Interpolation' in self.UI.keys() and self.UI['Interpolation'] in ['RBF']:

			Points=np.array(self.points)
			Points=Points[:,0]

			D=np.array(self.displacements)
			D=D[:,0]

			if 'BasisFunction' in self.UI.keys():
				f=self.UI['BasisFunction']
			else:
				f='thin-plate'
			Interpolator=Interpolation.setRBF(Points,D,self.UI['dot_distance'],f)
			print ('Using RBF interpolation with a '+f+' basis function')
		else:
			print ('Using linear triangluar interpolation')

		Nlist=[]
		Nlist_free=[]
		k=1

		# Define Boundary Conditions for all nodes withing the QD Array
		for node in self.CAE.getByBoundingBox( 'nodes', x_min-tol,y_min-tol,self.UI['substrate_thickness']-tol,x_max+tol,y_max+tol,self.UI['substrate_thickness']+tol):
			xyz=node.coordinates
			# Check if the point is inside a ignored region
			if self.UI['IgnoreSwitchedTriangles']==True:
				use=True
				for T in IgnoredRegions:
					if self.PointInTriangle ([xyz[0],xyz[1]], T[0], T[1],T[2])==True:
						use=False
						break
			else:
				use=True

			# Check that the node is on the surface and shoudl not be ignored
			if abs(xyz[2]-self.UI['substrate_thickness'])<tol and use==True:

				# Get interpolated displacement vector for the node
				if 'Interpolation' in self.UI.keys() and self.UI['Interpolation'] in ['RBF']:
					u=Interpolation.RBFInterp(xyz,Interpolator)
				else:
					u=Interpolation.TriangularLinear(self.triangles,self.points,self.displacements,np.matrix(xyz))

				# If Everything is fine, apply the nodal BC to the CAE Model
				if u!=False:
					N=node.label
					interp.append([xyz[0],xyz[1],xyz[2],u[0], u[1], u[2]])
					self.CAE.NodalBC(N,'Disp'+str(i),u,self.UI['z_constr'])
					Nlist.append('Disp'+str(i))
					if i%1000==0:
						print('\t'+str(i)+ ' BCs...')
					i+=1
				else:
					Nlist_free.append('Free'+str(k))
					self.CAE.CreateSet(node.label,'Free'+str(k))
					k=k+1

		# Create CAE Sets with the constrined Nodes
		self.CAE.CreateNodalDofSet(Nlist,Nlist_free)

		print(str(i-1)+ ' Displacement BCs have been applied in total!')
		self.CAE.Regenerate()

		# Store interpolated nodal displacement field in file
		if not os.path.exists('../AppliedData'):
			os.makedirs('../AppliedData')
		Functions.CSVWrite( '../AppliedData/'+self.UI['job_name']+'_AppliedBC.txt',interp)





	def sign ( self,p1, p2,  p3):
		# Sign Method needed for the PointInTriangle Method
		return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1]- p3[1])

	def PointInTriangle (self, pt,  v1,  v2, v3):
		# This method checks if a point is inside a triangle
		b1 = self.sign(pt, v1, v2) < 0.0
		b2 = self.sign(pt, v2, v3) < 0.0
		b3 = self.sign(pt, v3, v1) < 0.0
		return ((b1 == b2) and (b2 == b3))



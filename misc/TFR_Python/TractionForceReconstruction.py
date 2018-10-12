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
# TractionForceReconstruction Class													#
# Main Class for the Nonlinear FEM part of the cTFM package							#
# Reads the displacement vector ad discrete location and returns the data for the   #
# reconstruction of the traction field												#
# The adaptive meshing algorithm and the displacement field interpolation are		#
# defined in subclasses																#
# Required: ABAQUS 6.10EF with additional python packages (see Instructions)		#
#####################################################################################

import Functions
import Classes
import Interpolation


class TractionForceReconstruction():

	def __init__(self,UI):
		# Class initialization
		self.UI=UI
		print 'Initializing Traction Force Reconstruction...'
		print('')
		print('Job Name: '+self.UI['job_name'])

		# Save variable describing the path of the displacement vector data prepared
		# with the reference configuration reconstruction tool
		self.UI['Path_PreparedData']='../PreparedData/'

		# initialize an instance of the CAE class
		self.CAE=Classes.CAE_Container('../'+self.UI['BaseModelPath'])

		# set additional variables and the define predefined options (if not already defined in the user input UI)
		self.UI['substrate_thickness']=self.CAE.thickness
		self.UI['mesh_growt']=0
		self.UI['quadratic_formulation']=0
		self.UI['BasisFunction']='thin-plate'
		self.UI['IgnoreSwitchedTriangles']=False
		UI['Interpolation']='RBF'


	def ReadInputData(self):
		# Read and store the displacement vector data
		Data=Functions.ReadInput(self.UI)
		self.points=Data[0]
		self.displacements=Data[1]
		self.triangles=Data[2]

		# Check if the Basemodel is large enough compared to the position of the displacement vectors
		Check=Functions.CheckModelSize(self.points,self.CAE.XYLim)
		if Check[0]==False:

			print('The CAE model is too small!')
			print('CAE Size:')
			print('x_min: '+str(self.CAE.XYLim[0])+'um')
			print('y_min: '+str(self.CAE.XYLim[1])+'um')
			print('x_max: '+str(self.CAE.XYLim[2])+'um')
			print('x_max: '+str(self.CAE.XYLim[3])+'um')

			print('Measurement Field:')
			print('x_min: '+str(Check[1][0])+'um')
			print('y_min: '+str(Check[1][1])+'um')
			print('x_max: '+str(Check[1][2])+'um')
			print('x_max: '+str(Check[1][3])+'um')

			raise Exception('ERROR: The Basemodel is too small!')


	def AdaptiveMesh(self):
		# Run adaptive meshing scheme
		AM=Classes.AdaptiveMeshing(self.UI,self.CAE,self.points,self.displacements,self.triangles)
		AM.FindBoundaries()
		AM.Partition()
		AM.Mesh()

	def BCApplication(self):
		# Apply displacement boundary conditions on the surface
		BC=Classes.BC_Application(self.UI,self.CAE,self.points,self.displacements,self.triangles)
		BC.Apply()

	def ExportMeshData(self):
		# Export the meshdata (needed for post processing)
		self.CAE.ExportMeshData(self.UI['job_name'],self.UI['substrate_thickness'],self.UI['quadratic_formulation'])

	def CreateJob(self):
		# Create FEM solver job
		self.CAE.CreateJob(self.UI['job_name'],self.UI['CPUs'])

	def Run(self,wait=True):
		# Run finite elelement solver
		self.CAE.Run(self.UI['job_name'],wait)

	def Save(self):
		# Save the prepared CAE model
		self.CAE.Save(self.UI['job_name'])

	def ExportResults(self,frame=-1):
		# Export the results needed for the post-processing of the traction fields
		Post=Classes.ExportResults(self.UI)
		Post.ExportUF(frame)

	def Alert(self):
		# Alert user that the FE part is finished
		raise ValueError('The Nonlinear FEM Computation has finished!\n\n Please close Abaqus CAE and use the Matlab Postprocessing script for the visualization of the traction forces!')



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
# CAE_Container Class																#
# This Class is a main wrapper for the interaction of the cTFM tool with 			#
# ABAQUS 6.10EF for the generation of the CAE file									#
#####################################################################################

from abaqusConstants import *
from odbAccess import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from __main__ import *
from abaqus import *

import Functions

class CAE_Container():
	
	def __init__(self,Basepath):
		# Class Initialization
		
		# Load Basemodel and assign internal variables
		modelDB=openMdb(Basepath)
		self.model=modelDB.models['Model-1']
		self.part=self.model.parts['Part-1']
		self.assembly=self.model.rootAssembly
		self.instance=self.assembly.instances['Part-1-1']
		session.viewports['Viewport: 1'].setValues(displayedObject=self.assembly)
		session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
		self.assembly.regenerate()
		
		# Measure the size of the Basemodel 
		x_max=x_min=y_max=y_min=z_max=0
		for v in self.part.vertices:
			if v.pointOn[0][0]>x_max:
				x_max=v.pointOn[0][0]				
			if v.pointOn[0][0]<x_min:
				x_min=v.pointOn[0][0]				
			if v.pointOn[0][1]>y_max:
				y_max=v.pointOn[0][1]				
			if v.pointOn[0][1]<y_min:
				y_min=v.pointOn[0][1]
			if v.pointOn[0][2]>z_max:
				z_max=v.pointOn[0][2]
		
		self.thickness=z_max
		self.XYLim=[x_min,y_min,x_max,y_max]		
		print ('The basemodel substrate thickness is '+str(z_max)+'um.')
		
		
		
	def Mesh(self,mesh_growth,E_type):
		# Mesh generation method
		
		# Pick Region to mesh (everything)
		pickedRegions = self.part.cells.getSequenceFromMask(mask=('[#1 ]', ), )
		
		# Delete old mesh (if the substrate is meshed)
		self.part.deleteMesh(regions=pickedRegions)

		# Assing Mesh controls,element type and mesh growth parameters
		self.part.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
		if E_type==0:
			print('Hybrid Linear Element Formulation has been used!')
			self.part.setElementType(regions=(self.part.cells, ), elemTypes=(ElemType(elemCode=C3D4H),))
		else:
			print('Hybrid Quadratic Element Formulation has been used!')
			self.part.setElementType(regions=(self.part.cells, ), elemTypes=(ElemType(elemCode=C3D10H),))
			
		if mesh_growth==1:
			self.part.setMeshControls(regions=pickedRegions, sizeGrowth=MODERATE)
		elif mesh_growth==2:
			self.part.setMeshControls(regions=pickedRegions, sizeGrowth=MAXIMUM)
		else:
			self.part.setMeshControls(regions=pickedRegions, sizeGrowth=None)
		
		# Generate and show mesh
		self.part.generateMesh()
		self.assembly.regenerate()
		session.viewports['Viewport: 1'].setValues(displayedObject=self.assembly)
		session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
		
		
	def PartitionSketch(self,coords,FN,elementsize):
		# This method partitions the substrate surface with the line defined in coords
		# and seed the newly created edges with  element size (elementsize)
		
		# Find sketch plane and the edge of the sketch orientation
		z_max=0
		for fac in self.part.faces:
			if fac.pointOn[0][2]>z_max:
				sketch_plane=fac
				z_max=fac.pointOn[0][2]
		
		sketch_edge=self.part.edges[0]
		x_max=0
		for ed in self.part.edges:
			if ed.pointOn[0][2]==z_max and ed.pointOn[0][0]>x_max:
				sketch_edge=ed
				x_max=ed.pointOn[0][0]


		# create a sketch and draw the contour give in coords
		tt = self.part.MakeSketchTransform(sketchPlane=sketch_plane, sketchUpEdge=sketch_edge, sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.5))
		s = self.model.ConstrainedSketch(name='__profile__', sheetSize=169.7, gridSpacing=4.24, transform=tt)
		g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
		s.setPrimaryObject(option=SUPERIMPOSE)

		self.part.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

		for i in range(len(coords)-1):
			s.Line(point1=(coords[i].item(0), coords[i].item(1)), point2=(coords[i+1].item(0), coords[i+1].item(1)))

		# Find all faces on the upper substrate surface
		Part_face=[]
		for fac in self.part.faces:
			if fac.pointOn[0][2]==z_max:
				Part_face.append(fac)
		
		# partittion upper face
		self.part.PartitionFaceBySketch(sketchUpEdge=sketch_edge, faces=Part_face, sketch=s)
		s.unsetPrimaryObject()

		# Seed the newly generated edges
		pickedEdges=[]
		for ed in self.part.edges:
			if ed.featureName==FN:
				pickedEdges.append(ed)
		self.part.seedEdgeBySize(edges=pickedEdges, size=elementsize, deviationFactor=0.1, constraint=FINER)

	
	
	def getByBoundingBox( self,what,xmin,ymin,zmin,xmax,ymax,zmax):
		# Wrapper for the getByBoundingBox method defined in ABAQUS
		if what=='nodes':
			return self.instance.nodes.getByBoundingBox( xmin,ymin,zmin,xmax,ymax,zmax)
		if what=='elements':
			return self.instance.elements.getByBoundingBox( xmin,ymin,zmin,xmax,ymax,zmax)
			
			
	def NodalBC(self,N,BCName,u,z_constr):
		# Define displacement boundary condition for node N
		
		SetName=BCName
		self.part.SetFromNodeLabels(SetName,(N,))
		# depending if the out of plane component is avalayble or not, define the constrain with or without it respectively
		if z_constr==False:
			self.model.DisplacementBC(name=BCName, createStepName='Step-1', region=self.instance.sets[SetName], u1=u[0], u2=u[1], u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,  fieldName='', localCsys=None)
		else:
			self.model.DisplacementBC(name=BCName, createStepName='Step-1', region=self.instance.sets[SetName], u1=u[0], u2=u[1], u3=u[2], ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,  fieldName='', localCsys=None)
		
	
	def CreateSet(self, label,SetName):
		# Wrapper for the creation of nodal sets
		self.part.SetFromNodeLabels(SetName,(label,))
	
	
	def CreateNodalDofSet(self,SetList,SetListFree):
		# Create a set of all nodes inside the perimeter of defined displacement vectors
		# Create a special list for nodes without any displacement BCs
		AsetList=[]
		for Set in SetList:
			AsetList.append(self.part.sets[Set])
		
		self.part.SetByMerge(name='Constrained', sets=tuple(AsetList) )
		
		if len(SetListFree)>0:
			AsetList_free=[]
			for Set in SetListFree:
				AsetList_free.append(self.part.sets[Set])
			
			self.part.SetByMerge(name='Free', sets=tuple(AsetList_free) )
			
			self.part.SetByMerge(name='AllNodes', sets=(self.part.sets['Constrained'],self.part.sets['Free']) )
		else:
			self.part.SetByMerge(name='AllNodes', sets=tuple(AsetList) )
		
			
		
	def Regenerate(self):
		# Wrapper for the regenration of the assembly
		self.assembly.regenerate()
		
	def ExportMeshData(self,job_name,thickness,E_type):
		# This method exports the mesh data relative to the constrained nodes on the surface
		# Besides the postion in the reference configuration of the nodes themselves, the triangles of the surface elements are also exported
		
		print 'Exporting Mesh Data....'
		
		# define lists
		points_ex=[]
		nod=[]
		node_numbers=[]
		surf_triangles=[]
		Elements=[]
		
		# regenerate assembly to be sure that the geometry is up to date
		self.assembly.regenerate()
		
		# Tolerance for the control of the z-position of the nodes
		tol=1e-4
		
		# Go through all the nodes on the substrate surface for which boundary conditions have been defined
		for node in self.instance.sets['AllNodes'].nodes:
			xyz=node.coordinates
			nod.append(node)
			node_numbers.append(node.label)
			points_ex.append([node.label,xyz[0],xyz[1],xyz[2]])
			
			# get and store the label of Element which are sharing this node
			for E in node.getElements():
				Elements.append(E.label)
		
		# Since Elements are connected to multiple nodes, the duplicates in the list have to be removed
		Elements_unique=list(set(Elements))
		
		# counter and percentage steps for the progress visualization
		counter=0
		stp=5
		
		
		
		# Go through the elements
		for e in Elements_unique:
			
			# get the element object and get the connected node objects
			el=self.instance.elements.getFromLabel(e)
			el_nodes=el.getNodes()
			
			n_in_list=[]
			
			# create label list of of the elements which are on the substrates surface
			n_list=[]
			for n in el_nodes:
				if abs(n.coordinates[2]-thickness)<tol:
					n_list.append(n.label)

			# depending on the element type, the element needs too have 3 or 6 laying nodes on the surface
			if ( E_type==0 and len(n_list)>=3) or ( E_type==1 and len(n_list)>=6):
					
				
				# intersect the list of element nodes with the list of constrained nodes and count them
				n_in_list=list(set(n_list) & set(node_numbers))
				if ( E_type==0 and len(n_in_list)==3) or ( E_type==1 and len(n_in_list)==6):
					# if the intersection of the lists returns the right number, add element face to the list
					surf_triangles.append(n_in_list)
				
			# increment counter and if needed print status
			counter+=1	
			if float(counter)/len(Elements_unique)>=stp/100.:
				print('\t'+str(stp)+ '%...')
				stp+=5
		
		# if not existent, create a folder for the meshdata
		if not os.path.exists('../MeshData'):
			os.makedirs('../MeshData')
		# save meshdata in CSV formatted text files
		Functions.CSVWrite( '../MeshData/'+job_name+'_Nodes.txt',points_ex)
		Functions.CSVWrite( '../MeshData/'+job_name+'_Elements.txt',surf_triangles)
		self.assembly.regenerate()
		
		
	def CreateJob(self,job_name,CPUs):
		# Create Job for the solution of the nonlinear problem
		mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
		explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
		memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
		multiprocessingMode=DEFAULT, name=job_name, nodalOutputPrecision=SINGLE, 
		numCpus=CPUs, numDomains=CPUs, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
		waitHours=0, waitMinutes=0)
		
	def Run(self,job_name,wait):
		# Start FE solver
		# if wait==True the script waits until the solver has finished
		print ('Run Computation...')
		mdb.jobs[job_name].submit(consistencyChecking=OFF)
		if wait==True:
			mdb.jobs[job_name].waitForCompletion()

	def Save(self,job_name):
		# Save CAE file
		mdb.saveAs(pathName=job_name)
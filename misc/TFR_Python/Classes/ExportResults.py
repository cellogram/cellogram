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
# ExportResults Class																#
# This Class contains the method for the export of the needed data from the Result  # 
# database of the NLFEM computation.											    #
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

import Functions
import csv

class ExportResults():
	
	def __init__(self,UI):
		# Initialize
		time.sleep(10)
		self.UI=UI
		
		# Load ODB 
		ODBFile=self.UI['job_name']+'.odb'
		self.odb=session.openOdb(ODBFile)
		self.part=self.odb.parts['PART-1']
		self.assembly=self.odb.rootAssembly
		self.instance=self.assembly.instances['PART-1-1']

		session.viewports['Viewport: 1'].setValues(displayedObject=self.odb)
		session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 'Magnitude'), )
		session.viewports['Viewport: 1'].view.setValues(nearPlane=293.479,  farPlane=570.218,  cameraPosition=( 0, 0, 500), cameraUpVector=(0, 1, 0))
		session.viewports['Viewport: 1'].view.fitView()
		session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=( CONTOURS_ON_DEF, ))
		session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(visibleEdges=FEATURE)
		session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
		print('Processing job: '+self.UI['job_name'])
		
	def ExportUF(self,frame=-1):
		# Export the nodal displacement and reaction forces from the ODB
		
		# Define frame and step for result export
		step='Step-1'
		appendix=''
		if frame !=-1:
			appendix='_'+str(frame)
		# Get applied Node BCs
		try:
			file=open('../AppliedData/'+self.UI['job_name']+'_AppliedBC.txt','r')
		except:
			raise Exception('../AppliedData/'+self.UI['job_name']+'_AppliedBC.txt has not been found!')
		reader=csv.reader(file,delimiter='\t')
		points=[]
	
		x_max=0
		y_max=0
		y_min=0
		x_min=0
		
		tol=1e-4
		
		# Get geometric boundaries of the region with boundary conditions
		for row in reader:
			if float(row[0])>x_max:
				x_max=float(row[0])				
			if float(row[0])<x_min:
				x_min=float(row[0])				
			if float(row[1])>y_max:
				y_max=float(row[1])			
			if float(row[1])<y_min:
				y_min=float(row[1])
		t= self.UI['substrate_thickness']
		
		# Extract Strain Energy
		StrainEnergy=self.odb.steps[step].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLSE'].data[frame][1]
		
		# Get Displacement and reaction force Field
		Ffield=self.odb.steps[step].frames[frame].fieldOutputs['RT'].getSubset(region=self.instance.nodeSets['ALLNODES'])
		Ufield=self.odb.steps[step].frames[frame].fieldOutputs['U'].getSubset(region=self.instance.nodeSets['ALLNODES'])
		Ffield_FE=[]
		
		# Combine displacement and force Data
		for u in range(len(Ffield.values)):
			
			F=Ffield.values[u]
			U=Ufield.values[u]
			node=self.instance.getNodeFromLabel(F.nodeLabel)
			xyz=node.coordinates
			
			if F.instance==self.instance and abs(xyz[2]-t)<tol and  xyz[0]>x_min-tol and xyz[0]<x_max+tol and xyz[1]>y_min-tol and xyz[1]<y_max+tol :
				Ffield_FE.append([node.label,U.data[0],U.data[1],U.data[2],F.data[0],F.data[1],F.data[2]])

		
		if not os.path.exists('../ReconstructedData'):
			os.makedirs('../ReconstructedData')
		if not os.path.exists('../ReconstructedData/ReactionForces'):
			os.makedirs('../ReconstructedData/ReactionForces')
		if not os.path.exists('../ReconstructedData/HistoryOutputs'):
			os.makedirs('../ReconstructedData/HistoryOutputs')
			
		Functions.CSVWrite( '../ReconstructedData/ReactionForces/'+self.UI['job_name']+appendix+'.txt',Ffield_FE)
		
		StrainEnergy=self.odb.steps[step].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLSE'].data[-1][1]
		Functions.CSVWrite( '../ReconstructedData/HistoryOutputs/'+self.UI['job_name']+appendix+'.txt',[[StrainEnergy,'% Strain Energy']])
		
		time.sleep(1)
		print('Result Export finished')
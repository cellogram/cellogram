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
# Read all the necessary inputs for the NLFEM part from the reference position 		#
# estimation step of the cTFM package												#
#####################################################################################

import csv
import numpy as np

	
def ReadInput(UI):

	# Read the data for the QDs position, the relative displacement vectors and the triangulation data
	
	try:
		file=open(UI['Path_PreparedData']+UI['job_name']+'_points.txt','r')
	except:
		raise Exception(UI['Path_PreparedData']+UI['job_name']+'_points.txt has not been found!')
	reader=csv.reader(file,delimiter='\t')
	points=[]
	for row in reader:
		points.append(np.matrix([float(row[0]),float(row[1])]))
	
	try:
		file=open(UI['Path_PreparedData']+UI['job_name']+'_displacement.txt','r')
	except:
		raise Exception(UI['Path_PreparedData']+UI['job_name']+'_displacement.txt has not been found!')
	reader=csv.reader(file,delimiter='\t')
	U=[]
	for row in reader:
		if len(row)<3 and UI['z_constr']==True:
			UI['z_constr']=False
			print('Warning: Since the z-component of displacement is not given, the parameter z_constr will be set to false!')
	
		if UI['z_constr']==False:
			U.append(np.matrix([float(row[0]),float(row[1]),float(0)]))
		else:
			U.append(np.matrix([float(row[0]),float(row[1]),float(row[2])]))


	try:
		file=open(UI['Path_PreparedData']+UI['job_name']+'_triangles.txt','r')
	except:
		raise Exception(UI['Path_PreparedData']+UI['job_name']+'_triangles.txt has not been found!')	
	triangles=[]
	reader=csv.reader(file,delimiter='\t')
	for row in reader:
		triangles.append(np.matrix([int(row[0]),int(row[1]),int(row[2])]))


	print('Found ' + str(len(U))+ ' displacement boundary conditions for QDs')
	
	# Get additional informations from the optimization step
	file=open(UI['Path_PreparedData']+UI['job_name']+'_UserDataABQ.txt','r')
	reader=csv.reader(file,delimiter='\t')
	for row in reader:
		UI['image_scaling']=float(row[0])
		UI['dot_distance']=float(row[1])
		UI['imageSize']=[float(row[3]),float(row[2])]
		break # File has only one row

	return [points,U,triangles]	
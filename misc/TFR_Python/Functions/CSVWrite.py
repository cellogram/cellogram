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

# CSV File Writer Method

import csv
	
def CSVWrite( filename,data ):
	file=open(filename,'wb')
	writer=csv.writer(file,delimiter='\t')
	for d in data:
		writer.writerow(d)
	file.close()
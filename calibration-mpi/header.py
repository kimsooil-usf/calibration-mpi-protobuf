import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import os
import time
import scipy.stats as stats
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import scipy.integrate as integrate
import datetime as dt
import math
from pyDOE import *
import argparse
import random

np.set_printoptions(threshold=np.inf)
DATE_TO_RUN="1110-50k-960days"
# DATE_TO_RUN="1110-final-test-for-960d-onepiece"
INPUT_DIR = "/largedisk/mpi-folder2/calibration-mpi-protobuf/staticInst/output/hills-1.45m"
# INPUT_DIR = "/largedisk/mpi-folder2/calibration-mpi-protobuf/staticInst/output/hills-150k" # small pop for testing
OUTPUT_DIR =  "/largedisk2/mpi-folder/output-mpi-"+DATE_TO_RUN+"/"  # should be same as MPI_DIR # also use absolute path
MPI_DIR = "/largedisk2/mpi-folder/output-mpi-"+DATE_TO_RUN+"/"  # should be same as OUTPUT_DIR
PLOT_DIR = "/largedisk2/mpi-folder/output-mpi-"+DATE_TO_RUN+"/plots/"

if not os.path.exists(OUTPUT_DIR):
        os.system("mkdir -p " + OUTPUT_DIR)

PIECE = 1 # 32
# NUM_DAYS = 30
NUM_DAYS = 960

NRMSE=500 # How many RMSE indices to be used in the next piece (Top NRMSE smallest RMSE values)

# NPARAMS = 400
# NPROCESSORS = 360
#NPARAMS = 15

NPARAMS = 50000 # 20250=810x25 # Must be even number and multiple of NPROCESSORS, no error. (SK & JK 10/19)
NPROCESSORS = 1500
NPERNODE = 90 ###### NOT using at the most recent simulations 

# NPARAMS = 4 # 20250=810x25 # Must be even number and multiple of NPROCESSORS, no error. (SK & JK 10/19)
# NPROCESSORS = 4
# NPERNODE = 1 ###### NOT using at the most recent simulations 

# Initialize minValues, maxValues
# Later used as the parameter for drive_simulator in run_simulator.py
#region
#-----------------------------------------------------------------------------#
minValues={}
maxValues={}
#-----------------------------------------------------------------------------#

#-------------Initial fraction of infected/fixed number of infected-----------#
minValues['INIT_FIXED_NUMBER_INFECTED']=1
maxValues['INIT_FIXED_NUMBER_INFECTED']=500
#-----------------------------------------------------------------------------#

#-------------------INFECTIOUSNESS range is set in this section---------------#
minValues['INFECTIOUSNESS_ALPHA']=1
maxValues['INFECTIOUSNESS_ALPHA']=2

minValues['INFECTIOUSNESS_DELTA']=1.5
maxValues['INFECTIOUSNESS_DELTA']=4

minValues['INFECTIOUSNESS_OMICRON']=2
maxValues['INFECTIOUSNESS_OMICRON']=5


minValues['INFECTIOUSNESS_OMICRON_NEW']=2
maxValues['INFECTIOUSNESS_OMICRON_NEW']=6
#-----------------------------------------------------------------------------#

#---------------VIRULENCE values for strains----------------------------------#
minValues['VIRULENT_NEW_ALPHA']=1
maxValues['VIRULENT_NEW_ALPHA']=1.5

minValues['VIRULENT_NEW_DELTA']=1.5
maxValues['VIRULENT_NEW_DELTA']=3

minValues['VIRULENT_NEW_OMICRON']=0.3
maxValues['VIRULENT_NEW_OMICRON']=1

minValues['VIRULENT_NEW_OMICRON_NEW']=.3
maxValues['VIRULENT_NEW_OMICRON_NEW']=1
#-----------------------------------------------------------------------------#

#-----------Proportion of reinfections from different strains-----------------#
# minValues['REINFECTION_ALPHA']=0.001
# maxValues['REINFECTION_ALPHA']=0.1

# minValues['REINFECTION_DELTA']=0.0001
# maxValues['REINFECTION_DELTA']=0.2

# minValues['REINFECTION_OMICRON']=0.0001
# maxValues['REINFECTION_OMICRON']=0.3

# minValues['REINFECTION_OMICRON_NEW']=0.0001
# maxValues['REINFECTION_OMICRON_NEW']=0.4
#-----------------------------------------------------------------------------#

#--------------Fraction of new strains----------------------------------------#
minValues['FRACTION_NEW_ALPHA']=0.0001
maxValues['FRACTION_NEW_ALPHA']=0.1

minValues['FRACTION_NEW_DELTA']=0.0001
maxValues['FRACTION_NEW_DELTA']=0.1

minValues['FRACTION_NEW_OMICRON']=0.0001
maxValues['FRACTION_NEW_OMICRON']=0.1

minValues['FRACTION_NEW_OMICRON_NEW']=0.0001
maxValues['FRACTION_NEW_OMICRON_NEW']=0.1
#-----------------------------------------------------------------------------#

#----------------Fraction of susceptibles-------------------------------------#
# minValues['FRACTION_SUSCEPTIBLE_ALPHA']=0.0001
# maxValues['FRACTION_SUSCEPTIBLE_ALPHA']=1

# minValues['FRACTION_SUSCEPTIBLE_DELTA']=0.0001
# maxValues['FRACTION_SUSCEPTIBLE_DELTA']=1

# minValues['FRACTION_SUSCEPTIBLE_OMICRON']=0.0001
# maxValues['FRACTION_SUSCEPTIBLE_OMICRON']=1

# minValues['FRACTION_SUSCEPTIBLE_OMICRON_NEW']=0.0001
# maxValues['FRACTION_SUSCEPTIBLE_OMICRON_NEW']=1
#-----------------------------------------------------------------------------#

#----------------Vaccination Effectiveness------------------------------------#
minValues['VACCINATION_EFFECTIVENESS1']=0.5
maxValues['VACCINATION_EFFECTIVENESS1']=.9

minValues['VACCINATION_EFFECTIVENESS2']=0.5
maxValues['VACCINATION_EFFECTIVENESS2']=1

minValues['VACCINATION_EFFECTIVENESS_WANING']=0.5
maxValues['VACCINATION_EFFECTIVENESS_WANING']=1

minValues['VACCINATION_EFFECTIVENESS_BOOSTED']=0.8
maxValues['VACCINATION_EFFECTIVENESS_BOOSTED']=1
#-----------------------------------------------------------------------------#

#------------------Beta: home, school, work, community------------------------#
minValues['BETA_SCALE']=1
maxValues['BETA_SCALE']=9

minValues['BETA_H']=0.0001
maxValues['BETA_H']=1#9

minValues["BETA_RANDOM_COMMUNITY"]=0.0001#BETA_NBR_CELLS
maxValues["BETA_RANDOM_COMMUNITY"]=1#BETA_NBR_CELLS

minValues["BETA_W"]=0.0001
maxValues["BETA_W"]=1

minValues["BETA_S"]=0.001
maxValues["BETA_S"]=1.5

minValues["BETA_C"]=0.001
maxValues["BETA_C"]=1

minValues['BETA_PROJECT']=0.0001
maxValues['BETA_PROJECT']=1

minValues['BETA_NBR_CELLS']=0.0001
maxValues['BETA_NBR_CELLS']=1

minValues['BETA_CLASS']=0.0001
maxValues['BETA_CLASS']=1

minValues['BETA_TRAVEL']=0.0001
maxValues['BETA_TRAVEL']=1
#-----------------------------------------------------------------------------#

#---------------Other Parameters----------------------------------------------#
minValues['CALIBRATION_DELAY']=1
maxValues['CALIBRATION_DELAY']=13

minValues['DAYS_BEFORE_LOCKDOWN']=1
maxValues['DAYS_BEFORE_LOCKDOWN']=13
#-----------------------------------------------------------------------------#

#----------------Time of new variants-----------------------------------------#
minValues['TIME_ALPHA']=(dt.date(2020,12,30)-dt.date(2020,3,1)).days-45
maxValues['TIME_ALPHA']=(dt.date(2020,12,30)-dt.date(2020,3,1)).days+10

minValues['TIME_DELTA']=(dt.date(2021,5,2)-dt.date(2020,3,1)).days-45
maxValues['TIME_DELTA']=(dt.date(2021,5,2)-dt.date(2020,3,1)).days+10

minValues['TIME_OMICRON']=(dt.date(2021,12,11)-dt.date(2020,3,1)).days-45
maxValues['TIME_OMICRON']=(dt.date(2021,12,11)-dt.date(2020,3,1)).days+10

minValues['TIME_OMICRON_NEW']=(dt.date(2022,4,22)-dt.date(2020,3,1)).days-45
maxValues['TIME_OMICRON_NEW']=(dt.date(2022,4,22)-dt.date(2020,3,1)).days+10
#-----------------------------------------------------------------------------#

# Change
#----------------Random_seeds-----------------------------------------#
minValues['PROVIDE_INITIAL_SEED_GRAPH']=1
maxValues['PROVIDE_INITIAL_SEED_GRAPH']=NPARAMS

minValues['PROVIDE_INITIAL_SEED']=1
maxValues['PROVIDE_INITIAL_SEED']=NPARAMS
#-----------------------------------------------------------------------------#
#endregion

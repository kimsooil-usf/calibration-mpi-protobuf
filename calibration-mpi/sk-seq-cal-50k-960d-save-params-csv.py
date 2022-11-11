import os
import numpy as np
import pandas as pd
import time
import scipy.stats as stats
import datetime as dt
import datetime as dt
import math
from pyDOE import *
import argparse
import random

from multiprocessing import Process, Queue, Manager
import multiprocessing as mp

output_directory_base0 =  "../../temp/seq-cal-output-1110"
input_directory = "../staticInst/output/hills-1.45m/"
NUM_DAYS=960
NPARAMS=50000

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

############################################        

def priors(nParams):
    rangeParam = {}
    for ld in minValues:
        rangeParam[ld] = maxValues[ld] - minValues[ld]

    np.random.seed(2)#fixing  
    paraLHS=lhs(len(minValues), samples=nParams, criterion='corr')

    i = 0
    VarParamsLHS=pd.DataFrame()
    for ld in minValues:
        VarParamsLHS[ld] = minValues[ld]+rangeParam[ld]*paraLHS[:,i]
        i+=1

    VarParamsLHS['PROVIDE_INITIAL_SEED_GRAPH']=random.sample(range(1, nParams+1), nParams) # Change
    VarParamsLHS['PROVIDE_INITIAL_SEED']=random.sample(range(1, nParams+1), nParams) # Change
    # VarParamsLHS['PROVIDE_INITIAL_SEED_GRAPH']=1
    # VarParamsLHS['PROVIDE_INITIAL_SEED']=1
    return VarParamsLHS

def run_parallel_sim(parameters):
    # print("in run_parallel_sim")
    # print(len(parameters))
    #index=[x for x in range(len(parameters))]

    pool = mp.Pool(4)#processes=min((os.cpu_count() - 1), PROC_TO_RUN))

    for i in range(len(parameters)):
        pool.apply_async(run_parameter_simulator, (parameters,i))
    pool.close()
    pool.join()

def getCommand(params):
    # print("in getCommand")
    comms=[k for k in params]
    
    comms.pop(0)
    comms.pop(0)

    command = " "
    command += params["execDir"] + "drive_simulator" + " "
    command += params["SEED_FIXED_NUMBER"]

    for cm in comms:
        command+=" "+" --" + cm +" " + str(params[cm])

    return command

def run_parameter_simulator(params, rank):
    #print("in run_param...", rank)

    output_directory_base =  params["output_directory_base"][rank]
    output_directory =  output_directory_base +str(rank)+"/"
    # print("output_directory", output_directory)
    # if not os.path.exists(output_directory):
    #     os.system("mkdir -p " + output_directory)

    params_to_send = {
        "execDir": "../cpp-simulator/",

        "SEED_FIXED_NUMBER": "--SEED_FIXED_NUMBER",
        "NUM_DAYS": NUM_DAYS,
        "F_KERNEL_A": 10.751,
        "F_KERNEL_B": 5.384,

        "INTERVENTION": 18, #----Intervention 18 inlcudes hillsborough specific interventions

        "INIT_FIXED_NUMBER_INFECTED": int(float(params["INIT_FIXED_NUMBER_INFECTED"][rank])),

        "INFECTIOUSNESS_ALPHA":params["INFECTIOUSNESS_ALPHA"][rank],
        "INFECTIOUSNESS_DELTA":params["INFECTIOUSNESS_DELTA"][rank],
        "INFECTIOUSNESS_OMICRON":params["INFECTIOUSNESS_OMICRON"][rank],
        "INFECTIOUSNESS_OMICRON_NEW":params["INFECTIOUSNESS_OMICRON_NEW"][rank],

        "VIRULENT_NEW_ALPHA":params["VIRULENT_NEW_ALPHA"][rank],
        "VIRULENT_NEW_DELTA":params["VIRULENT_NEW_DELTA"][rank],
        "VIRULENT_NEW_OMICRON":params["VIRULENT_NEW_OMICRON"][rank],
        "VIRULENT_NEW_OMICRON_NEW":params["VIRULENT_NEW_OMICRON_NEW"][rank],

        # "REINFECTION_ALPHA":params["REINFECTION_ALPHA"][rank],
        # "REINFECTION_DELTA":params["REINFECTION_DELTA"][rank],
        # "REINFECTION_OMICRON":params["REINFECTION_OMICRON"][rank],
        # "REINFECTION_OMICRON_NEW":params["REINFECTION_OMICRON_NEW"][rank],

        "FRACTION_NEW_ALPHA":params["FRACTION_NEW_ALPHA"][rank],
        "FRACTION_NEW_DELTA":params["FRACTION_NEW_DELTA"][rank],
        "FRACTION_NEW_OMICRON":params["FRACTION_NEW_OMICRON"][rank],
        "FRACTION_NEW_OMICRON_NEW":params["FRACTION_NEW_OMICRON_NEW"][rank],
        
        # "FRACTION_SUSCEPTIBLE_ALPHA":params["FRACTION_SUSCEPTIBLE_ALPHA"][rank],
        # "FRACTION_SUSCEPTIBLE_DELTA":params["FRACTION_SUSCEPTIBLE_DELTA"][rank],
        # "FRACTION_SUSCEPTIBLE_OMICRON":params["FRACTION_SUSCEPTIBLE_OMICRON"][rank],
        # "FRACTION_SUSCEPTIBLE_OMICRON_NEW":params["FRACTION_SUSCEPTIBLE_OMICRON_NEW"][rank],

        "VACCINATION_EFFECTIVENESS1":params["VACCINATION_EFFECTIVENESS1"][rank],
        "VACCINATION_EFFECTIVENESS2":params["VACCINATION_EFFECTIVENESS2"][rank],
        "VACCINATION_EFFECTIVENESS_WANING":params["VACCINATION_EFFECTIVENESS_WANING"][rank],
        "VACCINATION_EFFECTIVENESS_BOOSTED":params["VACCINATION_EFFECTIVENESS_BOOSTED"][rank],

        "BETA_H":params["BETA_H"][rank],
        "BETA_RANDOM_COMMUNITY":params["BETA_RANDOM_COMMUNITY"][rank],
        "BETA_W":params[ "BETA_W"][rank],
        "BETA_S":params["BETA_S"][rank],
        "BETA_C":params["BETA_C"][rank],
        "BETA_PROJECT":params["BETA_PROJECT"][rank],
        "BETA_NBR_CELLS":params["BETA_NBR_CELLS"][rank],
        "BETA_CLASS":params["BETA_CLASS"][rank],
        "BETA_TRAVEL":params["BETA_TRAVEL"][rank],

        "CALIBRATION_DELAY": int(float(params["CALIBRATION_DELAY"][rank])),
        "DAYS_BEFORE_LOCKDOWN": int(float(params["DAYS_BEFORE_LOCKDOWN"][rank])),

        "TIME_ALPHA":int(float(params["TIME_DELTA"][rank])),
        "TIME_DELTA":int(float(params["TIME_DELTA"][rank])),
        "TIME_OMICRON":int(float(params["TIME_OMICRON"][rank])),
        "TIME_OMICRON_NEW":int(float(params["TIME_OMICRON_NEW"][rank])),

        "PROVIDE_INITIAL_SEED_GRAPH":params["PROVIDE_INITIAL_SEED_GRAPH"][rank],
        "PROVIDE_INITIAL_SEED":params["PROVIDE_INITIAL_SEED"][rank],
        # "PROVIDE_INITIAL_SEED_GRAPH":int(seed2),
        # "PROVIDE_INITIAL_SEED":int(seed1),

        "START_DAY": params["START_DAY"][rank],
        "output_directory": output_directory,
        "input_directory": params["input_directory"][rank],
        # "START_DAY": int(start_day),
        # "output_directory": out_dir,
        # "input_directory": in_dir,           
    }
    # print(params_to_send)
    command=getCommand(params_to_send)

    # print(command)

    ####os.system(command)
    
    #return command

############### main() ####################

for nsd in range(1,2):
    # create folder if not exist
    output_directory_base =  output_directory_base0 +"/piece_"+str(nsd)+"/"
    if not os.path.exists(output_directory_base):
        os.system("mkdir -p " + output_directory_base)
#---------Selecting priors sequentially-------------------------
#--------In the first piece of the simulation, the priors are selected from a fixed range as determined from literature/guess work.
#------- While the subsequent pieces the priors are informed from the best fitting parameters obtained in the previous step.
    n_Params = NPARAMS
    # if nsd>1:
    #     posterior0,posterior_index,prior=Priors_history(n_Params,nsd,output_directory_base0)
    if nsd==1:
        parameters = priors(n_Params)
    # else:
        # par1=Priors(int(n_Params/2))
        # par2=Priors_history(int(n_Params/2), nsd, output_directory_base0,ChiCrit)
        # parameters = pd.concat([par1,par2],ignore_index=True)#this function uses the information on current simulation piece, 
        # #and based on that it digs up the previous simulation results, calcuated best fitting parameters, 
        # #checks their range and then generates a new prior based on this new range.
        

    parameters["START_DAY"]=0
    parameters["output_directory_base"]=output_directory_base
    parameters["input_directory"]=input_directory
    #print("parameters length", len(parameters))

    parameters.to_csv(os.path.join(output_directory_base,"prior_parameters_sequential_"+str(nsd)+".csv"))
    
    run_parallel_sim(parameters)
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

NUM_DAYS=960
NUM_CORES=90

def run_parallel_sim(parameters):

    pool = mp.Pool(NUM_CORES)#processes=min((os.cpu_count() - 1), PROC_TO_RUN))

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
    print("output_directory", output_directory)
    if not os.path.exists(output_directory):
        os.system("mkdir -p " + output_directory)

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

    #print(command)

    #os.system(command)
    
    #return command

############### main() ####################
import sys

def main(argv):
    #print(argv) read full path of parameter csv file
    parameters=pd.read_csv(argv)
    #print(parameters)

    run_parallel_sim(parameters)

if __name__ == "__main__":

   main(sys.argv[1])
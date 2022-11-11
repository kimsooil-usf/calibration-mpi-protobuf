from header import *

def getCommand(params):
    comms=[k for k in params]
    
    comms.pop(0)
    comms.pop(0)

    command = " "
    command += params["execDir"] + "drive_simulator" + " "
    command += params["SEED_FIXED_NUMBER"]

    for cm in comms:
        command+=" "+" --" + cm +" " + str(params[cm])

    return command

def run_parameter_simulator(params, rank, seed1, seed2, start_day, in_dir, out_dir):
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
        "INFECTIOUSNESS_OMICRON_BA4":params["INFECTIOUSNESS_OMICRON_BA4"][rank],
        "INFECTIOUSNESS_OMICRON_BA5":params["INFECTIOUSNESS_OMICRON_BA5"][rank],
        "INFECTIOUSNESS_OMICRON_BA6":params["INFECTIOUSNESS_OMICRON_BA6"][rank],
        "INFECTIOUSNESS_OMICRON_BA7":params["INFECTIOUSNESS_OMICRON_BA7"][rank],
        "INFECTIOUSNESS_OMICRON_BA8":params["INFECTIOUSNESS_OMICRON_BA8"][rank],

        "VIRULENT_NEW_ALPHA":params["VIRULENT_NEW_ALPHA"][rank],
        "VIRULENT_NEW_DELTA":params["VIRULENT_NEW_DELTA"][rank],
        "VIRULENT_NEW_OMICRON":params["VIRULENT_NEW_OMICRON"][rank],
        "VIRULENT_NEW_OMICRON_NEW":params["VIRULENT_NEW_OMICRON_NEW"][rank],
        "VIRULENT_NEW_OMICRON_BA4":params["VIRULENT_NEW_OMICRON_BA4"][rank],
        "VIRULENT_NEW_OMICRON_BA5":params["VIRULENT_NEW_OMICRON_BA5"][rank],
        "VIRULENT_NEW_OMICRON_BA6":params["VIRULENT_NEW_OMICRON_BA6"][rank],
        "VIRULENT_NEW_OMICRON_BA7":params["VIRULENT_NEW_OMICRON_BA7"][rank],
        "VIRULENT_NEW_OMICRON_BA8":params["VIRULENT_NEW_OMICRON_BA8"][rank],

        # "REINFECTION_ALPHA":params["REINFECTION_ALPHA"][rank],
        # "REINFECTION_DELTA":params["REINFECTION_DELTA"][rank],
        # "REINFECTION_OMICRON":params["REINFECTION_OMICRON"][rank],
        # "REINFECTION_OMICRON_NEW":params["REINFECTION_OMICRON_NEW"][rank],

        "FRACTION_NEW_ALPHA":params["FRACTION_NEW_ALPHA"][rank],
        "FRACTION_NEW_DELTA":params["FRACTION_NEW_DELTA"][rank],
        "FRACTION_NEW_OMICRON":params["FRACTION_NEW_OMICRON"][rank],
        "FRACTION_NEW_OMICRON_NEW":params["FRACTION_NEW_OMICRON_NEW"][rank],
        "FRACTION_NEW_OMICRON_BA4":params["FRACTION_NEW_OMICRON_BA4"][rank],
        "FRACTION_NEW_OMICRON_BA5":params["FRACTION_NEW_OMICRON_BA5"][rank],
        "FRACTION_NEW_OMICRON_BA6":params["FRACTION_NEW_OMICRON_BA6"][rank],
        "FRACTION_NEW_OMICRON_BA7":params["FRACTION_NEW_OMICRON_BA7"][rank],
        "FRACTION_NEW_OMICRON_BA8":params["FRACTION_NEW_OMICRON_BA8"][rank],
        
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
        "TIME_OMICRON_BA4":int(float(params["TIME_OMICRON_BA4"][rank])),
        "TIME_OMICRON_BA5":int(float(params["TIME_OMICRON_BA5"][rank])),
        "TIME_OMICRON_BA6":int(float(params["TIME_OMICRON_BA6"][rank])),
        "TIME_OMICRON_BA7":int(float(params["TIME_OMICRON_BA7"][rank])),
        "TIME_OMICRON_BA8":int(float(params["TIME_OMICRON_BA8"][rank])),

        # "PROVIDE_INITIAL_SEED_GRAPH":params[rank][41],
        # "PROVIDE_INITIAL_SEED":params[rank][42],
        "PROVIDE_INITIAL_SEED_GRAPH":int(seed2),
        "PROVIDE_INITIAL_SEED":int(seed1),

        # "START_DAY": params[rank][43],
        # "output_directory": str(params[rank][44]) + str(rank) + "/",
        # "input_directory": params[rank][45],
        "START_DAY": int(start_day),
        "output_directory": out_dir,
        "input_directory": in_dir,           
    }
    if start_day>0:
        params_to_send["CALIBRATION_DELAY"]=0
        params_to_send["DAYS_BEFORE_LOCKDOWN"]=0
        
    command=getCommand(params_to_send)
    
    return command

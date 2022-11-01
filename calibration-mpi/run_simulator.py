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

def run_parameter_simulator(parameters, rank, seed1, seed2, start_day, in_dir, out_dir):
    params = {
        "execDir": "../cpp-simulator/",

        "SEED_FIXED_NUMBER": "--SEED_FIXED_NUMBER",
        "NUM_DAYS": NUM_DAYS,
        "F_KERNEL_A": 10.751,
        "F_KERNEL_B": 5.384,

        "INTERVENTION": 18, #----Intervention 18 inlcudes hillsborough specific interventions

        "INIT_FIXED_NUMBER_INFECTED": int(float(parameters[rank][0])),

        "INFECTIOUSNESS_ALPHA":parameters[rank][1],
        "INFECTIOUSNESS_DELTA":parameters[rank][2],
        "INFECTIOUSNESS_OMICRON":parameters[rank][3],
        "INFECTIOUSNESS_OMICRON_NEW":parameters[rank][4],

        "VIRULENT_NEW_ALPHA":parameters[rank][5],
        "VIRULENT_NEW_DELTA":parameters[rank][6],
        "VIRULENT_NEW_OMICRON":parameters[rank][7],
        "VIRULENT_NEW_OMICRON_NEW":parameters[rank][8],

        "REINFECTION_ALPHA":parameters[rank][9],
        "REINFECTION_DELTA":parameters[rank][10],
        "REINFECTION_OMICRON":parameters[rank][11],
        "REINFECTION_OMICRON_NEW":parameters[rank][12],

        "FRACTION_NEW_ALPHA":parameters[rank][13],
        "FRACTION_NEW_DELTA":parameters[rank][14],
        "FRACTION_NEW_OMICRON":parameters[rank][15],
        "FRACTION_NEW_OMICRON_NEW":parameters[rank][16],
        
        "FRACTION_SUSCEPTIBLE_ALPHA":parameters[rank][17],
        "FRACTION_SUSCEPTIBLE_DELTA":parameters[rank][18],
        "FRACTION_SUSCEPTIBLE_OMICRON":parameters[rank][19],
        "FRACTION_SUSCEPTIBLE_OMICRON_NEW":parameters[rank][20],

        "VACCINATION_EFFECTIVENESS1":parameters[rank][21],
        "VACCINATION_EFFECTIVENESS2":parameters[rank][22],
        "VACCINATION_EFFECTIVENESS_WANING":parameters[rank][23],
        "VACCINATION_EFFECTIVENESS_BOOSTED":parameters[rank][24],

        "BETA_H":parameters[rank][26],
        "BETA_RANDOM_COMMUNITY":parameters[rank][27],
        "BETA_W":parameters[rank][28],
        "BETA_S":parameters[rank][29],
        "BETA_C":parameters[rank][30],
        "BETA_PROJECT":parameters[rank][31],
        "BETA_NBR_CELLS":parameters[rank][32],
        "BETA_CLASS":parameters[rank][33],
        "BETA_TRAVEL":parameters[rank][34],

        "CALIBRATION_DELAY": int(float(parameters[rank][35])),
        "DAYS_BEFORE_LOCKDOWN": int(float(parameters[rank][36])),

        "TIME_ALPHA":int(float(parameters[rank][37])),
        "TIME_DELTA":int(float(parameters[rank][38])),
        "TIME_OMICRON":int(float(parameters[rank][39])),
        "TIME_OMICRON_NEW":int(float(parameters[rank][40])),

        # "PROVIDE_INITIAL_SEED_GRAPH":parameters[rank][41],
        # "PROVIDE_INITIAL_SEED":parameters[rank][42],
        "PROVIDE_INITIAL_SEED_GRAPH":int(seed2),
        "PROVIDE_INITIAL_SEED":int(seed1),

        # "START_DAY": params[rank][43],
        # "output_directory": str(params[rank][44]) + str(rank) + "/",
        # "input_directory": params[rank][45],
        "START_DAY": int(start_day),
        "output_directory": out_dir,
        "input_directory": in_dir,   
    }
    
    command=getCommand(params)

    return command

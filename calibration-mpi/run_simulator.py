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

        "INIT_FIXED_NUMBER_INFECTED": int(float(params[rank][0])),

        "INFECTIOUSNESS_ALPHA":params[rank][1],
        "INFECTIOUSNESS_DELTA":params[rank][2],
        "INFECTIOUSNESS_OMICRON":params[rank][3],
        "INFECTIOUSNESS_OMICRON_NEW":params[rank][4],

        "VIRULENT_NEW_ALPHA":params[rank][5],
        "VIRULENT_NEW_DELTA":params[rank][6],
        "VIRULENT_NEW_OMICRON":params[rank][7],
        "VIRULENT_NEW_OMICRON_NEW":params[rank][8],

        "REINFECTION_ALPHA":params[rank][9],
        "REINFECTION_DELTA":params[rank][10],
        "REINFECTION_OMICRON":params[rank][11],
        "REINFECTION_OMICRON_NEW":params[rank][12],

        "FRACTION_NEW_ALPHA":params[rank][13],
        "FRACTION_NEW_DELTA":params[rank][14],
        "FRACTION_NEW_OMICRON":params[rank][15],
        "FRACTION_NEW_OMICRON_NEW":params[rank][16],
        
        "FRACTION_SUSCEPTIBLE_ALPHA":params[rank][17],
        "FRACTION_SUSCEPTIBLE_DELTA":params[rank][18],
        "FRACTION_SUSCEPTIBLE_OMICRON":params[rank][19],
        "FRACTION_SUSCEPTIBLE_OMICRON_NEW":params[rank][20],

        "VACCINATION_EFFECTIVENESS1":params[rank][21],
        "VACCINATION_EFFECTIVENESS2":params[rank][22],
        "VACCINATION_EFFECTIVENESS_WANING":params[rank][23],
        "VACCINATION_EFFECTIVENESS_BOOSTED":params[rank][24],

        "BETA_H":params[rank][26],
        "BETA_RANDOM_COMMUNITY":params[rank][27],
        "BETA_W":params[rank][28],
        "BETA_S":params[rank][29],
        "BETA_C":params[rank][30],
        "BETA_PROJECT":params[rank][31],
        "BETA_NBR_CELLS":params[rank][32],
        "BETA_CLASS":params[rank][33],
        "BETA_TRAVEL":params[rank][34],

        "CALIBRATION_DELAY": int(float(params[rank][35])),
        "DAYS_BEFORE_LOCKDOWN": int(float(params[rank][36])),

        "TIME_ALPHA":int(float(params[rank][37])),
        "TIME_DELTA":int(float(params[rank][38])),
        "TIME_OMICRON":int(float(params[rank][39])),
        "TIME_OMICRON_NEW":int(float(params[rank][40])),

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
    
    command=getCommand(params_to_send)
    
    return command

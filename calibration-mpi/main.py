from header import *
from calculations import rmse_save_posterior

def generateCase():
    for nsd in range(1,PIECE+1):
        piece_directory =  OUTPUT_DIR +"piece_"+str(nsd)+"/"
        if not os.path.exists(piece_directory):
            os.system("mkdir -p " + piece_directory)

        timerange=list(np.arange(1+nsd*30, 1+nsd*30 + NUM_DAYS, 1))

        dataOriginal = pd.read_csv(os.path.join(INPUT_DIR,'HillsCases_deaths_hospitals_March1_Octover26.csv'),parse_dates=['date'])

        data = dataOriginal.loc[dataOriginal.time_step.isin(timerange)].reset_index()
        data = data.drop(['index'],axis=1)
        data.to_csv(os.path.join(piece_directory,'casesfile_'+str(nsd)+'.csv'),index=False)

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

def priors_history(nParams, nsd):
    piece_directory =  OUTPUT_DIR + "piece_"+str(nsd-1) + "/"
    #parameters_ran = [0] * nParams * 2
    parameters_ran = [0] * nParams

    rmse_save_posterior(parameters_ran, piece_directory, nsd-1)
    rmse_table = pd.read_csv(os.path.join(piece_directory,"rmse_priors_"+str(nsd-1)+".csv"))

    chicrit=stats.chi2.ppf(1-0.05, df=NUM_DAYS)
    print(chicrit)
    posterior_index=[]

# use big number for saving all/// relax conditions
    # for i in rmse_table['index']:
    #     if rmse_table['RMSE'][i] < chicrit and rmse_table['RMSE'][i]>0:
    #         posterior_index.append(i)
    # if len(posterior_index)==0: # make sure... posterior is not empty
    #     posterior_index=list(rmse_table.index.values)

    # Fix this: RMSE condition need to be relaxed to get enough parameter sets from prior
    # MIN_NUM_index=20
    # RMSE_THRESHOLD=1000

    MIN_NUM_index=20
    RMSE_THRESHOLD=1000000
    r_th=10
    while(len(list(set(posterior_index)))<MIN_NUM_index and r_th<RMSE_THRESHOLD): # sk & Shakir (10/20): relax conditions based on total # of simulations per 30 days
        for i in rmse_table['index']:
            if rmse_table['RMSE'][i]<r_th and rmse_table['RMSE'][i]>0:
                posterior_index.append(i)
        r_th+=1

    posterior_index=list(set(posterior_index))
    prior=pd.read_csv(os.path.join(piece_directory,"prior_parameters_sequential_"+str(nsd-1)+".csv"))
    posterior0=prior.iloc[posterior_index]
    infectedcases=[]
    for jk in posterior_index:
        output_directory =  piece_directory + str(int(jk)) + "/"
        START_DAY = (nsd-1) * NUM_DAYS
        filename="infections_from_new_strain"+str(START_DAY)+"_"+str(NUM_DAYS)+".csv"
        if(os.path.exists(filename)):
            inf=pd.read_csv(filename)
            infectedcases.append(inf['total_new_infections'][len(inf)-1])
    # sk 10/20
    minValues['INIT_FIXED_NUMBER_INFECTED']=min(infectedcases+[1]);
    maxValues['INIT_FIXED_NUMBER_INFECTED']=max(infectedcases+[1]);
#------------------Beta: home, school, work, community------------------------#
    minValues['BETA_SCALE']=min(posterior0['BETA_SCALE']);
    maxValues['BETA_SCALE']=max(posterior0['BETA_SCALE']);
    
    minValues['BETA_H']=min(posterior0['BETA_H'])#0.0001;
    maxValues['BETA_H']=max(posterior0['BETA_H']);#9
    
    minValues['BETA_PROJECT']=min(posterior0['BETA_PROJECT']);
    maxValues['BETA_PROJECT']=max(posterior0['BETA_PROJECT']);
    
    minValues['BETA_NBR_CELLS']=min(posterior0['BETA_NBR_CELLS']);
    maxValues['BETA_NBR_CELLS']=max(posterior0['BETA_NBR_CELLS']);
    
    minValues['BETA_CLASS']=min(posterior0['BETA_CLASS']);
    maxValues['BETA_CLASS']=max(posterior0['BETA_CLASS']);
    
    minValues['BETA_TRAVEL']=min(posterior0['BETA_TRAVEL']);
    maxValues['BETA_TRAVEL']=max(posterior0['BETA_TRAVEL']);
#-----------------------------------------------------------------------------#
    #priors(nParams)
    rangeParam={}
    for ld in minValues:
        rangeParam[ld]=maxValues[ld]-minValues[ld]
#-----------------------------------------------------------------------------#

    paraLHS=lhs(len(minValues), samples=nParams, criterion='corr');#to use the correlation feature, the dimension of the matrix should be more than 1.
    
#-----------------------------------------------------------------------------#
    i=0;
    VarParamsLHS=pd.DataFrame()
    for ld in minValues:
        VarParamsLHS[ld]=minValues[ld]+rangeParam[ld]*paraLHS[:,i]
        i=i+1

    VarParamsLHS['PROVIDE_INITIAL_SEED_GRAPH']=random.sample(range(1, nParams+1), nParams) # Change
    VarParamsLHS['PROVIDE_INITIAL_SEED']=random.sample(range(1, nParams+1), nParams) # Change
    # VarParamsLHS['PROVIDE_INITIAL_SEED_GRAPH']=1
    # VarParamsLHS['PROVIDE_INITIAL_SEED']=1
    return VarParamsLHS

def main():
    # Used for getting program completion time
    startTime = time.time()

    # Generate data
    generateCase()

    # Repeat for each piece
    for nsd in range(1,PIECE+1):
        #nsd=3 # for debugging... only (comment out for production)
        print("\nRunning MPI for piece ", nsd, " ////////////////////////////////////////////////////////////////////////////////////////////////////\n")

        # Create directories for each parameter
        for rank in range(NPARAMS):
            mpi_directory = MPI_DIR + "piece_"+str(nsd) + "/" + str(rank) + "/"
            if not os.path.exists(mpi_directory):
                os.system("mkdir -p " + mpi_directory)
                print("made output folders")

        # Calibration piece connection logic
        parameters = ''
        if nsd==1:
            parameters = priors(NPARAMS)
        else:
            # par1=priors(int(NPARAMS/2))
            # par2=priors_history(int(NPARAMS/2), nsd)
            # # print('par1',par1)
            # # print('par2',par2)
            # parameters = pd.concat([par1,par2],ignore_index=True)
            par2=priors_history(NPARAMS, nsd)
            # print('par1',par1)
            # print('par2',par2)
            parameters = par2
            
        #print('nsd',nsd,'parameters', parameters)
        #parameters["START_DAY"] = 1 + (nsd-1)*30
        parameters["START_DAY"] = (nsd-1)*30
        parameters["output_directory"] = MPI_DIR + "piece_"+str(nsd) + "/"
        parameters["input_directory"] = INPUT_DIR
        parameters.to_csv(os.path.join(MPI_DIR,"piece_"+str(nsd),"prior_parameters_sequential_"+str(nsd)+".csv"))

        # remain_params is the parameter remaining at the end
        # ex) remain_params=40 when NPARAMS is 400 and NPROCESSORS is 360
        # num_connect is the number of times we need to run the parallel simulations
        remain_params = 0 # sk: initializing
        if NPARAMS%NPROCESSORS != 0:
            remain_params = NPARAMS-NPROCESSORS
            num_connect = int(NPARAMS/NPROCESSORS) + 1
        else:
            remain_params = 0
            num_connect = int(NPARAMS/NPROCESSORS)

        # Run parallel simulations num_connect times
        for connect in range(num_connect):
            startTime_onebatch = time.time()
            # Running the last remaining simulations
            if remain_params != 0 and connect == num_connect-1:
                #command = "mpirun --host worker2:90 -np " + str(remain_params) + " \
                #command = "mpirun --hostfile /largedisk/mpi-test/host_file4 -npernode 1 -np " + str(remain_params) + " \
                command = "mpirun -map-by node --hostfile /largedisk/mpi-test/host_file_v2_12vm_1000core -np " + str(remain_params) + " \
                python run_parallel_simulations.py -c " + str(connect) + " -piece "+ str(nsd)+" -out_dir '"+str(MPI_DIR + "piece_"+str(nsd) + "/")+"'"
                
            # All other simulations
            else:
                #command = "mpirun --host worker2:90 -npernode " + str(NPERNODE) + " -np " + str(NPROCESSORS) + " \
                #command = "mpirun --hostfile /largedisk/mpi-test/host_file4 -npernode 1 -np " + str(NPROCESSORS) + " \
                command = "mpirun -map-by node --hostfile /largedisk/mpi-test/host_file_v2_12vm_1000core -np " + str(NPROCESSORS) + " \
                python run_parallel_simulations.py -c " + str(connect) + " -piece "+ str(nsd)+" -outdir '"+str(MPI_DIR + "piece_"+str(nsd) + "/")+"'"

            # Run command
            print(command)
            os.system(command)
            print("\n\n\n\n\nTime to complete a batch (", connect, " out of ", num_connect, "): ",time.time()-startTime_onebatch, " seconds")

    # Print program completion time
    print("\n\n\n\n\nProgram took total of ", time.time()-startTime, " seconds")

if __name__ == "__main__":
    main()

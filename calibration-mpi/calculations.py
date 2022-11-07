from header import *

def datamodel_Gap(model,data):
    rmse=0
    for i in range(len(model)):
        rmse+=(model[i]-data[i])*(model[i]-data[i])
    return np.sqrt(rmse)/len(model)

def rmse_save_posterior(parameters, piece_directory, nsd):
    nParams = len(parameters)
    RMSE_df=pd.DataFrame({"index":[0.0]*nParams,"RMSE":[0.0]*nParams,"p_value":[0.0]*nParams})
    RMSE_df=RMSE_df[:]*0.0
    csv_path=os.path.join(piece_directory,"casesfile_"+str(nsd)+".csv")
    data_temp=pd.read_csv(csv_path,parse_dates=['date']) ### fix....
    #dataCols = data_temp[['total Cases']]
    #data = dataCols['total Cases'].to_list()
    dataCols = data_temp[['total cases']]
    data = dataCols['total cases'].to_list()

    # print("i am in rmse\n")
    nSimPerDay = 4 #number of simulations per-day
#infections_from_new_strain0_150
    for i in range(len(parameters)):
        output_directory = piece_directory +str(int(i))+"/"
        START_DAY = (nsd-1) * NUM_DAYS # It is not passed from cpp_simulator(drive_simulator)'s command line. So be careful (SK 1101)
        filename="num_infected.csv"
        filename="infections_from_new_strain"+str(START_DAY)+"_"+str(NUM_DAYS)+".csv"
        #print("infections_from_new_strain: ", output_directory+filename)
        if(os.path.exists(output_directory+filename)):
            modeldata = pd.read_csv(output_directory+filename)
            modeldata = modeldata['total_new_infections']
            # modeldata = modeldata['num_infected']
            modeldata = modeldata.iloc[0::nSimPerDay].to_list()

            RMSE_df['RMSE'][i]=datamodel_Gap(modeldata,data)

            RMSE_df['p_value'][i] = 1000
            RMSE_df['index'][i] = i

    RMSE_df.to_csv(os.path.join(piece_directory, "rmse_priors_"+str(nsd)+".csv"),index=False)
    # print("i am done with rmse\n")

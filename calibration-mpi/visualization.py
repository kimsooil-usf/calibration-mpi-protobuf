from header import *
from calculations import rmse_save_posterior

def get_95CI_posterior(posterior_index, nsd):
    nParams = len(posterior_index)

    dataCols = DATA[['total Cases']]
    totpop=1450000
    data['variance']=np.sqrt((data['total Cases']/totpop)*(1-data['total Cases']/totpop))
    data['UL']=data['total Cases']+1.96*data['variance']
    data['LL']=data['total Cases']-1.96*data['variance']

    data = data['total Cases'].to_numpy()
    
    lendata = len(data)
    simdata = np.zeros(lendata)
    allmodels = pd.DataFrame()

    for i in posterior_index:
        output_directory =  OUTPUT_DIR +str(int(i))+"/"

        if(os.path.exists(os.path.join(output_directory,"num_infected.csv"))):
            modeldata=pd.read_csv(os.path.join(output_directory,"num_infected.csv"))

            modeldata=modeldata['num_infected']
            nSimPerDay = 4 #number of simulations per-day
            
            modeldata=modeldata.iloc[0::nSimPerDay]
            
            modeldata=modeldata.to_list()
            
            modeldata=np.array(modeldata)
            modeldata=modeldata[:lendata]
            allmodels[str(i)]=modeldata
            simdata=np.add(modeldata,simdata)
        
    mean=simdata/len(posterior_index)
    stdev=np.std(allmodels,axis=1)
    print(len(stdev),len(mean),len(allmodels))
    CI95=pd.DataFrame({"Mean":mean,"LL":mean-1.96*stdev/np.sqrt(len(mean)),"UL":mean+1.96*stdev/np.sqrt(len(mean))})
    return mean,stdev,allmodels,CI95


def makeDataFrame():
    CI95_full=pd.DataFrame()
    data_full=pd.DataFrame()

    for nsd in range(1,2):
        piece_directory =  OUTPUT_DIR + "piece_"+str(nsd) + "/"

        data = pd.read_csv(os.path.join(piece_directory, "casesfile_" + str(nsd) + ".csv"), parse_dates=['date'])
        data_full = pd.concat([data_full,data],ignore_index=True)
        parameters_ran = [0]*401#parameters.head(6853)
        rmse_save_posterior(parameters_ran, piece_directory, nsd)
        
        rmse_table = pd.read_csv(os.path.join(piece_directory, "rmse_priors_" + str(nsd) + ".csv"))
        posterior_index = []

        r_th = 100 * CHICRIT
        critlrmse = rmse_table[rmse_table['RMSE'] < r_th]

        if len(critlrmse > 0):
            for i in rmse_table['index']:
                if rmse_table['RMSE'][i]<r_th and rmse_table['RMSE'][i]>0:
                    posterior_index.append(i)

        posterior_index=list(set(posterior_index))
        if(len(posterior_index)>0):
            mn,sd,modl,CI95=get_95CI_posterior(posterior_index,piece_directory,nsd) 
            startdate1=data['date'][0]
            CI95['dates']=[startdate1 + dt.timedelta(days=int(x)) for x in range(1,len(CI95)+1)]
            CI95_full=pd.concat([CI95_full,CI95],ignore_index=True)

    return CI95_full, data_full

def plot(CI95_full, data_full):
    plt.figure(figsize=(10,6))

    plt.plot(CI95_full['dates'],CI95_full['Mean'])
    plt.plot(CI95_full['dates'],CI95_full['LL'])
    plt.plot(CI95_full['dates'],CI95_full['UL'])
    plt.plot(data_full['date'],data_full['total Cases'])
    plt.title("Infections")
    ax=plt.gca()
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
    plt.gcf().autofmt_xdate()

    plot_dir = PLOT_DIR
    plt.savefig(plot_dir + "plot" + '.png')

def main():
    CI95, data = makeDataFrame()
    plot(CI95, data)


if __name__ == "__main__":
    main()
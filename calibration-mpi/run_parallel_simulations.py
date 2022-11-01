from header import *
from mpi4py import MPI
import re

import multiprocessing as mp
from run_simulator import run_parameter_simulator

def getElements(data, n):
    for i in range(0, len(data), n):
        yield data[i:i+n]

def main():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-c')
    my_parser.add_argument('-piece')
    my_parser.add_argument('-outdir')
    args = my_parser.parse_args()

    connect = int(args.c)
    piece = int(args.piece)
    outdir = str(args.outdir)

    parameters = pd.read_csv(outdir+"prior_parameters_sequential_"+str(piece)+".csv")

    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()

    #comm.Barrier() ##############################################################

    print("Rank: ", rank, "\t-\t", time.ctime(time.time()))

    # Start parameter from where previous simulation ended
    rank += NPROCESSORS*connect
    #print('rank', rank, 'len(parameterArgument)', len(parameterArgument))
    command = run_parameter_simulator(parameters.values, rank)

    if piece==1:
        store_time_step = piece*NUM_DAYS*4 -1
        #command += " --STORE_STATE_TIME_STEP 119"
        command += " --STORE_STATE_TIME_STEP "+str(store_time_step)
        os.system(command)
        #os.system("gzip "+outdir+str(rank)+"/agentStore.pbstore") # jk 10/24
    else:
        load_time_step = (piece-1)*NUM_DAYS*4
        store_time_step = piece*NUM_DAYS*4 -1
        command += " --LOAD_STATE_TIME_STEP "+str(load_time_step)
        prev_out_dir=re.sub('piece_\d+', 'piece_'+str(piece-1), outdir)
        command += " --agent_load_file "+prev_out_dir+str(rank)+"/agentStore.pbstore"
        command += " --STORE_STATE_TIME_STEP "+str(store_time_step)
        #print("-----command", command)

        #os.system("gunzip "+prev_out_dir+str(rank)+"/agentStore.pbstore") # uncompress previous simulation's pbstore
        os.system(command)
        #comm.Barrier() ##############################################################
        #os.system("gzip "+prev_out_dir+str(rank)+"/agentStore.pbstore") # jk 10/24

    comm.Barrier() ##############################################################
    #MPI.Finalize() # sk: make sure... 

if __name__ == "__main__":
    main()
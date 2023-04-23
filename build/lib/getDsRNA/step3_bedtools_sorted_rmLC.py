import os,time,shutil
from multiprocessing import Pool
import subprocess
import re
import argparse


##----------------------------------------------------------------------------------------------------------
def sort_bed(ARG):
    file_path=ARG[0]
    sorted_path=ARG[1]
    #subprocess.Popen("sed -i '1d' "+file_path,shell=True).wait()
    subprocess.Popen(bedtools_path + " sort -i "+file_path+" > "+sorted_path,shell=True).wait()

def rmLC_bed(ARG):
    sorted_path=ARG[0]
    rmLC_path=ARG[1]
    subprocess.Popen(bedtools_path + " pairtobed -a "+sorted_path +" -b "+lc_path+" -f 0.5 -type neither > "+rmLC_path,shell=True).wait()
##----------------------------------------------------------------------------------------------------------

def main(args):    
    parser = argparse.ArgumentParser(description='python3 step3_bedtools_sorted_rmLC.py')
    parser.add_argument('-o', '--OUT_DIR',default = str(os.getcwd()),  help='path to Output Directory (default:Working Directory)')
    parser.add_argument('-B', '--Bedtools_path',default = str(os.getcwd())+"/bedtools",  help='path to BEDTOOLS (default:~/bedtools)')
    parser.add_argument('-lc', '--lc_path',default = str(os.getcwd())+"/Low_complexity.txt",  help='path to Low_complexity file (default:~/Low_complexity.txt)')
    parser.add_argument('-p', '--cpu',default=1,  help='CPU number (default=1)')

    try:
        args = parser.parse_args(args.split())

    except argparse.ArgumentError as e:
        print(e)
        exit(1)
    global CPU,base_path,bedtools_path,lc_path,rev_path,sort_path,rmlc_path
    CPU=int(int(args.cpu)/10)+1
    base_path=str(args.OUT_DIR)+"/"
    bedtools_path=str(args.Bedtools_path)
    lc_path=str(args.lc_path)

    rev_path=base_path+"ds_files_rev/"
    sort_path=base_path+"sorted_rev/"
    rmlc_path=base_path+"rmLC_rev/"

    if not os.path.isdir(sort_path):
        os.mkdir(sort_path)

    if not os.path.isdir(rmlc_path):
        os.mkdir(rmlc_path)

    ##----------------------------------------------------------------------------------------------------------
    print("Start step3 sort bedpe files and remove low complex sequences...")
    time_start = int(time.perf_counter())

    fn=os.listdir(rev_path)
    pool_size = CPU
    pool1 = Pool(pool_size)
    for f in fn:
        file_path=rev_path+f
        sorted_path=sort_path+f
        ARG=[file_path,sorted_path]
        pool1.apply_async(sort_bed, (ARG,)) 

    pool1.close() 
    pool1.join()  
    ##----------------------------------------------------------------------------------------------------------


    ##----------------------------------------------------------------------------------------------------------
    fn=os.listdir(sort_path)
    pool_size = CPU
    pool1 = Pool(pool_size)
    for f in fn:
        sorted_path=sort_path+f
        rmLC_path=rmlc_path+f
        ARG=[sorted_path,rmLC_path]
        pool1.apply_async(rmLC_bed, (ARG,)) 

    pool1.close() 
    pool1.join()  

    time_end = int(time.perf_counter())
    run_time = str(time_end-time_start)
    print("\t\tprocess running time is:" + run_time + " seconds" + "\n")        
##----------------------------------------------------------------------------------------------------------

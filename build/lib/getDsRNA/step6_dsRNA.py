import subprocess
import time
import os
from multiprocessing import Pool
import argparse

def dspair_add(ARG):
    ds_filter_hg19 = ARG[0]    
    ds_add_hg19 = ARG[1]
    add_num = int(ARG[2])
    
    with open(ds_filter_hg19) as fi, open(ds_add_hg19, "w") as fo:
        while True:
            lines = fi.readlines(100000)
            if not lines:
                break
            for line in lines:
                this_line = line.strip("\n")
                this_info = this_line.split("\t")
                this_chr = this_info[0]
                this_start1 = str(int(this_info[1]) - add_num)
                this_start2 = str(int(this_info[4]) - add_num)
                this_end1 = str(int(this_info[2]) + add_num)
                this_end2 = str(int(this_info[5]) + add_num)
                this_write = "\t".join([this_chr, this_start1, this_end1, this_chr, this_start2, this_end2])
                fo.write(this_write + "\n")

def main(args):
    parser = argparse.ArgumentParser(description='python3 step6_dsRNA.py')
    parser.add_argument('-o', '--OUT_DIR', default=os.getcwd(), help='path to Output Directory (default:Working Directory)')
    parser.add_argument('-p', '--cpu', default=1, help='CPU number (default=1)')
    args = parser.parse_args(args.split())
    global CPU,base_path,ds_dir,outdir
    CPU = int(args.cpu)
    base_path = str(args.OUT_DIR)+"/"
    ds_dir = base_path + "merge2_rev/"
    outdir = base_path + "dsRNA_file/"
    
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    print("Start step6 get candidate dsRNA...")
    time_start = int(time.perf_counter())
    
    fn = os.listdir(ds_dir)
    add_num = 300
    pool_size = CPU
    pool1 = Pool(pool_size)  
    for f in fn:
        ds_path = ds_dir + f
        file_out = outdir + f
        ARG = [ds_path, file_out, add_num]
        pool1.apply_async(dspair_add, (ARG,))

    pool1.close() 
    pool1.join()     

    subprocess.Popen("rm -r " + base_path + "*rev/", shell=True).wait()
    subprocess.Popen("rm -r " + base_path + "BLAT", shell=True).wait()
    subprocess.Popen("rm -r " + base_path + "query", shell=True).wait()
    subprocess.Popen("rm -r " + base_path + "reference", shell=True).wait()
    subprocess.Popen("rm -r " + base_path + "ds_files/", shell=True).wait()

    time_end = int(time.perf_counter())
    run_time = str(time_end - time_start)
    print("Processing running time is: " + run_time + " seconds" + "\n")


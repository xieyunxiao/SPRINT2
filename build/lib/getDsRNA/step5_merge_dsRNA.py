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
##----------------------------------------------------------------------------------------------------------

#ARG=[path_rev,path_merge,f]
def merge_ds(ARG):    
    path_rev=ARG[0]
    path_merge=ARG[1]
    f=ARG[2]

    this_sort=path_rev+f
    this_merge=path_merge+f
    this_tmp=path_merge+"tmp_"+f
    this_tmp_merge=path_merge+"merge_tmp_"+f
    
    #merge
    fi=open(this_sort)
    fo=open(this_merge,"w")
    fo_temp=open(this_tmp,"w")
    this_line=fi.readline()
    this_line=fi.readline()
    this_line=this_line.strip("\n")
    this_info=this_line.split("\t")

    if len(this_info)>5:
        this_chr=this_info[0]
        base_ds1_start=int(this_info[1])
        base_ds1_end=int(this_info[2])
        base_ds2_start=int(this_info[4])
        base_ds2_end=int(this_info[5])        
        #this_gene=this_info[6].split("-")[0]
        i=1        
        while 1 :
            lines=fi.readlines(100000000)
            if not lines:
                break
            for this_line in lines:
                if not this_line:
                    break            
                this_line=this_line.strip("\n")
                this_info=this_line.split("\t")
                
                i=i+1
                if len(this_info)<7:
                    continue
                    
                this_ds1_start=int(this_info[1])
                this_ds1_end=int(this_info[2])
                this_ds2_start=int(this_info[4])
                this_ds2_end=int(this_info[5])                
                
                if this_ds1_start >= base_ds1_start and this_ds1_start <= base_ds1_end:
                    base_ds1_end=max(this_ds1_end,base_ds1_end)
                    fo_temp.write("\t".join(this_info[3:6])+"\n")   

                else:
                    fo_temp.close()
                    subprocess.Popen("bedtools sort -i "+this_tmp+"  > "+this_tmp+"_sorted",shell=True).wait()
                    subprocess.Popen("bedtools merge -i "+this_tmp+"_sorted"+"  > "+this_tmp_merge,shell=True).wait()
                    fa=open(this_tmp_merge)
                    while 1 :
                        this_fa = fa.readline()
                        if not this_fa:
                            break
                        this_fa = this_fa.strip("\n")
                        #this_write=this_chr+"\t"+str(base_ds1_start)+"\t"+str(base_ds1_end)+"\t"+this_fa+"\t"+this_gene+"\n"
                        this_write=this_chr+"\t"+str(base_ds1_start)+"\t"+str(base_ds1_end)+"\t"+this_fa+"\n"
                        fo.write(this_write)
                    fa.close()
                    fo_temp=open(this_tmp,"w")
                    base_ds1_start=this_ds1_start
                    base_ds1_end=this_ds1_end
                    base_ds2_start=this_ds2_start
                    base_ds2_end=this_ds2_end       
    fi.close()
    fo.close()
##----------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------
def main(args):
    parser = argparse.ArgumentParser(description='python3 step5_merge_dsRNA.py')
    parser.add_argument('-o', '--OUT_DIR',default = str(os.getcwd()),  help='path to Output Directory (default:Working Directory)')
    parser.add_argument('-B', '--Bedtools_path',default = str(os.getcwd())+"/bedtools",  help='path to BEDTOOLS (default:~/bedtools)')
    parser.add_argument('-p', '--cpu',default=1,  help='CPU number (default=1)')

    try:
        args = parser.parse_args(args.split())

    except argparse.ArgumentError as e:
        print(e)
        exit(1)
    global CPU, base_path, bedtools_path, com_path, sort_path
    CPU=int(int(args.cpu)/5)+1
    base_path=str(args.OUT_DIR)+"/"
    bedtools_path=str(args.Bedtools_path)
    com_path=base_path+"combine_rev/"
    sort_path=base_path+"sort2_rev/"

    if not os.path.isdir(sort_path):
        os.mkdir(sort_path)

    ##----------------------------------------------------------------------------------------------------------
    print("Start step5 merge_dsRNA...")   
    fn=os.listdir(com_path)
    time_start = int(time.perf_counter())
    pool_size = CPU
    pool1 = Pool(pool_size)
    for f in fn:
        file_path=com_path+f
        sorted_path=sort_path+f
        ARG=[file_path,sorted_path]
        pool1.apply_async(sort_bed, (ARG,)) 

    pool1.close() 
    pool1.join()  
    ##----------------------------------------------------------------------------------------------------------


    ##----------------------------------------------------------------------------------------------------------
    CPU=int(args.cpu)
    path_merge=base_path+"merge2_rev/"
    fn=os.listdir(sort_path)
    if not os.path.isdir(path_merge):
        os.mkdir(path_merge)


    pool_size = CPU
    pool1 = Pool(pool_size)

    for f in fn:
        ARG=[sort_path,path_merge,f]
        pool1.apply_async(merge_ds, (ARG,)) 

    pool1.close() 
    pool1.join()  
    subprocess.Popen("rm "+path_merge+"*tmp*",shell=True).wait()
    time_end = int(time.perf_counter())
    run_time = str(time_end-time_start)
    print("\t\tprocess running time is:" + run_time + " seconds" + "\n")     
    ##----------------------------------------------------------------------------------------------------------














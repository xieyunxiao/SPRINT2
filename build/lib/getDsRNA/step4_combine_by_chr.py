import os,time
from multiprocessing import Pool
import subprocess
import argparse


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

    if len(this_info)>7:
        this_chr=this_info[0]
        base_ds1_start=int(this_info[1])
        base_ds1_end=int(this_info[2])
        base_ds2_start=int(this_info[4])
        base_ds2_end=int(this_info[5])        
        this_gene=this_info[6].split("-")[0]
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
                    subprocess.Popen(bedtools_path+" sort -i "+this_tmp+"  > "+this_tmp+"_sorted",shell=True).wait()
                    subprocess.Popen(bedtools_path+" merge -i "+this_tmp+"_sorted"+"  > "+this_tmp_merge,shell=True).wait()
                    fa=open(this_tmp_merge)
                    while 1 :
                        this_fa = fa.readline()
                        if not this_fa:
                            break
                        this_fa = this_fa.strip("\n")
                        this_write=this_chr+"\t"+str(base_ds1_start)+"\t"+str(base_ds1_end)+"\t"+this_fa+"\t"+this_gene+"\n"
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
#ARG=[rev_file,CHR,this_chr_path]
def combine_by_chr(ARG):
    rev_file=ARG[0]
    this_chr=ARG[1]
    this_chr_path=ARG[2]
    fi=open(rev_file)
    fo=open(this_chr_path,'a+')
    while 1:
        lines=fi.readlines(10000)
        if not lines:
            break
        for line in lines:
            if not line:
                break
            this_line=line.strip("\n")
            this_chr0=this_line.split("\t")[0]        
            if this_chr0==this_chr:
                fo.write(this_line+"\n")
    fi.close()
    fo.close()
##----------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------
def main(args):    
    print("Start step4 combine_by_chr...")     
    parser = argparse.ArgumentParser(description='python3 step4_combine_by_chr.py')
    parser.add_argument('-o', '--OUT_DIR',default = str(os.getcwd()),  help='path to Output Directory (default:Working Directory)')
    parser.add_argument('-t','--Transcript_file',default = str(os.getcwd())+"/transcript.bed", help='path to Transctript Annotation File (BED format) (default:~/transcript.bed)')
    parser.add_argument('-B', '--Bedtools_path',default = str(os.getcwd())+"/bedtools",  help='path to BEDTOOLS (default:~/bedtools)')
    parser.add_argument('-p', '--cpu',default=1,  help='CPU number (default=1)')

    try:
        args = parser.parse_args(args.split())
    except argparse.ArgumentError as e:
        print(e)
        exit(1)
    global CPU,base_path,trans_path,bedtools_path,path_rev,path_merge
    CPU=int(args.cpu)
    base_path=str(args.OUT_DIR)+"/"
    trans_path=str(args.Transcript_file)
    bedtools_path=str(args.Bedtools_path)

    path_rev=base_path+"rmLC_rev/"
    path_merge=base_path+"merge_rev/"

    if not os.path.isdir(path_merge):
        os.mkdir(path_merge)

    fn=os.listdir(path_rev)
    time_start = int(time.perf_counter())
    pool_size = CPU
    pool1 = Pool(pool_size)

    for f in fn:
        ARG=[path_rev,path_merge,f]
        pool1.apply_async(merge_ds, (ARG,)) 

    pool1.close() 
    pool1.join()  
    subprocess.Popen("rm "+path_merge+"*tmp*",shell=True).wait()
    print("\t\t start combine...")     
    ##----------------------------------------------------------------------------------------------------------


    ##----------------------------------------------------------------------------------------------------------
    rev_path=base_path+"merge_rev/"
    if not os.path.isdir(rev_path):
        os.mkdir(rev_path)

    com_path=base_path+"combine_rev/"
    if not os.path.isdir(com_path):
        os.mkdir(com_path)
    pool_size = CPU
    pool1 = Pool(pool_size)
    fi=open(trans_path)
    while 1:
        lines=fi.readlines(10000)
        if not lines:
            break
        for line in lines:
            if not line:
                break
            this_line=line.strip("\n")
            this_info=this_line.split("\t")
            this_chr0=this_info[0]
            this_name=this_info[5]+"-"+str(this_info[1])+"-"+str(this_info[2])
            rev_file=rev_path+"ds_"+this_name+".bedpe"
            #ds_chr1_id_2.bedpairs
            ch=this_chr0
            this_chr_path=com_path+ch+"_rev.bedpe"
            ARG=[rev_file,ch,this_chr_path]
            pool1.apply_async(combine_by_chr, (ARG,))

    fi.close()
    pool1.close() 
    pool1.join()  

    time_end = int(time.perf_counter())
    run_time = str(time_end-time_start)
    print("\t\tprocess running time is:" + run_time + " seconds" + "\n")     

##----------------------------------------------------------------------------------------------------------










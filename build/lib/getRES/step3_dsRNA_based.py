import os,time,shutil
from multiprocessing import Pool
import subprocess
import re
import argparse



##----------------------------------------------------------------------------------------------------------

#ARG=[bedpe_file,snv_file,outpath]
def get_dsRES(ARG):
    ds_filter=ARG[0]
    allsnv_path=ARG[1]
    outpath=ARG[2]
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    allsnv_sorted_path=outpath+"fa_genome_all_snv_sorted.bed"
    subprocess.Popen(bedtools_path + " sort -i "+allsnv_path+" > "+allsnv_sorted_path,shell=True).wait()
    
    ##----------------------------------------------------------------------------------------------------------
    ##fb_snv_duplets.bed
    ##fc_ds_snv_filter1.bedpe
    # 找duplets
    allsnv_sorted_path=outpath+"fa_genome_all_snv_sorted.bed"
    allsnv_duplets_path=outpath+"fb_snv_duplets.bed"

    fi=open(allsnv_sorted_path)
    snv_duplets={}
    snv_duplets_num={}

    this_chr1="chr0"
    this_pos1=int(0)
    this_type1="AA"
    while 1 :
        lines=fi.readlines(10000)
        if not lines:
            break
        for line in lines:
            this_snv=line.strip("\n")
            this_chr=this_snv.split("\t")[0]
            this_pos=int(this_snv.split("\t")[2])
            this_type=this_snv.split("\t")[3]

            if this_chr==this_chr1 and this_pos-this_pos1<=200 and this_type==this_type1:
                snv_duplets[this_snv1].append(this_snv)
                snv_duplets_num[this_snv1]=snv_duplets_num[this_snv1]+1
            else:
                this_chr1=this_chr
                this_pos1=this_pos
                this_type1=this_type
                this_snv1=this_snv
                snv_duplets[this_snv1]=[this_snv1]
                snv_duplets_num[this_snv1]=1
    fi.close()

    all_snv=list(snv_duplets.keys())
    fo=open(allsnv_duplets_path,"w")
    i=0
    for snv in all_snv:
        this_num=int(snv_duplets_num[snv])
        if this_num>1:
            i=i+1        
            this_cluster=snv_duplets[snv]
            for this_snv in this_cluster:
                this_write=this_snv+"\tC"+str(i)+"E\t"+str(this_num)+"\n"
                fo.write(this_write)    
    fo.close()

    del snv_duplets
    del snv_duplets_num
    
    ##----------------------------------------------------------------------------------------------------------
    #同时落在两段dsRNA上的SNV
    ds_snv_path=outpath+"fc_ds_snv_filter1.bedpe"
    subprocess.Popen(bedtools_path+" pairtobed -a "+ds_filter+" -b "+allsnv_duplets_path+"  > "+ds_snv_path,shell=True).wait()
    
    ##----------------------------------------------------------------------------------------------------------
    ##fd_ds_snv_pair_double_strands.bedpe
    ##fd_ds_snv_pair_single_strand.bedpe
    ##fe_snv_ds_filter_double_strands.bed
    ##fe_snv_ds_filter_double_strands.bed
    #dsRNA两条链上同时包含SNV的保留下来
    ds_snv_path=outpath+"fc_ds_snv_filter1.bedpe"
    ds_snv_pair=outpath+"fd_ds_snv_pair_double_strands.bedpe"
    ds_snv_pair2=outpath+"fd_ds_snv_pair_single_strand.bedpe"
    ds_snv_filter=outpath+"fe_snv_ds_filter_double_strands.bed"
    ds_snv_filter2=outpath+"fe_snv_ds_filter_single_strand.bed"
    
    fi=open(ds_snv_path)
    fo1=open(ds_snv_pair,'w')
    fo11=open(ds_snv_pair2,'w')
    fo2=open(ds_snv_filter,'w')
    fo3=open(ds_snv_filter2,'w')

    
    a=["A","T","C","G"]
    b=["T","A","G","C"]
    all_types={}
    for i in range(4):
        all_types[a[i]]=b[i]
    
    this_ds0="chr0\t0\t0\tchr0\t0\t0"
    this_ds_snv_num=[0,0]
    this_ds_snv=[[],[]]
    while 1 :
        lines=fi.readlines(100000)
        if not lines:
            break
        for line in lines:
            this_pair=line.strip("\n")
            this_pair=this_pair.split("\t")
            this_ds="\t".join(this_pair[0:6])

            this_snv="\t".join(this_pair[6:])
            ds1_start=int(this_pair[1])
            ds1_end=int(this_pair[2])
            ds2_start=int(this_pair[4])
            ds2_end=int(this_pair[5])    
            snv_pos=int(this_pair[8])
            
            if this_ds == this_ds0:
                if snv_pos>=ds1_start and snv_pos<= ds1_end:
                    this_ds_snv_num[0] = this_ds_snv_num[0]+1
                    this_ds_snv[0].append(this_snv)
                else:
                    this_ds_snv_num[1] = this_ds_snv_num[1]+1
                    this_ds_snv[1].append(this_snv)                   
                    
            else:
                if this_ds_snv_num[0]>0 and this_ds_snv_num[1]>0:
                    this_ds_seq1_snv=list(set(this_ds_snv[0]))
                    this_ds_seq2_snv=list(set(this_ds_snv[1]))
                    for snv1 in this_ds_seq1_snv:
                        snv1_info=snv1.split("\t")
                        snv1_type=snv1_info[3]                               
                        for snv2 in this_ds_seq2_snv:
                            snv2_info=snv2.split("\t")
                            snv2_type=snv2_info[3]                    
                            if snv1_type[0] in [snv2_type[0],all_types[snv2_type[0]]] and  snv1_type[1] in [snv2_type[1],all_types[snv2_type[1]]]:
                                this_cl=int(snv1_info[6])+int(snv2_info[6])
                                this_write1=this_ds+"\t"+snv1+"\t"+snv2+"\t"+str(this_cl)+"\n"
                                fo1.write(this_write1)
                                this_write="\t".join(snv1_info[0:6])+"\t"+str(this_cl)+"\n"                        
                                fo2.write(this_write)  
                                this_write="\t".join(snv2_info[0:6])+"\t"+str(this_cl)+"\n"       
                                fo2.write(this_write)                           

                            else:
                                this_write="\t".join(snv1_info[0:6])+"\t"+str(snv1_info[6])+"\n"          
                                fo3.write(this_write)
                                this_write="\t".join(snv2_info[0:6])+"\t"+str(snv2_info[6])+"\n"
                                fo3.write(this_write)
                                this_write=this_ds+"\t"+"\t".join(snv1_info[0:6])+"\t"+str(snv1_info[6])+"\n"
                                fo11.write(this_write)
                                this_write=this_ds+"\t"+"\t".join(snv2_info[0:6])+"\t"+str(snv2_info[6])+"\n"
                                fo11.write(this_write)
                                                       
                        
                else:
                    this_ds_seq1_snv=list(set(this_ds_snv[0]))
                    this_ds_seq2_snv=list(set(this_ds_snv[1]))
                    ds1_snv_num=len(this_ds_seq1_snv)
                    ds2_snv_num=len(this_ds_seq2_snv)
                    if ds1_snv_num !=0:
                        for snv1 in this_ds_seq1_snv:
                            snv1_info=snv1.split("\t")
                            this_write="\t".join(snv1_info[0:6])+"\t"+str(snv1_info[6])+"\n"
                            fo3.write(this_write)
                            this_write=this_ds+"\t"+"\t".join(snv1_info[0:6])+"\t"+str(snv1_info[6])+"\n"
                            fo11.write(this_write)
                    if ds2_snv_num !=0:
                        for snv2 in this_ds_seq2_snv:
                            snv2_info=snv2.split("\t")
                            this_write="\t".join(snv2_info[0:6])+"\t"+str(snv2_info[6])+"\n"
                            fo3.write(this_write)
                            this_write=this_ds+"\t"+"\t".join(snv2_info[0:6])+"\t"+str(snv2_info[6])+"\n"
                            fo11.write(this_write)
          
                                           
                this_ds0=this_ds
                this_ds_snv_num=[0,0]
                this_ds_snv=[[],[]]
                if snv_pos>=ds1_start and snv_pos<= ds1_end:
                    this_ds_snv_num[0] = this_ds_snv_num[0]+1
                    this_ds_snv[0].append(this_snv)
                else:
                    this_ds_snv_num[1] = this_ds_snv_num[1]+1
                    this_ds_snv[1].append(this_snv)
    fi.close()   
    fo1.close()
    fo2.close()
    fo3.close()
    fo11.close()

    subprocess.Popen("sort "+ds_snv_filter+" | uniq > "+ds_snv_filter+".uniq",shell=True).wait()
    subprocess.Popen("sort "+ds_snv_filter2+" | uniq > "+ds_snv_filter2+".uniq",shell=True).wait()
    subprocess.Popen("rm "+ds_snv_filter,shell=True).wait()
    subprocess.Popen("rm "+ds_snv_filter2,shell=True).wait()

    ##----------------------------------------------------------------------------------------------------------
    ##ff_snv_filter_ds
    ##ff_snv_filter_ss        
    ##filter RES double strands
    snv_ds=outpath+"fe_snv_ds_filter_double_strands.bed.uniq"
    snv_filter=outpath+"ff_res_filter_ds.bed"
    fi1=open(snv_ds)
    fsnv=open(snv_filter,'w')
    while 1:
        lines=fi1.readlines(10000)
        if not lines:
            break
        for line in lines:
            if not line:
                break
            this_line=line.strip("\n")          
            this_info=this_line.split("\t")
            this_ad=int(this_info[4])
            this_cl=int(this_info[6])
            if this_ad >=1 and this_cl>=4:                    
                fsnv.write("\t".join(this_info[0:5])+"\n")
    fi1.close()
    fsnv.close()

    subprocess.Popen(bedtools_path + " sort -i "+snv_filter+" | uniq "+" > "+snv_filter+".uniq",shell=True).wait()
    subprocess.Popen("rm "+snv_filter,shell=True).wait()

    ##----------------------------------------------------------------------------------------------------------
    ##filter RES single strand
    snv_ds=outpath+"fe_snv_ds_filter_single_strand.bed.uniq"
    snv_filter=outpath+"ff_res_filter_ss.bed"
    fi1=open(snv_ds)
    fsnv=open(snv_filter,'w')
    while 1:
        lines=fi1.readlines(10000)
        if not lines:
            break
        for line in lines:
            if not line:
                break
            this_line=line.strip("\n")           
            this_info=this_line.split("\t")
            this_ad=int(this_info[4])
            this_cl=int(this_info[6])
            if this_ad >=1 and this_cl>=4:                    
                fsnv.write("\t".join(this_info[0:5])+"\n")
    fi1.close()
    fsnv.close()

    subprocess.Popen(bedtools_path + " sort -i "+snv_filter+" | uniq "+" > "+snv_filter+".uniq",shell=True).wait()
    subprocess.Popen("rm "+snv_filter,shell=True).wait()    
    subprocess.Popen("cat "+outpath+"ff_res_filter_ds.bed.uniq "+outpath+"ff_res_filter_ss.bed.uniq | "+bedtools_path+" sort -i | uniq > "+outpath+"ff_res_filter_all.res",shell=True).wait()

    ##----------------------------------------------------------------------------------------------------------
    ##fg_RES_summery
    #统计
    snv_path=outpath+"ff_res_filter_all.res"
    all_type=[]
    all_type_num={}

    fi=open(snv_path)
    while 1 :
        lines=fi.readlines(10000)
        if not lines:
            break
        for line in lines:
            if not line:
                break
            this_line=line.strip("\n")
            this_snv=this_line.split("\t")
            this_type=this_snv[3]
            if this_type not in all_type:
                all_type.append(this_type)
                all_type_num[this_type]=1
            else:
                all_type_num[this_type]=int(all_type_num[this_type])+1
    fi.close()

    all_snv_num=int(sum(list(all_type_num.values())))

    all_type=list(all_type_num.keys())
    normal_type=["AG","TC","CT","GA"]
    for tp in normal_type:
        if tp not in all_type:
            all_type_num[tp]=0
    
    all_type_num["AG+TC"]=int(all_type_num["AG"])+int(all_type_num["TC"])
    all_type_num["CT+GA"]=int(all_type_num["CT"])+int(all_type_num["GA"])
    all_type_num["AG+TC+CT+GA"]=int(all_type_num["AG+TC"])+int(all_type_num["CT+GA"])
    all_type_num["ALL_RES"]=all_snv_num

    all_type=list(all_type_num.keys())
    all_num=list(all_type_num.values())

    if all_snv_num == 0:
        all_ratio=[round(n/(all_snv_num+1),3) for n in all_num]
    else:
        all_ratio=[round(n/all_snv_num,3) for n in all_num]

    summary_path=outpath+"fg_RES_summery.txt"
    fo=open(summary_path, 'w')
    this_line="Type\tNumber\tRatio\n"
    fo.write(this_line)
    for i in range(len(all_num)):
        this_line=str(all_type[i])+"\t"+str(all_num[i])+"\t"+str(all_ratio[i])+"\n"
        fo.write(this_line)
    fo.close()
##----------------------------------------------------------------------------------------------------------

#ARG=[res_file,stat_file]
def stat_res(ARG):
    snv_path=ARG[0]
    summary_path=ARG[1]
    all_type=[]
    all_type_num={}

    fi=open(snv_path)
    while 1 :
        lines=fi.readlines(10000)
        if not lines:
            break
        for line in lines:
            if not line:
                break
            this_line=line.strip("\n")
            this_snv=this_line.split("\t")
            this_type=this_snv[3]
            if this_type not in all_type:
                all_type.append(this_type)
                all_type_num[this_type]=1
            else:
                all_type_num[this_type]=int(all_type_num[this_type])+1
    fi.close()

    all_snv_num=int(sum(list(all_type_num.values())))

    all_type=list(all_type_num.keys())
    normal_type=["AG","TC","CT","GA"]
    for tp in normal_type:
        if tp not in all_type:
            all_type_num[tp]=0
    
    all_type_num["AG+TC"]=int(all_type_num["AG"])+int(all_type_num["TC"])
    all_type_num["CT+GA"]=int(all_type_num["CT"])+int(all_type_num["GA"])
    all_type_num["AG+TC+CT+GA"]=int(all_type_num["AG+TC"])+int(all_type_num["CT+GA"])
    all_type_num["ALL_RES"]=all_snv_num

    all_type=list(all_type_num.keys())
    all_num=list(all_type_num.values())

    if all_snv_num == 0:
        all_ratio=[round(n/(all_snv_num+1),3) for n in all_num]
    else:
        all_ratio=[round(n/all_snv_num,3) for n in all_num]

    fo=open(summary_path, 'w')
    this_line="Type\tNumber\tRatio\n"
    fo.write(this_line)
    for i in range(len(all_num)):
        this_line=str(all_type[i])+"\t"+str(all_num[i])+"\t"+str(all_ratio[i])+"\n"
        fo.write(this_line)
    fo.close()
##----------------------------------------------------------------------------------------------------------

def getBedDP(bed_in_path=0, bam_in_path=0, bed_out_path=0, bedtools_path=0):
    this_step =subprocess.Popen(' '.join([bedtools_path,'multicov','-bams',bam_in_path,'-bed',bed_in_path,'>',bed_out_path]),shell=True)
    this_step.wait()

##----------------------------------------------------------------------------------------------------------
def main(args):

    parser = argparse.ArgumentParser(description="python3 step3_dsRNA_based.py")
    parser.add_argument("-o", "--output",default = str(os.getcwd()),  help="path to Output Directory (defalt:Working Directory)")
    parser.add_argument("-s", "--snv",default = str(os.getcwd()),  help="path to SNV Directory (defalt:Working Directory)")
    parser.add_argument("-ds", "--dsrna",default = str(os.getcwd())+"/dsRNA_file/",  help="path to candidate dsRNA Directory (defalt:~/dsRNA_file/)")
    parser.add_argument("-b", "--bedtools",default = str(os.getcwd())+"/bedtools",  help="path to bedtools (defalt:~/bedtools)")
    parser.add_argument("-p", "--cpu",default = 1,  help="CPU number (defalt=1)")

    try:
        args = parser.parse_args(args.split())
    except argparse.ArgumentError as e:
        print(e)
        exit(1)

    global allsnv_path,outpath1,bedtools_path,ds_path,snv_path,outpath0
    outpath0=str(args.output)+"/"
    snv_path=str(args.snv)+"/"
    ds_path=str(args.dsrna)+"/"
    bedtools_path=str(args.bedtools)
    CPU=int(args.cpu)

    #print("Step3 Identifying RESs based dsRNA... ") 
    allsnv_path=snv_path+"/regular.snv"
    outpath1=outpath0+"ds_regular/"
    if not os.path.isdir(outpath0):
        os.mkdir(outpath0)
    if not os.path.isdir(outpath1):
        os.mkdir(outpath1)
    ##----------------------------------------------------------------------------------------------------------

    print("Step3 Identifying RESs based double stranded RNAs...")

    time_start = int(time.perf_counter())
    fn=os.listdir(ds_path)

    pool_size = CPU
    pool1 = Pool(pool_size)  
    for f in fn:    
        ch=f.split("_")[0]
        outpath2=outpath1+ch+"/"
        if not os.path.isdir(outpath2):
            os.mkdir(outpath2)
        ds_filter=ds_path+f
        ARG=[ds_filter,allsnv_path,outpath2]
        pool1.apply_async(get_dsRES, (ARG,))

    pool1.close() 
    pool1.join()     
    ##----------------------------------------------------------------------------------------------------------


    ##----------------------------------------------------------------------------------------------------------

    subprocess.Popen("cat "+outpath1+"/chr*/ff_res_filter_all.res | "+bedtools_path+" sort -i | uniq > "+outpath0+"ds_regular.res",shell=True).wait()
    ARG=[outpath0+"ds_regular.res",outpath0+"ds_regular_res.stat"]
    stat_res(ARG)

    time_end = int(time.perf_counter())
    run_time = str(time_end-time_start)
    print("dsRNA-based RES calling running time is:" + run_time + " seconds" + "\n")     

##----------------------------------------------------------------------------------------------------------




##cisNAT

import os,time,shutil
from multiprocessing import Pool
import subprocess

bedtools_path="~/bin/bedtools"
BASE_path="~/cisNAT/"
ds1_path="~/cisNAT/"

cis_NAT=BASE_path+"gene_ds_cisNAT.bed"
#fa_snv="fa_genome_all_snv_sorted.bed"
fb_duplets="2_RES_calling/ds_regular/chr21/fb_snv_duplets.bed"


def cisNAT_RES(ARG):

    ds_path=ds1_path
    base_path=BASE_path
    ##snv_filter_by cisNAT
    file_a=ds_path+fb_duplets
    file_b=cis_NAT
    file_c=base_path+"fc_snv_ds_filter.bedpe"
    subprocess.Popen(bedtools_path+" intersect -a "+file_a+" -b "+file_b+" -wao > "+file_c,shell=True).wait()

    ##filter RES
    snv_ds=base_path+"fc_snv_ds_filter.bedpe"

    snv_filter1=base_path+"fd_snv_filter.bed"
    
    fi1=open(snv_ds)

    fsnv=open(snv_filter1,'w')
    
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
            this_cl=int(this_info[7])
            if int(this_info[11])==1:                 
                if this_ad >=1 and this_cl>=3:                    
                    fsnv.write("\t".join(this_info[0:4])+"\tcis_NAT\t"+"\n")
                    
                if this_ad >=2 and this_cl>=2:
                    fsnv.write("\t".join(this_info[0:4])+"\tcis_NAT\t"+"\n")

    fi1.close()
    fsnv.close()
    
    snv_filter1=base_path+"fd_snv_filter.bed"
    snv_filter3=base_path+"fd_snv_filter_uniq.bed"
    subprocess.Popen(bedtools_path + " sort -i "+snv_filter1+" | uniq "+" > "+snv_filter3,shell=True).wait()
    
    #################################################
    ##stat
    snv_path=base_path+"fd_snv_filter_uniq.bed"
    all_type=["AG","TC","CT","GA"]
    all_type_num={"AG":0,"TC":0,"CT":0,"GA":0}

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

    all_type_num["AG+TC"]=int(all_type_num["AG"])+int(all_type_num["TC"])
    all_type_num["CT+GA"]=int(all_type_num["CT"])+int(all_type_num["GA"])
    all_type_num["AG+TC+CT+GA"]=int(all_type_num["AG+TC"])+int(all_type_num["CT+GA"])
    all_type_num["ALL"]=all_snv_num
    
    all_type=list(all_type_num.keys())
    all_num=list(all_type_num.values())

    all_ratio=[n/all_snv_num for n in all_num]

    summary_path=base_path+"fd_summary.txt"
    fo=open(summary_path, 'a+')
    this_line="Type\tNumber\tRatio\n"
    fo.write(this_line)
    for i in range(len(all_num)):
        this_line=str(all_type[i])+"\t"+str(all_num[i])+"\t"+str(all_ratio[i])+"\n"
        fo.write(this_line)

    fo.write("\n")
    fo.write("\n")
    fo.write("\n")
    fo.close()
    
    #############################################################################################################


cisNAT_RES()


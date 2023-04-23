##calling RES
##annotation-based
import subprocess,os,sys,time
import random
import argparse

def annotate(bed_in_dir=0,bed_anno_dir=0,bed_out_dir=0):
    if bed_in_dir ==0 or bed_anno_dir ==0 or bed_out_dir ==0:
        print("Please check the input directory! annotate(bed_in_dir,bed_anno_dir,bed_out_dir)")
        return("Please check the input directory! annotate(bed_in_dir,bed_anno_dir,bed_out_dir)")
    else:
        fi=open(bed_in_dir) #must by sorted
        fa=open(bed_anno_dir)
        fo=open(bed_out_dir,"w")
        anno={}
        for line in fa:
            seq=line[0:-1].split("\t")
            try:
                anno[seq[0]].append([int(seq[1])+1,int(seq[2]),seq[3],seq[4],seq[5]])    
            except Exception as e:
                anno[seq[0]]=[ [int(seq[1])+1,int(seq[2]),seq[3],seq[4],seq[5]]  ]
        for a in anno:
            anno[a].sort()
        top=0
        point=0
        lastchr=""
        for line in fi:
        
            seq = line.rstrip().split("\t")
            if 1==1:
                output=line.replace("\n","")
                if seq[0]==lastchr:
                    1+1
                else:
                    top=0
                    lastchr=seq[0]
                if seq[0] not in anno:
                    anno[seq[0]]=[ [0,0,0,0,0]  ]
                if top<len(anno[seq[0]]):            
                    while anno[seq[0]][top][1] < int(seq[2]) :
        
                        top=top+1
                        if top >= len(anno[seq[0]]):
                            break
                    point=top
                    if point < len(anno[seq[0]]):
                        while anno[seq[0]][point][0] <= int(seq[2])     and     anno[seq[0]][point][1] >= int(seq[2]) :
                                if anno[seq[0]][point][2]+"\t"+anno[seq[0]][point][3]+"\t"+anno[seq[0]][point][4] not in output:
                                    output=output+"\t"+anno[seq[0]][point][2]+"\t"+anno[seq[0]][point][3]+"\t"+anno[seq[0]][point][4]
                                point=point+1
                                if point >= len(anno[seq[0]]):
                                    break
        
                fo.write(output+"\n")
        fi.close()
        fa.close()
        fo.close()

##----------------------------------------------------------------------------------------------------------

def seperate(bed_in_dir,flag_out_dir,rp_out_dir,nonrp_out_dir,flag):
    fi=open(bed_in_dir)
    fo_flag=open(flag_out_dir,"w")
    fo_rp=open(rp_out_dir,"w")
    fo_nonrp=open(nonrp_out_dir,"w")
    for line in fi:
        if "Simple_repeat" in line or "Low_complexity" in line:
            next
        elif flag in line:
            fo_flag.write(line)
        elif "Repeat_region" in line:
            fo_rp.write(line)
        else:
            fo_nonrp.write(line)
    fi.close()
    fo_flag.close()
    fo_rp.close()
    fo_nonrp.close()

##----------------------------------------------------------------------------------------------------------

def get_snv_with_ad(snv_in_dir=0,snv_out_dir=0,flag=0):
    fi=open(snv_in_dir)
    fo=open(snv_out_dir,"w")
    
    for line in fi:
        seq=line.split("\t")
        try:
            if int(seq[4])>=int(flag):
                fo.write(line)
        except Exception as e:
            print(seq[0]+"\t"+seq[2]+"\twithout\tAD flag\n")        
    fi.close()
    fo.close()

##----------------------------------------------------------------------------------------------------------

def snv_cluster(bed_in_dir=0,bed_out_dir=0,cluster_distance=-1,cluster_size=-1):
    fi=open(bed_in_dir)
    fo=open(bed_out_dir,"w")
    tmp="chr0:0:AA"
    limitdistance=int(cluster_distance)
    limitnum=int(cluster_size)
    lst=[]
    for line in fi:
            seq=line.split("\t")
            tmpseq=tmp.split(":")
            if seq[0]==tmpseq[0] and int(seq[2])-int(tmpseq[1])<=limitdistance and seq[3]==tmpseq[2]:
                lst.append(line)
            else:
                if len(lst)>=limitnum:
                        begin=float(lst[0].split("\t")[1])
                        end=float(lst[-1].split("\t")[2])
                        density=len(lst)/(end-begin)
                        for one in lst:
                            fo.write(one[0:-1]+"\t"+str(len(lst))+"\t"+str(density)+"\n")
                lst=[]
                lst.append(line)
            tmp=seq[0]+":"+seq[2]+":"+seq[3]
    if len(lst)>=limitnum:
            begin=float(lst[0].split("\t")[1])
            end=float(lst[-1].split("\t")[2])
            density=len(lst)/(end-begin)
            for one in lst:
                fo.write(one[0:-1]+"\t"+str(len(lst))+"\t"+str(density)+"\n")
    fi.close()
    fo.close()

##----------------------------------------------------------------------------------------------------------

def bed_or(bed_in_dir1,bed_in_dir2,bed_out_dir):
    f1=open(bed_in_dir1)
    f2=open(bed_in_dir2)
    fo=open(bed_out_dir,"w")
    whole=[]
    for line in f1:
        seq=line.replace("\n","").split("\t")
        whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5]])
    for line in f2:
        seq=line.replace("\n","").split("\t")
        whole.append([seq[0],int(seq[2]),seq[3],seq[4],seq[5]])
    whole.sort()
    old=set()
    for one in whole:
        if one[0]+":"+str(one[1]) not in old:
            fo.write(one[0]+"\t"+str(one[1]-1)+"\t"+str(one[1])+"\t"+one[2]+"\t"+one[3]+"\t"+one[4]+"\n")
            old.add(one[0]+":"+str(one[1]))
    f1.close()
    f2.close()
    fo.close()
    
##----------------------------------------------------------------------------------------------------------

def combine_res(bed_in_dir1,bed_in_dir2,bed_in_dir3,bed_out_dir):
    f1=open(bed_in_dir1)
    f2=open(bed_in_dir2)
    f3=open(bed_in_dir3)
    fo=open(bed_out_dir,"w")
    whole=[]
    for line in f1:
        seq=line.replace("\n","").split("\t")
        whole.append([seq[0],int(seq[2]),seq[3],seq[4]])
    for line in f2:
        seq=line.replace("\n","").split("\t")
        whole.append([seq[0],int(seq[2]),seq[3],seq[4]])
    for line in f3:
        seq=line.replace("\n","").split("\t")
        whole.append([seq[0],int(seq[2]),seq[3],seq[4]])
    whole.sort()
    for one in whole:
        fo.write(one[0]+"\t"+str(one[1]-1)+"\t"+str(one[1])+"\t"+one[2]+"\t"+one[3]+"\n")
    f1.close()
    f2.close()
    f3.close()
    fo.close()
    
##----------------------------------------------------------------------------------------------------------

def res_or(bed_in_dir1,bed_in_dir2,bed_out_dir):
    f1=open(bed_in_dir1)
    f2=open(bed_in_dir2)
    fo=open(bed_out_dir,"w")
    whole={}
    for line in f1:
        seq=line.rstrip().split("\t")
        whole[":".join(seq[0:4])]=int(seq[4])
    for line in f2:
        seq=line.rstrip().split("\t")
        if ":".join(seq[0:4]) in whole:
            pass
            #whole[":".join(seq[0:4])+":"+seq[5]] +=int(seq[4])
        else:
            whole[":".join(seq[0:4])]=int(seq[4])
    lst=[]
    for one in whole:
        seq=one.split(":")
        lst.append([seq[0],int(seq[1]),int(seq[2]),seq[3],str(whole[one])])
    lst.sort()
    for one in lst:
        out=[]
        for i in one:
            out.append(str(i))
        fo.write("\t".join(out)+"\n")

    f1.close()
    f2.close()
    fo.close()
    
##----------------------------------------------------------------------------------------------------------

def o2b(bed_in,bed_out):
    fi=open(bed_in)
    fo=open(bed_out,"w")
    for line in fi:
        seq=line.rstrip().split("\t")
        #if float(seq[7])<0.2:
        fo.write("\t".join(seq[0:6])+"\n")
        
        
def get_depth(zz_in_dir=0,bed_in_dir=0,bed_out_dir=0):

    fread=open(zz_in_dir)# "./zz_folder/all.zz")
    fsnv=open(bed_in_dir) #"../bed_folder/try_new.bed") #   Hap_SRR521447.bed")
    fo=open(bed_out_dir,"w")#"./tmp/readspersite_new.zer","w")

    class Read:
        def __init__(self,read):
            
            self.snv=read.split("\t")[4].split(";")
            self.inter=read.split("\t")[3].split(";")
            self.direct=read.split("\t")[1]
        def locisin(self,loc):
            isin=0
            for inter in self.inter:
                inter=inter.split(":")
                if int(loc)<=int(inter[1]) and int(loc)>=int(inter[0]):
                    isin =1
                    break
            if isin ==0:
                return 0
            elif isin ==1:
                return 1
        def snvisin(self,snv):
            if snv in self.snv:
                return 1
            else:
                return 0
        def getmin(self):
            return int(self.inter[0].split(":")[0])
        def getmax(self):
            return int(self.inter[-1].split(":")[1])


    reads={}
    for line in fread:
        seq=line.split("\t")
        try:
                reads[seq[0]].append(line[0:-1])
        except Exception as e:
                
                reads[seq[0]]=[line[0:-1]]


    top=0
    chrr=""
    for line in fsnv:
        seq=line.rstrip().split("\t")
        deep=0
        altdeep=0
        try:
            snv=seq[3]+":"+seq[2]
            if seq[0] != chrr:
                top=0
                chrr=seq[0]
            if seq[0] not in reads:
                reads[seq[0]]=[]
            if top < len(reads[seq[0]]):
                while seq[0]==chrr and top < len(reads[seq[0]]) and Read(reads[seq[0]][top]).getmax() < int(seq[2]):
                    top=top+1
                point=top
                while seq[0]==chrr and point < len(reads[seq[0]]) and Read(reads[seq[0]][point]).getmin() <= int(seq[2]):
                    if Read(reads[seq[0]][point]).locisin(seq[2]) ==1:
                        deep=deep+1
                    
                    if Read(reads[seq[0]][point]).snvisin(snv)==1:
                        altdeep=altdeep+1
                
                    point=point+1
            fo.write(line[0:-1]+"\t"+str(altdeep)+":"+str(deep)+"\n")
        except Exception as e:
            print(line)
    fread.close()
    fsnv.close()
    fo.close()

##----------------------------------------------------------------------------------------------------------

def snv_or(bed_in_dir1,bed_in_dir2,bed_out_dir):
    f1=open(bed_in_dir1)
    f2=open(bed_in_dir2)
    fo=open(bed_out_dir,"w")
    whole={}
    for line in f1:
        seq=line.rstrip().split("\t")
        whole[":".join(seq[0:4])]=int(seq[4])
    for line in f2:
        seq=line.rstrip().split("\t")
        if ":".join(seq[0:4]) in whole:
            whole[":".join(seq[0:4])] +=int(seq[4])
        else:
            whole[":".join(seq[0:4])]=int(seq[4])
    lst=[]
    for one in whole:
        seq=one.split(":")
        lst.append([seq[0],int(seq[1]),int(seq[2]),seq[3],str(whole[one])])
    lst.sort()
    for one in lst:
        out=[]
        for i in one:
            out.append(str(i))
        fo.write("\t".join(out)+"\n")

    f1.close()
    f2.close()
    fo.close()

##----------------------------------------------------------------------------------------------------------

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

    all_type=list(all_type_num.keys())
    all_num=list(all_type_num.values())
    all_type_num["ALL_RES"]=all_snv_num

    all_ratio=[round(n/all_snv_num,3) for n in all_num]

    fo=open(summary_path, 'w')
    this_line="Type\tNumber\tRatio\n"
    fo.write(this_line)
    for i in range(len(all_num)):
        this_line=str(all_type[i])+"\t"+str(all_num[i])+"\t"+str(all_ratio[i])+"\n"
        fo.write(this_line)
    fo.close()

##----------------------------------------------------------------------------------------------------------

def getRES1():
    subprocess.Popen("cp " + ori_tmp+"regular.snv" + " " + tmp+"regular.snv",shell=True).wait() 
    subprocess.Popen("cp " + ori_tmp+"hyper_mskAG.snv" + " " + tmp+"hyper_mskAG.snv",shell=True).wait() 
    subprocess.Popen("cp " + ori_tmp+"hyper_mskTC.snv" + " " + tmp+"hyper_mskTC.snv",shell=True).wait() 
    annotate(tmp+"regular.snv",repeat,tmp+"regular.snv.anno")    
    seperate(tmp+"regular.snv.anno",tmp+"regular.snv.anno.alu",tmp+"regular.snv.anno.nalurp",tmp+"regular.snv.anno.nrp","Alu")
    get_snv_with_ad(tmp+"regular.snv.anno.alu",tmp+"regular.snv.anno.alu.ad2",2)
    snv_cluster(tmp+"regular.snv.anno.alu",tmp+"regular_alu.res.ad1", cluster_distance, cluster_size_alu_ad1)
    snv_cluster(tmp+"regular.snv.anno.alu.ad2",tmp+"regular_alu.res.ad2", cluster_distance, cluster_size_alu_ad2)
    bed_or(tmp+"regular_alu.res.ad1",tmp+"regular_alu.res.ad2",tmp+"regular_alu.res")
    snv_cluster(tmp+"regular.snv.anno.nalurp",tmp+"regular_nalurp.res", cluster_distance, cluster_size_nalurp)
    snv_cluster(tmp+"regular.snv.anno.nrp",tmp+"regular_nrp.res", cluster_distance, cluster_size_nrp)
    combine_res(tmp+"regular_alu.res",tmp+"regular_nalurp.res",tmp+"regular_nrp.res",tmp+"regular_split.res")
    cluster_size_regular_max=max([cluster_size_alu_ad1,cluster_size_alu_ad2,cluster_size_nalurp,cluster_size_nrp])
    combine_res(tmp+"regular.snv.anno.alu",tmp+"regular.snv.anno.nalurp",tmp+"regular.snv.anno.nrp",tmp+"regular.snv.anno.rmsrp")
    snv_cluster(tmp+"regular.snv.anno.rmsrp",tmp+"regular_overall.res", cluster_distance, cluster_size_regular_max)
    res_or(tmp+"regular_split.res",tmp+"regular_overall.res",tmp+"regular.res")
    

    annotate(tmp+"hyper_mskTC.snv",repeat,tmp+"hyper_mskTC.snv.anno")    
    seperate(tmp+"hyper_mskTC.snv.anno",tmp+"hyper_mskTC.snv.anno.alu",tmp+"hyper_mskTC.snv.anno.nalurp",tmp+"hyper_mskTC.snv.anno.nrp","Alu")
    snv_cluster(tmp+"hyper_mskTC.snv.anno.alu",tmp+"hyper_mskTC_alu.res", cluster_distance, cluster_size_alu_hp)
    snv_cluster(tmp+"hyper_mskTC.snv.anno.nalurp",tmp+"hyper_mskTC_nalurp.res", cluster_distance, cluster_size_nalurp_hp)
    snv_cluster(tmp+"hyper_mskTC.snv.anno.nrp",tmp+"hyper_mskTC_nrp.res", cluster_distance, cluster_size_nrp_hp)
    combine_res(tmp+"hyper_mskTC_alu.res",tmp+"hyper_mskTC_nalurp.res",tmp+"hyper_mskTC_nrp.res",tmp+"hyper_mskTC_split.res")            
    cluster_size_hyper_max=max([cluster_size_alu_hp,cluster_size_nalurp_hp,cluster_size_nrp_hp])
    combine_res(tmp+"hyper_mskTC.snv.anno.alu",tmp+"hyper_mskTC.snv.anno.nalurp",tmp+"hyper_mskTC.snv.anno.nrp",tmp+"hyper_mskTC.snv.anno.rmsrp")
    snv_cluster(tmp+"hyper_mskTC.snv.anno.rmsrp",tmp+"hyper_mskTC_overall.res", cluster_distance, cluster_size_hyper_max)
    res_or(tmp+"hyper_mskTC_split.res",tmp+"hyper_mskTC_overall.res",tmp+"hyper_mskTC.res")


    annotate(tmp+"hyper_mskAG.snv",repeat,tmp+"hyper_mskAG.snv.anno")    
    seperate(tmp+"hyper_mskAG.snv.anno",tmp+"hyper_mskAG.snv.anno.alu",tmp+"hyper_mskAG.snv.anno.nalurp",tmp+"hyper_mskAG.snv.anno.nrp","Alu")
    snv_cluster(tmp+"hyper_mskAG.snv.anno.alu",tmp+"hyper_mskAG_alu.res", cluster_distance, cluster_size_alu_hp)
    snv_cluster(tmp+"hyper_mskAG.snv.anno.nalurp",tmp+"hyper_mskAG_nalurp.res", cluster_distance, cluster_size_nalurp_hp)
    snv_cluster(tmp+"hyper_mskAG.snv.anno.nrp",tmp+"hyper_mskAG_nrp.res", cluster_distance, cluster_size_nrp_hp)
    combine_res(tmp+"hyper_mskAG_alu.res",tmp+"hyper_mskAG_nalurp.res",tmp+"hyper_mskAG_nrp.res",tmp+"hyper_mskAG_split.res")            
    cluster_size_hyper_max=max([cluster_size_alu_hp,cluster_size_nalurp_hp,cluster_size_nrp_hp])
    combine_res(tmp+"hyper_mskAG.snv.anno.alu",tmp+"hyper_mskAG.snv.anno.nalurp",tmp+"hyper_mskAG.snv.anno.nrp",tmp+"hyper_mskAG.snv.anno.rmsrp")
    snv_cluster(tmp+"hyper_mskAG.snv.anno.rmsrp",tmp+"hyper_mskAG_overall.res", cluster_distance, cluster_size_hyper_max)
    res_or(tmp+"hyper_mskAG_split.res",tmp+"hyper_mskAG_overall.res",tmp+"hyper_mskAG.res")

    snv_or(tmp+"hyper_mskTC.res",tmp+"hyper_mskAG.res",tmp+"hyper.res")

    snv_cluster(tmp+"regular.snv",tmp+"regular.res_tmp",cluster_distance,cluster_size_rg) 
    o2b(tmp+"regular.res_tmp",tmp+"regular.res") 

    snv_cluster(tmp+"hyper_mskTC.snv",tmp+"hyper_mskTC.res",cluster_distance,cluster_size_hp)
    snv_cluster(tmp+"hyper_mskAG.snv",tmp+"hyper_mskAG.res",cluster_distance,cluster_size_hp)
    snv_or(tmp+"hyper_mskTC.res",tmp+"hyper_mskAG.res",tmp+"hyper.res")

##----------------------------------------------------------------------------------------------------------
def getRES2():
    subprocess.Popen("cp " + ori_tmp+"hyper_mskAG.snv" + " " + tmp+"hyper_mskAG.snv",shell=True).wait() 
    subprocess.Popen("cp " + ori_tmp+"hyper_mskTC.snv" + " " + tmp+"hyper_mskTC.snv",shell=True).wait() 
    snv_cluster(tmp+"hyper_mskTC.snv",tmp+"hyper_mskTC.res", cluster_distance, cluster_size_nrp_hp)
    snv_cluster(tmp+"hyper_mskAG.snv",tmp+"hyper_mskAG.res", cluster_distance, cluster_size_nrp_hp)
    snv_or(tmp+"hyper_mskTC.res",tmp+"hyper_mskAG.res",tmp+"hyper.res")


##----------------------------------------------------------------------------------------------------------

def main(args):

    parser = argparse.ArgumentParser(description="python3 step2_annotation_based.py")
    parser.add_argument("-o", "--output",default = str(os.getcwd())+"/2_RES_calling/",  help="path to Output Directory (defalt:Working Directory)")
    parser.add_argument("-s", "--snv",default = str(os.getcwd())+"/1_SNV_calling/",  help="path to SNV Directory (defalt:Working Directory)")
    parser.add_argument("-rp", "--repeat",help="path to repeat file #Optional")

    try:
        args = parser.parse_args(args.split())

    except argparse.ArgumentError as e:
        print(e)
        exit(1)

    global ori_tmp,tmp,repeat
    ori_tmp=str(args.snv)+"/"
    tmp=str(args.output)+"/"
    repeat=str(args.repeat)

    if args.repeat:
        repeat=str(args.repeat)
    else:
        repeat=False

    global cluster_distance,cluster_size_alu_ad1,cluster_size_alu_ad2,cluster_size_nalurp,cluster_size_nrp
    global cluster_size_rg,cluster_size_hp,cluster_size_alu_hp,cluster_size_nalurp_hp,cluster_size_nrp_hp
    global strand_specify,var_limit,poly_limit,rm_multi
    cluster_distance=200
    cluster_size_alu_ad1 = 3
    cluster_size_alu_ad2 = 2
    cluster_size_nalurp = 5
    cluster_size_nrp = 7
    cluster_size_rg = 5
    cluster_size_hp = 5
    cluster_size_alu_hp = 5
    cluster_size_nalurp_hp = 5
    cluster_size_nrp_hp = 5
    strand_specify=0
    var_limit=20
    poly_limit=10
    rm_multi=0

    ##----------------------------------------------------------------------------------------------------------

    print("Step2 Identifying RESs based repeat annotation...")
    time_start = int(time.perf_counter())
    if not os.path.isdir(tmp):
        os.mkdir(tmp)
    if repeat != False:
        getRES1()    
        if not os.path.isdir(tmp+"anno_regular"):
            os.mkdir(tmp+"anno_regular")
        subprocess.Popen("mv " + tmp+"regular* " + tmp+"anno_regular/",shell=True).wait() 
        cmd="awk -F '\t' '{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5}' "
        subprocess.Popen(cmd + tmp+"anno_regular/regular.res > " + tmp+"anno_regular.res",shell=True).wait() 
        ARG=[tmp+"anno_regular.res",tmp+"anno_regular_res.stat"]
        stat_res(ARG)

    else:
        getRES2()

    if not os.path.isdir(tmp+"hyper"):
        os.mkdir(tmp+"hyper")
    
    subprocess.Popen("mv " + tmp+"hyper_* " + tmp+"hyper/",shell=True).wait() 
    subprocess.Popen("mv " + tmp+"hyper.* " + tmp+"hyper/",shell=True).wait() 
    subprocess.Popen("cp " + tmp+"hyper/hyper.res " + tmp+"hyper.res",shell=True).wait() 

    ARG=[tmp+"hyper.res",tmp+"SPRINT2_Hyper_RES.stat"]
    stat_res(ARG)
    time_end = int(time.perf_counter())
    run_time = str(time_end-time_start)

    print("Repeat-based RES calling running time is:" + run_time + " seconds" + "\n")     

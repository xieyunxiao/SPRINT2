##SNV calling
import subprocess,os,sys,pysam
import time, pysam,re
from Bio import SeqIO 
import argparse

def description():
    print('')
    print("##############################################################################################")
    print('')
    print("   SPRINT-2.0: An enhanced tool to identify RNA editing sites")
    print('')
    print("   https://github.com/xieyunxiao/SPRINT2/")
    print('')
    print("   Please contact 20210700107@fudan.edu.cn when questions arise.")
    print('')
    print("##############################################################################################")

def get_baseq_cutoff(fq_in_dir=0,cutoff_out_dir=0):
    fi=open(fq_in_dir)
    fo=open(cutoff_out_dir,"w")
    line1=fi.readline()
    line2=fi.readline()
    line3=fi.readline()
    line4=fi.readline()
    did=0
    while line1 !="":
        if did==1:
            break
        qua=line4[0:-1]
        for i in qua:
            if ord(i) > 76:
                fo.write("89") #64+25
                did=1
                break
            if ord(i) < 60:
                fo.write("58") #33+25
                did=1
                break

        line1=fi.readline()
        line2=fi.readline()
        line3=fi.readline()
        line4=fi.readline()
    fi.close()
    fo.close()    

##----------------------------------------------------------------------------------------------------------

def cut(fq_in_dir=0,fq_out_dir=0,cutnum=0,name="read1"):
    fi=open(fq_in_dir)
    fo=open(fq_out_dir,"w")
    cutnum=int(cutnum)
    line1=fi.readline()
    line2=fi.readline()
    line3=fi.readline()
    line4=fi.readline()
    idd=1
    while line1 !="":
        CELL_TAG=""
        if "XC:Z:" in line1:
            seq=line1.split("_")
            for one in seq:
                if one[:5]=="XC:Z:":
                    CELL_TAG=one        
        if CELL_TAG !="":
            fo.write("@id_"+str(idd)+"_"+CELL_TAG+"_"+name+"\n")
        else:
            fo.write("@id_"+str(idd)+"_"+name+"\n")
        fo.write(line2[cutnum:])
        fo.write(line3)
        fo.write(line4[cutnum:])
        line1=fi.readline()
        line2=fi.readline()
        line3=fi.readline()
        line4=fi.readline()
        idd=idd+1
    fi.close()
    fo.close()    

##----------------------------------------------------------------------------------------------------------

def fq2bam(out_dir,paired_end,TAG,read1,read2,refgenome,bwa,samtools,mapcpu):
    if paired_end==True:
        mapcpu=max([int(int(mapcpu)/2.0),1])
    
    step1_1=subprocess.Popen(bwa+" aln -t "+str(mapcpu)+" "+refgenome+" "+read1+" > "+out_dir+"read1.sai",shell=True)
    step1_1.wait()
    if paired_end==True:
        step1_2=subprocess.Popen(bwa+" aln  -t "+str(mapcpu)+" "+refgenome+" "+read2+" > "+out_dir+"read2.sai",shell=True)
        step1_2.wait()
        
    step1_3=subprocess.Popen(bwa+" samse -n4 "+refgenome+" "+out_dir+"read1.sai "+read1+" > "+out_dir+"name_read1.sam",shell=True)
    step1_3.wait()
    os.remove(out_dir+"read1.sai")
    if paired_end==True:
        step1_4=subprocess.Popen(bwa+" samse -n4 "+refgenome+" "+out_dir+"read2.sai "+read2+" > "+out_dir+"name_read2.sam",shell=True)
        step1_4.wait()
        os.remove(out_dir+"read2.sai")        
  
    step1_7=subprocess.Popen(samtools+" view -@ "+str(mapcpu)+" -bS "+out_dir+"name_read1.sam >"+out_dir+"name_read1.bam",shell=True)
    step1_7.wait()
    if paired_end==True:
        step1_8=subprocess.Popen(samtools+" view -@ "+str(mapcpu)+" -bS "+out_dir+"name_read2.sam >"+out_dir+"name_read2.bam",shell=True)
        step1_8.wait()
        step1_9=subprocess.Popen(samtools+" sort -@ "+str(mapcpu)+" "+out_dir+"name_read1.bam "+out_dir+"name_read1_sorted",shell=True)
        step1_10=subprocess.Popen(samtools+" sort -@ "+str(mapcpu)+" "+out_dir+"name_read2.bam "+out_dir+"name_read2_sorted",shell=True)
        step1_9.wait()
        step1_10.wait()
        step1_11=subprocess.Popen(samtools+" merge -f "+out_dir+"all.bam "+out_dir+"name_read1_sorted.bam "+out_dir+"name_read2_sorted.bam",shell=True)
        step1_11.wait()

        if os.path.exists(out_dir+"all.bam"):
                if os.path.exists(out_dir+"name_read1.sam"):
                        os.remove(out_dir+"name_read1.sam")
                if os.path.exists(out_dir+"name_read1.bam"):
                    os.remove(out_dir+"name_read1.bam")
                if os.path.exists(out_dir+"name_read1_sorted.bam"):
                        os.remove(out_dir+"name_read1_sorted.bam")
                if os.path.exists(out_dir+"name_read2.sam"):
                        os.remove(out_dir+"name_read2.sam")
                if os.path.exists(out_dir+"name_read2.bam"):
                        os.remove(out_dir+"name_read2.bam")
                if os.path.exists(out_dir+"name_read2_sorted.bam"):
                        os.remove(out_dir+"name_read2_sorted.bam")
    
    else:
        step1_9=subprocess.Popen(samtools+" sort -@ "+str(mapcpu)+" "+out_dir+"name_read1.bam "+out_dir+"all",shell=True)
        step1_9.wait()
        if os.path.exists(out_dir+"all.bam"):
            if os.path.exists(tmp+"name_read1.sam"):
                os.remove(tmp+"name_read1.sam")
            if os.path.exists(tmp+"name_read1.bam"):
                os.remove(tmp+"name_read1.bam")

    step2_2=subprocess.Popen(samtools+" view -h -o "+out_dir+"all.sam "+out_dir+"all.bam",shell=True)
    step2_2.wait()
    subprocess.Popen("cp "+out_dir+"/all.sam "+tmp+"/"+TAG+"_all.sam",shell=True).wait()
    if os.path.exists(out_dir+"all.sam"):
        os.remove(out_dir+"all.sam")

##----------------------------------------------------------------------------------------------------------

def fq2sam(TAG,paired_end,read1,read2,tmp,refgenome,bwa,samtools,mapcpu,):
    if paired_end==True:
        mapcpu=max([int(int(mapcpu)/2.0),1])
    ori_tmp=tmp
    tmp=tmp+"/"+TAG+"/"
    if os.path.exists(tmp)==False:
        os.mkdir(tmp)
        
    step1_1=subprocess.Popen(bwa+" aln -t "+str(mapcpu)+" "+refgenome+" "+read1+" > "+tmp+"read1.sai",shell=True)
    if paired_end==True:
        step1_2=subprocess.Popen(bwa+" aln  -t "+str(mapcpu)+" "+refgenome+" "+read2+" > "+tmp+"read2.sai",shell=True)

    step1_1.wait()
    if paired_end==True:
        step1_2.wait()
    step1_3=subprocess.Popen(bwa+" samse -n4 "+refgenome+" "+tmp+"read1.sai "+read1+" > "+tmp+"name_read1.sam",shell=True)
    if paired_end==True:
        step1_4=subprocess.Popen(bwa+" samse -n4 "+refgenome+" "+tmp+"read2.sai "+read2+" > "+tmp+"name_read2.sam",shell=True)
    step1_3.wait()
    if paired_end==True:
        step1_4.wait()
    if os.path.exists(tmp+"name_read1.sam"):
        if os.path.exists(tmp+"read1.sai"):
            os.remove(tmp+"read1.sai")
        if os.path.exists(ori_tmp+"cut_read1.fastq"):
            os.remove(ori_tmp+"cut_read1.fastq")
    if os.path.exists(tmp+"name_read2.sam"):
        if os.path.exists(tmp+"read2.sai"):
            os.remove(tmp+"read2.sai")
        if os.path.exists(ori_tmp+"cut_read2.fastq"):
            os.remove(ori_tmp+"cut_read2.fastq")
    
    step1_7=subprocess.Popen(samtools+" view -bS "+tmp+"name_read1.sam >"+tmp+"name_read1.bam",shell=True)
    if paired_end==True:
        step1_8=subprocess.Popen(samtools+" view -bS "+tmp+"name_read2.sam >"+tmp+"name_read2.bam",shell=True)
    step1_7.wait()
    if paired_end==True:
        step1_8.wait()
    if paired_end==True:
        step1_9=subprocess.Popen(samtools+" sort "+tmp+"name_read1.bam "+tmp+"name_read1_sorted",shell=True)
        step1_10=subprocess.Popen(samtools+" sort "+tmp+"name_read2.bam "+tmp+"name_read2_sorted",shell=True)
        step1_9.wait()
        step1_10.wait()
        step1_11=subprocess.Popen(samtools+" merge -f "+tmp+"all.bam "+tmp+"name_read1_sorted.bam "+tmp+"name_read2_sorted.bam",shell=True)
        step1_11.wait()
        if os.path.exists(tmp+"all.bam"):
            if os.path.exists(tmp+"name_read1.sam"):
                os.remove(tmp+"name_read1.sam")
            if os.path.exists(tmp+"name_read1.bam"):
                os.remove(tmp+"name_read1.bam")
            if os.path.exists(tmp+"name_read1_sorted.bam"):
                os.remove(tmp+"name_read1_sorted.bam")
            if os.path.exists(tmp+"name_read2.sam"):
                os.remove(tmp+"name_read2.sam")
            if os.path.exists(tmp+"name_read2.bam"):
                os.remove(tmp+"name_read2.bam")
            if os.path.exists(tmp+"name_read2_sorted.bam"):
                os.remove(tmp+"name_read2_sorted.bam")

    else:
        step1_9=subprocess.Popen(samtools+" sort "+tmp+"name_read1.bam "+tmp+"all",shell=True)
        step1_9.wait()
        if os.path.exists(tmp+"all.bam"):
            if os.path.exists(tmp+"name_read1.sam"):
                os.remove(tmp+"name_read1.sam")
            if os.path.exists(tmp+"name_read1.bam"):
                os.remove(tmp+"name_read1.bam")
    step2_2=subprocess.Popen(samtools+" view -h -o "+tmp+"all.sam "+tmp+"all.bam",shell=True)
    step2_2.wait()
    subprocess.Popen("cp "+tmp+"./all.sam "+ori_tmp+"/"+TAG+"_all.sam",shell=True).wait()
    if os.path.exists(tmp+"all.sam"):
        os.remove(tmp+"all.sam")
    
##----------------------------------------------------------------------------------------------------------

def bam2SNV(bam_in_path=0, fa_in_path=0, zz_out_path=0,bed_out_path=0):
    fi = pysam.AlignmentFile(bam_in_path, "rb")  # bam
    fo = open(zz_out_path, "w")  # zz
    fo.write("#chr\tsam_flag\tsam_mapq\tread_interval\tmismatch_list\tmismatch_quality\tmismatch_read_pos\tseq\tread_name\tfragment_end_dist\n")
    fo2=open(bed_out_path,"w")
    fo2.write("#chr\tstart_0\tend_1\tmismatch\tquality\tfragment_end_dist\n")    
    for line in fi:
        if line.alen==None:
            continue
        this_name = str(line.query_name)  # read_name
        this_chr = str(line.reference_name)  # chr
        this_CG = str(line.cigarstring)  #        
        this_pos = str(line.pos + 1)  #
        this_seq = str(line.seq)  # seq
        this_qseq = str(line.qual)  #
        this_mapq = str(line.mapq)  # sam_mapq
        this_flag = str(line.flag)  # sam_flag
        this_seqlen = int(line.alen) 
        mismatch_list = ""
        mismatch_quality = ""
        mismatch_read_pos = ""
        fragment_end_dist = ""
        this_interval_list=[]
        for a in line.aligned_pairs:
            seq_ind = a[0]
            ref_ind = a[1]
            if ref_ind != None:
                this_interval_list.append(int(ref_ind))
            if seq_ind != None and ref_ind != None:
                seq_a = this_seq[seq_ind].upper()
                ref_a = chrom[this_chr][ref_ind].upper()
                if seq_a != ref_a and ref_a != "N" and seq_a != "N":
                    mis_a = ref_a + seq_a + ":" + str(int(ref_ind) + 1)
                    mismatch_list = mismatch_list + mis_a + ";"
                    mismatch_quality = mismatch_quality + str(ord(this_qseq[seq_ind])) + ","
                    mismatch_read_pos = mismatch_read_pos + str(int(seq_ind) + 1) + ","
                    fragment_end_dist1=min(int(this_seqlen) - int(seq_ind),int(seq_ind) + 1)
                    fragment_end_dist = fragment_end_dist + str(fragment_end_dist1) + ","
                    this_out2=this_chr + "\t" + str(int(ref_ind)) + "\t" + str(int(ref_ind) + 1) + "\t" + ref_a + seq_a + "\t" + str(ord(this_qseq[seq_ind])) + "\t" + str(fragment_end_dist1) + "\n"
                    fo2.write(this_out2)
        this_interval = str(min(this_interval_list)+1) + ":" + str(max(this_interval_list)+1)# read_interval
        this_out = this_chr + "\t" + this_flag + "\t" + this_mapq + "\t" + this_interval + "\t" + mismatch_list[0:-1] + "\t" + mismatch_quality[0:-1] + "\t" + mismatch_read_pos[0:-1] + "\t" + this_seq + "\t" + this_name + "\t" + fragment_end_dist[0:-1] + "\n"
        if mismatch_list != "":
            fo.write(this_out)
    fi.close()
    fo.close()
    fo2.close()

##----------------------------------------------------------------------------------------------------------

def umsam2fq(sam_in_dir=0,fq_out_dir=0):

    fi=open(sam_in_dir)
    fo=open(fq_out_dir,"w")
    for line in fi:
        seq=line.rstrip().split("\t")
        if line[0] !="@" and len(bin(int(seq[1])))>=5 and bin(int(seq[1]))[-3]=="1":
            if len(bin(int(seq[1])))>=9 and bin(int(seq[1]))[-7]=="1":
                seq[0]=seq[0][0:-2]+"_1"
            elif len(bin(int(seq[1])))>=10 and bin(int(seq[1]))[-8]=="1":
                seq[0]=seq[0][0:-2]+"_2"
            elif line[-1]=="1":
                seq[0]=seq[0][0:-2]+"_1"
            elif line[-1]=="2":
                seq[0]=seq[0][0:-2]+"_2"
            fo.write("@"+seq[0]+"\n"+seq[9]+"\n+\n"+seq[10]+"\n")
    fo.close()
    fi.close()

##----------------------------------------------------------------------------------------------------------

def antisense_reverse(read):
    read=read.upper()
    read_change_base=""
    for one in read:
        if one == "A":
            read_change_base += "T"
        elif one == "C":
            read_change_base += "G"
        elif one == "G":
            read_change_base += "C"
        elif one == "T":
            read_change_base += "A"
        else:
            read_change_base += "N"
    read_reverse=read_change_base[::-1]
    return read_reverse

##----------------------------------------------------------------------------------------------------------

def maskfq(fq_in_dir,mask_from,mask_to):
    mask_from=mask_from.upper()
    mask_to=mask_to.upper()
    fi=open(fq_in_dir)
    fo=open(fq_in_dir[0:-3]+"_"+mask_from+"_to_"+mask_to+".fq","w")
    line1=fi.readline().replace("\n","")
    line2=fi.readline().upper().replace("\n","")
    line3=fi.readline().replace("\n","")
    line4=fi.readline().replace("\n","")
    while line1 !="":
        if line1[-1]=="1":
            line2=antisense_reverse(line2)
            line4=line4[::-1]
        record="1"
        line2_new=""
        for one in line2.replace("\n",""):
            if one==mask_from:
                record=record+"1"
            else:
                record=record+"0"
        
        fo.write(line1+"_|_"+mask_from+"_to_"+mask_to+"_|_"+str(int(record,2))+"_|_read2"+"\n")
        fo.write(line2.replace(mask_from,mask_to)+"\n")            
        fo.write("+\n")
        fo.write(line4+"\n")
        line1=fi.readline().replace("\n","")
        line2=fi.readline().upper().replace("\n","")
        line3=fi.readline().replace("\n","")
        line4=fi.readline().replace("\n","")
    fi.close()
    fo.close()

##----------------------------------------------------------------------------------------------------------

def poly_check(seq,poly_limit):
    if "A"*poly_limit not in seq and "T"*poly_limit not in seq and "G"*poly_limit not in seq and "C"*poly_limit not in seq:
        return True
    else:
        return False

##----------------------------------------------------------------------------------------------------------

def var_check(change_from, change_to, seq,var_limit):
    ALL=["A","T","C","G"]
    tmp=[]
    for one in ALL:
        if one !=change_from.upper() and one !=change_to.upper():
            tmp.append(one)
    flag=1
    for one in tmp:
        if seq.count(one) < var_limit/(float(len(tmp))+2):
            flag=0
    if flag==1:
        return True
    else:
        return False
    
##----------------------------------------------------------------------------------------------------------

def reverse_base(base):
    base=base.upper()
    if base=="A":
        return "T"
    elif base=="C":
        return "G"
    elif base=="G":
        return "C"
    elif base=="T":
        return "A"
    else:
        return "N"    
        
##----------------------------------------------------------------------------------------------------------

def recover_sam(sam_in_dir,sam_out_dir, var_limit=20,poly_limit=10,rm_multi=0):
    fi=open(sam_in_dir)
    fo=open(sam_out_dir,"w")
    for line in fi:
        seq=line.split("\t")
        if line[0]=="@":
            fo.write(line)
        elif seq[1]=="4" and seq[2]=="*":
            break
        elif seq[1]!="4" and len(seq)>=9:
            seq=line.split("\t")
            seq[9]=seq[9].upper()
            seq[1]=int(seq[1])
            if len(bin(seq[1]))>=7:
                if bin(seq[1])[-3]!="1":
                    if bin(seq[1])[-5]=="1":
                        seq[1]="16"
                    else:
                        seq[1]="0"
            seq[1]=str(seq[1])
            record=bin(int(seq[0].split("_|_")[2]))[3:]
            change_from=seq[0].split("_|_")[1].split("_")[0]
            change_to=seq[0].split("_|_")[1].split("_")[2]
            
            
            if len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=="1":   #seq[1]=="16":
                change_from=reverse_base(change_from)
                change_to=reverse_base(change_to)
                record=record[::-1]
            else:
                record=record
            changed_read=seq[9]
            i=0
            recovered_read=""
            while i<len(seq[9]):
                if record[i]=="1" and seq[9][i]==change_to:
                    recovered_read += change_from
                elif record[i]=="1" and seq[9][i]!=change_to:
                    #print "recover error in "+seq[0]
                    recovered_read += "N"
                else:
                    recovered_read += seq[9][i]
                i=i+1
            seq[9]=recovered_read
            #fo.write(seq[0])
            #if len(record)==len(seq[9]) and "I" not in seq[5] and "D" not in seq[5] and len(record)-changed_read.count(change_to) > 25 and poly_check(seq[9],poly_limit):
            if len(record)==len(seq[9]) and len(record)-changed_read.count(change_to) > var_limit and poly_check(seq[9],poly_limit): #and var_check(change_from,change_to,seq[9],var_limit):
                if (rm_multi==1 and "XA:Z:" not in line) or rm_multi==0:
                    fo.write(seq[0])
                    j=1
                    while j<len(seq):
            
                        fo.write("\t"+seq[j])
                        j=j+1
    fo.close()
    fi.close()
    
##----------------------------------------------------------------------------------------------------------

def sam2zz(sam_in_dir=0,fa_in_dir=0,zz_out_dir=0):
    #-----------------------------------------------------------
    #reading refgenome
    fref=open(fa_in_dir)
    chrom={}
    chrr=""
    line=fref.read()
    line=line.split(">")
    for seq in line:
            if " " in seq:
                chrr=seq[0:seq.find(" ")]
            else:
                chrr=seq[0:seq.find("\n")]
            chrom[chrr]=seq[seq.find("\n"):].replace("\n","")

    #0 base

    fref.close()
    #------------------------------------------------------------
    #------------------------------------------------------------
    #donee:dolist;lst:donumlist
    def doCG(a):
        donee=[]
        lst=re.findall( "(\d+|\+|-|\*|/)", a )
        for i in a:
            if i == "I" or i == "D" or i== "M" or i=="S" or i=="P" or i=="N" :
                donee.append(i)
        return donee,lst
    #donefunction

    def doneCG(CG,chrr,pos,seq,qseq):#pos is 1 base
        donee,lst=doCG(CG)
        errorsite=""
        intersite=""
        quasite=""
        locsite=""
        pieceloc=""
        refseq=""
        seqseq=""
        refpos=int(pos)-1
        seqpos=0
        step=0
        while step<len(donee):
            if donee[step]=="I":
                seqpos=seqpos+int(lst[step])
            elif donee[step]=="D":
                refpos=refpos+int(lst[step])
            elif donee[step]=="N":
                refpos=refpos+int(lst[step])
            elif donee[step]=="S":
                seqpos=seqpos+int(lst[step])
            elif donee[step]=="M":
                
                refseq=refseq+chrom[chrr][refpos:refpos+int(lst[step])]
                seqseq=seqseq+seq[seqpos:seqpos+int(lst[step])]
                j=refpos
                jj=seqpos
                while j<refpos+int(lst[step]):
                    try:
                        if chrom[chrr][j].upper() != seq[jj].upper() and chrom[chrr][j].upper() !="N" and seq[jj].upper() != "N":
                            errorsite=errorsite+chrom[chrr][j].upper()+seq[jj].upper()+":"+str(j+1)+";"
                            quasite=quasite+","+str(ord(qseq[jj]))
                            locsite=locsite+","+str(jj+1)
                            pieceloc=pieceloc+","+str( min( jj+1-seqpos,seqpos+int(lst[step])-jj  ) )
                                     
                    except Exception as e:
                        pass
                        #print "error with",chrr,pos,e
                    j=j+1
                    jj=jj+1
                intersite=intersite+str(refpos+1)+":"+str(refpos+int(lst[step]))+";"
                refpos=refpos+int(lst[step])
                seqpos=seqpos+int(lst[step])
            step=step+1
        refseq=refseq.upper()
        seqseq=seqseq.upper()
        return refseq,seqseq,errorsite,intersite,quasite,locsite,pieceloc
    #------------------------------------------------------------
    #------------------------------------------------------------            
    fi=open(sam_in_dir) #sam
    fo=open(zz_out_dir,"w") #zz
    fo.write("#chr\tsam_flag\tsam_mapq\tread_interval\tmismatch_list\tmismatch_quality\tmismatch_read_pos\tseq\tread_name\tfragment_end_dist\n") 
    for line in fi:
        seq=line.split("\t")
        if line[0]!="@" and len(seq)>5:
            name=seq[0].split("_|_")[0]
            if seq[0][0]!="@" and seq[2]!="*" and seq[5]!="*" :# and whole[name]<2:
                refseq,seqseq,errorsite,intersite,quasite,locsite,pieceloc=doneCG(seq[5],seq[2],seq[3],seq[9],seq[10])
                quasite=quasite[1:]
                locsite=locsite[1:]
                pieceloc=pieceloc[1:]
                
                if len(intersite[0:-1])==0:
                    intersite="*;"
                if len(errorsite[0:-1])==0:
                    errorsite="*;"
                if len(quasite)==0:
                    quasite="*"
                if len(locsite)==0:
                    locsite="*"
                if len(pieceloc)==0:
                    pieceloc="*"
                
                if len(bin(int(seq[1])))>=9 and bin(int(seq[1]))[-7]=="1":
                    fo.write(seq[2]+"\t"+seq[1]+"\t"+seq[4]+"\t"+intersite[0:-1]+"\t"+errorsite[0:-1]+"\t"+quasite+"\t"+locsite+"\t"+seq[9]+"\t"+seq[0]+"_1"+"\t"+pieceloc+"\n")
                elif len(bin(int(seq[1])))>=10 and bin(int(seq[1]))[-8]=="1":
                    fo.write(seq[2]+"\t"+seq[1]+"\t"+seq[4]+"\t"+intersite[0:-1]+"\t"+errorsite[0:-1]+"\t"+quasite+"\t"+locsite+"\t"+seq[9]+"\t"+seq[0]+"_2"+"\t"+pieceloc+"\n")
                else:
                    fo.write(seq[2]+"\t"+seq[1]+"\t"+seq[4]+"\t"+intersite[0:-1]+"\t"+errorsite[0:-1]+"\t"+quasite+"\t"+locsite+"\t"+seq[9]+"\t"+seq[0]+"\t"+pieceloc+"\n")
                    
    fi.close()    
    fo.close()
    
##----------------------------------------------------------------------------------------------------------

def zz2snv(zz_in_dir=0,bed_out_dir=0,baseq_cutoff_dir=0):        
    fi=open(zz_in_dir)
    fo=open(bed_out_dir,"w")

    fqua=open(baseq_cutoff_dir)
    limitbasequa=int(fqua.readline().replace("\n",""))
    fqua.close()
    limitad=1
    limitloc=5
    limitmpqua=0
    allsnv={}
    for line in fi:
        if line[0]=="#":
            continue
        truesnv=[]
        seq=line.rstrip().split("\t")
        mismatch=seq[4].split(";")
        basequa=seq[5].split(",")
        loc=seq[9].split(",") #fragment-loc
        mpqua=int(seq[2])
        ####################################################################
         #change the sam flag "seq[1]" when you didn"t use "bwa -aln" as mapper
        seq[1]=int(seq[1])
        if len(bin(seq[1]))>=7:
            if bin(seq[1])[-3]!="1":
                if bin(seq[1])[-5]=="1":
                    seq[1]="16"
                else:
                    seq[1]="0"
        elif len(bin(seq[1]))>=5:
            if bin(seq[1])[-3]!="1":
                seq[1]="0"
                             
        else: 
            seq[1]="0"
        #####################################################################        
        if basequa[0]!="*" and mpqua >= limitmpqua and mpqua < 200:
            i=0
            baseqlst=[]
            mistype={}
            while i < len(basequa):
                baseqlst.append(int(basequa[i]))
                try:
                    mistype[mismatch[i].split(":")[0]] += 1
                except Exception as e:
                    mistype[mismatch[i].split(":")[0]] = 1
                if int(basequa[i]) >= limitbasequa  and  int(loc[i]) > limitloc  :

                    truesnv.append([mismatch[i],seq[1],seq[8]])
                    
                i=i+1

            #fflag=1
            #masktype_tmp=seq[8].split("_|_")[1].split("_to_")
            #masktype=masktype_tmp[0]+masktype_tmp[1]
            miss=[]
            for mis in mistype:
                miss.append(mistype[mis])
            miss.sort()
            if len(miss)>=2:
                missnum=sum(miss[:-1])
            else:
                missnum=0
            
            if len(baseqlst)>0 and sum(baseqlst)/float(len(baseqlst)) >= limitbasequa: #and missnum <= mismatch_num(len(seq[7])):
                for snv in truesnv:
                    try:
                        allsnv[seq[0]+"\t"+snv[0].split(":")[0]+"\t"+snv[0].split(":")[1]][0]=allsnv[seq[0]+"\t"+snv[0].split(":")[0]+"\t"+snv[0].split(":")[1]][0]+1
                        if (len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=="1" and snv[2][-1] == "1" ) or ((len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-5]=="0") and snv[2][-1] == "2" ):
                            allsnv[seq[0]+"\t"+snv[0].split(":")[0]+"\t"+snv[0].split(":")[1]][1]=allsnv[seq[0]+"\t"+snv[0].split(":")[0]+"\t"+snv[0].split(":")[1]][1]+1
                        elif (len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=="1" and snv[2][-1] == "2" ) or ((len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-5]=="0") and snv[2][-1] == "1" ):
                            allsnv[seq[0]+"\t"+snv[0].split(":")[0]+"\t"+snv[0].split(":")[1]][2]=allsnv[seq[0]+"\t"+snv[0].split(":")[0]+"\t"+snv[0].split(":")[1]][2]+1
                    except Exception as e:
                        if (len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=="1" and snv[2][-1] == "1" ) or ((len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-5]=="0") and snv[2][-1] == "2" ):
                            allsnv[seq[0]+"\t"+snv[0].split(":")[0]+"\t"+snv[0].split(":")[1]]=[1,1,0]
                        elif (len(bin(int(seq[1]))) > 5 and bin(int(seq[1]))[-5]=="1" and snv[2][-1] == "2" ) or ((len(bin(int(seq[1]))) < 5 or bin(int(seq[1]))[-5]=="0") and snv[2][-1] == "1" ):
                            allsnv[seq[0]+"\t"+snv[0].split(":")[0]+"\t"+snv[0].split(":")[1]]=[1,0,1]

    snv_bed=[]
    for snv in allsnv:
        seq=snv.split("\t")
        if allsnv[snv][0]>=limitad:
                if allsnv[snv][1] > allsnv[snv][2]: 
                    snv_bed.append([seq[0],int(seq[2]),seq[1],"+",allsnv[snv][0]])
                elif allsnv[snv][2] > allsnv[snv][1]:
                    snv_bed.append([seq[0],int(seq[2]),seq[1],"-",allsnv[snv][0]])
                else:
                    snv_bed.append([seq[0],int(seq[2]),seq[1],".",allsnv[snv][0]])

    snv_bed.sort()
    for one in snv_bed:
        fo.write(one[0]+"\t"+str(one[1]-1)+"\t"+str(one[1])+"\t"+one[2]+"\t"+str(one[4])+"\t"+one[3]+"\n")
    fi.close()
    fo.close()
    
##----------------------------------------------------------------------------------------------------------

def tsnv2gsnv(trans_loc_file_in , transcript_snv_in , genome_snv_out):
    fa=open(trans_loc_file_in)
    l1=fa.readline()
    l2=fa.readline()
    TRANS={}
    while l1 !="":
        trans=l1[1:].rstrip()
        seq=l2.split(";")[:-1]
        TRANS[trans]=seq
        l1=fa.readline()
        l2=fa.readline()

    fa.close()

    def loc_t2g(tCHR,tLOC):
        gCHR=tCHR.split("_|_")[1]
        seq=TRANS[tCHR]
        tLOC=int(tLOC)
        flag=1
        tmp=0
        i=0
        while i<len(seq) and flag==1:
            end=int(seq[i].split(",")[1])
            start=int(seq[i].split(",")[0])
            tmp +=  end-start+1
            if tmp >= tLOC :
                flag=0
            i += 1
        j=i-1
        dis2end  = tmp-tLOC
        gLOC = int(seq[j].split(",")[1]) - dis2end
        return gLOC



    fi=open(transcript_snv_in)
    fo=open(genome_snv_out,"w")
    for line in fi:
        if line[0]=="#":
            fo.write(line)
            continue
        seq=line.split("\t")
        
        tCHR = seq[0]
        gCHR=tCHR.split("_|_")[1]
        out=loc_t2g(tCHR,seq[2])
        seq[0]=gCHR
        seq[1]=str(int(out)-1)
        seq[2]=str(int(out))
        this_write="\t".join(seq)

        fo.write(this_write)
    fi.close()
    fo.close()

##----------------------------------------------------------------------------------------------------------

def getBedAD(bed_in_path=0, allsnv_in_path=0, bed_out_path=0):
    AD={}

    fi=open(allsnv_in_path)
    for line in fi:
        line=line.rstrip()
        if line[0]=="#":
            continue
        seq=line.split("\t")
        this_snv=":".join([seq[0],seq[2],seq[3]])
        AD[this_snv]=0
    fi.close()

    fi=open(allsnv_in_path)
    for line in fi:
        line=line.rstrip()
        if line[0]=="#":
            continue
        seq=line.split("\t")
        this_snv=":".join([seq[0],seq[2],seq[3]])
        AD[this_snv]=AD[this_snv]+1
    fi.close()


    fi=open(bed_in_path)
    fo=open(bed_out_path,"w")
    for line in fi:
        line=line.rstrip()
        if line[0]=="#":
            continue
        seq=line.rstrip().split("\t")
        this_snv=":".join([seq[0],seq[2],seq[3]])
        this_ad=str(AD[this_snv])
        fo.write("\t".join(seq[0:4])+"\t"+this_ad+"\n")
    fo.close()
##----------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------
def main(args):

    parser = argparse.ArgumentParser(description="python3 step1_SNVcalling.py")
    parser.add_argument("-o", "--output",default = str(os.getcwd()),  help="path to Output Directory (defalt:Working Directory)")
    parser.add_argument("-r1", "--read1",default = str(os.getcwd())+"/test_1.fastq",  help="path to read1.fastq (defalt:~/test_1.fastq)")
    parser.add_argument("-r2", "--read2",  help="path to read2.fastq ##optional")
    parser.add_argument("-R", "--reference",default = str(os.getcwd())+"/reference.fa",  help="path to reference FASTA file (defalt:~/reference.fa)")
    parser.add_argument("-b", "--bwa",default = str(os.getcwd())+"/bwa",  help="path to BWA (defalt:~/bwa)")
    parser.add_argument("-B", "--bedtools",default = str(os.getcwd())+"/bedtools",  help="path to BEDTOOLS ((defalt:~/bedtools)")
    parser.add_argument("-s", "--samtools",default = str(os.getcwd())+"/samtools",  help="path to SAMTOOLS  ((defalt:~/samtools)")
    parser.add_argument("-p", "--cpu",default=1,  help="CPU number (defalt=1)")

    try:
        args = parser.parse_args(args.split())

    except argparse.ArgumentError as e:
        print(e)
        exit(1)

    ##----------------------------------------------------------------------------------------------------------

    global cutbp,paired_end,read0,tmp,chrom
    cutbp=6
    paired_end=True
    read0=""
    tmp=str(args.output)+"/"
    read1=str(args.read1)

    refgenome=str(args.reference)
    bwa=str(args.bwa)
    samtools=str(args.samtools)
    bedtools=str(args.bedtools)
    mapcpu=int(args.cpu)

    if args.read2:
        read2=str(args.read2)
        paired_end=True
    else:
        read2=read0
        paired_end=False


    ##----------------------------------------------------------------------------------------------------------
    description()
    print("\n\n\nStep1:SNV calling...")

    time_start = int(time.perf_counter())

    if not os.path.isdir(tmp):
        os.mkdir(tmp)
    ##read_cut
    cut(read1,tmp+"cut_read1.fastq",cutbp,"read1")
    if paired_end==True:
        cut(read2,tmp+"cut_read2.fastq",cutbp,"read2")
    get_baseq_cutoff(read1,tmp+"baseq.cutoff")

    ##----------------------------------------------------------------------------------------------------------

    ##----------------------------------------------------------------------------------------------------------
    ##genome
    print("\tStart mapping genome ...")
    if not os.path.isdir(tmp+"genome/"):
        os.mkdir(tmp+"genome/")
    TAG="genome"
    fq2bam(tmp+"genome/",paired_end,TAG,read1,read2,refgenome,bwa,samtools,mapcpu)
    subprocess.Popen(samtools+" index "+tmp+"genome/all.bam",shell=True).wait()

    os.remove(tmp+"cut_read1.fastq")
    if os.path.exists(tmp+"cut_read2.fastq"):
        os.remove(tmp+"cut_read2.fastq")

    subprocess.Popen(samtools+" view -f4 "+tmp+"genome/all.bam > "+tmp+"genome_unmapped.sam",shell=True).wait()
    umsam2fq(tmp+"genome_unmapped.sam",tmp+"genome_unmapped.fq")

    chrom=SeqIO.to_dict(SeqIO.parse(refgenome, "fasta"))
    bam2SNV(bam_in_path=tmp+"genome/all.bam", fa_in_path=refgenome, zz_out_path=tmp+"genome_all.zz",bed_out_path=tmp+"genome_all.snv")
    subprocess.Popen(bedtools+" sort -i "+tmp+"genome_all.snv"+" | uniq >"+tmp+"genome_all.snv.dedup.snv",shell=True).wait()                
    getBedAD(bed_in_path=tmp+"genome_all.snv.dedup.snv", allsnv_in_path=tmp+"genome_all.snv", bed_out_path=tmp+"genome_all_ad.snv")
    subprocess.Popen(bedtools+" sort -i "+tmp+"genome_all_ad.snv"+" | uniq >"+tmp+"genome_all_ad.snv.dedup",shell=True).wait()                

    ##----------------------------------------------------------------------------------------------------------

    ##transcript
    if os.path.exists(refgenome+".trans.fa"):
        print("\tStart mapping transcript ...")
        if not os.path.isdir(tmp+"transcript/"):
            os.mkdir(tmp+"transcript/")
        TAG="transcript"
        fq2bam(tmp+"transcript/",False,TAG,tmp+"/genome_unmapped.fq",read0,refgenome+".trans.fa",bwa,samtools,mapcpu)
        subprocess.Popen(samtools+" index "+tmp+"transcript/all.bam",shell=True).wait()
        chrom=SeqIO.to_dict(SeqIO.parse(refgenome+".trans.fa", "fasta"))
        bam2SNV(bam_in_path=tmp+"transcript/all.bam", fa_in_path=refgenome+".trans.fa", zz_out_path=tmp+"transcript_all.zz",bed_out_path=tmp+"transcript_all.snv")
        subprocess.Popen(bedtools+" sort -i "+tmp+"transcript_all.snv"+" | uniq >"+tmp+"transcript_all.snv.dedup.snv",shell=True).wait()                
        subprocess.Popen(samtools+" view -f4 "+tmp+"transcript/all.bam > "+tmp+"transcript_unmapped.sam",shell=True).wait()
        umsam2fq(tmp+"transcript_unmapped.sam",tmp+"/regular_unmapped.fq")    
        maskfq(tmp+"/regular_unmapped.fq","A","G")
        maskfq(tmp+"/regular_unmapped.fq","T","C")

        getBedAD(bed_in_path=tmp+"transcript_all.snv", allsnv_in_path=tmp+"transcript_all.snv", bed_out_path=tmp+"transcript_all_ad.snv")
        subprocess.Popen(bedtools+" sort -i "+tmp+"transcript_all_ad.snv"+" | uniq >"+tmp+"transcript_all_ad.snv.dedup",shell=True).wait()                
        tsnv2gsnv(refgenome+".trans.fa.loc", tmp+"transcript_all_ad.snv.dedup", tmp+"transcript_all_ad.snv.dedup.genome.snv")
        subprocess.Popen(bedtools+" sort -i "+tmp+"transcript_all_ad.snv.dedup.genome.snv"+" | uniq >"+tmp+"transcript_all_ad.snv.genome.snv.dedup.snv",shell=True).wait()                

    else:    
        subprocess.Popen("cp " + tmp+"genome_unmapped.fq" + " " + tmp+"/regular_unmapped.fq",shell=True).wait()    
        maskfq(tmp+"/regular_unmapped.fq","A","G")
        maskfq(tmp+"/regular_unmapped.fq","T","C")

    ##----------------------------------------------------------------------------------------------------------

    #genome_maskAG  
    if os.path.exists(refgenome+".mskAG.fa"):
        print("\tStart re-mapping genome ...")
        if not os.path.isdir(tmp+"genome_mskAG/"):
            os.mkdir(tmp+"genome_mskAG/")
        TAG="genome_mskAG"
        fq2sam(TAG,False,tmp+"/regular_unmapped_A_to_G.fq",read0,tmp,refgenome+".mskAG.fa",bwa,samtools,mapcpu)   
        subprocess.Popen(samtools+" index "+tmp+"genome_mskAG/all.bam",shell=True).wait()
        subprocess.Popen(samtools+" view -f4 "+tmp+"genome_mskAG/all.bam > "+tmp+"genome_mskAG_unmapped.sam",shell=True).wait()
        umsam2fq(tmp+"genome_mskAG_unmapped.sam",tmp+"genome_mskAG_unmapped.fq")    
        recover_sam(tmp+"genome_mskAG_all.sam",tmp+"genome_mskAG_all.sam.rcv", 20, 10, 0)
        sam2zz(tmp+"genome_mskAG_all.sam.rcv",refgenome,tmp+"genome_mskAG_all.zz")
        zz2snv(zz_in_dir=tmp+"genome_mskAG_all.zz",bed_out_dir=tmp+"genome_mskAG_all.snv",baseq_cutoff_dir=tmp+"baseq.cutoff")
        subprocess.Popen(bedtools+" sort -i "+tmp+"genome_mskAG_all.snv"+" | uniq >"+tmp+"genome_mskAG_all.snv.dedup.snv",shell=True).wait()                
        getBedAD(bed_in_path=tmp+"genome_mskAG_all.snv.dedup.snv", allsnv_in_path=tmp+"genome_mskAG_all.snv", bed_out_path=tmp+"genome_mskAG_all_ad.snv")
        subprocess.Popen(bedtools+" sort -i "+tmp+"genome_mskAG_all_ad.snv"+" | uniq >"+tmp+"genome_mskAG_all_ad.snv.dedup",shell=True).wait()                

    ##----------------------------------------------------------------------------------------------------------

    #genome_maskTC      
    if os.path.exists(refgenome+".mskTC.fa"):
        if not os.path.isdir(tmp+"genome_mskTC/"):
            os.mkdir(tmp+"genome_mskTC/")
        TAG="genome_mskTC"
        fq2sam(TAG,False,tmp+"/regular_unmapped_T_to_C.fq",read0,tmp,refgenome+".mskTC.fa",bwa,samtools,mapcpu)   
        subprocess.Popen(samtools+" index "+tmp+"genome_mskTC/all.bam",shell=True).wait()
        subprocess.Popen(samtools+" view -f4 "+tmp+"genome_mskTC/all.bam > "+tmp+"genome_mskTC_unmapped.sam",shell=True).wait()
        umsam2fq(tmp+"genome_mskTC_unmapped.sam",tmp+"genome_mskTC_unmapped.fq")    
        recover_sam(tmp+"genome_mskTC_all.sam",tmp+"genome_mskTC_all.sam.rcv", 20, 10, 0)
        sam2zz(tmp+"genome_mskTC_all.sam.rcv",refgenome,tmp+"genome_mskTC_all.zz")
        zz2snv(zz_in_dir=tmp+"genome_mskTC_all.zz",bed_out_dir=tmp+"genome_mskTC_all.snv",baseq_cutoff_dir=tmp+"baseq.cutoff")
        subprocess.Popen(bedtools+" sort -i "+tmp+"genome_mskTC_all.snv"+" | uniq >"+tmp+"genome_mskTC_all.snv.dedup.snv",shell=True).wait()                
        getBedAD(bed_in_path=tmp+"genome_mskTC_all.snv.dedup.snv", allsnv_in_path=tmp+"genome_mskTC_all.snv", bed_out_path=tmp+"genome_mskTC_all_ad.snv")
        subprocess.Popen(bedtools+" sort -i "+tmp+"genome_mskTC_all_ad.snv"+" | uniq >"+tmp+"genome_mskTC_all_ad.snv.dedup",shell=True).wait()                

    ##----------------------------------------------------------------------------------------------------------

    #transcript_maskAG   
    if os.path.exists(refgenome+".trans.fa.mskAG.fa"):
        print("\tStart re-mapping transcript ...")
        if os.path.exists(refgenome):
            if not os.path.isdir(tmp+"transcript_mskAG/"):
                os.mkdir(tmp+"transcript_mskAG/")
            TAG="transcript_mskAG"
            fq2sam(TAG,False,tmp+"/genome_mskAG_unmapped.fq",read0,tmp,refgenome+".trans.fa.mskAG.fa",bwa,samtools,mapcpu)   
            subprocess.Popen(samtools+" index "+tmp+"transcript_mskAG/all.bam",shell=True).wait()
            recover_sam(tmp+"transcript_mskAG_all.sam",tmp+"transcript_mskAG_all.sam.rcv", 20, 10, 0)
            sam2zz(tmp+"transcript_mskAG_all.sam.rcv",refgenome+".trans.fa",tmp+"transcript_mskAG_all.zz")
            zz2snv(zz_in_dir=tmp+"transcript_mskAG_all.zz",bed_out_dir=tmp+"transcript_mskAG_all.snv",baseq_cutoff_dir=tmp+"baseq.cutoff")
            getBedAD(bed_in_path=tmp+"transcript_mskAG_all.snv", allsnv_in_path=tmp+"transcript_mskAG_all.snv", bed_out_path=tmp+"transcript_mskAG_all_ad.snv")
            subprocess.Popen(bedtools+" sort -i "+tmp+"transcript_mskAG_all_ad.snv"+" | uniq >"+tmp+"transcript_mskAG_all_ad.snv.dedup",shell=True).wait()                
            tsnv2gsnv(refgenome+".trans.fa.loc", tmp+"transcript_mskAG_all_ad.snv.dedup", tmp+"transcript_mskAG_all_ad.snv.dedup.genome.snv")
            subprocess.Popen(bedtools+" sort -i "+tmp+"transcript_mskAG_all_ad.snv.dedup.genome.snv"+" | uniq >"+tmp+"transcript_mskAG_all_ad.snv.genome.snv.dedup.snv",shell=True).wait()                

    ##----------------------------------------------------------------------------------------------------------

    #transcript_maskTC
    if os.path.exists(refgenome+".trans.fa.mskTC.fa"):
        if os.path.exists(refgenome):
            if not os.path.isdir(tmp+"transcript_mskTC/"):
                os.mkdir(tmp+"transcript_mskTC/")
            TAG="transcript_mskTC"
            fq2sam(TAG,False,tmp+"/genome_mskTC_unmapped.fq",read0,tmp,refgenome+".trans.fa.mskTC.fa",bwa,samtools,mapcpu)   
            subprocess.Popen(samtools+" index "+tmp+"transcript_mskTC/all.bam",shell=True).wait()
            recover_sam(tmp+"transcript_mskTC_all.sam",tmp+"transcript_mskTC_all.sam.rcv", 20, 10, 0)
            sam2zz(tmp+"transcript_mskTC_all.sam.rcv",refgenome+".trans.fa",tmp+"transcript_mskTC_all.zz")
            zz2snv(zz_in_dir=tmp+"transcript_mskTC_all.zz",bed_out_dir=tmp+"transcript_mskTC_all.snv",baseq_cutoff_dir=tmp+"baseq.cutoff")
            getBedAD(bed_in_path=tmp+"transcript_mskTC_all.snv", allsnv_in_path=tmp+"transcript_mskTC_all.snv", bed_out_path=tmp+"transcript_mskTC_all_ad.snv")
            subprocess.Popen(bedtools+" sort -i "+tmp+"transcript_mskTC_all_ad.snv"+" | uniq >"+tmp+"transcript_mskTC_all_ad.snv.dedup",shell=True).wait()                
            tsnv2gsnv(refgenome+".trans.fa.loc", tmp+"transcript_mskTC_all_ad.snv.dedup", tmp+"transcript_mskTC_all_ad.snv.dedup.genome.snv")
            subprocess.Popen(bedtools+" sort -i "+tmp+"transcript_mskTC_all_ad.snv.dedup.genome.snv"+" | uniq >"+tmp+"transcript_mskTC_all_ad.snv.genome.snv.dedup.snv",shell=True).wait()                

    time_end = int(time.perf_counter())
    run_time = str(time_end-time_start)
    print("SNV calling running time is:" + run_time + " seconds" + "\n")        

    ##----------------------------------------------------------------------------------------------------------

    if os.path.exists(refgenome+".trans.fa"):
        subprocess.Popen("cat "+tmp+"/transcript_all_ad.snv.genome.snv.dedup.snv "+tmp+"/genome_all_ad.snv.dedup | uniq | bedtools sort -i > "+tmp+"/regular.snv",shell=True).wait()
        subprocess.Popen("cat "+tmp+"/transcript_mskTC_all_ad.snv.genome.snv.dedup.snv "+tmp+"/genome_mskTC_all_ad.snv.dedup | uniq | bedtools sort -i > "+tmp+"/hyper_mskTC.snv",shell=True).wait()
        subprocess.Popen("cat "+tmp+"/transcript_mskAG_all_ad.snv.genome.snv.dedup.snv "+tmp+"/genome_mskAG_all_ad.snv.dedup | uniq | bedtools sort -i > "+tmp+"/hyper_mskAG.snv",shell=True).wait()    
        subprocess.Popen("mv "+tmp+"/transcript_mskAG_* "+tmp+"/transcript_mskAG/",shell=True).wait()    
        subprocess.Popen("mv "+tmp+"/transcript_mskTC_* "+tmp+"/transcript_mskTC/",shell=True).wait()    
        subprocess.Popen("mv "+tmp+"/regular_unmap* "+tmp+"/transcript/",shell=True).wait()        
        subprocess.Popen("mv "+tmp+"/transcript_*.* "+tmp+"/transcript/",shell=True).wait()

    else:
        subprocess.Popen("cp "+tmp+"/genome_all_ad.snv.dedup "+tmp+"/regular.snv",shell=True).wait()
        subprocess.Popen("cp "+tmp+"/genome_mskTC_all_ad.snv.dedup "+tmp+"/hyper_mskTC.snv",shell=True).wait()
        subprocess.Popen("cp "+tmp+"/genome_mskAG_all_ad.snv.dedup "+tmp+"/hyper_mskAG.snv",shell=True).wait()
        subprocess.Popen("mv "+tmp+"/genome_mskAG_* "+tmp+"/genome_mskAG/",shell=True).wait()    
        subprocess.Popen("mv "+tmp+"/genome_mskTC_* "+tmp+"/genome_mskTC/",shell=True).wait()
        subprocess.Popen("mv "+tmp+"/genome_*.* "+tmp+"/genome/",shell=True).wait()
##----------------------------------------------------------------------------------------------------------



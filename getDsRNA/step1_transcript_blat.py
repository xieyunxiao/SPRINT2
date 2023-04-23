##BLAT all reference
##pip install biopython
##pip install multiprocessing
##pip install subprocess
##pip install argparse


from Bio import SeqIO 
from multiprocessing import Pool
import subprocess 
import os,time
import shutil
import re
import argparse

##----------------------------------------------------------------------------------------------------------
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
    
def oneChrRef(CHR="chr1"):
    this_chr=str(CHR)
    this_ref_path=ref_out_dir+"/"+this_chr+".fa"
    this_seq=">"+this_chr+"\n"+str(chrom[this_chr].upper())+"\n"
    fo=open(this_ref_path,'w')
    fo.write(this_seq)
    fo.close()

##----------------------------------------------------------------------------------------------------------
##find all Query sequences in all hg38 reference file with multiprocessing pool
##running time is:334 seconds (5.6 minutes)

def OneChrQueryFinder(CHR="chr1",query_add=20,blat_out_dir=0):
    this_chr=str(CHR)
    this_chr_seq=chrom[this_chr]
    out_path=blat_out_dir+"/"    
    chr_len=len(this_chr_seq)
    
    ##the first not N site
    noN_num_start=0
    for i in range(chr_len):
        this_site=this_chr_seq[i].upper()
        if this_site != "N":
            noN_num_start=i
            break
        if i==(chr_len-1):
            noN_num_start=(chr_len-1)
    ##the last not N site        
    noN_num_end=(chr_len-1)
    if noN_num_start < (chr_len-1):        
        for i in range((chr_len-1),-1,-1):
            this_site=this_chr_seq[i].upper()
            if this_site != "N":
                noN_num_end=i
                break    
    
        #screening        
        Query_seqs=[]    
        for this_point in range(noN_num_start+query_add,noN_num_end-query_add):
            query_add=int(query_add)
            query_start=int(max(noN_num_start,int(this_point-query_add)))
            query_end=int(min(noN_num_end,int(this_point + query_add)))
            query_seq=chrom[CHR][query_start:query_end].upper()
            query_seq=str(query_seq)
            Query_seqs.append(query_seq)
        Query_seqs=list(set(Query_seqs))       

        #writing in fasta file
        fa_out_path=out_path+CHR+"_Query.fa"
        fo=open(fa_out_path,'w') 
        for seqid in range(len(Query_seqs)):
            seq=Query_seqs[seqid]
            line = ">"+CHR+",seq"+str(seqid)+"\n"+str(seq)+"\n"
            fo.write(line)            
        fo.close()

##----------------------------------------------------------------------------------------------------------
##BLAT by chromosome
##running time is : ~ 3 days
#ARG=[reference_path,query_path,blat_out_path,blat_path,this_chr_id]
def blatWorker(ARG):
    ref_in_path=ARG[0]
    query_in_path=ARG[1]
    map_out_path=ARG[2]+"o.blast9"+"_"+str(ARG[4])
    blat_path=ARG[3]
    this_step =subprocess.Popen(' '.join([blat_path,'-out=blast9','-minScore=25','-fastMap','-minIdentity=70', '-noHead', ref_in_path, query_in_path, map_out_path+'.tmpOut', '>', map_out_path+'.tmpInfo' ]),shell=True)
    this_step.wait()
    fmap=open(map_out_path,'w')
    ftmp=open(map_out_path+'.tmpOut')
    for l_tmp in ftmp:
        if l_tmp[0]!='#':
            fmap.write(l_tmp)
    ftmp.close()
    fmap.close()


##----------------------------------------------------------------------------------------------------------
def main(args):
    parser = argparse.ArgumentParser(description='python3 step1_transcript_blat.py')
    parser.add_argument('-o', '--OUT_DIR',default = str(os.getcwd()),  help='path to Output Directory (default:Working Directory)')
    parser.add_argument('-r', '--Reference_file',default = str(os.getcwd())+"/reference.fa", help='path to Genome Reference File (FASTA format) (default:~/reference.fa)')
    parser.add_argument('-t','--Transcript_file',default = str(os.getcwd())+"/transcript.bed", help='path to Transctript Annotation File (BED format) (default:~/transcript.bed)')
    parser.add_argument('-b', '--Blat_file',default = str(os.getcwd())+"/blat",  help='path to blat (default:~/blat)')
    parser.add_argument('-p', '--cpu',default=1,  help='CPU number (default=1)')

    try:
        args = parser.parse_args(args.split())

    except argparse.ArgumentError as e:
        print(e)
        exit(1)
        
    global OUT_DIR_BASE,OUT_DIR,blat_path,trans_path,ref_in_path,CPU,SITE_ADD_AROUND,trans_add,add2
    OUT_DIR_BASE=str(args.OUT_DIR)+"/"
    blat_path = str(args.Blat_file)
    trans_path=str(args.Transcript_file)
    ref_in_path = str(args.Reference_file)
    CPU=int(args.cpu)

    SITE_ADD_AROUND=40
    trans_add=2
    ##trans_add:...kb

    SITE_ADD_AROUND=int(SITE_ADD_AROUND)
    trans_add=int(trans_add)*1000
    add2=2*SITE_ADD_AROUND

    ##----------------------------------------------------------------------------------------------------------
    description()
    print("\n\n\nStart step1 BLAT...")

    OUT_DIR = OUT_DIR_BASE+"/"
    if not os.path.isdir(OUT_DIR):
        os.mkdir(OUT_DIR)
    ##----------------------------------------------------------------------------------------------------------
    ##load reference genome fasta file and compute time
    ##running time is:40 seconds
    print("1.Loading reference genome file...")
    global chrom
    chrom=SeqIO.to_dict(SeqIO.parse(ref_in_path, "fasta"))
    #print(list(chrom.keys())[0:50])
    #print("\t\tReference genome file is OK!")
    global all_chroms,all_chroms_new
    ##print number of bases in each chrom
    all_chroms=list(chrom.keys())###640 sequences

    ##----------------------------------------------------------------------------------------------------------

    ##----------------------------------------------------------------------------------------------------------
    ##split each chrom seqence
    all_chroms_new={}
    fi=open(trans_path)
    while 1 :
        lines=fi.readlines(10000)
        if not lines:
            break
        for line in lines:
            this_line=line.strip("\n")
            this_info=this_line.split("\t")
            if True:
                this_chr=str(this_info[0])
                this_name=this_info[5]+"-"+str(this_info[1])+"-"+str(this_info[2])
                this_start=max(int(this_info[1])-trans_add,0)
                this_end=min(int(this_info[2])+trans_add,len(chrom[this_chr].seq))
                this_seq=this_seq=chrom[this_chr].seq[this_start-1:this_end+1].upper()
                all_chroms_new.update({this_name:this_seq})
    fi.close()
    all_chroms=list(all_chroms_new.keys())
    chrom=all_chroms_new

    ##----------------------------------------------------------------------------------------------------------

    ##----------------------------------------------------------------------------------------------------------
    ##split reference sequence by chromosome
    ##running time is:4 seconds
    global ref_out_dir
    ref_out_dir=OUT_DIR+"/reference"
    if not os.path.isdir(ref_out_dir):
        os.mkdir(ref_out_dir)

    pool_size = CPU
    pool1 = Pool(pool_size)  
    for this_chr in all_chroms:
        chrr=str(this_chr)
        pool1.apply_async(oneChrRef, (chrr,)) 
    pool1.close() 
    pool1.join()  

    ##----------------------------------------------------------------------------------------------------------

    ##----------------------------------------------------------------------------------------------------------
    ##find all Query sequences in all hg38 reference file with multiprocessing pool
    ##running time is:334 seconds (5.6 minutes)
    ## all sequence query finder 
    print("2.Findind query sequences....")
    global fa_out_dir
    fa_out_dir=OUT_DIR+"/query"
    if not os.path.isdir(fa_out_dir):
        os.mkdir(fa_out_dir)

    pool_size = int(CPU)
    pool1 = Pool(pool_size)  

    for this_chr in all_chroms:
        chr_len=len(str(chrom[this_chr]))
        pool1.apply_async(OneChrQueryFinder, (this_chr, SITE_ADD_AROUND, fa_out_dir))

    pool1.close() 
    pool1.join()

    ##----------------------------------------------------------------------------------------------------------

    ##----------------------------------------------------------------------------------------------------------
    ##BLAT by chromosome
    ##running time is : ~ 3 days
    ## all BLAT
    global blat_out_dir
    blat_out_dir=OUT_DIR+"/BLAT"
    CPU=CPU
    if not os.path.isdir(blat_out_dir):
        os.mkdir(blat_out_dir)
    print("3.Running BLAT....")     
    time_start = int(time.perf_counter())

    querys=os.listdir(OUT_DIR+"/query")
    pool_size = CPU
    pool1 = Pool(pool_size)  
    for qu in querys:
        this_chr=qu.split("_id_")[0]
        this_chr_id=qu.split("_Query")[0]
        reference_path=OUT_DIR+"/reference/"+this_chr_id+".fa"
        query_path=OUT_DIR+"/query/"+qu
        blat_out_path=blat_out_dir+"/"
        params=[reference_path,query_path,blat_out_path,blat_path,this_chr_id]
        pool1.apply_async(blatWorker, (params,)) 
    pool1.close()  
    pool1.join()  

    time_end = int(time.perf_counter())
    run_time = str(time_end-time_start)
    print("\t\tprocess BLAT running time is:" + run_time + " seconds" + "\n")     

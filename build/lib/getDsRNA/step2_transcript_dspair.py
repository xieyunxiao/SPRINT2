import os,time,shutil
from multiprocessing import Pool
import argparse


##----------------------------------------------------------------------------------------------------------

#ARG=[blat_in_file, ds_pair_out_file, CHR, start_pos, SITE_ADD_AROUND]    
def ds_pair(ARG):
    #o.blast9_SNX10-ENSG00000086300.16_9
    this_in_path=ARG[0]
    this_out_path=ARG[1]
    this_chr=ARG[2]
    this_start_pos=ARG[3]
    SITE_ADD_AROUND=ARG[4]     
    fi = open(this_in_path) 
    #5S_rRNA-ENSG00000201285.1,seq0  5S_rRNA-ENSG00000201285.1       100.00  80      0       0       1       80      166     245     7.7e-38 155.0
    fo=open(this_out_path,'w')
    this_write="seq1_chr\tseq1_start\tseq1_end\tseq2_chr\tseq2_start\tseq2_end\tquery_id\tseq1_identity\tseq1_score\tseq2_identity\tseq2_score\tdirection\tdistance"
    fo.write(this_write+'\n')        

    i = -1
    identity=[]
    score=[]
    start=[]
    end=[]
    this_query1="chr0_id_0,seq0"
    this_query2="chr0_id_0,seq0"
    while 1:
        lines = fi.readlines(100000)
        if not lines:
            break
        for line in lines:
            i=i+1
            if not line:
                break
            this_line=line.strip("\n")
            this_seq=this_line.split("\t")[0]
            if True:
                this_query2=this_seq
                
                if this_query2 == this_query1:                    
                    this_query1=this_query2                    
                    identity.append(this_line.split("\t")[2])
                    score.append(this_line.split("\t")[11])
                    start.append(this_line.split("\t")[8])
                    end.append(this_line.split("\t")[9])

                if this_query2 != this_query1:
                    #print(score)
                    #print(identity)
                    diff=[int(end[j])-int(start[j]) for j in range(len(start))]
                    #print(diff)
                    direction=['reverse' if j < 0 else 'forward' for j in diff]
                    i_reverse=[j for j in range(len(direction)) if direction[j]=='reverse']
                    #print(i_reverse)
                    i_forward=[j for j in range(len(direction)) if direction[j]=='forward']
                    if len(i_reverse)>0 and len(i_forward)>0:
                        i1_100=[j for j in range(len(identity)) if identity[j]=="100.00"]
                        i2_100=[j for j in range(len(identity)) if diff[j]==((2*SITE_ADD_AROUND)-1)]
                        i_100 = list(set(i1_100).intersection(set(i2_100)))
                        
                        for j_100 in i_100:
                            seq1_start=int(start[j_100])+this_start_pos
                            seq1_end=int(end[j_100])+this_start_pos
                            seq1_identity=identity[j_100]
                            seq1_score=score[j_100]
                            for j_forward in i_forward:
                                if j_forward != j_100:                        
                                    seq2_start=int(start[j_forward])+this_start_pos
                                    seq2_end=int(end[j_forward])+this_start_pos
                                    seq2_identity=identity[j_forward]
                                    seq2_score=score[j_forward]
                                    this_direction='forward'
                                    this_distance=max(abs(seq1_start-seq2_start),abs(seq1_start-seq2_end),abs(seq1_end-seq2_start),abs(seq1_end-seq2_end))
                                    
                                    if True:
                                        this_write=str(this_chr)+"\t"+str(seq1_start)+"\t"+str(seq1_end)+"\t"+str(this_chr)+"\t"+str(seq2_start)+"\t"+str(seq2_end)+"\t"+str(this_query1)+"\t"+str(seq1_identity)+"\t"+str(seq1_score)+"\t"+str(seq2_identity)+"\t"+str(seq2_score)+"\t"+str(this_direction)+"\t"+str(this_distance)
                                        fo.write(this_write+'\n')
                            for j_reverse in i_reverse:
                                if j_reverse != j_100:                        
                                    seq2_start=int(start[j_reverse])+this_start_pos
                                    seq2_end=int(end[j_reverse])+this_start_pos
                                    seq2_identity=identity[j_reverse]
                                    seq2_score=score[j_reverse]                                
                                    this_direction='reverse'
                                    this_distance=max(abs(seq1_start-seq2_start),abs(seq1_start-seq2_end),abs(seq1_end-seq2_start),abs(seq1_end-seq2_end))
                                    
                                    if True:
                                        this_write=str(this_chr)+"\t"+str(seq1_start)+"\t"+str(seq1_end)+"\t"+str(this_chr)+"\t"+str(seq2_end)+"\t"+str(seq2_start)+"\t"+str(this_query1)+"\t"+str(seq1_identity)+"\t"+str(seq1_score)+"\t"+str(seq2_identity)+"\t"+str(seq2_score)+"\t"+str(this_direction)+"\t"+str(this_distance)
                                        fo.write(this_write+'\n')    

                    score=[]
                    identity=[]
                    start=[]
                    end=[]
                    identity.append(this_line.split("\t")[2])
                    score.append(this_line.split("\t")[11])
                    start.append(this_line.split("\t")[8])
                    end.append(this_line.split("\t")[9])  
                    this_query1=this_query2


    fi.close()
    fo.close()

##----------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------

#ARG=[bed_in,bed_out]
def rev_ds(ARG):   
    bedpe_in=ARG[0]
    bedpe_out=ARG[1]
    fi = open(bedpe_in)
    fo = open(bedpe_out,"w")
    line=fi.readline()
    while 1:
        lines=fi.readlines(100000000)
        if not lines:
            break
        for this_line in lines:
            if not this_line:
                break       
            if "reverse" in this_line:
                fo.write(this_line)
            
    fi.close()
    fo.close()

##----------------------------------------------------------------------------------------------------------

def main(args):  
    parser = argparse.ArgumentParser(description='python3 step2_transcript_dspair.py')
    parser.add_argument('-o', '--OUT_DIR',default = str(os.getcwd()),  help='path to Output Directory (default:Working Directory)')
    parser.add_argument('-t','--Transcript_file',default = str(os.getcwd())+"/transcript.bed", help='path to Transctript Annotation File (BED format) (default:~/transcript.bed)')
    parser.add_argument('-p', '--cpu',default=1,  help='CPU number (default=1)')
    try:
        args = parser.parse_args(args.split())

    except argparse.ArgumentError as e:
        print(e)
        exit(1)
    global OUT_DIR,trans_path,CPU,SITE_ADD_AROUND,blat,ds_path,rev_path
    OUT_DIR_BASE=str(args.OUT_DIR)+"/"
    trans_path=str(args.Transcript_file)
    CPU=int(args.cpu)

    SITE_ADD_AROUND=40
    blat=OUT_DIR_BASE+"BLAT/"
    ds_path=OUT_DIR_BASE+"ds_files/"
    rev_path=OUT_DIR_BASE+"ds_files_rev/"

    if not os.path.isdir(ds_path):
        os.mkdir(ds_path)
    if not os.path.isdir(rev_path):
        os.mkdir(rev_path)

    ##----------------------------------------------------------------------------------------------------------

    print("Start step2 get dspairs...")

    time_start = int(time.perf_counter())
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
            this_info=line.strip("\n")
            this_info=this_info.split("\t")
            if True:
                CHR=this_info[0]
                this_name=this_info[5]+"-"+str(this_info[1])+"-"+str(this_info[2])
                start_pos=max(int(this_info[1])-2001,0)
                blat_in_file=blat+"o.blast9_"+this_name
                ds_pair_out_file=ds_path+"ds_"+this_name+".bedpe"
                #ARG=[file_in, file_out, CHR, SITE_ADD_AROUND,seq_id1,seq_id2]
                ARG=[blat_in_file, ds_pair_out_file, CHR, start_pos, SITE_ADD_AROUND]
                pool1.apply_async(ds_pair, (ARG,)) 

    fi.close()

    pool1.close() 
    pool1.join()  

    ##----------------------------------------------------------------------------------------------------------

    ##----------------------------------------------------------------------------------------------------------

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
            this_info=line.strip("\n")
            this_info=this_info.split("\t")
            if True:
                this_name=this_info[5]+"-"+str(this_info[1])+"-"+str(this_info[2])
                #this_gene_id=this_info[6]
                start_pos=max(int(this_info[1])-2001,0)
                file_in=ds_path+"ds_"+this_name+".bedpe"
                file_out=rev_path+"ds_"+this_name+".bedpe"
                ARG=[file_in,file_out]
                pool1.apply_async(rev_ds, (ARG,)) 

    fi.close()

    pool1.close() 
    pool1.join()  

    time_end = int(time.perf_counter())
    run_time = str(time_end-time_start)
    print("\t\tprocess running time is:" + run_time + " seconds" + "\n")        



##----------------------------------------------------------------------------------------------------------








import subprocess,os,sys
from multiprocessing import Pool
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
    
##----------------------------------------------------------------------------------------------------------

def transcript_assembler(fref_in_dir=0, fgtf_in_dir=0, ftrans_out_dir=0):
    if fref_in_dir==0 or fgtf_in_dir==0:
        print('fref_in_dir\tfgtf_in_dir\tftrans_out_dir')
        return 0
    fref=open(fref_in_dir)
    chrom={}
    chrr=''
    line=fref.read()
    line=line.split('>')
    for seq in line:
        if ' ' in seq:
            chrr=seq[0:seq.find(' ')]
        else:
            chrr=seq[0:seq.find('\n')]
        chrom[chrr]=seq[seq.find('\n'):].replace('\n','').upper()
    fref.close()

    fgtf=open(fgtf_in_dir)
    transcript={}
    trc=[]

    for line in fgtf:
        if '#' != line[0]:
            seq=line.rstrip().split('\t')
            if seq[0] in chrom and seq[2]=='exon' and 'transcript_id ' in seq[8] and int(seq[4]) <= len(chrom[seq[0]]):
                chrr=seq[0]
                begin=int(seq[3])
                end=int(seq[4])
                transcript_id = seq[8].split('transcript_id ')[1].split(';')[0].replace('"','')
                strand = seq[6]
                if transcript_id not in transcript:
                    trc.append(transcript_id)
                    transcript[transcript_id] = [chrr,strand,[begin,end]]
                else:
                    transcript[transcript_id].append([begin,end])
    
    
    fgtf.close()

    for transcript_id in transcript:
        tmp=transcript[transcript_id][2:]
        tmp.sort()
        transcript[transcript_id][2:]=tmp
        #if transcript[transcript_id][1]=='-':
        #    transcript[transcript_id][2:]=transcript[transcript_id][2:][::-1]

    ftrans=open(ftrans_out_dir,'w')
    ftransloc=open(ftrans_out_dir+'.loc','w')


    for transcript_id in trc:
        trans_name='>'+transcript_id+'_|_'+transcript[transcript_id][0]+'_|_'+transcript[transcript_id][1]+'\n'
        chrr=transcript[transcript_id][0]
        strand=transcript[transcript_id][1]
        loc=transcript[transcript_id][2:]
        trans_seq=''
        ftransloc.write(trans_name)
        for one in loc:
            ftransloc.write(str(one[0])+','+str(one[1])+';')
            trans_seq += chrom[chrr][one[0]-1:one[1]]
        ftransloc.write('\n')
        ftrans.write(trans_name)
        i=0
        while i < len(trans_seq):
            ftrans.write(trans_seq[i])
            i += 1
            if i%50==0:
                ftrans.write('\n')
        if i%50!=0:
            ftrans.write('\n')
##----------------------------------------------------------------------------------------------------------

def maskAwithG(fa_in_dir=0,fa_out_dir=0):
    if fa_in_dir==0 or fa_out_dir==0:
        print('fa_in_dir\tfa_our_dir')
        return 0
    fi=open(fa_in_dir)
    fo=open(fa_out_dir,'w')

    for line in fi:
        if line[0]=='>':
            fo.write(line)
        else:
            fo.write(line.replace('A','G').replace('a','g'))

##----------------------------------------------------------------------------------------------------------

def maskTwithC(fa_in_dir=0,fa_out_dir=0):
    if fa_in_dir==0 or fa_out_dir==0:
        print('fa_in_dir\tfa_our_dir')
        return 0
    fi=open(fa_in_dir)
    fo=open(fa_out_dir,'w')

    for line in fi:
        if line[0]=='>':
            fo.write(line)
        else:
            fo.write(line.replace('T','C').replace('t','c'))

##----------------------------------------------------------------------------------------------------------

def bwa_worker(fa_in_dir=0):
    step=subprocess.Popen(bwa+' index -a bwtsw '+fa_in_dir,shell=True)
    step.wait()

##----------------------------------------------------------------------------------------------------------

def main(): 
    parser = argparse.ArgumentParser(description='Mask A with G on reference FASTA file.')
    parser.add_argument('-r', '--reference',default = str(os.getcwd())+"/reference.fa",  help='path to reference FASTA file (default:~/reference.fa)')
    parser.add_argument('-o', '--output',default = str(os.getcwd()),  help='path to output directory (default:~/)')
    parser.add_argument('-t', '--gtf',  help='path to transcript annotation GTF file #Optional')
    parser.add_argument('-b', '--bwa',default = str(os.getcwd())+"/bwa",  help='path to bwa (default:~/bwa)')
    parser.add_argument('-p', '--cpu',default=1,  help='CPU number (default=1)')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        print(e)
        exit(1)

    global bwa,refgenome,CPU, gtf_file
    refgenome=str(args.reference)
    base_path=str(args.output)+"/"
    bwa=str(args.bwa)
    CPU=int(args.cpu)

    if args.gtf:
        gtf_file=str(args.gtf)
    else:
        gtf_file=False

    global read2
    read2=''
    paired_end=True
    description()
    pool_size = CPU
    pool1 = Pool(pool_size)

    print('Masking A with G in reference genome...')
    maskAwithG(refgenome,refgenome+'.mskAG.fa')
    print('Masking T with C in reference genome...')
    maskTwithC(refgenome,refgenome+'.mskTC.fa')
    
    if gtf_file != False:
        print('Assembling transcripts...')
        transcript_assembler(refgenome,gtf_file,refgenome+'.trans.fa')
        print('Masking A with G in transcripts...')
        maskAwithG(refgenome+'.trans.fa',refgenome+'.trans.fa.mskAG.fa')
        print('Masking T with C in transcripts...')
        maskTwithC(refgenome+'.trans.fa',refgenome+'.trans.fa.mskTC.fa')
    
    
    print('Building BWA index...')
    pool1.apply_async(bwa_worker, (refgenome,)) 
    pool1.apply_async(bwa_worker, (refgenome+'.mskAG.fa',)) 
    pool1.apply_async(bwa_worker, (refgenome+'.mskTC.fa',)) 

    if gtf_file !=False:
        pool1.apply_async(bwa_worker, (refgenome+'.trans.fa',)) 
        pool1.apply_async(bwa_worker, (refgenome+'.trans.fa.mskAG.fa',)) 
        pool1.apply_async(bwa_worker, (refgenome+'.trans.fa.mskTC.fa',)) 
    pool1.close() 
    pool1.join()  
##----------------------------------------------------------------------------------------------------------

    
    


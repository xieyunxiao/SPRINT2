import argparse,os
import subprocess
from .step1_SNVcalling import *
from .step2_annotation_based import *
from .step3_dsRNA_based import *

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
    this_step1 = subprocess.Popen(' '.join([bedtools_path,'multicov','-bams',bam_in_path,'-bed',bed_in_path,' >',bed_out_path+"_hp"]),shell=True)
    this_step1.wait()
    this_step2 = subprocess.Popen('| sed -i "s/\b0\b/1/g" '+bed_out_path+"_hp"+' > '+bed_out_path,shell=True)
    this_step2.wait()
    this_step3 = subprocess.Popen('rm '+bed_out_path+"_hp",shell=True)
    this_step3.wait()
def run():
    parser = argparse.ArgumentParser(description='Get candidate double-stranded RNA')
    parser.add_argument('-s', '--SNV_PATH',default = str(os.getcwd())+"/1_SNV_calling/",  help='path to SNV directory (default:~/1_SNV_calling/)')
    parser.add_argument('-o', '--RES_PATH',default = str(os.getcwd())+"/2_RES_calling/",  help='path to output directory (default:~/2_RES_calling/)')
    parser.add_argument('-r1', '--read1',default = str(os.getcwd())+"/test_1.fastq", help='path to sample FASTQ_1 file (default:~/test_1.fastq)')
    parser.add_argument('-r2', '--read2', help='path to reference FASTQ_2 file #Optional')
    parser.add_argument('-r', '--Reference_file',default = str(os.getcwd())+"/reference.fa", help='path to reference FASTA file (default:~/reference.fa)')
    parser.add_argument('-rp', '--Repeat_file',  help='path to repeat BED file #Optional')
    parser.add_argument('-b', '--bwa',default = str(os.getcwd())+"/bwa",  help='path to BWA (default:~/bwa)')
    parser.add_argument('-B', '--bedtools',default = str(os.getcwd())+"/bedtools",  help='path to BEDTOOLS (default:~/bedtools)')
    parser.add_argument('-S', '--samtools',default = str(os.getcwd())+"/samtools",  help='path to SAMTOOLS (default:~/samtools)')
    parser.add_argument('-ds', '--Candidate_dsRNA',  help='path to candidate dsRNA directory #Optional')
    parser.add_argument('-p', '--cpu',default=1,  help='CPU number (default=1)')
    args = parser.parse_args()

    if args.read2:
        step1_SNVcalling.main("-o "+args.SNV_PATH+" -r1 "+args.read1+" -r2 "+args.read2+" -R "+args.Reference_file+" -b "+args.bwa+" -B "+args.bedtools+" -s "+args.samtools+" -p "+args.cpu)
    else:
        step1_SNVcalling.main("-o "+args.SNV_PATH+" -r1 "+args.read1+" -R "+args.Reference_file+" -b "+args.bwa+" -B "+args.bedtools+" -s "+args.samtools+" -p "+args.cpu)

    if args.Repeat_file:
        step2_annotation_based.main("-o "+args.RES_PATH+" -s "+args.SNV_PATH+" -rp "+args.Repeat_file)
    else:
        step2_annotation_based.main("-o "+args.RES_PATH+" -s "+args.SNV_PATH)

    if args.Candidate_dsRNA:
        step3_dsRNA_based.main("-o "+args.RES_PATH+" -s "+args.SNV_PATH+" -ds "+args.Candidate_dsRNA+" -b "+args.bedtools+" -p "+args.cpu)

    ##============================================================================================================
    outpath0=args.RES_PATH
    snv_path=args.SNV_PATH
    bedtools_path=args.bedtools
    if os.path.exists(outpath0+"ds_regular.res "):
        if os.path.exists(outpath0+"anno_regular.res "):
            subprocess.Popen("cat "+outpath0+"ds_regular.res "+outpath0+"anno_regular.res | "+bedtools_path+" sort -i | uniq > "+outpath0+"regular.res",shell=True).wait()
        else:
            subprocess.Popen("cat "+outpath0+"ds_regular.res | "+bedtools_path+" sort -i | uniq > "+outpath0+"regular.res",shell=True).wait()
    else:
        if os.path.exists(outpath0+"anno_regular.res "):
            subprocess.Popen("cat "+outpath0+"anno_regular.res | "+bedtools_path+" sort -i | uniq > "+outpath0+"regular.res",shell=True).wait()
        else:
            subprocess.Popen("touch "+outpath0+"ds_regular.res | "+bedtools_path+" sort -i | uniq > "+outpath0+"regular.res",shell=True).wait()

    ARG=[outpath0+"regular.res",outpath0+"SPRINT2_Regular_RES.stat"]
    stat_res(ARG)

    subprocess.Popen("cat "+outpath0+"regular.res "+outpath0+"hyper.res | "+bedtools_path+" sort -i | uniq > "+outpath0+"ALL.res",shell=True).wait()
    ARG=[outpath0+"ALL.res",outpath0+"SPRINT2_ALL_RES.stat"]
    stat_res(ARG)

    getBedDP(bed_in_path=outpath0+"ALL.res", bam_in_path=snv_path+"genome/all.bam", bed_out_path=outpath0+"ALL.res.dp", bedtools_path=bedtools_path)
    getBedDP(bed_in_path=outpath0+"regular.res", bam_in_path=snv_path+"genome/all.bam", bed_out_path=outpath0+"regular.res.dp", bedtools_path=bedtools_path)
    getBedDP(bed_in_path=outpath0+"hyper.res", bam_in_path=snv_path+"genome/all.bam", bed_out_path=outpath0+"hyper.res.dp", bedtools_path=bedtools_path)

    subprocess.Popen('echo "#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tDepth" | cat - '+outpath0 +'/ALL.res.dp   > '+outpath0+'/SPRINT2_identified_all.res',shell=True).wait()
    subprocess.Popen('echo "#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tDepth" | cat - '+outpath0 +'/regular.res.dp   > '+outpath0+'/SPRINT2_identified_regular.res',shell=True).wait()
    subprocess.Popen('echo "#Chrom\tStart(0base)\tEnd(1base)\tType\tSupporting_reads\tDepth" | cat - '+outpath0 +'/hyper.res.dp   > '+outpath0+'/SPRINT2_identified_hyper.res',shell=True).wait()

    if not os.path.isdir(outpath0+"RES"):
        os.mkdir(outpath0+"RES")
    if not os.path.isdir(outpath0+"anno_regular"):
        os.mkdir(outpath0+"anno_regular")
    subprocess.Popen("mv " + outpath0+"regular.res* " + outpath0+"anno_regular/",shell=True).wait() 
    subprocess.Popen("mv " + outpath0+"hyper.res* " + outpath0+"RES/",shell=True).wait() 
    subprocess.Popen("mv " + outpath0+"ALL.res* " + outpath0+"RES/",shell=True).wait() 
    if os.path.exists(outpath0+"anno_regular.res"):
        subprocess.Popen("mv " + outpath0+"anno*res* " + outpath0+"RES/",shell=True).wait() 
    if os.path.exists(outpath0+"ds_regular.res"):
        subprocess.Popen("mv " + outpath0+"ds*res* " + outpath0+"RES/",shell=True).wait() 


if __name__ == '__main__':
    run()


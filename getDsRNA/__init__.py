import argparse,os
import subprocess
from .step1_transcript_blat import *
from .step2_transcript_dspair import *
from .step3_bedtools_sorted_rmLC import *
from .step4_combine_by_chr import *
from .step5_merge_dsRNA import *
from .step6_dsRNA import *


def run():
    parser = argparse.ArgumentParser(description='Get candidate double-stranded RNA')
    parser.add_argument('-o', '--OUT_DIR',default = str(os.getcwd()),  help='path to output directory (default:~/)')
    parser.add_argument('-r', '--Reference_file', help='path to reference FASTA file (default:~/reference.fa)')
    parser.add_argument('-t', '--Transcript_file',  help='path to transcript annotation bed file (default:~/transcript.bed)')
    parser.add_argument('-b', '--Blat_file',default = str(os.getcwd())+"/blat",  help='path to blat (default:~/blat)')
    parser.add_argument('-B', '--Bedtools_path',default = str(os.getcwd())+"/bedtools",  help='path to bedtools (default:~/bedtools)')
    parser.add_argument('-lc', '--lc_path',default = str(os.getcwd())+"/Low_complex.txt",  help='path to output Low_complex.txt (default:~/Low_complex.txt)')
    parser.add_argument('-p', '--cpu',default=1,  help='CPU number (default=1)')
    args = parser.parse_args()

    step1_transcript_blat.main("-o "+args.OUT_DIR+" -r "+args.Reference_file+" -t "+args.Transcript_file+" -b "+args.Blat_file+" -p "+args.cpu)
    step2_transcript_dspair.main("-o "+args.OUT_DIR+" -t "+args.Transcript_file+" -p "+args.cpu)
    step3_bedtools_sorted_rmLC.main("-o "+args.OUT_DIR+" -B "+args.Bedtools_path+" -lc "+args.lc_path+" -p "+args.cpu)
    step4_combine_by_chr.main("-o "+args.OUT_DIR+" -B "+ args.Bedtools_path+" -t "+args.Transcript_file+" -p "+args.cpu)
    step5_merge_dsRNA.main("-o "+args.OUT_DIR+" -B "+args.Bedtools_path+" -p "+args.cpu)
    step6_dsRNA.main("-o "+args.OUT_DIR+" -p "+args.cpu)

if __name__ == '__main__':
    run()

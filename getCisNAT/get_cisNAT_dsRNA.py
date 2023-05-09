import subprocess 
bedtools_path="~/bin/bedtools"
gene_bed="~/gene.bed"
gene_string1_bed="~/gene_string1.bed"
gene_string2_bed="~/gene_string2.bed"
gene_string12_overlap_bed="~/gene_string12_overlap.bed"
gene_ds_cisNAT= "~/gene_string12_overlap_cisNAT.bedpe"

fi = open(gene_bed)
fo1=open(gene_string1_bed,"w")
fo2=open(gene_string2_bed,"w")
while 1:
    lines=fi.readlines(10000000)
    if not lines:
        break
    for line in lines:
        line=line.strip("\n")
        if not line:
            break
        if "+" in line:
            fo1.write(line+"\n")
        if "-\t" in line:
            fo2.write(line+"\n")
fo1.close()
fo2.close()
subprocess.Popen(bedtools_path+" intersect -a "+gene_string1_bed+" -b "+gene_string2_bed+" -wa -wb > "+gene_string12_overlap_bed,shell=True).wait()

fi=open(gene_string12_overlap_bed)
fo=open(gene_ds_cisNAT,"w")
while 1:
    lines=fi.readlines(1000000)
    if not lines:
        break
    for line in lines:
        if not line:
            break
        line=line.strip("\n")
        this_info=line.split("\t")

        this_chr=this_info[0]
        this_start1=int(this_info[1])
        this_end1=int(this_info[2])
        this_start2=int(this_info[9])
        this_end2=int(this_info[10])

        if this_start1 < this_end2 and this_start1 > this_start2:
            this_w=[this_chr,str(this_start1),str(this_end2)]
            this_w.extend(this_info)
            this_write="\t".join(this_w)
            fo.write(this_write+"\n")

fi.close()
fo.close()

gene_pair= "~/gene_string12_overlap_cisNAT.bedpe"
gene_ds_cisNAT= "~/gene_ds_cisNAT.bed"
fi=open(gene_pair)
fo=open(gene_ds_cisNAT,"w")
while 1:
    lines=fi.readlines(1000000)
    if not lines:
        break
    for line in lines:
        if not line:
            break
        line=line.strip("\n")
        this_info=line.split("\t")
        this_write="\t".join(this_info[0:3])
        fo.write(this_write+"\n")


fi.close()
fo.close()

gene_string1_bed="~/gene_string1.bed"
gene_string2_bed="~/gene_string2.bed"
gene_string12_overlap_bed="~/gene_string12_overlap.bed"
gene_ds_cisNAT= "~/gene_string12_overlap_cisNAT.bedpe"
gene_pair= "~/gene_string12_overlap_cisNAT.bedpe"
gene_ds_cisNAT= "~/gene_ds_cisNAT.bed"
subprocess.Popen("rm "+gene_string1_bed,shell=True).wait()
subprocess.Popen("rm "+gene_string2_bed,shell=True).wait()
subprocess.Popen("rm "+gene_string12_overlap_bed,shell=True).wait()

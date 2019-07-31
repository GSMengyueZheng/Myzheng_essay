#!/usr/bin/env python
# File Name: WGSporadic_Run_202_BAM2vcf.py
# Author: Wuy
# Created Time: Sat 03 Sep 2016 05:47:34 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import os

def conf2dic(config):
    from polangxin import ParseConf
    conf_dic     = ParseConf(config)
    return conf_dic

def id2bam(id):
    return "%s/HQData/Sample_%s/%s.sorted.rmdup.recal.bam" % (os.getcwd(), id, id)

#hc
def hc():
    cmd = """%s HaplotypeCaller \
             -R %s \
             -I %s \
             --genotyping-mode DISCOVERY \
             -stand-call-conf 30 \
             --min-base-quality-score 20 \
             -L %s \
             -O %s/raw/%s.hc.vcf \
             """.strip()
    return cmd

#def vt():
#    cmd = """vt decompose -s %s/raw/%s.hc.vcf | \
#        vt normalize -r %s - > %s/%s.hc.vcf"""
#    return cmd

def filtered():
    cmd = "WGSporadic_Script_200_filterVCF.py %s/raw/%s.hc.vcf \
            > %s/%s.hc.vcf 2> TMP/%s.hc.failed.vcf"
    return cmd


# pileup
def pileup():
    cmd = """python /GPFS02/zhengmy/WGS_Pip/WGS2019_hg19/Vcf2Bed.py \
             %s/raw/%s.hc.vcf > TMP/%s.hc.bed&&%s merge -i TMP/%s.hc.bed > \
             TMP/%s.hc.merge.bed&&%s mpileup \
            -A -d 100000 -B -q 0 -Q 20  \
            -l TMP/%s.hc.merge.bed \
            -f %s \
            %s \
            | gzip > %s/%s.pileup.gz"""
    return cmd

def annovarcmd():
    cmd1 = "perl %s/convert2annovar.pl\
            -format vcf4  %s/%s.hc.vcf > TMP/%s.hc.anoinput "
    cmd2 = "%s/table_annovar.pl \
            TMP/%s.hc.anoinput \
            %s \
            -buildver hg19 -protocol \
            refGene,snp138,snp138NonFlagged,cosmic70,clinvar_20160302,popfreq_all_20150413,ljb26_all,gnomad_genome,gnomad_exome,rmsk  \
            -operation g,f,f,f,f,f,f,f,f,r         -nastring .    \
            -remove -out TMP/%s.hc"
    cmd3 =  "mv TMP/%s.hc.hg19_multianno.txt %s/%s.hc.tsv "
    return " && ".join([cmd1, cmd2, cmd3])
def vep():
    cmd1 = "WGSporadic_Script_201_Vcf2Vep.py %s/%s.hc.vcf \
            > TMP/%s.hc.vepinput"
#    cmd2 = "source ~njsh/bashrc.vep"
    cmd3 = "%s \
            -i TMP/%s.hc.vepinput \
            -o TMP/%s.hc.vepoutput    \
            --verbose --no_progress  \
            --shift_hgvs 1 --force_overwrite   --everything  \
            --fasta %s \
            --refseq --dir %s \
            --offline  --buffer_size 250000 --species homo_sapiens"
    cmd4 = "WGSporadic_Script_202_VEParse.py TMP/%s.hc.vepoutput \
            > %s/%s.hc.vep.csv"
    return " && ".join([cmd1, cmd3, cmd4])

def oncotator():
    cmd = "%s \
            -v --input_format=VCF --output_format=TCGAMAF \
            --db-dir %s \
            %s/%s.hc.vcf \
            %s/%s.hc.maf hg19"
    return cmd

def main(config, caseid, outdir="Data_TMP"):
    dic      = conf2dic(config)
    ref      = dic['ref']
    tmp_dir  = dic['tmp_dir']
    gatk     = dic['gatk']
    target   = dic['target']
    java18   = dic['java']
    samtools = dic['samtools']
    vep83    = dic["vep83"]
    onco= dic["oncotator"]
    onco_db  = dic["oncotator_db"]

    annovar = dic["annovar"]
    annovar_db = dic["annovar_db"]
    vep83db = dic["vep83db"]
    ensembl_fa = dic["ensembl_fa"]
    bedtools = dic["bedtools"]

    print 'echo "== Start hc call =="'
    print 'date'
    print hc()%(gatk, ref,"%s"%(id2bam(caseid)),
           target, outdir, caseid)

    print 'echo "== Start filtered VCF =="'
    print 'date'
    print filtered()%(outdir, caseid, outdir, caseid, caseid)

    print 'echo "== Start call peliup =="'
    print 'date'
    print pileup()%(outdir,caseid,caseid,bedtools,caseid,caseid,
           samtools, caseid, ref, id2bam(caseid),
           outdir, caseid)

    print 'echo "== Start Annovar =="'
    print 'date'
    print annovarcmd()%(annovar, outdir, caseid, caseid,
            annovar, caseid, annovar_db, caseid,
            caseid, outdir, caseid)

    print 'echo "== Start VEP annotate =="'
    print 'date'
    print vep()%(outdir, caseid, caseid,
            vep83, caseid, caseid,ensembl_fa, vep83db,
            caseid, outdir, caseid)

    print 'echo "== Start Oncotator annotate =="'
    print 'date'
    print oncotator()%(onco, onco_db, outdir, caseid, outdir, caseid)
    print 'date'

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print "py config sampleid"
        exit(1)
    t = sys.argv
    main(t[1], t[2])


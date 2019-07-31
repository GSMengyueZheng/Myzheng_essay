#!/usr/bin/env python
# File Name: WGSporadic_Run_102_RawData2BAM.py
# Author:
# Created Time: Thu 01 Sep 2016 02:20:04 PM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import sys
import os


def conf2dic(config):
    from polangxin import ParseConf
    conf_dic     = ParseConf(config)
    return conf_dic

def prepare(key):
    os.system('mkdir -p TMP HQData OUT/CovDep OUT/InsertSize OUT/VCF run_scripts Infos_TMP')

def QC(config, key, rawdir = "RawData", hqdir = "HQData"):
    dic = conf2dic(config)
    trimmomatic = dic['trimmomatic']
    tmp_dir = dic["tmp_dir"]
    qc_fa = dic["qc_fa"]
    java18 = dic["java"]

    os.system("mkdir -p %s/Sample_%s" % (hqdir, key))
    r1 = "%s/%s_R1.fastq.gz"%(rawdir, key)
    r2 = "%s/%s_R2.fastq.gz"%(rawdir, key)
    h1 = "%s/Sample_%s/%s_R1.fastq.gz"%(hqdir, key, key)
    h2 = "%s/Sample_%s/%s_R2.fastq.gz"%(hqdir, key, key)

    qc_raw = "FastqInfos.py %s/%s_R1.fastq.gz >> Infos_TMP/RawInfos.csv && \
               FastqInfos.py %s/%s_R2.fastq.gz >> Infos_TMP/RawInfos.csv & \
               " % (rawdir, key, rawdir, key)
    cmd1 = """%s  -Djava.io.tmpdir=%s \
            -XX:ParallelGCThreads=8 \
            -jar %s PE \
            -threads  8 \
            %s \
            %s \
            %s \
            TMP/%s_R1.unpaired.fastq.gz \
            %s \
            TMP/%s_R2.unpaired.fastq.gz \
            ILLUMINACLIP:%s:2:20:10:1:true \
            LEADING:15 \
            TRAILING:15 \
            SLIDINGWINDOW:5:20 \
            AVGQUAL:20 \
            MINLEN:36 \
            """ % (java18, tmp_dir, trimmomatic, r1, r2, h1, key, h2, key, qc_fa)
    qc_hq = "FastqInfos.py %s/Sample_%s/%s_R1.fastq.gz >> Infos_TMP/CleanInfos.csv && \
             FastqInfos.py %s/Sample_%s/%s_R2.fastq.gz >> Infos_TMP/CleanInfos.csv & \
            " % (hqdir, key, key, hqdir, key, key)
    print "date"
    print 'echo "== Start QC =="'
    print qc_raw
    print cmd1
    print qc_hq

def Mapping(config, key, hqdir = "HQData"):
    dic = conf2dic(config)
    bwa = dic["bwa"]
    ref = dic["ref"]
    picard25 = dic["picard25"]
    java18 = dic["java"]
    tmp_dir = dic["tmp_dir"]
    cmd1 = """%s mem -M -t 8 \
            -R "@RG\\tID:%s\\tSM:%s\\tLB:%s\\tPL:illumina" \
            %s \
            %s/Sample_%s/%s_R1.fastq.gz \
            %s/Sample_%s/%s_R2.fastq.gz \
            | gzip > %s/Sample_%s/%s.sam.gz \
            """ % (bwa, key, key, key, ref, hqdir, key, key,
                    hqdir, key, key, hqdir, key, key)
    cmd2 = """/GPFS02/zhengmy/WGS_Pip/WGS2019_hg19/pigz-2.3.4/unpigz \
              -c -p 10 %s/Sample_%s/%s.sam.gz | /GPFS01/softwares/sambamba/sambamba-0.6.9-linux-static view -t 10 -f bam -l 0 -S /dev/stdin | /GPFS01/softwares/sambamba/sambamba-0.6.9-linux-static sort -t 10 -m 3G --tmpdir TMP -o %s/Sample_%s/%s.sorted.bam /dev/stdin \
           """ % (hqdir,key,key,hqdir,key,key) 
    print "date"
    print 'echo "== Start Mapping =="'
    print cmd1
    print "date"
    print 'echo "== Start Sorted BAM =="'
    print cmd2

def Rmdup2BQSR(config, key, hqdir = "HQData"):
    dic = conf2dic(config)
    ref = dic["ref"]
    picard25 = dic["picard25"]
    tmp_dir = dic["tmp_dir"]
    gatk = dic["gatk"]
    indel_gsdb = dic["indel_gsdb"]
    indel_1000g = dic["indel_1000g"]
    indel_mills = dic["indel_mills"]
    target = dic["target"]
    snp_gsdb = dic["snp_gsdb"]
    java18 = dic["java"]
    snp_dbsnp = dic["snp_dbsnp"]

    # rmdup
    cmd1 = """%s -Djava.io.tmpdir=%s \
            -jar %s MarkDuplicates \
            INPUT=%s/Sample_%s/%s.sorted.bam \
            OUTPUT=%s/Sample_%s/%s.sorted.rmdup.bam \
            METRICS_FILE=%s/Sample_%s/%s.sorted.rmdup.metrics \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            CREATE_INDEX=true \
            """ % (java18, tmp_dir, picard25, hqdir, key, key,
                    hqdir, key, key, hqdir, key, key)

    # BaseRecalibrator
    cmd2 = """%s BaseRecalibrator -I %s/Sample_%s/%s.sorted.rmdup.bam \
              -R %s \
              --known-sites %s \
              --known-sites %s \
              --known-sites %s \
              --known-sites %s \
              --known-sites %s \
              -O %s/Sample_%s/%s.sorted.rmdup.recal_data.table \
            """ % (gatk,hqdir,key,key,ref,indel_gsdb,indel_1000g,
                   indel_mills,snp_gsdb,snp_dbsnp,hqdir,key,key)
    
    # ApplyBQSR
    cmd3 = """%s ApplyBQSR -I %s/Sample_%s/%s.sorted.rmdup.bam \
            -bqsr %s/Sample_%s/%s.sorted.rmdup.recal_data.table \
            -O %s/Sample_%s/%s.sorted.rmdup.recal.bam
            """ % ( gatk, hqdir, key, key,
                    hqdir, key, key, hqdir, key, key)

    # SUMMARY
    cmd4 = """python /GPFS02/zhengmy/WGS_Pip/WGS2019_hg19/SUMMARY.WGS.hg19.py \
              -Id %s -iP Infos_TMP/ -bP %s/Sample_%s/ -oP %s/Sample_%s/ \
           """ % (key,hqdir,key,hqdir,key)
    
    print 'date'
    print 'echo "== Start MarkDuplicates =="'
    print cmd1

    print 'date'
    print 'echo "== Start BaseRecalibrator =="'
    print cmd2

    print 'date'
    print 'echo "== Start ApplyBQSR =="'
    print cmd3

    print 'date'
    print 'echo "== Start SUMMARY =="'
    print cmd4


def main(config, key):
    QC(config, key)
    Mapping(config, key)
    Rmdup2BQSR(config, key)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "%s config in.id" % os.path.basename(__file__)
        exit(1)
    main(sys.argv[1], sys.argv[2])


#!/usr/bin/env python
#-*-coding:utf-8
import os
import sys
import argparse

def fil2dic(fil):
  with open(fil) as fh:
    dicR1 = {}
    dicR2 = {}
    for i in fh:
        i = i.rstrip()
        if i.startswith("#SAMPLE"):
            next
        t = i.split(",")
        name = t[0]
        na = t[0].split("_")[0]
        Reads = t[1]
        Length = t[2]
        Bases = t[3]
        GC = t[4]
        Q20 = t[5]
        Q30 = t[6]
        key = ",".join([Length,GC,Q20,Q30,Reads,Bases])
        if "R1" in name:
            dicR1[na] = key
        elif "R2" in name:
            dicR2[na] = key
    return dicR1,dicR2

def sumdic(dicR1,dicR2):
    dic = {}
    for sample in dicR1:
        if sample in dicR2:
            A = [float(x) for x in dicR1[sample].split(",")]
            B = [float(x) for x in dicR2[sample].split(",")]
            Length = str(round((A[0]+B[0])/2,0))
            GC = str(round((A[1]+B[1])/2,2))
            Q20 = str(round((A[2]+B[2])/2,2))
            Q30 = str(round((A[3]+B[3])/2,2))
            Reads = str(round((A[4]+B[4])/2,0))
            Bases = str(round((A[5]+B[5])/2,1))
            dic[sample] = ",".join([Length,GC,Q20,Q30,Reads,Bases])
    return dic

def pickmapped(sample,BamPath):
    with open(BamPath+sample+".flagstat") as fh:
        n = 0
        for lin in fh:
            if n == 4:
                t = lin.split("(")[1]
                m = t.split(":")[0]
                n += 1
                MappedRate = str(m).replace("%","")
            else:
                n += 1
                MappedRate = "NoFound"
    return MappedRate

def pickdepth(sample,typ,BamPath):
    with open(BamPath+sample+typ+"_wgs_metrics.txt") as fh:
        for lin in fh:
            if lin.startswith("GENOME_TERRITORY"):
                lin2 = fh.next()
                lin2 = lin2.rstrip()
                line = lin2.split()
                MEAN_depth = str(round(float(line[1]),2))
                MEDIAN_depth = str(round(float(line[3]),2))
                COV1X = str(round(float(line[12]),2))
                COV5X = str(round(float(line[13]),2))
                COV10X = str(round(float(line[14]),2))
                COV20X = str(round(float(line[16]),2))
                COV30X = str(round(float(line[18]),2))
                res = ",".join([MEAN_depth,MEDIAN_depth,COV1X,COV5X,COV10X,COV20X,COV30X])
                break
    return res

def duprate(sample,BamPath):
    with open(BamPath+sample+".sorted.rmdup.metrics") as fh:
        for lin in fh:
            if lin.startswith("LIBRARY"):
                lin2 = fh.next()
                lin2 = lin2.rstrip()
                line = lin2.split()
                dup = round(float(line[-2])*100,2)
                dup2 = str(dup)
                break
    return dup2

def insertsize(sample,BamPath):
    with open(BamPath+sample+".insert_size_metrics.txt") as fh:
        for lin in fh:
            if lin.startswith("MEDIAN_INSERT_SIZE"):
                lin2 = fh.next()
                lin2 = lin2.rstrip()
                line = lin2.split()
                size = line[0]
                break
    return str(size)

def judge_fil_exist(sample,RawDic,CleanDic,BamPath): 
    dic = {}
    cmd1 = "/bin/java -Djava.io.tmpdir=/GPFS01/JavaTemp -jar /GPFS01/softwares/picard-tools-2.5.0/picard.jar CollectRawWgsMetrics \
            I=%s%s.%s.bam \
            O=%s%s.%s_wgs_metrics.txt \
            R= /GPFS01/databases/GSCAP/hg19.fa \
            INCLUDE_BQ_HISTOGRAM=true 2> %s%s.%s_wgs_metrics.err"
    cmd2 = "/bin/java -Djava.io.tmpdir=/GPFS01/JavaTemp -jar /GPFS01/softwares/picard-tools-2.5.0/picard.jar MarkDuplicates \
            INPUT=%s%s.sorted.bam \
            OUTPUT=%s%s.sorted.rmdup.Re.bam \
            METRICS_FILE=%s%s.sorted.rmdup.metrics \
            REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true 2> %s%s.rmdup.err"
    cmd3 = "/bin/java -Djava.io.tmpdir=/GPFS01/JavaTemp -jar /GPFS01/softwares/picard-tools-2.5.0/picard.jar CollectInsertSizeMetrics \
            I=%s%s.sorted.bam \
            O=%s%s.insert_size_metrics.txt \
            H=%s%s.insert_size_histogram.pdf \
            M=0.5 2> %s%s.insertsize.err"
    cmd4 = "samtools flagstat %s%s.sorted.bam > %s%s.flagstat 2> %s%s.flagstat.err"
    if os.path.exists(BamPath+sample+".sorted.bam"):
        dic["Bam"] = "yes"
    else:
        print("There is no "+sample+".sorted.bam!")
        exit
    if os.path.exists(BamPath+sample+".sorted.rmdup.bam"):
        dic["DupBam"] = "yes"
    else:
        print("There is no "+sample+".sorted.rmdup.bam!")
        exit
    if os.path.exists(BamPath+sample+".raw_wgs_metrics.txt"):
        dic["RawDepth"] = "yes"
    else:
        os.system(cmd1 % (BamPath,sample,"sorted",BamPath,sample,"raw",BamPath,sample,"raw"))
    if os.path.exists(BamPath+sample+".dup_wgs_metrics.txt"):
        dic["Dupdepth"] = "yes"
    else:
        os.system(cmd1 % (BamPath,sample,"sorted.rmdup",BamPath,sample,"dup",BamPath,sample,"dup"))
    if os.path.exists(BamPath+sample+".sorted.rmdup.metrics"):
        dic["DupRate"] = "yes"
    else:
        os.system(cmd2 % (BamPath,sample,BamPath,sample,BamPath,sample,BamPath,sample))
    if os.path.exists(BamPath+sample+".insert_size_metrics.txt"):
        dic["InsertSize"] = "yes"
    else:
        os.system(cmd3 % (BamPath,sample,BamPath,sample,BamPath,sample,BamPath,sample))
    if os.path.exists(BamPath+sample+".flagstat"):
        dic["MappedRate"] = "yes"
    else:
        os.system(cmd4 % (BamPath,sample,BamPath,sample,BamPath,sample))
    if sample in RawDic:
        dic["Raw"] = "yes"
    else:
        print("There is no "+sample+" RawData Info!")
        exit        
    if sample in CleanDic:
        dic["Clean"] = "yes"
    else:
        print("There is no "+sample+" CleanData Info!")
        exit

def main(argv=None):
    """
    Command Line Usage.
    > python SUMMARY.WGS.SINGLE.py -h
    usage: SUMMARY.WGS.SINGLE.py [-h] [-v] [-V] -ID sample \
    -inPath PathInfo -bamPath PathBam -o PathOut
    
      Created by Mengyue.Zheng on 2019-02-20.
      All rights reserved.
    
    USAGE
    
    optional arguments:
        -h,  --help            Show this help message and exit   
        -v,  --verbose         Set verbosity level [default: True]
        -Id sample, --caseId CaseName
                               Id of the case to be analyzed              
        -iP PathToInfo,  --infoPath PathToInfo
                               Full path to the Infos file [default: ./] 
        -bP PathToBam, --bamPath PathToBam
                               Full path to the Bam file [default: ./]
        -oP PathToOut, --outputPath PathToOut
                               Full path to the output directory [default: ./]
    :Example:
    .. code-block::sh
      python SUMMARY.WGS.SINGLE.py -ID sample \
      -inPath /path/to/info/ -bamPath /path/to/bam -o /path/to/output/ 
    
    """    
    # Command line options.
    # Setup argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", \
                        help="set verbosity level [default: %(default)s]")
    parser.add_argument("-Id","--CaseId",help="Id of the case to be analyzed",required=True)
    parser.add_argument("-iP","--infoPath",help="Full path to the Infos file",default="./")
    parser.add_argument("-bP","--bamPath",help="Full path to the Bam file",default="./")
    parser.add_argument("-oP","--outputPath",help="Full path to the output directory",default="./")
    # Process arguments
    args = parser.parse_args()
    if len(sys.argv) < 1:
        print("Command is not right!")
    else:
        header = ['#SAMPLE','LENGTH','GC(%)','Q20(%)','Q30(%)','PF_READS','PF_BASES',
                  'CLEAN_READS','CLEAN_BASES','RATIO_OF_READS(%)','RATIO_OF_BASES(%)',
                  'INSERT_SIZE','DUPLICATE(%)','RATIO_OF_MAPPED(%)','MEAN_DEPTH',
                  'MEDIAN_DEPTH','1X_COVERAGE(%)','5X_COVERAGE(%)','10X_COVERAGE(%)',
                  '20X_COVERAGE(%)','30X_COVERAGE(%)','MEAN_DEPTH_DEDUP','MEDIAN_DEPTH_DEDUP',
                  '1X_COVERAGE_DEDUP(%)','5X_COVERAGE(%)','10X_COVERAGE(%)',
                  '20X_COVERAGE_DEDUP(%)','30X_COVERAGE(%)']
        dicR1,dicR2 = fil2dic(args.infoPath+"RawInfos.csv")
        dicR3,dicR4 = fil2dic(args.infoPath+"CleanInfos.csv")
        RawDic = sumdic(dicR1,dicR2)
        CleanDic = sumdic(dicR3,dicR4)
        judge_fil_exist(args.CaseId,RawDic,CleanDic,args.bamPath)
        if pickdepth(args.CaseId,".raw",args.bamPath):
            RawDepth = pickdepth(args.CaseId,".raw",args.bamPath)
        else:
            RawDepth = "NA,NA,NA,NA,NA,NA,NA"
        if pickdepth(args.CaseId,".dup",args.bamPath):
            DupDepth = pickdepth(args.CaseId,".dup",args.bamPath)
        else:
            DupDepth = "NA,NA,NA,NA,NA,NA,NA"
        if duprate(args.CaseId,args.bamPath):
            DupRate = duprate(args.CaseId,args.bamPath)
        else:
            DupRate = "NA"
        if pickmapped(args.CaseId,args.bamPath):
            MapRate = pickmapped(args.CaseId,args.bamPath)
        else:
            MapRate = "NA"
        if insertsize(args.CaseId,args.bamPath):
            InsertSize = insertsize(args.CaseId,args.bamPath)
        else:
            InsertSize = "NA"
        with open(args.outputPath+args.CaseId+".SUMMARY.csv","w") as fw: 
            fw.write(",".join(header)+"\n")
            rawres = RawDic[args.CaseId].split(",")
            cleanres = CleanDic[args.CaseId].split(",")
            readrate = str(round((float(cleanres[-2])/float(rawres[-2]))*100,2))
            baserate = str(round((float(cleanres[-1])/float(rawres[-1]))*100,2))
            cleand = ",".join([cleanres[-2],cleanres[-1],readrate,baserate])
            fw.write(args.CaseId+","+RawDic[args.CaseId]+","+cleand+","+InsertSize+","+DupRate+","+MapRate+","+RawDepth+","+DupDepth+"\n")
if __name__ == '__main__':
    main()

#!/usr/bin/env python
#-*-coding:utf-8
import sys
import pyfaidx
import argparse
import operator

global reference,movelst,prtlst
reference = pyfaidx.Fasta("/GPFS02/zhengmy/code/hs37d5.fa")
movelst = []
prtlst = []

# import CPX file， create diccpx\diclst
def importCPX(fil):
    diclst = {}
    diccpx = {}
    diccount = {}
    dicend = {}
    with open(fil) as fh:
        for i in fh:
            i = i.rstrip()
            line = i
            if i.startswith("id"):
                continue
            if line in prtlst:
                next
            else:
                prtlst.append(line)
            t = i.split(",")
            info = t[5],t[6],t[7],t[8]
            key1 = t[5]
            key2 = t[6]           
            pos1 = info[0].split(":")
            pos2 = info[1].split(":")
            cpxchr,cpxsta = pos1[0],int(pos1[1])
            cpxend = int(pos2[1])
            tupcpx = [cpxchr,cpxsta,cpxend,info[2],info[3]]
            diclst[line] = tupcpx
            endlst = dicend.values()
            if key1 in diccpx or key1 in endlst:
                diccount[key1] = 2
                diccount[key2] = 2
            elif key2 in diccpx or key2 in endlst:
                diccount[key1] = 2
                diccount[key2] = 2
            else:
                diccpx[key1] = tupcpx
                dicend[key1] = key2
                diccount[key1] = 1
    rmlst = []
    for keym in diccpx:
            if diccount[keym] == 1:
                next
            else:
                rmlst.append(keym)
    for rmid in rmlst:
        del diccpx[rmid]
    return diclst,diccpx,dicend

# import raw and mut result
def importMut(fil,mutlst):
    with open(fil) as fh:
        for i in fh:
            i = i.rstrip()
            if i.startswith("id"):
                continue
            t = i.split(",")
            info = t[5],t[6],t[7],t[8]
            pos1 = info[0].split(":")
            pos2 = info[1].split(":")
            cpxchr,cpxsta = pos1[0],int(pos1[1])
            cpxend = int(pos2[1])
            tupmut =  (cpxchr,cpxsta,cpxend,info[2],info[3])
            mutlst.append(tupmut)            
    return mutlst

# sort chr of raw and mut result
def sortMutchr(mutlst):
    dicchr = {}
    chrlst = ["1","2","3","4","5","6",\
              "7","8","9","10","11","12",\
              "13","14","15","16","17",\
              "18","19","20","21","22","X","Y"]
    for ch in chrlst:
        dicchr[ch] = []
    for i in mutlst:
        CHR = i[0].replace("chr","")
        if CHR in dicchr:
            dicchr[CHR].append(i)
        else:
            next           
    return dicchr

# sort possta of raw and mut result
def sortMutpos(dicchr):
    SORTMUT = []
    sorted_dicchr = {}
    for  i in dicchr:
        dicpos = {}
        for n in dicchr[i]:
            pos = n[1]
            dicpos[n] = pos
        sorted_dicchr[i] = sorted(dicpos.items(),key=operator.itemgetter(1))
    for n in sorted_dicchr:
        for m in sorted_dicchr[n]:
            SORTMUT.append(m[0])
    return SORTMUT

#grouing sub_mutation
def SubMut(SORTMUT,diccpx,dicend):
    dicSubMut = {}
    for i in SORTMUT:
        line = i
        for n in range(-5,5):
            pos1 = i[1]+n
            pos2 = i[2]+n
            key1 = str(i[0])+":"+str(pos1)
            key2 = str(i[0])+":"+str(pos2)
            for m in diccpx:
                if key1 == m or key1 == dicend[m]:
                    if m in dicSubMut:
                        if line in dicSubMut[m]:
                            next
                        else:
                            dicSubMut[m].append(line)
                    else:
                        dicSubMut[m] = [line]
                if key2 == m or key2 == dicend[m]:
                    if m in dicSubMut:
                        if line in dicSubMut[m]:
                            next
                        else:
                            dicSubMut[m].append(line)
                    else:
                        dicSubMut[m] = [line]
    return dicSubMut

# LeftAigned and Deladjust
def Deladjust(MUT,CPX,tag):
    CHR = CPX[0].replace("chr","")
    n = 0
    if MUT[3] == "-":
        if tag == "Before":
           STA2 = MUT[1]-1
           END2 = MUT[2]-1
           REF2 = MUT[3]
           ALT2 = MUT[4]
           n = 1
        elif tag == "After":
           STA2 = MUT[1]+1
           END2 = MUT[2]+1
           REF2 = MUT[3]
           ALT2 = MUT[4]
           n = 1
        else:
           STA2 = MUT[1]
           END2 = MUT[2]
           REF2 = MUT[3]
           ALT2 = MUT[4]
           n = 1
    else:
        STA2 = MUT[1]   
        END2 = MUT[2]
        REF2 = MUT[3]
        ALT2 = MUT[4]
    if MUT[1] < CPX[1]:
        if n == 1:
            if MUT[1] < (CPX[1]-1):
                ref1 = reference[CHR][MUT[1]:CPX[1]-1].seq.upper()
                ref2 = reference[CHR][MUT[2]:(MUT[2]+(CPX[1]-MUT[1])-1)].seq.upper()
                if ref1 == ref2:                 
                    ref = reference[CHR][MUT[2]:(MUT[2]+(CPX[1]-MUT[1])-1)].seq.upper()
                    REF = "-"
                    END = END2+(CPX[1]-STA2)-2
                    STA = STA2+(CPX[1]-STA2)-2
                    ALT_tmp = ALT2.lstrip(ref)
                    ALT = ALT_tmp+reference[CHR][STA:END].seq.upper()
                else:
                    STA = STA2
                    REF = REF2
                    END = END2
                    ALT = ALT2
            else:
                STA = STA2
                REF = REF2
                END = END2
                ALT = ALT2
                
        else:
            ref1 = reference[CHR][MUT[1]-1:CPX[1]-1].seq.upper()
            ref2 = reference[CHR][MUT[2]:(MUT[2]+(CPX[1]-MUT[1]))].seq.upper()
            if ref1 == ref2:                 
                ref = reference[CHR][(CPX[1]-1):(MUT[2]+(CPX[1]-MUT[1]))].seq.upper()          
                REF = ref
                END = END2+(CPX[1]-STA2)
                STA = STA2+(CPX[1]-STA2)
                ALT = ALT2
            else:
                STA = STA2
                REF = REF2
                END = END2
                ALT = ALT2
    else:
        STA = STA2
        END = END2
        REF = REF2
        ALT = ALT2
    if MUT[-1] == 1:
        n = 1 
    MUTadj = [MUT[0],STA,END,REF,ALT,n]
    return MUTadj

# If Have_gap
def MergeGap(MUT1,MUT2):
    CHR = MUT1[0].replace("chr","")
    n = 0
    if MUT1[-1] == 1:
       n = 1
    elif MUT2[-1] == 1:
       n = -1
    else:
       n = 0
    if MUT1[2] >= MUT2[1]:
        MUT3 = "NoMerge"
    elif (MUT1[2]+1) == MUT2[1]:
        ref = MUT1[3].replace("-","")+MUT2[3].replace("-","")
        alt = MUT1[4].replace("-","")+MUT2[4].replace("-","")
        MUT3 = [MUT1[0],MUT1[1],MUT2[2],ref,alt,n]
    else:
        refgap = reference[CHR][MUT1[2]:MUT2[1]-1].seq.upper()
        ref = MUT1[3].replace("-","")+refgap+MUT2[3].replace("-","")
        alt = MUT1[4].replace("-","")+refgap+MUT2[4].replace("-","")
        MUT3 = [MUT1[0],MUT1[1],MUT2[2],ref,alt,n]
    return MUT3

# MergeMut    
def Merge2Mut(SubMut,CPX):        
    MUT1 = Deladjust(SubMut[0],CPX,"Before")
    MUT2 = Deladjust(SubMut[1],CPX,"After")
    MergeMut = MergeGap(MUT1,MUT2)
    CPX2 = CPXMove(MUT1,MUT2,CPX)
    if MergeMut == CPX2:
        movelst.append(SubMut[0])
        movelst.append(SubMut[1])

def CPXMove(MUT1,MUT2,CPX):
    if MUT1[3] == "-" or MUT1[-1] == 1:
        STA = CPX[1]-2
    else:
        STA = CPX[1]
    if MUT2[3] == "-" or MUT2[-1] == 1:
        END = CPX[2]+2
    else:
        END = CPX[2]
    CPX2  = [CPX[0],STA,END,CPX[3],CPX[4],MUT1[-1]]
    return CPX2
       
def Merge3Mut(SubMut,CPX,l):
  m = (l-2)
  for i in range(m):
    MUT1b = Deladjust(SubMut[i+0],CPX,"Before")
    MUT2a = Deladjust(SubMut[i+1],CPX,"After") 
    CPX2 = CPXMove(MUT1b,MUT2a,CPX)
    MergeMut1 = MergeGap(MUT1b,MUT2a)
    if MergeMut1 == CPX2:
        movelst.append(SubMut[i+0])
        movelst.append(SubMut[i+1])   
    else:
        MUT3a = Deladjust(SubMut[i+2],CPX,"After")
        MergeMut2 = MergeGap(MUT1b,MUT3a)
        CPX2 = CPXMove(MUT1b,MUT3a,CPX)
        if MergeMut2 == CPX2:
            movelst.append(SubMut[i+0])
            movelst.append(SubMut[i+2])        
        else: 
             MUT2b = Deladjust(SubMut[1],CPX,"Before")   
             MergeMut3 = MergeGap(MUT2b,MUT3a)
             CPX2 = CPXMove(MUT2b,MUT3a,CPX)
             if MergeMut3 == CPX2:
                 movelst.append(SubMut[1])
                 movelst.append(SubMut[2])
             else:
                 MUT1 = Deladjust(SubMut[0],CPX,"NA")
                 MUT2 = Deladjust(SubMut[1],CPX,"NA")
                 MergeMut1re = MergeGap(MUT1,MUT2)
                 MergeMut1re2 = Deladjust(MergeMut1re,CPX,"Before")
                 MergeMut4 = MergeGap(MergeMut1re2,MUT3a)
                 CPX2 = CPXMove(MergeMut1re2,MUT3a,CPX)
                 if MergeMut4 == CPX2:
                     movelst.append(SubMut[0])
                     movelst.append(SubMut[1])
                     movelst.append(SubMut[2])
                 else:
                     next

#def Merge4Mut(SubMut,CPX,l):
#    m = (l-3)
#    for i in range(m):
#        MUT1b = Deladjust(SubMut[i+0],CPX,"Before")
#        MUT2a = Deladjust(SubMut[i+1],CPX,"After")
#        MUT2b = Deladjust(SubMut[i+1],CPX,"Before")
#        MUT3a = Deladjust(SubMut[i+2],CPX,"After")
#        MUT3b = Deladjust(SubMut[i+2],CPX,"Before")
#        MUT4a = Deladjust(SubMut[i+3],CPX,"After")
#        MergeMut1 = MergeGap(MUT1b,MUT2a)
#        MergeMut2 = MergeGap(MUT1b,MUT3a)
#        MergeMut3 = MergeGap(MUT1b,MUT4a)       
#        MergeMut4 = MergeGap(MUT2b,MUT3a)
#        MergeMut5 = MergeGap(MUT2b,MUT4a)
#        MergeMut1re = MergeGap(MUT1b,MUT2b)
#        MergeMut6 = MergeGap(MergeMut1re,MUT3a)
#        MergeMut7 = MergeGap(MergeMut1re,MUT4a)
#        MergeMut2re = MergeGap(MUT1b,MUT3b)
#        MergeMut8 = MergeGap(MergeMut2re,MUT4a)
#        MergeMut4re = MergeGap(MUT2b,MUT3b)
#        MergeMut9 = MergeGap(MergeMut4re,MUT4a)
#        MergeMut6re = MergeGap(MergeMut1re,MUT3b)
#        MergeMut10 = MergeGap(MergeMut6re,MUT4a)
#        if MergeMut1 == CPX:
#            movelst.append(SubMut[i+0])
#            movelst.append(SubMut[i+1])
#        elif MergeMut2 == CPX:
#            movelst.append(SubMut[i+0])
#            movelst.append(SubMut[i+2])
#        elif MergeMut3 == CPX:
#            movelst.append(SubMut[i+0])
#            movelst.append(SubMut[i+3])
#        elif MergeMut4 == CPX:
#            movelst.append(SubMut[i+1])
#            movelst.append(SubMut[i+2])
#        elif MergeMut5 == CPX:
#            movelst.append(SubMut[i+1])
#            movelst.append(SubMut[i+3])
#        elif MergeMut6 == CPX:
#            movelst.append(SubMut[i+0])
#            movelst.append(SubMut[i+1])
#            movelst.append(SubMut[i+2])
#        elif MergeMut7 == CPX:
#            movelst.append(SubMut[i+0])
#            movelst.append(SubMut[i+1])
#            movelst.append(SubMut[i+3])
#        elif MergeMut8 == CPX:
#            movelst.append(SubMut[i+0])
#            movelst.append(SubMut[i+2])
#            movelst.append(SubMut[i+3])
#        elif MergeMut9 == CPX:
#            movelst.append(SubMut[i+1])
#            movelst.append(SubMut[i+2])
#            movelst.append(SubMut[i+3])
#        elif MergeMut10 == CPX:
#            movelst.append(SubMut[i+1])
#            movelst.append(SubMut[i+2])
#            movelst.append(SubMut[i+3])
#            movelst.append(SubMut[i+4])
#        else:
#            next

def SubMerge(dicSubMut,diccpx):
    for i in dicSubMut:
        l = len(dicSubMut[i])
        if l < 2:
            next
        elif l == 2:
            Merge2Mut(dicSubMut[i],diccpx[i])
        elif l >= 3:
            Merge3Mut(dicSubMut[i],diccpx[i],l)
#        elif l >= 4:
#            Merge4Mut(dicSubMut[i],diccpx[i],l)
        else:
            next

def FilterMut(fil,movelst):
    with open(fil) as fh:
        for i in fh:
            i = i.rstrip()
            line = i
            if i.startswith("id"):
                continue
            t = i.split(",")
            info = t[5],t[6],t[7],t[8]
            pos1 = info[0].split(":")
            pos2 = info[1].split(":")
            cpxchr,cpxsta = pos1[0],int(pos1[1])
            cpxend = int(pos2[1])
            tupmut =  (cpxchr,cpxsta,cpxend,info[2],info[3])
            if tupmut in movelst:
                next
            else:
                if line in prtlst:
                    next
                else:
                    prtlst.append(line)

def main(argv=None):
    """
    Command Line Usuage.
    > python CPXDetection_hs37d5.py -h
    Usuage: python CPXDetection_hs37d5.py [-h] [-v] -id sample \
    -cpxpath CPX.csv -mutpath MUT.csv -rawpath mut_raw.csv -logpath \
    log.file -op RES.csv
    
        Created by Mengyue.Zheng on 2019-06-25
        All right reserved.
    USUAGE
    
    optional arguments:
        -h,  --help    show this help message and exit
        -v,  --verbose    Set verbosity level [default: True]
        -id sample,  --caseid  CaseName    id of the case to be analyzed
        -cp  pathtocpx,  --cpxpath  pathtocpx    Full path to the CPX.csv [default:./]
        -mp  pathtomut,  --mutpath  pathtomut    Full path to the MUT.csv [default:./]
        -rp  pathtoraw,  --rawpath  pathtoraw    Full path to the mut_raw.csv [default:./]
        -lp  pathtolog,  --logpath  pathtolog    Full path to the log file [default:./]
        -op  pathtores,  --respath  pathtores    Full path to the output directory [default:./]  
    :Example:
    ..code-block::sh
        python CPXDetection_hs37d5.py -id sample -cp /path/CPX.csv \
        -mp /path/MUT.csv -rp /path/mut_raw.csv -lp /path/log \
        -o /path/RES.csv
    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-v","--verbose",dest="verbose",action="store_true",\
                        help="set verbosity level [default:%(default)s]")
    parser.add_argument("-id","--caseid",help="id of the case to be analyzed",required=True)
    parser.add_argument("-cp","--cpxpath",help="Full path to the CPX.csv",default="./")
    parser.add_argument("-mp","--mutpath",help="Full path to the MUT.csv",default="./")
    parser.add_argument("-rp","--rawpath",help="Full path to the mut_raw.csv",default="./")
    parser.add_argument("-lp","--logpath",help="Full path to the log file",default="./")
    parser.add_argument("-op","--respath",help="Full path to the output directory",default="./")
    args = parser.parse_args()
    if len(sys.argv) <1:
        print("Command is not right!")
    else:
        header = ["id","Type","Gene","Gene.ID","AAChange","Chr.start","Chr.end","Ref","Alt","Hom.Het",\
                  "ExonicFunc","rsID","AF","copynumber","DP(ref:alt)","NormalID","NormalAF",\
                  "NormalDP(ref:alt)","Hotspot","CLNSIG","CLNREVSTAT","CLNDNINCL","CLINID",\
                  "CosmicID","Cosmic.Occurence","HGMD_class","Variant_effect_LOVD","InterVar_Class",\
                  "1000g","1000gEAS","ExAC_ALL","ExAC_EAS","gnomAD_exome_ALL","gnomAD_exome_EAS",\
                  "gnomAD_genome_ALL","gnomAD_genome_EAS","GS_ALL","SIFT.score","SIFT.pred",\
                  "PolyPhen.score","PolyPhen.pred","CADD_phred","GERP++_RS","OMIM_Phenotype",\
                  "OMIM_Inheritance","OMIM_Link","InterVar_ACMG","OrphaNumber","CANONICAL",\
                  "AAChange(Annovar)","BrkpFusion","BrkptType","TumorReferenceCount",\
                  "TumorSplitReferenceCount","TumorVariantCount","TumorSplitVariantCount",\
                  "Cosmic_Fusion_Counts","CC_Tumour_Types(Somatic)","Orientation"]
        mutlst = []
        MUTLST_tmp = importMut(args.mutpath+args.caseid+".MUT.csv",mutlst)
        MUTLST = importMut(args.rawpath+args.caseid+".mut_raw.csv",MUTLST_tmp)
        dicchr = sortMutchr(MUTLST)
        SORTMUT = sortMutpos(dicchr)
        diclst,diccpx,dicend = importCPX(args.cpxpath+args.caseid+".CPX.csv")
        dicSubMut = SubMut(SORTMUT,diccpx,dicend)
        SubMerge(dicSubMut,diccpx)
        FilterMut(args.mutpath+args.caseid+".MUT.csv",movelst)
        with open(args.respath+args.caseid+".RES.csv","w") as fw:
            fw.write(",".join(header)+"\n")
            for i in prtlst:
                fw.write(i+"\n")
        with open(args.logpath+args.caseid+".log","w") as fw1:
            for i in movelst:
                mvl = list(map(str,i))
                fw1.write(",".join(mvl)+"\n") 
if __name__ == '__main__':
    main()

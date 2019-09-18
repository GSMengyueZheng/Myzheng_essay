#!/usr/bin/env python
#-*-coding:utf-8
import sys
import pyfaidx
import argparse
import operator
import os

global reference,movelst,prtlst,prtcpxlst,cpxmglst,prtlstrm,numlst
reference = pyfaidx.Fasta("/GPFS02/zhengmy/code/hs37d5.fa")
movelst,prtlst,prtcpxlst,cpxmglst,prtlstrm,numlst = [],[],[],[],[],[]


#import CPX file， create diccpx\diclst
def importCPX(fil):
    diclst,diccpx,diccount,dicend,cpxlst  = {},{},{},{},[]
    with open(fil) as fh:
        for i in fh:
          if i.startswith("id"):
              continue
          line = i.rstrip()
          t = i.rstrip().split(",")
          key1,key2 = t[5],t[6]
          if "%" in t[12]:
              AF = float(t[12].replace("%",""))
          else:
              AF = "None"
          if t[5] == "." or t[6] == ".":
              next
          else:
            pos1 = t[5].split(":")
            pos2 = t[6].split(":")
            cpxchr,cpxsta = pos1[0],int(pos1[1])
            cpxend = int(pos2[1])
            tupcpx = [cpxchr,cpxsta,cpxend,t[7],t[8],AF,"CPX"]
            cpxlst.append(tupcpx)
            diclst[line] = tupcpx
            endlst = dicend.values()
            if key1 in diccpx or key1 in endlst:
                diccount[key1],diccount[key2] = 2,2
            elif key2 in diccpx or key2 in endlst:
                diccount[key1],diccount[key2] = 2,2
            else:
                diccpx[key1],dicend[key1],diccount[key1] = tupcpx,key2,1
    rmlst = []
    for keym in diccpx:
            if diccount[keym] == 1:
                next
            else:
                rmlst.append(keym)
    for rmid in rmlst:
        del diccpx[rmid]
    return diclst,diccpx,dicend,cpxlst

# import raw and mut result
def importMut(fil,mutlst,cpxlst,frm):
    with open(fil) as fh:
        for i in fh:
          if i.startswith("id"):
              continue
          t = i.split(",")
          if "%" in t[12]:
              AF = float(t[12].replace("%",""))
          else:
              AF = "None"
          if t[5] == "." or t[6] == ".":
	      next
          else:
            refl = len(t[7])
            altl = len(t[8])
            pos1 = t[5].split(":")
            pos2 = t[6].split(":")
            cpxchr,cpxsta = pos1[0],int(pos1[1])
            cpxend = int(pos2[1])
            if frm == "MUT":
              if refl >= 50 or altl >= 50:
                if refl >= 50 and altl >= 50:
                    tupcpx2 = [cpxchr,cpxsta,cpxend,str(refl),str(altl),"CPX"]
                    if tupcpx2 in cpxlst:
                        numlst.append(tupcpx2)
                    else:
                        next
                elif refl >= 50:
                    tupcpx2 = [cpxchr,cpxsta,cpxend,str(refl),t[8],"CPX"]
                    if tupcpx2 in cpxlst:
                        numlst.append(tupcpx2)
                    else:
                        next
                else:
                    tupcpx2 = [cpxchr,cpxsta,cpxend,t[7],str(altl),"CPX"]
                    if tupcpx2 in cpxlst:
                        numlst.append(tupcpx2)
                    else:
                        next
            else:
              next
            tupmut =  (cpxchr,cpxsta,cpxend,t[7],t[8],AF,"MUT")
            if tupmut in mutlst:
                next
            else:
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
        for n in range(-10,10):
            pos1 = i[1]+n
            pos2 = i[2]+n
            key1 = str(i[0])+":"+str(pos1)
            key2 = str(i[0])+":"+str(pos2)
            ref,alt = i[3],i[4]
            refl,altl = len(i[3]),len(i[4])
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
  if MUT == "NoMerge":
        MUTadj = "NoMerge"
  else:  
    CHR = CPX[0].replace("chr","")
    n = 0
    if MUT[-1] == -2:
       STA2 = MUT[1]
       END2 = MUT[2]-2
       REF2 = MUT[3]
       ALT2 = MUT[4]
       n = 2 
    elif MUT[3] == "-" or MUT[-1] == -1:
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
        elif tag == "Forth":
           STA2 = MUT[1]-2
           END2 = MUT[2]-2
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
                ref3 = reference[CHR][MUT[1]:MUT[2]].seq.upper()
                if ref1 == ref2:                 
                    ref = reference[CHR][MUT[2]:(MUT[2]+(CPX[1]-MUT[1])-1)].seq.upper()
                    REF = "-"
                    END = END2+(CPX[1]-STA2)-2
                    STA = STA2+(CPX[1]-STA2)-2
                    ALT_tmp = ALT2.lstrip(ref)
                    ALT = ALT_tmp+reference[CHR][STA:END].seq.upper()
                elif ref3 == ALT2:
                    STA = MUT[1]+(MUT[2]-MUT[1])
                    REF = "-"
                    END = MUT[2]
                    ALT = ALT2
                    n = 2
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
    if tag == "Third" and MUT[-1] == -1:
        n = -3
    #if tag == "Forth" and MUT[-1] == -1:
    #    n = -2
    MUTadj = [MUT[0],STA,END,REF,ALT,MUT[5],MUT[6],n]
  return MUTadj

# If Have_gap
def MergeGap(MUT1,MUT2): 
  if MUT1 == "NoMerge" or MUT2 == "NoMerge":
     MUT3 = "NoMerge"
  else:
    CHR = MUT1[0].replace("chr","")
    if MUT1[-1] == 1:
       n = 1
    elif MUT2[-1] == 1:
       n = -1
    else:
        n = MUT1[-1]
    if MUT1[2] >= MUT2[1]:
        MUT3 = "NoMerge"
    elif (MUT1[2]+1) == MUT2[1]:
        if n == -3:
            refgap = reference[CHR][MUT1[2]-2:MUT2[1]-1].seq.upper()
            ref = MUT1[3].replace("-","")+refgap+MUT2[3].replace("-","")
            alt = MUT1[4].replace("-","")+refgap+MUT2[4].replace("-","")
            MUT3 = [MUT1[0],MUT1[1],MUT2[2],ref,alt,"AF","NA",n]
        else:
           ref = MUT1[3].replace("-","")+MUT2[3].replace("-","")
           alt = MUT1[4].replace("-","")+MUT2[4].replace("-","")
           MUT3 = [MUT1[0],MUT1[1],MUT2[2],ref,alt,"AF","NA",n]
    else:
        if n == -3:
            refgap = reference[CHR][MUT1[2]-2:MUT2[1]-1].seq.upper()
            ref = MUT1[3].replace("-","")+refgap+MUT2[3].replace("-","")
            alt = MUT1[4].replace("-","")+refgap+MUT2[4].replace("-","")
            MUT3 = [MUT1[0],MUT1[1],MUT2[2],ref,alt,"AF","NA",n]
        else:
            refgap = reference[CHR][MUT1[2]:MUT2[1]-1].seq.upper()
            ref = MUT1[3].replace("-","")+refgap+MUT2[3].replace("-","")
            alt = MUT1[4].replace("-","")+refgap+MUT2[4].replace("-","")
            MUT3 = [MUT1[0],MUT1[1],MUT2[2],ref,alt,"AF","NA",n]
  return MUT3


###分3组的：
def getRealSub(l):
    realList = []
    for i in range(0,l):
        for j in range(0,l):
            for z in range(0,l):
                for k in range(0,l):
                    if i == j or j == z or i == z or i == k or j == k or z == k:
                        next
                    else:
                        if i <= j <= z <= k:
                            realList.append([i,j,z,k])
                        else:
                            next
    return realList

def checkAF(MUT1,MUT2):
    if MUT1[5] == "None" or MUT2[5] == "None":
        AF = "NoSuit"
    elif (MUT1[5]-MUT2[5]) > 30 and MUT1[6] == "MUT":
        AF = "NoSuit"
    elif (MUT2[5]-MUT1[5]) > 30 and MUT2[6] == "MUT":
        AF = "NoSuit"
    else:
        AF = "Suit"
    return AF

# MergeMut    
def Merge2Mut(SubMut,CPX):
    MUT1 = Deladjust(SubMut[0],CPX,"Before")
    MUT2 = Deladjust(SubMut[1],CPX,"After")
    AFStat = checkAF(MUT1,MUT2)
    if AFStat == "Suit":
       MergeMut = MergeGap(MUT1,MUT2)
       CPX2 = CPXMove(MUT1,MUT2,CPX)
       if MergeMut[0:5] == CPX2[0:5]:
           movelst.append(SubMut[0])
           movelst.append(SubMut[1])
           cpxmglst.append(CPX)
    else:
        next

def CPXMove(MUT1,MUT2,CPX):
    if MUT1[-1] == 2:
       STA = CPX[1]-1
    elif MUT1[3] == "-" or MUT1[-1] == 1:
        STA = CPX[1]-2
    else:
        STA = CPX[1]
    if MUT2[3] == "-" or MUT2[-1] == 1:
        END = CPX[2]+2
    else:
        END = CPX[2]
    CPX2  = [CPX[0],STA,END,CPX[3],CPX[4],CPX[5],CPX[6],MUT1[-1]]
    return CPX2

def Merge3Mut(SubMut,CPX):
    State = ""
    i,j,k = 0,1,2
    MUT1b = Deladjust(SubMut[i],CPX,"Before")
    MUT2a = Deladjust(SubMut[j],CPX,"After")    
    CPX2 = CPXMove(MUT1b,MUT2a,CPX)
    AFStat = checkAF(MUT1b,MUT2a)
    MergeMut1 = MergeGap(MUT1b,MUT2a)
    if AFStat == "Suit" and MergeMut1[0:5] == CPX2[0:5]:
            movelst.append(SubMut[i])
            movelst.append(SubMut[j])
            cpxmglst.append(CPX)
            State = "Merge"
    else:
        MUT3a = Deladjust(SubMut[k],CPX,"After")
        AFStat = checkAF(MUT1b,MUT3a)
        MergeMut2 = MergeGap(MUT1b,MUT3a)
        CPX2 = CPXMove(MUT1b,MUT3a,CPX)
        if AFStat == "Suit" and MergeMut2[0:5] == CPX2[0:5]:
            movelst.append(SubMut[i])
            movelst.append(SubMut[k])
            cpxmglst.append(CPX)
            State = "Merge"
        else: 
             MUT2b = Deladjust(SubMut[j],CPX,"Before")
             AFStat = checkAF(MUT2b,MUT3a)
             MergeMut3 = MergeGap(MUT2b,MUT3a)
             CPX2 = CPXMove(MUT2b,MUT3a,CPX)
             if AFStat == "Suit" and MergeMut3[0:5] == CPX2[0:5]:
                 movelst.append(SubMut[j])
                 movelst.append(SubMut[k])
                 cpxmglst.append(CPX)
                 State = "Merge"
             else:
                 MUT1 = Deladjust(SubMut[i],CPX,"Before")
                 MUT2 = Deladjust(SubMut[j],CPX,"NA")                 
                 MergeMut1re = MergeGap(MUT1,MUT2)
                 MergeMut1re2 = Deladjust(MergeMut1re,CPX,"Before")
                 AFStat = checkAF(MUT1,MUT2)
                 AFStat2 = checkAF(MUT1,MUT3a)
                 AFStat3 = checkAF(MUT2,MUT3a)
                 MergeMut4 = MergeGap(MergeMut1re2,MUT3a)
                 CPX2 = CPXMove(MergeMut1re2,MUT3a,CPX)
                 if AFStat == "Suit" and AFStat2 == "Suit" and AFStat3 == "Suit" and MergeMut4[0:5] == CPX2[0:5]:
                     movelst.append(SubMut[i])
                     movelst.append(SubMut[j])
                     movelst.append(SubMut[k])
                     cpxmglst.append(CPX)
                     State = "Merge"
                 else:
                    next
    return State

def Merge4Mut(SubMut,CPX,l):
  realList = getRealSub(l)
  State = ""
  for lst in realList:
    i,j,k,v = lst
    State = Merge3Mut([SubMut[i],SubMut[j],SubMut[k]],CPX)
    if State == "Merge":
        next
    else:
         State == Merge3Mut([SubMut[j],SubMut[k],SubMut[v]],CPX)
    if State == "Merge":
        next
    else:
        MUT1b = Deladjust(SubMut[i],CPX,"Before")
        MUT2b = Deladjust(SubMut[j],CPX,"After")
        MUT3b = Deladjust(SubMut[k],CPX,"After")
        MUT4a = Deladjust(SubMut[v],CPX,"After")
        MergeMut1 = MergeGap(MUT1b,MUT2b)       
        MergeMut1re = Deladjust(MergeMut1,CPX,"Third")
        MergeMut2 = MergeGap(MergeMut1re,MUT3b)
        MergeMut2re = Deladjust(MergeMut2,CPX,"Forth")
        MergeMut3 = MergeGap(MergeMut2re,MUT4a)
        CPX2 = CPXMove(MergeMut2re,MUT4a,CPX)
        if MergeMut3[0:5] == CPX2[0:5]:
            movelst.append(SubMut[i])
            movelst.append(SubMut[j])
            movelst.append(SubMut[k])
            movelst.append(SubMut[v])
            cpxmglst.append(CPX)
            State = "Merge"
        else:
            next

def SubMerge(dicSubMut,diccpx):
    for i in dicSubMut:
        l = len(dicSubMut[i])
        if l < 2:
            next
        elif l == 2:
            Merge2Mut(dicSubMut[i],diccpx[i])
        elif l == 3:
            Merge3Mut(dicSubMut[i],diccpx[i])
        elif l >= 4:
            Merge4Mut(dicSubMut[i],diccpx[i],l)
        else:
            next

def is_int(str):
    try:
        int(str)
        return True
    except ValueError:
        return False

def FilterMut(fil,movelst,TMBState):
    with open(fil) as fh:
        for i in fh:
          i = i.rstrip()
          line = i
          if i.startswith("id"):
                continue
          t = i.split(",")
          if "%" in t[12]:
                AF = float(t[12].replace("%",""))
          else:
                AF = "None"
          if t[5] == "." or t[6] == ".":
                if line in prtlst:
                    next
                else:
                    prtlst.append(line)
          else:
            pos1 = t[5].split(":")
            pos2 = t[6].split(":")
            cpxchr,cpxsta = pos1[0],int(pos1[1])
            cpxend = int(pos2[1])
            tupmut =  (cpxchr,cpxsta,cpxend,t[7],t[8],AF,"MUT")
            if TMBState == "TMB_ori":
               if is_int(t[7]) == True  or is_int(t[8]) == True or len(t[7]) >= 50 or len(t[8]) >= 50:
                   FLG = "TMB"
               else:
                   FLG = "No"
            else:
               FLG = "No"
            if FLG == "TMB":
               next
            else:
              if tupmut in movelst:
                prtlstrm.append("MUT"+","+line)
                #prtlst.append(line+","+"MUTFilter")
              else:
                if line in prtlst:
                    next
                else:
                    #prtlst.append(line+","+"MUTOK")
                    prtlst.append(line)

def importCpx2(fil,TMBState):
    with open(fil) as fh:
        for i in fh:
          i = i.rstrip()
          line = i
          if i.startswith("id"):
                continue
          t = i.split(",")
          if "%" in t[12]:
                AF = float(t[12].replace("%",""))
          else:
                AF = "None"
          if t[5] == "." or t[6] == ".":
                if line in prtcpxlst:
                    next
                else:
                    prtcpxlst.append(line)
          else:
            pos1 = t[5].split(":")
            pos2 = t[6].split(":")
            cpxchr,cpxsta = pos1[0],int(pos1[1])
            cpxend = int(pos2[1])
            tupmut =  [cpxchr,cpxsta,cpxend,t[7],t[8],AF,"CPX"]
            if TMBState == "TMB_ori":
               if is_int(t[7]) == True  or is_int(t[8]) == True or len(t[7]) >= 50 or len(t[8]) >= 50:
                   FLG = "TMB"
               else:
                   FLG = "No"
            else:
               FLG = "No"
            if FLG == "TMB":
               next
            else:
              if tupmut in numlst:
                #prtlst.append(line+","+"CPXFilter")
                 next
              else:
                if line in prtcpxlst:
                    next
                else:
                    #prtlst.append(line+","+"CPXOK")
                     prtcpxlst.append(line)
              if tupmut in cpxmglst:
                  prtlstrm.append("CPX"+","+line)
              else:
                next

# import raw and mut result
def importRaw(fil):
    with open(fil) as fh:
        for i in fh:
          i = i.rstrip()
          line = i
          if i.startswith("id"):
                continue
          t = i.split(",")
          if "%" in t[12]:
                AF = float(t[12].replace("%",""))
          else:
                AF = "None"
          if t[5] == "." or t[6] == ".":
                next
          else:
            pos1 = t[5].split(":")
            pos2 = t[6].split(":")
            cpxchr,cpxsta = pos1[0],int(pos1[1])
            cpxend = int(pos2[1])
            tupmut =  (cpxchr,cpxsta,cpxend,t[7],t[8],AF,"MUT")
            if tupmut in movelst:
                    prtlstrm.append("RAW"+","+line)
            else:
                next

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
        -tp type,    --type    type    TMB.csv or MUT.csv
        -cp  pathtocpx,  --cpx  cpx.csv     CPX.csv
        -mp  pathtomut,  --mut  mut.csv     MUT.csv or TMB.csv
        -rp  pathtoraw,  --mut_raw  mut_raw.csv    mut_raw.csv
        -lp  pathtolog,  --logpath  pathtolog    Full path to the log file [default:./]
        -op  pathtores,  --respath  pathtores    Full path to the output directory [default:./]  
    :Example:
    ..code-block::sh
        python CPXDetection_hs37d5.py -id sample -type MUT_ori/TMB_ori -cp /path/CPX.csv \
        -mp /path/MUT.csv -rp /path/mut_raw.csv -lp /path/log \
        -o /path/RES.csv
    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-v","--verbose",dest="verbose",action="store_true",\
                        help="set verbosity level [default:%(default)s]")
    parser.add_argument("-id","--caseid",help="id of the case to be analyzed",required=True)
    parser.add_argument("-typ","--typ",help="Type of TMB_ori or MUT_ori",required=True)
    parser.add_argument("-cp","--cpx",help="File of CPX.csv",required=True)
    parser.add_argument("-mp","--mut",help="File of MUT.csv or TMB.csv",required=True)
    parser.add_argument("-rp","--mut_raw",help="File of mut_raw.csv",required=True)
    parser.add_argument("-lp","--logpath",help="Path to the log file",default="./")
    parser.add_argument("-op","--respath",help="Path to the output directory",default="./")
    args = parser.parse_args()
    TMBState = args.typ
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
    if len(sys.argv) <1:
        print("Command is not right!")
    elif os.path.exists(args.mut) == False:
        if TMBState == "TMB_ori":
           with open(args.respath+args.caseid+".TMB.csv","w") as fw:
                HEAD = ",".join(header)
                fw.write(HEAD+","+"TmbTag"+"\n")
        else:
           with open(args.respath+args.caseid+".MUT.csv","w") as fw:
                HEAD = ",".join(header)
                fw.write(HEAD+"\n")
        with open(args.logpath+args.caseid+".log","w") as fw1:
                HEAD = ",".join(header)
                fw1.write(HEAD+"\n")
    elif os.path.exists(args.cpx) == False or os.path.exists(args.mut_raw) == False:
           if TMBState == "TMB_ori":
              with open(args.respath+args.caseid+".TMB.csv","w") as fw:
                  HEAD = ",".join(header)
                  fw.write(HEAD+","+"TmbTag"+"\n")
                  with open(args.mut) as fh:
                     for i in fh:
                        i = i.rstrip()
                        fw.write(i+",."+"\n")
           else:
              with open(args.respath+args.caseid+".MUT.csv","w") as fw:
                HEAD = ",".join(header)
                fw.write(HEAD+"\n")
                with open(args.mut) as fh:
                    for i in fh:
                        i = i.rstrip()
                        fw.write(i+"\n")
           with open(args.logpath+args.caseid+".log","w") as fw1:
                HEAD = ",".join(header)
                fw1.write(HEAD+","+"TmbTag"+"\n")
    else:
        diclst,diccpx,dicend,cpxlst = importCPX(args.cpx)
        mutlst = []
        cpxlst_tmp = []
        MUTLST_tmp = importMut(args.mut,mutlst,cpxlst,"MUT")
        MUTLST = importMut(args.mut_raw,MUTLST_tmp,cpxlst_tmp,"RAW")
        dicchr = sortMutchr(MUTLST)
        SORTMUT = sortMutpos(dicchr)
        dicSubMut = SubMut(SORTMUT,diccpx,dicend)
        SubMerge(dicSubMut,diccpx)
        FilterMut(args.mut,movelst,TMBState)
        importCpx2(args.cpx,TMBState)
        importRaw(args.mut_raw)
        if TMBState == "TMB_ori":
          with open(os.path.join(args.respath,args.caseid+".TMB.csv"),"w") as fw:
             HEAD = ",".join(header)
             fw.write(HEAD+","+"TmbTag"+"\n")
             for i in prtlst:
                fw.write(i+"\n")
             for i in prtcpxlst:
                val = i+",."
                if val in prtlst:
                   next
                else:
                   fw.write(i+",."+"\n")
        else:
          with open(os.path.join(args.respath,args.caseid+".MUT.csv"),"w") as fw:
             HEAD = ",".join(header)
             fw.write(HEAD+"\n")
             for i in prtlst:
                fw.write(i+"\n")
             for i in prtcpxlst:
                if i in prtlst:
                   next
                else:
                   fw.write(i+"\n")
        with open(os.path.join(args.logpath,args.caseid+".log"),"w") as fw1:
            for n in prtlstrm:
                fw1.write(n+"\n")
if __name__ == '__main__':
    main()


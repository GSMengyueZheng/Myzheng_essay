#!/usr/bin/env python
#-*-coding:utf-8

import sys
import os
import pyfaidx
reference = pyfaidx.Fasta(r"D:\Zhengmy\software\spyder_py3\python_test\wuyong_test1\clinical_QC_revise\ComplexMutation\Test2\hs37d5.fa")
path = r"D:\Zhengmy\software\spyder_py3\python_test\wuyong_test1\PolProject\MutationSignature\克隆性造血项目"
os.chdir(path)

###############################################
###########1. 统计三联核苷酸的，考虑互补，不区分癌肿

#############从reference中提取前后碱基
with open("All_noindel.txt") as fh:
  with open("All_noindel_2.txt","w") as fw:
    for i in fh:
        i = i.rstrip()
        line = i
        if i.startswith("Chr"):
            print(line+"\t"+"REF"+"\t"+"ALT",file=fw)
            continue
        t = i.split("\t")
        pos1,pos2,ref,alt = t[0],t[1],t[2],t[3]
        CHR = pos2.split(":")[0].replace("chr","")
        STA = int(pos1.split(":")[1])
        ref1 = str(reference[CHR][STA-2:STA-1])
        ref2 = str(reference[CHR][STA:STA+1])
        REF = ref1+ref+ref2
        ALT = ref1+alt+ref2
        print(line+"\t"+REF+"\t"+ALT,file=fw)       


################统计个数  
a = [x+y1+z + "->"+ x+y2+z for x in "ATCG" for y1 in "CT" for z in "ATCG" for y2 in "ATCG" if y1!=y2]
b = [x+y1+z + "->"+ x+y2+z for x in "TAGC" for y1 in "GA" for z in "TAGC" for y2 in "TAGC" if y1!=y2]
l = len(a)

dic2 = {}
dicVal = {}
for i in range(l):
    dic2[a[i]] = b[i]

for val in a:
         key = val
         dicVal[key] = 0

with open("All_noindel_2.txt") as fh:
     for i in fh:
        i = i.rstrip()
        if i.startswith("id"):
            continue
        t = i.split("\t")
        ref,alt = t[-2],t[-1]
        ALT = ref+"->"+alt
        for n in dicVal:
            if ALT == n or ALT == dic2[n]:
                dicVal[n] += 1
                    
with open("All_noindel.txt.sta.txt","w") as fw:
     for i in dicVal:
         print(i+"\t"+str(dicVal[i]),file=fw)




######################################################
##############################
########2. 统计单个核苷酸变化情况
def Subtype(Ref,Alt,lst):
    if Ref == "C" and Alt == "A":
        lst[0] += 1
    elif Ref == "G" and Alt == "T":
        lst[0] += 1
    elif Ref == "C" and Alt == "G":
        lst[1] += 1
    elif Ref == "G" and Alt == "C":
        lst[1] += 1
    elif Ref == "C" and Alt == "T":
        lst[2] += 1
    elif Ref == "G" and Alt == "A":
        lst[2] += 1
    elif Ref == "T" and Alt == "A":
        lst[3] += 1
    elif Ref == "A" and Alt == "T":
        lst[3] += 1
    elif Ref == "T" and Alt == "C":
        lst[4] += 1
    elif Ref == "A" and Alt == "G":
        lst[4] += 1 
    elif Ref == "T" and Alt == "G":
        lst[5] += 1
    elif Ref == "A" and Alt == "C":
        lst[5] += 1
    else:
        next
    return lst
 
    
with open() as fh:
    dic = {}
    dic2 = {}
    for i in fh:
        i = i.rstrip()
        if i.startswith("id"):
            continue
        t = i.split("\t")
        Sex,Cancer,Ref,Alt = t[3],t[6],t[16],t[17]
        if Cancer in dic:
                dic2[Cancer] = Subtype(Ref,Alt,dic[Cancer])
                dic[Cancer] = dic2[Cancer]
        else:
                dic[Cancer] = [0,0,0,0,0,0]
                dic2[Cancer] = Subtype(Ref,Alt,dic[Cancer])
                dic[Cancer] = dic2[Cancer]


with open() as fw:               
  for i in dic:
    print(str(i)+ "\t"+str(dic[i]),file = fw)


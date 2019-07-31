#!/usr/bin/env python
#-*- coding;utf-8 -*-
import sys
with open(sys.argv[1]) as fh:
    for i in fh:
        i = i.rstrip()
        if i.startswith("#"):
            continue
        t = i.split()
        CHR,POS,REF,ALT = t[0],int(t[1]),t[3],[4]
        REFL = len(REF)
        ALTL = len(ALT)
        if REFL == 1 and ALTL == 1:
            print(CHR+"\t"+str(POS)+"\t"+str(POS))
        else:
            print(CHR+"\t"+str(POS-1)+"\t"+str(POS+1))

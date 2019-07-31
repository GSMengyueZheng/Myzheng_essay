#!/usr/bin/env python
# File Name: WGSporadic_Main.py
# Author:
# Created Time: Fri 01 Sep 2017 10:14:21 AM CST
# Copyright: GENESEEQ Technology Inc. All rights reserved.
# Description:
#########################################################################
import sys
import os


if len(sys.argv) < 3:
    print "%s config id.dic" % os.path.basename(__file__)
    exit(1)

os.system('mkdir -p TMP HQData OUT/CovDep OUT/InsertSize OUT/VCF run_scripts Infos_TMP Data_TMP/raw')

for i in open(sys.argv[2]):
    i = i.rstrip()
    t = i.split()

    # Rename and Comnine Data
    os.system("WGSporadic_Run_101_RenameCombineData.py %s %s \
            > run_scripts/%s.sh" % (t[0], t[1], t[1]))
    # Mapping and Infos
    os.system("/GPFS02/zhengmy/WGS_Pip/Version1/Pip/WGSporadic_Run_102_RawData2BAM_2019.py %s %s \
            >> run_scripts/%s.sh" % (sys.argv[1], t[1], t[1]))
    os.system("/GPFS02/zhengmy/WGS_Pip/Version1/Pip/WGSporadic_Run_202_BAM2vcf_2019.py %s %s \
            >> run_scripts/%s.sh" % (sys.argv[1], t[1], t[1]))

    os.system("WGSporadic_S3_Single.py %s %s \
            >> run_scripts/%s.sh" % (sys.argv[1], t[1], t[1]))

    print "sh run_scripts/%s.sh >> run_scripts/%s.sh.log \
            2>> run_scripts/%s.sh.log"%(t[1], t[1], t[1])

# Myzheng_essay
Some essays were written by ZhengMy, and we hope one can give support for you!

## 1 ComplexMutation Detection
```python
python /GPFS02/zhengmy/code/CPXDetection_hs37d5.py -h
usage: CPXDetection_hs37d5.py [-h] [-v] -id CASEID [-cp CPXPATH] [-mp MUTPATH]
                              [-rp RAWPATH] [-lp LOGPATH] [-op RESPATH]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         set verbosity level [default:False]
  -id CASEID, --caseid CASEID
                        id of the case to be analyzed
  -cp CPXPATH, --cpxpath CPXPATH
                        Full path to the CPX.csv
  -mp MUTPATH, --mutpath MUTPATH
                        Full path to the MUT.csv
  -rp RAWPATH, --rawpath RAWPATH
                        Full path to the mut_raw.csv
  -lp LOGPATH, --logpath LOGPATH
                        Full path to the log file
  -op RESPATH, --respath RESPATH
                        Full path to the output directory
                        
example:
path:  /GPFS02/zhengmy/Clincal_Pipeline/ComplexMutation/Test3/Try2
eg:  python /GPFS02/zhengmy/code/CPXDetection_hs37d5.py -id ct-B190321194866-Y189-CLN-KY295-T -lp OUT/ -op OUT/

```

## 2 WGSporadic_2019 for Rare Genetic Disease
```
example: python WGSporadic_Main_2019.py PUM_2019.WGS id.dic > dic.sh
         bb -c 9 dic.sh -o big
```

## 3 Mutation Signature Plot
### signature PLOT and trinucleotide PLOT and nucleotide PLOT
```
python Calculate.py All_noindel.txt All_noindel.out1 All_noindel.out2 All_noindel
```

## 4 ... 

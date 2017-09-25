import math
import numpy
import sys


skip1File = sys.argv[1]
skipkFile = sys.argv[2]

Peak_Path = sys.argv[3]
Pnts_Path = sys.argv[4]

estLst = open(skip1File, 'r').readlines()[1:]
LnLst  = map(lambda x: x[:-1].split(), estLst)
loc1 = [int(i[8]) for i in LnLst]
loc2 = [int(i[9]) for i in LnLst]
rs2  = [i[5] for i in LnLst]


inputCsv = open(skipkFile, 'r').readlines()[1:]

box = {}

for lineNum in range(len(inputCsv)):
    line = inputCsv[lineNum][:-1].split()
    if loc1[lineNum] == int(line[8]):
        end_pnt = int(line[9])
    else:
        end_pnt = 'NA'
    cMset = line[16:]
    cMset = [25*math.log((1+2*float(j))/(1-2*float(j))) for j in cMset if j != 'NA' and float(j) <= 0.48] # <0 means 0cM
    for i in range(lineNum, loc2.index(end_pnt)+1):
        ratio=(loc2[i] - loc1[i])*1./(loc2[loc2.index(end_pnt)]-loc1[lineNum])
        if rs2[i] not in box:
            if line[0] == 'NA' or 'rs' in line[0]: 
                box[rs2[i]]=[('NA',)]
            else:
                box[rs2[i]]=[(str(ratio*float(line[0])),)]
        else:
            if line[0] == 'NA' or 'rs' in line[0]: 
                box[rs2[i]][0] += ('NA',)
            else:
                box[rs2[i]][0] += (str(ratio*float(line[0])),)
        box[rs2[i]].extend([str(ratio*j) for j in cMset])




PeakFile = open(Peak_Path, 'w')
for i in range(len(rs2)):
    PeakFile.write('chr\t'+rs2[i]+'\t'+str(loc1[i])+'\t'+str(loc2[i])+'\t')
    PeakFile.write('\t'.join(box[rs2[i]][0]))
    PeakFile.write('\n')

PeakFile.close()


# PntsFile = open(Pnts_Path, 'w')
# for rs in rs2:
#     PntsFile.write('chr22\t'+rs+'\t')
#     PntsFile.write('\t'.join(box[rs][1:])+'\n')
# PntsFile.close()


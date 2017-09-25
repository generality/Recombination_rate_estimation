## pre_step4.compnts.py

import sys
import numpy as np

input_path = sys.argv[1]
output_path= sys.argv[2]




col1_4 = []
rs = []
for i in range(10):
    lineNum = 0
    with open(input_path, 'r') as f:
        for line in f:
            if i==0:
                rs.append(line[:-1].split()[4:])
                col1_4.append('\t'.join([line[:-1].split()[0],line[:-1].split()[1],line[:-1].split()[3]])+'\t')
            else:
                rs[lineNum].extend(line[:-1].split()[4:])
            lineNum += 1

outfile = open(output_path, 'w')
for i in range(len(rs)):
    tempLst = np.array([float(j) for j in rs[i] if j != 'NA'])
    outfile.write(col1_4[i])
     if len(tempLst) > 3:
         outfile.write(str(np.median(tempLst))+'\n')
     else:
         outfile.write('NA\n')


outfile.close()
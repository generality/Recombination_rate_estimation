import sys

def g2r(gmapFilePath, inwidth, binstart, binend):
    rmapLoc = [i for i in range(binstart, binend+1, inwidth)]
    gmapLst = open(gmapFilePath, 'r').readlines()
    gmapLst = map(lambda x: x[:-1].split(), gmapLst)
    gmapLoc = [int(i[2]) for i in gmapLst]

    r = 0
    g = 0
    cMLst = []
    infoLst=[]
    temp_ind = []
    while r < len(rmapLoc):
        if rmapLoc[r] - inwidth/2 >= gmapLoc[g] and rmapLoc[r] + inwidth/2 <= gmapLoc[g+1]:
            if gmapLst[g+1][3] !='NA':
                cMLst.append([rmapLoc[r], float(gmapLst[g+1][3]) * inwidth/(gmapLoc[g+1]-gmapLoc[g]), float(gmapLst[g+1][4]) * inwidth/(gmapLoc[g+1]-gmapLoc[g])])
            else:
                cMLst.append([rmapLoc[r], 'NA'])
            r += 1
            continue

        # cross
        # ========= temp_interval doesn't mean interval length
        #
        if gmapLoc[g] <= rmapLoc[r] - inwidth/2 and rmapLoc[r] - inwidth/2 <= gmapLoc[g+1]:
            temp_ind.append(g)  # start point
            g+=1
        elif gmapLoc[g] <=  rmapLoc[r] + inwidth/2 and rmapLoc[r] + inwidth/2 <= gmapLoc[g+1]:
            temp_ind.append(g)  # before end point
            temp_interval = gmapLoc[temp_ind[0]:(temp_ind[-1]+2)]
            temp_interval_cM  = [gmapLst[i+1][3] for i in range(temp_ind[0], temp_ind[-1]+1)] # n+1
            temp_interval_num = [gmapLst[i+1][4] for i in range(temp_ind[0], temp_ind[-1]+1)]
            temp_interval_num  = map(float, temp_interval_num)
            temp_interval_num[0] = temp_interval_num[0] * (temp_interval[1] -(rmapLoc[r]-inwidth/2))/(temp_interval[1]-temp_interval[0])
            temp_interval_num[-1]= temp_interval_num[-1]* (rmapLoc[r] +inwidth/2 -temp_interval[-2])/(temp_interval[-1]-temp_interval[-2])

            if 'NA' not in temp_interval_cM:
                temp_interval_cM  = map(float, temp_interval_cM)
                temp_interval_cM[0] = temp_interval_cM[0] * (temp_interval[1] -(rmapLoc[r]-inwidth/2))/(temp_interval[1]-temp_interval[0])
                temp_interval_cM[-1]= temp_interval_cM[-1]* (rmapLoc[r] +inwidth/2 -temp_interval[-2])/(temp_interval[-1]-temp_interval[-2])
   
                cMLst.append([rmapLoc[r], sum(temp_interval_cM), sum(temp_interval_num)])
            
            else:
                cMLst.append([rmapLoc[r], 'NA', sum(temp_interval_num)])
            
            r += 1
        else:
            g += 1
    return cMLst



binstart = 20417698
binend   = 44407698
inwidth  = 10000


gmapFilePath = sys.argv[1]
res = g2r(gmapFilePath, inwidth, binstart, binend)

print '\n'.join(map(lambda x: '\t*\t'.join(map(str,x)), res))


# res = {}
# for i in range(10):
#     gmapFilePath = './gmap/20170905test.F.'+str(i+1)+'.gmap'
#     res['res.'+str(i+1)] = g2r(gmapFilePath, inwidth, binstart, binend)


# for i in range(len(res['res.1'])):
#     DLst = [res['res.'+str(j+1)][i][2] for j in range(len(res))]
#     RLst = [res['res.'+str(j+1)][i][3] for j in range(len(res))]

#     maxDcol = DLst.index(max(DLst))
#     maxRcol = RLst.index(max(RLst))

#     print '\t'.join(map(str, [res['res.'+str(maxDcol+1)][i][1], res['res.'+str(maxRcol+1)][i][1], res['res.'+str(maxRcol+1)][i][0], res['res.'+str(maxDcol+1)][i][0]]))
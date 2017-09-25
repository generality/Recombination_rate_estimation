require(tmvtnorm)

args<-commandArgs(T)
InputPath=args[1]
OutputPath=args[2]

est = as.matrix(read.csv(InputPath, stringsAsFactors = F, header = F))

starow = data.frame()

for(j in 1:nrow(est)){
  print(j)
  dat = est[j,]
  b = dat[13:length(dat)]
  b[b>=0.48] = NA
  b = b[!is.na(b)]
    if(length(b) > 250){
      b = as.numeric(b)
      b[b<0] = 0
      if(quantile(b,0.95)-quantile(b, 0.05)<= 1e-4){
        res = c(quantile(b,0.5),0,0,0)               # zero as the estiamte value
      }else{
        gmmfitl = gmm.tmvnorm(as.matrix(b))
        summaryLst = summary(gmmfitl)
        if(summaryLst$coefficients[1,4]<0.05){
          ker = summaryLst$coefficients[1,1]
          cM = 25*log((1+2*ker)/(1-2*ker))
          res = c(cM,summaryLst$coefficients[2,1],summaryLst$coefficients[1,2],summaryLst$coefficients[2,2])
        }else{
          res = c(NA,NA,NA,NA)
        }
      }
    }else{
      res = c(NA,NA,NA,NA)
    }
  starow = rbind(starow, res)
  }

alldat = cbind(starow, est)

write.table(alldat, OutputPath, row.names = F, col.names = F, quote = F)




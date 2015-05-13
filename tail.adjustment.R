rm(list=ls())
### The first part of this code fits distributions to each percentile group from the WTID
### The second part adjusts the CPS

### foreign package is needed to load state data sets
library("foreign")

### This file is produced from the WTID and has the lower bound (qmin), upper bound (qmax), and average (qave) for each percentile group
data <- read.csv("income.tail.csv")

### This sets up the data objects that hold the fitted exponential distribution for each percentile group. “a” and “c” are parameters for the exponential distribution and “error” is a measure of how far off the estimates are.
cnames <- c("year","qmin","qmax","qave","a","c","error")
d90 <- as.data.frame(matrix(data = NA,nrow = length(unique(data$Year)),ncol = length(cnames)))
row.names(d90) <- data$Year
names(d90) <- cnames
d90$year <- data$Year

### This is a lame fix to prevent a problem with losing significant digits when working with big integers and exponentials.
d90$qmin <- data$top10/1e5
d90$qmax <- data$top5/1e5
d90$qave <- data$ave10/1e5

d95 <- d90
d95$qmin <- data$top5/1e5
d95$qmax <- data$top1/1e5
d95$qave <- data$ave5/1e5

### The WTID data split the top 1%-0.1% into the 1%-0.5% and the top 0.5%-0.1%. For aesthetic reasons I merged them.
d99 <- d90
d99$qmin <- data$top1/1e5
d99$qmax <- data$top.1/1e5
d99$qave <- ((.5*data$ave1+.4*data$ave.5)/.9)/1e5

d999 <- d90
d999$qmin <- data$top.1/1e5
d999$qmax <- data$top.01/1e5
d999$qave <- data$ave.1/1e5

### I assume that the “maximum” taxable income (excluding capital gains) is one billion dollars a year. The results are not very sensitive to this assumption and it was really just made for convince. For the very top, other functional forms could be used (e.g. Pareto).
d9999 <- d90
d9999$qmin <- data$top.01/1e5
d9999$qmax <- 1e4
d9999$qave <- data$ave.01/1e5

### Here I use a fixed point method to find the parameters for the exponential distribution. I found that -2*qmin seems to be a good place to start. error is the absolute error of the average for the group. In most cases, it is less than a cent.
for (y in 1:nrow(d90)) {
  q0 <- d90$qmin[y]
  q1 <- d90$qmax[y]
  qm <- d90$qave[y]
  
  c <- -2*q0
  for (i in 1:500) {
    a = c/(exp(c*q1) - exp(c*q0))
    cnext = a*(exp(c*q1)*(c*q1 - 1) - exp(c*q0)*(c*q0 - 1))/(c*qm)
    c <- cnext
  }
  d90$a[y] <- a
  d90$c[y] <- c
  d90$error[y] <- abs(a*(exp(c*q1)*(c*q1 - 1)/c^2 - exp(c*q0)*(c*q0 - 1)/c^2) - qm)
}

for (y in 1:nrow(d95)) {
  q0 <- d95$qmin[y]
  q1 <- d95$qmax[y]
  qm <- d95$qave[y]
  
  c <- -2*q0
  for (i in 1:1000) {
    a = c/(exp(c*q1) - exp(c*q0))
    cnext = a*(exp(c*q1)*(c*q1 - 1) - exp(c*q0)*(c*q0 - 1))/(c*qm)
    c <- cnext
  }
  d95$a[y] <- a
  d95$c[y] <- c
  d95$error[y] <- abs(a*(exp(c*q1)*(c*q1 - 1)/c^2 - exp(c*q0)*(c*q0 - 1)/c^2) - qm)
}

for (y in 1:nrow(d99)) {
  q0 <- d99$qmin[y]
  q1 <- d99$qmax[y]
  qm <- d99$qave[y]
  
  c <- -2*q0
  for (i in 1:1000) {
    a = c/(exp(c*q1) - exp(c*q0))
    cnext = a*(exp(c*q1)*(c*q1 - 1) - exp(c*q0)*(c*q0 - 1))/(c*qm)
    c <- cnext
  }
  d99$a[y] <- a
  d99$c[y] <- c
  d99$error[y] <- abs(a*(exp(c*q1)*(c*q1 - 1)/c^2 - exp(c*q0)*(c*q0 - 1)/c^2) - qm)
}

for (y in 1:nrow(d999)) {
  q0 <- d999$qmin[y]
  qm <- d999$qave[y]
  q1 <- d999$qmax[y]
  
  c <- -.01
  for (i in 1:1000) {
    a = c/(exp(c*q1) - exp(c*q0))
    cnext = a*(exp(c*q1)*(c*q1 - 1) - exp(c*q0)*(c*q0 - 1))/(c*qm)
    c <- .5*cnext+.5*c
  }
  d999$a[y] <- a
  d999$c[y] <- c
  d999$error[y] <- abs(a*(exp(c*q1)*(c*q1 - 1)/c^2 - exp(c*q0)*(c*q0 - 1)/c^2) - qm)
}

for (y in 1:nrow(d9999)) {
  q0 <- d9999$qmin[y]
  qm <- d9999$qave[y]
  q1 <- d9999$qmax[y]
  
  c <- -.01
  for (i in 1:1000) {
    a = c/(exp(c*q1) - exp(c*q0))
    cnext = a*(exp(c*q1)*(c*q1 - 1) - exp(c*q0)*(c*q0 - 1))/(c*qm)
    c <- .5*cnext+.5*c
  }
  d9999$a[y] <- a
  d9999$c[y] <- c
  d9999$error[y] <- abs(a*(exp(c*q1)*(c*q1 - 1)/c^2 - exp(c*q0)*(c*q0 - 1)/c^2) - qm)
}

####################################
year.list <- as.character(c(1979:2012))

### This is a file with the CPS-U-RS.
cpi <- read.csv("CPIURS.csv")

for (year in year.list) {
### This loads the appropriate CPS. For this work I used the CEPR extracts of the CPS. I’ve used the UNICON extract to do this too—you need to change all of the variable names, likewise with the IPUMS extracts.
  print(paste("Beginning: ",year,sep=""))
  str.year <- as.character(as.numeric(year)+1)
  fname <- paste("cepr_march_",str.year,".dta",sep="")
  cps <- read.dta(fname)

### These are the variables that I use in the analysis. These are the CEPR extract names and do not all match those in IPUMS or UNICON.  
  cps <- cps[,c("id","hhseq","hjrid","famwgt","perno","incp_all","incp_ern","hrslyr","wgt","child",
                "female","married","educ","age","wkslyr","parid","spouseid","peridh","incp_wag",
                "incp_sefrm","incp_senf","incp_div","incp_ret","incp_int","incp_rnt")]
 
### Here I define taxable income for a person based on the various types of income. This includes things like wages, self employment (farm and non-farm), dividends, interest, rent, and some retirement income. 
  cps$txinc <- 0
  for(tx in c("incp_wag","incp_sefrm","incp_senf","incp_div","incp_int","incp_rnt","incp_ret")) {
    cps[!is.na(cps[,tx]),"txinc"] <- cps[!is.na(cps[,tx]),"txinc"] + cps[!is.na(cps[,tx]),tx]
  }

### These lines construct tax units in a similar fashion to Piketty & Saez as well as Burkhauser, et al. A Tax unit is any married person with his or her spouse as well as any single person age twenty or over.  
  cps$tid <- NA
  #singles
  cps$tid[cps$age >= 20 & cps$married == 0] <- cps$peridh[cps$age >= 20 & cps$married == 0]
  #married
  cps$spouseid <- as.integer(cps$spouseid)
  cps$peridh <- as.integer(cps$peridh)
  cps$tid[!is.na(cps$spouseid) & cps$married == 1] <- pmin.int(cps[!is.na(cps$spouseid) & cps$married == 1,c("spouseid")],
                                                               cps[!is.na(cps$spouseid) & cps$married == 1,c("peridh")])
  cps$spouseid <- as.character(cps$spouseid)
  cps$peridh <- as.character(cps$peridh)
  #married no spouseid--treat as single
  cps$tid[is.na(cps$spouseid) & cps$married == 1] <- cps$peridh[is.na(cps$spouseid) & cps$married == 1]
  
  #dependents
  #Chlidren 
  cps$tid[cps$age < 20 & cps$parid %in% cps$tid & !is.na(cps$parid)] <- cps$parid[cps$age < 20 & cps$parid %in% cps$tid & !is.na(cps$parid)]  
  
  if (length(cps$spouseid) == length(cps$spouseid[is.na(cps$spouseid)])) {
    #married
    cps$tid[cps$married == 1] <- as.numeric(cps$hhseq[cps$married == 1])*10+1
    cps$tid[cps$married == 0 & cps$age < 20] <- as.numeric(cps$hhseq[cps$married == 0 & cps$age < 20])*10+1
    #single
    cps$tid[cps$married == 0 & cps$age >= 20] <- as.numeric(cps$hhseq[cps$married == 0 & cps$age >= 20])*10+
      cps$perno[cps$married == 0 & cps$age >= 20]+1
  }
  
### I want to focus on prime working age people so as to avoid discussion of changes in retirement age and other factors
  cps$prime <- 0
  cps$prime[cps$age > 24 & cps$age < 55] <- 1
  
  cps$tid <- as.character(cps$tid)

### This is for the output file. Here we put together the tax unit incomes, etc.
  out.names <- c("id","w","t","prime")
  output <- as.data.frame(matrix(data = NA,nrow = length(unique(cps$tid[!is.na(cps$tid)])),ncol = length(out.names)))
  names(output) <- out.names
  
### Here I build the output file. “tid” is the tax unit id, “t” is the taxable income, “w” is the weight for the tax unit (there are a handful of tax units where mean != unique and so I chose to just average them all.), and “prime” is a marker as to whether or not there is a prime working age person in the tax unit.
  output$id <- sapply(split(cps[!is.na(cps$tid),"tid"],cps$tid[!is.na(cps$tid)]),function(z) unique(z))
  output$t <- sapply(split(cps$txinc[!is.na(cps$tid)],cps$tid[!is.na(cps$tid)]),function(z) sum(z))
  output$w <-  sapply(split(cps$famwgt[!is.na(cps$tid)],cps$tid[!is.na(cps$tid)]),function(z) mean(z))
  output$prime <-  sapply(split(cps$prime[!is.na(cps$tid)],cps$tid[!is.na(cps$tid)]),function(z) max(z))

### This is an adjustment for inflation  
  inflation <- cpi$cpi[cpi$year == 2012]/cpi$cpi[cpi$year == as.numeric(year)]
  output$t <- output$t*inflation

### Cutting out the misfits. Most of these are people under 20 that are not tied to adults. Let me know if there are better ways of dealing with this population.
  output <- output[!is.na(output$id),]
  
  print(paste("Beginning adjustment: ",year,sep=""))
  #Tail Adjustments
### Ordering the tax units by taxable income
  output <- output[order(output$t),]
  output$trnk <- cumsum(output$w)
  output$trnk <- 100*output$trnk/max(output$trnk)
  
### Adjust the taxable income based on the formulas discussed in the technical appendix. “t.alt” is the adjusted taxable income.
  #Adjust people between 90% and 95%
  output$t.alt <- output$t
  output$fw <- output$w
  samp <- nrow(output[output$trnk > 90 & output$trnk < 95,])
  q0 <- d90[as.character(year),"qmin"]
  a <- d90[as.character(year),"a"]
  c <- d90[as.character(year),"c"]
  r.samp <- (1:samp)/samp-1e-8
  output$t.alt[output$trnk > 90 & output$trnk < 95] <- 1e5*sort(log(r.samp*c/a + exp(c*q0))/c)
  
  #Adjust people between 95% and 99%
  samp <- nrow(output[output$trnk > 95 & output$trnk < 99,])
  q0 <- d95[as.character(year),"qmin"]
  a <- d95[as.character(year),"a"]
  c <- d95[as.character(year),"c"]
  r.samp <- (1:samp)/samp-1e-8
  output$t.alt[output$trnk > 95 & output$trnk < 99] <- 1e5*sort(log(r.samp*c/a + exp(c*q0))/c)
  
  #Duplicates to get the top 1% to 0.1%
### To keep enough records to do reasonable analysis, I duplicate records in the top 1%. This lets me keep enough records at the top of the distribution for additional analysis. The weights are scaled down to preserve the overall population size.
  dups <- output[output$trnk > 99,]
  dups$fw <- dups$w*0.09
  dups$id <- as.character(as.numeric(dups$id) + 1)
  samp <- nrow(dups)
  q0 <- d999[as.character(year),"qmin"]
  a <- d999[as.character(year),"a"]
  c <- d999[as.character(year),"c"]
  r.samp <- (1:samp)/samp-1e-8
  dups$t.alt <- 1e5*sort(log(r.samp*c/a + exp(c*q0))/c)
  
  #Duplicates to get the top 1% to 0.01%
  dups2 <- output[output$trnk > 99,]
  dups2$fw <- dups2$w*0.01
  dups2$id <- as.character(as.numeric(dups2$id) + 2)
  samp <- nrow(dups2)
  q0 <- d9999[as.character(year),"qmin"]
  a <- d9999[as.character(year),"a"]
  c <- d9999[as.character(year),"c"]
  r.samp <- (1:samp)/samp-1e-8
  dups2$t.alt <- 1e5*sort(log(r.samp*c/a + exp(c*q0))/c)
  
  #Adjust people between 99% and 99.9%
  samp <- nrow(output[output$trnk > 99,])
  q0 <- d99[as.character(year),"qmin"]
  a <- d99[as.character(year),"a"]
  c <- d99[as.character(year),"c"]
  r.samp <- (1:samp)/samp-1e-8
  output$t.alt[output$trnk > 99] <- 1e5*sort(log(r.samp*c/a + exp(c*q0))/c)
  #Scale down weights for the top 1% because we will be adding the duplicates for the top 0.1% and 0.01% later
  output$fw[output$trnk > 99] <- output$w[output$trnk > 99]*.9
  
  #Add the duplicated groups to the rest of the population
  output <- rbind(output,dups,dups2)
  output <- output[order(output$t.alt),]
  output$frnk <- cumsum(output$fw)
  output$frnk <- 100*output$frnk/max(output$frnk)
  write.csv(output,file=paste(“adjusted.cps.”,year,”.csv”,sep=“”))
}

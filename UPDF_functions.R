require(compiler)
require(msm)
require(snowfall)
require(parallel)
require(doParallel)

## This script contains two functions. The first "updf" calculates the uniform probability density function for the
## supplied data and the second "updf.plot" plots the results.

# Multicore setup
cl <- makeCluster(detectCores())
registerDoParallel(cl)

###############################################################
###############################################################
###############################################################

updf <- function(site,cer.type,ct,start,end,chron,interval=5,cutoff=0.1,min.period=25) {
  
  ## Remove any type with a count of 0 from all vectors
  trim <- which(ct>0)
  cer.type <- as.vector(cer.type[trim])
  site <- as.vector(site[1])
  start <- start[trim]
  end <- end[trim]
  ct <- ct[trim]
  chron <- chron[trim]
  
  
  ## define function for rounding by base defined in "interval" argument and round start and end dates for each type
  mround <- function(x,base){base*round(x/base)}
  start <- mround(start,interval)
  end <- mround(end,interval)
  
  ## find minimal ceramic intervals based on rounded type date overlaps
  years <- unique(sort(c(start,end)))
  period <- do.call(rbind,foreach(i=1:(length(years)-1)) %dopar% rbind(c(years[i],years[i+1]))) 
  # Ensure that no period is shorter than the specificed min.period argument and combine periods to satisfy this requirement
  if (length(period)>2) {
    per.short <- foreach(i=1:(nrow(period)-1)) %dopar% (period[i+1,1]-period[i,1])
    comb.per <- which(per.short<min.period)
    if (length(comb.per)>0){ period <- period[-comb.per,]
    for (i in 1:(nrow(period)-1)) {period[i,2] <- period[i+1,1]}}
  }  
  
  ## Define lists of ceramic period lengths and site period lengths
  cer.per <- foreach(i=1:length(start)) %dopar% (seq((start[i]:end[i]))+(start[i]-1))
  per.len <- foreach(i=1:nrow(period)) %dopar% (seq(period[i,1]:period[i,2])+(period[i,1])-1)
  
  ## Create Uniform Distance dataframe via parallel function
  ucalc <- function(a,b) {
    out <- (length(intersect(a,b))-1)/(length(a)-1)
    if (out<0) {out <- 0}
    return(out)}
  
  udist <- foreach(a=cer.per,.combine='rbind') %:% foreach(b=per.len,.combine='c') %dopar% (ucalc(a,b))
  if (length(udist)==1) {udist <- as.matrix(udist)}
  colnames(udist) <- paste('AD',period[,1],'-',period[,2],sep='')
  
  ## Create dataframe of prior values and sum to calculate prior probabilities
  prior.tab <- udist*ct
  if (ncol(prior.tab)>1) {
    if (length(prior.tab[which(chron==1),])==ncol(prior.tab)) {prior <- prior.tab[which(chron==1),]/sum(prior.tab[which(chron==1),])} else
    {prior <- colSums(prior.tab[which(chron==1),])/sum(prior.tab[which(chron==1),])}} else {
      prior <- udist[1,]
    }
  prior[is.na(prior)] <- 0
  
  ## Create dataframe of count probabilities
  pij <- sweep(prior.tab,2,colSums(prior.tab),'/')
  pij[is.na(pij)] <- 0
  
  ## Create dataframe of probabilities based on uniform distribution
  uij <- sweep(udist,2,colSums(udist),'/')
  uij[is.na(uij)] <- 0
  
  ## Create dataframe of standard deviations of uniform distribution probabilities
  sd.t <- sqrt(sweep((uij*(1-uij)),2,colSums(ceiling(udist)),'/'))
  
  ## Create dataframe of conditionals and calculate conditional probabilities
  cij <- dnorm(pij,uij,sd.t)
  cij[which(uij-sd.t==1)] <- dnorm(1,1,1) # For intervals with only one ceramic type, set to standard normal
  is.na(cij) <- !is.finite(cij) 
  
  # calculate conditional probability and rescale from 0-1
  conditional <- apply(cij,2,mean,na.rm=T)
  conditional[is.na(conditional)] <- 0
  if (sum(conditional)>0) {conditional <- conditional/sum(conditional)}
  
  ## Calculate posterior proportions and remove any generated NA
  posterior <- foreach(i=1:length(prior),.combine='c') %dopar% ((prior[i]*conditional[i]))
  posterior <- posterior/sum(posterior)
  posterior[is.na(posterior)] <- 0
  
  ## Deal with edge cases where sum of posterior probabilities = 0 or conditional = prior
  if(sum(posterior)==0) {posterior <- prior}
  if (identical(conditional,prior)) {posterior <- prior}
  
  # calculated beginning (lwr) and ending (upr) dates based on the first and last period with posterior probability about the selected threshold
  lwr <- period[min(which(posterior>cutoff*max(posterior))),1]
  upr <- period[max(which(posterior>cutoff*max(posterior))),2]
  
    ## Create output list and return
  out.list <- list(site=site,prior=prior,posterior=posterior,conditional=conditional,period=period,samp.size=sum(ct),occ=cbind(lwr,upr))
  return(out.list)
}



##############################################################################
##############################################################################
##############################################################################

## Function for plotting the output of the updf function
## Plots prior, conditional, and posterior probabilities by period

updf.plot <- function(updf.out, beg.date, end.date, plot.interval) {
  prior <- updf.out$prior
  conditional <- updf.out$conditional
  posterior <- updf.out$posterior
  samp.size <- updf.out$samp.size
  period <- updf.out$period
  site <- updf.out$site
  occ <- updf.out$occ
  
  g.per <- as.vector(t(period))
  g.per <- c(g.per[1],g.per,g.per[length(g.per)])
  g.post <- as.vector(t(cbind(posterior,posterior)))
  g.post <- c(0,g.post,0)
  g.prior <- as.vector(t(cbind(prior,prior)))
  g.prior <- c(0,g.prior,0)
  g.con <- as.vector(t(cbind(conditional,conditional)))
  g.con <- c(0,g.con,0)
  
  y.lim <- ceiling(max(c(prior,conditional,posterior))*10)/10
  plot(g.per,g.post,type='n',ylim=c(0,y.lim),xlim=c(beg.date, end.date),ylab='Probability of Occupation',xlab='Date',xaxt='n',main=paste(site," N=",samp.size,', AD',occ[1],'-',occ[2],sep=''))
  axis(1, at = seq(beg.date, end.date, by = plot.interval), las=2)
  polygon(g.per,g.post,col='lightgray')
  lines(g.per,g.prior,type='l',col='red',lwd=2)
  lines(g.per,g.con,lty=2,lwd=2,col='blue')
  abline(v=occ[1],lty=2)
  abline(v=occ[2],lty=2)
  legend('top',bty='n',cex=0.75,inset=-0.07,ncol=3,c('Prior','Conditional','Posterior'),col=c('red','blue','lightgray'),lty=c(1,2,1),lwd=c(2,2,5),xpd=T)
}

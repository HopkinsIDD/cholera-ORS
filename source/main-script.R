#########################
## Plots for ORS Paper ##
#########################
library(fields) # for some image.plot stuff
library(MASS) # for neg. binom model (but not using at the moment)
library(TTR) # for time series smoothing
library(magrittr)
library(dplyr)
library(irr)
library(lme4)
source("Code/R/bangpharma-functions.R")

remake_pdfs <- FALSE

## ----------------- ##
## Bring in key data ##
## ----------------- ##

## dates
date.seq <- seq.Date(as.Date("2013-04-17"),as.Date("2013-10-30"),by=1)

## cholera
trimmed.cholera.daily <- readRDS("Generated_Data/trimmed_cholera_daily.rds")
trimmed.cholera.weekly<-readRDS("Generated_Data/trimmed_cholera_weekly.rds")
## diarrhea
trimmed.all.diarrhea.daily <- readRDS("Generated_Data/trimmed_all_diarrhea_daily.rds")
trimmed.nonc.diarrhea.daily <- readRDS("Generated_Data/trimmed_nonc_diarrhea_daily.rds")
trimmed.all.diarrhea.weekly <- readRDS("Generated_Data/trimmed_all_diarrhea_weekly.rds")
trimmed.diarrhea.weekly <- readRDS("Generated_Data/trimmed_nonc_diarrhea_weekly.rds")

## ors
packets.indmat.daily <- readRDS("Generated_Data/packets_indmat_daily.rds")
packets.daily <- readRDS("Generated_Data/packets_daily.rds")
cust.daily <- readRDS("Generated_Data/cust_daily.rds")
cust.indmat.daily<-readRDS("Generated_Data/cust_indmat_daily.rds")
cust.indmat.weekly<-readRDS("Generated_Data/cust_indmat_weekly.rds")
packets.indmat.weekly<-readRDS("Generated_Data/packets_indmat_weekly.rds")
packets.weekly<-readRDS("Generated_Data/packets_weekly.rds")
cust.weekly<-readRDS("Generated_Data/cust_weekly.rds")
cust.weekly[is.nan(cust.weekly)] <- NA # NaN's are a result of Missing data

## ----------------- ##
## ORS summary stats ##
## ----------------- ##

## how many packaets and customers
sum(cust.daily,na.rm=T)
sum(packets.daily,na.rm=T)
## range of daily customers
range(cust.daily,na.rm=T)
quantile(cust.daily,na.rm=T)
range(packets.daily,na.rm=T)

tmp <- c(packets.daily/cust.daily )
range(ifelse(is.nan(tmp),NA,tmp),na.rm=T)
mean(ifelse(is.nan(tmp),NA,tmp),na.rm=T)
quantile(ifelse(is.nan(tmp),NA,tmp),na.rm=T)
quantile(cust.daily,na.rm=T)

## per pharmacy mean
plot(ecdf(apply(packets.daily,1,function(x) mean(x,na.rm=T))))
plot(ecdf(apply(cust.daily,1,function(x) mean(x,na.rm=T))))
#hist()
range(apply(cust.daily,1,function(x) sum(x==0,na.rm=T)/sum(!is.na(x))))

## missing data days
## note that we don't want to count 'missing days' before the
## pharmacy was enrolled

## let's get the first non-NA for each row
tail(cust.daily,35)[,1:10]
first.days <- apply(cust.daily,1,function(x) min(which(!is.na(x))))
cust.daily[-c(cbind(1:length(first.days),c(first.days)))]

## missing days
missing.pharm.days <- sapply(1:nrow(cust.daily),function(x) sum(is.na(cust.daily[x,first.days[x]:ncol(cust.daily)])))
quantile(missing.pharm.days)
mean(missing.pharm.days)

##########################
## load bangpharma data ##
##########################

bp.raw <- read.csv("Data/bangpharmaData/bangpharma_sales_2014-04-30.csv")
## do a little pre-processing of the data
bp <- clean.new.data(bp.raw,TRUE)
combined.ors <- ddply(bp,.(sale.created_at),function(x) sum(x$sale.ors))
colnames(combined.ors) <- c("date","ors")
bp.red <- bp[bp$sale.created_at %within% interval(min(date.seq),max(date.seq)),]
out <- ddply(bp.red,"sale.pharmacy_id",function(x) {
    merge(x[,c("sale.created_at","sale.ors")],data.frame(sale.created_at=date.seq),
          all.x=F,all.y=T)
},.drop=FALSE)
out <- dcast(out,sale.created_at~sale.pharmacy_id,fun.aggregate = sum)
rownames(out) <- out[,1]
out <- out[,-1]


## let's look by week
wn <- week(as.Date(rownames(out),format="%Y-%m-%d"))
out <- data.frame(out)
out$week <- wn
out.by.week <- ddply(out,'week',function(x){
    tmp <- x[,-c(ncol(x))]
    rc <- colSums(tmp,na.rm=T)
})


## first let's check out reporting volume by pharmacy overtime
## by week

pp <- bp %>% group_by(sale.pharmacy_id,week) %>% summarise(reports=n()) %>%
  ggplot() + geom_point(aes(x=week,y=reports,group=1)) +
  geom_smooth(aes(x=week,y=reports,group=1),method=lm) +
  facet_wrap('sale.pharmacy_id',scales = 'free') +xlab('week') + ylab('number of sales reports to system')
pp + theme(axis.text.x = element_text(size = 6),axis.ticks.length = unit(.1, "cm")) +
  theme(axis.text.y = element_text(size = 6),axis.ticks.length = unit(.1, "cm")) +
  scale_x_continuous(breaks=pretty_breaks(n=5)) -> pp
if(remake_pdfs){
    to.pdf(print(pp),filename="Figures/bp_reports_overtime_bypharm.pdf",
           width=15,height=10)
}

mod.dat <- bp %>% group_by(sale.pharmacy_id,week) %>% summarise(reports=n())
colnames(mod.dat) <- c('pid','week','n')
fit <- lmer(n ~ factor(week) + (1|pid),data=mod.dat)

## let's first compare weeks to weeks
paper.week <- t(cust.weekly)
bp.week <-out.by.week[,-1]

## comparing all pharmacies each week
week.cor <- sapply(1:29,function(x) cor(unlist(paper.week[x,]),unlist(bp.week[x,]),method="pearson",use="pairwise.complete.obs"))
## compariing individual pharamcies to themselves
week.cor2 <- sapply(1:50,function(x) cor(unlist(paper.week[,x]),unlist(bp.week[,x]),method="pearson",use="pairwise.complete.obs"))

pdf("Figures/community_weekly_correlation.pdf")
hist(week.cor,breaks="fd",col="grey",border="white",main="",xlab="Pearson Correlation")
abline(v=median(week.cor),col="orange",lwd=2)
text(0.6,9,sprintf("Median (comunity-wide) \n Correlation= %.2f",median(week.cor)))
dev.off()

pdf("Figures/withinpharmacy_weekly_correlation.pdf")
hist(week.cor2,breaks="fd",col="grey",border="white",main="",xlab="Pearson Correlation")
abline(v=median(week.cor2),col="orange",lwd=2)
text(0.4,15,sprintf("Median (within pharmacy) \n Correlation= %.2f",median(week.cor2)))
dev.off()

## now look at community-wide ORS per week correlation
## not so impressive
cor(rowSums(paper.week,na.rm=T),rowSums(bp.week,na.rm=T),method="pearson")
cbind(diff(rowSums(paper.week,na.rm=T))/mean(rowSums(paper.week,na.rm=T)),
    diff(rowSums(bp.week,na.rm=T))/mean(rowSums(bp.week,na.rm=T)),
    method="pearson")

## how well correlated are the sales from bangpharma and paper?
daily.sales <- t(out[,-ncol(out)])
quantile(
    sapply(1:50,function(x)
        cor(unlist(daily.sales[x,]),unlist(cust.daily[x,]),use="pairwise.complete.obs"))
)

median(sapply(1:197,function(x)
            cor(daily.sales[,x],cust.daily[,x],use="pairwise.complete.obs")))


###########################################################
## compare the cumaltive customers between bp and papers ##
###########################################################
phar.ids <- unique(bp.raw$sale.pharmacy_id)
phar.ids <- phar.ids[order(phar.ids)]

## get cumulatives
cum.bp <- apply(out,2,function(x) {x <- ifelse(is.na(x),0,x);cumsum(x)})[,-ncol(out)]
cum.papier <- apply(cust.daily,1,function(x) {x <- ifelse(is.na(x),0,x);cumsum(x)})
rownames(cum.papier) <- rownames(cum.bp)
colnames(cum.bp)<- colnames(cum.papier)
cum.bp.m <- cbind(melt(cum.bp),system="phone")
cum.papier.m <- cbind(melt(cum.papier),system="paper")

colnames(cum.bp.m) <- colnames(cum.papier.m) <- c("date","pid","customers","system")
cum.bp.m[,1] <- as.Date(cum.bp.m[,1])
cum.papier.m[,1] <- as.Date(cum.papier.m[,1])
cum.comb <- rbind(cum.bp.m,cum.papier.m)
cum.comb$pid<-factor(cum.comb$pid)
gg <- ggplot(cum.comb) +
  geom_line(aes(x=date,y=customers,color=system)) +
  facet_wrap('pid',scales = 'free') + xlab('day') + ylab('cumulative customers') +
  theme(axis.text.x =  element_text(size = rel(0)),axis.ticks.length = unit(.1, "cm"),
        axis.text.y = element_text(size = rel(.5)),axis.ticks.length = unit(.1, "cm"))
if(remake_pdfs){
    to.pdf(print(gg),"Figures/cdfs-bangpharma-vs-paper.pdf",width=8,height=6)
}
## let's use a concordance statistic to look at how well the paper and phone
## systems predict (the same) high and low sales excess sales
paper <- ifelse(apply(cust.daily,1,function(x) x-mean(x,na.rm=T)/sd(x,na.rm=T))>0,1,0)
phone <- ifelse(apply(t(out[,1:50]),1,function(x) x-mean(x,na.rm=T)/sd(x,na.rm=T))>0,1,0)
ratings <- cbind(c(paper),c(phone))
miss <- rowSums(ratings) %>% is.na %>% which
table(ratings[-miss,1],ratings[-miss,2])
mean(ratings[-miss,1]==ratings[-miss,2])
kappa2(ratings[-miss,])

## -------------------------------------- ##
## Plot of pharmacy locations in arichpur ##
## -------------------------------------- ##


## ---------------------------------------------- ##
## Plot of cholera and non-cholera diarrhea cases ##
## ---------------------------------------------- ##

my.dates <- seq.Date(as.Date("2013-04-17"),as.Date("2013-10-30"),by=1)

pdf("Figures/ORS-Fig2.pdf",height=5,width=10)
palette(brewer.pal(8,"Dark2"))
par(mfrow=c(2,1),mar=c(1,0,0,0),oma=c(2,3,2,1),mgp=c(2,.2,0),tck=-.01)
plot(trimmed.all.diarrhea.daily,type="h",col="grey",xlab="Day",ylab="Non-cholera Diarrhea",xaxt = "n")
axis(3,at=seq(0,length(my.dates),length=10),labels=format(seq.Date(as.Date("2013-04-17"),as.Date("2013-10-20"),length=10),"%d-%b"),cex.axis=.8)
lines(1:length(my.dates),SMA(trimmed.all.diarrhea.daily,n=3),col=1)
#lines(smooth.spline(trimmed.nonc.diarrhea.daily),type="l",col=1,lty=1)
fit.all <- glm(trimmed.all.diarrhea.daily~1,family="quasipoisson")
#exp(confint(glm.nb(trimmed.all.diarrhea.daily~1)))
cis.all <- exp(confint(fit.all))

#abline(h=mean(trimmed.nonc.diarrhea.daily),col=AddAlpha(2,.5))
abline(h=cis.all ,col=AddAlpha(2,.5),lty=2)
#rug(which(trimmed.all.diarrhea.daily > cis.all[2]),col=AddAlpha(2,.5),lwd=1)
mtext("all diarrhea",side=2,outer=T,line=1.5,at = .76)
mtext("cholera diarrhea",side=2,outer=T,line=1.5,at = .26)
mtext("day of clinical presentation",side=1,outer=T,line=.5,at = .5)
text(1,8,"A",cex=1.2)

plot(trimmed.cholera.daily,type="h",col="grey",xlab="Day",ylab="Cholera",xaxt="n")
axis(1,at=seq(0,length(my.dates),length=10),labels=format(seq.Date(as.Date("2013-04-17"),as.Date("2013-10-20"),length=10),"%d-%b"),cex.axis=.8)
lines(1:length(my.dates),SMA(trimmed.cholera.daily,n=7),col=1)
#lines(smooth.spline(trimmed.cholera.daily,df=25),type="l",col=1,lty=1)
fit.c <- glm(trimmed.cholera.daily~1,family="quasipoisson")
cis.c <- exp(confint(fit.c))
#normal.cis.c <- mean(trimmed.cholera.daily) + c(-1,1)*qnorm(.975)*sqrt(mean(trimmed.cholera.daily)/sum(trimmed.cholera.daily))
#abline(h=mean(trimmed.cholera.daily),col=AddAlpha(2,.5))
abline(h=cis.c ,col=AddAlpha(2,.5),lty=2)
#rug(which(trimmed.cholera.daily > cis.c[2]),col=AddAlpha(2,.5),lwd=1)
legend("topright",c("raw data","3-day moving average","95% CI for incidence rate"),col=c("grey",1,2),lty=c(1,1,2),bty="n")
text(1,4.5,"B",cex=1.2)

fit.nc <- glm(trimmed.nonc.diarrhea.daily~1,family="quasipoisson")
cis.nc <- exp(confint(fit.nc))
mean(trimmed.nonc.diarrhea.daily > cis.nc[2])
sum(trimmed.nonc.diarrhea.daily > cis.nc[2])
## plot(trimmed.nonc.diarrhea.daily,type="h",col="grey",xlab="Day",ylab="Non-cholera Diarrhea",xaxt = "n")
## axis(1,at=seq(0,length(my.dates),length=10),labels=format(seq.Date(as.Date("2013-04-17"),as.Date("2013-10-12"),length=10),"%d-%b"),cex.axis=.8)
## lines(1:length(my.dates),SMA(trimmed.nonc.diarrhea.daily,n=3),col=1)
## #lines(smooth.spline(trimmed.nonc.diarrhea.daily),type="l",col=1,lty=1)
## normal.cis.nc <- mean(trimmed.nonc.diarrhea.daily) + c(-1,1)*qnorm(.975)*sqrt(mean(trimmed.nonc.diarrhea.daily))
## #abline(h=mean(trimmed.nonc.diarrhea.daily),col=AddAlpha(2,.5))
## abline(h=normal.cis.nc ,col=AddAlpha(2,.5),lty=2)
## abline(h=0 ,col=AddAlpha(2,.5),lty=2)
## rug(which(trimmed.nonc.diarrhea.daily > normal.cis.nc[2]),col=AddAlpha(2,.5),lwd=3)
## mtext("cholera",side=2,outer=T,line=1.5,at = .76)
## mtext("non-cholera diarrhea",side=2,outer=T,line=1.5,at = .26)
## mtext("day of clinical presentation",side=1,outer=T,line=.5,at = .5)
## text(1,7,"B",cex=1.2)
dev.off()


## ---------------------------------- ##
## Cholera and Diarrhea Summary Stats ##
## ---------------------------------- ##
sum(trimmed.all.diarrhea.daily) # number of hosptalized diarrhea cases from Arichpur
sum(trimmed.cholera.daily)
sum(trimmed.nonc.diarrhea.daily)

# divisions by where we found the case
# from debug(combine.and.anon.clinic.sample.data) call:
## tongi.culture.results
## Negative Positive
##      132       28
## Browse[2]> table(icddrb.culture.results)
## icddrb.culture.results
## Negative Positive
##      114       41

# days with and without diarrhea
mean(trimmed.all.diarrhea.daily==0)
mean(trimmed.cholera.daily==0)
mean(trimmed.nonc.diarrhea.daily==0)

range(trimmed.all.diarrhea.daily)
range(trimmed.cholera.daily)


## --------------------------- ##
## Figure 3: Heat map of cases ##
## --------------------------- ##
## heat map of ORS and chlera
zmax <- 3
zmin <- -3
quartz("",height=6,width=12)

make.fig3 <- function(){
    pdf("Figures/Fig3-ORS-heatmap3.pdf",height=6,width=12)
    par(oma=c(2,2,2,2))
    palette(brewer.pal(8,"Dark2"))
    cust.daily.norm <- apply(cust.daily[,],1,
                             function(x) (x - mean(x,na.rm=T))/sd(x,na.rm=T))
    pids <- colnames(cust.daily.norm)
                                        # sort by mmulti dimensional scaling coefficient
    mds <- cmdscale(dist(t(cust.daily.norm),method="manhattan"),k=1)
    ## or by missing days
                                        #mds <- apply(cust.daily.norm,2,function(x) sum(is.na(x)))
    cust.daily.norm <- cust.daily.norm[,order(mds,decreasing = T)]
    cust.daily.norm <- ifelse(cust.daily.norm>zmax,zmax,cust.daily.norm)
    cust.daily.norm <- ifelse(cust.daily.norm<zmin,zmin,cust.daily.norm)
    cust.daily.norm <- ifelse(is.na(cust.daily.norm),NaN,cust.daily.norm)
                                        #split.screen(rbind(c(0,0.9,0,.8),c(0.9,1,0,.8),c(0,.9,.79,1)))


    split.screen(rbind(c(0,.9,0,.2),c(0,0.9,0.2,.8),c(0.9,1,0,.8),c(0,.9,.79,1)))
    screen(1)
    palette(brewer.pal(4,"RdBu"))
    par(mar=c(2,4.1,0,2.1),oma=c(2,2,1,1),mgp=c(2,.2,0),tck=-.01)
    date.grid <- seq.Date(from=as.Date("2013-04-17"),to=as.Date("2013-10-20"),by="week")
    date.big <- seq.Date(from=as.Date("2013-04-17"),to=as.Date("2013-10-20"),by="month")
    plot(date.seq,
         apply(cust.daily.norm,1,function(x) mean(x,na.rm=T)),
         type="h",axes=FALSE,
         col=ifelse(rowSums(cust.daily.norm,na.rm=T)>0,1,4),
         xlab="",ylab="",
         xlim=c(as.Date("2013-04-17")+7,as.Date("2013-10-30")-7))
    axis(3,
         at=seq.Date(as.Date("2013-04-17"),
             as.Date("2013-10-30"),length=10),
         labels=rep("",10),
         cex.axis=.8)
    axis(1,
         at=seq.Date(as.Date("2013-04-17"),
             as.Date("2013-10-30"),length=10),
         labels=format(
             seq.Date(as.Date("2013-04-17"),
                      as.Date("2013-10-30"),
                      length=10),"%d-%b"),
         cex.axis=.8)
    axis(2,
         cex.axis=.8)


    ## adding to the lims to help them line up with other [plots
    ## this doesn't truncate the data
                                        #axis(1,at=date.big,labels=format(as.Date(date.big),"%b"))
    #abline(v=date.grid,col=AddAlpha("grey",.3),lty=1)
    abline(h=seq(-40,40,by=.25),col=AddAlpha("grey",0.5),lty=1)

    screen(2)
    par(mar=c(2,4.1,0,2.1))
    image(cust.daily.norm,
          col=colorRampPalette(rev(brewer.pal(8,"RdBu")))(1000),
          axes=F,main="",zlim=c(zmin,zmax))
    mtext(text=pids[order(mds)],#,1:length(pids),
          side=2, line=0.3,
          at=seq(0,1,length=length(1:50)), las=1, cex=0.6)
    axis(1,
         at=seq(0,1,length=10),
         labels=format(
             seq.Date(as.Date("2013-04-17"),
                      as.Date("2013-10-30"),
                      length=10),"%d-%b"),
         cex.axis=.8)
    ## axis(1,at=seq(0,1,length=12),
    ##      labels=seq(1,179,length=12))
                                        #mtext(text=seq(1,144,by=14), side=1, line=0.3, at=seq(0,1,length=11), las=2, cex=0.5)
screen(3)
image.plot(cust.daily.norm, legend.only=T,
           col=colorRampPalette(rev(brewer.pal(8,"RdBu")))(1000),
           smallplot=c(.1,.2, .3,.75))
screen(4)
par(mar=c(.5,2.4,2,.5))
plot(0:196,SMA(trimmed.all.diarrhea.daily,n=3),col="black",
     lty=1,type="l",axes=FALSE,xlab="",ylab="")
lines(0:196,SMA(trimmed.cholera.daily,n=3),col="black",lty=2)
#grid()
abline(h=1:8,col=AddAlpha("grey",.1))
abline(v=seq(0,196,by=7),col=AddAlpha("grey",.1))
#axis(3,at=seq(1,179,length=12))
axis(3,
     at=seq(0,length(my.dates),length=10),
     labels=format(seq.Date(as.Date("2013-04-17"),
         as.Date("2013-10-30"),length=10),"%d-%b"),
     cex.axis=.8)
axis(2,cex.axis=.8)
    close.screen(all = TRUE)
    dev.off()
}

palette(brewer.pal(8,"Dark2"))
pharm.locs <- read.csv("Data/PharmacyRecruitement/Pharmacy_Registration_Data.csv",as.is=T)
all.pharms <- read.csv("Data/All_Coordinate_of_Pharmacies.csv",as.is=TRUE)
all.pharms$sid <- as.numeric(gsub("AP","",all.pharms[,1]))


pharm.locs <- data.frame(pharm.locs[,-ncol(pharm.locs)][,c(1,3,4)])
colnames(pharm.locs) <- c("id","lat","lon")

warning("merge here is messed up, not sure what happened but need to get this
working again or delete")
new.locs.df <- merge(pharm.locs,all.pharms,by="sid")
new.locs.df$mds <- mds
rbPal <- colorRampPalette(c('red','blue'))
new.locs.df$col <- rbPal(100)[as.numeric(cut(mds,breaks = 100))]

library(ggmap)
mymap = get_map(location = c(lon = 90.407, lat = 23.887),
    zoom = 15, maptype = 'satellite',color="bw")
pharmmap = ggmap(mymap,legend="topleft") +
geom_point(data=new.locs.df, aes(x=Longitude, y=Latitude,col=mds)) +
labs(title = "Spatial Specficity")+ scale_colour_gradient2(low="red", high="blue",midpoint=0)
pdf("Figures/spatialspec.pdf")
print(pharmmap)
dev.off()
#ggplot_build(pharmmap)


#####################
## Cholera and ORS ##
#####################

pred <- apply(cust.daily.norm,1,function(x) mean(x,na.rm=T))

## first let's look at the cross correlations
pdf("Figures/ccf-dia-andchol.pdf",width=6,height=6)
par(mfrow=c(2,1),mar=c(.5,1,.5,0),oma=c(2,2,1,1),mgp=c(2,.2,0),tck=-.01)
dia.ccf <- ccf(trimmed.all.diarrhea.daily,pred,
               ci.type="ma",
               plot=F,
               lag.max=10)
upperCI <- qnorm((1 + 0.95)/2)/sqrt(dia.ccf$n.used)
lowerCI <- -qnorm((1 + 0.95)/2)/sqrt(dia.ccf$n.used)
plot(-10:10,dia.ccf$acf,type="h",ylim=c(-.2,.3),lwd=2,xaxt="n")
axis(3)
abline(h=0)
polygon(x=c(-100,100,100,-100),y=c(lowerCI,lowerCI,upperCI,upperCI),
        border=FALSE,col=AddAlpha("grey",.5))
abline(h=upperCI,col="grey",lty=2)
abline(h=lowerCI,col="grey",lty=2)
grid()
chol.ccf <- ccf(trimmed.cholera.daily,pred,lag=10,plot=F)#rowSums(cust.daily.norm,na.rm=T))
upperCI <- qnorm((1 + 0.95)/2)/sqrt(chol.ccf$n.used)
lowerCI <- -qnorm((1 + 0.95)/2)/sqrt(chol.ccf$n.used)
plot(-10:10,chol.ccf$acf,type="h",ylim=c(-.2,.25),lwd=2)
abline(h=0)
polygon(x=c(-100,100,100,-100),y=c(lowerCI,lowerCI,upperCI,upperCI),
        border=FALSE,col=AddAlpha("grey",.5))
abline(h=upperCI,col="grey",lty=2)
abline(h=lowerCI,col="grey",lty=2)
grid()
mtext("Cross Correlation \n (All Diarrhea)",side=2,at=0.55,line=1,cex=0.95)
mtext("Cross Correlation \n (Confirmed Cholera)",side=2,at=0,line=1,cex=0.95)
mtext("Lag (days)",side=1,line=1.25,cex=0.95)
dev.off()

#ccf(trimmed.nonc.diarrhea.daily,pred)#rowSums(cust.daily.norm,na.rm=T))

## now lets model differnt lags
aics <- matrix(nrow=11,ncol=3)
for (i in 0:10){
    fit.c <- glm(lag.and.pad(trimmed.cholera.daily,i)~pred,family="poisson")
    fit.all <- glm(lag.and.pad(trimmed.all.diarrhea.daily,i)~pred,family="poisson")
    fit.nonc <- glm(lag.and.pad(trimmed.nonc.diarrhea.daily,i)~pred,family="poisson")

    aics[i+1,1] <- AIC(fit.c)
    aics[i+1,2] <- AIC(fit.all)
    aics[i+1,3] <- AIC(fit.nonc)
}

(best.lags <- c(0:10)[apply(aics,2,which.min)])

par(mfrow=c(3,1))
plot(lag.and.pad(trimmed.cholera.daily,best.lags[1])[-c(1)],rowSums(cust.daily.norm,na.rm=T)[-c(1)])
lines(smooth.spline(lag.and.pad(trimmed.cholera.daily,best.lags[1])[-1],
                  rowSums(cust.daily.norm,na.rm=T)[-1],tol=.1),col="grey")
abline(h=0,lty=2)
plot(lag.and.pad(trimmed.all.diarrhea.daily,best.lags[2])[-c(1:7)],rowSums(cust.daily.norm,na.rm=T)[-c(1:7)])
lines(smooth.spline(lag.and.pad(trimmed.all.diarrhea.daily,best.lags[2])[-c(1:7)],
                    rowSums(cust.daily.norm,na.rm=T)[-c(1:7)]),col="grey")
abline(h=0,lty=2)

plot(lag.and.pad(trimmed.nonc.diarrhea.daily,best.lags[3])[-c(1:7)],rowSums(cust.daily.norm,na.rm=T)[-c(1:7)])
lines(smooth.spline(lag.and.pad(trimmed.nonc.diarrhea.daily,best.lags[3])[-c(1:7)],
                    rowSums(cust.daily.norm,na.rm=T)[-c(1:7)]),col="grey")
abline(h=0,lty=2)

## probably want quasi-poisson
fit.cholera <- glm(lag.and.pad(trimmed.cholera.daily,best.lags[1])~pred ,family="quasipoisson")
fit.all <- glm(lag.and.pad(trimmed.all.diarrhea.daily,best.lags[2])~pred ,family="quasipoisson")
fit.nonc <- glm(lag.and.pad(trimmed.nonc.diarrhea.daily,best.lags[3])~pred ,family="quasipoisson")

lag = 7
smoothcases=c(rep(NA,lag),smooth.spline((lag.and.pad(trimmed.all.diarrhea.daily,lag))[-c(1:lag)])$y)
glm(lag.and.pad(trimmed.all.diarrhea.daily,best.lags[2])~pred ,family="quasipoisson")

## ## let's try a little cross-validation


## xy <- data.frame(cases=lag.and.pad(trimmed.all.diarrhea.daily,best.lags[2]),ors=pred)
## xy <- xy[-which(is.na(rowSums(xy))),]

## fit.all <- glm(cases ~ ors,family="quasipoisson",data=xy)

## boot_me <- function(preddata,fit,data){
##     tmp_dat <- data[sample(nrow(data),size=nrow(data),replace=T),]
##     tmp_fit <- update(fit.all,data=tmp_dat)
##     tmp_pred <- predict(tmp_fit,type='response',newdata=preddata)
##     rpois(length(tmp_pred),lambda=tmp_pred)
## }

## simvals <- raply(500,boot_me(preddata=data.frame(ors=xy[,2]),fit=fit.all,data=xy))
## w = apply(simvals,2,function(x) quantile(x,c(.025,.975)))
## plot(xy[,1],ylim=c(0,10))
## polygon(x=c(1:190,190:1),y=c(w[1,],rev(w[2,])),col=AddAlpha('grey',.1),border="white",bty="n")

## choose a random month period

## months <- month(date.seq) %>% as.numeric
## unique.months <- unique(months)

## pred.list <- list()
## for (i in seq_along(unique.months)){
##     ## create temp dataset
##     rem.inds <- which(months==unique.months[i])
##     y.chol <- lag.and.pad(trimmed.cholera.daily[-rem.inds],best.lags[1])
##     x.sales <- pred[-rem.inds]
##     fit.tmp <- glm(y.chol~x.sales,family="quasipoisson")
##     pred.list[[i]] <-
##         cbind(pred=
##                 predict(fit.tmp,data.frame(x.sales=pred[rem.inds]),interval="prediction",type="response",se.fit=TRUE)
##              ,
##                             truth=lag.and.pad(trimmed.cholera.daily[rem.inds],best.lags[1]))
## }


#plot(date.seq,rowSums(cust.daily.norm,na.rm=T),type="h",col=ifelse(rowSums(cust.daily.norm,na.rm=T)>0,1,3))
pdf("Figures/meancommunitywideors.pdf",width=6,height=3)
palette(brewer.pal(4,"RdBu"))
par(mfrow=c(1,1),mar=c(.5,1,.5,0),oma=c(2,2,1,1),mgp=c(2,.2,0),tck=-.01)
date.grid <- seq.Date(from=as.Date("2013-04-1"),to=as.Date("2013-10-30"),by="week")
date.big <- seq.Date(from=as.Date("2013-04-1"),to=as.Date("2013-11-1"),by="month")
plot(date.seq,
     apply(cust.daily.norm,1,function(x) mean(x,na.rm=T)),
     type="h",
     col=ifelse(rowSums(cust.daily.norm,na.rm=T)>0,1,4),
     xaxt="n")
axis(1,at=date.big,labels=format(as.Date(date.big),"%b"))
abline(v=date.grid,col=AddAlpha("grey",.3),lty=2)
abline(h=seq(-40,40,by=.25),col=AddAlpha("grey",0.5),lty=2)
## plot(date.seq,trimmed.all.diarrhea.daily,type="h",col="brown")
## lines(date.seq,trimmed.cholera.daily,type="h",col="orange")
## abline(v=date.grid,col=AddAlpha("grey",.3),lty=2)
## abline(h=seq(-40,40,by=1),col=AddAlpha("grey",0.5),lty=2)
## legend("topright",c("All Diarrhea","Confirmed Cholera"),lty=1,col=c("brown","orange"),bty="n")
mtext("Community-Wide ORS Sales",side=2,outer=2,cex=.95,line=.2)
## mtext("Medically Attended Cases",side=2,outer=2,at = .25,cex=.95,line=.2)
dev.off()


n <- 10000
b0 <- .1
b1 <- .3
x <- runif(n,-1,1)
mu <- exp(b0+b1*x)
y <- rpois(n,mu)
predict(
    fit =  glm(y~x,family="poisson")
    range(y)
    range(predict(fit,type="response"))
   plot( predict(fit,data.frame(x=seq(-1,1,length=100)),type="response"))

   ,type='response')



#########################
## Temperature and ORS ##
#########################
require(gridExtra)
## bring in data from weather underground (VGHS station)
wdat <- read.csv("Data/dhaka_weather_wunder.csv",colClasses=c("Date","numeric","numeric","numeric"))

ggplot(wdat) + geom_point(aes(Date,Max_TemperatureC)) + geom_smooth(aes(Date,Max_TemperatureC)) -> maxT
ggplot(wdat) + geom_point(aes(Date,Min_TemperatureC)) + geom_smooth(aes(Date,Min_TemperatureC)) -> minT
ggplot(wdat) + geom_point(aes(Date,Mean_TemperatureC)) + geom_smooth(aes(Date,Mean_TemperatureC)) -> meanT

to.pdf(grid.arrange(maxT,meanT,minT,nrow=3),filename="Figures/minmaxtemp.pdf",width=10,height=6)
wdat$ors <- pred
glimpse(wdat)
lm(pred ~ Max_TemperatureC,data=wdat) %>% summary
lm(pred ~ Mean_TemperatureC,data=wdat) %>% confint


ggplot(wdat) + geom_smooth(aes(x=Max_TemperatureC,y=pred),method='glm') +
#  geom_smooth(aes(x=Max_TemperatureC,y=pred),col=2,se = FALSE) +
  geom_point(aes(x=Max_TemperatureC,y=pred)) -> gg
to.pdf(print(gg),filename="Figures/max_temp_ors.pdf")

ggplot(wdat) + geom_smooth(aes(x=Mean_TemperatureC,y=pred),method='glm') +
  geom_point(aes(x=Mean_TemperatureC,y=pred)) +
  ylab('Mean Community-Wide Excess ORS') +
  xlab("Mean Daily Temperature [C]") -> gg
to.pdf(print(gg),filename="Figures/mean_temp_ors.pdf")

## is temperature associated with disease in univariate analysis?
meantemp_cent <- (wdat$Mean_TemperatureC - mean(wdat$Mean_TemperatureC))
fit.c.temp <-   lm(trimmed.cholera.daily~meantemp_cent)
fit.all.temp <- lm(trimmed.all.diarrhea.daily~meantemp_cent)
fit.nonc <-  glm(trimmed.nonc.diarrhea.daily~lag.and.pad(meantemp_cent,best.lags[3]),family="poisson")

## probably want quasi-poisson
fit.cholera <- lm(trimmed.cholera.daily~meantemp_cent)
fit.all <- glm(lag.and.pad(trimmed.all.diarrhea.daily,best.lags[2])~pred + meantemp_cent,family="quasipoisson")
fit.nonc <- glm(lag.and.pad(trimmed.nonc.diarrhea.daily,best.lags[3])~pred + meantemp_cent,family="quasipoisson")

## ggplot(wdat) + geom_smooth(aes(x=Min_TemperatureC,y=pred),method='glm') +
##   geom_smooth(aes(x=Min_TemperatureC,y=pred),col=2,se = FALSE) +
##   geom_point(aes(x=Min_TemperatureC,y=pred))



###########################################
## Revised Temperature Adjusted Analyses ##
###########################################

## go back through model selection
aics <- matrix(nrow=11,ncol=3)
effects <- array(dim=c(11,3,3))

library(MASS)
for (i in 0:10){
    ## fit.c <- glm(lag.and.pad(trimmed.cholera.daily,i)~pred  ,family="poisson")
    ## fit.all <- glm(lag.and.pad(trimmed.all.diarrhea.daily,i)~pred ,family="poisson")
    ## fit.nonc <- glm(lag.and.pad(trimmed.nonc.diarrhea.daily,i)~pred ,family="poisson")

    fit.c <- glm.nb(lag.and.pad(trimmed.cholera.daily,i)~pred)
    fit.all <- glm.nb(lag.and.pad(trimmed.all.diarrhea.daily,i)~pred)
    fit.nonc <- glm.nb(lag.and.pad(trimmed.nonc.diarrhea.daily,i)~pred)

    ## fit.c <- glm(lag.and.pad(trimmed.cholera.daily,i)~rowSums(cust.daily.norm,na.rm=T),family="poisson")
    ## fit.all <- glm(lag.and.pad(trimmed.all.diarrhea.daily,i)~rowSums(cust.daily.norm,na.rm=T),family="poisson")
    ## fit.nonc <- glm(lag.and.pad(trimmed.nonc.diarrhea.daily,i)~rowSums(cust.daily.norm,na.rm=T),family="poisson")
    aics[i+1,1] <- AIC(fit.c)
    aics[i+1,2] <- AIC(fit.all)
    aics[i+1,3] <- AIC(fit.nonc)
    effects[i+1,1,] <- c(exp(coef(fit.c)[2]),exp(confint(fit.c)[2,]))
    effects[i+1,2,] <-  c(exp(coef(fit.all)[2]),exp(confint(fit.all)[2,]))
    effects[i+1,3,] <-  c(exp(coef(fit.nonc)[2]),exp(confint(fit.nonc)[2,]))

}

(best.lags_temp <- c(0:10)[apply(aics,2,which.min)])

round(apply(aics,2,function(x) x - min(x)),1)

## new estimates
library(mgcv)
fit.cholera <- glm(lag.and.pad(trimmed.cholera.daily,best.lags_temp[1])~pred + meantemp_cent,family="quasipoisson")
summary(fit.cholera)
exp(coef(fit.cholera))
exp(confint(fit.cholera))
fit.all <- glm(lag.and.pad(trimmed.all.diarrhea.daily,best.lags_temp[2])~ pred + meantemp_cent,family="quasipoisson")
summary(fit.all)
exp(confint(fit.all))
fit.all <- glm(lag.and.pad(trimmed.all.diarrhea.daily,7)~ pred + meantemp_cent,family="quasipoisson")
summary(fit.all)
exp(coef(fit.all))
exp(confint(fit.all))

fit.nonc <- glm(lag.and.pad(trimmed.nonc.diarrhea.daily,best.lags_temp[3])~pred +meantemp_cent,family="quasipoisson")
summary(fit.nonc)

## quick look at AGE or ORS recipients
library(lubridate)
dig_ors_dat = read.csv("Data/digitized_ors_records.csv")
dig_ors_dat$rdate <- as.Date(dig_ors_dat$rdate,format="%Y-%m-%d")
dig_ors_dat$bdate <- as.Date(dig_ors_dat$bdate,format="%Y-%m-%d")
date.seq <- seq.Date(as.Date("2013-04-17"),as.Date("2013-10-30"),by=1)
dig_ors_dat <- subset(dig_ors_dat,bdate %in% date.seq)
dig_ors_dat$week <- week(dig_ors_dat$bdate)
dig_ors_dat %>% group_by(week) %>%
  summarize(prop_lt15=mean(yr15le == 1,na.rm=T))%>%
  ggplot() + geom_point(aes(week,y=prop_lt15)) + stat_smooth(aes(week,prop_lt15)) +
  ylab("proportion of intended ORS \n receiptients under 15-years old") -> gg
if(remake_pdfs){
to.pdf(print(gg),filename = "Figures/age_by_week_ors.pdf",width=7,height=4)
}
mean(dig_ors_dat$yr15le == 1,na.rm=T)

load("Generated_Data/q2cln.Rda") #q2c1n object has
q2cln[which(q2cln[,2] == "0001003"),"indcc"] <- 0 # correction
index.ids <- q2cln[q2cln[,"indcc"] == 1,'datamemberid']
quantile(q2cln[which(q2cln[,'datamemberid'] %in% index.ids),]$q1_years)



##########################################
## Plotting Pharmacy Locations with OSM ##
##########################################
library(sp)
library(maptools)
library(ggmap)

## bring in pharmacy locations
pharm.locs.lookup <- read.csv("Data/PharmacyRecruitement/Pharmacy_Registration_Data.csv",as.is=T) #actual ones enrolled but the GPS on these doesn't seem correct
pharm.locs <- read.csv("Data/All_Coordinate_of_Pharmacies.csv")
pharm.locs$GPS.AP.ID <- as.numeric(sapply(as.character(pharm.locs[,1]),function(x) strsplit(x,split="AP")[[1]][2]))
pharm.locs <- join(pharm.locs.lookup[,c("Pharmacy.ID","GPS.AP.ID")],pharm.locs)
colnames(pharm.locs)[4:5] <- c("latitude","longitude")

## bring in Arichpur border
arichpur_border <- read.csv("Data/Boundary_of_Arichpur.csv")[,1:3]
ch <- chull(arichpur_border[,2:3])
coords <- arichpur_border[c(ch, ch[1]), ]  # closed polygon
plot(arichpur_border[,2:3])
lines(coords[,2:3],col="red")

## bring in case locations
dat=read.csv("Generated_Data/clinic_spatial_mlva_wide.csv")
dat = subset(dat, Sample.Loc!="pharm")
glimpse(dat)
dat[which(is.na(dat$longitude)),] %>% length
dat$Household.ID[which(duplicated(dat$Household.ID))] <- 9999

get_map(location = c(lat = 23.888, lon = 90.406), zoom = 16,
scale = "auto", maptype = "roadmap", source = "osm",
messaging = FALSE, urlonly = FALSE, filename = "ggmapTemp",
crop = TRUE, color = "bw", language = "en-EN") -> my_map
                                        #fortify(vac_areas) -> vac_areas_f
cols <- brewer.pal(3,"Dark2")

## now a little magic to make a scale bar
bb <- attr(my_map,"bb")

sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                   lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                   lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                   lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))


sbar$distance = distHaversine(long = c(sbar$lon.start,sbar$lon.end),
                              lat = c(sbar$lat.start,sbar$lat.end))

ptspermm <- 2.83464567  # need this because geom_text uses mm, and themes use pts. Urgh.

ggmap(my_map, extent = "normal",maprange = FALSE) +
  geom_point(data=pharm.locs,aes(x=longitude,y=latitude),alpha=.95,pch=17,col=cols[1]) +
  geom_point(data=dat,aes(x=longitude,y=latitude),col=cols[2]) +
  geom_segment(data = sbar,
               aes(x = lon.start,
                   xend = lon.end,
                   y = lat.start,
                   yend = lat.end)) +
                     geom_text(data = sbar,
            aes(x = (lon.start + lon.end)/2,
                y = lat.start + 0.025*(bb$ur.lat - bb$ll.lat),
                label = paste(format(distance,
                    digits = 3,
                    nsmall = 0),
                    'meters')),
            hjust = 0.5,
            vjust = 0,
            size = 8/ptspermm)  +
             coord_map(projection="mercator",
                       xlim=c(bb$ll.lon, bb$ur.lon),
                       ylim=c(bb$ll.lat, bb$ur.lat)) +
                         xlab('longitude') + ylab('latitude') -> gg
to.pdf(print(gg),filename = "Figures/pharmacies_and_cases_osm.pdf")
png(file="Figures/pharmacies_and_cases_osm.png",width=4000,height=3500,res=500)
print(gg)
dev.off()
bitmap("pharmacies_and_cases_osm.tiff", height = 4, width = 4, units = 'in', type="tifflzw", res=300)
print(gg)

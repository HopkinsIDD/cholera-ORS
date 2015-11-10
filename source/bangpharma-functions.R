# functions to be used on exported bangpahrma data
library(lubridate)
library(dplyr)
library(scales)

##' Processes data downloaded from bangpharma
##' TODO: assign day cleverly. If it is early early morning than it should be considered a report for the previous day
##' @param bp
##' @return
##' @author Andrew Azman
clean.new.data <- function(bp,adjust.probable.type.3s=T){
    ## drop blank sales
    bp <- bp[which(!is.na(bp[,1])),]
    ## drop our phamacy
    bp <- subset(bp,sale.pharmacy_id != 250)
    ## strip out time zone
    bp$sale.created_at <- gsub("UTC ",replacement="",x=bp$sale.created_at)
    ## adjut for dhaka time vs. UTC time
    bp$sale.created_at <- as.POSIXct(bp$sale.created_at,format="%a %b %d %X %Y",tzone="UTC") + 6*60*60

    h.str <- as.numeric(format(bp$sale.created_at, "%H")) +
               as.numeric(format(bp$sale.created_at, "%M"))/60


    ## deal with known mistakes reported by the field staff
    # may 15th sale by 271 for 133 ors should be 13
    bp$sale.ors[with(bp,which(sale.pharmacy_id == 271 & sale.ors == 133))] <- 13

    ## rule: if there are more than 3 type 1 entries in a 30 minute period it is probably an error. Take the last type 1.
    tmp.index <- which(bp$sale.created_at %within% new_interval(dmy("25-June-2013"),dmy("29-June-2013")) & bp$sale.pharmacy_id ==293)
    bp <- bp[-tmp.index[1:(length(tmp.index)-1)],]


    if (adjust.probable.type.3s){
        warning("Assuming that those suscpisously high reports are actually type 3")

        bp$sale.report_type[which(bp$sale.pharmacy_id == 275 & bp$sale.ors == 175)]<- 3
        bp$sale.report_type[which(bp$sale.pharmacy_id == 272 & bp$sale.ors == 109)] <- 3
        bp$sale.report_type[which(bp$sale.pharmacy_id == 272 & bp$sale.ors == 95)] <- 3

        ## 174 by 275 on 6-14-2013 seems like it was probably a type 3 report
        ## 109 by 272 on 6-15-2013 seems like it was probably a type 3 report
        ## 120 by 287 on 6-20-2013 seems like it was probably a type 3 report
    }


    ## if they are type 1 then we will add all sales from each day

    ## deal with sales type 2s - replace the previous value if there was one the same day or just keep it
    bp.rev <- ddply(bp,.(sale.pharmacy_id),function(x) {
        x <- x[order(x[,"sale.created_at"]),] # sort just to make sure
        type2s <- which(x$sale.report_type == 2)
        rem.ind <- c()
        if (length(type2s) > 0){
            for (i in seq_along(type2s)){
                if ((x[type2s[i]-1,"sale.created_at"] - x[type2s[i]-1,"sale.created_at"]) < 1){
                    rem.ind <- c(rem.ind,type2s[i]-1)
                }
            }
            if (length(rem.ind) > 0) x <- x[-rem.ind,]
        }
        x
    })


    print("section 3")
    ## now deal with sales type 3's
    bp.rev <- ddply(bp.rev,.(sale.pharmacy_id),function(x) {
        x <- x[order(x[,"sale.created_at"]),] # sort just to make sure
        type3s <- which(x$sale.report_type == 3)  # which ones are type 3s?
        rm.inds <- c()
        if (length(type3s) > 0){
            for (i in seq_along(type3s)){
                ## how many days has it been since they last reported?
                ## -1 b/c this counts the day of the last report
                days.since.last <- ceiling(x[type3s[i],"sale.created_at"] - x[type3s[i]-1,"sale.created_at"]) - 1

                if (length(days.since.last) != 0 && days.since.last > 0){
                    tmp.df <- x[rep(type3s[i],days.since.last),]

                    ## update the times for each
                    for (j in 1:nrow(tmp.df)){
                        tmp.df[j,"sale.created_at"] <- tmp.df[1,"sale.created_at"] - (60*60*24)*(j -1)
                                        #tmp.df[j,"week"] <- week(tmp.df[j,"sale.created_at"])
                                        #tmp.df[j,"day"] <- day(tmp.df[j,"sale.created_at"])
                                        #tmp.df[j,"yday"] <- yday(tmp.df[j,"sale.created_at"])
                    }

                    tmp.df[,"sale.report_type"] <- 4 # mark it as type 4
                    tmp.df[,"sale.ors"] <- x[type3s[i],"sale.ors"]/as.numeric(days.since.last)

                    ## make entries of type 4 that divide the sales evenly over the missed days
                    ## remove the original entry
                    rm.inds <- c(rm.inds,type3s[i])
                    x <- rbind(x,tmp.df)
                } else {
                    x[type3s[i],"sale.report_type"] <- 1 # assuming that is an additional sale if this happens
                }
            }
            if (length(rm.inds) > 0) x <- x[-rm.inds,]
        }

        x
    })

    ## add yearday, weeks, day
    bp.rev$sale.created_at <- as.Date(bp.rev$sale.created_at)
    bp.rev$week <- week(bp.rev$sale.created_at)
    bp.rev$day <- day(bp.rev$sale.created_at)
    bp.rev$yday <- yday(bp.rev$sale.created_at)

    ## now go through and add sales together from the same day

    bp.rev <- ddply(bp.rev,.(sale.pharmacy_id,yday),function(x){
        if (nrow(x) > 1){
            combined.ors <- sum(x[,"sale.ors"])
            x <- x[nrow(x),]
            x[,"sale.ors"] <- combined.ors
        }
        x
    })

    print("end")
    # convert to Date format (and drop the time bit by default)


    return(bp.rev)
}

make.image.plot <- function(aggregate.to="day"){
    require(fields)
#    bp <- CleanNewData(bp.raw)
    bp[,"week"] <- factor(bp[,"week"])
    bp[,"yday"] <- factor(bp[,"yday"])
                                        # get sales from each day for each pharmacy NA if no data
    ddply(bp,.(sale.pharmacy_id,yday),.drop=FALSE,function(x){
        if (nrow(x) > 0){
            return(x$sale.ors)
        } else {
            return(NA)
        }}
          ) -> ors.dat.per.pid.per.time

    ddply(ors.dat.per.pid.per.time,.(sale.pharmacy_id),function(x) {
        x[,3]
    }) -> ors.mat


    ors.mat <- as.matrix(ors.mat)
    rownames(ors.mat) <-ors.mat[,1]
    ors.mat <- ors.mat[,-1]
    colnames(ors.mat) <- levels(factor(bp[,"yday"]))

    ## now get the mean per day foor use in weeks
    daily.mean <- apply(ors.mat,1,function(x) mean(x,na.rm=T))

    ## normalize rows by the daily mean
    normalized.ors <- apply(ors.mat,1,function(x) {
        x/mean(x,na.rm=T)
    })


    if (aggregate.to == "week"){
        ddply(bp,.(sale.pharmacy_id,week),.drop=FALSE,function(x){
            if (nrow(x) > 0){
                return(sum(x$sale.ors))
            } else {
                return(NA)
            }}
              ) -> ors.dat.per.pid.per.time

        # how many reporting days in each week for each pharmacy?
        ddply(ddply(bp,.(sale.pharmacy_id,week),.drop=FALSE,nrow),.(sale.pharmacy_id),function(x){
            x[,3]
        }) -> days.reporting
        rownames(days.reporting) <- days.reporting[,1]
        days.reporting <- days.reporting[,-1]

        ddply(ors.dat.per.pid.per.time,.(sale.pharmacy_id),function(x) {
            x[,3]
        }) -> ors.mat
        ors.mat <- as.matrix(ors.mat)
        rownames(ors.mat) <-ors.mat[,1]
        ors.mat <- ors.mat[,-1]
        colnames(ors.mat) <- levels(factor(bp[,"week"]))
        normalized.ors <- sapply(1:nrow(ors.mat),function(x) ors.mat[x,]/(daily.mean[x]*as.matrix(days.reporting[x,])))
    }

    par(mar=c(4,4,4,6),mgp=c(2,3,0))
    image(log1p(normalized.ors),col=colorRampPalette(brewer.pal(9,"Blues"))(10000),axes=FALSE,xlab=sprintf("%s of year",aggregate.to),ylab="pharmacy id")
    mtext(text=rownames(ors.mat), side=2, line=0.3, at=seq(0,1,length=nrow(ors.mat)), las=1, cex=0.7)
    mtext(text=colnames(ors.mat), side=1, line=0.3, at=seq(0,1,length=ncol(ors.mat)), las=2, cex=0.5)
    image.plot(log1p(normalized.ors), legend.only=T,col=colorRampPalette(brewer.pal(9,"Blues"))(10000))
    title("ORS Sales per Pharmacy (normalized by pharmacy)")
}


## for each pharmacy
## total number of calls made to BP
## total number of dropped calls to BP
## total sales recorded (by type)
## Number of days with no sucsessful reports
## min median max sales
single.pharm.report <- function(data,past.n.days){
    if (!missing(past.n.days)) data <- subset(data,as.Date(sale.created_at) > today() - past.n.days)

    if (nrow(data) == 0) return(rep(0,11)) # just in case no data in this time range

    dropped.calls <- subset(data,sale.outcome == 0)     # number of dropped calls
    tmp <- subset(data,sale.outcome != 0)

    sales.count.by.type <- sapply(1:3,function(x) sum(tmp$sale.report_type == x))

                                        # number of days with no report
    possible.days <- seq.Date(from=as.Date(min(data$sale.created_at,na.rm=T)),to=as.Date(today()),by=1)
    days.reporting <- length(unique(as.Date(tmp$sale.created_at)))

    sales.quants <- quantile(tmp$sale.ors,c(0,0.25,.5,0.75,1))

    out <- round(c(nrow(dropped.calls),
             nrow(tmp),
             days.reporting/length(possible.days),
             sales.count.by.type,
             sales.quants),2)

    return(out)
}

all.pharm.report <- function(dat,...){
    rc <- ddply(dat,.(sale.pharmacy_id),function(x) {
        single.pharm.report(x,...)
    })

    colnames(rc) <- c("pid","dropped.calls.num",
                      "total.suc.calls.num",
                      "pct.days.reporting",
                      "sales.count.1",
                      "sales.count.2",
                      "sales.type.3",
                      "ors.cust.min",
                      "ors.cust.25pct",
                      "ors.cust.50pct",
                      "ors.cust.75pct",
                      "ors.cust.max")
    return(rc)
}



##'
##' @param tongi.file path to the csv version of Culture Tracking Sheet Template tongi sheet
##' @param icddrb.file path to the csv version of Culture Tracking Sheet Template icddrb sheet
##' @param pharm.hhl.file
##' @param q1.file
##' @param max.valid.date
##' @param aggregated
##' @param with.gps - if unaggregated, do we want the gps points too
##' @param hosptial.cases.only - if true willl exclude hospital cases
##' @return data frame of dates of sample collection, cholera result, and clinic from which the patient came from
##' @author Andrew Azman
combine.and.anon.clinic.sample.data <- function(
    tongi.file="Data/MainData/Current/TongiSamples_current.csv",
    icddrb.file="Data/MainData/Current/ICDDRBSamples_current.csv",
    pharm.hhl.file="Data/MainData/Current/HouseholdSamples_current.csv",
    q1.file = "Data/MainData/OldDataSets/householdques1_15_sep_2013.dta",
    max.valid.date=today(),
    aggregated=FALSE,
    with.gps=FALSE,
    hospital.cases.only=FALSE){

    require(foreign)
    require(lubridate)
    tongi.sample.data <- read.csv(tongi.file,as.is=T,colClasses=c("Study.ID"="character"))
    icddrb.sample.data <- read.csv(icddrb.file,as.is=T)
    hhl.sample.data <- read.csv(pharm.hhl.file,as.is=T)
    q1 <- read.dta(q1.file)

    tongi.dates <- as.Date(tongi.sample.data[,"Date.of.Sample.Collection"],format="%d-%b-%y")
    icddrb.dates <- as.Date(icddrb.sample.data[,"Date.of.Sample.Collection"],format="%d-%b-%y")
    ## some of the houshold samples are don't have dates since they didn't have sampels taen
    hhl.sample.data <- subset(hhl.sample.data,Sample.ID != "")
    hhl.dates <- as.Date(hhl.sample.data[,"Date.of.Sample.Collection"],format="%d-%b-%y")

    if (with.gps){
        warning("THis may not work as expected, check this before using!")
        q1.sub <- q1[,c("dataid","longitude","latitude")]
        q1.sub[,2] <- as.numeric(q1.sub[,2])
        q1.sub[,3] <- as.numeric(q1.sub[,3])

        tongi.sample.data$dataid <- sapply(strsplit(tongi.sample.data[,"Study.ID"],split="-"),function(x) x[2])
        tongi.sample.data <- merge(q1.sub,tongi.sample.data,by="dataid",all.y=T)

        icddrb.sample.data$dataid <- sapply(strsplit(icddrb.sample.data[,"Study.ID"],split="-"),function(x) x[2])
        icddrb.sample.data <- merge(q1.sub,icddrb.sample.data,by="dataid",all.y=T)

        hhl.sample.data$dataid <- sapply(strsplit(hhl.sample.data[,"Household.ID"],split="-"),function(x) x[1])
        hhl.sample.data <- merge(q1.sub,hhl.sample.data,by="dataid",all.y=T)
    }

    ## check that the dates seem to both in the expected format and dueing the correct time frame
    if (any(is.na(tongi.dates))) stop("some of the tongi dates are not formatted correctly")
    if (any(is.na(icddrb.dates))) stop("some of the icddrb dates are not formatted correctly")
    if (any(is.na(hhl.dates))) stop("some of the hhl dates are not formatted correctly")

    int <- new_interval(ymd("2013-03-25"),max.valid.date)
    if (any(!icddrb.dates %within% int)) stop("some icddrb dates don't seem to be in the correct range")
    if (any(!tongi.dates %within% int)) stop("some tongi dates don't seem to be in the correct range")
    if (any(!hhl.dates %within% int)) stop("some houshold/pharmacy dates don't seem to be in the correct range")

    ## grab culture results
    tongi.culture.results <- tongi.sample.data[,"Culture.Result"]
    icddrb.culture.results <- icddrb.sample.data[,"Culture.Result"]
    hhl.culture.results <- hhl.sample.data[,"Culture.Result"]

    ## check to make sure they are either positive or negative (or "" for results that haven't come in yet
    if (!all(tongi.culture.results %in% c("Negative","Positive",""))) stop("some tongi culture results are either mispelled or invalid")
    if (!all(icddrb.culture.results %in% c("Negative","Positive",""))) stop("some icddrb culture results are either mispelled or invalid")
    if (!all(hhl.culture.results %in% c("Negative","Positive",""))) stop("some houshold culture results are either mispelled or invalid")

    if (any(tongi.culture.results == "")) {
        drop.me <- which(tongi.culture.results == "")
        tongi.culture.results[drop.me] <- rep("Negative",length(drop.me))
        cat("Imputing Negatives for",length(drop.me),"entries from tongi due to incomplete culture results.\n")
        ## tongi.dates <- tongi.dates[-drop.me]
        ## tongi.culture.results <- tongi.culture.results[-drop.me]
        ## tongi.sample.data <- tongi.sample.data[-drop.me,]
        ## cat("Dropping",length(drop.me),"entries from tongi due to incomplete culture results.\n")
    }

    if (any(icddrb.culture.results == "")) {
        drop.me <- which(icddrb.culture.results == "")
        ## for now going to assume they are negative
        icddrb.culture.results[drop.me] <- rep("Negative",length(drop.me))
        cat("Imputing Negatives for",length(drop.me),"entries from icddrb due to incomplete culture results.\n")
        ## icddrb.dates <- icddrb.dates[-drop.me]
        ## icddrb.culture.results <- icddrb.culture.results[-drop.me]
        ## icddrb.sample.data <- icddrb.sample.data[-drop.me,]
        ## cat("Dropping",length(drop.me),"entries from icddrb due to incomplete culture results.\n")
    }

    if (any(hhl.culture.results == "")) {
            drop.me <- which(hhl.culture.results == "")
            hhl.dates <- hhl.dates[-drop.me]
            hhl.culture.results <- hhl.culture.results[-drop.me]
            hhl.sample.data <- hhl.sample.data[-drop.me,]
            cat("Dropping",length(drop.me),"entries from hhl due to incomplete culture results.\n")
    }

    if (hospital.cases.only){
        dat <- data.frame(date=c(tongi.dates,icddrb.dates))
        dat$cholera <- c(ifelse(tongi.culture.results == "Positive",1,0),
                         ifelse(icddrb.culture.results == "Positive",1,0)
                         )

        dat$location <- c(rep("tongi",length(tongi.dates)),
                          rep("icddrb",length(icddrb.dates)))
        if (with.gps){
            dat$dataid <-  c(tongi.sample.data$dataid,icddrb.sample.data$dataid)
            dat$longitude <- c(tongi.sample.data$longitude,icddrb.sample.data$longitude)
            dat$latitude <- c(tongi.sample.data$latitude,icddrb.sample.data$latitude)
        }

    } else {
        dat <- data.frame(date=c(tongi.dates,icddrb.dates,hhl.dates))
        dat$cholera <- c(ifelse(tongi.culture.results == "Positive",1,0),
                         ifelse(icddrb.culture.results == "Positive",1,0),
                         ifelse(hhl.culture.results == "Positive",1,0)
                         )

        dat$location <- c(rep("tongi",length(tongi.dates)),
                          rep("icddrb",length(icddrb.dates)),
                          rep("hhl.pharm",length(hhl.dates)))

        if (with.gps){
            dat$dataid <-  c(tongi.sample.data$dataid,icddrb.sample.data$dataid,hhl.sample.data$dataid)
            dat$longitude <- c(tongi.sample.data$longitude,icddrb.sample.data$longitude,hhl.sample.data$longitude)
            dat$latitude <- c(tongi.sample.data$latitude,icddrb.sample.data$latitude,hhl.sample.data$latitude)
        }
    }

    if (aggregated){
        ## first this stuff
        all.dates <- data.frame(date=seq.Date(from=min(dat$date),by="day",to=max(dat$date)))
        dat <- merge(all.dates,dat,all.x=T)
        dat[is.na(dat)] <- 2

        ## now aggregate by day
        dat <- ddply(dat,.(date),function(x) {
            c("non.c.diarrhea"=sum(x$cholera == 0),
              "cholera"=sum(x$cholera==1))
        },.drop=FALSE)

    }
    return(dat)
}


##' Creates a daily dataset with both ORS and Cholera counts
##' @param daily.ors.reports
##' @param comb.cultures
##' @param trim.by.dates only return dates for which we have both ors and cholera
##' @param day.lag days to lag by cholera by (ahead)
##' @return
##' @author Andrew Azman
create.ors.plus.cholera.dataset <- function(daily.ors.reports=comb.r,
                                            comb.cultures=comb.cultures,
                                            trim.by.dates=TRUE,
                                            day.lag=0){
    # lag cholera dates
    comb.cultures$date <- comb.cultures$date + day.lag

    if (trim.by.dates){
        rc <- merge(comb.cultures,daily.ors.reports)
    } else {
        rc <- merge(comb.cultures,daily.ors.reports,all.x=T,all.y=T)
    }
    return(rc)
}

##' helper function for lagging variabel
##' @title
##' @param x
##' @param k
##' @return
##' @author asa
lag.and.pad <- function(x, k) {
    c(rep(NA, k), x)[1 : length(x)]
}

## get dustance from lat long
distHaversine <- function(long, lat){
  dlong = (long[2] - long[1])*pi/180
  dlat  = (lat[2] - lat[1])*pi/180

  # Haversine formula:
  R = 6371;
  a = sin(dlat/2)*sin(dlat/2) + cos(lat[1])*cos(lat[2])*sin(dlong/2)*sin(dlong/2)
  c = 2 * atan2( sqrt(a), sqrt(1-a) )
  d = R * c
  return(d*1000) # in m
}

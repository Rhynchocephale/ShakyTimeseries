rm(list=ls())

options(error=recover)
#setwd("./timeseries/qrfit-data")
#library(xts)
library(changepoint)
library(zoo)

findJumps <- function(folder)
{

    fileList <- paste(folder,"residuals",list.files(paste(folder,"residuals",sep="/"), pattern = "*.dat"),sep="/")
    #fileList <- paste(folder,"smoothed",list.files(paste(folder,"smoothed",sep="/"), pattern = "*.dat"),sep="/") #with smoothed
    #fileList <- paste(folder,list.files(folder, pattern = "*.dat"),sep="/") #with raw

    for(file in fileList)
    {
        staName <- strsplit(file,"/")
        staName <- staName[[1]][length(staName[[1]])]
        staName <- strsplit(staName,"_")[[1]][1]

        print(paste("finding jumps: ",staName))

        mydata <- read.table(file, sep="\t", header=TRUE, stringsAsFactor=FALSE) #reads file

        xyz <- mydata[,3:5] #get columns

        splittedDates <- strsplit(mydata[,"epoch"],"T")
        splittedDates <- sapply(splittedDates, "[", 1)
        dates <- as.Date(splittedDates) #gets epochs string vector, splits to only get the dates, tranforms to date

        #dates2 <- duplicated(dates)
        #print(dates[dates2])

        #xyzCUSUM <- toCUSUM(xyz)

        #should make timeseries (or something very close to timeseries)
        for(i in 1:3)
        {
            timeseries <- zoo(xyz[,i],dates)
            b <- cpt.meanvar(timeseries, penalty="Asymptotic", pen.value=0.0001, method="BinSeg",Q=100, test.stat="Normal")
            a <- cpts(b) #b$cps[1,]
            #plot(dates,xyz[,i],col="blue",pch=20)
            #points(dates[a],xyz[a,i],col="red",pch=15)

            toText(dates[a],folder,staName)
        }
    }
}

toCUSUM <- function(xyz)
{
    xyzCUSUM <- xyz**2
    xyzMean = c(mean(xyz[,1]),mean(xyz[,2]),mean(xyz[,3]))
    for (i in 2:length(xyz[,1]))
    {
        xyzCUSUM[i,] <- xyzCUSUM[i-1,] + xyzCUSUM[i,]
    }

    xyzCUSUM <- xyzCUSUM - xyzMean
    return(xyzCUSUM)
}

toText <- function(dates,folder,staName)
{
    write(strftime(dates,"%Y-%m-%d"),paste(folder,'/offsets/',staName,"_jumps.dat",sep=""),append=TRUE)
}

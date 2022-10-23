#Yin-Phan Tsang 2022.07.21======================================================
# Preparing EflowState for 15 streams gauges in Lisi et al. 2022
# 16717000
# 16010000
# 16108000
# 16620000
# 16604500
# 16618000
# 16614000
# 16400000
# 16240500
# 16284200
# 16304200
# 16294100
# 16294900
# 16330000
#===============================================================================
#clean the working space
rm(list = ls(all = TRUE))
ls(all = TRUE)

#set-up working directory-------------------------------------------------------
#remeber use back slash "/" in directory
setwd('C:/Users/tsang/Dropbox/PC/Desktop/P_Otilith/R')


# install packages, only need to do it once-------------------------------------
#install.packages("dataRetrieval")
#install.packages("lmomco")
#install.packages("EcoHydRology")
#install.packages("ggplot2")
#To install the latest release of EflowStats see releases to find the latest version and use:
#remotes::install_github("USGS-R/EflowStats@v5.0.1")
#To install the latest and greatest build use the following code:
#remotes::install_github("USGS-R/EflowStats")


#[1. Set up site information for dataRetrieval]---------------------------------
#input USGS site information, reading parameters, start/end period
# 16717000
# 16010000
# 16108000<-- less than 4 otolith in 2011
# 16620000
# 16604500
# 16618000<-- less than 4 otolith in 2009
# 16614000
# 16400000
# 16240500<-- less than 4 otolith in 2011
# 16284200
# 16296500
# 16304200<-- less than 4 otolith in 2009
# 16294100
# 16294900
# 16330000<-- less than 4 otolith in 2009 and 2011
siteNumbers <- c('16717000','16010000','16108000','16620000','16604500','16618000','16614000','16400000',
                 '16240500','16284200','16296500','16304200','16294100','16294900','16330000')


#[2.calculate the metrics based on daily flow]------------------------------------------------------

library(dataRetrieval)
library(EflowStats)
library(stringr)
# # https://rdrr.io/github/USGS-R/EflowStats/f/vignettes/intro.Rmd


#----calculate CV
library('EnvStats')
#calendar year
allCV<-c()
for(year in 2001:2011){
  
  startDate <- paste((year),"-01-01",sep="")
  endDate <- paste(year,"-12-31",sep="")
  
  for(i in siteNumbers){
    print(i)
    #Get some data
    dailyQ <- readNWISdv(siteNumber = i,
                         parameterCd = "00060",
                         startDate = startDate,
                         endDate = endDate)
    print(min(dailyQ$Date))
    print(max(dailyQ$Date))
    
    #if any na, replaced with average of previous and next day
    i_na<-which(is.na(dailyQ[c("X_00060_00003")]))# found Makiki King stream NA on 2021-7-28
    if(length(i_na)>0){
      for(j in i_na){
        dailyQ[i_na,c("X_00060_00003")]<-mean(dailyQ[j-1,c("X_00060_00003")], dailyQ[j+1,c("X_00060_00003")])
      }
    }
    
    
    #Check data for completeness
    completedates<-seq(min(dailyQ$Date),max(dailyQ$Date),by=1)
    missingdates<-as.Date(setdiff(completedates,dailyQ$Date), origin="1970-01-01")
    print(paste("missing:",missingdates))
    
    makingQ<-dailyQ$X_00060_00003
    #lookup the previous dates of missingdates
    if(length(missingdates)==0){
      makingQ<-dailyQ$X_00060_00003
    }else if(length(missingdates)==1){
      i_previousdatesindailyQ<-which(dailyQ$Date==(missingdates-1))
      i_j_dailyQ<-i_previousdatesindailyQ
      makingQ<-c(dailyQ$X_00060_00003[1:i_j_dailyQ],dailyQ$X_00060_00003[i_j_dailyQ],dailyQ$X_00060_00003[(i_j_dailyQ+1):nrow(dailyQ)])
    }else if(length(missingdates)>1){
      i_previousdatesindailyQ<-which(dailyQ$Date==(missingdates[1]-1))
      i_j_dailyQ<-i_previousdatesindailyQ
      makingQ<-c(dailyQ$X_00060_00003[1:i_j_dailyQ],dailyQ$X_00060_00003[i_j_dailyQ],dailyQ$X_00060_00003[(i_j_dailyQ+1):nrow(dailyQ)])
      for(j in missingdates[2:length(missingdates)]){
        k<-1
        i_previousdatesindailyQ<-which(dailyQ$Date==(as.Date(j,origin="1970-01-01")-k))
        while(length(i_previousdatesindailyQ)==0){
          k<-k+1
          i_previousdatesindailyQ<-which(dailyQ$Date==(as.Date(j,origin="1970-01-01")-k))
        }
        i_j_dailyQ<-i_previousdatesindailyQ#previous value in dailyQ
        i_j_completedate<-which(completedates==(as.Date(j,origin="1970-01-01")-1))#previous date in completedates
        makingQ<-c(makingQ[1:i_j_completedate],dailyQ$X_00060_00003[i_j_dailyQ],makingQ[(i_j_completedate+1):length(makingQ)])
      }
    }
    print(paste("completenum",length(completedates)))
    print(paste("makingQnum",length(makingQ)))
    dailyQClean<-data.frame(date=completedates,discharge=makingQ)
    
    site_cv<-cv(dailyQClean$discharge)
    
    #preparing output table
    siteInfo <- readNWISsite(siteNumber = i)
    drainageArea <- siteInfo$drain_area_va
    temp<-c(i,drainageArea,year,site_cv)
    names(temp)<-c("USGSsiteID","Area_ac", "calendaryear","siteCV")
    allCV<-rbind(allCV,temp)
  }
}

write.csv(allCV,"allcvfor15gaugesinLisi2022_bycalendaryear.csv")


#format
#as numeric
library(ggplot2)
allCV<-as.data.frame(allCV)
for(i in c(2,4)){
  allCV[,i]<-as.numeric(allCV[,i])
}
jpeg(paste("CV of daily flow 2001-2011_15 gauges_",Sys.Date(),".jpeg",sep=""), width = 900, height = 450)
#par(mfrow = c(2, 1))
p1<-ggplot(allCV, aes(x=calendaryear, y=siteCV, group=calendaryear)) + 
        geom_boxplot()+theme_bw()+
        theme(text = element_text(size = 20))+
        geom_boxplot(data=allCV[allCV$calendaryear=="2009",],
                     aes(x = calendaryear, y=siteCV),fill="#F8644D")+
        geom_boxplot(data=allCV[allCV$calendaryear=="2011",],
                     aes(x = calendaryear, y=siteCV),fill="#00BFC4")+
        labs(title="a",x="Year", y = "CV of daily flow")
print(p1)
dev.off()


#---- recreate Fig 1a in Lisi et al. 2022, but with lines indicating the direction of flow changes
CV_FW<-read.csv('CVvsFW_takeoutlessthan4.csv')
CV_FW$calendaryear<-as.factor(CV_FW$calendaryear)
names(CV_FW)[3]<-"year"

jpeg(paste("CVvsFW_RecreateFig1a_color_",Sys.Date(),".jpeg",sep=""), width = 900, height = 575)
p2<-ggplot(CV_FW, aes(x=siteCV,y=FW,group=USGSsiteID)) + 
        geom_point(aes(color=year), size=4)+theme_bw()+
        geom_path(data=CV_FW, aes(x=siteCV,y=FW,group=USGSsiteID),
                  color="black",
                  arrow = arrow(angle = 30, length=unit(0.1, "in"), type="closed"),
                  show.legend = FALSE,
                  inherit.aes = FALSE)+
        theme_classic() +
        theme(text = element_text(size = 20))+
        #scale_shape_manual(values = c(19, 1))+
        labs(title="b",x="flow variability (CV of daily flow)", y = "Proportion of freshwater residence")
print(p2)
dev.off()

library(patchwork)
jpeg(paste("CV of daily flow 2001-2011_15 gauges and Fig1g",Sys.Date(),".jpeg",sep=""), width = 900, height = 1150)
p1/p2
dev.off()


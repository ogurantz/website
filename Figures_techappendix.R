rm(list=ls())

library(tidyverse)
library(ggplot2)
library(ggtext)
library(estimatr)
library(patchwork)
#library(causalTree)
library(fastDummies)
library(grf)
library(car)
library(aod)
library(haven)
library(interflex)
library(glue)
library(ggtext)
library(xtable)
library(dplyr)


## settings

Fig1 <- TRUE
Fig2 <- TRUE
Fig3 <- TRUE
Fig4 <- TRUE

col1 <- c("cornflowerblue")
col2 <- c("firebrick1")
col3 <- c("darkgreen")

path <- "G:/XXX"

# load data ---------------------------------------------------------------
setwd(paste(path,"Data/transfer240625",sep=""))

# treatment group
tr <- read.csv("TestSamp_24JUN2024.csv")
tr$tr <- 1

# control group
co <- list()
for(i in 1:4){
  co[[i]] <- read.csv(paste("CtrlSamp_24JUN2024_Part",i,".csv",sep=""))
}
co <- do.call(rbind, co)
co$tr <- 0

# combine
d <- bind_rows(tr,co)
rm(tr,co)


# preprocess --------------------------------------------------------------

# mailing delivered indicator
d$mdelivered <- 0
d$mdelivered[d$tr==1] <- 1
d$mdelivered[d$RETURN_DATE_1!=""] <- 0

lm(mdelivered~tr,data=d)

# outcome
d$NATZ_FILE_DT_d <- as.Date(d$NATZ_FILE_DT , format = "%m/%d/%Y")
d$n400 <- as.numeric(!is.na(d$NATZ_FILE_DT_d))*100

# rand block
d <- d %>% group_by(LPR_GROUP,COA_GROUP,ICOB_GROUP) %>% mutate(rblock=cur_group_id())
d$izip <- d$IZIP


# Fig. 1. Effect of the Letter on N-400 Submission ---------------------------------------------------------

  
  # set dates  
  dates <- seq(as.Date("2023-10-26"),as.Date("2024-06-24"),by="week")
  
  store <- matrix(NA,nrow=length(dates),ncol=9)
  for(i in 1:length(dates)){
    
    cat("Computing Results for Week:",i,"out of" ,length(dates) ,"\n")
    # compute cumulative indicator for having nat up to specific week  
    d$n400_cum <-   (as.numeric( (d$NATZ_FILE_DT_d < dates[i])==TRUE))*100
    d$n400_cum[is.na(d$NATZ_FILE_DT_d)] <- 0
    
    # compute means by group
    tout1 <- t.test(d$n400_cum[d$tr==1], conf.level = 0.95)
    tout0 <- t.test(d$n400_cum[d$tr==0], conf.level = 0.95)
    # ATE
    lout    <- lm_robust(n400_cum ~ tr,
                         fixed_effects=as.factor(rblock),
                         data = d,
                         se_type="HC2")
    # store
    store[i,] <- c(tout1$estimate,tout1$conf.int,
                   tout0$estimate,tout0$conf.int,
                   summary(lout)$coefficients[1,c(1,5,6)]
    )
  } # end loop across weeks
  
  # combine results by week
  st <- data.frame(store)
  names(st) <- c("tr_m","tr_lb","tr_ub",
                 "co_m","co_lb","co_ub",
                 "tau_m","tau_lb","tau_ub")
  st$week_no   <- (1:nrow(store))-1
  st$week_date <- dates
  
  # reshape
  st_stacked <- data.frame(m=c(st$tr_m,st$co_m),
                           lb=c(st$tr_lb,st$co_lb),
                           ub=c(st$tr_ub,st$co_ub),
                           tr=c(rep(1,nrow(st)),rep(0,nrow(st))),
                           week_no=(c(1:nrow(st),1:nrow(st)))-1)
  
  # plot
  st_stacked$tr <- factor(st_stacked$tr,levels=c(1,0),labels=c("Treatment","Control"))
  
  write.csv(st, "G:\\XXX\\fig1_R.csv", row.names=TRUE)
  
  p1 <- ggplot(st_stacked,
               aes(x=week_no,y=m,ymin=lb,ymax=ub,group=tr,fill=tr,color=tr))+
    geom_ribbon(alpha=.4,color=NA)+
    geom_line()+
    ylab("N-400 Applications Submitted (%)")+xlab("Weeks Since First Letters")+ylim(0,10)+
    scale_color_manual(name="", values = c(col2,col1))+
    scale_fill_manual(name="", values = c(col2,col1)) +  theme(legend.position = "top")+
    ggtitle("A: Cumulative N-400 Application Rate")
  
  p2 <- ggplot(st,
               aes(x=week_no,y=tau_m,ymin=tau_lb,ymax=tau_ub))+
    geom_ribbon(alpha=.4,color=NA,fill=col3)+
    geom_line(color=col3)+ylim(-0.05,1.25)+scale_y_continuous(breaks=seq(0,1.24,.25))+
    ylab("Effect of Letter on N-400 Application Rate (pp)")+xlab("Weeks Since First Letters")+
    ggtitle("B: Increase in N-400 Application Rate")
  
  p <- p1 + p2
  
  
  #####
  
  
  # merge in median income
  setwd(paste(path,"/Data",sep=""))
  zip <- read_dta("medinc.dta")
  
  d$izip <- d$IZIP
  d <- left_join(d,zip,by="izip")
  d <- data.frame(d)
  
  # run moderation effects
  d$trplot <- factor(d$tr,levels=c(0,1),labels=c("Control","Treatment"))
  d$n400 <- d$n400/100
  iout <-interflex(Y = "n400", D = "trplot", X = "medinc"
                   ,data = d, estimator = "kernel", bw=2000, na.rm = TRUE)
  
  # save
  setwd(paste(path,"/XXX/results",sep=""))
  
  p1 <- plot(iout,
             order = c("Treatment"),
             main = "B: Intention-to-treat Effects by Median Income", 
             ylab = "Effect of Letter (pp)",
             xlab="Median Income of Zip Code ($)",
             color=c(col3),
             Xdistr = "histogram", xlim = c(10000,150000), 
             ylim = c(-0.03,.04), cex.main = .8, 
             cex.lab = .8, cex.axis = .8,fill=col3,font.main = 2)
  p1
  
  dp <- data.frame(iout$est.kernel$Treatment)
  percent_no_sign <- function(x) {
    paste0(x * 100)
  }
  p1 <- ggplot(dp,aes(x=X,y=TE,ymin=lower.CI.95..,ymax=upper.CI.95..)) + 
    geom_hline(yintercept=0,col="white",size=2)+
    geom_ribbon(alpha=.4,color=NA,fill=col3)+
    geom_line(color=col3,size=1.2)+
    ylab("Effect of Letter on N-400 Application Rate (pp)")+xlab("Median Income of Zip Code ($)")+
    ggtitle("B: Intention-to-treat Effects by Median Income") + 
    coord_cartesian(xlim = c(10000,150000),ylim = c(-0.001,.028))+
    theme(plot.title = element_text(face = "bold",size=14),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11),
    )+scale_y_continuous(labels = percent_no_sign)
  
  
  write.csv(dp, "G:\\XXX\\fig4A.csv", row.names=TRUE)
  
  p2 <- predict(iout,order = c("Treatment","Control"),
                subtitle = c("Treatment","Control"), pool = T, color=c(col1,col2),
                main = "A: N400 Submission Rates by Median Income", 
                ylab = "N-400 Applications Submitted (%)",
                xlab="Median Income of Zip Code ($)",
                legend.title = "",
                xlim = c(10000,150000),ylim=c(0.05,.11),cex.main = .8,cex.axis=.9,cex.label=.9,
                Xdistr="hist")+scale_y_continuous(labels = percent_no_sign)
  
  test <- data.frame(iout$pred.kernel)
  test2 <- data.frame(iout$hist.out$mids)
  test3 <- data.frame(iout$count.tr)
  
  write.csv(test, "G:\\XXX\\fig41.csv", row.names=TRUE)
  write.csv(test2, "G:\\XXX\\fig42.csv", row.names=TRUE)
  write.csv(test3, "G:\\XXX\\fig43.csv", row.names=TRUE)
  
  p3 <- p2 + p1
  
  write.csv(dp, "G:\\XXX\\fig3A.csv", row.names=TRUE)
  
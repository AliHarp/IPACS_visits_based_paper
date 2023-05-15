######## ensure LOS params match simulation ###########

rm(list=ls())
require(doSNOW)
require(foreach)
require(tidyverse)
require(parallel)
require(readxl)
require(zoo)
require(stringr)
library(beepr)
library(moments)

rm(list=ls())

setwd("~/Dropbox/03_IPAC/P1_paper_scenarios")

##### params #####

avg_LOS <- 12.489
#avg_LOS <- 12.1304
#avg_LOS <- 12
LOS_sd <- 4
#LOS_sd <- 6


#function for custom prob dist for los of those already in service
`%>%` <- magrittr::`%>%`
freqcalc <- function(rem_los){
  rem_los2 <- data.frame(rem_los) %>% dplyr::group_by(rem_los) %>%
    dplyr::summarise(n=dplyr::n()) %>% dplyr::mutate(Freq=n/sum(n))
  remaining_los_dist <- rem_los2$Freq #new distribution
  return(remaining_los_dist)
}

#function to sample from discrete dist and create custom dist
sample_reps <- function(share, seq){
  #sample once from disc dist and from 1:sample
  rem_los_fn <- function(seq, share){
    samples <- sample(seq, 1, replace=TRUE, prob=share)#sample patient in service
    samples <- round(samples)
    los_so_far <- runif(1, 1, samples)
    return(los_so_far)
  }
  rem_los <- lapply(rep(1,10000), function(x) rem_los_fn(seq, share)) #distribution of remaining LoS
  rem_los <- unlist(rem_los)
  rem_los <- round(rem_los)
  #create custom freq dist
  remaining_los_dist_visit <- freqcalc(rem_los)
  return(remaining_los_dist_visit)
}


# LOS for patients already in service using custom prob dist #############

###### 1. los if normal 
dis_los2 <- function(){
  sample_init <- function(visit_los_mean, visit_los_sd){
    seq <- 1:(visit_los_mean + 3*visit_los_sd) #range to sample from
    seq_0.5 <- (seq[1]-0.5):(rev(seq)[1]+0.5) #sampling chance of random pt having los of..
    los_prob_0.5 <- pnorm(seq_0.5, mean=visit_los_mean, sd=visit_los_sd, log=FALSE) #prob dist
    los_prob_between <- zoo::rollapply(los_prob_0.5, width=2, FUN=function(q) q[2]-q[1])
    los_prob_between_sum <- sum(los_prob_between) #sum for normalising
    normalised_los_prob_between <- los_prob_between/los_prob_between_sum
    a_prob <- normalised_los_prob_between * seq
    adj_mean <- sum(los_prob_between * seq)
    share <- a_prob/adj_mean #los of those in service probabilities
    seq_div_2 <- seq/2 #
    newmean <- sum(share*seq_div_2)
    sample_reps(share, seq)
  }
  remaining_los_dist_visit <- sample_init(avg_LOS, LOS_sd)
  x<- sample(1:length(remaining_los_dist_visit), 100000, replace=TRUE,
             prob=remaining_los_dist_visit) #sample patient in service
  return(x)
}

######### 2. los if exp 
dis_los2e <- function(){
  sample_init <- function(visit_los_mean, visit_los_sd){  
    seq <- 1:(visit_los_mean + 3*visit_los_sd) #range to sample from
    seq_0.5 <- (seq[1]-0.5):(rev(seq)[1]+0.5) #sampling chance of random pt having los of..
    los_prob_0.5 <- pexp(seq_0.5, rate=1/visit_los_mean, lower.tail = TRUE, log.p = FALSE)
    #los_prob_0.5 <- pnorm(seq_0.5, mean=visit_los_mean, sd=visit_los_sd, log=FALSE) #prob dist
    los_prob_between <- zoo::rollapply(los_prob_0.5, width=2, FUN=function(q) q[2]-q[1])
    los_prob_between_sum <- sum(los_prob_between) #sum for normalising
    normalised_los_prob_between <- los_prob_between/los_prob_between_sum
    a_prob <- normalised_los_prob_between * seq
    adj_mean <- sum(los_prob_between * seq)
    share <- a_prob/adj_mean #los of those in service probabilities
    seq_div_2 <- seq/2 #
    newmean <- sum(share*seq_div_2)
    sample_reps(share, seq)
  }
  remaining_los_dist_visit <- sample_init(avg_LOS, LOS_sd)
  x<- sample(1:length(remaining_los_dist_visit), 100000, replace=TRUE,
             prob=remaining_los_dist_visit) #sample patient in service
  return(x)
}


# checking with plots
library(ggplot2)
y<-dis_los2e()
y <- as.data.frame(y)
png(filename="hist_exp.png", width =1200, height=700)
ggplot(y, aes(x=y)) + geom_histogram(binwidth=1)
dev.off()

n<-dis_los2()
n<- as.data.frame(n)
png(filename="hist_norm.png", width =1200, height=700)
ggplot(n, aes(x=n)) + geom_histogram(binwidth=1)
dev.off()
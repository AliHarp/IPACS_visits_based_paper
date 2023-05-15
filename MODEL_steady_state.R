######  all outputs - raw data, summary by run, summary statistics

library(doSNOW)
library(foreach)
library(tidyverse)
library(zoo)
library(beepr)
library(moments)

rm(list=ls())

setwd("~/Dropbox/03_IPAC/P1_paper_scenarios")
# all parameter combination options
SEED <-1
nruns <- 1000
sim_length <- 1000 #days
warmup <-2000 

#####input variables#####
lambda <- 5.5
#ISR <- 2.5
ISR <- 2
endSR <- 1.7372
#endSR <- 1
sd<-0.5
#avg_LOS <- 12.1304
#avg_LOS <- 12.489
avg_LOS <- 12
LOS_sd <- 4
#LOS_sd <- 6
#dist <- paste0("BRN(",avg_LOS,",", LOS_sd,")")
dist <- paste0("BRExp(",avg_LOS,")")

#n_patients = 70 #number of patients s ystem take
#n_slots  <- n_patients * mean(c(ISR, endSR)) #number of visit slots available per day
n_slots <- 125 #hours
n_patients <- n_slots/mean(c(ISR, endSR))

date_ts = data.frame(seq(as.Date(format(Sys.time(), "%Y-%m-%d")), 
                         by = "day", length.out = sim_length))

outputs<-NULL
######calc distributions for visits for P1 paper to compare with analytical method

seq <- seq(0.5, (ISR+4*sd), by=0.5)
seq_0.5 <- seq((seq[1]-0.25), (rev(seq)[1]+0.25),by=0.5) 

IVR_prob_0.5 <- pnorm(seq_0.5, mean=ISR, sd=sd, log=FALSE) #column B
los_prob_between <- zoo::rollapply(IVR_prob_0.5, width=2, FUN=function(z) z[2]-z[1]) #column C
los_prob_between_sum <- sum(los_prob_between) #sum for normalising
normalised_los_prob_between <- los_prob_between/los_prob_between_sum #column D
normalised_los_prob_between <- normalised_los_prob_between[1:(length(normalised_los_prob_between)-1)] #remove last digit = IVR
seqshort <- seq[1:length(seq)-1]
adj_mean <- sum(normalised_los_prob_between * seqshort) #mean IVR

endSR_prob_0.5 <- pnorm(seq_0.5, mean=endSR, sd=sd, log=FALSE) #column B
los_prob_between_end <- zoo::rollapply(endSR_prob_0.5, width=2, FUN=function(z) z[2]-z[1]) #column C
los_prob_between_sum_end <- sum(los_prob_between_end) #sum for normalising
normalised_los_prob_between_end <- los_prob_between_end/los_prob_between_sum_end #column D
normalised_los_prob_between_end <- normalised_los_prob_between_end[1:(length(normalised_los_prob_between_end)-1)] #remove last digit
seqshort <- seq[1:length(seq)-1]
adj_mean_end <- sum(normalised_los_prob_between_end * seqshort) #mean endVR before end<=initial vr

#sample IVR and a FVR, given FVR<=IVR
dis_slots <- function(){
  x<- sample(seqshort, 1, replace=TRUE, prob=normalised_los_prob_between)
  y<-sample(seqshort, 1, replace=TRUE, prob=normalised_los_prob_between_end)
  if (x<=y){
    y<-x
    return(c(x,y))
  }
  else {return(c(x,y))}
}

#1 million replications for calc prop distribution of endSR given constraints
visits = lapply(seq_len(100), FUN = function(x)  dis_slots())
# 
# #bind into two columsn: ivr and fvr
visitsdf <- do.call(rbind, visits)
colnames(visitsdf) <- c("initials", "ends")
visitsdf<-as.data.frame(visitsdf)
# 
# #check column means and compare with theoretical means in DW excel file
colMeans(visitsdf[sapply(visitsdf, is.numeric)])

#Calculate distribution for FVR to use in simulation
#for IVR, use normalised_los_prob_between
vis <- data.frame(visitsdf) %>% group_by(ends) %>% dplyr::summarise(n=n()) %>% mutate(Freq=n/sum(n))
visfreqEND <- vis$Freq #endSR


    ####runs##################################
    cl<-parallel::makeCluster(17) #, setup_strategy="sequential")
    registerDoSNOW(cl)
    RESULTS<-foreach(run=1:nruns,.combine="rbind") %dopar% {
      set.seed(nruns*(SEED-1)+run)
      
      #distributions
      #arrival distribution
      dis_arrival <- function(){
        x<- round(rpois(1,lambda=lambda))
        if (x<=0){ #if x=0 or negative, then x = 0, no one arrives that day
          x <-0
          return(x)
        } else {return(x)}
      }
      
      ###### SELECT APPROPRIATE LOS DISTRIBUTION
      # #LOS distribution - NOT ACCEPTING 0-DAY lOS
      #   dis_los <- function(){
      #   x<- round(rnorm(1,mean=avg_LOS,sd=LOS_sd))
      #   while (x<=0){
      #     x <-round(rnorm(1,mean=avg_LOS,sd=LOS_sd))
      #   }
      #   return(x)
      # }


      #exponential LOS distribution not accepting 0-DAY los
      dis_los <- function(){
        x<- round(rexp(1, 1/avg_LOS))
       while (x<=0){
        x <-round(rexp(1, 1/avg_LOS))
       }
        return(x)
      }
        
     # ISR distribution
      dis_init_slots <- function(){
        x<- sample(seqshort, 1, replace=TRUE, prob=normalised_los_prob_between) #ISR freq distribution
        return(x)
      }

      # endSR distribution
      dis_end_slots <- function(){
        x<- sample(seqshort, 1, replace=TRUE, prob=normalised_los_prob_between_end) #FVR freq distribution
        return(x)
      }
      # #end SR distribution 
      # dis_end_slots <- function(){
      #   x<- sample(seqshort[1:length(normalised_los_prob_between_end)], 1, replace=TRUE, 
      #              prob=normalised_los_prob_between_end)
      #   return(x)
      # }
      
      #####output variables#####
      ent_sys <- 0 # number of entities that entered the system
      left_sys <-0 # number of entities that left the system
      
      #output after warm up period
      output<-data.frame(RUNX=integer(sim_length), #run number x
                         day= integer(sim_length), #output per day
                         q_length = integer(sim_length), #number of patients in the queue
                         patients_in_service = numeric(sim_length),
                         res_used=numeric(sim_length), #used slots
                         res_idle=numeric(sim_length), #idle slots
                         in_sys=numeric(sim_length))#number of patinets in the system
      

      visits<-data.frame(visits=list(sim_length))
      #####creating necessary data structures#####
      
      #patient list
      patients<-data.frame(id=integer((sim_length+warmup)*2),            #patient id
                           los=integer((sim_length+warmup)*2),           #length of stay
                           arrival_time =integer((sim_length+warmup)*2), # day in the simulation the entity arrived
                           start_service=integer((sim_length+warmup)*2), # day actual service started
                           end_service=integer((sim_length+warmup)*2),   # day service ended
                           wait_time=integer((sim_length+warmup)*2),     # number of days spent in the queue
                           exit=logical((sim_length+warmup)*2),          # boolean variable, TRUE if the entity has left the systemvi
                           stringsAsFactors=FALSE)
      
     
       npat<-0 #initialising counter for patients dataframe
      
      #list with required visit vectors for each patient
      req_visits <- list()

      #resources
      resources <- matrix(nrow=(sim_length+warmup)*10, ncol = 1) #times 2 to make the calculations for resources work at the end of the simulation
      
      resources[,] <- n_slots
      
      #vector for storing waiting time, kept for each patient who left the system
      waittime_vec <- data.frame(RUNX=integer(),
                                 start_service= integer(),
                                 waittime = integer(),
                                 stringsAsFactors=FALSE)
      
      #####simulation#####
      id<-0

      for (t in 1:(sim_length+warmup)) {
        #arrivals to service
        narr<-dis_arrival()
        if(narr>0){
          ent_sys <- ent_sys + narr
          
          #for each arrived patient
          for (j in 1:narr) {
            id<-id+1
            npat<-npat+1
            los<- dis_los()
            arrival_time <- t
            exit <-FALSE
            patients[npat, ] <- c(id,los,arrival_time,NA, NA, 0, exit)
            
            #initial slots and creating required visits vector
            init_slots <- dis_init_slots()
            end_slots <- dis_end_slots()
            if(end_slots>init_slots){
              end_slots = init_slots
            }
            #to accept 0-day LoS
            visit_vector <- function(){
              v <- (seq(init_slots,end_slots,length.out = los)) 
              if (is.integer(v) && length(v) == 0L){
                v <- 0}
              return(v)}
            visit_vector <- visit_vector()
            #remove rounding bias
            for (i in 1:length(visit_vector)){
              sam <- sample(c(-0.01, 0.01), size=1, replace=TRUE)
              visit_vector[[i]] <- visit_vector[i]+sam
              visit_vector[i] <- round(visit_vector[i]/0.5)*0.5
            }
            req_visits[[id]] <- visit_vector
            
            #planning service, checking resources
            tt<-t #temporary t for incrementing when no resources available
  
            while (is.na(patients$start_service[npat])==TRUE){
              if (all((resources[((tt):((tt)+patients$los[npat]-1)),]>= req_visits[[id]])==TRUE)){
                patients$start_service[npat] <- tt
                patients$end_service[npat] <- patients$start_service[npat]+(patients$los[npat]-1)
                
                #decrease capacity
                resources[((tt):((tt)+patients$los[npat]-1)),] <- resources[((tt):((tt)+patients$los[npat]-1)),] - req_visits[[id]]
              } else {
                tt<-tt+1 #if no sufficient resources, check for starting on the next day
              } 
            }          
          }
        }
        
        #increase wait time for patients in the queue
        in_q<-which((patients$start_service>t)&(patients$id>0))
        if (length(in_q)>0){
          patients[in_q,6]<- patients[in_q,6]+1
        }
        
        #recording output from the day warm up period has finished
        if (t>warmup){ #only start recording after the warm up period
          #visits[t-warmup,]<-visits=req_visits 
          if (npat>0 & nrow(waittime_vec)>0) {
            output[t-warmup, ]<- c(RUNX=run, 
                                   day= t,
                                   q_length = length(in_q),
                                   patients_in_service=(n_slots-(resources[t,]))/(mean(c(ISR, endSR))),
                                   res_used= 1- (resources[t,]/n_slots),
                                   res_idle= resources[t,]/n_slots,
                                   in_sys = (ent_sys - left_sys))
          } else if (npat>0 & nrow(waittime_vec)==0) {
            output[t-warmup, ]<- c(RUNX=run,
                                   day= t,
                                   q_length = length(in_q),
                                   patients_in_service=(n_slots-(resources[t,]))/(mean(c(ISR, endSR))),
                                   res_used= 1- (resources[t,]/n_slots),
                                   res_idle= resources[t,]/n_slots,
                                   in_sys = (ent_sys - left_sys))
          } else {
            output[t-warmup, ]<- c(RUNX=run,
                                   day= t,
                                   q_length = length(in_q),
                                   patients_in_service=(n_slots-(resources[t,]))/(mean(c(ISR, endSR))),
                                   res_used= 1- (resources[t,]/n_slots),
                                   res_idle= resources[t,]/n_slots,
                                   in_sys = (ent_sys - left_sys))
          }
        }
        
        #remove patients whose service has ended from the patients table
        remove <- which(patients$end_service==t)
        if(length(remove)>0){
          if(t>=warmup){
            df<-data.frame(RUNX = run, start_service= patients$start_service[remove], waittime= patients[remove,6])
            waittime_vec <- rbind(waittime_vec,df) #keeping waiting time
          }  
          patients <- patients[-remove,] #remove from patient list
          npat<- npat - length(remove)
          left_sys <- left_sys + length(remove)
        }
        
      }
      
      list<-list(output, resources, waittime_vec)
      
      return(list)
      
    }
    stopCluster(cl)
    ###############################################
    
    #creating dataframe for summary info
    summary <- data.frame(LOS = integer(nruns),
                          ISR = integer(nruns),
                          nruns = integer(nruns),
                          sim_length = integer(nruns),
                          warm_up=integer(nruns),
                          capacity = integer(nruns),
                          mean_wait= numeric(nruns),
                          q_length = numeric(nruns),
                          patients_in_service = numeric(nruns),
                          sd_queue = numeric(nruns),
                          skew_queue = numeric(nruns),
                          var_queue = numeric(nruns),
                          res_used= numeric(nruns),
                          res_idle= numeric(nruns),
                          in_sys = numeric(nruns),
                          prop_delay = numeric(nruns))
    
    
    #splitting up RESULTS list in 2
    output<-RESULTS[,1]
    out<-do.call(rbind, output)
    #combining in one dataframe
    
    resources<-RESULTS[,2]
    res<-do.call(cbind, resources) 
    colnames(res)<- c(1:nruns)
    
    waittimes <- RESULTS[,3]
    wait<-do.call(rbind, waittimes)
    
   
    #summary of all runs
    for (k in 1:nruns){ 
      r.out <- which(out[,1]==k)
      k.wait <- which(wait[,1]==k)
      summary[k,]<- c(LOS = avg_LOS,
                      ISR = ISR,
                      nruns = nruns,
                      sim_length = sim_length,
                      warm_up=warmup,
                      capacity = n_slots,
                      mean_wait= round(mean(wait$waittime[k.wait]),8),
                      q_length = round(mean(out$q_length[r.out]),8),
                      patients_in_service = round(mean(out$patients_in_service[r.out]),8),
                      sd_queue = round(sd(out$q_length[r.out]),8),
                      skew_queue = skewness(out$q_length[r.out]),
                      var_queue = var(out$q_length[r.out]),
                      res_used= round(mean(out$res_used[r.out]),8),
                      res_idle= round(mean(out$res_idle[r.out]),8),
                      in_sys= round(mean(out$in_sys[r.out]),8),
                      prop_delay = length(which(wait$waittime[k.wait]!=0))/length(wait$waittime[k.wait]))
    }
    ########### output ##############
    
    CI_z <- function (x, ci = 0.95)
    {
      `%>%` <- magrittr::`%>%`
      standard_deviation <- sd(x)
      sample_size <- length(x)
      Margin_Error <- abs(qnorm((1-ci)/2))* standard_deviation/sqrt(sample_size)
      df_out <- data.frame( sample_size=length(x), Mean=mean(x), sd=sd(x),
                            Margin_Error=Margin_Error,
                            'CI lower limit'=(mean(x) - Margin_Error),
                            'CI Upper limit'=(mean(x) + Margin_Error)) %>%
        tidyr::pivot_longer(names_to = "Measurements", values_to ="values", 1:6 )
      return(df_out)
    }
    
    q<- CI_z(summary$q_length)
    q$values <- round(q$values,4)
    w<-CI_z(summary$mean_wait)
    w$values <- round(w$values,4)
    u<-CI_z(summary$res_used)
    u$values <- round(u$values, 4)
    d<-CI_z(summary$prop_delay)
    d$values <- round(d$values, 4)
    
    var <- var(out$q_length)
    sd <- sqrt(var)
    sk <- skewness(out$q_length)

    
    output2<- cbind(q,w[,2],u[,2],d[,2], sk, sd, var, dist, lambda, ISR, endSR)
    colnames<- cbind("measure", "E(q)", "E(w)", "Resource use","P(delayed)",
                     "skewness(q)","sd(q)","var(q)", "dist", "lambda", "IVR", "FVR")
    colnames(output2)<-colnames
    
    #summary by run
    write.csv(summary, paste0("csv_output/summary_by_run", dist, lambda, " ISR=", ISR, 
                             " endSR=", endSR, ".csv"), row.names=FALSE)
    #summary statistics
    write.csv(output2, paste0("csv_output/summary/dist_summary_", dist, lambda, "ISR= ", ISR,
                              " endSR ", endSR, ".csv"), row.names=FALSE)
    
    #raw output runs and runtime
    write.csv(out, paste0("csv_output/raw_output_", dist, lambda, "ISR= ",ISR,", endSR ", endSR, ".csv"))
    
    #write.csv(out, paste0("Example resources available per day, lambda=", lambda, ".csv", row.names=FALSE))
    #write.csv(wait, paste0("sample_waits, lambda=", lambda, ".csv"), row.names=FALSE)
    #write.csv(summary, paste0("csv_output/summary_by_run, BRExp(12)_L.csv"))
    
    
    beep(sound=100)
    

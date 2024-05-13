##install packages and libraries

install.packages("pscl")
install.packages("jtools")
install.packages("ggpubr")

require(ggplot2)
require(pscl)

library(nnet)
library(MASS)
library(jtools)
library(dplyr)


###Import Data###

##import control cohort data
controldata <- read.table("C:/Users/em286502/OneDrive - NHS Wales/Documents/PASS/ControlNoPassFreq.txt",sep="\t",col.names=c("Recruitment Cohort",              "GMC Number",             "MSRA score",  "1st 6m Extension",        "2nd 6m Extension",      "Exceptional 3rd Extension",              "Extension Count",        "Extensions Y/N",           "AKT 1st attempt",              "AKT 2nd attempt",        "AKT 3rd attempt",         "AKT 4th attempt",              "AKT excpetional 5th attempt", "AKT Attempt Count",    "AKT Overall Pass/Fail",         "CSA/RCA 1st attempt",             "CSA/RCA 2nd attempt",            "CSA/RCA 3rd attempt",             "CSA/RCA 4th attempt",              "CSA/RCA exceptional 5th attempt",    "CSA/RCA Attempt Count", "CSA/RCA Overall Pass/Fail",   "Released from training",              "Successfully completed training",        "PASS", "Notes", "ExtensionFreq"))

##import PASS cohort data
PASSdata <- read.table("C:/Users/em286502/OneDrive - NHS Wales/Documents/PASS/PASSCohortMSRANew.txt",sep="\t",col.names=c("Cohort",              "GMC Number",             "MSRA score",  "MSRA Grouped"))


###Data Cleaning for Control Cohort### 

#remove 0s from AKT attempt count, CSA attempt count, MSRA Score
ctrlcohort <- subset(controldata, controldata$MSRA.scwhatore != 0)
ctrlcohort <- subset(ctrlcohort, ctrlcohort$AKT.Attempt.Count != 0)
ctrlcohort <- subset(ctrlcohort, ctrlcohort$CSA.RCA.Attempt.Count != 0)


###Data Cleaning for PASS Cohort###

##remove 0s from MSRA score
PASSdata <- subset(PASSdata, PASSdata$MSRA.score != 0)



###Data Binning###

#PASS cohort bins for graph
PASSdata$binmidpoint <- round(PASSdata$MSRA.score, -1)

#control cohort bins for graph
ctrlcohort$binmidpoint <- round(ctrlcohort$MSRA.score, -1)



###Pearson Correlation Coefficient###
pearsons <- cor(ctrlcohort$MSRA.score, ctrlcohort$Extension.Count, method="pearson")

#results show a moderate correlation
pearsons




###Check the Mean and Variance for the Count Data used in Poisson Model###

dispersion_test <- function(x) 
{
  res <- 1-2 * abs((1 - pchisq((sum((x - mean(x))^2)/mean(x)), length(x) - 1))-0.5)
  
  cat("Dispersion test of count data:\n",
      length(x), " data points.\n",
      "Mean: ",mean(x),"\n",
      "Variance: ",var(x),"\n",
      "Probability of being drawn from Poisson distribution: ", 
      round(res, 3),"\n", sep = "")
  
  invisible(res) }
#results
dispersion_test(ctrlcohort$Extension.Countl)



  
  ###Poisson Model for Control Data###
  
  poissonmodel <- glm(data = ctrlcohort, Extension.Count ~ MSRA.score, family = poisson)
   
  ##gives p values
  
  summary(poissonmodel)
  
  #log likelihood
  loglik_poisson <- logLik(poissonmodel)
  loglik_poisson
  
 
  
  ###stores results of prediction###
  PASSdata$Poisson_Result <- predict(poissonmodel, PASSdata)
  
  PASSdata$Poisson_Prediction <- round(exp(PASSdata$Poisson_Result),0)
  
  
  
  
  #CREATING HISTOGRAM FOR PASS COHORT WITH PREDICTIONS OBTAINED
  
 
  #Gives the frequency of the predicted value per MSRA score (PASS Cohort)
  PASSdata$Poisson_Pred_Extension_Freq <- sapply(1:nrow(PASSdata), function(i){
    PASSdata %>% filter(MSRA.score == PASSdata$MSRA.score[i], Poisson_Prediction == PASSdata$Poisson_Prediction[i]) %>% nrow
  })

  #PASS cohort predicted extension frequency per bin
  PASSdata$Poisson_Pred_Extension_Freq_Bin <- sapply(1:nrow(PASSdata), function(i){
    PASSdata %>% filter(binmidpoint == PASSdata$binmidpoint[i], Poisson_Prediction == PASSdata$Poisson_Prediction[i]) %>% nrow
  })
  
  
  ###Creating the histogram of PASS cohort and predicted extension count###
  PASS_poisson_predicted_histogram <- ggplot(PASSdata, aes(x = binmidpoint)) + 
    scale_x_continuous(breaks = seq(375,485, by = 10)) +
    geom_histogram(binwidth = 10, fill = "#489cc9", colour = "#325083") +
    geom_point(data = PASSdata, aes(y = Poisson_Pred_Extension_Freq_Bin, colour = as.factor(Poisson_Prediction)), size = 3) + 
    scale_colour_manual("Number of Predicted Extensions", values = c("#575756", "#afca0b", "#e52e7d","#ed8900")) +
    scale_y_continuous(breaks = seq(0,32, by = 2), name = "MSRA Score Frequency", sec.axis = sec_axis(~.*1, breaks = seq(0,32, by = 2), name = "Predicted Extensions Frequency")) +
    theme_nice(legend.pos = "bottom") +
    labs(title = "PASS Cohort: MSRA Scores and Predicted Number of Extensions", x = "MSRA Score")
 
   #Output graph 
  PASS_poisson_predicted_histogram
  
  
  
  
  
#CREATING HISTOGRAM FOR CONTROL COHORT 
  

#control cohort bin for extension frequency
ctrlcohort$Extension_Freq_Bin <- sapply(1:nrow(ctrlcohort), function(i){
  ctrlcohort %>% filter(binmidpoint == ctrlcohort$binmidpoint[i], Extension.Count == ctrlcohort$Extension.Count[i]) %>% nrow
})

###creating the histogram of control cohort's extensions by MSRA###
ctrlcohort_extension_histogram <- ggplot(ctrlcohort, aes(x = binmidpoint)) + 
  scale_x_continuous(breaks = seq(375,700, by = 10)) +
  geom_histogram(binwidth = 10, fill = "#489cc9", colour = "#325083") +
  geom_point(data = ctrlcohort, aes(y = Extension_Freq_Bin, colour = as.factor(Extension.Count)), size = 3) + 
  scale_colour_manual("Number of Extensions", values = c("#575756", "#afca0b", "#e52e7d","#ed8900", "#78368c","#325083","#00917f")) +
  scale_y_continuous(breaks = seq(0,50, by = 2), name = "MSRA Score Frequency", sec.axis = sec_axis(~.*1, breaks = seq(0,50, by = 2), name = "Extensions Frequency")) +
  theme_nice(legend.pos = "bottom") +
  labs(title = "Control Cohort: MSRA Scores and Number of Extensions", x = "MSRA Score")

#output graph: use compare with PASS prediction histogram
ctrlcohort_extension_histogram


  
  
  ##save results to csv##
  write.csv(PASSdata,"C:/Users/em286502/OneDrive - NHS Wales/Documents/PASS/PASSPredictResults.csv")
  
  
  
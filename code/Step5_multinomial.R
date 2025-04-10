library(data.table)
library(bigsnpr)
library(brglm2)
library(dplyr)


load("../data//Clinical_maternal_babysexnotcomplete_10_20.RData")
colnames(meta_info)[1]<-"IID"
snp<-fread("../data/PE_maternal_control_region_filtered.raw")
full_data<-merge(meta_info,snp,by="IID")
full_data<-full_data[,-c(10:14)]

listname<-make.names(colnames(full_data))
colnames(full_data)<-listname

snplist<-colnames(full_data)[10:length(listname)]

length(snplist)

full_data$race_recode<-relevel(full_data$race_recode,ref=3)




control <- brglmControl(maxit = 200, epsilon = 1e-07, slowit = 0.5)
p_values<-list()
count<-lst()

for (i in 1:length(snplist)){
  snp_name<-snplist[i]
  formula<-as.formula(paste("race_recode ~",snp_name,"+MomAge+BMI",collapse=""))
  model<-tryCatch({brmultinom(formula,data=full_data,control=control)},error=function(e){
    message("Error fitting ",snp_name)
    cat("An error occurred: ", conditionMessage(e),"\n")
    return(NULL)
  })
  if(is.null(model))
    next
  coeff<-summary(model)$coefficients[,2]
  std_errors<-summary(model)$standard.errors[,2]
  z_score<-coeff/std_errors
  p_value<-2*(1-pnorm(abs(z_score)))
  #Get the p-value
  p_values[[snp_name]]<-p_value
  #Get the count
  snp_data<-full_data[c('race_recode',snp_name)]
  df_summary<-snp_data %>%
    group_by(race_recode) %>%
    summarize(total_count=sum(across(all_of(snp_name)),na.rm=TRUE))
  count[[snp_name]]<-df_summary
  print(snp_name)
  cat("\n")
  print(p_value)
  cat("\n")
  print(df_summary)
  cat("\n")
}

save(p_values,count,file="../output/multinomresult_Momage_BMI.RData")
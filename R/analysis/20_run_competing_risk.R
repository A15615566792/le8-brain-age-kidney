args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2)stop("Usage: Rscript 20_run_competing_risk.R full.csv out_dir")
full_path<-args[1];out_dir<-args[2];dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages(library(survival))
z<-function(x)as.numeric((as.numeric(x)-mean(as.numeric(x),na.rm=TRUE))/sd(as.numeric(x),na.rm=TRUE))
eth4<-function(x){x<-tolower(trimws(as.character(x)));w<-x%in%c("british","irish","white","any other white background");a<-x%in%c("indian","pakistani","bangladeshi","chinese","asian or asian british","any other asian background");b<-x%in%c("caribbean","african","black or black british","any other black background");factor(ifelse(w,"White",ifelse(a,"Asian",ifelse(b,"Black","Other/unknown"))),levels=c("White","Asian","Black","Other/unknown"))}
tidy_cox<-function(f,t){b<-coef(f)[t];se<-sqrt(vcov(f)[t,t]);data.frame(outcome="Incident_CKD",model="Fine_Gray_death_competing",term=t,beta=unname(b),se=unname(se),HR=unname(exp(b)),CI_low=unname(exp(b-1.96*se)),CI_high=unname(exp(b+1.96*se)),p=unname(2*pnorm(-abs(b/se))),n=f$n,events=f$nevent)}
f<-read.csv(full_path,check.names=FALSE,na.strings=c("","NA"))
f$bd<-as.Date(substr(f$baseline_date,1,10));f$t0<-as.Date(substr(f$time0_full,1,10));f$age_webq<-as.numeric(f$Age_at_recruitment)+as.numeric(f$t0-f$bd)/365.25
f$sexf<-relevel(factor(f$Sex),ref="Female");f$eth4<-eth4(f$Ethnic_0);f$tdi<-as.numeric(f$TDI_at_recruitment);f$edu<-as.numeric(f$area_education_deprivation_z);f$em<-factor(ifelse(is.na(f$edu),"Missing","Observed"));f$ei<-ifelse(is.na(f$edu),0,f$edu);f$ctr<-factor(f$baseline_center);f$en<-z(f$energy_kcal_primary);f$rnd<-as.numeric(f$webq_n_dates);f$LE8<-as.numeric(f$LE8)/10
f$status<-factor(ifelse(f$kidney_event==1,"CKD",ifelse(f$death_competing==1,"Death","Censor")),levels=c("Censor","CKD","Death"))
vars<-c("followup_years","status","LE8","age_webq","sexf","eth4","tdi","ei","em","ctr","en","rnd")
d<-f[complete.cases(f[,vars]),vars];d$id_num<-seq_len(nrow(d))
fg<-finegray(Surv(followup_years,status)~LE8+age_webq+sexf+eth4+tdi+ei+em+ctr+en+rnd+id_num,data=d,etype="CKD",id=id_num)
fit<-coxph(Surv(fgstart,fgstop,fgstatus)~LE8+age_webq+sexf+eth4+tdi+ei+em+ctr+en+rnd+cluster(id_num),data=fg,weights=fgwt)
res<-tidy_cox(fit,"LE8");res$source_n<-nrow(d);res$expanded_n<-nrow(fg)
write.csv(res,file.path(out_dir,"Table_P3_FineGray_total_LE8.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")

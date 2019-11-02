#Data Collector & Analysis

#Library Download --------------------------------------------------------------------

rm(list=ls())		#Clear Memory

library(dplyr)		#Data Management
library(e1071)		#Skewness Formula
library(ggplot2)		#Graphical Presentation
library(gridExtra)	#Graphical Presentation

#Setting -----------------------------------------------------------------------------

DirInK1000<-"C:/Users/Andrea/Desktop/PhD/PORTFOLIO OPTIMIZATION/OUTPUT/Trials_1000"
DirInK1500<-"C:/Users/Andrea/Desktop/PhD/PORTFOLIO OPTIMIZATION/OUTPUT/Trials_1500"
DirInK2000<-"C:/Users/Andrea/Desktop/PhD/PORTFOLIO OPTIMIZATION/OUTPUT/Trials_2000"
DirOut<-"C:/Users/Andrea/Desktop/PhD/PORTFOLIO OPTIMIZATION/OUTPUT"

Homogenity_Scaling<-4.5

#Function ----------------------------------------------------------------------------

PreProcessing<-function(DirIn, Trials){

	#Declare Working Directory
	setwd(DirIn)

	#Single Files Reference
	File_Name<-list.files(pattern="*.csv")
	File_Length<-length(File_Name)

	#Data Container
	Run<-NULL

	#Data Download and Standardization
	for(i in 1:File_Length){

		#Download csv file
		Temp<-read.csv(file=File_Name[i], header=TRUE, sep=",", row.names= NULL) %>%
			as.data.frame()

		#Single file Processing
		Temp_Name<-Temp[,2]

		Temp_Mtrx<-data.matrix(Temp[,-c(1,2)], rownames.force = NA) %>% 
			     	t

		row.names(Temp_Mtrx)<-NULL

		Temp_db<-as.data.frame(Temp_Mtrx)
		names(Temp_db)<-Temp_Name

		#Standardize Performance Measures
		Temp_db$UMCS_CardinalityStd<-(Temp_db$UMCS_Cardinality)/max(Temp_db$UMCS_Cardinality,Temp_db$DEMO_Cardinality)
		Temp_db$DEMO_CardinalityStd<-(Temp_db$DEMO_Cardinality)/max(Temp_db$UMCS_Cardinality,Temp_db$DEMO_Cardinality)

		Temp_db$UMCS_SpacingStd<-(Temp_db$UMCS_Spacing)/max(Temp_db$UMCS_Spacing,Temp_db$DEMO_Spacing)
		Temp_db$DEMO_SpacingStd<-(Temp_db$DEMO_Spacing)/max(Temp_db$UMCS_Spacing,Temp_db$DEMO_Spacing)

		Temp_db$UMCS_SpreadStd<-(Temp_db$UMCS_Spread)/max(Temp_db$UMCS_Spread,Temp_db$DEMO_Spread)
		Temp_db$DEMO_SpreadStd<-(Temp_db$DEMO_Spread)/max(Temp_db$UMCS_Spread,Temp_db$DEMO_Spread)

		Temp_db$UMCS_RangeStd<-(Temp_db$UMCS_Range)/max(Temp_db$UMCS_Range,Temp_db$DEMO_Range)
		Temp_db$DEMO_RangeStd<-(Temp_db$DEMO_Range)/max(Temp_db$UMCS_Range,Temp_db$DEMO_Range)

		Temp_db$UMCS_DominanceStd<-(Temp_db$UMCS_Dominance)/(Temp_db$DEMO_Dominance+Temp_db$UMCS_Dominance)
		Temp_db$DEMO_DominanceStd<-(Temp_db$DEMO_Dominance)/(Temp_db$DEMO_Dominance+Temp_db$UMCS_Dominance)

		#Add Run Reference
		Temp_db$Portfolio<-substr(File_Name[i],8,10)
		Temp_db$Market<-substr(File_Name[i],15,17)	
		Temp_db$Run<-i
		Temp_db$Trials<-Trials

		#Combine Runs
		Run<-rbind(Run,Temp_db)
	}

	#Prepare Output
	return(Run)
}

Summary<-function(x){		
	Stat<-c(min(x),quantile(x,0.05),quantile(x,0.25),quantile(x,0.50),quantile(x,0.75),quantile(x,0.95),max(x),
		  mean(x),sd(x),skewness(x),sum(x>0)/length(x))
}

Summary_DB<-function(db, Trials, Type){

	out<-as.data.frame(apply(db[,6:11], 2, Summary))
	rownames(out)<-c("Min","q_0.05","q_0.25","q_0.50","q_0.75","q_0.95","Max","Mean","Sd","Skew","Prob(DEMO>UMCS)")

	out$Trials<-Trials
	out$Type<-Type

	return(out)
}

#Analysis ----------------------------------------------------------------------------

#Download Single Run Results
Run_K1000<-PreProcessing(DirInK1000, Trials = 1000)
Run_K1500<-PreProcessing(DirInK1500, Trials = 1500)
Run_K2000<-PreProcessing(DirInK2000, Trials = 2000)

#Run Database
Run<-rbind(Run_K1000, Run_K1500, Run_K2000)

#Performance Database
DB<-data.frame(Trials = Run$Trials,
		   m = Run$m,
		   Prt_Dimension = Run$m,
		   Prt_Homogenity = (Run$Freq_Sd+Run$Sev_Mean_Sd+Run$Sev_Sd_Sd)/Homogenity_Scaling,
		   Market = Run$Market,
		   P_Hyper = Run$DEMO_Hyper - Run$UMCS_Hyper,
		   P_CardinalityStd= Run$DEMO_CardinalityStd- Run$UMCS_CardinalityStd,
		   P_SpacingStd= Run$DEMO_SpacingStd- Run$UMCS_SpacingStd,
		   P_SpreadStd= Run$DEMO_SpreadStd- Run$UMCS_SpreadStd,
		   P_RangeStd= Run$DEMO_RangeStd- Run$UMCS_RangeStd,
		   P_DominanceStd= Run$DEMO_DominanceStd- Run$UMCS_DominanceStd)

#Summary Statistics
DB_Stat<-rbind(Summary_DB(db=DB[DB$Trials==1000,],Trials=1000,Type="Total"),
		   Summary_DB(db=DB[DB$Trials==1000&DB$Prt_Dimension==100,],Trials=1000,Type="Portfolio Low"),
		   Summary_DB(db=DB[DB$Trials==1000&DB$Prt_Dimension==200,],Trials=1000,Type="Portfolio High"),
		   Summary_DB(db=DB[DB$Trials==1000&DB$Market=="Low",],Trials=1000,Type="Market Low"),
		   Summary_DB(db=DB[DB$Trials==1000&DB$Market=="Med",],Trials=1000,Type="Market Medium"),
		   Summary_DB(db=DB[DB$Trials==1000&DB$Market=="Hig",],Trials=1000,Type="Market High"),
		   Summary_DB(db=DB[DB$Trials==1500,],Trials=1500,Type="Total"),
		   Summary_DB(db=DB[DB$Trials==2000,],Trials=2000,Type="Total"))

#Portfolio Homogenity
hist(DB[DB$Trials==1000,"Prt_Homogenity"],xlab ="Standardize Level",main="Portfolio Homogenity")

par(mfrow=c(3,2))
plot(DB[DB$Trials==1000,"Prt_Homogenity"],DB[DB$Trials==1000,"P_Hyper"],xlab="Homogenity Level",ylab="Delta HyperVolume",main="Delta HyperVolume per Homogenity Level")
plot(DB[DB$Trials==1000,"Prt_Homogenity"],DB[DB$Trials==1000,"P_CardinalityStd"],xlab="Homogenity Level",ylab="Delta Cardinality",main="Delta Cardinality per Homogenity Level")
plot(DB[DB$Trials==1000,"Prt_Homogenity"],DB[DB$Trials==1000,"P_SpacingStd"],xlab="Homogenity Level",ylab="Delta Spacing",main="Delta Spacing per Homogenity Level")
plot(DB[DB$Trials==1000,"Prt_Homogenity"],DB[DB$Trials==1000,"P_SpreadSpd"],xlab="Homogenity Level",ylab="Delta Spread",main="Delta Spread per Homogenity Level")
plot(DB[DB$Trials==1000,"Prt_Homogenity"],DB[DB$Trials==1000,"P_RangeStd"],xlab="Homogenity Level",ylab="Delta Range",main="Delta Range per Homogenity Level")
plot(DB[DB$Trials==1000,"Prt_Homogenity"],DB[DB$Trials==1000,"P_DominanceStd"],xlab="Homogenity Level",ylab="Delta Dominance",main="Delta Dominance per Homogenity Level")

#Performance Histograms Complete
hist(DB[DB$Trials==1000,"P_Hyper"],		   main="HyperVolume Delta", xlab="Delta HyperVolume"); abline(v=0)
hist(DB[DB$Trials==1000,"P_CardinalityStd"], main="Cardinality Delta", xlab="Delta Cardinality"); abline(v=0)
hist(DB[DB$Trials==1000,"P_SpacingStd"],     main="Spacing Delta",     xlab="Delta Spacing");     abline(v=0)
hist(DB[DB$Trials==1000,"P_SpreadStd"],      main="Spread Delta",      xlab="Delta Spread");      abline(v=0)
hist(DB[DB$Trials==1000,"P_RangeStd"],       main="Range Delta",       xlab="Delta Range");       abline(v=0)
hist(DB[DB$Trials==1000,"P_DominanceStd"],   main="Dominance Delta",   xlab="Delta Dominance");   abline(v=0)

#Output ------------------------------------------------------------------------------

setwd(DirOut)

write.csv(Run,     file = "Run_Result.csv",      row.names=FALSE)
write.csv(DB,      file = "Run_Performance.csv", row.names=FALSE)
write.csv(DB_Stat, file = "Run_Statistics.csv",  row.names=TRUE)

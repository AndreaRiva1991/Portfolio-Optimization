#PORTFOLIO OPTIMIZATION - EFFICIENT FRONTIER SEARCH

#00.Library --------------------------------------------------------------------------------------

#Clear Memory
rm(list=ls())

#Library Download
library(emoa)		#HyperVolume Evaluation
library(e1071)		#Skewness Formula
library(foreach)		#Parallel Computing
library(doParallel)	#Parallel Computing

#01.Function Base --------------------------------------------------------------------------------

#Simulate Database
SyntheticDB<-function(s, m, F_mn, S_mn, S_sd, F_mn_sd, S_mn_sd, S_sd_sd, Seed){
  
  #Reproducibility
  set.seed(Seed)								
  
  #Portfolio Parameter
  P_Par<-as.data.frame(matrix(0,m,3))					#Create data container for Freq/Sev parameter
  names(P_Par)<-c("Freq_Mean","Sev_Mean","Sev_Sd")
  
  P_Par$Freq_Mean<-abs(F_mn+rnorm(m,0,F_mn_sd))			#Simulate random Frequency distribution
  P_Par$Sev_Mean<-abs(S_mn+rnorm(m,0,S_mn_sd))			#Simulate random Severity - mean
  P_Par$Sev_Sd<-abs(S_sd+rnorm(m,0,S_sd_sd))			#Simulate random Severity - sd
  
  #Full Portfolio - PREMIUM
  P<-exp(P_Par$Sev_Mean+P_Par$Sev_Sd^2/2)*P_Par$Freq_Mean	#Fair Premium
  
  #Full Portfolio - NUMBER OF CLAIMS	
  F<-matrix(rpois(s,P_Par$Freq_Mean),m,s)				#Number of claims
  F_max<-max(F)								#Max Number of Claims in single simulation
  F_bool<-array(0,c(m,s,F_max))					#Number of claims - booleand cube
  for(i in 1:F_max){F_bool[,,i]<-ifelse(F>=i,1,0)}		#Selection
  
  #Full Portfolio - SEVERITY
  S<-rowSums(array(rlnorm(n = m*s*F_max, mean = P_Par$Sev_Mean, sd = P_Par$Sev_Sd),
                   c(m,s,F_max))*F_bool, 
             dims = 2)
  
  #Prepare Output
  out<-list()
  out$P<-P
  out$S<-S
  out$R<-P_Par
  
  return(out)
}

#Market Dynamic
Market<-function(Beta0, Beta1, Beta2, Data){
  
  #GLM Score
  Score<-Beta0+Beta1*Data[1]+Beta2*Data[2]
  
  #Acceptance Probability
  Prob<-if(Data[2]!=0){
    1-exp(Score)/(1+exp(Score))	
  }else{
    0
  }
  
  #Prepare Output
  return(Prob)
}

#Evaluate Portfolio
Ptf_Eval<-function(S, P, V, alfa, Beta0, Beta1, Beta2, Scaling, Seed){
  
  #Data Container
  out<-NULL	
  
  #Reference
  v1<-V[,1]
  v2<-V[,2]
  m<-dim(V)[1]
  s<-dim(S)[2]
  Data<-cbind(P/max(P), v2)
  
  #Reproducibility
  set.seed(Seed)
  
  #Market Dynamic
  Prob<-if(Scaling == FALSE){
    apply(Data, 1, Market, Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2)
  }else{
    rep(1, m)
  }
  
  Acceptance<-t(apply(as.matrix(Prob), 1 , function(y){sample(x = 0:1, size = s, replace = TRUE, prob = c(1-y, y))}))
  
  #Quantile Calculation
  X<-as.numeric(t(v1)%*%(S*Acceptance))
  q<-quantile(X,alfa)
  
  #TVaR Calculation
  out[1]<-(mean(X[X>=q])-mean(X))
  
  #Premium Calculation	
  out[2]<-if(Scaling == FALSE){
    mean((P*v1*(1+v2))%*%Acceptance)
  }else{
    mean((P*v1)%*%Acceptance)
  }
  
  #Retain Calculation
  out[3]<-if(sum(v1)!=0){
    mean(colSums(Acceptance))/sum(v1)
  }else{
    0
  }
  
  #Prepare Output
  return(out)
}

#Calc Frontier
CalcFrontier<-function(Evl){
  
  #Reference
  Evl_n<-nrow(Evl)
  
  #Data Container
  Eff<-rep(0,Evl_n)	
  
  #Efficient Frontier Finding
  for(i in 1:Evl_n){
    
    #Select observation
    Ref_test<-Evl[i,]
    Ref_db<-Evl[-i,]
    
    #Efficency Conditions
    c_rsk<-Ref_db[,1]<=Ref_test[1]
    c_prm<-Ref_db[,2]>=Ref_test[2]
    c_ret<-Ref_db[,3]>=Ref_test[3]
    c_eff<-c_rsk&c_prm&c_ret
    
    #Flag if Efficient
    if(sum(c_eff)==0){Eff[i]<-1}
  }
  
  #Prepare Output
  out<-cbind(Evl,Eff)
  return(out)
}

#Clear Memory
ClearMemory<-function(Evl, Pop){
  
  #Data Container
  Memory_Evl<-Evl
  Memory_Pop<-Pop
  
  #Sort Solution
  c0<-order(Memory_Evl[,1])
  Memory_Evl<-Memory_Evl[c0,]
  Memory_Pop<-Memory_Pop[,c0,]	
  
  #Delate Duplicated Solution
  c1<-!duplicated(Memory_Evl)
  Memory_Evl<-Memory_Evl[c1,]	
  Memory_Pop<-Memory_Pop[,c1,]	
  
  #Check Updated View of Efficient Frontier
  Memory_Eval<-CalcFrontier(Evl = Memory_Evl[,1:3])
  c2<-Memory_Eval[,4]==1
  Memory_Evl<-Memory_Evl[c2,]
  Memory_Pop<-Memory_Pop[,c2,]
  
  #Prepare Output
  out<-list(Memory_Evl, Memory_Pop)
  names(out)<-c("Evl", "Pop")
  return(out)
}

#Random Search
SEARCH_RS<-function(m, S, P, alfa, Iteration, Beta0, Beta1, Beta2, C_Min, C_Max, Seed){
  
  #Data Container
  out<-list()
  out$Evl<-NULL
  out$Pop<-NULL
  Ptf<-array(0, c(m, Iteration, 2))
  
  #Reproducibility
  set.seed(Seed)
  
  #Simulate Portfolio Dimensions
  Ptf_Prob<-runif(n = Iteration, min = 0, max = 1)
  
  #Generate Random Portfolios
  Ptf[,,1]<-apply(as.matrix(Ptf_Prob), 1 , function(y){sample(x = 0:1, size = m, replace = TRUE, prob = c(y, 1-y))})
  Ptf[,,2]<-matrix(runif(n = m*Iteration, min = C_Min, max = C_Max), m, Iteration)*Ptf[,,1]
  
  #Total Portfolio Scaling Factor
  Scaling_Factor<-Ptf_Eval(S = S, P = P, V = matrix(1, m, Iteration), alfa = alfa, Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2, Scaling = TRUE, Seed = Seed)
  
  #Risk/Return/Retain Portfolios Evaluation
  Ptf_Evl<-t(apply(Ptf, 2, Ptf_Eval, S = S, P = P, alfa = alfa, Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2, Scaling = FALSE, Seed = Seed))
  Ptf_Evl[,1]<-Ptf_Evl[,1]/Scaling_Factor[1]
  Ptf_Evl[,2]<-Ptf_Evl[,2]/Scaling_Factor[2]	
  
  #Calculate Observed Efficient Frontier
  Ptf_Eff<-CalcFrontier(Evl = Ptf_Evl)
  c0<-Ptf_Eff[,4]==1
  Eff_Evl<-Ptf_Eff[c0,]
  Eff_Ptf<-Ptf[,c0,]
  
  EF<-ClearMemory(Evl = Eff_Evl, Pop = Eff_Ptf)
  
  #Prepare Output
  out$Evl<-EF$Evl[,1:3]
  out$Pop<-EF$Pop
  return(out)
}

#02.Function Evolutionary ------------------------------------------------------------------------

#Generate Initial Population
EO_Init<-function(m, Pop_N, C_Min, C_Max, Seed){
  
  #Data Container
  out<-array(0, c(m, Pop_N, 2))
  
  #Reproducibility
  set.seed(Seed)
  
  #Simulate Portfolio Dimensions
  Prob<-runif(n = Pop_N, min = 0, max = 1)
  
  #Generate Random Portfolios
  out[,,1]<-apply(as.matrix(Prob), 1 , function(y){sample(x = 0:1, size = m, replace = TRUE, prob = c(y, 1-y))})
  out[,,2]<-matrix(runif(n = m*Pop_N, min = C_Min, max = C_Max), m, Pop_N)*out[,,1]
  
  #Output
  return(out)
}

#Evaluate TVaR/Return
EO_Eval<-function(m, S, P, Pop, alfa, Beta0, Beta1, Beta2, Seed){
  
  #Data Container
  out<-matrix(0, dim(Pop)[2], 5)
  colnames(out)<-c("Risk", "Return", "Retain", "Pareto", "Crowding")
  
  #Reproducibility
  set.seed(Seed)
  
  #Total Portfolio Scaling Factor	
  Scaling_Factor<-Ptf_Eval(S = S, P = P, V = matrix(1, m, Iteration), alfa = alfa, Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2, Scaling = TRUE, Seed = Seed)
  
  #Risk/Return Path Evaluation
  Ptf_Evl<-t(apply(Pop, 2, Ptf_Eval, S = S, P = P, alfa = alfa, Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2, Scaling = FALSE, Seed = Seed))
  out[,1]<-Ptf_Evl[,1]/Scaling_Factor[1]
  out[,2]<-Ptf_Evl[,2]/Scaling_Factor[2]	
  out[,3]<-Ptf_Evl[,3]
  
  #Prepare Output
  return(out)
}

#Pareto Sorting
ParetoSort<-function(Pop_Evl){
  
  #Reference
  Temp_Evl<-Pop_Evl$Evl
  Temp_Pop<-Pop_Evl$Pop	
  Temp_N<-nrow(Temp_Evl)
  
  options(warn=-1)
  
  #Data Container
  out<-list()
  
  #Pareto Sorting
  c0<-Temp_Evl[,4]==0
  i<-1
  
  while(sum(c0)!=0){
    
    #Turn Off Points Already Analyzed
    Temp<-Temp_Evl[,1:3]
    Temp[!c0,]<-0
    
    if(sum(Temp)!=0){
      
      #Evaluate Frontier
      Temp_EF<-CalcFrontier(Evl = Temp)
      
      #Pareto Order
      Temp_Evl[c0&Temp_EF[,4]==1,4]<-i
      
    }else{
      #Pareto Order
      Temp_Evl[Temp_Evl[,4]==0,4]<-max(Temp_Evl[,4])+1	
    }
    
    #Update Conditions
    c0<-Temp_Evl[,4]==0
    i<-i+1
  }
  
  #Prepare Output
  out<-list(Temp_Evl, Temp_Pop)
  names(out)<-c("Evl", "Pop")
  return(out)
}

#Crowding Distance
CrowdingDist<-function(Pop_Evl){
  
  #Reference
  Temp_Evl<-Pop_Evl$Evl
  Temp_Pop<-Pop_Evl$Pop	
  Temp_N<-nrow(Temp_Evl)
  
  options(warn=-1)
  
  #Data Container
  out<-list()
  
  #Crowding Distance
  Pareto_EF<-max(Temp_Evl[,4])
  
  for(i in 1:Pareto_EF){
    
    #Check Condition
    c0<-Temp_Evl[,4]==i
    
    if(sum(c0)>1){
      
      #Filter Pareto Level
      Temp<-Temp_Evl[c0,]
      
      #Reference
      Temp_n<-nrow(Temp)		
      c1<-order(Temp[,1])
      c2<-order(Temp[,2])
      c3<-order(Temp[,3])
      
      #Data Container		
      Dist1<-rep(0,Temp_n)
      Dist2<-rep(0,Temp_n)
      Dist3<-rep(0,Temp_n)
      
      #Crowding Distance - Risk
      Dist1[2:(Temp_n-1)]<-Temp[c1,1][3:Temp_n]-Temp[c1,1][1:(Temp_n-2)]
      Dist1[Dist1==0]<-max(Dist1[Dist1!=0])
      
      #Crowding Ditsance - Return
      Dist2[2:(Temp_n-1)]<-Temp[c2,2][3:Temp_n]-Temp[c2,2][1:(Temp_n-2)]
      Dist2[Dist2==0]<-max(Dist2[Dist2!=0])
      
      #Crowding Ditsance - Retain
      Dist3[2:(Temp_n-1)]<-Temp[c3,3][3:Temp_n]-Temp[c3,3][1:(Temp_n-2)]
      Dist3[Dist3==0]<-max(Dist3[Dist3!=0])
      
      #Total Crowding Distance
      Temp[c1,5]<-Dist1
      Temp[c2,5]<-Temp[c2,5]+Dist2
      Temp[c3,5]<-Temp[c3,5]+Dist3
      
      #Update
      Temp_Evl[Temp_Evl[,4]==i,5]<-Temp[,5]	
    }
  }
  
  #Adjustment
  Temp_Evl[abs(Temp_Evl[,5])==Inf,5]<-0
  Temp_Evl[is.na(Temp_Evl[,5]),5]<-0
  
  #Prepare Output
  out<-list(Temp_Evl, Temp_Pop)
  names(out)<-c("Evl", "Pop")
  return(out)
}

#Reproduction
EO_Reproduction<-function(Evl, Pop, Archive, prob_m, C_Min, C_Max, Seed){
  
  #Reference
  Pop_row<-nrow(Pop)
  Pop_col<-ncol(Pop)
  Arc_col<-ncol(Archive$Pop)
  
  #Data Container
  NP<-Pop
  AR<-Archive
  
  #Reproducibility
  set.seed(Seed)
  
  #Archive Truncate
  AR<-ParetoSort(Pop_Evl = AR)
  AR<-CrowdingDist(Pop_Evl = AR)
  c0<-order(-rank(AR$Evl[,3]),AR$Evl[,4], decreasing = TRUE)
  AR$Evl<-AR$Evl[c0[1:min(Arc_col,Pop_N)],]
  AR$Pop<-AR$Pop[,c0[1:min(Arc_col,Pop_N)],]
  Arc_col<-ncol(AR$Pop)
  
  #Cross Over	
  Parent_A<-sample(x = Arc_col, size = Pop_col/2, replace = TRUE)
  Parent_B<-sample(x = Pop_col, size = Pop_col/2, replace = TRUE)			
  Select_A<-matrix(sample(x = 0:1, size = Pop_row*Pop_col/2, replace = TRUE), nrow = Pop_row, ncol = Pop_col/2, byrow = TRUE)
  Select_B<-ifelse(Select_A==0,1,0)
  NP[,,1]<-AR$Pop[,,1][,Parent_A]*Select_A+NP[,,1][,Parent_B]*Select_B
  NP[,,2]<-AR$Pop[,,2][,Parent_A]*Select_A+NP[,,2][,Parent_B]*Select_B  
  
  #Mutation
  Select_M<-matrix(sample(x = 0:1, size = Pop_row*Pop_col, replace = TRUE, prob = c(1-prob_m, prob_m)), nrow = Pop_row, ncol = Pop_col, byrow = TRUE)
  Delta_M<-matrix(runif(n = Pop_row*Pop_col, min = C_Min, max = C_Max), nrow = Pop_row, ncol = Pop_col, byrow = TRUE)
  NP[,,1]<-ifelse(Select_M==1,ifelse(NP[,,1]==0,1,0),NP[,,1])
  NP[,,2]<-NP[,,2]+Delta_M*NP[,,1]
  NP[,,2]<-ifelse(NP[,,2]<C_Min,C_Min,ifelse(NP[,,2]>C_Max,C_Max,NP[,,2]))
  
  #Prepare Output
  return(NP)
}

#Merge Population
EO_Merge<-function(Pop, Chl, Pop_Evl, Chl_Evl){
  
  #Reference
  Pop_col<-ncol(Pop)
  Pop_row<-nrow(Pop)
  
  #Data Container
  out<-list()
  out$Evl<-Pop_Evl
  Pop_Temp<-Pop
  
  #Dominance Conditions
  c_p1<-Pop_Evl[,1]< Chl_Evl[,1] & Pop_Evl[,2]>=Chl_Evl[,2] & Pop_Evl[,3]>=Chl_Evl[,3]
  c_p2<-Pop_Evl[,1]<=Chl_Evl[,1] & Pop_Evl[,2]> Chl_Evl[,2] & Pop_Evl[,3]>=Chl_Evl[,3]
  c_p3<-Pop_Evl[,1]<=Chl_Evl[,1] & Pop_Evl[,2]>=Chl_Evl[,2] & Pop_Evl[,3]> Chl_Evl[,3]
  c_c1<-Chl_Evl[,1]< Pop_Evl[,1] & Chl_Evl[,2]>=Pop_Evl[,2] & Chl_Evl[,3]>=Pop_Evl[,3]
  c_c2<-Chl_Evl[,1]<=Pop_Evl[,1] & Chl_Evl[,2]> Pop_Evl[,2] & Chl_Evl[,3]>=Pop_Evl[,3]
  c_c3<-Chl_Evl[,1]<=Pop_Evl[,1] & Chl_Evl[,2]>=Pop_Evl[,2] & Chl_Evl[,3]> Pop_Evl[,3]
  
  #Dominance Positions
  p_p<-which(c_p1&c_p2&c_p3)
  p_c<-which(c_c1&c_c2&c_c1)
  p_o<-(1:Pop_col)[!(1:Pop_col)%in%c(p_p,p_c)]
  
  #Substitute Dominant Children
  out$Evl[p_c,]<-Chl_Evl[p_c,]
  Pop_Temp[,p_c,]<-Chl[,p_c,]
  
  #Add Non Dominant Children
  out$Evl<-rbind(out$Evl, Chl_Evl[p_o,])
  out$Pop<-array(0,c(Pop_row, nrow(out$Evl), 2))
  out$Pop[,,1]<-cbind(Pop_Temp[,,1], Chl[,p_o,1])
  out$Pop[,,2]<-cbind(Pop_Temp[,,2], Chl[,p_o,2])
  
  #Preapre Output
  return(out)
}

#Truncate Population
EO_Truncate<-function(Pop_Mrg, Pop_N){
  
  #Reference
  Temp<-Pop_Mrg
  
  #Data Container
  out<-list()
  
  #Pareto Sorting
  Temp<-ParetoSort(Pop_Evl = Temp)
  
  #Crowding Distance
  Temp<-CrowdingDist(Pop_Evl = Temp)
  
  #Truncate Population
  c0<-order(-rank(Temp$Evl[,4]),Temp$Evl[,5], decreasing = TRUE)
  out$Evl<-Temp$Evl[c0[1:Pop_N],]
  out$Pop<-Temp$Pop[,c0[1:Pop_N],]
  
  #Prepare Output	
  out$Evl[,4:5]<-0
  return(out)
}

#Evolutionary Search
SEARCH_GA<-function(m, S, P, alfa, Beta0, Beta1, Beta2, C_Min, C_Max, Generation, Pop_N, prob_m, Seed){
  
  #Data Container
  out<-list()
  Archive<-list()
  
  #Reproducibility
  set.seed(Seed)
  
  #Population Initialization
  Pop<-EO_Init(m = m, Pop_N = Pop_N, C_Min = C_Min, C_Max = C_Max, Seed = Seed)
  
  #Risk/Return/Retain Evaluation
  Pop_Evl<-EO_Eval(m = m, S = S, P = P, Pop = Pop, alfa = alfa, Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2, Seed = Seed)
  
  #Update Archive
  Archive<-ClearMemory(Evl = Pop_Evl, Pop = Pop)
  
  #Evolutionary Optimization
  for(i in 1:Generation){
    
    #Candidate Generation
    Chl<-EO_Reproduction(Evl = Pop_Evl, Pop = Pop, Archive = Archive, prob_m = prob_m, C_Min = C_Min, C_Max = C_Max, Seed = Seed+i)
    
    #Candidate Risk/Return Evaluation
    Chl_Evl<-EO_Eval(m = m, S = S, P = P, Pop = Chl, alfa = alfa, Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2, Seed = Seed+i)
    
    #Merge Population
    Pop_Mrg<-EO_Merge(Pop = Pop, Chl = Chl, Pop_Evl = Pop_Evl, Chl_Evl = Chl_Evl)
    
    #Selection
    Pop_Trc<-EO_Truncate(Pop_Mrg = Pop_Mrg, Pop_N = Pop_N)
    Pop<-Pop_Trc$Pop
    Pop_Evl<-Pop_Trc$Evl
    
    #Update EF Memory
    Archive_Evl<-rbind(Archive$Evl, Pop_Evl)
    Archive_Pop<-array(0,c(m,nrow(Archive_Evl),2))
    Archive_Pop[,,1]<-cbind(Archive$Pop[,,1], Pop[,,1])
    Archive_Pop[,,2]<-cbind(Archive$Pop[,,2], Pop[,,2])
    Archive<-ClearMemory(Evl = Archive_Evl, Pop = Archive_Pop)
    
    #Message
    #message("Generation: ", i)
  }
  
  #Prepare Output
  out$Evl<-Archive$Evl
  out$Pop<-Archive$Pop
  
  return(out)
}

#03.Performance Measures -------------------------------------------------------------------------

Performance<-function(EF_RS, EF_GA){
  
  #Library Upload
  library(emoa)
  
  #Data Container
  out<-as.data.frame(matrix(0, 6, 3))
  names(out)<-c("Metric", "Random", "Evolutionary")
  
  #M1_Hyper Volume
  out[1,1]<-"HyperVolume"
  out[1,2]<-dominated_hypervolume(points = t(EF_RS), ref = c(1,1,1))
  out[1,3]<-dominated_hypervolume(points = t(EF_GA), ref = c(1,1,1))
  
  #M2_Cardinality
  out[2,1]<-"Cardinality"
  out[2,2]<-nrow(EF_RS)
  out[2,3]<-nrow(EF_GA)
  
  #M3_Spacing
  Measure_Spacing<-function(x){
    D<-as.matrix(dist(x))
    D_Min<-apply(D, 2, function(x){sort(x)[2]})
    sum((D_Min-mean(D_Min))^2)/nrow(x)*1000
  }
  
  out[3,1]<-"Spacing"
  out[3,2]<-Measure_Spacing(EF_RS)
  out[3,3]<-Measure_Spacing(EF_GA)
  
  #M4_Spread
  Measure_Spread<-function(x){
    D<-as.matrix(dist(x))		
    D_e<-sum(apply(x, 2, function(x){max(x)-min(x)}))
    D_Min<-apply(D, 2, function(x){sort(x)[2]})
    D_Mean<-mean(D)
    (D_e+sum(abs(D_Min-D_Mean)))/(D_e+nrow(x)*D_Mean)
  }
  
  out[4,1]<-"Spread"
  out[4,2]<-Measure_Spread(EF_RS)
  out[4,3]<-Measure_Spread(EF_GA)
  
  #M5_Range
  out[5,1]<-"Range"
  out[5,2]<-sum(apply(EF_RS, 2, function(x){max(x)-min(x)}))
  out[5,3]<-sum(apply(EF_GA, 2, function(x){max(x)-min(x)}))
  
  #M6_Dominance
  out[6,1]<-"Dominance"
  EF_MR<-CalcFrontier(Evl = rbind(EF_RS,EF_GA))
  EF_MR<-EF_MR[EF_MR[,4]==1,1:3]
  out[6,2]<-sum(apply(EF_RS, 1, function(x){ifelse(x[1]%in%EF_MR[,1]&x[2]%in%EF_MR[,2]&x[3]%in%EF_MR[,3],1,0)}))
  out[6,3]<-sum(apply(EF_GA, 1, function(x){ifelse(x[1]%in%EF_MR[,1]&x[2]%in%EF_MR[,2]&x[3]%in%EF_MR[,3],1,0)}))
  
  #Prepare Output
  return(out)
}

#04.Master ---------------------------------------------------------------------------------------

Run_Single<-function(s, m, F_mn, S_mn, S_sd, F_mn_sd, S_mn_sd, S_sd_sd, alfa, Iteration, Generation, Pop_N, prob_m, Beta0, Beta1, Beta2, C_Min, C_Max, Seed){
  
  #Data Container
  out<-list()
  
  #Generate Synthetic Portfolio
  DB<-SyntheticDB(s = s, m = m, F_mn = F_mn, S_mn = S_mn, S_sd = S_sd, F_mn_sd = F_mn_sd, S_mn_sd = S_mn_sd, S_sd_sd = S_sd_sd, Seed = Seed)
  
  #Pareto Frontier Search
  PTF_RS<-SEARCH_RS(m = m, S = DB$S, P = DB$P, alfa = alfa, Iteration = Iteration,                                   Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2, C_Min = C_Min, C_Max = C_Max, Seed = Seed)
  PTF_GA<-SEARCH_GA(m = m, S = DB$S, P = DB$P, alfa = alfa, Generation = Generation, Pop_N = Pop_N, prob_m = prob_m, Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2, C_Min = C_Min, C_Max = C_Max, Seed = Seed)
  
  #Performance Evaluation 
  PERFORMANCE<-Performance(EF_RS = PTF_RS$Evl[,1:3], EF_GA = PTF_GA$Evl[,1:3])
  
  #Prepare Output
  out$Portfolio<-DB$R
  out$RS_Evl<-PTF_RS$Evl
  out$RS_Pop<-PTF_RS$Pop
  out$GA_Evl<-PTF_GA$Evl
  out$GA_Pop<-PTF_GA$Pop
  out$Performance<-PERFORMANCE
  return(out)
}

Run_Statistic<-function(Run){
  
  #Data Container
  out<-data.frame(matrix(0,1,48))
  names(out)<-c("m",  
                "Freq_Mean", "Freq_Sd", "Freq_Skew", "Freq_Min", "Freq_Max",
                "Sev_Mean_Mean", "Sev_Mean_Sd", "Sev_Mean_Skew", "Sev_Mean_Min", "Sev_Mean_Max",
                "Sev_Sd_Mean", "Sev_Sd_Sd", "Sev_Sd_Skew", "Sev_Sd_Min", "Sev_Sd_Max",
                "UMCS_Hyper", "UMCS_Cardinality", "UMCS_Spacing", "UMCS_Spread", "UMCS_Range", "UMCS_Dominance",
                "UMCS_Selection_Mean", "UMCS_Selection_Sd", "UMCS_Selection_Skew", "UMCS_Selection_Min", "UMCS_Selection_Max",
                "UMCS_Renewal_Mean", "UMCS_Renewal_Sd", "UMCS_Renewal_Skew", "UMCS_Renewal_Min", "UMCS_Renewal_Max",
                "DEMO_Hyper", "DEMO_Cardinality", "DEMO_Spacing", "DEMO_Spread", "DEMO_Range", "DEMO_Dominance",
                "DEMO_Selection_Mean", "DEMO_Selection_Sd", "DEMO_Selection_Skew", "DEMO_Selection_Min", "DEMO_Selection_Max",
                "DEMO_Renewal_Mean", "DEMO_Renewal_Sd", "DEMO_Renewal_Skew", "DEMO_Renewal_Min", "DEMO_Renewal_Max")
  
  
  #Selections Statistics
  UMCS_Selec<-rowSums(Run$RS_Pop[,,1])
  DEMO_Selec<-rowSums(Run$GA_Pop[,,1])
  
  UMCS_Renew<-apply(Run$RS_Pop[,,2],2,function(x){mean(x[x!=0])})
  DEMO_Renew<-apply(Run$GA_Pop[,,2],2,function(x){mean(x[x!=0])})
  
  UMCS_Renew<-ifelse(is.na(UMCS_Renew),0,UMCS_Renew)
  DEMO_Renew<-ifelse(is.na(DEMO_Renew),0,DEMO_Renew)
  
  #Summary Statistic
  out$m<-nrow(Run$Portfolio)
  
  out$Freq_Mean<-mean(Run$Portfolio[,1])
  out$Freq_Sd<-sd(Run$Portfolio[,1])
  out$Freq_Skew<-skewness(Run$Portfolio[,1])
  out$Freq_Min<-min(Run$Portfolio[,1])
  out$Freq_Max<-max(Run$Portfolio[,1])
  
  out$Sev_Mean_Mean<-mean(Run$Portfolio[,2])
  out$Sev_Mean_Sd<-sd(Run$Portfolio[,2])
  out$Sev_Mean_Skew<-skewness(Run$Portfolio[,2])
  out$Sev_Mean_Min<-min(Run$Portfolio[,2])
  out$Sev_Mean_Max<-max(Run$Portfolio[,2])
  
  out$Sev_Sd_Mean<-mean(Run$Portfolio[,3])
  out$Sev_Sd_Sd<-sd(Run$Portfolio[,3])
  out$Sev_Sd_Skew<-skewness(Run$Portfolio[,3])
  out$Sev_Sd_Min<-min(Run$Portfolio[,3])
  out$Sev_Sd_Max<-max(Run$Portfolio[,3])
  
  out$UMCS_Hyper<-Run$Performance$Random[1]
  out$UMCS_Cardinality<-Run$Performance$Random[2]
  out$UMCS_Spacing<-Run$Performance$Random[3]
  out$UMCS_Spread<-Run$Performance$Random[4]
  out$UMCS_Range<-Run$Performance$Random[5]
  out$UMCS_Dominance<-Run$Performance$Random[6]
  
  out$UMCS_Selection_Mean<-mean(UMCS_Selec)
  out$UMCS_Selection_Sd<-sd(UMCS_Selec)
  out$UMCS_Selection_Skew<-skewness(UMCS_Selec)
  out$UMCS_Selection_Min<-min(UMCS_Selec)
  out$UMCS_Selection_Max<-max(UMCS_Selec)
  
  out$UMCS_Renewal_Mean<-mean(UMCS_Renew)
  out$UMCS_Renewal_Sd<-sd(UMCS_Renew)
  out$UMCS_Renewal_Skew<-skewness(UMCS_Renew)
  out$UMCS_Renewal_Min<-min(UMCS_Renew)
  out$UMCS_Renewal_Max<-max(UMCS_Renew)
  
  out$DEMO_Hyper<-Run$Performance$Evolutionary[1]
  out$DEMO_Cardinality<-Run$Performance$Evolutionary[2]
  out$DEMO_Spacing<-Run$Performance$Evolutionary[3]
  out$DEMO_Spread<-Run$Performance$Evolutionary[4]
  out$DEMO_Range<-Run$Performance$Evolutionary[5]
  out$DEMO_Dominance<-Run$Performance$Evolutionary[6]
  
  out$DEMO_Selection_Mean<-mean(DEMO_Selec)
  out$DEMO_Selection_Sd<-sd(DEMO_Selec)
  out$DEMO_Selection_Skew<-skewness(DEMO_Selec)
  out$DEMO_Selection_Min<-min(DEMO_Selec)
  out$DEMO_Selection_Max<-max(DEMO_Selec)
  
  out$DEMO_Renewal_Mean<-mean(DEMO_Renew)
  out$DEMO_Renewal_Sd<-sd(DEMO_Renew)
  out$DEMO_Renewal_Skew<-skewness(DEMO_Renew)
  out$DEMO_Renewal_Min<-min(DEMO_Renew)
  out$DEMO_Renewal_Max<-max(DEMO_Renew)
  
  #Prepare Output
  return(out)
}

Run_All<-function(s, m, F_mn, S_mn, S_sd, F_mn_sd_Max, S_mn_sd_Max, S_sd_sd_Max, alfa, Iteration, Generation, Pop_N, prob_m, Beta0, Beta1, Beta2, C_Min, C_Max, Seed, Macro_Iter, Macro_Cluster){
  
  #Reproducibility
  set.seed(Seed)
  
  #Simulate Individual Variability levels
  DB_Var<-data.frame(F_mn_sd = runif(n = Macro_Iter, min = 0, max = F_mn_sd_Max),
                     S_mn_sd = runif(n = Macro_Iter, min = 0, max = S_mn_sd_Max),
                     S_sd_sd = runif(n = Macro_Iter, min = 0, max = S_sd_sd_Max))
  
  DB_Var<-DB_Var[order(rowSums(DB_Var)),]
  
  #Parallel Run - Set Up
  Cluster<-makeCluster(Macro_Cluster)
  registerDoParallel(Cluster)
  
  Function_Exp<-c("CalcFrontier","ClearMemory","CrowdingDist","EO_Eval","EO_Init","EO_Merge","EO_Reproduction","EO_Truncate","Market","ParetoSort","Performance","Ptf_Eval","Run_Single","SEARCH_GA","SEARCH_RS","SyntheticDB","Run_Single")
  
  #Parallel Run - Execution
  Run<-foreach(i=1:Macro_Iter, .export=Function_Exp)%dopar%{
    
    Run_Single(s = s, m = m, F_mn = F_mn, S_mn = S_mn, S_sd = S_sd, 
               F_mn_sd = DB_Var$F_mn_sd[i], 
               S_mn_sd = DB_Var$S_mn_sd[i], 
               S_sd_sd = DB_Var$S_sd_sd[i], 
               alfa = alfa, 
               Iteration = Iteration, 
               Generation = Generation, 
               Pop_N = Pop_N, prob_m = prob_m, 
               Beta0 = Beta0, Beta1 = Beta1, Beta2 = Beta2, 
               C_Min = C_Min, C_Max = C_Max, 
               Seed = Seed)    
  }
  
  #Parallel Run - Close
  stopCluster(Cluster)
  
  #Prepare Output
  out<-as.data.frame(sapply(Run, Run_Statistic))
  names(out)<-paste("Simulation_",1:ncol(out),sep="")
  
  #Prepare Output
  return(out)
}

#05.Setting --------------------------------------------------------------------------------------

#Macro Setting
alfa<-0.99					#TVaR Level of confidence
C_Min<- -0.50				#Premium Change Floor
C_Max<- +0.50				#Premium Change Ceiling
Macro_Iter<-36				#Number of Algorithms Restart
Macro_Cluster<-3				#Number of Available Core

#Synthetic Data Generator
s<-10000					#Number of simulation
F_mn<-1					#Portfolio Frequency Mean
S_mn<-5					#Portfolio Severity Mean
S_sd<-1					#Portfolio Severity Standard Deviation
F_mn_sd_Max<-2.5				#Max Portfolio Frequency Mean Variability
S_mn_sd_Max<-1.5				#Max Portfolio Severity Mean Variability
S_sd_sd_Max<-0.5				#Max Portfolio Severity Standard Deviation Varaibility

#Random Search
Iteration<-1500				  #Number of Simulated Portfolio			

#Multi Objective GA
Generation<-10				  #Number of Algorithm Iteration
Pop_N<-150					    #Population Dimension
prob_m<-0.01				    #Mutation Frequency

#06.Experiment -----------------------------------------------------------------------------------

#Upload Console
setwd("C:/Users/Andrea/Desktop/PhD/PORTFOLIO OPTIMIZATION")
Console<-as.data.frame(read.csv("RunConsole.csv", sep=";"))

#Console Reference
Console_N<-nrow(Console)

#Experiment
for(i in 1:Console_N){
  
  #Run Simulation
  Console_Temp<-Run_All(Seed = Console$Seed[i], 
                        m = Console$m[i], 
                        Beta0 = Console$Beta0[i] , 
                        Beta1 = Console$Beta1[i], 
                        Beta2 = Console$Beta2[i], 
                        Macro_Iter = Macro_Iter, 
                        Macro_Cluster = Macro_Cluster, 
                        s = s, 
                        F_mn = F_mn, 
                        S_mn = S_mn, 
                        S_sd = S_sd, 
                        F_mn_sd_Max = F_mn_sd_Max, 
                        S_mn_sd_Max = S_mn_sd_Max, 
                        S_sd_sd_Max = S_sd_sd_Max, 
                        alfa = alfa, 
                        Iteration = Iteration, 
                        Generation = Generation, 
                        Pop_N = Pop_N, 
                        prob_m = prob_m, 
                        C_Min = C_Min, 
                        C_Max = C_Max)
  
  #Download Result
  write.csv(cbind(rownames(Console_Temp),
                  apply(Console_Temp,2,as.character)), 
            paste(Console$Name[i],".csv",sep=""), 
            row.names = TRUE)
}


rm(list = ls())

dfweb.surv <- data.table::data.table(haven::read_sav(
        'Y:/3rd_yr/DCE_missing/All_participants_No_exclusion/datasets/dfwebpref.sav'))

# 1st keep those who chose the optout less than once
dfweb.surv<- dfweb.surv[nb_optout<=1, ]

# New PID for the personID after exclusion of those with >1 optout choices
dfweb.surv[, `:=`(PID = match(personID, table = unique(personID)))]

I<- length(unique(dfweb.surv[,PID]))
S<- length(unique(dfweb.surv[,setID]))
M<- length(unique(dfweb.surv[,altID]))

# choice data
cdat.MSI <- array(rep(0,M*S*I), c(M, S, I))
for(i in 1:I){  for(s in 1:S){
        cdat.MSI[,s,i] <- matrix(as.numeric(
                dfweb.surv[PID==i & setID==s, c(select)]), nrow = M)      
} }

# design data
des <- unique(dfweb.surv[,.(setID, altID, optout, trustlabel, delivery_time, delivery_price, 
                       returning, web_usability, discount, headqtrs, distance_shipment100, 
                       I_info_trustlabel, I_info_headqtrs, I_info_distance_shipment)])

P <- ncol(des[,-(1:2)])
des.PMS<- array(rep(0, P*M*S), c(P, M, S))
for(s in 1:S){ for(m in 1:M) {
        des.PMS[,m,s] <- matrix(as.numeric( des[setID==s & altID==m, 3:14]), nrow = P)
}}

# choice sets to be included in the analysis for those that dropped out
fdat.SI <- array(rep(0,S*I), c(S,I))
for(i in 1:I){ fdat.SI[,i] <- matrix(1- as.numeric(unique(
        dfweb.surv[PID==i, .(PID, setID, stop_experiment)])[,stop_experiment]), nrow = S)}

# choice complexity data
xwd <- 'C:/KU_Leuven/DCE/PhD.projects/3rd_yr/DCE_missing/'
xwd <-paste0(xwd,'datasets/dfsim_design_and_metrics.txt')
setComplexity <- data.table::data.table(read.table(xwd))

setComplexity[, `:=`(av_NmP = (NmP12+NmP13+NmP23)/3, 
                     av_NmP_lag1 = (av_NmP12_lag1+av_NmP13_lag1+av_NmP23_lag1)/3,
                     av_NmP_prevs = (av_NmP12_prevs+av_NmP13_prevs+av_NmP23_prevs)/3)]

setComplexity<- setComplexity[,.(setID, aSD_j, dSD_j,av_NmP)]
setComplexity[, `:=`(I= 1)]
Pcomplex <- 4

Xcomp.PS <- array(rep(0, Pcomplex*S), c(Pcomplex,S))

for(s in 1:S){  Xcomp.PS[,s] <- matrix(as.numeric(
        setComplexity[setID==s, .(I, av_NmP, aSD_j, dSD_j)]), nrow = Pcomplex)
} 

Xcomp.PSI <- Xcomp.PS%o%rep(1, I)

drop.SI<-array(rep(0,S*I),c(S,I))
for(i in 1:I){ drop.SI[,i] <- matrix(as.numeric(unique(
        dfweb.surv[PID==i, .(PID, setID, stop_experiment)])[,stop_experiment]), nrow = S)}

# use only relevant pieces in each respondent's array
INmodel.SI <- array(rep(1, I*S), c(S,I))
dfweb.survTemp <- unique(dfweb.surv[,.(PID, setID, stop_experiment)])
dfweb.survTemp[, `:=`(ToDrop = cumsum(stop_experiment)), by = .(PID)]
for (i in 1:I) {
        for(s in 1:S){
        INmodel.SI[s,i] <- ifelse(s %in% dfweb.survTemp[PID==i & ToDrop<=1, setID], 1, 0)
        }}
rm(dfweb.survTemp)
# MNL + SURV function ####
logist<-function(x){1/(1+exp(-x))} # sigmoid function

fdce.surv<-function(param){
        
        #prior
        # prior<-sum(log(dnorm(param)))
        
        beta<-param[1:P]
        
        beta.PMS<-((beta%o%rep(1,M))%o%rep(1,S))
         
        # DCE: utilities and choice probabilities
        
        u.MS<-apply(des.PMS*beta.PMS,c(2,3),sum)
        pi.MSI<-(exp(u.MS)/rep(1,M)%o%apply(exp(u.MS),c(2),sum))%o%rep(1,I)
        
        tempdce.SI<- INmodel.SI*(apply(cdat.MSI*log(pi.MSI),c(2,3),sum))
        
        # dropout behavior
        theta <- param[(P+1):(P+Pcomplex)]
        theta.PSI <- ((theta%o%rep(1,S))%o%rep(1,I))
        
        pdrop.SI <- logist(apply(Xcomp.PSI*theta.PSI, c(2,3), sum))
        
        # dropout probabilities
        tempdrop.SI <- INmodel.SI*(drop.SI*log(pdrop.SI)+(1-drop.SI)*log(1-pdrop.SI))
        
        temp <- exp(apply(tempdce.SI + tempdrop.SI,c(2),sum))
        # f <- -2*(prior + sum(log(temp)))
        f <- -sum(log(temp))
}
beta.s <- c(-3.1, 0.7,-0.07,-0.1,-0.60,0.12,0.01,0.36,-0.03, 0.48, -0.01, -0.04)
theta.s <- c(-0.10, 0.109, 0.11, -1.46)
param.s <- c(beta.s, theta.s)

mdce.surv<-optim(param.s, fdce.surv,hessian=TRUE,method="BFGS")
mdce.surv$par
mdce.surv$value

cbind('Est' = mdce.surv$par, 'Std. Error' = sqrt(diag(solve(mdce.surv$hessian))), 
      'z-value' = mdce.surv$par/ sqrt(diag(solve(mdce.surv$hessian))),
      'Pr(>|z|)' = 2*pnorm(mdce.surv$par/ sqrt(diag(solve(mdce.surv$hessian)))))





############# END ####
# ==== DCE only code === ####
f<-function(param){
  
  #prior
  # prior<-sum(log(dnorm(param)))
  
  beta<-param[1:P]
  
  beta.PMS<-((beta%o%rep(1,M))%o%rep(1,S))
  
  #utilities and choice probabilities
  
  u.MS<-apply(des.PMS*beta.PMS,c(2,3),sum)
  pi.MSI<-(exp(u.MS)/rep(1,M)%o%apply(exp(u.MS),c(2),sum))%o%rep(1,I)
  
  temp.SI<-exp(apply(cdat.MSI*log(pi.MSI),c(2,3),sum))
  # f<--2*(prior+sum(log(temp)))
  f<--2*(sum(apply(fdat.SI*log(temp.SI), c(2), sum)))
}
param.s<-c(beta.s)
m<-optim(param.s,f,hessian=TRUE,method="BFGS")

round(m$par,3) # reproduced successfully
cbind('Est' = m$par, 'Std. Error' = sqrt(diag(solve(m$hessian))), 
      'z-value' = m$par/ sqrt(diag(solve(m$hessian))),
      'Pr(>|z|)' = 2*pnorm(m$par/ sqrt(diag(solve(m$hessian)))))


# surv/dropout ONLY CODE ####
fcomplex <- function(param){
  
  # prior
  # prior <- sum(log(dnorm(param)))
  
  beta <- param[1:Pcomplex]
  beta.PSI <- ((beta%o%rep(1,S))%o%rep(1,I))
  
  # 
  pin.SI <- logist(apply(Xcomp.PSI*beta.PSI, c(2,3), sum))
  
  #probabilities
  temp <- exp(apply(INmodel.SI*(drop.SI*log(pin.SI)+(1-drop.SI)*log(1-pin.SI)),2, sum))
  # f <- -2*(prior + sum(log(temp)))
  f <- -sum(log(temp))
}
beta.complex <- c(-1, 1.0, 5, -2)
mf.complex <-optim(beta.complex, fcomplex, hessian=TRUE, method = "BFGS")
mf.complex$par

# glm for comparison
dfweb.survTemp <- unique(dfweb.surv[,.(PID, setID, stop_experiment)])
dfweb.survTemp[, `:=`(ToDrop = cumsum(stop_experiment)), by = .(PID)]
dfweb.survTemp <- merge(dfweb.survTemp, setComplexity, by = 'setID')
mcomplex <- glm(stop_experiment~ av_NmP + aSD_j + dSD_j, family = 'binomial', 
                data = dfweb.survTemp[ToDrop<=1,])

cbind('Est' = mf.complex$par, 'Std. Error' = sqrt(diag(solve(mf.complex$hessian))), 
      'z-value' = mf.complex$par/ sqrt(diag(solve(mf.complex$hessian))), 
      'Pr(>|z|)' = 2*pnorm(mf.complex$par/ sqrt(diag(solve(mf.complex$hessian)))))
summary(mcomplex)$coefficients
logLik(mcomplex); mf.complex$value


rm(list = ls())

setwd('C:/KU_Leuven/DCE/PhD.projects/3rd_yr/DCE_missing/All_participants_No_exclusion')

# data used in the article
dfweb <- data.table::data.table(haven::read_sav(
        "Y:/2nd_yr/dynamic MNL ANA/project_website_preference/datasets/dfwebpref.sav"))

# MNL code
dfwebmlogit<- mlogit::mlogit.data(dfweb, choice = "select2", shape = 'long', 
                                  alt.var = "altID", id.var = "personID")

#### MNL ####
# mlogit ####
mnl <- gmnl::gmnl(select ~ optout + trustlabel + delivery_time + delivery_price 
                  + returning + web_usability + discount + headqtrs +
                          distance_shipment100 + I_info_trustlabel + 
                          I_info_headqtrs + I_info_distance_shipment| 0, 
                  data = dfwebmlogit, model = 'mnl')
summary(mnl)
BIC(mnl)

# LEO code ####
I<- length(unique(dfweb[,personID]))
S<- length(unique(dfweb[,setID]))
M<- length(unique(dfweb[,altID]))

# choice data
cdat.MSI <- array(rep(0,M*S*I), c(M, S, I))
for(i in 1:I){  for(s in 1:S){
        cdat.MSI[,s,i] <- matrix(as.numeric(
                dfweb[personID==i & setID==s, c(select)]), nrow = M)      
        } }

# design data
des <- unique(dfweb[,.(setID, altID, optout, trustlabel, delivery_time, delivery_price, 
                returning, web_usability, discount, headqtrs, distance_shipment100, 
                I_info_trustlabel, I_info_headqtrs, I_info_distance_shipment)])

P <- ncol(des[,-(1:2)])
des.PMS<- array(rep(0, P*M*S), c(P, M, S))
for(s in 1:S){ for(m in 1:M) {
        des.PMS[,m,s] <- matrix(as.numeric( des[setID==s & altID==m, 3:14]), nrow = P)
}}

f<-function(param){
        
        #prior
        # prior<-sum(log(dnorm(param)))
        
        beta<-param[1:P]
        
        beta.PMS<-((beta%o%rep(1,M))%o%rep(1,S))
        
        #utilities and choice probabilities
        
        u.MS<-apply(des.PMS*beta.PMS,c(2,3),sum)
        pi.MSI<-(exp(u.MS)/rep(1,M)%o%apply(exp(u.MS),c(2),sum))%o%rep(1,I)
        
        temp<-exp(apply(cdat.MSI*log(pi.MSI),c(3),sum))
        # f<--2*(prior+sum(log(temp)))
        f<--2*(sum(log(temp)))
}
beta.s <- c(-3.1, 0.7,-0.07,-0.1,-0.60,0.12,0.01,0.36,-0.03, 0.48, -0.01, -0.04)
param.s<-c(beta.s)
m<-optim(param.s,f,hessian=TRUE,method="BFGS")

cbind(round(m$par,3), round(coef(mnl),3))

# dropout data ####
dfsurv <- data.table::data.table(haven::read_sav("datasets/dfweb_survdata.sav"))

Isurv <- max(dfsurv[,PID]) ; S <- max(dfsurv[,setID]); Psurv <- length(5:16) #6
drop.SI<-array(rep(0,S*Isurv),c(S,Isurv))
for (i in 1:Isurv){
        drop.SI[,i]<- matrix(as.numeric(dfsurv[PID==i , stop_experiment]), nrow = S)       
}

# model design data
X.PSI <- array(rep(0, Psurv*S*Isurv), c(Psurv,S,Isurv))

for (i in 1:Isurv) { 
        for(s in 1:S){ 
                X.PSI[,s,i] <- 
                        matrix(as.numeric(dfsurv[PID==i & setID==s, 5:16]), nrow = Psurv)
        } 
}

# use only relevant pieces in each respondent's array
INmodel.SI <- array(rep(1, Isurv*S), c(S,Isurv))
for (i in 1:Isurv) {
        for(s in 1:S){
                INmodel.SI[s,i] <- ifelse(s %in% dfsurv[PID==i & ToDrop<=1, setID], 1, 0)
        }}

logist<-function(x){1/(1+exp(-x))} # sigmoid function
fsurv <- function(param){
        
        # prior
        # prior <- sum(log(dnorm(param)))
        
        beta <- param[1:Psurv]
        beta.PSI <- ((beta%o%rep(1,S))%o%rep(1,Isurv))
        
        # 
        pin.SI <- logist(apply(X.PSI*beta.PSI, c(2,3), sum))
        
        #probabilities
        temp <- exp(apply(INmodel.SI*(drop.SI*log(pin.SI)+(1-drop.SI)*log(1-pin.SI)),2, sum))
        # f <- -2*(prior + sum(log(temp)))
        f <- -sum(log(temp))
}

tbeta<-c(-0.1, 0.05,0.013,0.04,0.05,0.051,0.081,0.081,0.082, 0.082, 0.083, 0.083)
msurv <- glm(stop_experiment~ as.factor(setID), family = 'binomial', 
             data = dfsurv[ToDrop <=1,]) # 1's after the first event dropped
mf.surv <-optim(tbeta, fsurv, hessian=TRUE, method = "BFGS")
mf.surv$par
cbind('Est' = mf.surv$par, 'Std. Error' = sqrt(diag(solve(mf.surv$hessian))), 
      'z-value' = mf.surv$par/ sqrt(diag(solve(mf.surv$hessian))), 
      'Pr(>|z|)' = 2*pnorm(mf.surv$par/ sqrt(diag(solve(mf.surv$hessian)))))

summary(msurv)$coefficients

cbind(mf.surv$par, coef(msurv))
cbind(sqrt(diag(solve(mf.surv$hessian))),summary(msurv)$coefficients[,2])

mf.surv$value
logLik(msurv)

# drop-out on choice set complexity ####

# design data+choice complexity measures ####
# On how the metrics were calculated and the concerned
# reference article, refer to Design_data_&_metrics.R script

xwd <- 'C:/KU_Leuven/DCE/PhD.projects/3rd_yr/DCE_missing/'
xwd <-paste0(xwd,'datasets/dfsim_design_and_metrics.txt')
setComplexity <- data.table::data.table(read.table(xwd))

setComplexity[, `:=`(av_NmP = (NmP12+NmP13+NmP23)/3, 
                  av_NmP_lag1 = (av_NmP12_lag1+av_NmP13_lag1+av_NmP23_lag1)/3,
                  av_NmP_prevs = (av_NmP12_prevs+av_NmP13_prevs+av_NmP23_prevs)/3)]
dfsurv <- merge(dfsurv, setComplexity, by = "setID")

mcomplex <- glm(stop_experiment~ av_NmP + aSD_j + dSD_j, family = 'binomial', 
                data = dfsurv[ToDrop<=1,])

dfsurv[, `:=`(I= 1)]
Pcomplex <- 4
Xcomp.PSI <- array(rep(0, Pcomplex*S*Isurv), c(Pcomplex,S,Isurv))

for (i in 1:Isurv) { 
        for(s in 1:S){ 
                Xcomp.PSI[,s,i] <- matrix(as.numeric(
                        dfsurv[PID==i & setID==s, .(I, av_NmP, aSD_j, dSD_j)]), 
                        nrow = Pcomplex)
        } 
}

fcomplex <- function(param){
        
        # prior
        # prior <- sum(log(dnorm(param)))
        
        beta <- param[1:Pcomplex]
        beta.PSI <- ((beta%o%rep(1,S))%o%rep(1,Isurv))
        
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
cbind('Est' = mf.complex$par, 'Std. Error' = sqrt(diag(solve(mf.complex$hessian))), 
      'z-value' = mf.complex$par/ sqrt(diag(solve(mf.complex$hessian))), 
      'Pr(>|z|)' = 2*pnorm(mf.complex$par/ sqrt(diag(solve(mf.complex$hessian)))))
summary(mcomplex)$coefficients

logLik(mcomplex); mf.complex$value


# combine MNL and survival codes ####
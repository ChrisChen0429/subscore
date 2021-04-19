library(haven)
library(tidyverse)
library(cocor)
data = read_sav("MSAUSAM3.sav")
data = data %>% filter(IDBOOK==5)

## select the students
data <- data %>% dplyr::select(-names(which(sapply(data, function(x) sum(is.na(x))) == nrow(data))))
data <- data %>% dplyr::select(starts_with("MA"))
names(data)

## score the items
data[is.na(data)] <- 0
data$MA23144 <- ifelse(data$MA23144 == 3, 1, 0)
data$MA23185 <- ifelse(data$MA23185 == 1, 1, 0)
data$MA23054 <- ifelse(data$MA23054 %in% c(20, 21), 1, 0)
data$MA23064 <- ifelse(data$MA23064 == 1, 1, 0)
data$MA23131A <- ifelse(data$MA23131A == 10, 1, 0)
data$MA23131B <- ifelse(data$MA23131B %in% c(10, 11), 1, 0)
data$MA23157 <- ifelse(data$MA23157==20, 1, 0)
data$MA23045 <- ifelse(data$MA23045==3,1,0)
data$MA23082 <- ifelse(data$MA23082==1,1,0)
data$MA23020 <- ifelse(data$MA23020==4,1,0)
data$MA23094 <- ifelse(data$MA23094 %in% c(20,21,22),1,0)
data$MA33027 <- ifelse(data$MA33027==3,1,0)
data$MA33091 <- ifelse(data$MA33091==10,1,0)
data$MA33106 <- ifelse(data$MA33106==2,1,0)
data$MA33090 <- ifelse(data$MA33090==20,1,0)
data$MA33126 <- ifelse(data$MA33126==2,1,0)
data$MA33118A <- ifelse(data$MA33118A==1,1,0)
data$MA33118B <- ifelse(data$MA33118B==3,1,0)
data$MA33118C <- ifelse(data$MA33118C==2,1,0)
data$MA33118 <- ifelse(data$MA33118==20,1,0)
data$MA33243 <- ifelse(data$MA33243%in%c(20,21),1,0)
data$MA33229 <- ifelse(data$MA33229==10,1,0)
data$MA33011 <- ifelse(data$MA33011==1,1,0)
data$MA33159 <- ifelse(data$MA33159==10,1,0)
data$MA33054 <- ifelse(data$MA33054==4,1,0)
data$MA33085 <- ifelse(data$MA33085==4,1,0)
data$MA33190 <- ifelse(data$MA33190==10,1,0)
data$MA33115 <- ifelse(data$MA33115==3,1,0)
data$MA33237 <- ifelse(data$MA33237==10,1,0)
data$MA33077 <- ifelse(data$MA33077==1,1,0)
data$MA33132 <- ifelse(data$MA33132==1,1,0)
data$MA33218 <- ifelse(data$MA33218==20,1,0)
data$MA33236A <- ifelse(data$MA33236A==3,1,0)
data$MA33236B <- ifelse(data$MA33236B==2,1,0)
data$MA33236C <- ifelse(data$MA33236C==1,1,0)
data$MA33236D <- ifelse(data$MA33236D==2,1,0)
data$MA33236 <- ifelse(data$MA33236==20,1,0)
data$MA33181 <- ifelse(data$MA33181==3,1,0)
data$MA33002 <- ifelse(data$MA33002==10,1,0)
data$MA33169 <- ifelse(data$MA33169==2,1,0)
data$MA33235 <- ifelse(data$MA33235==10,1,0)

library(readxl)
ItemInformation <- read_excel("TA15_MAT_ItemInformation.xlsx")
items <- names(data)
domain <- c()
for (i in items){
  domain <- c(domain, ItemInformation %>% filter(`Item ID`==i) %>% dplyr::select(`Content Domain`) %>% pull())
}

alg_index <- which(domain=="Algebra")
cal_index <- which(domain=="Calculus")
geo_index <- which(domain=="Geometry")

final <- cbind(data[,alg_index],data[,cal_index],data[,geo_index])

## subscore 

library(subscore)
test.data<-data.prep(scored.data = final,
                     subtest.infor = c(length(table(domain)),as.numeric(table(domain))),
                     subtest.names = as.vector(names(table(domain))))
ns <- length(table(domain))

## define different types of reliabilities 
# 1. Spearman-Brown
spearman_brown <- function(raw_data){
  total_a <- rep(0,nrow(raw_data))
  total_b <- rep(0,nrow(raw_data))
  for (i in 1:ncol(raw_data)){
    if (i %% 2 == 1){
      total_a <- total_a + raw_data[,i]
    }else{
      total_b <- total_b + raw_data[,i]
    }
  }
  r_ab <- cor(total_a,total_b)
  ans <- (2 * r_ab) / (1 + r_ab)
  return(ans)
}
# 2. Raju 
raju <- function(raw_data){
  total_a <- rep(0,nrow(raw_data))
  total_b <- rep(0,nrow(raw_data))
  lambda1 <- 0
  lambda2 <- 0
  for (i in 1:ncol(raw_data)){
    if (i %% 2 == 1){
      total_a <- total_a + raw_data[,i]
      lambda1 <- lambda1 + 1
    }else{
      total_b <- total_b + raw_data[,i]
      lambda2 <- lambda2 + 1
    }
  }
  cov_ab <- cov(total_a,total_b) ### covariance not correlation
  v_x <- var(rowSums(raw_data)) ### variance of the total score
  lambda1 <- lambda1 / ncol(raw_data)
  lambda2 <- lambda2/ ncol(raw_data)
  ans <- cov_ab / (lambda1 * lambda2 * v_x)
  return(ans)
}
# 3. Angoff
angoff <- function(raw_data){
  total_a <- rep(0,nrow(raw_data))
  total_b <- rep(0,nrow(raw_data))
  for (i in 1:ncol(raw_data)){
    if (i %% 2 == 1){
      total_a <- total_a + raw_data[,i]
    }else{
      total_b <- total_b + raw_data[,i]
    }
  }
  cov_ab <- cov(total_a,total_b)
  v_a <- var(total_a)
  v_b <- var(total_b)
  v_x <- var(rowSums(raw_data))
  ans <- 4 * cov_ab / (v_x - ((v_a - v_b)/sqrt(v_x))^2)
  return(ans)
}
# 4. KR21
kr21 <- function(raw_data){
  s_x <- var(rowSums(raw_data))
  x_bar <- mean(rowSums(raw_data))
  n <- ncol(raw_data)
  ans <- ((n/(n - 1))*(1 - x_bar*(n - x_bar)/(n*s_x)))
  return(ans)
}


# redefine the CTTsub function
CTTsub<-function (test.data, reliability_type = "Alpha") {
  n.tests<-length(test.data)
  n.subtests<-n.tests-1
  n.items<-rep(NA,n.tests)
  n.cases<-rep(NA,n.tests)
  
  for (t in 1:n.tests) {
    n.items[t]<-dim(test.data[[t]])[2] 
    n.cases[t]<-dim(test.data[[t]])[1] 
  } 
  n.items.total<-n.items[n.tests]
  reliability.alpha<-rep(NA, (n.tests))  
  
  subscore.list <- as.list(rep(NA, n.tests))
  names(subscore.list) <- names(test.data)
  for (t in 1 : (n.tests))  {
    subscore.list[[t]]<- rowSums(test.data[[t]],na.rm = T)
  }  
  
  subscore.original.matrix<-do.call(cbind, subscore.list) 
  corr<-cor(subscore.original.matrix)
  itemstrata<-cbind(colnames(test.data[[n.tests]]),rep(1:(n.tests-1),lengths(test.data)[1:(n.tests-1)]))
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  str.alpha<-quiet(stratified.cronbach.alpha(test.data[[n.tests]], 
                                             itemstrata=itemstrata)$alpha.stratified)
  stratefied.alpha<-c(str.alpha[-1],str.alpha[1])
  
  
  if (reliability_type == "Alpha"){
    for (r in 1:(n.tests)) {
      reliability.alpha[r]<-itemAnalysis(test.data[[r]],NA.Delete=T, itemReport=F)$alpha
    }   
  }
  if (reliability_type == "spearman_brown"){
    for (r in 1:(n.tests)) {
      reliability.alpha[r]<-spearman_brown(test.data[[r]])
    }   
  }
  if (reliability_type == "raju"){
    for (r in 1:(n.tests)) {
      reliability.alpha[r]<-raju(test.data[[r]])
    }   
  }
  if (reliability_type == "angoff"){
    for (r in 1:(n.tests)) {
      reliability.alpha[r]<-angoff(test.data[[r]])
    }   
  }
  if (reliability_type == "kr21"){
    for (r in 1:(n.tests)) {
      reliability.alpha[r]<-kr21(test.data[[r]])
    }   
  }
  
  Reliabilities<-cbind(reliability.alpha, stratefied.alpha)
  disattenuated.corr<-disattenuated.cor(corr, reliability.alpha)[-n.tests,-n.tests]
  
  sigma.obs<-rep(NA,n.tests)
  for (t in 1:n.tests) {
    sigma.obs[t]<-sd(subscore.list[[t]],na.rm = TRUE)
  }
  
  var.obs<-sigma.obs^2
  CovMat.Obs<-cov(subscore.original.matrix)
  var.true<-var.obs*reliability.alpha
  sigma.true<-sqrt(var.true)
  CovMat.true<-CovMat.Obs
  for (t in 1:n.tests) {
    CovMat.true[t,t]<-var.true[t]
  }
  
  mean<-rep(NA,n.tests)
  SD<-rep(NA,n.tests)
  for (t in 1:n.tests) {
    mean[t]<-mean(subscore.list[[t]],na.rm = TRUE)
    SD[t]<-sd(subscore.list[[t]],na.rm = TRUE)
  }
  mylist.names <- c(paste('Subscore.s.',rep(names(test.data)[-length(test.data)]),sep=''))
  subscore.list.RegOnSub <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnSub) <- mylist.names
  subscore.dataframe<-as.data.frame(subscore.original.matrix)
  for (t in 1: n.subtests) {
    subscore.list.RegOnSub[[t]]<-mean[t]+reliability.alpha[t]*(subscore.dataframe[,t]-mean[t])
  } 
  PRMSE.s<-rep(NA,n.tests)
  PRMSE.s[1:n.subtests]<-reliability.alpha[1:n.subtests]
  
  PRMSE.x<-rep(NA,n.tests)
  r.StXt<-rep(NA,n.tests)
  
  cov.rowsum<-rowSums(CovMat.true[,1:n.subtests],na.rm = TRUE)
  
  for (t in 1:n.subtests) {    
    r.StXt[t]<-cov.rowsum[t]^2/(var.true[t]*var.true[n.tests])
    PRMSE.x[t]<-r.StXt[t]*reliability.alpha[n.tests]
  } 
  mylist.names <- c(paste('Subscore.x.',rep(names(test.data)[-length(test.data)]),sep=''))
  subscore.list.RegOnTot <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnTot) <- mylist.names
  
  for (t in 1:n.subtests) { 
    subscore.list.RegOnTot[[t]]<-mean[t]+sqrt(PRMSE.x[t])*(sigma.true[t]/(sigma.obs[n.tests])*(subscore.dataframe[,n.tests]-mean[n.tests]))
  } 
  
  tao<-rep(NA,n.tests)
  beta<-rep(NA,n.tests)
  gamma<-rep(NA,n.tests)
  
  for (t in 1:n.tests) {
    tao[t]<-(sqrt(reliability.alpha[n.tests])*sqrt(r.StXt[t])-corr[t,n.tests]*sqrt(reliability.alpha[t]))/(1-corr[t,n.tests]^2)
    beta[t]<- sqrt(reliability.alpha[t])*(sqrt(reliability.alpha[t])-corr[t,n.tests]*tao[t])
    gamma[t]<-sqrt(reliability.alpha[t])*tao[t]*(sigma.obs[t]/sigma.obs[n.tests])
  } 
  
  PRMSE.sx<-rep(NA, n.tests)
  for (t in 1:n.subtests) { 
    PRMSE.sx[t]<-reliability.alpha[t]+tao[t]^2*(1-corr[t,n.tests]^2)
  } 
  
  mylist.names <- c(paste('Subscore.sx.',rep(names(test.data)[-length(test.data)]),sep=''))
  subscore.list.RegOnTotSub <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnTotSub) <- mylist.names
  
  for (t in 1: n.subtests) { 
    subscore.list.RegOnTotSub[[t]]<-mean[t]+beta[t]*(subscore.dataframe[,t]-mean[t])+gamma[t]*(subscore.dataframe[,n.tests]-mean[n.tests])
  } 
  
  added.value.s<-PRMSE.s>PRMSE.x
  added.value.sx<-(PRMSE.sx-pmax(PRMSE.s,PRMSE.x))>(.1*(1-pmax(PRMSE.s,PRMSE.x)))
  
  subscore.information.list<-list(Alpha=reliability.alpha, 
                                  PRMSE.s=PRMSE.s, PRMSE.x=PRMSE.x, PRMSE.sx=PRMSE.sx,
                                  added.value.s=added.value.s,added.value.sx=added.value.sx)
  subscore.information<-do.call(cbind,subscore.information.list)
  rownames(subscore.information)<-names(test.data)
  
  subscore.original<-do.call(cbind,subscore.list)
  subscore.s<-do.call(cbind,subscore.list.RegOnSub)
  subscore.x<-do.call(cbind,subscore.list.RegOnTot)
  subscore.sx<-do.call(cbind,subscore.list.RegOnTotSub)
  
  Orig.mean<-mean[1:n.subtests]
  Orig.sd<-SD[1:n.subtests]
  subscore.s.mean<-colMeans(subscore.s,na.rm=T)
  subscore.x.mean<-colMeans(subscore.x,na.rm=T)
  subscore.sx.mean<-colMeans(subscore.sx,na.rm=T)
  subscore.s.sd<-apply(subscore.s, 2, sd,na.rm=T)
  subscore.x.sd<-apply(subscore.x, 2, sd,na.rm=T)
  subscore.sx.sd<-apply(subscore.sx, 2, sd,na.rm=T)
  
  summary.list<-list(Orig.mean=Orig.mean, Orig.sd=Orig.sd,subscore.s.mean=subscore.s.mean,subscore.s.sd=subscore.s.sd,
                     subscore.x.mean=subscore.x.mean,subscore.x.sd=subscore.x.sd,subscore.sx.mean=subscore.sx.mean,
                     subscore.sx.sd=subscore.sx.sd)
  summary<-do.call(cbind,summary.list)
  rownames(summary)<-names(test.data)[1:n.subtests]
  
  #below code from Sinharay (2019)
  # Compute PRMSEs suggested by Sinharay (2013) from Haberman's PRMSEs
  Olkin.Z<-rep(NA,n.tests)
  Williams.t<-rep(NA,n.tests)
  Hedges.Olkin.Z<-rep(NA,n.tests)
  
  PRs<-PRMSE.s[1:n.subtests]*PRMSE.s[1:n.subtests]
  PRx<-PRMSE.x[1:n.subtests]*PRMSE.s[1:n.subtests]
  PRsx<-PRMSE.sx[1:n.subtests]*PRMSE.s[1:n.subtests]
  rsx<-cor(subscore.original)[1:n.subtests,(n.subtests+1)]
  olk<-rep(0,n.subtests)
  wil<-rep(0,n.subtests)
  n<-n.cases[length(test.data)]
  for (j in 1:n.subtests)
  {olk[j]=cocor.dep.groups.overlap(sqrt(PRs[j]),sqrt(PRx[j]),
                                   rsx[j],n)@olkin1967$statistic
  wil[j]=cocor.dep.groups.overlap(sqrt(PRs[j]),
                                  sqrt(PRx[j]),rsx[j],n)@williams1959$statistic}
  
  compsd<-function(r01,r02,r12,r012,ns,n) {
    a2=2*(r02-r12*r01)/(1-r12*r12)
    a1=-r12*a2
    a3=2*(r12*r01*r01+r12*r02*r02-r01*r02*(1+r12^2))/((1-r12^2)**2)
    V=matrix(0,3,3)
    v11s=(1-r01^2)^2/n
    v22s=(1-r02^2)^2/n
    v33s=(1-r12^2)^2/n
    z=rep(0,4)
    for (j in 1:ns) {
      V[1,2]=(0.5*(2*r12[j]-r01[j]*r02[j])*(1-r12[j]^2-r01[j]^2-r02[j]^2)+r12[j]^3)/n
      V[2,3]=(0.5*(2*r01[j]-r12[j]*r02[j])*(1-r12[j]^2-r01[j]^2-r02[j]^2)+r01[j]^3)/n
      V[1,3]=(0.5*(2*r02[j]-r12[j]*r01[j])*(1-r12[j]^2-r01[j]^2-r02[j]^2)+r02[j]^3)/n
      V[1,1]=v11s[j]
      V[2,2]=v22s[j]
      V[3,3]=v33s[j]
      for (i in 1:2) {for (k in (i+1):3){V[k,i]=V[i,k]}}
      vec=c(a1[j],a2[j],a3[j])
      SD=sqrt(t(vec)%*%V%*%vec)
      z[j]=(r012[j]*r012[j]-r01[j]^2)/SD}
    return(z)}
  
  zs<-compsd(sqrt(PRs),sqrt(PRx),rsx,sqrt(PRsx),n.subtests,n)
  zx<-compsd(sqrt(PRx),sqrt(PRs),rsx,sqrt(PRsx),n.subtests,n)
  hedgesolkin<-ifelse(PRs>PRx,zs,zx)
  
  Olkin.Z[1:n.subtests]<-olk
  Williams.t[1:n.subtests]<-wil
  Hedges.Olkin.Z[1:n.subtests]<-hedgesolkin
  
  subscore.information.list<-list(Alpha=reliability.alpha, 
                                  PRMSE.s=PRMSE.s, PRMSE.x=PRMSE.x, PRMSE.sx=PRMSE.sx,
                                  added.value.s=as.numeric(added.value.s),added.value.sx=as.numeric(added.value.sx),
                                  Olkin.Z=Olkin.Z, Williams.t=Williams.t, 
                                  Hedges.Olkin.Z=Hedges.Olkin.Z)
  subscore.information<-do.call(cbind,subscore.information.list)
  rownames(subscore.information)<-names(test.data)
  
  if (sum(PRMSE.s[PRMSE.s>1],na.rm=T)>=1 | sum(PRMSE.x[PRMSE.x>1],na.rm=T)>=1 | sum(PRMSE.sx[PRMSE.sx>1],na.rm=T)>1) {
    warning ("PRMSE value(s) exceeds 1. The corresponding (augmented) subscore does not have added value.",
             call. = FALSE)
  }
  
  return (list(summary=summary,
               Correlation=corr,
               Disattenuated.correlation=disattenuated.corr, 
               PRMSE=subscore.information[,1:6], 
               PRMSE.test=subscore.information,
               subscore.original=subscore.original,
               subscore.s=subscore.s,
               subscore.x=subscore.x,
               subscore.sx=subscore.sx)) 
}  

Subscores_Alpha = CTTsub(test.data = test.data,reliability_type = 'Alpha')
Subscores_Alpha$PRMSE.test

Subscores_SB = CTTsub(test.data = test.data,reliability_type = "spearman_brown")
Subscores_SB$PRMSE.test

Subscores_raju = CTTsub(test.data = test.data,reliability_type = "raju")
Subscores_raju$PRMSE.test

Subscores_angoff = CTTsub(test.data = test.data,reliability_type = "angoff")
Subscores_angoff$PRMSE.test

Subscores_kr21 = CTTsub(test.data = test.data,reliability_type = "kr21")
Subscores_kr21$PRMSE.test






# Alpha results
dat_Alpha <- data.frame(dens = c(Subscores_Alpha$subscore.original, Subscores_Alpha$subscore.s, Subscores_Alpha$subscore.x, Subscores_Alpha$subscore.sx)
                        , lines = c(rep("subscore.original", each = 2004), rep("subscore.s", each = 1503), rep("subscore.x", each = 1503), rep("subscore.sx", each = 1503)))
# Plot.
plot_1 <- ggplot(dat_Alpha, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + ggtitle("Subscore Augmentation Distribution with Coefficient Alpha Reliability Index")+ 
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))


# Spearman_Brown results
dat_SB <- data.frame(dens = c(Subscores_SB$subscore.original, Subscores_SB$subscore.s, Subscores_SB$subscore.x, Subscores_SB$subscore.sx)
                     , lines = c(rep("subscore.original", each = 2004), rep("subscore.s", each = 1503), rep("subscore.x", each = 1503), rep("subscore.sx", each = 1503)))
# Plot.
plot_2 <- ggplot(dat_SB, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + ggtitle("Subscore Augmentation Distribution with Spearman_Brown Reliability Index")+ 
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))


# Raju results
dat_raju<- data.frame(dens = c(Subscores_raju$subscore.original, Subscores_raju$subscore.s, Subscores_raju$subscore.x, Subscores_raju$subscore.sx)
                      , lines = c(rep("subscore.original", each = 2004), rep("subscore.s", each = 1503), rep("subscore.x", each = 1503), rep("subscore.sx", each = 1503)))
# Plot.
plot_3 <- ggplot(dat_raju, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + ggtitle("Subscore Augmentation Distribution with Raju Reliability Index") + 
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))

# Angoff results
dat_angoff <- data.frame(dens = c(Subscores_angoff$subscore.original, Subscores_angoff$subscore.s, Subscores_angoff$subscore.x, Subscores_angoff$subscore.sx)
                         , lines = c(rep("subscore.original", each = 2004), rep("subscore.s", each = 1503), rep("subscore.x", each = 1503), rep("subscore.sx", each = 1503)))
  
# Plot.
plot_4 <- ggplot(dat_angoff, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + ggtitle("Subscore Augmentation Distribution with Angoff Reliability Index")+ 
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))

# KR21 results
dat_kr21 <- data.frame(dens = c(Subscores_kr21$subscore.original, Subscores_kr21$subscore.s, Subscores_kr21$subscore.x, Subscores_kr21$subscore.sx)
                       , lines = c(rep("subscore.original", each = 2004), rep("subscore.s", each = 1503), rep("subscore.x", each = 1503), rep("subscore.sx", each = 1503)))
# Plot.
plot_5 <- ggplot(dat_kr21, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + ggtitle("Subscore Augmentation Distribution with KR21 Reliability Index")+ 
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))

library(ggpubr)
ggarrange(plot_1, plot_2, plot_3, plot_4,plot_5, ncol=2, nrow=3, common.legend = TRUE, legend="bottom",label.x = seq(0,3,0.5),align = "hv")











data <- data.frame(observed = Subscores_Alpha$subscore.original[,"Algebra"],
                   alpha = Subscores_Alpha$subscore.s[,1],
                   spearman_brown = Subscores_SB$subscore.s[,1],
                   kr21 = Subscores_kr21$subscore.s[,1],
                   raju = Subscores_raju$subscore.s[,1],
                   angoff = Subscores_angoff$subscore.s[,1]
)

p1 <- ggplot(data) + ggtitle("Observed v.s. Predicted Subscore (s) for Algebra")+ 
  geom_point(aes(x=observed,y=alpha,col="alpha")) + geom_line(aes(x=observed,y=alpha,col="alpha")) + 
  geom_point(aes(x=observed,y=spearman_brown,col="spearman_brown")) + geom_line(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  geom_point(aes(x=observed,y=kr21,col="kr21")) + geom_line(aes(x=observed,y=kr21,col="kr21")) + 
  geom_point(aes(x=observed,y=raju,col="raju")) + geom_line(aes(x=observed,y=raju,col="raju")) + 
  geom_point(aes(x=observed,y=angoff,col="angoff")) + geom_line(aes(x=observed,y=angoff,col="angoff")) + 
  ylab("Predicted Subscore") +
  xlab("Observed Subscore") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))



data <- data.frame(observed = Subscores_Alpha$subscore.original[,"Calculus"],
                   alpha = Subscores_Alpha$subscore.s[,2],
                   spearman_brown = Subscores_SB$subscore.s[,2],
                   kr21 = Subscores_kr21$subscore.s[,2],
                   raju = Subscores_raju$subscore.s[,2],
                   angoff = Subscores_angoff$subscore.s[,2]
)

p2 <- ggplot(data) + ggtitle("Observed v.s. Predicted Subscore (s) for Calculus")+ 
  geom_point(aes(x=observed,y=alpha,col="alpha")) + geom_line(aes(x=observed,y=alpha,col="alpha")) + 
  geom_point(aes(x=observed,y=spearman_brown,col="spearman_brown")) + geom_line(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  geom_point(aes(x=observed,y=kr21,col="kr21")) + geom_line(aes(x=observed,y=kr21,col="kr21")) + 
  geom_point(aes(x=observed,y=raju,col="raju")) + geom_line(aes(x=observed,y=raju,col="raju")) + 
  geom_point(aes(x=observed,y=angoff,col="angoff")) + geom_line(aes(x=observed,y=angoff,col="angoff")) + 
  ylab("Predicted Subscore") +
  xlab("Observed Subscore") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))



data <- data.frame(observed = Subscores_Alpha$subscore.original[,"Geometry"],
                   alpha = Subscores_Alpha$subscore.s[,3],
                   spearman_brown = Subscores_SB$subscore.s[,3],
                   kr21 = Subscores_kr21$subscore.s[,3],
                   raju = Subscores_raju$subscore.s[,3],
                   angoff = Subscores_angoff$subscore.s[,3]
)

p3 <- ggplot(data) + ggtitle("Observed v.s. Predicted Subscore (s) for Geometry")+ 
  geom_point(aes(x=observed,y=alpha,col="alpha")) + geom_line(aes(x=observed,y=alpha,col="alpha")) + 
  geom_point(aes(x=observed,y=spearman_brown,col="spearman_brown")) + geom_line(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  geom_point(aes(x=observed,y=kr21,col="kr21")) + geom_line(aes(x=observed,y=kr21,col="kr21")) + 
  geom_point(aes(x=observed,y=raju,col="raju")) + geom_line(aes(x=observed,y=raju,col="raju")) + 
  geom_point(aes(x=observed,y=angoff,col="angoff")) + geom_line(aes(x=observed,y=angoff,col="angoff")) + 
  ylab("Predicted Subscore") +
  xlab("Observed Subscore") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))

library(ggpubr)
ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="bottom",label.x = seq(0,3,0.5),align = "hv")










data <- data.frame(observed = Subscores_Alpha$subscore.original[,"Algebra"],
                   alpha = Subscores_Alpha$subscore.x[,1],
                   spearman_brown = Subscores_SB$subscore.x[,1],
                   kr21 = Subscores_kr21$subscore.x[,1],
                   raju = Subscores_raju$subscore.x[,1],
                   angoff = Subscores_angoff$subscore.x[,1]
)

p1 <- ggplot(data) + ggtitle("Observed v.s. Predicted Subscore (x) for Algebra")+ 
  geom_point(aes(x=observed,y=alpha,col="alpha")) + 
  #geom_line(aes(x=observed,y=alpha,col="alpha")) + 
  geom_point(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  #geom_line(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  geom_point(aes(x=observed,y=kr21,col="kr21")) + 
  #geom_line(aes(x=observed,y=kr21,col="kr21")) + 
  geom_point(aes(x=observed,y=raju,col="raju")) + 
  #geom_line(aes(x=observed,y=raju,col="raju")) + 
  geom_point(aes(x=observed,y=angoff,col="angoff")) + 
  #geom_line(aes(x=observed,y=angoff,col="angoff")) + 
  ylab("Predicted Subscore") +
  xlab("Observed Subscore") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))




data <- data.frame(observed = Subscores_Alpha$subscore.original[,"Calculus"],
                   alpha = Subscores_Alpha$subscore.x[,2],
                   spearman_brown = Subscores_SB$subscore.x[,2],
                   kr21 = Subscores_kr21$subscore.x[,2],
                   raju = Subscores_raju$subscore.x[,2],
                   angoff = Subscores_angoff$subscore.x[,2]
)

p2 <- ggplot(data) + ggtitle("Observed v.s. Predicted Subscore (x) for Calculus")+ 
  geom_point(aes(x=observed,y=alpha,col="alpha")) + 
  #geom_line(aes(x=observed,y=alpha,col="alpha")) + 
  geom_point(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  #geom_line(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  geom_point(aes(x=observed,y=kr21,col="kr21")) + 
  #geom_line(aes(x=observed,y=kr21,col="kr21")) + 
  geom_point(aes(x=observed,y=raju,col="raju")) + 
  #geom_line(aes(x=observed,y=raju,col="raju")) + 
  geom_point(aes(x=observed,y=angoff,col="angoff")) + 
  #geom_line(aes(x=observed,y=angoff,col="angoff")) + 
  ylab("Predicted Subscore") +
  xlab("Observed Subscore") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))



data <- data.frame(observed = Subscores_Alpha$subscore.original[,"Geometry"],
                   alpha = Subscores_Alpha$subscore.x[,3],
                   spearman_brown = Subscores_SB$subscore.x[,3],
                   kr21 = Subscores_kr21$subscore.x[,3],
                   raju = Subscores_raju$subscore.x[,3],
                   angoff = Subscores_angoff$subscore.x[,3]
)

p3 <- ggplot(data) + ggtitle("Observed v.s. Predicted Subscore (x) for Geometry")+ 
  geom_point(aes(x=observed,y=alpha,col="alpha")) + 
  #geom_line(aes(x=observed,y=alpha,col="alpha")) + 
  geom_point(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  #geom_line(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  geom_point(aes(x=observed,y=kr21,col="kr21")) + 
  #geom_line(aes(x=observed,y=kr21,col="kr21")) + 
  geom_point(aes(x=observed,y=raju,col="raju")) + 
  #geom_line(aes(x=observed,y=raju,col="raju")) + 
  geom_point(aes(x=observed,y=angoff,col="angoff")) + 
  #geom_line(aes(x=observed,y=angoff,col="angoff")) + 
  ylab("Predicted Subscore") +
  xlab("Observed Subscore") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))

ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="bottom",label.x = seq(0,3,0.5),align = "hv")











data <- data.frame(observed = Subscores_Alpha$subscore.original[,"Algebra"],
                   alpha = Subscores_Alpha$subscore.sx[,1],
                   spearman_brown = Subscores_SB$subscore.sx[,1],
                   kr21 = Subscores_kr21$subscore.sx[,1],
                   raju = Subscores_raju$subscore.sx[,1],
                   angoff = Subscores_angoff$subscore.sx[,1]
)

p1 <- ggplot(data) + ggtitle("Observed v.s. Predicted Subscore (sx) for Algebra")+ 
  geom_point(aes(x=observed,y=alpha,col="alpha")) + 
  #geom_line(aes(x=observed,y=alpha,col="alpha")) + 
  geom_point(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  #geom_line(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  geom_point(aes(x=observed,y=kr21,col="kr21")) + 
  #geom_line(aes(x=observed,y=kr21,col="kr21")) + 
  geom_point(aes(x=observed,y=raju,col="raju")) + 
  #geom_line(aes(x=observed,y=raju,col="raju")) + 
  geom_point(aes(x=observed,y=angoff,col="angoff")) + 
  #geom_line(aes(x=observed,y=angoff,col="angoff")) + 
  ylab("Predicted Subscore") +
  xlab("Observed Subscore") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))




data <- data.frame(observed = Subscores_Alpha$subscore.original[,"Calculus"],
                   alpha = Subscores_Alpha$subscore.sx[,2],
                   spearman_brown = Subscores_SB$subscore.sx[,2],
                   kr21 = Subscores_kr21$subscore.sx[,2],
                   raju = Subscores_raju$subscore.sx[,2],
                   angoff = Subscores_angoff$subscore.sx[,2]
)

p2 <- ggplot(data) + ggtitle("Observed v.s. Predicted Subscore (sx) for Calculus")+ 
  geom_point(aes(x=observed,y=alpha,col="alpha")) + 
  #geom_line(aes(x=observed,y=alpha,col="alpha")) + 
  geom_point(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  #geom_line(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  geom_point(aes(x=observed,y=kr21,col="kr21")) + 
  #geom_line(aes(x=observed,y=kr21,col="kr21")) + 
  geom_point(aes(x=observed,y=raju,col="raju")) + 
  #geom_line(aes(x=observed,y=raju,col="raju")) + 
  geom_point(aes(x=observed,y=angoff,col="angoff")) + 
  #geom_line(aes(x=observed,y=angoff,col="angoff")) + 
  ylab("Predicted Subscore") +
  xlab("Observed Subscore") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))



data <- data.frame(observed = Subscores_Alpha$subscore.original[,"Geometry"],
                   alpha = Subscores_Alpha$subscore.sx[,3],
                   spearman_brown = Subscores_SB$subscore.sx[,3],
                   kr21 = Subscores_kr21$subscore.sx[,3],
                   raju = Subscores_raju$subscore.sx[,3],
                   angoff = Subscores_angoff$subscore.sx[,3]
)

p3 <- ggplot(data) + ggtitle("Observed v.s. Predicted Subscore (sx) for Geometry")+ 
  geom_point(aes(x=observed,y=alpha,col="alpha")) + 
  #geom_line(aes(x=observed,y=alpha,col="alpha")) + 
  geom_point(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  #geom_line(aes(x=observed,y=spearman_brown,col="spearman_brown")) + 
  geom_point(aes(x=observed,y=kr21,col="kr21")) + 
  #geom_line(aes(x=observed,y=kr21,col="kr21")) + 
  geom_point(aes(x=observed,y=raju,col="raju")) + 
  #geom_line(aes(x=observed,y=raju,col="raju")) + 
  geom_point(aes(x=observed,y=angoff,col="angoff")) + 
  #geom_line(aes(x=observed,y=angoff,col="angoff")) + 
  ylab("Predicted Subscore") +
  xlab("Observed Subscore") +
  theme_minimal()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text =element_text(size=20),
        legend.title =element_text(size=20),
        legend.key.size = unit(1, "cm"))

ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="bottom",label.x = seq(0,3,0.5),align = "hv")





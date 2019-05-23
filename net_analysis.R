setwd("/github/gillnet_analysis-RM-KO")

install.packages("TropFishR")
library(TropFishR)


gillnet_k<-read.csv(file.choose(),stringsAsFactors = T, header = T)
colnames(gillnet_k)<- c("0.5",  "2.0",  "3.0",  "4.0",  "6.0",  "8.0", "10.0")

#unlist the gillnet dataset into a vector - to read it in as matrix
CatchPerNet_mat<- unlist(gillnet_k, use.names = FALSE)
CatchPerNet_mat

midLengths<-c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205, 215, 225, 235, 245, 255, 265, 275, 285)
midLengths
meshSizes<-c(0.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0)
meshSizes

#i developed the CatchPerNet_mat manually but im sure there's a code to auotomate this... Will still search
CatchPerNet_mat = matrix(
  c(34, 143, 78, 16, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 2, 150, 176, 36, 4, 4, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 37, 121, 17, 9, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 13, 21, 18, 29, 22, 7, 10, 17, 8, 2, 6, 10, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 15, 57, 25, 15, 27, 31, 10, 3, 3, 2, 1, 1, 1, 0, 0, 0,
    0, 4, 6, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3, 5, 4, 3, 2, 3, 4, 3, 1, 1, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 6, 12, 5, 8, 2, 4, 4,
    3, 2, 3, 1, 1, 0, 1, 0, 0, 2, 0, 1, 1, 0),
  nrow = 29, #number of rows
  ncol = 7, #number of columns
  byrow = TRUE) #fill matrix by rows

millar = list(midLengths = midLengths, meshSizes = meshSizes, CatchPerNet_mat = CatchPerNet_mat)

# millar_Select Output section 2 ---------------------------------------------------------
output<- select_Millar(millar, x0 = NULL, rtype = "lognorm") #you can change the graph types between normal curve and lognormal
ncolor<- length(output$meshSizes)
plot(output, plotlens=seq(40,90,0.1), deviance_plot = FALSE,
     lty=1, col=rainbow(ncolor))
legend("topright", col = rainbow(ncolor), legend = output$meshSizes,
       lty = 1, title = "Mesh_size [cm]")


# gillnet fit section -----------------------------------------------------
#trial


gillnetfit <- function(data, meshsizes,
                       rtype="norm.loc",
                       rel.power=NULL,
                       plotlens=NULL,
                       details=FALSE){if(sum(sort(meshsizes)==meshsizes)!=length(meshsizes))
                         stop("Mesh sizes must be ascending order")
  lens=rep(data[,1],ncol(data[,-1]))
  msizes=rep(meshsizes,rep(nrow(data),ncol(data[,-1])))
  msize1=msizes[1]
  dat=as.vector(data[,-1])
  var1=lens*msizes; var2=msizes^2; var3=(lens/msizes)^2
  var4=lens/msizes; var5=-log(msizes); var6=log(msizes/msizes[1])
  var7=var6*log(lens) - 0.5*var6*var6; var8=lens*lens
  var9=msizes/lens
  
  if(is.null(plotlens)) plotlens=data[,1]
  if(is.null(rel.power)) os=0
  else os=rep(log(rel.power),rep(nrow(data),ncol(data[,-1])))
  switch(rtype,
         "norm.loc"={
           if(missing(rel.power))
             fit=glm(dat ~ -1 + var1 + var2 + as.factor(lens),family=poisson)
           else
             fit=glm(dat ~ -1 + var1 + var2 + as.factor(lens) + offset(os),family=poisson)
           x=coef(fit)[c("var1","var2")]
           varx=summary(fit)$cov.unscaled[1:2,1:2]
           k=-2*x[2]/x[1]; sigma=sqrt(-2*x[2]/(x[1]^2))
           vartemp=msm::deltamethod(list(~-2*x2/x1,~sqrt(-2*x2/(x1^2))),x,varx,ses=F)
           pars=c(k,sigma,k*msizes[1],sigma)
           form1=as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
           varpars=msm::deltamethod(list(~x1,~x2,form1,~x2),c(k,sigma),vartemp,ses=F)
           gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
           rownames(gear.pars)=c("k","sigma","mode(mesh1)","std_dev(all meshes)") },
         "norm.sca"={
           if(missing(rel.power))
             fit=glm(dat ~ -1 + var3 + var4 + as.factor(lens),family=poisson)
           else
             fit=glm(dat ~ -1 + var3 + var4 + as.factor(lens) + offset(os),family=poisson)
           x=coef(fit)[c("var3","var4")]
           varx=summary(fit)$cov.unscaled[1:2,1:2]
           k1=-x[2]/(2*x[1]); k2=-1/(2*x[1])
           vartemp=msm::deltamethod(list(~-x2/(2*x1),~-1/(2*x1)),x,varx,ses=F)
           pars=c(k1,k2,k1*msizes[1],sqrt(k2*msizes[1]^2))
           form1=as.formula(sprintf("~x1*%f",msize1)) #Deltamethod quirk
           form2=as.formula(sprintf("~sqrt(x2*%f^2)",msize1)) #Deltamethod quirk
           varpars=msm::deltamethod(list(~x1,~x2,form1,form2),c(k1,k2),vartemp,ses=F)
           gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
           rownames(gear.pars)=c("k1","k2","mode(mesh1)","std_dev(mesh1)") },
         "gamma"={
           if(missing(rel.power))
             fit=glm(dat ~ -1 + var4 + var5 + as.factor(lens),family=poisson)
           else
             fit=glm(dat ~ -1 + var4 + var5 + as.factor(lens) + offset(os),family=poisson)
           x=coef(fit)[c("var4","var5")]
           varx=summary(fit)$cov.unscaled[1:2,1:2]
           alpha=x[2]+1; k=-1/x[1]
           vartemp=msm::deltamethod(list(~x2+1,~-1/x1),x,varx,ses=F)
           pars=c(alpha,k,(alpha-1)*k*msizes[1],sqrt(alpha*(k*msizes[1])^2))
           form1=as.formula(sprintf("~(x1-1)*x2*%f",msize1)) #Deltamethod quirk
           form2=as.formula(sprintf("~sqrt(x1*(x2*%f)^2)",msize1)) #Deltamethod quirk
           varpars=msm::deltamethod(list(~x1,~x2,form1,form2),c(alpha,k),vartemp,ses=F)
           gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
           rownames(gear.pars)=c("alpha","k","mode(mesh1)","std_dev(mesh1)")  },
         "lognorm"={
           if(missing(rel.power))
             fit=glm(dat ~ -1 + var6 + var7 + as.factor(lens),family=poisson)
           else
             fit=glm(dat ~ -1 + var6 + var7 + as.factor(lens) + offset(os),family=poisson)
           x=coef(fit)[c("var6","var7")]
           varx=summary(fit)$cov.unscaled[1:2,1:2]
           mu1=-(x[1]-1)/x[2]; sigma=sqrt(1/x[2])
           vartemp=msm::deltamethod(list(~-(x1-1)/x2,~sqrt(1/x2)),x,varx,ses=F)
           pars=c(mu1,sigma,exp(mu1-sigma^2),sqrt(exp(2*mu1+sigma^2)*(exp(sigma^2)-1)))
           varpars=msm::deltamethod(list(~x1,~x2,~exp(x1-x2^2),
                                         ~sqrt(exp(2*x1+x2^2)*(exp(x2^2)-1))),c(mu1,sigma),vartemp,ses=F)
           gear.pars=cbind(estimate=pars,s.e.=sqrt(diag(varpars)))
           rownames(gear.pars)=c("mu1(mode log-scale, mesh1)","sigma(std_dev log scale)",
                                 "mode(mesh1)","std_dev(mesh1)")  },
         stop(paste("\n",rtype, "not recognised, possible curve types are ",
                    "\"norm.loc\", \"norm.sca\", \"gamma\", and \"lognorm\"")))
  rselect=rcurves_Millar(rtype,meshsizes,rel.power,pars,plotlens)
  devres=matrix(resid(fit,type="deviance"),nrow(data),ncol(data[,-1]))
  # if(plots[1]) plot.curves(type,plotlens,rselect)
  # if(plots[2]) plot.resids(devres,meshsizes,data[,1])
  g.o.f=c(deviance(fit),sum(resid(fit,type="pearson")^2),fit$df.res,fit$null)
  names(g.o.f)=c("model_dev","Pearson chi-sq","dof","null_dev")
  fit.type=paste(paste(rtype,ifelse(is.null(rel.power),"",": with unequal mesh efficiencies")))
  if(details==FALSE)
    return(list(fit.type=fit.type,gear.pars=gear.pars,fit.stats=g.o.f))
  else
    return(list(
      fit.type=fit.type,gear.pars=gear.pars,fit.stats=g.o.f,devres=devres,rselect=rselect))
  
}

data("gillnet_k")


gillnetfit(data = millar,meshSizes,type="norm.loc",rel.power = NULL,
           plots=c(T,T),plotlens = NULL,details = F)



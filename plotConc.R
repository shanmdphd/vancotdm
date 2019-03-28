  
  # Chapter 7. plotCONC ----
  output$plotCONC <- renderPlot({
    
    data1=sim.data()
    
    plot(data1$TIME,data1$DV,ylab="concentration", xlab="time",type="l", pch=4,main="Concentration Curve of Vancomycin",col="#009ACD", ylim=c(0,50),xlim=c(0,data1$maxtime[1]))
    
    if (input$Observations=='1'){
      points(data1$pointtime,data1$observedConc, col="red", pch=6,cex=1.0)	#y point
      points(data1$pointtime[1],subset(data1$DV,data1$time==data1$pointtime), col="green", pch=6,cex=1.0)	#predicted DV point
      legend(data1$maxtime[1]-11,50,legend=c("observed conc.","predicted conc."),col=c("red","green"),pch=6:6, cex=0.8)
    }
    if(input$Observations=='2'){
      points(data1$pointtime1,data1$observedConc1, col="red", pch=6,cex=1.0)	#y point
      points(data1$pointtime1[1],subset(data1$DV,data1$time==data1$pointtime1), col="green", pch=6,cex=1.0)	#predicted DV point
      
      points(data1$pointtime2,data1$observedConc2, col="blue", pch=1,cex=1.0)   #y2 point
      points(data1$pointtime2[1],subset(data1$DV,data1$time==data1$pointtime2), col="black", pch=1,cex=1.0)   #predicted DV point 1
      legend(data1$maxtime[1]-11,50,legend=c("observed conc.","predicted conc."),col=c("red","green"),pch=6:6, cex=0.8)
    }
  })
  
  sim2.data <- reactive({
    
    obs1conc= input$obsc
    obs1time=input$obst
    obs1dat=input$obsd
    
    obs11conc= input$obsc1
    obs11time=input$obst1
    obs11dat=input$obsd1
    
    obs2conc= input$obsc2
    obs2time=input$obst2
    obs2dat=input$obsd2
    
    inFile <- input$file1
    WEIGHT <- input$weight
    AGE<-input$age
    SCR<-input$scr
    
    CLCR <- calculate_crcl(AGE, WEIGHT, input$sex, SCR)
    
    #	    Dose Data1 -------------------
    
    rawdata=read.csv(inFile$datapath,header=T)
    rawdata
    
    rawdata2=rawdata[complete.cases(rawdata), ]
    rawdata2
    
    infTime3=as.numeric(difftime(strptime(rawdata2[,3],"%H:%M"),strptime(rawdata2[,2],"%H:%M"),units="hours"))
    infTime3
    
    
    for(i in 1:length(infTime3)){
      if(infTime3[i]<0){
        infTime3[i]=24+infTime3[i]
      }
    }
    
    
    RATE3=rawdata2[,4]/infTime3
    RATE3
    
    cat("RATE3:",RATE3,"\n")
    
    before=strptime(paste(rawdata[,1],rawdata[,2]), "%y.%m.%d %H:%M")
    before
    after=strptime(paste(rawdata[,1],rawdata[,3]), "%y.%m.%d %H:%M")
    after
    
    
    
    if (input$Observations=='1')
    {
      
      Observeddate=paste(input$obsDate, substr(input$obsTime, 12, 20)) 
      pointtime=abs(as.numeric(difftime(before[1],Observeddate,units="hours")))	#pointtime
      
      dose2 = input$newdose
      
      infTime2=input$newinf					#infTime2=1	#infusion time
      infTime2
      
      
      
      # population INFORM. --------------------------------
      
      CLPOP = 2.82			#Creatine Level Compartment1
      CLPOP2 = 0.837			#Creatine Level Compartment2
      V1POP = 31.7			#Volume Compartment1
      QPOP = 11.8
      V2POP = 75.9			#Volume Compartment1
      
      # calculation	 --------------------------------------
      
      nd1=sum(cumprod(rawdata2[,4]>0))
      B1T=0						#Time of first dosing
      RATE2=dose2/infTime2				#dose by 1 time
      nd2=nd1*6
      newtau2=input$newtau
      
      result222=c()
      
      for (i in 1:nd1){
        result222[i]=abs(as.numeric(difftime(before[i],before[1],units="hours")))
        # print(result222)
        
      }
      
      tau=result222
      
      maxtime2=max(tau)+newtau2			#nd*tau
      maxtime2
      # cat("maxtime2:",maxtime2)
      
      i=c()
      v = RATE3
      
      result=c()
      
      for(i in 1:nd1){
        # print(i)
        result[i] = list(append(c(v[i]),0))
        #print(result)
        if(i==nd1)break;
      }
      
      
      
      
      
      result = c(do.call("cbind",result))
      result
      
      cat("result:", result,"\n")
      
      RATEinf<-append(result,0)		#RATE by point
      RATEinf
      flag<-complete.cases(RATEinf)
      RATEinf=RATEinf[flag]
      
      
      resultf=c()
      result2=c()
      result22=c()
      
      i=c()
      num=1
      a=append(tau,c(outer(max(tau),c(abs(newtau2)*(1:nd1)),"+")))
      b=infTime3
      
      cat("result:",result,"\n")
      cat("a:",a)
      
      
      for (i in 1:(nd1+1)){
        result22[i]=a[i]
        #cat("i=",i,"num=",num,"\n")
        #print(result22[i])
        resultf[num]=result22[i]
        
        
        if(i==nd1+1)break;
        
        num=num+1
        
        result2[num]=a[i]+b[i]
        #cat("i=",i,"num=",num,"\n")
        #print(result2[num])
        resultf[num]=result2[num]
        
        i=i+1
        num=num+1
        
      }
      
      
      
      TIMEinf=resultf
      flag<-complete.cases(TIMEinf)
      TIMEinf=TIMEinf[flag]
      TIMEinf
      
      
      TIMEinf2<-append(TIMEinf,outer(c(maxtime2+infTime2,maxtime2+newtau2), newtau2 * 0:ceiling((nd2/2)-1),"+")[-(length(outer(c(maxtime2+infTime2,maxtime2+newtau2), newtau2 * 0:ceiling((nd2/2)-1),"+")))])
      TIMEinf2
      
      RATEinf2<-append(c(rep(c(RATE2,0),nd2/2)),result,0)		#RATE by point
      RATEinf2
      flag<-complete.cases(RATEinf2)
      RATEinf2=RATEinf2[flag]
      
      length(TIMEinf2)
      length(RATEinf2)
      
      cat("a:",a,"\n")
      cat("newtau2:",newtau2,"\n")
      cat("TIMEinf2:",TIMEinf2,length(TIMEinf),"\n")
      cat("RATEinf2:",RATEinf2,length(RATEinf),"\n")
      cat("infTime3:",infTime3,length(infTime3),"\n")
      
      Cstepdoseinf <- approxfun(TIMEinf, RATEinf, method = "const")
      Cstepdoseinf(0:max(TIMEinf))
      length(Cstepdoseinf(0:max(TIMEinf)))
      
      
      Cstepdoseinf2 <- approxfun(TIMEinf2, RATEinf2, method = "const")
      Cstepdoseinf2(0:max(TIMEinf2))
      length(Cstepdoseinf2(0:max(TIMEinf2)))
      
      
      n=1
      ID = seq(from = 1, to = n, by = 1)
      TVCL=CLPOP*(CLCR/72)**CLPOP2
      TVV2 <- V2POP*(WEIGHT/60)
      TVV1 <- V1POP
      Q  <- QPOP
      
      B1T=0
      TIME <- seq(from = 0, to = maxtime2, by =0.1)
      TIME <- sort(unique(c(TIME,B1T)))
      TIMElast <- max(TIME)
      TIME3 <- seq(from = 0, to = max(TIMEinf2), by =0.1)
      
      ndoses=TIMElast/tau
    }
    
    
    RATE3=rawdata2[,4]/infTime3
    RATE3
    
    if (input$Observations=='2')
    {
      infTime2=input$newinf
      
      Observeddate1=paste(input$obsd1, substr(input$obst1, 12, 20))
      Observeddate2=paste(input$obsd2, substr(input$obst2, 12, 20))
      
      pointtime1=abs(as.numeric(difftime(before[1],Observeddate1,units="hours")))	#pointtime
      pointtime2=abs(as.numeric(difftime(before[1],Observeddate2,units="hours")))	#pointtime
      
      
      dose2 = input$newdose
      
      ##################################
      ######  population INFORM ----
      ##################################
      
      CLPOP = 2.82			#Creatine Level Compartment1
      CLPOP2 = 0.837			#Creatine Level Compartment2
      V1POP = 31.7			#Volume Compartment1
      QPOP = 11.8
      V2POP = 75.9			#Volume Compartment1
      
      
      # calculation	         ----
      
      nd1=sum(cumprod(rawdata2[,4]>0))	#ndoses=number of dosing
      # nd=sum(cumprod(rawdata2[,4]>0))*2	#ndoses=number of dosing
      RATE2=dose2/infTime2				#dose by 1 time
      B1T=0						#Time of first dosing
      nd2=nd1*6
      
      vec=c(Observeddate1,Observeddate2)
      
      Observeddatef=c()
      
      fun.obsdate<-function(vec){
        if(vec[1]>vec[2]) Observeddatef=vec[1]
        if(vec[1]<vec[2]) Observeddatef=vec[2]
        return(Observeddatef)
      }
      
      fun.obsdate(vec)
      
      pointtime3=abs(as.numeric(difftime(before[1],fun.obsdate(vec),units="hours")))   #pointtime
      
      #newtau2=ceiling(as.numeric(difftime(fun.obsdate(vec),before[nd1],units="hours")))	 
      newtau2=input$newtau
      
      result222=c()
      
      for (i in 1:nd1){
        result222[i]=abs(as.numeric(difftime(before[i],before[1],units="hours")))
        #print(result222)
      }
      
      tau=result222
      
      maxtime2=max(tau)+newtau2			#nd*tau
      maxtime2
      
      cat("maxtau=",max(tau))
      cat("maxtime2=",maxtime2)
      
      i=c()
      v = RATE3
      result=c()
      
      for(i in 1:nd1){
        print(i)
        result[i] = list(append(c(v[i]),0))
        print(result)
        if(i==nd1)break;
      }
      
      
      result = c(do.call("cbind",result))
      result
      
      RATEinf<-append(result,0)	#RATE by point
      RATEinf
      flag<-complete.cases(RATEinf)
      RATEinf=RATEinf[flag]
      RATEinf
      
      
      
      resultf=c()
      result2=c()
      result22=c()
      
      a=append(tau,c(outer(max(tau),c(abs(newtau2)*(1:nd1)),"+")))
      a
      i=c()
      num=1
      b=infTime3
      
      
      for (i in 1:(nd1+1)){
        result22[i]=a[i]
        #cat("i=",i,"num=",num,"\n")
        #print(result22[i])
        resultf[num]=result22[i]
        
        
        if(i==nd1+1)break;
        
        num=num+1
        
        result2[num]=a[i]+b[i]
        #cat("i=",i,"num=",num,"\n")
        #print(result2[num])
        resultf[num]=result2[num]
        
        i=i+1
        num=num+1
        
      }
      
      TIMEinf=resultf
      flag<-complete.cases(TIMEinf)
      TIMEinf=TIMEinf[flag]
      TIMEinf
      
      
      
      TIMEinf2<-append(TIMEinf,outer(c(maxtime2+infTime2,maxtime2+newtau2), newtau2 * 0:ceiling((nd2/2)-1),"+")[-(length(outer(c(maxtime2+infTime2,maxtime2+newtau2), newtau2 * 0:ceiling((nd2/2)-1),"+")))])
      TIMEinf2
      cat("nd2=",nd2,"\n")
      cat("maxtime2=", maxtime2,"\n")
      
      RATEinf2<-append(c(rep(c(RATE2,0),nd2/2)),result,0)		#RATE by point
      RATEinf2
      flag<-complete.cases(RATEinf2)
      RATEinf2=RATEinf2[flag]
      
      length(TIMEinf2)
      length(RATEinf2)
      
      cat("newtau2:",newtau2,length(newtau2),"\n")
      cat("TIMEinf2:",TIMEinf2,length(TIMEinf2),"\n")
      cat("RATEinf2:",RATEinf2,length(RATEinf2),"\n")
      cat("infTime32:",infTime3,length(infTime3),"\n")
      
      
      cat("TIMEinf:",TIMEinf,length(TIMEinf),"\n")
      cat("RATEinf:",RATEinf,length(RATEinf),"\n")
      
      
      Cstepdoseinf <- approxfun(TIMEinf, RATEinf, method = "const")
      Cstepdoseinf(0:max(TIMEinf))
      length(Cstepdoseinf(0:max(TIMEinf)))
      
      
      Cstepdoseinf2 <- approxfun(TIMEinf2, RATEinf2, method = "const")
      Cstepdoseinf2(0:max(TIMEinf2))
      length(Cstepdoseinf2(0:max(TIMEinf2)))
      
      n=1
      ID = seq(from = 1, to = n, by = 1)
      TVCL=CLPOP*(CLCR/72)**CLPOP2
      TVV2 <- V2POP*(WEIGHT/60)
      TVV1 <- V1POP
      Q  <- QPOP
      
      B1T=0
      TIME <- seq(from = 0, to = maxtime2, by =0.1)
      TIME <- sort(unique(c(TIME,B1T)))
      TIMElast <- max(TIME)
      #TIME2 <- seq(from = maxtime2+1, to=tau*10, by=0.1)
      TIME3 <- seq(from = 0, to = max(TIMEinf2), by =0.1)
      
      ndoses=TIMElast/tau
    }
    
    #Cstepdoseinf(T) function : please refer to  http://onlinelibrary.wiley.com/doi/10.1002/psp4.21/full#footer-support-info
    
    DOSEdata <- data.frame(var    = rep(1, times = nd1),
                           time   = tau, #seq(0,TIMElast-tau,tau),
                           value  = rep(0, times = nd1),
                           method = rep("add", times = nd1))
    
    n=1
    ID = seq(from = 1, to = n, by = 1)
    
    #ETA EPS
    EPS1SD <- 1.00E-04      #Additive residual error   
    EPS2SD <- 0.253     #Proportional residual error
    EPS2SDsq=(EPS2SD)^2
    ETA1SD <- 0.257    # Clearance eta 
    ETA2SD <- 0.959    # Volume eta
    
    cmat=function(vec){
      LENGTH=length(vec)
      DIM=as.integer(round((sqrt(8*LENGTH+1)-1)/2,0))
      if(DIM*(DIM+1)/2!=LENGTH) return(NULL)
      mat=matrix(nrow=DIM, ncol=DIM)
      mat[upper.tri(mat, diag=TRUE)]=vec		
      mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
      return(mat)
    }
    
    omega <- cmat(c(ETA1SD ,0,ETA2SD))	#c(Clearance eta,0,volume eta)
    omega.inv <- solve(omega)
    omega.inv
    
    # model 
    model <- function(T,A,eta){
      RateC <-Cstepdoseinf(T)
      K1=QPOP/V1POP
      K2=QPOP/(TVV2*exp(eta[2]))
      K3=(TVCL*exp(eta[1])) / TVV1
      
      dA <- vector(length = 2)
      dA[1] = (RateC) - (K1*A[1]) + (K2*A[2]) - (K3*A[1])  #Central compartment 
      dA[2] = (K1*A[1]) - (K2*A[2])				#Peripheral compartment 
      list(dA)
    }
    
    mod.cmp =cmpfun(model)
    
    # model 
    model2 <- function(T,A,eta){
      RateC <-Cstepdoseinf2(T)
      K1=QPOP/V1POP
      K2=QPOP/(TVV2*exp(eta[2]))
      K3=(TVCL*exp(eta[1])) / TVV1
      
      dA <- vector(length = 2)
      dA[1] = (RateC) - (K1*A[1]) + (K2*A[2]) - (K3*A[1])  #Central compartment 
      dA[2] = (K1*A[1]) - (K2*A[2])				#Peripheral compartment 
      list(dA)
    }
    
    mod2.cmp = cmpfun(model2)
    
    A_0=c(A1 = 0, A2 = 0)
    
    if (input$Observations=='1')
    {
      
      y=c(obs1conc)
      
      mapb22 <- function(eta){
        etamat=matrix(unlist(eta))
        out <- lsoda(c(A1=0 ,A2=0), TIME, mod.cmp, eta, events=list(data=DOSEdata))
        out <- cbind(out, DV=out[,"A1"]/TVV1)
        
        out31 <- lsoda(c(A1=0 ,A2=0), TIME3, mod2.cmp, eta, events=list(data=DOSEdata))
        out31 <- cbind(out31, DV=out31[,"A1"]/TVV1)
        
        eta=c(eta[1],eta[2])
        eta_m=unlist(matrix(eta,nrow = 2))
        sig2=EPS2SDsq
        sig2j <- subset(out[,4],out[,1]==pointtime)^2*sig2
        sqwres <- log(sig2j) + (1/sig2j)*(y-subset(out[,4],out[,1]==pointtime))^2
        nOn <- diag(t(eta_m) %*% omega.inv %*% eta_m)
        return(sum(sqwres)+ nOn)
      }
      
    }
    
    if (input$Observations=='2')
    {
      y=c(input$obsc1,input$obsc2)
      
      mapb22 <- function(eta){
        etamat=matrix(unlist(eta))
        out <- lsoda(c(A1=0 ,A2=0), TIME, mod.cmp, eta, events=list(data=DOSEdata))
        out <- cbind(out, DV=out[,"A1"]/TVV1)
        
        out31 <- lsoda(c(A1=0 ,A2=0), TIME3, mod2.cmp, eta, events=list(data=DOSEdata))
        out31 <- cbind(out31, DV=out31[,"A1"]/TVV1)
        
        eta=c(eta[1],eta[2])
        eta_m=unlist(matrix(eta,nrow = 2))
        sig2=EPS2SDsq
        sig2j <- subset(out[,4],out[,1]==pointtime1)^2*sig2
        sqwres <- log(sig2j) + (1/sig2j)*(y[1]-subset(out[,4],out[,1]==pointtime1))^2 + (1/sig2j)*(y[2]-subset(out[,4],out[,1]==pointtime2))^2
        nOn <- diag(t(eta_m) %*% omega.inv %*% eta_m)
        return(sum(sqwres)+ nOn)
      }
    }
    
    mapb22.cmp = cmpfun(mapb22)
    
    ini=c(0.1,0.1)
    
    
    withProgress(message = 'minimization in progress', min = 0, max = 100,value = 99,
                 
                 
                 {
                   
                   FIT=optim(ini,mapb22.cmp,control = list(trace=1,REPORT=1) ,method="L-BFGS-B")
                   FIT$par
                 }
    )
    
    
    # cat("FIT$par:",FIT$par,length(FIT$par),"\n")
    # cat("omega.inv:",omega.inv,length(omega.inv),"\n")
    
    troughLV=input$ll
    peakLV=input$ul
    
    
    out <- lsoda(c(A1=0 ,A2=0), TIME, mod.cmp, FIT$par, events=list(data=DOSEdata))
    out <- cbind(out, DV2=out[,"A1"]/TVV1)
    
    out31 <- lsoda(c(A1=0 ,A2=0), TIME3, mod2.cmp, FIT$par, events=list(data=DOSEdata))
    out31 <- cbind(out31, DV3=out31[,"A1"]/TVV1)
    
    outs<-data.frame(out)
    outs2<-data.frame(out31)
    
    # DV2 <- out[,"A1"]/TVV1
    # cat("is vector?",is.vector(DV2),"\n")
    # outs2$DV2 <- DV2
    
    outs$TIME=TIME
    outs2$TIME3=TIME3
    
    
    outs2$tau2=input$newtau
    outs2$maxtau=max(tau)
    outs2$maxtime2=maxtime2
    
    outs2$troughLV=troughLV
    outs2$peakLV=peakLV
    
    
    outs3=merge(x=outs,y=DOSEdata, by="time",all.x=TRUE)
    outs4=merge(x=outs2,y=DOSEdata, by="time",all.x=TRUE)
    outs5=merge(x=outs4,y=outs3, by="time",all.x=TRUE)
  })
  
  # Chapter 8. plotCONC2 ----
  output$plotCONC2 <- renderPlot({
    data2 = sim2.data()
    
    cat("data2$TIME=",data2$TIME,"\n")
    cat("data2$maxtau=",data2$maxtau,"\n")
    cat("data2$TIME3=",data2$TIME3,"\n")
    cat("data2$DV2=",data2$DV2,"\n")
    cat("data2$DV3=",data2$DV3,"\n")
    cat("data2$maxtime2=",data2$maxtime2,"\n")
    
    aaa=cbind(subset(data2$TIME,data2$TIME>data2$maxtau[1]), subset(data2$DV2,data2$TIME>data2$maxtau[1]))
    bbb=cbind(subset(data2$TIME3,data2$TIME3>data2$maxtime2[1]),subset(data2$DV3,data2$TIME3>data2$maxtime2[1]))
    fff=rbind(aaa,bbb)
    
    plot(fff[,1],fff[,2],ylab="concentration", xlab='time',type="l", pch=4,main="Concentration Curve of Vancomycin",col="#009ACD", ylim=c(0,50),lwd=2.0)#,xlim=c(0,200))
    
    abline(h=data2$troughLV, col='blue', lwd=3)
    abline(h=data2$peakLV,col='red',lwd=3, lty=2)
  })  

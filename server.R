# setup ----

default_dose_example_csv <- '"Date","Dosing_Time","Route","Dose"
"17.05.03","10:30","PO","500"
"17.05.04","10:30","PO","500"'

input <- list( # These values should be defined in ui.R.
  sex = 'male',
  obsc = 50, obsDate = "2017-05-05", obsTime = strptime("10:30", "%R"),
  weight = 60,
  Observations = '1',
  total_bilirubin = 1,
  scr = 0.9,
  post_op_date = 3,
  newdose = 1000,
  newtau = 48,
  newinf = 1,
  ll = 15,
  ul = 40,
  age = 20
)

library(deSolve)
library(plyr)
library(grid)
library(compiler)
library(shinyTime)
library(lubridate)
library(TeachingDemos)
library(rmarkdown)
library(knitr)
library(DT)
library(rsconnect)
library(tidyverse)

# dose in mg kg당 6-8mg : 8*60

calculate_crcl <- function(age, weight, sex, scr){
  crcl <- ((140-age) * weight * ifelse(sex == 'Female', 0.85, 1)) / (72*scr)
  return(crcl)
}

## ltv2mat copy right:: Prof. Bae 

cmat=function(vec){
  LENGTH=length(vec)
  DIM=as.integer(round((sqrt(8*LENGTH+1)-1)/2,0))
  if(DIM*(DIM+1)/2!=LENGTH) return(NULL)
  mat=matrix(nrow=DIM, ncol=DIM)
  mat[upper.tri(mat, diag=TRUE)]=vec
  mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
  return(mat)
}

# main ----

shiny::shinyServer(function(input, output) {
  
  # chapter 1. downloadData ----
  
  output$downloadData <- downloadHandler(
    filename <- function() { paste("dose_example", '.csv', sep='') },
    content <- function(file) {
      write.csv(datasetInput(), file, row.names = F)
    }
  )
  
  # Chapter 2. dosing_history_contents ----
  output$dosing_history_contents <- renderTable({
    if (is.null(input$file1)) return(read_csv(default_dose_example_csv,
                                              col_types = 'cccd'))
    return(read_csv(input$file1$datapath,
                    col_types = 'cccd'))
  })
  
  # Chapter 3. creatinine_clearance ----
  # crcl: https://www.mdcalc.com/creatinine-clearance-cockcroft-gault-equation
  output$creatinine_clearance <- renderText({
    return(calculate_crcl(input$age, input$weight, input$sex, input$scr) %>% round(digits = 2))
  })
  
  # Chapter 4. output_table1_time_predicted_concentration ----
  output$output_table1_time_predicted_concentration <- renderTable({
    prt1 <- sim.data()
    prt2 <- prt1[complete.cases(prt1),]
    if (input$Observations=='1') {
      prtx_predicted_concentration <- prt2 %>% 
        slice(2) %>% 
        select(pointtime, 
               `observed conc.`=observedConc, 
               `predicted conc.`=predictedConc)
    } 
    if (input$Observations=='2') {
      prtx_predicted_concentration <- prt2 %>% 
        slice(2) %>% 
        select(pointtime1, 
               observedConc1, 
               predictedConc1)
    }
    return(prtx_predicted_concentration)
  })
  
  # Chapter 5. outputtable2 ----
  output$outputtable2 <- renderTable({
    prt1=sim.data()
    prt2=prt1[complete.cases(prt1),]
    if (input$Observations=='1')
    {
      prtx=prt2[2,c("CL","V1","V2")]
    }
    
    if (input$Observations=='2')
    {
      prtx=prt2[2,c("pointtime2","observedConc2","predictedConc2")]
    }
    prtx
  })
  
  # Chapter 6. outputtable3 ----
  output$outputtable3 <- renderTable({
    prt1=sim.data()
    prt2=prt1[complete.cases(prt1),]
    if(input$Observations=='1'){
    }
    if (input$Observations=='2')
    {
      prtx=prt2[2,c("CL","V1","V2")]
      return(prtx)
    }
  })
  
  # prelude 1. datasetInput ----
  
  datasetInput <- reactive({
    coln=c("Date","Inf_st_Time","Inf_ed_Time","Dose" )
    dat1=c("17.05.03", "10:30", "11:30","500")
    dat2=c("17.05.03", "22:30", "23:30","500")
    dat3=c("17.05.04", "10:30", "11:30","1000")
    dat4=c("17.05.04", "22:30", "23:30","1000")
    
    suppl <- data.frame (rbind(dat1,dat2,dat3,dat4))
    colnames(suppl)<-coln
    rownames(suppl)<-NULL
    suppl
  })
  
  # prelude 2. dose.data ----
  
  dose.data <- reactive({
    inFile <- input$file1 
    if (is.null(inFile)) return(NULL) 
    a=read.csv(inFile$datapath, header=T, stringsAsFactors = T) 
    b=a[complete.cases(a), ]  
    b$paste=paste(b$Date,b$Time) 
    b
    #Time calculation code is copyrighted
  })
  
  sim.data <- reactive({
    
    # prelude 3. sim.data ----
    
    # input: Observation
    obs1conc <- input$obsc
    obs1time <- input$obst
    obs1dat  <- input$obsd
    
    # input: Demog
    WEIGHT   <- input$weight
    AGE      <- input$age
    TBIL     <- input$total_bilirubin
    POD_week <- input$post_op_date/7
    
    # Typical Values
    TVCL <- 28.5 - 1.24*POD_week - 0.252*(TBIL-10) + 0.188*(WEIGHT-60) - 0.191*(AGE - 40)
    TVV1 <- 133.00
    TVKA <-   1.28
    
    # Eta (Omega)
    ETA1SD <- 0.04     # Vancomycin Clearance eta 
    ETA2SD <- 0.04     # Vancomycin Volume eta
    omega <- cmat(c(ETA1SD ,0,ETA2SD)) %>%  # Prof. Bae's function
      print() #c(Clearance eta,0,volume eta)
    omega.inv <- solve(omega) %>% 
      print()
    
    # Eps (Sigma)
    
    EPS1SD <- 40      # Vancomycin Additive residual error   
    EPS2SD <- 0.1     # Vancomycin Proportional residual error
    EPS2SDsq=(EPS2SD)^2 # Vancomycin square
    
    WEIGHT;AGE;TBIL;POD_week;TVCL;TVV1;TVKA;omega;omega.inv
    
    input_file_text <- ifelse(is.null(input$file1), 
                              yes = default_dose_example_csv, 
                              no = input$file1$datapath)
    
    rawdata <- read_csv(input_file_text, col_types = 'cccd') %>% 
      as.data.frame() %>% 
      print()
    
    rawdata2 <- rawdata[complete.cases(rawdata), ]
    
    # infusion duration
    infTime3 <- as.numeric(difftime(strptime(rawdata2[,3],"%H:%M"), # end
                                    strptime(rawdata2[,2],"%H:%M"), # start
                                    units="hours")) %>% 
      print()
    
    for(i in 1:length(infTime3)){
      if(infTime3[i]<0){
        infTime3[i]=24+infTime3[i]
      }
    }
    
    RATE3=rawdata2[,4]/infTime3
    RATE3
    
    before=strptime(paste(rawdata[,1],rawdata[,2]), "%y.%m.%d %H:%M")
    before
    after=strptime(paste(rawdata[,1],rawdata[,3]), "%y.%m.%d %H:%M")
    after
    
    if (input$Observations=='1')
    {
      Observeddate=paste(input$obsDate, substr(input$obsTime, 12, 20))
      pointtime=abs(as.numeric(difftime(before[1],Observeddate,units="hours")))  #pointtime
      
      dose2 = input$newdose
      
      # calculation
      
      nd1=sum(cumprod(rawdata2$Dose>0))  #ndoses=number of dosing
      nd=sum(cumprod(rawdata2$Dose>0))*2  # number of dosing x 2
      B1T=0            #Time of first dosing 
      nd2=nd*3 # number of dosing x 3
      
      # first dose to the last dose
      newtau=ceiling(as.numeric(difftime(Observeddate,before[nd1],units="hours")))   #time duration
      
      result222=c()
      
      for (i in 1:nd1){
        result222[i]=abs(as.numeric(difftime(before[i],before[1],units="hours")))
        #print(result222)
      }
      
      tau=result222
      
      maxtime=c(max(tau)+newtau)      #nd*tau
      maxtime
      
      i=c()
      v = RATE3
      result=c()
      
      for(i in 1:(nd)){
        print(i)
        result[i] = list(append(c(v[i]),0))
        #print(result)
        if(i==nd)break;
      }
       
      result = c(do.call("cbind",result))
      result
      
      RATEinf <- result  #RATE by point
      RATEinf
      
      flag <- complete.cases(RATEinf)
      RATEinf=RATEinf[flag]
      RATEinf
      
      resultf=c()
      result2=c()
      result22=c()
      
      a=append(tau,c(outer(max(tau),c(abs(newtau)*(1:nd1)),"+")))
      a
      i=c()
      num=1
      b=infTime3
      
      for (i in 1:(nd+1)){
        result22[i]=a[i]
        #cat("i=",i,"num=",num,"\n")
        #print(result22[i])
        resultf[num]=result22[i]
        
        if(i==nd+1)break;
        num=num+1
        result2[num]=a[i]+b[i]
        #cat("i=",i,"num=",num,"\n")
        #print(result2[num])
        resultf[num]=result2[num]
        
        i=i+1
        num=num+1 
      }
      
      TIMEinf=resultf
      TIMEinf
      
      flag<-complete.cases(TIMEinf)
      TIMEinf=TIMEinf[flag]
      list(TIMEinf, RATEinf) %>% map(length)
      
      Cstepdoseinf <- approxfun(TIMEinf, RATEinf, method = "constant")
      Cstepdoseinf(0:max(TIMEinf))
      Cstepdoseinf(0)
      length(Cstepdoseinf(0:max(TIMEinf)))
      
      Cstepdoseinf(c(0.1,0.2,1))
       
      n=1
      ID = seq(from = 1, to = n, by = 1)
      # TVCL=CLPOP*(CLCR/72)**CLPOP2
      # TVV2 <- V2POP*(WEIGHT/60)
      # TVV1 <- V1POP
      # Q  <- QPOP
      
      B1T=0
      TIME <- seq(from = 0, to = pointtime, by = 0.1)
      TIME <- sort(unique(c(TIME,B1T)))
      TIMElast <- max(TIME)
    }
    
    # Cstepdoseinf(T) function : 
    # please refer to  http://onlinelibrary.wiley.com/doi/10.1002/psp4.21/full#footer-support-info
    
    DOSEdata <- data.frame(var    = rep(1, times = nd1),
                           time   = tau,  #seq(0,TIMElast-tau,tau),
                           value  = rep(500, times = nd1),
                           method = rep("add", times = nd1))
    
    # model ----
    
    model <- function(Time, A, eta){
      Ka <- TVKA
      Cl <- TVCL * exp(eta[1])
      Vd <- TVV1 * exp(eta[2])
      
      Ke <-  Cl/Vd  # Elimination
      
      dA <- vector(length = 2)
      dA[1] <- -Ka * A[1]              # Depot compartment - intestine 
      dA[2] <-  Ka * A[1] - Ke * A[2]   # Central compartment - plasma
      return(list(dA))
    }
    
    mod.cmp <-  compiler::cmpfun(model)
    initial_amount <- c(A1 = 0, A2 = 0) 
    if (input$Observations=='1') {
      y <- c(obs1conc)
      #' @example mapb2(c(0.1, 0.1))
      mapb2 <- function(eta){ # eta is a list of 2. eta <- c(0.4321,0.4321)
        etamat=matrix(unlist(eta))
        out <- deSolve::lsoda(y = initial_amount, 
                              times = TIME, 
                              func = model, 
                              parms = eta, 
                              events=list(data=DOSEdata))
        out <- cbind(out, DV=out[,"A2"]/TVV1)
        head(out)
        
        eta <- c(eta[1],eta[2])
        eta_m <- unlist(matrix(eta,nrow = 2))
        sig2 <- EPS2SDsq
        sig2j <- subset(out[,3],out[,1]==pointtime)^2*sig2
        sqwres <- log(sig2j) + (1/sig2j) * (y[1]-subset(out[,3],out[,1]==pointtime))^2
        nOn <- diag(t(eta_m) %*% omega.inv %*% eta_m)
        return(sum(sqwres)+ nOn)
      } 
    }
    
    mapb2.cmp <- cmpfun(mapb2)
    ini <- c(0.04,0.04)
    mapb2.cmp(ini)
    
    # Method "L-BFGS-B" is that of Byrd et. al. (1995) which allows box constraints, that is each variable can be given a lower and/or upper bound. The initial value must satisfy the constraints. This uses a limited-memory modification of the BFGS quasi-Newton method. If non-trivial bounds are supplied, this method will be selected, with a warning.
    
    #shiny::withProgress(
    #  message = 'Minimization in progress', 
    #  min = 0, 
    #  max = 100,
    #  value = 99, 
    #  {
    FIT <- stats::optim(par = ini, # CL, V => fitting
                        fn = mapb2.cmp, # A function to be minimized (or maximized)
                        method="L-BFGS-B", 
                        control = list(trace=TRUE,REPORT=TRUE))
    print(FIT$par)
    
    #  }
    #)
    #FIT <- list(par=c(-0.15379907 , 0.08570668))
    # cat("FIT$par=",FIT$par)
    # cat("omega.inv=",omega.inv)
    
    outs <- lsoda(y = initial_amount, 
                  times = TIME, 
                  func = mod.cmp,  
                  parms = FIT$par,
                  events = list(data = DOSEdata )) %>% 
      as.data.frame() %>% 
      mutate(DV = A2/TVV1) %>% 
      mutate(maxtime = maxtime) %>% 
      print()
    
    #out <- cbind(out, DV=out[,"A1"]/TVV1)
    #outs <- data.frame(out)
    #outs$maxtime=maxtime
    
    if (input$Observations=='1'){
      outs$pointtime=pointtime  
      outs$Observeddate=Observeddate
      outs$predictedConc=subset(outs[,3],outs[,1]==pointtime)
      outs$observedConc=obs1conc
    }
    
    outs$TIME <- TIME
    outs$CL <- TVCL*(exp(FIT$par[1]))    #predicted CL
    outs$V1 <- TVV1*(exp(FIT$par[2]))    #predicted V2
    
    # 
    #     if (input$Observations=='1')
    #     {
    #       
    #     }
    #     if (input$Observations=='2')
    #     {
    #       outs$obs1=obs11conc
    #       outs$obs2=obs2conc
    #     }
    
    outs2=merge(x=outs,y=DOSEdata, by="time",all.x=TRUE)
    head(outs2);tail(outs2)
    outs2 %>% as_tibble()
    outs2 %>% as_tibble() %>% 
      ggplot(aes(time, DV)) + geom_point() + geom_hline(yintercept = c(0.050, 0.200))
  })
  # end ----  
})

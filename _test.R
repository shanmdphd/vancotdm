
library(shinydashboard)
library(shinyTime)
library(lubridate)
library(shinythemes)
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
library(shiny)  


default_dose_example_csv <- '"Date","Inf_st_Time","Inf_ed_Time","Dose"
"17.05.03","10:30","11:30","500"
"17.05.03","22:30","23:30","750"
"17.05.04","10:30","11:30","1000"
"17.05.04","22:30","23:30","1000"'

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

input <- list(
  sex = 'male',
  obsc = 5,
  obst = '2300',
  obsd = '2017-05-06',
  weight = 60,
  Observations = '1',
  obsDate = "2017-05-06",
  obsTime = strptime("23:00", "%R"),
  scr = 0.9,
  total_bilirubin = 1,
  post_op_date = 3,
  age = 20
)

input$file1$datapath <- 'dose_example.csv'

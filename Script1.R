#SCRIPT 1. Leverage calculation

#This script performs a standard leverage analysis for the three bioregions that
#cover the Southern Forests of SW Australia. Leverage is the inverse of the slope
#for a simple linear regression of wildfire area as a function of the area of
#recently burnt land. As per Boer et al (2009), all wildfire recorded for each
#region is used.

library(dplyr)

#Define start and end years for analysis
start <- 1954
end <- 2018
young <- 6

#Include miscellaneous fires as WF? (Y/N)
misc <- "N"

####################################
SW <- read.csv("SW1950.csv")
Year <- seq(start, end, 1)
tbl <- as.data.frame(Year)

names <- as.data.frame(distinct(SW,Bioregion))
n <- as.numeric(nrow(names))

#Create empty summary dataframe
res <- data.frame( 'Bioregion' = character(0), 'Leverage' = numeric(0), 'Rsquared' = numeric(0), 'p'= numeric(0))

#Loop through regions
for(N in 1:n) {
  
  #Isolate records to the bioregion
  region <- as.character(names$Bioregion[N])
  x <- SW[SW$Bioregion == region, ]
  
  #Find area of WF per year
  if (misc == "Y") {
    wf <- x[x$Binary == "WF", ]
  } else {
    wf <- x[x$Type == "WF", ] 
  }

  wf <- wf[ which( wf$Year >= start & wf$Year <= end) , ]%>%
    group_by(Year) %>%
    summarize_if(is.numeric, sum)%>%
    mutate(WF = Area_ha)%>%
    select(Year, WF)%>%
    left_join(tbl)
  
  #Find area of young per year
  y <- x[ which( x$Year >= start & x$Year <= end) , ]%>%
    group_by(Year) %>%
    summarize_if(is.numeric, sum)%>%
    mutate(Tot = Area_ha)%>%
    select(Year, Tot)%>%
    left_join(tbl)%>%
    mutate(Recent = RcppRoll::roll_sum(Tot, n = young, fill=NA, align="right"))%>%
    mutate(Year = Year+1)

  #Join tables
  vals <- y %>%
    left_join(wf)%>%
    select(Year, Recent, WF)
  vals[is.na(vals)] <- 0
  vals <- vals[ which( vals$Year >= start+young & vals$Year <= end) , ]
  
  #Fit linear
  Linm<-lm(vals$WF ~ vals$Recent)
  test <- cor.test(vals$Recent,vals$WF)
  leverage <- round(-(lm(vals$WF ~ vals$Recent)$coefficients[2]),3)
  p <- round(test$p.value, 3)
  Rsq <- round(test$estimate^2, 3)
  a <- round((lm(vals$WF ~ vals$Recent)$coefficients[2]),4) 
  b <- round((lm(vals$WF ~ vals$Recent)$coefficients[1]),0)
  out <- as.data.frame(list('Bioregion' = region, 'Leverage' = leverage, 'Rsquared' = Rsq, 'p'=p, 'a'=a, 'b'=b))
  
  #Add to table
  res <- rbind(res, out)
  #Name table of raw bioregion res
  assign(region,vals)
}

write.csv(res, "Table_S2.csv")

#SCRIPT 2. Climatic drivers and leverage recalculation
#Part 1 Searches for correlations between wildfire and potential climatic drivers
#Constructs simple models for the most likely drivers

library(ggplot2)
library(patchwork)
library(tidyverse)
library(splines)

#CLIMATIC DATA
SA <- read.csv("climFireDat.csv")%>%
  mutate(synVar = HdA+LdA,
         wWF = Warren)

#####################################################################
testSet <- c("Year", "wWF", "SW", "SWf", "SA", "SAf")
colNames <- c("Year", "wWF", "SW", "SWf", "SA", "SAf", "Year_p", "wWF_p", "SW_p", "SWf_p", "SA_p", "SAf_p")
n <- names(SA)
Climatic <- n[10:length(n)]
res2 <- as.data.frame(Climatic)
for(i in colNames){
  res2[,i] <- NA
}

#Loop through columns finding correlations
cola <- 2
colb <- 8
for (test in 1:length(testSet)) {
  refN <- which( colnames(SA)==testSet[test] )
 for (row in 1:length(Climatic)) {
  r1 <- (cor.test(SA[ ,refN],SA[,(row+9)]))
  p1 <- round(as.numeric(r1$p.value), 3)
  if (p1 <= 0.05) {
    c1 <- round(as.numeric(r1$estimate),3)
  } else {
    c1 <- ""
    p1 <- ""
  }
  res2[row,cola] <- c1
  res2[row,colb] <- p1
 }
  cola <- cola+1
  colb <- colb+1
}

write.csv(res2, "Table_S3.csv")

#Create a subset of Warren predictors
library(dplyr)
subWarren <- SA %>%
  select("wWF", "WarrenY", "TmeanA", "TmaxA", "TminS", "synVar", "LiS")

cortab <- as.data.frame(round(cor(subWarren),3))
write.csv(cortab, "Table_S4.csv")

# Create Fig. S4
windows(10,10)
pairs(subWarren)


##############################
# Part 2. Calculate leverage for Warren Bioregion, including climatic drivers
# NB: Basic leverage result will differ slightly from that in script 1 as the starting year is 1957 rather than 1954


# Build a function to extract a p value from linear model object
# Source: https://gettinggeneticsdone.blogspot.com/2011/01/rstats-function-for-extracting-f-test-p.html
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

colNames <- c("p", "Rsq", "Leverage", "Age_p", "AIC")
Model <- c("Y", "synVar", "TmaxA", "LiS", 
           "Y+synVar", "Y+TmaxA", "Y+LiS", 
           "synVar+TmaxA", "synVar+LiS", "TmaxA+LiS", 
           "Y+synVar+TmaxA", "Y+synVar+LiS", "Y+TmaxA+LiS",
           "synVar+TmaxA+LiS", "Y+synVar+TmaxA+LiS")
res3 <- as.data.frame(Model)
for(i in colNames){
  res3[,i] <- NA
}

# Y
Mod<-lm(SA$wWF ~ SA$WarrenY)
Ma <- summary(Mod)
Ma_r <- round(as.numeric(Ma$r.squared),3)
Ma_Yp <- round(as.numeric(Ma$coefficients[8]),3)
if (Ma_Yp < 0.05) {
  Ma_YL <- -1*round(as.numeric(Ma$coefficients[2]),3)  
} else {
  Ma_YL <- "N.S."
}
Ma_p <- round(lmp(Mod),3)
Ma_AIC <- round(AIC(Mod),1)
res3[1,2] <- Ma_p
res3[1,3] <- Ma_r
res3[1,4] <- Ma_YL
res3[1,5] <- Ma_Yp
res3[1,6] <- Ma_AIC

# synVar
Mod<-lm(SA$wWF ~ SA$synVar)
Mb <- summary(Mod)
Mb_r <- round(as.numeric(Mb$r.squared),3)
Mb_p <- round(lmp(Mod),3)
Mb_AIC <- round(AIC(Mod),1)
res3[2,2] <- Mb_p
res3[2,3] <- Mb_r
res3[2,6] <- Mb_AIC


# TmaxA
Mod<-lm(SA$wWF ~ SA$TmaxA)
Mc <- summary(Mod)
Mc_r <- round(as.numeric(Mc$r.squared),3)
Mc_p <- round(lmp(Mod),3)
Mc_AIC <- round(AIC(Mod),1)
res3[3,2] <- Mc_p
res3[3,3] <- Mc_r
res3[3,6] <- Mc_AIC


# LiS
Mod<-lm(SA$wWF ~ SA$LiS)
Md <- summary(Mod)
Md_r <- round(as.numeric(Md$r.squared),3)
Md_p <- round(lmp(Mod),3)
Md_AIC <- round(AIC(Mod),1)
res3[4,2] <- Md_p
res3[4,3] <- Md_r
res3[4,6] <- Md_AIC

# Y+synVar
Mod<-lm(SA$wWF ~ SA$WarrenY + SA$synVar)
Mab <- summary(Mod)
Mab_r <- round(as.numeric(Mab$r.squared),3)
Mab_Yp <- round(as.numeric(Mab$coefficients[11]),3)
if (Mab_Yp < 0.05) {
  Mab_YL <- -1*round(as.numeric(Mab$coefficients[2]),3)  
} else {
  Mab_YL <- "N.S."
}
Mab_p <- round(lmp(Mod),3)
Mab_AIC <- round(AIC(Mod),1)
res3[5,2] <- Mab_p
res3[5,3] <- Mab_r
res3[5,4] <- Mab_YL
res3[5,5] <- Mab_Yp
res3[5,6] <- Mab_AIC

# Y+TmaxA
Mod<-lm(SA$wWF ~ SA$WarrenY + SA$TmaxA)
Mac <- summary(Mod)
Mac_r <- round(as.numeric(Mac$r.squared),3)
Mac_Yp <- round(as.numeric(Mac$coefficients[11]),3)
if (Mac_Yp < 0.05) {
  Mac_YL <- -1*round(as.numeric(Mac$coefficients[2]),3)  
} else {
  Mac_YL <- "N.S."
}
Mac_p <- round(lmp(Mod),3)
Mac_AIC <- round(AIC(Mod),1)
res3[6,2] <- Mac_p
res3[6,3] <- Mac_r
res3[6,4] <- Mac_YL
res3[6,5] <- Mac_Yp
res3[6,6] <- Mac_AIC

# Y+LiS
Mod<-lm(SA$wWF ~ SA$WarrenY + SA$LiS)
Mad <- summary(Mod)
Mad_r <- round(as.numeric(Mad$r.squared),3)
Mad_Yp <- round(as.numeric(Mad$coefficients[11]),3)
if (Mad_Yp < 0.05) {
  Mad_YL <- -1*round(as.numeric(Mad$coefficients[2]),3)  
} else {
  Mad_YL <- "N.S."
}
Mad_p <- round(lmp(Mod),3)
Mad_AIC <- round(AIC(Mod),1)
res3[7,2] <- Mad_p
res3[7,3] <- Mad_r
res3[7,4] <- Mad_YL
res3[7,5] <- Mad_Yp
res3[7,6] <- Mad_AIC

# synVar+TmaxA
Mod<-lm(SA$wWF ~ SA$synVar + SA$TmaxA)
Mbc <- summary(Mod)
Mbc_r <- round(as.numeric(Mbc$r.squared),3)
Mbc_p <- round(lmp(Mod),3)
Mbc_AIC <- round(AIC(Mod),1)
res3[8,2] <- Mbc_p
res3[8,3] <- Mbc_r
res3[8,6] <- Mbc_AIC

# synVar+LiS
Mod<-lm(SA$wWF ~ SA$synVar + SA$LiS)
Mbd <- summary(Mod)
Mbd_r <- round(as.numeric(Mbd$r.squared),3)
Mbd_p <- round(lmp(Mod),3)
Mbd_AIC <- round(AIC(Mod),1)
res3[9,2] <- Mbd_p
res3[9,3] <- Mbd_r
res3[9,6] <- Mbd_AIC

# TmaxA+LiS
Mod<-lm(SA$wWF ~ SA$TmaxA + SA$LiS)
Mcd <- summary(Mod)
Mcd_r <- round(as.numeric(Mcd$r.squared),3)
Mcd_p <- round(lmp(Mod),3)
Mcd_AIC <- round(AIC(Mod),1)
res3[10,2] <- Mcd_p
res3[10,3] <- Mcd_r
res3[10,6] <- Mcd_AIC

# Y+synVar+TmaxA
Mod<-lm(SA$wWF ~ SA$WarrenY + SA$synVar + SA$TmaxA)
Mabc <- summary(Mod)
Mabc_r <- round(as.numeric(Mabc$r.squared),3)
Mabc_Yp <- round(as.numeric(Mabc$coefficients[14]),3)
if (Mabc_Yp < 0.05) {
  Mabc_YL <- -1*round(as.numeric(Mab$coefficients[2]),3)  
} else {
  Mabc_YL <- "N.S."
}
Mabc_p <- round(lmp(Mod),3)
Mabc_AIC <- round(AIC(Mod),1)
res3[11,2] <- Mabc_p
res3[11,3] <- Mabc_r
res3[11,4] <- Mabc_YL
res3[11,5] <- Mabc_Yp
res3[11,6] <- Mabc_AIC

# Y+synVar+LiS
Mod<-lm(SA$wWF ~ SA$WarrenY + SA$synVar + SA$LiS)
Mabd <- summary(Mod)
Mabd_r <- round(as.numeric(Mabd$r.squared),3)
Mabd_Yp <- round(as.numeric(Mabd$coefficients[14]),3)
if (Mabd_Yp < 0.05) {
  Mabd_YL <- -1*round(as.numeric(Mabd$coefficients[2]),3)  
} else {
  Mabd_YL <- "N.S."
}
Mabd_p <- round(lmp(Mod),3)
Mabd_AIC <- round(AIC(Mod),1)
res3[12,2] <- Mabd_p
res3[12,3] <- Mabd_r
res3[12,4] <- Mabd_YL
res3[12,5] <- Mabd_Yp
res3[12,6] <- Mabd_AIC

# Y+TmaxA+LiS
Mod<-lm(SA$wWF ~ SA$WarrenY + SA$TmaxA + SA$LiS)
Macd <- summary(Mod)
Macd_r <- round(as.numeric(Macd$r.squared),3)
Macd_Yp <- round(as.numeric(Macd$coefficients[14]),3)
if (Macd_Yp < 0.05) {
  Macd_YL <- -1*round(as.numeric(Macd$coefficients[2]),3)  
} else {
  Macd_YL <- "N.S."
}
Macd_p <- round(lmp(Mod),3)
Macd_AIC <- round(AIC(Mod),1)
res3[13,2] <- Macd_p
res3[13,3] <- Macd_r
res3[13,4] <- Macd_YL
res3[13,5] <- Macd_Yp
res3[13,6] <- Macd_AIC

# synVar+TmaxA+LiS
Mod<-lm(SA$wWF ~ SA$synVar + SA$TmaxA + SA$LiS)
Mbcd <- summary(Mod)
Mbcd_r <- round(as.numeric(Mbcd$r.squared),3)
Mbcd_p <- round(lmp(Mod),3)
Mbcd_AIC <- round(AIC(Mod),1)
res3[14,2] <- Mbcd_p
res3[14,3] <- Mbcd_r
res3[14,6] <- Mbcd_AIC

# Y+synVar+TmaxA+LiS
Mod<-lm(SA$wWF ~ SA$WarrenY + SA$synVar + SA$TmaxA + SA$LiS)
Mabcd <- summary(Mod)
Mabcd_r <- round(as.numeric(Mabcd$r.squared),3)
Mabcd_Yp <- round(as.numeric(Mabcd$coefficients[17]),3)
if (Mabcd_Yp < 0.05) {
  Mabcd_YL <- -1*round(as.numeric(Mabcd$coefficients[2]),3)  
} else {
  Mabcd_YL <- "N.S."
}
Mabcd_p <- round(lmp(Mod),3)
Mabcd_AIC <- round(AIC(Mod),1)
res3[15,2] <- Mabcd_p
res3[15,3] <- Mabcd_r
res3[15,4] <- Mabcd_YL
res3[15,5] <- Mabcd_Yp
res3[15,6] <- Mabcd_AIC

write.csv(res3, "Table_S5.csv")

# Graph the climatic model for Warren bioregion

SA <- SA%>%
  mutate(warrenMod = 8337*synVar+2616*LiS-2713795,
         warRes = wWF - warrenMod)

names(SA)

# Create Fig. 1c
windows(4.7,4.5)
ggplot(SA , aes(y = warRes , x = WarrenY))+
  geom_point(alpha = 0.4, size = 2, colour = "blue")+ 
  ggtitle("c.") +
  geom_smooth(method = "lm", formula = y ~ bs(x, degree=1), color= "firebrick",
              fill = "firebrick", alpha = 0.4) +
  labs(x = "Area of young forest (ha)", y = "Residual WF area (ha)" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=11, colour = "black"),
        axis.text.y  = element_text(size=11, colour = "black"),
        axis.title.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=14),
        plot.title = element_text(vjust=1.5, size=14, colour = "black"))

#Find Rsq for leverage against climate-corrected data
cor <- cor.test(SA$warRes, SA$WarrenY)
Rsq <- as.numeric(cor$statistic)^2


# Fig S5
windows(5,4)
ggplot(SA , aes(y = wWF , x = warrenMod, colour = Year))+
  geom_point(alpha = 0.6, size = 2)+ 
  geom_smooth(method = "lm", formula = y ~ bs(x, degree=3), color= "firebrick",
              fill = "firebrick", alpha = 0.4) +
  labs(x = "Modelled burnt area (ha)", y = "Measured burnt area (ha)" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=11, colour = "black"),
        axis.text.y  = element_text(size=11, colour = "black"),
        axis.title.y = element_text(face="bold", size=11),
        axis.title.x = element_text(face="bold", size=11))


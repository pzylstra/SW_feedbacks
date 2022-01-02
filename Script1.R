# SCRIPT 1

# Identifies climatic correlations & autocorrelations (Fig. S1, Table S2)
# Builds and tests climatic models (Table S3)
# Groups years into tertile divisions and creates Fig. 2

library(ggplot2)
library(patchwork)
library(tidyverse)
library(splines)

#LOAD CLIMATIC DATA
SA <- read.csv("climFireDat.csv")%>%
  mutate(synVar = HdA+LdA)

#####################################################################
# Part 1 - Identify climatic drivers

testSet <- c("Year", "WF")
colNames <- c("Year", "WF", "Year_p", "WF_p")
n <- names(SA)
Climatic <- n[3:length(n)]
res1 <- as.data.frame(Climatic)
for(i in colNames){
  res1[,i] <- NA
}

#Loop through columns finding correlations
cola <- 2
colb <- 4
for (test in 1:length(testSet)) {
  refN <- which( colnames(SA)==testSet[test] )
  for (row in 1:length(Climatic)) {
    r1 <- (cor.test(SA[ ,refN],SA[,(row+2)]))
    p1 <- round(as.numeric(r1$p.value), 3)
    if (p1 <= 0.05) {
      c1 <- round(as.numeric(r1$estimate),3)
    } else {
      c1 <- ""
      p1 <- ""
    }
    res1[row,cola] <- c1
    res1[row,colb] <- p1
  }
  cola <- cola+1
  colb <- colb+1
}

write.csv(res1, "Table_S2.csv")

# Examine autocorrelations

SA <- read.csv("climFireDat.csv") %>%
  mutate(synVar = HdA+LdA)%>%
  select("Year", "WF", "TmeanA", "TmaxA", "HdA", "HdS", "LdA", "LdS", "synVar")

cortabSA <- as.data.frame(round(cor(SA),3))
write.csv(cortabSA, "Table_S2.csv")

# Fig. S1
windows(10,10)
pairs(SA)


#_____________________________________________________________________
# Part 2. Build climatic models

# Build a function to extract a p value from linear model object
# Source: https://gettinggeneticsdone.blogspot.com/2011/01/rstats-function-for-extracting-f-test-p.html
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

colNames <- c("p", "Rsq", "AIC")
Model <- c("synVar", "TmeanA", "TmaxA",  
           "synVar+TmeanA", "synVar+TmaxA", "TmaxA+TmeanA",
           "synVar+TmeanA+TmaxA")
res2 <- as.data.frame(Model)
for(i in colNames){
  res2[,i] <- NA
}

# synVar
Mod<-lm(SA$WF ~ SA$synVar)
Ma <- summary(Mod)
Ma_r <- round(as.numeric(Ma$r.squared),3)
Ma_p <- round(lmp(Mod),3)
Ma_AIC <- round(AIC(Mod),1)
res2[1,2] <- Ma_p
res2[1,3] <- Ma_r
res2[1,4] <- Ma_AIC

# TmeanA
Mod<-lm(SA$WF ~ SA$TmeanA)
Mb <- summary(Mod)
Mb_r <- round(as.numeric(Mb$r.squared),3)
Mb_p <- round(lmp(Mod),3)
Mb_AIC <- round(AIC(Mod),1)
res2[2,2] <- Mb_p
res2[2,3] <- Mb_r
res2[2,4] <- Mb_AIC


# TmaxA
Mod<-lm(SA$WF ~ SA$TmaxA)
Mc <- summary(Mod)
Mc_r <- round(as.numeric(Mc$r.squared),3)
Mc_p <- round(lmp(Mod),3)
Mc_AIC <- round(AIC(Mod),1)
res2[3,2] <- Mc_p
res2[3,3] <- Mc_r
res2[3,4] <- Mc_AIC

# synVar+TmeanA
Mod<-lm(SA$WF ~ SA$synVar + SA$TmeanA)
Mab <- summary(Mod)
Mab_r <- round(as.numeric(Mab$r.squared),3)
Mab_p <- round(lmp(Mod),3)
Mab_AIC <- round(AIC(Mod),1)
res2[4,2] <- Mab_p
res2[4,3] <- Mab_r
res2[4,4] <- Mab_AIC

# synVar+TmaxA
Mod<-lm(SA$WF ~ SA$synVar + SA$TmaxA)
Mac <- summary(Mod)
Mac_r <- round(as.numeric(Mac$r.squared),3)
Mac_p <- round(lmp(Mod),3)
Mac_AIC <- round(AIC(Mod),1)
res2[5,2] <- Mac_p
res2[5,3] <- Mac_r
res2[5,4] <- Mac_AIC

# TmeanA+TmaxA
Mod<-lm(SA$WF ~ SA$TmeanA + SA$TmaxA)
Mbc <- summary(Mod)
Mbc_r <- round(as.numeric(Mbc$r.squared),3)
Mbc_p <- round(lmp(Mod),3)
Mbc_AIC <- round(AIC(Mod),1)
res2[6,2] <- Mbc_p
res2[6,3] <- Mbc_r
res2[6,4] <- Mbc_AIC

# synVar+TmeanA+TmaxA
Mod<-lm(SA$WF ~ SA$synVar + SA$TmeanA +  SA$TmaxA)
Mabc <- summary(Mod)
Mabc_r <- round(as.numeric(Mabc$r.squared),3)
Mabc_p <- round(lmp(Mod),3)
Mabc_AIC <- round(AIC(Mod),1)
res2[7,2] <- Mabc_p
res2[7,3] <- Mabc_r
res2[7,4] <- Mabc_AIC

write.csv(res2, "Table_S3.csv")

#_____________________________________________________________________
# Part 3. Group years into tertile divisions

Div <- filter(SA, Year > 1963)%>%
  arrange(synVar) %>%
  mutate(Risk = NA)
LM <- Div[c(1:37),]
All <- Div[c(1:55),]
Div <- rbind(Div, LM)
Div <- rbind(Div, All)
Div[c(1:18),10] <- "L"
Div[c(19:37),10] <- "M"
Div[c(38:55),10] <- "H"
Div[c(56:92),10] <- "LM"
Div[c(93:147),10] <- "All"


#_____________________________________________________________________
# Part 4. Collect climate-vegetation-fire data
library(reshape2)
library(dplyr)
library(stringr)
library(outliers)
library(readr)
library(berryFunctions)
library(nls.multstart)
library(nlstools)
library(ggplot2)
library(extraDistr)

# STEP 1. Function to sort and group likelihood values by climatic tertiles
climFireTab <- function(Div, IN, tr = "L", outs = "F") {
  
  #Extract group names and build dataframe
  ifelse(outs == "T", IN <- filter(IN, outs == "TRUE"), IN <- filter(IN, outs == "FALSE"))
  study <- filter(Div, Risk == tr)%>%
    mutate(n = paste("n",Year, sep=""),)
  subs <- IN[, study$n]
  subs$X <- 1:nrow(subs) 
  studyx <- study %>% group_by(Risk) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE)
  study[length(study$n)+1,ncol(study)-1]<- "X"
  long <- melt(subs, id.vars = "X", na.rm = TRUE)
  long <- long %>% 
    mutate(gp = as.character(5*ceiling(X/5)),)%>%
    select(X, value, gp)
  tab <- long %>%
    mutate(clim = tr)
  return(tab)
}

#Load climate data and empty output table
IN <- read.csv("Pip.csv")

tab <- data.frame(X = integer(),
                  value = double(),
                  gp = double(),
                  clim = character())

for (tr in c("L", "M", "LM", "H", "All")) {
  res <- climFireTab(Div, IN=IN, tr = tr, outs = "F")
  tab <- rbind(tab, res)
}

# Summarize data for graphing and analysis
forStats <- tab %>%
  select(clim, value, X, gp) %>%
  mutate(gp = as.numeric(gp))%>%
  group_by(clim, X) %>%
  summarize_all(mean) %>%
  mutate(gp = as.character(gp))
forStats$names <- factor(forStats$gp, levels=c("5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65"))
forStats$cols <- factor(forStats$clim, levels=c("L", "M", "LM", "H", "All"))
supp.labs <- c("Low synoptic variability", "Moderate synoptic variability", "Low-mod synoptic variability", "High synoptic variability", "All")
names(supp.labs) <- c("L", "M", "LM", "H", "All")
write.csv(forStats, "forStats.csv")

means <- forStats %>%
  group_by(names) %>%
  summarize_if(is.numeric, mean)

# Fig. 2
names(SA)

Figa <- ggplot(SA , aes(y = WF, x = synVar))+
  geom_point(size = 2)+ 
  ggtitle("(a)") +
  geom_smooth(method = "lm", formula = y ~ bs(x, degree=1), color= "firebrick",
              fill = "firebrick", alpha = 0.4) +
  labs(x = "Synoptic variability (changes/season)", y = "Wildfire area (ha)" ) +
  theme_bw() +
  ylim(0,NA) +
  theme(axis.text.x  = element_text(vjust=1.5, size=11, colour = "black"),
        axis.text.y  = element_text(size=11, colour = "black"),
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))

Figb <- ggplot(SA , aes(y = synVar , x = Year))+
  geom_line(size = 0.5)+ 
  ggtitle("(b)") +
  geom_smooth(method = "lm", formula = y ~ bs(x, degree=2), color= "firebrick",
              fill = "firebrick", alpha = 0.4) +
  labs(x = "Year", y = "Synoptic variability (changes/season)" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=11, colour = "black"),
        axis.text.y  = element_text(size=11, colour = "black"),
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))

windows(8,4.5)
Figa+Figb
  plot_layout(ncol = 2)
  
  forTwo <- filter(forStats, cols == "LM" || cols == "H")
  library(ggpubr)
  windows(8,5)
  ggplot(data = forTwo, aes(names, value)) +
    ggtitle("(c)") +
    geom_boxplot(fill = "sienna4", alpha = 0.7, varwidth = TRUE) +
    geom_smooth(method="lm", formula = y ~ poly(x, 2), size = 1, color= "black",
                fill = "burlywood", alpha = 0.4, aes(group=1)) +
    stat_regline_equation(label.x = 10, label.y = 0.1) +
    ylim(0,NA) +
    stat_summary(fun=mean, geom="point", shape=1, size=7, color="black") +
    labs(y = expression(paste("Fire likelihood (ha"^"-1"*"year"^"-1"*")")), x = "Years since fire") + 
    theme_bw() +
    theme(axis.text.x  = element_text(vjust=1.5, size=11),
          axis.text.y  = element_text(size=11),
          axis.title.y = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=14),
          plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black")) +
    facet_wrap(~ cols, 
               labeller = labeller(NA,cols = supp.labs))+
    theme(
      strip.text.x = element_text(
        size = 12, face = "bold"
      )
    )

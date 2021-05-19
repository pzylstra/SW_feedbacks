#SCRIPT 4. Analysis and hypothesis-testing

library(dplyr)

set <- read.csv("forStats.csv")
setA <- filter(set, clim != "LM")
setL <- filter(set, clim == "L")
setLM <- filter(set, clim == "LM")
setM <- filter(set, clim == "M")
setH <- filter(set, clim == "H")

#___________________________________________________________________
# MANUALLY CHOOSE THE DATASET TO EXAMINE
set <- setLM
#___________________________________________________________________

out <- data.frame(matrix(ncol=0,nrow=as.numeric(max(set$X)-2)))
out$age <- seq(from = 2, to = (max(set$X)-1))
l <- as.numeric(min(set$X)+1)
u <- as.numeric(max(set$X))

for (ageDelin in 2:(u-1)) {
  a <- mean(as.numeric((filter(set, X <= ageDelin))$value))
  b <- mean(as.numeric((filter(set, X > ageDelin))$value))
  c <- mean(as.numeric(set$value))
  if (length((filter(set, X > ageDelin))$value)<2) {
    p <- NA
  } else {
  test <- t.test((as.numeric((filter(set, X <= ageDelin))$value)), 
                 (as.numeric((filter(set, X > ageDelin))$value)))
  testc <- t.test((as.numeric((filter(set, X <= ageDelin))$value)), 
                 (as.numeric(set$value)))
  testm <- t.test((as.numeric((filter(set, X > ageDelin))$value)), 
                  (as.numeric(set$value)))
  }
  p <- test$p.value
  pY <- testc$p.value
  pM <- testm$p.value
  FS <- a/b
  Y <- a/c
  M <- b/c
  out$a[ageDelin-1] <- a
  out$b[ageDelin-1] <- b
  out$FS[ageDelin-1] <- FS
  out$p[ageDelin-1] <- round(p,3)
  out$Y[ageDelin-1] <- Y
  out$pY[ageDelin-1] <- round(pY,3)
  out$M[ageDelin-1] <- M
  out$pM[ageDelin-1] <- round(pM,3)
}
#___________________________________________________________________
# EXAMINE THE OUTPUTS TO IDENTIFY AGE CLASS DELINEATIONS
View(out)
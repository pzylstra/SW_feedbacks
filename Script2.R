#SCRIPT 2. Analysis and hypothesis-testing

library(dplyr)

set <- read.csv("forStats.csv")

rep <- 1
outSum <- data.frame(matrix(ncol=13,nrow=5))
colnames(outSum) <- c('Climate', 'Y_age', 'Y_likelihood', 'y_likelihood', 
                      'R_age', 'R_likelihood', 'r_likelihood', 'D_likelihood', 'd_likelihood',
                      'M_likelihood', 'm_likelihood', 'FS', 'FSs')
for (clim in c("L", "M", "LM", "H", "All")) {
  climate <- clim
  Set <- filter(set, clim == climate)
  
  out <- data.frame(matrix(ncol=0,nrow=as.numeric(max(Set$X)-2)))
  out$age <- seq(from = 2, to = (max(Set$X)-1))
  l <- as.numeric(min(Set$X)+1)
  u <- as.numeric(max(Set$X))
  
  for (ageDelin in 2:(u-1)) {
    a <- mean(as.numeric((filter(Set, X <= ageDelin))$value))
    b <- mean(as.numeric((filter(Set, X > ageDelin))$value))
    c <- mean(as.numeric(Set$value))
    if (length((filter(Set, X > ageDelin))$value)<2) {
      p <- NA
    } else {
      test <- t.test((as.numeric((filter(Set, X <= ageDelin))$value)), 
                     (as.numeric((filter(Set, X > ageDelin))$value)))
      testc <- t.test((as.numeric((filter(Set, X <= ageDelin))$value)), 
                      (as.numeric(Set$value)))
      testm <- t.test((as.numeric((filter(Set, X > ageDelin))$value)), 
                      (as.numeric(Set$value)))
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
  #Summarise stats
  
  for (r in 1:nrow(out)) {
    out$Yage[r] = if (out$pY[r] < 0.05) {
      out$age[r]} else {0}
    out$Rage[r] = if (out$pM[r] < 0.05) {
      out$age[r]} else {100}
  }
  Ya <- max(out$Yage)
  Ra <- min(out$Rage)
  Yl <- mean(filter(Set, X <= Ya)$value)
  ly <- mean(filter(Set, X <= 5)$value)
  Rl <- mean(filter(Set, X > Ya & X < Ra)$value)
  lr <- mean(filter(Set, X > 5 & X < 51)$value)
  lD <- mean(filter(Set, X <= Ra)$value)
  ld <- mean(filter(Set, X <= 50)$value)
  Ml <- mean(filter(Set, X > Ra)$value)
  l_m <- mean(filter(Set, X > 50)$value)
  FSa <- round(out[(Ra-1),4],1) 
  FSs <- round(out[50,4],1) 
  
  outSum$Climate[rep] <- clim
  outSum$Y_age[rep] <- Ya
  outSum$Y_likelihood[rep] <- signif(Yl,2)
  outSum$y_likelihood[rep] <- signif(ly,2)
  outSum$R_age[rep] <- Ra
  outSum$R_likelihood[rep] <- signif(Rl,2)
  outSum$r_likelihood[rep] <- signif(lr,2)
  outSum$D_likelihood[rep] <- signif(lD,2)
  outSum$d_likelihood[rep] <- signif(ld,2)
  outSum$M_likelihood[rep] <- signif(Ml,2)
  outSum$m_likelihood[rep] <- signif(l_m,2)
  outSum$FS[rep] <- FSa
  outSum$FSs[rep] <- FSs
  rep <- rep+1
}

View(outSum)
write.csv(outSum, "FS.csv")

# SCRIPT 5. Attribute treatment vs climatic effects

library(dplyr)

set <- read.csv("forStats.csv")
year <- read.csv("treat_test.csv")

burnLM <- data.frame()
ageLM <- data.frame()
burnH <- data.frame()
ageH <- data.frame()
LM <- data.frame()
H <- data.frame()

for (n in 1:length(year$Year)) {
 testset <- filter(set, clim == year$synVar[n])
 climate <- year$synVar[n]
 bL <- year$Bl[n]
 bU <- year$Bu[n]
 aL <- year$Al[n]
 aU <- year$Au[n]
 burn <- as.numeric((filter(testset, X > bL & X <= bU))$value)
 datB <- as.data.frame(burn)
 age <- as.numeric((filter(testset, X > aL & X <= aU))$value)
 datA <- as.data.frame(age)
 if (climate == "LM") {
         burnLM <- rbind(burnLM, datB)
         ageLM <- rbind(ageLM, datA)
 } else {
         burnH <- rbind(burnH, datB)
         ageH <- rbind(ageH, datA)
 }
}

#Test Low-mod difference
bLMmean <- round(mean(burnLM$burn),4)
aLMmean <- round(mean(ageLM$age),4)
diffLM <- round(bLMmean/aLMmean,4)
testLM <- t.test(burnLM$burn, ageLM$age)
pLM <- round(testLM$p.value,4)

#Test High difference
bHmean <- round(mean(burnH$burn),4)
aHmean <- round(mean(ageH$age),4)
diffH <- round(bHmean/aHmean,4)
testH <- t.test(burnH$burn, ageH$age)
pH <- round(testH$p.value,4)

#Test all
bAll <- rbind(burnLM, burnH)
aAll <- rbind(ageLM, ageH)
bmean <- round(mean(bAll$burn),4)
amean <- round(mean(aAll$age),4)
diff <- round(bmean/amean,4)
test <- t.test(bAll$burn, aAll$age)
p <- round(test$p.value,4)

# Tabulate data
comp <- data.frame(synVar = c("LM", "H", "All"), 
                   Burn = c(bLMmean, bHmean, bmean),
                   NoBurn = c(aLMmean, aHmean,amean),
                   Difference = c(diffLM, diffH, diff),
                   p = c(pLM, pH, p))
write.csv(comp, "Table2.csv")

# Arrange and plot results
burnLM <- burnLM %>%
        mutate(Likelihood = burn,
               Treatment = "Burn",
               Climate = "LM") %>%
        select(Treatment, Climate, Likelihood)
burnH <- burnH %>%
        mutate(Likelihood = burn,
               Treatment = "Burn",
               Climate = "H") %>%
        select(Treatment, Climate, Likelihood)
ageLM <- ageLM %>%
        mutate(Likelihood = age,
               Treatment = "No-burn",
               Climate = "LM") %>%
        select(Treatment,Climate, Likelihood)
ageH <- ageH %>%
        mutate(Likelihood = age,
               Treatment = "No-burn",
               Climate = "H") %>%
        select(Treatment, Climate, Likelihood)
all <- rbind(burnLM, burnH, ageLM, ageH)
all$cols <- factor(all$Climate, levels=c("LM", "H"))
supp.labs <- c("Low-mod synoptic variability", "High synoptic variability")
names(supp.labs) <- c("LM", "H")

# Fig. 1b

library(ggplot2)
windows(8,4.5)
ggplot(data = all, aes(Treatment, Likelihood)) +
        ggtitle("b.") +
        geom_violin(fill = "darkolivegreen") +
        geom_boxplot(width = 0.1, notch=TRUE)+
        stat_summary(fun=mean, geom="point", shape=1, size=7, color="black") +
        labs(y = expression(paste("Fire likelihood (ha"^"-1"*"year"^"-1"*")")), x = "Treatment") + 
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

# Climate effect

cBurn <- filter(all, Treatment == "Burn")%>%
        group_by(Climate) %>%
        summarise_all(mean) %>%
        mutate(Burn = Likelihood) %>%
        select(Climate, Burn)
LM <- filter(all, Climate == "LM" & Treatment == "Burn")
H <- filter(all, Climate == "H" & Treatment == "Burn")
dBurn <- round(cBurn$Burn[1] / cBurn$Burn[2],2)
test <- t.test(LM$Likelihood, H$Likelihood)
pBurn <- round(test$p.value,4)

cNB <- filter(all, Treatment == "No-burn")%>%
        group_by(Climate) %>%
        summarise_all(mean) %>%
        mutate(NB = Likelihood) %>%
        select(Climate, NB)
LM <- filter(all, Climate == "LM" & Treatment == "No-burn")
H <- filter(all, Climate == "H" & Treatment == "No-burn")
dNB <- round(cNB$NB[1] / cNB$NB[2],2)
test <- t.test(LM$Likelihood, H$Likelihood)
pNB <- round(test$p.value,4)

climEff <- left_join(cBurn, cNB, by = "Climate")
analysis <- data.frame(Climate = c("Difference", "p"), Burn = c(dBurn, pBurn), NB = c(dNB, pNB))


# Fig. S8

windows(4.5,4.5)
climBurn <- filter(all, Treatment == "Burn")
ggplot(data = climBurn, aes(Climate, Likelihood)) +
        geom_violin(fill = "darkolivegreen") +
        geom_boxplot(width = 0.1, notch=TRUE)+
        stat_summary(fun=mean, geom="point", shape=1, size=7, color="black") +
        labs(y = expression(paste("Fire likelihood (ha"^"-1"*"year"^"-1"*")")), x = "Synoptic variability class") + 
        theme_bw() +
        theme(axis.text.x  = element_text(vjust=1.5, size=11),
              axis.text.y  = element_text(size=11),
              axis.title.y = element_text(face="bold", size=14),
              axis.title.x = element_text(face="bold", size=14))

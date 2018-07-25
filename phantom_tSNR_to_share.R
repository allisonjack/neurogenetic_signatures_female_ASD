setwd("/Users/YourUsername/WhereverYouPutTheFile")

# -------(Install and) Load Packages-----------
# Uncomment if necessary to install packages

#install.packages(c("nlme", "multcomp"))

library(nlme)
library(multcomp)

# ------ANOVA-------------

df <- read.csv("phantom_tSNR_data_to_share.csv", header = TRUE)

summary(A <- lme(tSNR.mean ~ Site, data = df, random = ~1|Name))
anova(A)
summary(glht(A, linfct = mcp(Site = "Tukey")))

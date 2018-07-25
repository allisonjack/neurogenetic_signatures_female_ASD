#don't forget to change the path!
#setwd("/Users/YourUsername/WhereverYouPutTheFile")

# --------Description of Variables--------

#anon.ID: anonymized ID numbers
#rms.abs.max: Maximum absolute root mean squared head motion (in mm) across the run
#rms.rel.max: Maximum relative root mean squared head motion (in mm) across the run
#Group:
  #ASD: Autism Spectrum Disorder
  #TD: Typically developing
#scan.age.mos: Age (in months) at which participant was scanned
#Sex:
  #M: Male
  #F: Female
#SexGroup:
  #TDf: TD female
  #ASDf: ASD female
  #TDm: TD male
  #ASDm: ASD male
#DAS.GCA.SS: Differential Ability Scales 2nd edition General Conceptual Ability Standard Score

# -------(Install and) Load Packages & Data-----------

# Uncomment if necessary to install packages
#install.packages(c("MatchIt", "optmatch", "tidyverse", "stringr"))

library(MatchIt)
library(optmatch)
library(tidyverse)
#library(stringr)

#read in main file
df <- read_csv("dataframe_for_matching_to_share.csv")

# --------Visualize head motion distribution and remove some outliers------------

ggplot(df, aes(x = rms.abs.max, y = rms.rel.max, color = SexGroup)) +
  geom_point() +
  geom_text(aes(y = rms.rel.max + .15,
                label = anon.ID),
            size = 2) +
  geom_hline(yintercept = 3) + 
  geom_vline(xintercept = 3)
# only n = 3 fall above 3mm maximum in both absolute & relative motion
# probably best to remove as outliers

motion_outliers <- c(209, 208, 210)

df <- df %>% 
  filter(df$anon.ID %in% motion_outliers == FALSE)

# ---Test whether the 4 groups differ on key metrics using ANOVA and follow up pairwise comparisons-----------

#do sexgroups differ on head motion, IQ, or age?
(aov.models <- tibble(x = c('rms.abs.max', 'rms.rel.max', 'DAS.GCA.SS', 'scan.age.mos')) %>%
    mutate(
      aovmod  = map(x, ~ aov(as.formula(paste(., '~ SexGroup')), data = df)),
      tuk     = map(aovmod, ~ TukeyHSD(.)),
      aovsumm = map(aovmod, ~ summary(.)),
      fval     = map_dbl(aovsumm, ~ unlist(.)[["F value1"]]),
      pval     = map_dbl(aovsumm, ~ unlist(.)[["Pr(>F)1"]]),
      df1      = map_dbl(aovsumm, ~ unlist(.)[["Df1"]]),
      df2      = map_dbl(aovsumm, ~ unlist(.)[["Df2"]])
    )
)

#view the pairwise comparisons for metrics with significant p value
aov.models.sig <- aov.models %>%
  #scan age doesn't differ, so no need to see those pairwise comparisons
  filter(x != "scan.age.mos") %>%
  #just look at Tukey HSD
  select(x, tuk)

#create custom function to only return pairwise comparisons where p < .10
get.tuk <- function(y){
  tuk.df <- as.data.frame(aov.models.sig$tuk[aov.models.sig$x == y][[1]][[1]])
  return(subset(tuk.df, tuk.df[4] < 0.10)[,c(1,4)])
}

get.tuk("rms.abs.max")
get.tuk("rms.rel.max")
get.tuk("DAS.GCA.SS")


#---------Examine differences specifically between TDf & ASDf and TDm & ASDm using t-tests----------

by_sex <- df %>%
  group_by(Sex) %>%
  nest()

t_fx <- function(x, df){
  t.test(as.formula(paste(x, '~ Group')), data = df)
}

t.df <- by_sex %>%
  mutate(
    t.abs   = map(data, t_fx, x = 'rms.abs.max'),
    #ABS.t   = map_dbl(t.abs, "statistic"),
    #ABS.df  = map_dbl(t.abs, 'parameter'),
    ABS.p   = map_dbl(t.abs, 'p.value'),
    
    t.rel   = map(data, t_fx, x = 'rms.rel.max'),
    #REL.t   = map_dbl(t.rel, 'statistic'),
    #REL.df  = map_dbl(t.rel, 'parameter'),
    REL.p   = map_dbl(t.rel, 'p.value'),
    
    t.IQ    = map(data, t_fx, x = 'DAS.GCA.SS'),
    #IQ.t    = map_dbl(t.IQ, 'statistic'),
    #IQ.df   = map_dbl(t.IQ, 'parameter'),
    IQ.p    = map_dbl(t.IQ, 'p.value')) %>%
  select(-data, -t.abs, -t.rel, -t.IQ)

t.df
# for girls, IQ differs but not head motion
# for boys, IQ and motion both differ.

# ----------------Add binarized grouping variable for matching------------

df <- df %>%
  mutate(dx.bin = ifelse(Group == "ASD", 1, 0))

# -----------Find TD girl matches for ASD girls--------------
df.f <- as.data.frame(subset(df, Sex == "F"))

#having previously run the matching algorithm, I remove anon.ID == 24 because
#distance from the best match (0.7253730) was too high/couldn't find a good match
df.f <- subset(df.f, anon.ID != 24) 

rownames(df.f) <- df.f$anon.ID
df.f <- subset(df.f, select = c(anon.ID, Group, dx.bin, DAS.GCA.SS, rms.rel.max, rms.abs.max))
m.outF <- matchit(dx.bin ~ DAS.GCA.SS,
                  data = df.f,
                  method = "optimal")
match.F <- match.data(m.outF, group = "treat")
match.F$match.ID <- as.vector(m.outF$match.matrix)

ASDf.matched <- match.F$anon.ID
TDf.matched <-match.F$match.ID
f.matched <- as.factor(c(ASDf.matched, TDf.matched))
f.matched <- data.frame(f.matched)
f.matched <- rename(f.matched, anon.ID = f.matched)

#Do ASDf and TDf differ on key metrics after matching?
df.f.matched <- merge(df.f, f.matched, by = "anon.ID")
t.test(DAS.GCA.SS ~ Group, data = df.f.matched)
t.test(rms.abs.max ~ Group, data = df.f.matched)
t.test(rms.rel.max ~ Group, data = df.f.matched)

# -----------Find TD boy matches for ASD boys--------------
df.m <- as.data.frame(subset(df, Sex == "M"))

#having previously run the matching algorithm, I remove anon.ID == 160 because
#distance from the best match (0.8001177) was too high/couldn't find a good match
df.m <- subset(df.m, anon.ID != 160) 

rownames(df.m) <- df.m$anon.ID

df.m <- subset(df.m, select = c(anon.ID, Group, dx.bin, DAS.GCA.SS, rms.rel.max, rms.abs.max))
m.outM <- matchit(dx.bin ~ DAS.GCA.SS + rms.abs.max,
                  data = df.m,
                  method = "optimal")
match.M <- match.data(m.outM, group = "treat")
match.M$match.ID <- as.vector(m.outM$match.matrix)

ASDm.matched <- match.M$anon.ID
TDm.matched <-match.M$match.ID
m.matched <- as.factor(c(ASDm.matched, TDm.matched))
m.matched <- data.frame(m.matched)
m.matched <- rename(m.matched, anon.ID = m.matched)

df.m.matched <- merge(df.m, m.matched, by = "anon.ID")
t.test(DAS.GCA.SS ~ Group, data = df.m.matched)
t.test(rms.rel.max ~ Group, data = df.m.matched) 
t.test(rms.abs.max ~ Group, data = df.m.matched) 

# -----------Make a new dataframe of the matched sample ------
identical(names(df.m.matched), names(df.f.matched))
match.IDs <- rbind(df.f.matched, df.m.matched)
match.IDs <- subset(match.IDs, select = anon.ID)
df.match <- merge(df, match.IDs)

# ---------------Do groups differ on key metrics after matching?---------

(aov.models <- tibble(x = c('rms.abs.max', 'rms.rel.max', 'DAS.GCA.SS')) %>%
    mutate(
      aovmod  = map(x, ~ aov(as.formula(paste(., '~ SexGroup')), data = df.match)),
      tuk     = map(aovmod, ~ TukeyHSD(.)),
      aovsumm = map(aovmod, ~ summary(.)),
      fval     = map_dbl(aovsumm, ~ unlist(.)[["F value1"]]),
      pval     = map_dbl(aovsumm, ~ unlist(.)[["Pr(>F)1"]]),
      df1      = map_dbl(aovsumm, ~ unlist(.)[["Df1"]]),
      df2      = map_dbl(aovsumm, ~ unlist(.)[["Df2"]])
    )
)
aov.models.sig <- aov.models %>%
  select(x, tuk)

#IQ no longer differs across groups

get.tuk("rms.abs.max")
get.tuk("rms.rel.max")

#head motion metrics do have a significant p-value on the omnibus test but
#pairwise comparisons reveal that the only pairs with significant differences
#are not contrasts of interest
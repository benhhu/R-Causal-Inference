#Causal Inference
library(sas7bdat)
library(car)

setwd("//file.phhp.ufl.edu/home/huihu/files/RWD/R-Causal-Inference")
nhefs <- read.sas7bdat("nhefs_book.sas7bdat")
nhefs$cens <- as.numeric(is.na(nhefs$wt82))
nhefs$older <- as.numeric(nhefs$age > 50 & !is.na(nhefs$age))


nhefs$education.code <- recode(nhefs$school, " 0:8 = 1; 9:11 = 2; 12 = 3; 13:15 = 4; 16:hi = 5; NA = 6 ")
nhefs$education <- recode(nhefs$education.code, " 1 = '1. 8th grade or less'; 2 = '2. HS dropout'; 3 = '3. HS'; 4 = '4. College dropout'; 5 = '5. College or more'; 6 = 'Unknown' ")

#IP weighting
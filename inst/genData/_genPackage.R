#------------ generating package tw.DEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")

.tmp.f <- function(){
	library(twMisc)
	source("R/logitnorm.R")
	twUtestF()
}

#R CMD build logitnorm
#R CMD check --as-cran logitnorm_XXX | tee tmp.txt | grep -v OK$
#R CMD check --no-vignettes --no-latex --no-codoc logitnorm
#R CMD INSTALL --html logitnorm

.tmp.f <- function(){
	library(sos)
	fres1 <- findFn("logitnormal") 
} 

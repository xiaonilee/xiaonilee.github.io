pkgname <- "npar"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('npar')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("oneway")
### * oneway

flush(stderr()); flush(stdout())

### Name: oneway
### Title: Nonparametric group comparisons
### Aliases: oneway

### ** Examples

results <- oneway(hlef ~ region, life)
summary(results)
plot(results, col="lightblue", main="Multiple Comparisons",
     xlab="US Region", ylab="Healthy Life Expectancy at Age 65")



cleanEx()
nameEx("plot.oneway")
### * plot.oneway

flush(stderr()); flush(stdout())

### Name: plot.oneway
### Title: Plot nonparametric group comparisons
### Aliases: plot.oneway

### ** Examples

results <- oneway(hlef ~ region, life)
plot(results, col="lightblue", main="Multiple Comparisons",
     xlab="US Region", ylab="Healthy Life Expectancy at Age 65")



cleanEx()
nameEx("print.oneway")
### * print.oneway

flush(stderr()); flush(stdout())

### Name: print.oneway
### Title: Print multiple comparisons
### Aliases: print.oneway

### ** Examples

results <- oneway(hlef ~ region, life)
print(results)



cleanEx()
nameEx("summary.oneway")
### * summary.oneway

flush(stderr()); flush(stdout())

### Name: summary.oneway
### Title: Summarize oneway nonparametric analyses
### Aliases: summary.oneway

### ** Examples

results <- oneway(hlef ~ region, life)
summary(results)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')


library("STACMR")

## load data from Exp. 1 of Dunn, Newell, & Kalish (2012)
data(delay)
head(delay)

y <- gen2list(data=delay)
str(y, 2)

delaystats <- staSTATS(delay)
str(delaystats)

delaystats[[1]]$means

delaystats[[2]]$means

staPLOT(data = delay, 
        groups = list(c(1:4), c(5:8)), 
        grouplabels = list("No delay", "Delay"), 
        axislabels = list("RB","II"))


### Fits Partial Order Model to Data
# create partial order
E <- list(c(1:4),c(5:8),c(5,1),c(6,2),c(7,3),c(8,4))
out2 <- staMR(data=delay, partial=E)

cbind(out2$x[[1]],out2$x[[2]]) # simplify presentation
out2$fval
out2$shrinkage

### Conduct CMR State-Trace Analysis
out1 <- staCMR(data=delay, partial=E)

cbind(out1$x[[1]],out1$x[[2]])
out1$fval
out1$shrinkage

staPLOT(data = delay, 
        groups = list(c(1:4), c(5:8)), 
        grouplabels = list("No delay", "Delay"), 
        axislabels = list("RB","II"),
        pred=out1$x)

### Test fit of Partial Order Model

out3 <- staMRFIT(delay, partial=E, nsample=10000)
out3$p
out3$datafit

### p-Value of Difference in Fit of Conjoint Monotonic and Fit of Partial Order Model
out4 <- staCMRFIT(delay, partial=E, nsample=10000)
out4$p
out4$datafit


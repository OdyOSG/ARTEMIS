library(reticulate)
library(oncoRegimens)

g <- 0.5
Tfac <- 0.2
local_align = 1
verbose = 1
mem = as.integer(10)
removeOverlap = 1

#QCAA
regimen <- encode("0Q.2C.2A.1A")
#CA*QCAA*AQQ*QCAA*AQQ*QCAA*
drugRecord <- encode("0C.0A.9Q.2C.2A.1A.0A.0Q.0Q.0Q.2C.2A.1A.0A.0Q.0Q.0Q.3C.2A.1A")

output <- align(regimen,drugRecord,g,Tfac,NA,local_align,verbose,mem,removeOverlap)



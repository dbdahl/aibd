# Close all R sessions.
# In a new R session, run aibd2:::installJARs()
# In RStudio menu, select Build -> Install and Restart

# What about exaluating multiple Zs? ----------------------------------------------------
library(aibd2)
ibp1 <- ibp(1,4)
samples <- sampleFeatureAllocation(5e3, ibp1) # This is a list

# R
system.time(y <- prFeatureAllocation(samples, distribution = ibp1))

# Scala
system.time(y <- prFeatureAllocation(samples,distribution = ibp1, implementation='scala'))

# Scala a second time
system.time(y <- prFeatureAllocation(samples,distribution = ibp1, implementation='scala'))

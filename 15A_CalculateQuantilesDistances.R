variants <- 97690
genome.length <- 103014376

#Draw random positions with uniform distribution
positions <- runif(variants, 1, genome.length)

#Sort and round
positions.sorted <- sort(positions)
positions.rounded <- round(positions.sorted)

#Calculate distances to subsequent positions
distances <- positions.rounded[2:97690] - positions.rounded[1:97689]

#Plot distances
hist(distances, breaks=200)

#Get the percentiles for the distribution
quantile(distances, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))

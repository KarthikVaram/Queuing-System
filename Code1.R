simulate_mm1 <- function(lambda, mu, n_customers){
  
  # Generate interarrival and service times
  interarrival <- rexp(n_customers, rate = lambda)
  service <- rexp(n_customers, rate = mu)
  
  # Arrival times
  arrival <- cumsum(interarrival)
  
  # Initialize vectors
  start <- numeric(n_customers)
  departure <- numeric(n_customers)
  waiting <- numeric(n_customers)
  
  # First customer
  start[1] <- arrival[1]
  departure[1] <- start[1] + service[1]
  waiting[1] <- 0
  
  # Remaining customers
  for(i in 2:n_customers){
    
    start[i] <- max(arrival[i], departure[i-1])
    departure[i] <- start[i] + service[i]
    waiting[i] <- start[i] - arrival[i]
  }
  
  # System time = waiting + service
  system_time <- waiting + service
  
  return(list(
    arrival = arrival,
    service = service,
    start = start,
    departure = departure,
    waiting = waiting,
    system_time = system_time
  ))
}


result <- simulate_mm1(lambda = 2, mu = 3, n_customers = 10000)

mean(result$waiting)
mean(result$system_time)
lambda <- 2
mu <- 3

rho <- lambda / mu

theoretical_wait <- rho / (mu - lambda)
simulated_wait <- mean(result$waiting)

theoretical_wait
simulated_wait

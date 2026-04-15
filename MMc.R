simulate_mmc <- function(lambda, mu, c, max_events = 100000){
  
  t <- 0
  n <- 0
  
  times <- numeric(max_events)
  states <- numeric(max_events)
  
  for(i in 1:max_events){
    
    # rates
    arrival_rate <- lambda
    service_rate <- min(n, c) * mu
    
    total_rate <- arrival_rate + service_rate
    
    # time to next event
    dt <- rexp(1, total_rate)
    t <- t + dt
    
    times[i] <- t
    states[i] <- n
    
    # event type
    if(runif(1) < arrival_rate / total_rate){
      n <- n + 1
    } else {
      n <- max(0, n - 1)
    }
  }
  
  list(times = times, states = states)
}

estimate_pj_mmc <- function(times, states, max_state = 15){
  
  durations <- diff(c(0, times))
  total_time <- sum(durations)
  
  p_hat <- sapply(0:max_state, function(j)
    sum(durations[states == j]) / total_time)
  
  return(p_hat)
}

estimate_metrics_mmc <- function(times, states, c){
  
  durations <- diff(c(0, times))
  total_time <- sum(durations)
  
  # L
  L_hat <- sum(states * durations) / total_time
  
  # Lq
  queue <- pmax(states - c, 0)
  Lq_hat <- sum(queue * durations) / total_time
  
  # P_wait
  P_wait_hat <- sum(durations[states >= c]) / total_time
  
  list(L = L_hat, Lq = Lq_hat, P_wait = P_wait_hat)
}

theoretical_mmc <- function(lambda, mu, c){
  
  rho <- lambda / (c * mu)
  
  # compute p0
  sum1 <- sum(sapply(0:(c-1), function(j) (lambda/mu)^j / factorial(j)))
  sum2 <- (lambda/mu)^c / factorial(c) * (1 / (1 - rho))
  
  p0 <- 1 / (sum1 + sum2)
  
  # P_wait
  P_wait <- ((lambda/mu)^c / factorial(c)) * (1 / (1 - rho)) * p0
  
  # Lq
  Lq <- (P_wait * rho) / (1 - rho)
  
  # L
  L <- Lq + lambda / mu
  
  list(L = L, Lq = Lq, P_wait = P_wait)
}

theoretical_pj_mmc <- function(lambda, mu, c, max_state = 15){
  
  rho <- lambda / (c * mu)
  
  # compute p0
  sum1 <- sum(sapply(0:(c-1), function(j) (lambda/mu)^j / factorial(j)))
  sum2 <- (lambda/mu)^c / factorial(c) * (1 / (1 - rho))
  
  p0 <- 1 / (sum1 + sum2)
  
  pj <- numeric(max_state + 1)
  
  for(j in 0:max_state){
    if(j < c){
      pj[j+1] <- ((lambda/mu)^j / factorial(j)) * p0
    } else {
      pj[j+1] <- ((lambda/mu)^c / factorial(c)) *
        (lambda/(c*mu))^(j-c) * p0
    }
  }
  
  return(pj)
}

set.seed(42)

lambda <- 4
mu <- 2
c <- 3
max_state <- 15

# -------------------------------
# Run simulation
# -------------------------------
sim <- simulate_mmc(lambda, mu, c, 100000)

# -------------------------------
# Estimate metrics
# -------------------------------
est <- estimate_metrics_mmc(sim$times, sim$states, c)
theory <- theoretical_mmc(lambda, mu, c)

print(est)
print(theory)

# -------------------------------
# Step 3: Verify p_j
# -------------------------------
p_hat <- estimate_pj_mmc(sim$times, sim$states, max_state)
p_true <- theoretical_pj_mmc(lambda, mu, c, max_state)

pj_results <- data.frame(
  j = 0:max_state,
  simulated = round(p_hat, 4),
  theoretical = round(p_true, 4)
)

print(pj_results)

cat("Max error:", max(abs(p_hat - p_true)), "\n")

# -------------------------------
# Step 4: Plot comparison
# -------------------------------
plot(0:max_state, p_hat, type="b", pch=16,
     xlab="State j", ylab="Probability",
     main="M/M/c: Stationary Distribution")

lines(0:max_state, p_true, col="red", type="b", pch=17)

legend("topright",
       legend=c("Simulated", "Theoretical"),
       col=c("black", "red"),
       pch=c(16,17))


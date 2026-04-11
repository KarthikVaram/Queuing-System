# -------------------------------
# M/M/1 Queue Simulation
# Limiting Probabilities Verification
# -------------------------------

# -------------------------------
# 1. Simulation (CTMC approach)
# -------------------------------
simulate_mm1_states <- function(lambda, mu, max_events = 100000){
  
  t <- 0
  n <- 0  # current state
  
  times <- numeric(max_events)
  states <- numeric(max_events)
  
  for(i in 1:max_events){
    
    # Total rate depends on state
    rate <- if(n == 0) lambda else (lambda + mu)
    
    # Time to next event
    dt <- rexp(1, rate)
    t <- t + dt
    
    # Record state BEFORE jump
    times[i] <- t
    states[i] <- n
    
    # Determine event type
    if(n == 0){
      n <- n + 1
    } else {
      if(runif(1) < lambda / (lambda + mu)){
        n <- n + 1   # arrival
      } else {
        n <- n - 1   # departure
      }
    }
  }
  
  return(list(times = times, states = states))
}

# -------------------------------
# 2. Estimate limiting probabilities
# -------------------------------
estimate_pj <- function(times, states, max_state = 10){
  
  durations <- diff(c(0, times))
  total_time <- sum(durations)
  
  p_hat <- numeric(max_state + 1)
  
  for(j in 0:max_state){
    p_hat[j+1] <- sum(durations[states == j]) / total_time
  }
  
  return(p_hat)
}

# -------------------------------
# 3. Theoretical probabilities
# -------------------------------
theoretical_pj <- function(lambda, mu, max_state = 10){
  rho <- lambda / mu
  p <- sapply(0:max_state, function(j) (1 - rho) * rho^j)
  return(p)
}

# -------------------------------
# 4. Run Simulation
# -------------------------------
set.seed(42)

lambda <- 2
mu <- 3
max_state <- 10

sim <- simulate_mm1_states(lambda, mu, max_events = 100000)

p_hat <- estimate_pj(sim$times, sim$states, max_state)
p_true <- theoretical_pj(lambda, mu, max_state)

# -------------------------------
# 5. Compare Results
# -------------------------------
results <- data.frame(
  j = 0:max_state,
  simulated = round(p_hat, 4),
  theoretical = round(p_true, 4)
)

print(results)

# -------------------------------
# 6. Plot
# -------------------------------
plot(0:max_state, p_hat, type="b", pch=16,
     xlab="State j", ylab="Probability",
     main="Limiting Distribution: Simulation vs Theory")

lines(0:max_state, p_true, type="b", col="red", pch=17)

legend("topright",
       legend=c("Simulated", "Theoretical"),
       col=c("black", "red"),
       pch=c(16,17))

# -------------------------------
# 7. Error Check 
# -------------------------------
cat("\nMax absolute error:",
    max(abs(p_hat - p_true)), "\n")
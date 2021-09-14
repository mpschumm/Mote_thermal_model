# Runtime in years
runtime <- 320
# Longevity of the fish in years
longevity = 25
# Setting up 3d list that tracks cohort variables through time
timepoints <<- vector(mode = "list", length = runtime)
timepoints[[1]] <- matrix(data=0, nrow=longevity, ncol=3)
timepoints[[1]][1,1] <- 1000
timepoints[[1]][1,2] <- 0.1
timepoints[[1]][1,3] <- 0
# Creating temperature trajectories
temps <- c(rep(10,200), rep(10,120))
temps <- c(rep(10,200), seq(10,13, by=(3/(120-1))))
temps <- c(rep(10,200), seq(10,16, by=(6/(120-1))))
# What type of recruitment to use?
recruitment_toggle <- BH
recruitment_toggle <- Ricker

# Run the model
mapply( function(x) {
  timepoints[[x]] <<- model_run(timepoints[[x-1]], x)
}, 2:runtime)

# Function that runs the model for a year
model_run <- function(stages, x) {
  temp <- temps[x]
  lengths<-stages[,2]
  popn <- stages[,1]
  lengths <- lengths + growth(lengths)*TPC(temp)
  # Apply mortality from non-fishing causes
  popn <- mortality(lengths, popn, m_p_max, m, 0)*TPC(temp)
  realizable_catch <- sum(popn)
  # Apply fishing mortality
  popn <- mortality(lengths, popn, 0, 0, F_M_init)
  realizable_catch <- realizable_catch - sum(popn)
  FSSB <- sum(fecundity(popn, lengths))
  popn[longevity] <- popn[longevity] + popn[(longevity-1)]
  popn[2:(longevity-1)] <- popn[1:(longevity-2)] 
  if (recruitment_toggle==BH) {
  popn[1] <- 1000*(BH(FSSB/(1000^2)))/2
  } else {
  popn[1] <- 1000*(Ricker(FSSB/(1000^2)))/2 
  }
  lengths[2:longevity] <- lengths[1:(longevity-1)]
  lengths[1] <- 0.1
  return(as.matrix(cbind(popn, lengths, realizable_catch)))
}

# Function for thermal performance
omega_max <- 1
T_H <- 17.5
T_max <- 11
T_L <- 2.5
TPC <- function(temp) {
  return((omega_max*(temp-T_H)*(temp-T_L)^2)/((T_max-T_L)*((T_max-T_L)*(temp-T_max)-(T_max-T_H)*(T_max+T_L-2*temp))))
}

# Length-dependent growth function
growth <- function(l) {
  return(ifelse(l>0, 0.27*(1.19-l), 0))
}

# Predation mortality (for higher mortality on smaller individuals): roughly 0.5 year-based instantaneous mortality at ~1 yrs., Koster et al 2003
# exp(-0.5) = 0.6065307 = 60.65307% remaining after 1 yr.
# 0.6065307^(1/365) = 0.9986
m_p_max <- log(0.9986)*(-1)
# Set fishing mortality in absence of restriction due to temperature 
# e.g., 0.2 for F=0.2
F_M_init <- exp(-(0.2))^(1/365)
F_M_init <- log(0.997264)*(-1)
# Base natural mortality
m<-log(0.999452)*(-1)
# Function for applying mortality
mortality <- function(l, no, m_p_max, m, F_M) {
  if (!(is.numeric(F_M))) {
    return(0)
  } else {
  return(no*(exp(-m -(m_p_max-m)*exp(-8*l)) -F_M*selective_mortality(l))^365)
}
}

# Selective harvest function
mincatchsize = 0.5
selective_mortality <- function(length) {
  1 / (1 + exp(-(length - mincatchsize)*100))
}

# Fecundity function
fecundity <- function(no, l) {
  return(no*5787*l^3)
}

# Density dependent recruitment, Ricker
Alpha <- 5.6
beta <- -3e-06
Ricker <- function(biomass) {
  return(Alpha*biomass*exp(beta*biomass))
}

# Density dependent recruitment, Beverton Holt
a <- 8
b <- 1e-05
BH <- function(biomass) {
  return(a*biomass/(1+b*biomass))
}

# Functions for summarizing and plotting model output
recruit_sum<-function(point) {
  sum(point[1:longevity,1])
}

egg_sum<-function(point) {
  return((unricker(point[1,1]*2/1000))*(1000^2))
}

unricker <- function(x) {
  return(lambertW0((x*beta)/Alpha)/beta)
}

harvest_sum<-function(point) {
  sum(point[1:longevity,1]*(1-exp(-F_M_init*selective_mortality(point[1:longevity,2]))^365)*point[1:longevity,2])
}

harvest_adaptive_sum<-function(point) {
  sum(point[1:longevity,1]*(1-exp(-F_M_init*TPC(point[1:longevity,3])*selective_mortality(point[1:longevity,2]))^365)*point[1:longevity,2])
}

sum_biomass<-function(point) {
  sum(point[1:longevity,1]*point[1:longevity,2])
}

sum_biomass<-function(point) {
  sum(point[1:longevity,1]*point[1:longevity,2])
}

plot(unlist(lapply(timepoints, recruit_sum)))

plot(unlist(lapply(timepoints, harvest_sum)))

plot(unlist(lapply(timepoints, harvest_adaptive_sum)))

plot(unlist(lapply(timepoints, sum_biomass)))


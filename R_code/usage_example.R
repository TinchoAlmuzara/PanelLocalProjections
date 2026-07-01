### USAGE EXAMPLE #########################################################
# This script illustrates the usage of panel_LP, the function to implement
# estimation and inference as recommended in "Micro Responses to Macro Shocks"
# by M. Almuzara and V. Sancibrian.
#
# It mirrors usage_example.m (the Matlab version). Because R and Matlab use
# different random number generators, the simulated data -- and hence the
# numbers in the plot -- will not match the Matlab script exactly, but the
# data-generating process and the call to panel_LP are identical.
#
# Version: 2026 June 30 - R 4.4.2
##########################################################################

# Clear memory 
rm(list = ls())
source("panel_LP.R")
set.seed(1918)

# Simulate a simple dataset
T_dim   <- 30
N       <- 1000
macro   <- rnorm(T_dim)                 
X       <- rep(macro, each = N)         
firm    <- 1 + rnorm(N)                 
s       <- rep(firm, times = T_dim)
b_true  <- 1
y       <- (s * X) * b_true +
           as.vector(outer(0.5*s[1:N] + rnorm(N), rnorm(T_dim))) +  
           rnorm(N * T_dim)                            
t_index <- rep(1:T_dim, each = N)
i_index <- rep(1:N, times = T_dim)
H       <- 5

# Call panel local projections function
# In R the inputs are passed as named arguments. Optional arguments and their defaults:
#   s            : heterogeneity characteristic(s); omit for a homogeneous shock
#   W            : controls; omit (NULL) if no controls are to be used
#   FE           : unit and time fixed effects; omit (NULL) if no fixed effects are to be used
#   p_max        : number of lags of regressand and shock to be added as controls
#   small_sample : Imbens-Kolesar (2016, REStat) small-sample refinement; defaults to FALSE
#   cumulative   : report cumulative impulse responses; defaults to FALSE
LP_out <- panel_LP(y            = y,
                   X            = X,
                   s            = s,
                   i_index      = i_index,
                   t_index      = t_index,
                   W            = NULL,                  
                   FE           = cbind(i_index, t_index), 
                   H            = H,                     # impulse response horizon
                   p_max        = 0,                     
                   small_sample = TRUE,                 
                   cumulative   = TRUE)                  

# Plot estimates with 90% confidence bands
horizons <- 0:H
blue     <- rgb(0, 85, 164, maxColorValue = 255)
plot_df  <- data.frame(h        = horizons,
                       estimate = LP_out$estimate[, 1],
                       lower    = LP_out$CI90[, 1, 1],
                       upper    = LP_out$CI90[, 1, 2])

fig0 <- ggplot(plot_df, aes(x = h)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "90% CI"), alpha = 0.2) +
  geom_line(aes(y = lower), colour = blue, linetype = "dotted", linewidth = 1) +
  geom_line(aes(y = upper), colour = blue, linetype = "dotted", linewidth = 1) +
  geom_line(aes(y = estimate, colour = "Point estimate"), linewidth = 1.5) +
  geom_point(aes(y = estimate), colour = blue, size = 3) +
  scale_colour_manual(name = NULL, values = c("Point estimate" = blue)) +
  scale_fill_manual(name = NULL, values = c("90% CI" = blue)) +
  scale_x_continuous(breaks = horizons, limits = c(0, H)) +
  scale_y_continuous(breaks = seq(-0.5, 2.5, by = 0.5), limits = c(-0.7, 2.7)) +
  labs(x = "h", y = "LP estimates of cumulative impulse responses") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(fig0)



###############################################################################
# Purpose: Methane & oxygen dynamics in Eutrophic Lake Model (MODEL2)
# Author:  Daniel Reed (dan.reed@wsu.edu)
# Date:    July 4th 2017
###############################################################################

library(ReacTran)
library(magrittr)
library(ggplot2)
library(cowplot)

# LOAD DATA ----------------------------------------------------------------
# Load profile data
Methane <- read.csv("../data/Methane.csv")
Oxygen <-  read.csv("../data/Oxygen.csv")

# Load rate data
Rates <- read.csv("../data/Rates.csv")

# Load area data
Area <- read.csv("../data/Depth_Area.csv",  stringsAsFactors=FALSE)

# Load mixing data
Mixing <- read.csv("../data/Kz.csv")

# MODEL2 DEFINITION & SOLUTION ----------------------------------------------------------------
#Define the model domain ranging from z=5 to z=15 with a resolution of 1 cm
domain <- setup.grid.1D(x.up=5, x.down=15, dx.1=0.01)

# Create area function to convert depths into areas
Area.fun <- approxfun(Area$Depth..m., Area$Area..m2.)

# Create mixing function
Mixing.func <- approxfun(Mixing$Depth, Mixing$K.z)

# Parameters and boundary conditions
O2.0 <- Oxygen$Oxygen[5] * 1e-3
CH4.L <- tail(Methane$Methane, 1) * 1e-3
D.mix <- Mixing.func(domain$x.int)
Area <- Area.fun(domain$x.int)
Vol <- Area.fun(domain$x.mid) * domain$dx
BDC <- Area * D.mix / domain$dx[1]
k.ox <- 0.1205

#Define model function
Model<-function(t,C,p){
  # Divide up state variables
  O2 <- C[1:domain$N]
  CH4 <- C[(domain$N+1):(domain$N*2)]
  
  # Methane oxidation rate 
  R <- k.ox * CH4 
  
  # Transport processes + reactions
  dO2 <- tran.volume.1D(C = O2, Disp = BDC, C.up = O2.0, 
                        flow.lat = 0, V = Vol, flow = 0)$dC - 2 * R * (O2 / (O2 + 3e-6)) 
  dCH4 <- tran.volume.1D(C = CH4, Disp = BDC, C.down = CH4.L,
                         flow.lat = 0, V = Vol, flow = 0)$dC - R
  
  #Return results
  return(list(c(dO2, dCH4), Rate = R, Rate.O2 = R * O2 / (O2 + 3e-6)))
}

# Solve the model
soln <- steady.1D(y = c(rep(O2.0, domain$N), rep(CH4.L, domain$N)), func=Model,
                  nspec=2, parms = NULL, method = "stode", positive = TRUE)

# Reformat data into data frame
res <- data.frame(O2 = soln$y[1:domain$N], CH4 = soln$y[(domain$N + 1):(2 * domain$N)],
                  Depth  = domain$x.mid, Rate = soln$Rate, Rate.O2 = soln$Rate.O2)

# Plot the results
p1 <- ggplot(res, aes(y = -Depth)) + geom_line(aes(x = O2), lwd = 1.5, colour = "red") 
p1 <- p1 + geom_line(aes(x = CH4), lwd = 1.5, col = "blue")
p1 <- p1 + geom_point(data = Oxygen, aes(y = -Depth, x = Oxygen * 1e-3), col = "red", size = 5)
p1 <- p1 + geom_point(data = Methane, aes(y = -Depth, x = Methane * 1e-3), col = "blue", size = 5)
p1 <- p1 + xlab(expression(CH[4]*","~O[2]~"concentrations (mmol"~L^{-1}*")"))
p1 <- p1 + ylab("Depth (m)") + theme_bw()
p1 <- p1 + scale_y_continuous(breaks = seq(0, -15, by =-5), 
                              labels = seq(0, 15, by = 5), limits = c(-15, -5))

p2 <- ggplot(res, aes(x = -Depth)) + geom_line(aes(y = Rate), lwd = 1.5, colour = "black") 
p2 <- p2 + geom_line(aes(y = Rate.O2), colour = "red", lwd = 1.5, linetype="dashed")
p2 <- p2 + geom_point(data = Rates, aes(x = -Depth, y = Oxidation.rate * 1e-3), 
                      colour = "black", size = 5)
p2 <- p2 + coord_flip() + xlab("") + theme_bw()
p2 <- p2 + ylab(expression(CH[4]~"oxidation rate (mmol"~L^{-1}~d^{-1}*")"))
p2 <- p2 + scale_x_continuous(breaks = seq(0, -15, by =-5), 
                              labels = NULL, limits = c(-15, -5))

plot_grid(p1, p2)
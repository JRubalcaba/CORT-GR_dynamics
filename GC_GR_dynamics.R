############### R Code for the manuscript ###############
# "Modelling the role of glucocorticoid receptor as mediator of endocrine responses to environmental challenge":
# Jimeno B, Rubalcaba JG. 2024. Phil. Trans. R. Soc. B 20220501. https://doi.org/10.1098/rstb.2022.0501

##### Juan G. Rubalcaba: jg.rubalcaba@gmail.com #####

library(deSolve)
require(ggplot2)
require(reshape2)

CORT_dyn <- function(time, state, parameters, signal) {
  p <- signal(time)
  with(as.list(c(state, parameters)),{
    # HORMONES
    dC = -bc * C + p - kcg * O_G - kcm * O_M # CRH
    dA = -ba * A + ka * C - kag * O_G - kam * O_M # ACTH
    dO = -bo * O + ko * A + kmb * O_M + kgb * O_G - kmf * O^n_M * M - kgf * O^n_G * G # CORT
    # RECEPTORS
    dO_G = kgf * O^n_G * G - kgb * O_G
    dG = -kgf * O^n_G * G + kgb * O_G
    dO_M = kmf * O^n_M * M - kmb * O_M
    dM = -kmf * O^n_M * M + kmb * O_M
    list(c(dC,dA,dO,dO_G,dG,dO_M,dM))
  })
}

fit_pars <- c(
  bc=0.73305999, ba=0.44988769, ka=0.31825634, bo=0.01696430, ko=1.19808705, kgb=0.66700166,
  kgf=0.08371887, kmb=0.36783814, kmf=0.43707626, G=1.26916088, M=0.57658257,
  sd_t1=0.84875326, sd_t2=0.66977003, sd_t3=0.68020784, sd_t4=0.41854675
)

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypotalamus
  kcm = 0.1, # prop MR in hypotalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

parameters[12]/parameters[13] # Michaelis K for GR
parameters[14]/parameters[15] # Michaelis K for MR

state <- c(
  C = 1,
  A = 1,
  O = 1,
  O_G = 0,
  fit_pars[10],
  O_M = 0,
  fit_pars[11]
)

signal0 <- function(time, b5=0, alpha=0.02, B=0.3) b5 * alpha * exp(-alpha * time) + B
signal <- function(time, b5=50, alpha=0.02, B=0.3) b5 * alpha * exp(-alpha * time) + B

curve(signal(x, b5=10, alpha = 0.03, B=0.3), 0, 300) # Input stressor (eq. 2.1)

out0 <- ode(y = state, times = seq(0,1000, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
eq <- out0[nrow(out0),2:8]
out0 <- ode(y = eq, times = seq(0,200, by = 1), func = CORT_dyn, parms = parameters, signal = signal)

# Stressor (500x400)
require(ggplot2)
ggplot() + theme_classic() + ylim(0,2) + ylab("") + xlab("Time (min)") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16, colour="black")) +
  geom_line(mapping=aes(y=signal(1:200)-0.3, x=1:200), cex=1) + geom_hline(yintercept = 0.3, cex=1, linetype="dashed")

# CRH
ggplot() + theme_classic() + ylim(0,2) + ylab("[CRH]") + xlab("Time (min)") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16, colour="black")) +
  geom_line(mapping=aes(y=out0[,2], x=0:200), cex=1) 

# ACTH
ggplot() + theme_classic() + ylim(0,2) + ylab("[ACTH]") + xlab("Time (min)") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16, colour="black")) +
  geom_line(mapping=aes(y=out0[,3], x=0:200), cex=1) 

# CORT
ggplot() + theme_classic() + ylim(0,20) + ylab("[GC]") + xlab("Time (min)") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16, colour="black")) +
  geom_line(mapping=aes(y=out0[,4], x=0:200), cex=1) 

# Receptor binding
ggplot() + theme_classic() + ylim(0,1) + ylab("GC - GR Bounds") + xlab("Time (min)") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16, colour="black")) +
  geom_line(mapping=aes(y=out0[,5]/fit_pars[10], x=0:200), cex=1) +
  geom_line(mapping=aes(y=out0[,7]/fit_pars[11], x=0:200), cex=1, linetype="dashed") 

######## Zebra finch data (Jimeno et al. 2019) ######## 

# parameter fitting

fit_pars <- c(
  bc=0.73305999, ba=0.44988769, ka=0.31825634, bo=0.01696430, ko=1.19808705, kgb=0.66700166,
  kgf=0.08371887, kmb=0.36783814, kmf=0.43707626, G=1.26916088, M=0.57658257,
  sd_t1=0.84875326, sd_t2=0.66977003, sd_t3=0.68020784, sd_t4=0.41854675
) # Estimated with ML

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypotalamus
  kcm = 0.1, # prop MR in hypotalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

fit_pars[6]/fit_pars[7] # affinity of GR
fit_pars[8]/fit_pars[9] # affinity of MR

state <- c(
  C = 1,
  A = 1,
  O = 1,
  O_G = 0,
  fit_pars[10],
  O_M = 0,
  fit_pars[11]
)

signal <- function(time, b5=50, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B
signal0 <- function(time, b5=0, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B

out0 <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
eq <- out0[nrow(out0),2:8]
injectevents <- data.frame(var = "A",
                           time = 83,
                           value = 10,
                           method = "add")
out <- ode(y = eq, times = seq(0, 110, by = 1), func = CORT_dyn, parms = parameters, 
           signal = signal, events = list(data = injectevents))

require(ggplot2)
mod <- as.data.frame(out)

# Predicted GC concentration (Fig. 2A)
ggplot() + theme_classic() +
  ylab("log CORT level (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  geom_line(mod, mapping=aes(y=log(O), x=time), cex=1) 

######## _Varying GR levels ####

fit_pars <- c(
  bc=0.73305999, ba=0.44988769, ka=0.31825634, bo=0.01696430, ko=1.19808705, kgb=0.66700166,
  kgf=0.08371887, kmb=0.36783814, kmf=0.43707626, G=1.26916088, M=0.57658257,
  sd_t1=0.84875326, sd_t2=0.66977003, sd_t3=0.68020784, sd_t4=0.41854675
)

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypotalamus
  kcm = 0.1, # prop MR in hypotalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

signal <- function(time, b5=50, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B
signal0 <- function(time, b5=0, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B

injectevents <- data.frame(var = "A",
                           time = 83,
                           value = 10,
                           method = "add") # Simulated ACTH injection

GR_levels <- seq(0.5,2,length.out=10)

out_O <- out_O_G <- array(NA, dim=c(101, length(GR_levels)))
for(i in 1:length(GR_levels)){
  state <- c(
    C = 1,
    A = 1,
    O = 1,
    O_G = 0,
    G = GR_levels[i],
    O_M = 0,
    M = 0.57658257
  )
  # equilibrium
  out <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
  eq <- out[nrow(out),2:8]
  
  # test
  out <- ode(y = eq, times = seq(0, 100, by = 1), func = CORT_dyn, 
             parms = parameters, signal = signal, events = list(data = injectevents))
  
  
  out_O[,i] <- out[,4]
  out_O_G[,i] <- out[,5]
}

out_GR_sample <- out_O
out_GR_sample[-c(1,20,80,100),] <- NA

require(reshape)
plot_GR <- melt(out_O, varnames = c("time","GR"))
out_GR_sample <- out_GR_sample[complete.cases(out_GR_sample),]
out_GR_sample <- melt(out_GR_sample, varnames = c("time","GR"))
out_GR_sample$time[which(out_GR_sample$time==2)] <- 20
out_GR_sample$time[which(out_GR_sample$time==3)] <- 80
out_GR_sample$time[which(out_GR_sample$time==4)] <- 100

t1 <- out_GR_sample[out_GR_sample$time==1,3]
t2 <- out_GR_sample[out_GR_sample$time==20,3]
t3 <- out_GR_sample[out_GR_sample$time==80,3]
t4 <- out_GR_sample[out_GR_sample$time==100,3]

response_plot <- data.frame(CORT=c(t1, t2, t3, t4),
                            GR=c(rep(GR_levels, 4)),
                            response=c(rep("t1",10), rep("t2",10), rep("t3",10), rep("t4",10)))
ggplot() + theme_classic() +
  ylim(-2,60) +
  ylab("% Reduction in CORT level") +
  xlab("GR level") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  geom_line(size = 1) +
  geom_line(response_plot[response_plot$response == "t1",], mapping=aes(y=CORT, x=GR), size=1, colour="red") +
  geom_line(response_plot[response_plot$response == "t2",], mapping=aes(y=CORT, x=GR), size=1, colour="blue") +
  geom_line(response_plot[response_plot$response == "t3",], mapping=aes(y=CORT, x=GR), size=1, colour="green") +
  geom_line(response_plot[response_plot$response == "t4",], mapping=aes(y=CORT, x=GR), size=1, colour="black")

######## House Sparrow data (Zimmer et al. 2021) ######## 

# parameter fitting
fit_pars <- c(bc=0.32629691, ba=0.21617635, ka=0.86906774, bo=0.18375305, ko=0.79438481,
              kgb=0.97259699, kgf=0.09941762, kmb=0.28568222, kmf=1.30058438, G=0.97055172, M=0.31364320,
              sd_t1=0.63780192, sd_t2=0.46723044,
              sd_t3=0.45672276, sd_t4=0.92089914) #Estimated with ML

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypotalamus
  kcm = 0.1, # prop MR in hypotalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

fit_pars[6]/fit_pars[7] # affinity of GR
fit_pars[8]/fit_pars[9] # affinity of MR

state <- c(
  C = 1,
  A = 1,
  O = 1,
  O_G = 0,
  fit_pars[10],
  O_M = 0,
  fit_pars[11]
)

signal <- function(time, b5=50, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B
signal0 <- function(time, b5=0, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B

out0 <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
eq <- out0[nrow(out0),2:8]
out <- ode(y = eq, times = seq(0, 70, by = 1), func = CORT_dyn, parms = parameters, signal = signal)

require(ggplot2)
mod <- as.data.frame(out)

# Predicted GC concentration (Fig. 2B)
ggplot() + theme_classic() +
  ylab("log CORT level (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  geom_line(mod, mapping=aes(y=log(O), x=time), cex=1) 

######## _Varying GR levels ####

GR_levels <- seq(0,1.25,length.out=10)

out_O <- out_O_G <- array(NA, dim=c(71, length(GR_levels)))
for(i in 1:length(GR_levels)){
  state <- c(
    C = 1,
    A = 1,
    O = 1,
    O_G = 0,
    G = GR_levels[i],
    O_M = 0,
    M = 0.31364320
  )
  # equilibrium
  out <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
  eq <- out[nrow(out),2:8]
  
  # test
  out <- ode(y = eq, times = seq(0, 70, by = 1), func = CORT_dyn, 
             parms = parameters, signal = signal)
  
  
  out_O[,i] <- out[,4]
  out_O_G[,i] <- out[,5]
}

out_GR_sample <- out_O
out_GR_sample[-c(1,30,60),] <- NA

require(reshape)
plot_GR <- melt(out_O, varnames = c("time","GR"))
out_GR_sample <- out_GR_sample[complete.cases(out_GR_sample),]
out_GR_sample <- melt(out_GR_sample, varnames = c("time","GR"))
out_GR_sample$time[which(out_GR_sample$time==1)] <- 1
out_GR_sample$time[which(out_GR_sample$time==2)] <- 30
out_GR_sample$time[which(out_GR_sample$time==3)] <- 60

# Predicted
t1 <- out_GR_sample[out_GR_sample$time==1,3]
t2 <- out_GR_sample[out_GR_sample$time==30,3]
t3 <- out_GR_sample[out_GR_sample$time==60,3]

response_plot <- data.frame(CORT=c(t1, t2, t3, t4),
                            GR=c(rep(GR_levels, 4)),
                            response=c(rep("t1",10), rep("t2",10), rep("t3",10), rep("t4",10)))
ggplot() + theme_classic() +
  ylim(0,80) +
  ylab("CORT level") +
  xlab("GR level") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  geom_line(size = 1) +
  geom_line(response_plot[response_plot$response == "t1",], mapping=aes(y=CORT, x=GR), size=1, colour="red") +
  geom_line(response_plot[response_plot$response == "t2",], mapping=aes(y=CORT, x=GR), size=1, colour="green") +
  geom_line(response_plot[response_plot$response == "t3",], mapping=aes(y=CORT, x=GR), size=1, colour="blue") 

######## Tree swallow (Zimmer et al. 2023) ######## 

# parameter fitting
fit_pars <- c(bc=0.38476305, ba=0.18627488, ka=0.41333127, bo=0.17241886, ko=1.89064829,       
              kgb=0.45048750, kgf=0.09156401, kmb=0.21403514, kmf=0.66983096, G=1.24534487,         
              M=0.80172800, sd_t1=1.13742490, sd_t2=0.40118879,
              sd_t3=0.85694958, sd_t4=1.03177214) # Estimated with ML

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypothalamus
  kcm = 0.1, # prop MR in hypothalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

fit_pars[7]/fit_pars[6] # affinity of GR
fit_pars[9]/fit_pars[8] # affinity of MR

state <- c(
  C = 1,
  A = 1,
  O = 1,
  O_G = 0,
  fit_pars[10],
  O_M = 0,
  fit_pars[11]
)

signal <- function(time, b5=50, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B
signal0 <- function(time, b5=0, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B

out0 <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
eq <- out0[nrow(out0),2:8]
out <- ode(y = eq, times = seq(0, 70, by = 1), func = CORT_dyn, parms = parameters, signal = signal)

require(ggplot2)
mod <- as.data.frame(out)

# Predicted GC concentration (Fig. 2C)
ggplot() + theme_classic() +
  ylab("log CORT level (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  geom_line(mod, mapping=aes(y=log(O), x=time), cex=1) 

######## _Varying GR levels ####

GR_levels <- seq(0,1.25,length.out=10)

out_O <- out_O_G <- array(NA, dim=c(71, length(GR_levels)))
for(i in 1:length(GR_levels)){
  state <- c(
    C = 1,
    A = 1,
    O = 1,
    O_G = 0,
    G = GR_levels[i],
    O_M = 0,
    M = 0.80172800
  )
  # equilibrium
  out <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
  eq <- out[nrow(out),2:8]
  
  # test
  out <- ode(y = eq, times = seq(0, 70, by = 1), func = CORT_dyn, 
             parms = parameters, signal = signal)
  
  
  out_O[,i] <- out[,4]
  out_O_G[,i] <- out[,5]
}

out_GR_sample <- out_O
out_GR_sample[-c(1,30,60),] <- NA

require(reshape)
plot_GR <- melt(out_O, varnames = c("time","GR"))
out_GR_sample <- out_GR_sample[complete.cases(out_GR_sample),]
out_GR_sample <- melt(out_GR_sample, varnames = c("time","GR"))
out_GR_sample$time[which(out_GR_sample$time==1)] <- 1
out_GR_sample$time[which(out_GR_sample$time==2)] <- 30
out_GR_sample$time[which(out_GR_sample$time==3)] <- 60

# Predicted
t1 <- out_GR_sample[out_GR_sample$time==1,3]
t2 <- out_GR_sample[out_GR_sample$time==30,3]
t3 <- out_GR_sample[out_GR_sample$time==60,3]

response_plot <- data.frame(CORT=c(t1, t2, t3, t4),
                            GR=c(rep(GR_levels, 4)),
                            response=c(rep("t1",10), rep("t2",10), rep("t3",10), rep("t4",10)))
ggplot() + theme_classic() +
  ylim(-4,60) +
  ylab("% Reduction in CORT level") +
  xlab("GR level") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  geom_line(size = 1) +
  geom_line(response_plot[response_plot$response == "t1",], mapping=aes(y=CORT, x=GR), size=1, colour="red") +
  geom_line(response_plot[response_plot$response == "t2",], mapping=aes(y=CORT, x=GR), size=1, colour="blue") +
  geom_line(response_plot[response_plot$response == "t3",], mapping=aes(y=CORT, x=GR), size=1, colour="green") 

######## _Varying MR levels ####

MR_levels <- seq(0,1.25,length.out=10)

out_O <- out_O_G <- array(NA, dim=c(71, length(GR_levels)))
for(i in 1:length(GR_levels)){
  state <- c(
    C = 1,
    A = 1,
    O = 1,
    O_G = 0,
    G = 1.2453448,
    O_M = 0,
    M = MR_levels[i]
  )
  # equilibrium
  out <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
  eq <- out[nrow(out),2:8]
  
  # test
  out <- ode(y = eq, times = seq(0, 70, by = 1), func = CORT_dyn, 
             parms = parameters, signal = signal)
  
  
  out_O[,i] <- out[,4]
  out_O_G[,i] <- out[,5]
}

out_MR_sample <- out_O
out_MR_sample[-c(1,30,60),] <- NA

require(reshape)
plot_MR <- melt(out_O, varnames = c("time","MR"))
out_MR_sample <- out_MR_sample[complete.cases(out_MR_sample),]
out_MR_sample <- melt(out_MR_sample, varnames = c("time","MR"))
out_MR_sample$time[which(out_MR_sample$time==1)] <- 1
out_MR_sample$time[which(out_MR_sample$time==2)] <- 30
out_MR_sample$time[which(out_MR_sample$time==3)] <- 60

# Predicted
t1 <- out_MR_sample[out_MR_sample$time==1,3]
t2 <- out_MR_sample[out_MR_sample$time==30,3]
t3 <- out_MR_sample[out_MR_sample$time==60,3]

response_plot <- data.frame(CORT=c(t1, t2, t3, t4),
                            MR=c(rep(MR_levels, 4)),
                            response=c(rep("t1",10), rep("t2",10), rep("t3",10), rep("t4",10)))
ggplot() + theme_classic() +
  # ylim(0,30) +
  ylab("CORT level (ng / mL)") +
  xlab("MR level") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  geom_line(size = 1) +
  geom_line(response_plot[response_plot$response == "t1",], mapping=aes(y=CORT, x=MR), size=1, colour="red") +
  geom_line(response_plot[response_plot$response == "t2",], mapping=aes(y=CORT, x=MR), size=1, colour="blue") +
  geom_line(response_plot[response_plot$response == "t3",], mapping=aes(y=CORT, x=MR), size=1, colour="green") 

######## Varying stressor I: 10 stress levels + 2 GR levels ######## 

fit_pars <- c(
  bc=0.73305999, ba=0.44988769, ka=0.31825634, bo=0.01696430, ko=1.19808705, kgb=0.66700166,
  kgf=0.08371887, kmb=0.36783814, kmf=0.43707626, G=1.26916088, M=0.57658257,
  sd_t1=0.84875326, sd_t2=0.66977003, sd_t3=0.68020784, sd_t4=0.41854675
)

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypotalamus
  kcm = 0.1, # prop MR in hypotalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

signal0 <- function(time, b5=0, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B

stressor <- seq(0, 100, length.out=10)

out_O_high <- out_O_low <- out_O_G_high <- out_O_G_low <- array(NA, dim=c(201, length(stressor)))
for(i in 1:length(stressor)){
  state_high <- c(
    C = 1,
    A = 1,
    O = 1,
    O_G = 0,
    G = 2,
    O_M = 0,
    M = 0.57658257
  )
  state_low <- c(
    C = 1,
    A = 1,
    O = 1,
    O_G = 0,
    G = 0.5,
    O_M = 0,
    M = 0.57658257
  )
  
  # equilibrium
  out0_high <- ode(y = state_high, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
  eq_high <- out0_high[nrow(out0_high),2:8]
  
  out0_low <- ode(y = state_low, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
  eq_low <- out0_low[nrow(out0_low),2:8]
  
  signal <- function(time, b5=stressor[i], alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B
  
  # test
  out_high <- ode(y = eq_high, times = seq(0, 200, by = 1), func = CORT_dyn, 
             parms = parameters, signal = signal)
  
  out_low <- ode(y = eq_low, times = seq(0, 200, by = 1), func = CORT_dyn, 
                  parms = parameters, signal = signal)
  
  
  out_O_high[,i] <- out_high[,4]
  out_O_G_high[,i] <- out_high[,5]
  
  out_O_low[,i] <- out_low[,4]
  out_O_G_low[,i] <- out_low[,5]
}

plot_GR_high <- melt(out_O_high, varnames = c("time","GR"))
p1 <- ggplot() + theme_classic() +
  ylim(0, 50) +
  ylab("log CORT level (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_high, mapping=aes(y=value, x=time, colour=GR, group=GR), size=1)

plot_GR_low <- melt(out_O_low, varnames = c("time","GR"))
p2 <- ggplot() + theme_classic() +
  ylim(0, 50) +
  ylab("log CORT level (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_low, mapping=aes(y=value, x=time, colour=GR, group=GR), size=1) 


plot_GR_high <- melt(out_O_high, varnames = c("time","GR"))
plot_GR_high$response <- plot_GR_high$value * state_high[5] / (parameters[12]/parameters[13] + plot_GR_high$value)
p3 <- ggplot() + theme_classic() +
  ylim(0, 2) +
  ylab("Physiologcal response (GR-CORT bounds)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_high, mapping=aes(y=response, x=time, colour=GR, group=GR), size=1)

plot_GR_low <- melt(out_O_low, varnames = c("time","GR"))
plot_GR_low$response <- plot_GR_low$value * state_low[5] / (parameters[12]/parameters[13] + plot_GR_low$value)
p4 <- ggplot() + theme_classic() +
  ylim(0, 0.5) +
  ylab("Physiologcal response (GR-CORT bounds)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_low, mapping=aes(y=response, x=time, colour=GR, group=GR), size=1) 


require(ggpubr)
ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)

######## Varying stressor II: 2 stress levels + 10 GR levels ######## 

fit_pars <- c(
  bc=0.73305999, ba=0.44988769, ka=0.31825634, bo=0.01696430, ko=1.19808705, kgb=0.66700166,
  kgf=0.08371887, kmb=0.36783814, kmf=0.43707626, G=1.26916088, M=0.57658257,
  sd_t1=0.84875326, sd_t2=0.66977003, sd_t3=0.68020784, sd_t4=0.41854675
)

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypotalamus
  kcm = 0.1, # prop MR in hypotalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

signal0 <- function(time, b5=0, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B
signal_high <- function(time, b5=100, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B
signal_low <- function(time, b5=20, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B

GR_levels <- seq(0,4,length.out=10)

out_O_high <- out_O_low <- out_O_G_high <- out_O_G_low <- array(NA, dim=c(201, length(GR_levels)))
for(i in 1:length(GR_levels)){
  state <- c(
    C = 1,
    A = 1,
    O = 1,
    O_G = 0,
    G = GR_levels[i],
    O_M = 0,
    M = 0.57658257
  )

  # equilibrium
  out0 <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
  eq <- out0[nrow(out0),2:8]

  # test
  out_high <- ode(y = eq, times = seq(0, 200, by = 1), func = CORT_dyn, 
                  parms = parameters, signal = signal_high)
  
  out_low <- ode(y = eq, times = seq(0, 200, by = 1), func = CORT_dyn, 
                 parms = parameters, signal = signal_low)
  
  
  out_O_high[,i] <- out_high[,4]
  out_O_G_high[,i] <- out_high[,5]
  
  out_O_low[,i] <- out_low[,4]
  out_O_G_low[,i] <- out_low[,5]
}

plot_GR_high <- melt(out_O_high, varnames = c("time","GR"))
p1 <- ggplot() + theme_classic() +
  ylim(0, 50) +
  ylab("log CORT level (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_high, mapping=aes(y=value, x=time, colour=GR, group=GR), size=1)

plot_GR_low <- melt(out_O_low, varnames = c("time","GR"))
p2 <- ggplot() + theme_classic() +
  ylim(0, 50) +
  ylab("log CORT level (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_low, mapping=aes(y=value, x=time, colour=GR, group=GR), size=1) 


plot_GR_high <- melt(out_O_high, varnames = c("time","GR"))
plot_GR_high$response <- plot_GR_high$value * sort(rep(GR_levels, 201)) / (parameters[12]/parameters[13] + plot_GR_high$value)
p3 <- ggplot() + theme_classic() +
  ylim(0, 3) +
  ylab("Physiologcal response (GR-CORT bounds)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_high, mapping=aes(y=response, x=time, colour=GR, group=GR), size=1)

plot_GR_low <- melt(out_O_low, varnames = c("time","GR"))
plot_GR_low$response <- plot_GR_low$value * sort(rep(GR_levels, 201))  / (parameters[12]/parameters[13] + plot_GR_low$value)
p4 <- ggplot() + theme_classic() +
  ylim(0, 3) +
  ylab("Physiologcal response (GR-CORT bounds)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_low, mapping=aes(y=response, x=time, colour=GR, group=GR), size=1) 

require(ggpubr)
ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)

y1 <- signal_high(1:200)
y2 <- signal_low(1:200)

ggplot() + theme_classic() + ylim(0,2) + geom_line(mapping=aes(y=y1, x=1:200), size=1)
ggplot() + theme_classic() + ylim(0,2) + geom_line(mapping=aes(y=y2, x=1:200), size=1)

######## Varying stressor III: 2 stress durations + 10 GR levels ######## 

fit_pars <- c(
  bc=0.73305999, ba=0.44988769, ka=0.31825634, bo=0.01696430, ko=1.19808705, kgb=0.66700166,
  kgf=0.08371887, kmb=0.36783814, kmf=0.43707626, G=1.26916088, M=0.57658257,
  sd_t1=0.84875326, sd_t2=0.66977003, sd_t3=0.68020784, sd_t4=0.41854675
)

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypotalamus
  kcm = 0.1, # prop MR in hypotalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

signal0 <- function(time, b5=0, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B
signal_short <- function(time, b5=3, alpha=0.04, B=0.2) b5 * exp(-alpha * time) + B
signal_long <- function(time, b5=3, alpha=0.01, B=0.2) b5 * exp(-alpha * time) + B

GR_levels <- seq(0,2,length.out=10)

out_O_long <- out_O_short <- out_O_G_long <- out_O_G_short <- array(NA, dim=c(201, length(GR_levels)))
for(i in 1:length(GR_levels)){
  state <- c(
    C = 1,
    A = 1,
    O = 1,
    O_G = 0,
    G = GR_levels[i],
    O_M = 0,
    M = 0.57658257
  )
  
  # equilibrium
  out0 <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
  eq <- out0[nrow(out0),2:8]
  
  # test
  out_long <- ode(y = eq, times = seq(0, 200, by = 1), func = CORT_dyn, 
                  parms = parameters, signal = signal_long)
  
  out_short <- ode(y = eq, times = seq(0, 200, by = 1), func = CORT_dyn, 
                 parms = parameters, signal = signal_short)
  
  
  out_O_long[,i] <- out_long[,4]
  out_O_G_long[,i] <- out_long[,5]
  
  out_O_short[,i] <- out_short[,4]
  out_O_G_short[,i] <- out_short[,5]
}

plot_GR_long <- melt(out_O_long, varnames = c("time","GR"))
p1 <- ggplot() + theme_classic() +
  ylim(0, 100) +
  ylab("log CORT level (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_long, mapping=aes(y=value, x=time, colour=GR, group=GR), size=1)

plot_GR_short <- melt(out_O_short, varnames = c("time","GR"))
p2 <- ggplot() + theme_classic() +
  ylim(0, 100) +
  ylab("log CORT level (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_short, mapping=aes(y=value, x=time, colour=GR, group=GR), size=1) 


plot_GR_long <- melt(out_O_long, varnames = c("time","GR"))
plot_GR_long$response <- plot_GR_long$value * sort(rep(GR_levels, 201)) / (parameters[12]/parameters[13] + plot_GR_long$value)
p3 <- ggplot() + theme_classic() +
  ylim(0, 2) +
  ylab("Physiologcal response (GR-CORT bounds)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_long, mapping=aes(y=response, x=time, colour=GR, group=GR), size=1)

plot_GR_short <- melt(out_O_short, varnames = c("time","GR"))
plot_GR_short$response <- plot_GR_short$value * sort(rep(GR_levels, 201))  / (parameters[12]/parameters[13] + plot_GR_short$value)
p4 <- ggplot() + theme_classic() +
  ylim(0, 2) +
  ylab("Physiologcal response (GR-CORT bounds)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_short, mapping=aes(y=response, x=time, colour=GR, group=GR), size=1) 

require(ggpubr)
ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)

y1 <- signal_long(1:200)
y2 <- signal_short(1:200)

ggplot() + theme_classic() + ylim(0,4) + geom_line(mapping=aes(y=y1, x=1:200), size=1)
ggplot() + theme_classic() + ylim(0,4) + geom_line(mapping=aes(y=y2, x=1:200), size=1)

######## Varying stressor with senstivity analysis ######## 

fit_pars <- c(
  bc=0.73305999, ba=0.44988769, ka=0.31825634, bo=0.01696430, ko=1.19808705, kgb=0.66700166,
  kgf=0.08371887, kmb=0.36783814, kmf=0.43707626, G=1.26916088, M=0.57658257,
  sd_t1=0.84875326, sd_t2=0.66977003, sd_t3=0.68020784, sd_t4=0.41854675
)

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypotalamus
  kcm = 0.1, # prop MR in hypotalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

signal0 <- function(time, b5=0, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B

GR_levels <- seq(0,4,length.out=10)

stressor <- array(NA, dim=c(50,2))
stressor[1:50,1] <- seq(0, 200, length.out=50) # parameter b5
stressor[1:50,2] <- 0.02 # parameter alpha

out_O <- out_response <- array(NA, dim=c(50, length(GR_levels)))
for(i in 1:50){
  for(j in 1:length(GR_levels)){
    state <- c(
      C = 1,
      A = 1,
      O = 1,
      O_G = 0,
      G = GR_levels[j],
      O_M = 0,
      M = 0.57658257
    )
    
    # equilibrium
    out0 <- ode(y = state, times = seq(0,500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
    eq <- out0[nrow(out0),2:8]
    
    # test
    signal <- function(time, b5=stressor[i,1], alpha=stressor[i,2], B=0.2) b5 * alpha * exp(-alpha * time) + B
    
    
    out <- ode(y = eq, times = seq(0, 200, by = 1), func = CORT_dyn, 
                    parms = parameters, signal = signal)

    out_O[i,j] <- max(out[,4]) - min(out[,4]) # SET indicator from GC response
    response <- state[5] * out[,5] / (parameters[12]/parameters[13] + out[,5])
    out_response[i,j] <- max(response) - min(response)
  }
}

plot_GR_b5 <- melt(out_O[1:50,], varnames = c("b5","GR"))
plot_GR_b5$b5 <- rep(stressor[,1],10)
p1b <- ggplot() + theme_classic() +
  ylim(0, 120) +
  ylab("Range in CORT level (ng / mL)") +
  xlab("Stressor intensity") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_b5, mapping=aes(y=value, x=b5, colour=GR, group=GR), size=1)

# plot_GR_alpha <- melt(out_O[51:100,], varnames = c("alpha","GR"))
# p2 <- ggplot() + theme_classic() +
#   ylim(0, 140) +
#   ylab("Range in CORT level (ng / mL)") +
#   xlab("Stressor duration") +
#   theme(axis.title = element_text(size=14),
#         axis.text = element_text(size=14, colour="black"), 
#         legend.position = "none") +
#   scale_colour_continuous(low="grey90", high="black")+
#   geom_line(plot_GR_alpha, mapping=aes(y=value, x=alpha, colour=GR, group=GR), size=1)

plot_GR_b5_response <- melt(out_response[1:50,], varnames = c("b5","GR"))
plot_GR_b5_response$b5 <- rep(stressor[,1],10)
p3b <- ggplot() + theme_classic() +
  ylim(0, 1.2) +
  ylab("Range in Physiological response") +
  xlab("Stressor intensity") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(plot_GR_b5_response, mapping=aes(y=value, x=b5, colour=GR, group=GR), size=1)

# plot_GR_alpha_response <- melt(out_response[51:100,], varnames = c("alpha","GR"))
# p4 <- ggplot() + theme_classic() +
#   ylim(0, 0.5) +
#   ylab("Range in CORT level (ng / mL)") +
#   xlab("Stressor duration") +
#   theme(axis.title = element_text(size=14),
#         axis.text = element_text(size=14, colour="black"), 
#         legend.position = "none") +
#   scale_colour_continuous(low="grey90", high="black")+
#   geom_line(plot_GR_alpha_response, mapping=aes(y=value, x=alpha, colour=GR, group=GR), size=1)

require(ggpubr)
ggarrange(p1,p2,p1b,p3,p4,p3b,nrow=2,ncol=3) # 1500 x 650

########  Physiological response vs CORT ######## 

fit_pars <- c(
  bc=0.73305999, ba=0.44988769, ka=0.31825634, bo=0.01696430, ko=1.19808705, kgb=0.66700166,
  kgf=0.08371887, kmb=0.36783814, kmf=0.43707626, G=1.26916088, M=0.57658257,
  sd_t1=0.84875326, sd_t2=0.66977003, sd_t3=0.68020784, sd_t4=0.41854675
)

parameters <- c(
  fit_pars[1], # degradation rate CRH
  kcg = 0.1, # prop GR in hypotalamus
  kcm = 0.1, # prop MR in hypotalamus
  
  fit_pars[2], # degradation rate ACTH
  fit_pars[3], # production of ACTH
  kag = 0.1, # prop GR in pituitary
  kam = 0.1, # prop MR in pituitary
  
  fit_pars[4], # degradation rate CORT 
  fit_pars[5], # prod CORT
  
  n_G = 1, # number of binding sites GR
  n_M = 1, # number of binding sites MR
  fit_pars[6], # releasing (backwards) for GR
  fit_pars[7], # binding (forward) for GR
  fit_pars[8], # releasing (backwards) for MR
  fit_pars[9] # binding (forward) for MR
)

signal <- function(time, b5=50, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B
signal0 <- function(time, b5=0, alpha=0.02, B=0.2) b5 * alpha * exp(-alpha * time) + B

GR_levels <- seq(0.5,2,length.out=10)

out_O <- out_O_G <- array(NA, dim=c(201, length(GR_levels)))
for(i in 1:length(GR_levels)){
  state <- c(
    C = 1,
    A = 1,
    O = 1,
    O_G = 0,
    G = GR_levels[i],
    O_M = 0,
    M = 0.57658257
  )
  # equilibrium
  out <- ode(y = state, times = seq(0, 500, by = 1), func = CORT_dyn, parms = parameters, signal = signal0)
  eq <- out[nrow(out),2:8]
  
  # test
  out <- ode(y = eq, times = seq(0, 200, by = 1), func = CORT_dyn, 
             parms = parameters, signal = signal)
  
  
  out_O[,i] <- out[,4]
  out_O_G[,i] <- out[,5]
}

out_GR_sample <- out_O
out_GR_sample[-c(1,30,120),] <- NA

require(reshape)
plot_GR <- melt(out_O, varnames = c("time","GR"))
out_GR_sample <- out_GR_sample[complete.cases(out_GR_sample),]
out_GR_sample <- melt(out_GR_sample, varnames = c("time","GR"))
out_GR_sample$time[which(out_GR_sample$time==1)] <- 0
out_GR_sample$time[which(out_GR_sample$time==2)] <- 30
out_GR_sample$time[which(out_GR_sample$time==3)] <- 120

require(ggplot2)
require(RColorBrewer) # todas 500 x 450
time1cols <- colorRampPalette(brewer.pal(8, "Blues"))
time20cols <- colorRampPalette(brewer.pal(8, "Greens"))
time80cols <- colorRampPalette(brewer.pal(8, "Oranges"))
# time100cols <- colorRampPalette(brewer.pal(8, "Reds"))
ggplot(plot_GR, aes(y=value, x=time, colour=GR, group=GR)) + theme_classic() +
  ylab("GC Concentration (ng / mL)") +
  xlab("Time (min)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") +
  scale_colour_continuous(low="grey90", high="black")+
  geom_line(size=1) +
  geom_point(out_GR_sample[out_GR_sample$time==0,], mapping=aes(y=value, x=time), colour=time1cols(10),cex=7) +
  geom_point(out_GR_sample[out_GR_sample$time==30,], mapping=aes(y=value, x=time), colour=time20cols(10),cex=7) +
  geom_point(out_GR_sample[out_GR_sample$time==120,], mapping=aes(y=value, x=time), colour=time80cols(10),cex=7) 

require(reshape)
x <- seq(0, 30, length.out=100)
GR_bounds <- array(NA, dim=c(100, 10))
for(i in 1:10){
  GR_bounds[,i] <- GR_levels[i] * x /(parameters[12]/parameters[13] + x)
}
GR_bounds <- melt(GR_bounds, varnames = c("NA","GR"))  
GR_bounds$CORT <- rep(x,10)

# 500 x 450
ggplot() + theme_classic() +
  ylim(0,2) +
  ylab("Physiological response (GC - GR bounds)") +
  xlab("GC concentration (ng / mL)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14, colour="black"), 
        legend.position = "none") + 
  scale_colour_continuous(low="white", high="black") +
  geom_line(GR_bounds, mapping=aes(y=value, x=(CORT), colour=GR, group=GR), size=1) +
  geom_point(mapping=aes(y=out_O_G[1,], x=(out_O[1,])), colour=time1cols(10), cex=7) +
  geom_point(mapping=aes(y=out_O_G[30,], x=(out_O[30,])), colour=time20cols(10), cex=7) +
  geom_point(mapping=aes(y=out_O_G[120,], x=(out_O[120,])), colour=time80cols(10), cex=7) 





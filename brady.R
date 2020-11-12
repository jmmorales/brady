
library(ape)

# prepare data and fit the model ----
filop <- read.nexus("Phylip_tree_with_distances_12Dec2019.nexus")

# cophenetic distance among plants
coP <- cophenetic.phylo(filop)

# rearrange the matrix so that order of names matches the data
tmp <- dimnames(coP)
idsp <- as.numeric(as.factor(tmp[[1]]))

DistP <- matrix(NA, ncol(coP), ncol(coP))
for(i in 1:ncol(coP)){
  for(j in 1:ncol(coP)){
    DistP[idsp[i],idsp[j]] <- coP[i,j]
  }
}


#datos <- read.csv("data/up_data_genus.csv", header = TRUE)
datos <- read.csv("BradyD1.csv", header = TRUE)
head(datos)

pred <- (datos$Pred) # depredated seeds
density <- datos$new_density
lat <- datos$latitude
hem <- as.numeric(as.factor(datos$hemisphere)) - 1 
log_seed_mass <- log(datos$Seed_Mass)
site <- datos$Site
total <- datos$Tot_Seed

pspp <- datos$plant_sp
psp <- unique(pspp)

ns_p <- length(psp) # number of plant species
# nsites <- length(unique(site))


N <- length(pred)
Y <- pred
#K <- 4
#X <- matrix(c(rep(1, N), density, lat, hem), N, K )
X <- as.matrix(cbind(scale(log(density)), scale(lat), hem)) #, scale(log_seed_mass)))
colnames(X) = c("log_density", "lat", "hem") #, "log_seed_m")
K <- dim(X)[2]

J_2 <- as.numeric(as.factor(datos$plant_sp))

N_2 <- max(J_2)


library(rstan)
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tmp = unique(datos$plant_sp)
TT = numeric(length(tmp))
for(i in 1:length(tmp)){
  tmp1 = which(datos$plant_sp == tmp[i])
  TT[i] = log(datos$Seed_Mass[tmp1[1]])
}

pr_dat <- list(N = N,
               Y = Y,
               total = total,
               X = X,
               K = K,
               M_2 = 1,
               J_2 = J_2,
               N_2 = N_2,
               log_seed_mass = as.numeric(scale(TT)),
               DistP = DistP/100)

fit <- stan(file = 'zibin_theta_res.stan', 
            data = pr_dat,
            iter = 50000, 
            thin = 25, 
            chains = 3, 
            control = list(max_treedepth = 15)
            )



fit_summary <- summary(fit)$summary



# figure options
res = 300
he = 17
wi = 17


# check R_hat
hist(fit_summary[,10], 100)

# check n_eff
hist(fit_summary[,9], 100)

tmp <- which(fit_summary[,9] < 1000)

print(fit, pars = c("a", "b",  "beta","betaz","rhosqp", "deltap", "etasqp", "az", "bz"))

res <- summary(fit, pars = c("a", "b", "beta", "rhosqp", "deltap", "etasqp", "zi"))$summary
write.csv(res, file = "reszi.csv")


library(coda)
new_y <- extract(fit, pars = "y_pred")
length(which(rowMeans(new_y$y_pred) > mean(Y)))/nrow(new_y$y_pred)

sds = apply(new_y$y_pred, 1, sd)
length(which(sds > sd(Y)))/nrow(new_y$y_pred)

preds = extract(fit, pars = c("sr", "sr_rep"))
plot(preds$sr, preds$sr_rep)
abline(0,1)



pred <- apply(new_y[[1]], 2, quantile, probs=c(0.025,0.5,0.975)) #the median line with 95% 


tmp <- sort(lat, index.return = TRUE)

plot(lat, Y)
lines(lat[tmp$ix], pred[1, tmp$ix])
lines(lat[tmp$ix], pred[3, tmp$ix])

curve(plogis(res[1,1] + x * res[2,1]), xlim = c(-3,3))

# plant random effects --------------------------------------------

# The covariance between plants decreases exponentially with the squared 
# cophenetic distance among them. The rate of decline is guiven by rho.
# The shape of this decline is a half-gaussian
# etasqp is the maximum covariance among plants
# deltap is the extra covariance of a plant species with itself. 

rq <- extract(fitz, pars = c("rhosqp", "deltap", "etasqp"))

# let's see how is the decrease in covariance with cophenetic distance

xs = seq(0, 4, length.out = 200)
op <- par( las = 1, bty = "n", cex = 1, cex.lab = 1)
plot(xs * 100,  exp(- mean(rq$rhosq) * xs^2) , type = "l", lwd = 3, 
     xlab = "cophenetic distance", ylab = "correlation")
rug(DistP[], lwd = 1)
par(op)



# curve(mean(rq$etasq) * exp(- mean(rq$rhosq) * x^2) + mean(rq$delta), 
#      xlim = c(0,3), lwd = 3, ylim = c(0,7))

nsim = 1500
xs = seq(0, 4, length.out = 200)

co <- matrix(NA, nsim, length(xs))

for(i in 1:nsim){
  co[i, ] <- rq$etasq[i] * exp(- rq$rhosq[i] * xs^2) + rq$delta[i]
} 

library(coda)

lo <- HPDinterval(as.mcmc(co)) 

#tiff(filename="crep.tif", res=res, unit="cm", height=10, width=10)
op <- par( las = 1, bty = "n", cex = 1, cex.lab = 1)
plot(xs * 100, lo[,1], type = "l", ylim = c(0,10), xlab = "cophenetic distance", ylab = "covariance")
hist(DistP, 30, freq = FALSE, main = "", xlab = "cophenetic distance", ylim = c(0,1))
lines(xs * 100, lo[,2])
lines(xs * 100,  exp(- mean(rq$rhosq) * xs^2) , lwd = 3)
par(op)
#dev.off()
#--------


fitf <- stan(file = 'zibin_theta_pp_res.stan', 
            data = pr_dat,
            iter = 50000, 
            thin = 25, 
            chains = 3, 
            control = list(max_treedepth = 15)
)

fit_summary <- summary(fitf)$summary

# check R_hat
hist(fit_summary[,10], 100)

# check n_eff
hist(fit_summary[,9], 100)

tmp <- which(fit_summary[,9] < 1000)

print(fitf, pars = c("a", "b",  "beta","betaz","rhosqp", "deltap", "etasqp", "az", "bz",
                     "rhosqs", "deltas", "etasqs"))

res <- summary(fit, pars = c("a", "b", "beta", "rhosqp", "deltap", "etasqp", "zi"))$summary
write.csv(res, file = "reszi.csv")


library(coda)
new_y <- extract(fitf, pars = "y_pred")
length(which(rowMeans(new_y$y_pred) > mean(Y)))/nrow(new_y$y_pred)

sds = apply(new_y$y_pred, 1, sd)
length(which(sds > sd(Y)))/nrow(new_y$y_pred)

preds = extract(fitf, pars = c("sr", "sr_rep"))
plot(preds$sr, preds$sr_rep)
abline(0,1)

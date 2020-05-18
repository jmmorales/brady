library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

C = as.matrix(read.csv('https://github.com/jmmorales/brady/raw/master/C.csv'))
Ys = read.csv('https://github.com/jmmorales/brady/raw/master/Y.csv')
TT = as.matrix(read.csv('https://github.com/jmmorales/brady/raw/master/TT.csv'))
X = as.matrix(read.csv('https://github.com/jmmorales/brady/raw/master/X.csv'))

dat <- list(N = 249,
                Y = Y$Y,
                total = Y$total,
                X = X,
                K = 4,
                J = 36,
                jj = Y$J_2,
                L = 2,
                TT = TT,
                C = C,
                ones = numeric(36) + 1)

fit <- stan(file='stancode.stan', data = dat, iter = 10, chains = 1)
save(fit, file="fitb.RData")

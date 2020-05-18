library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

C = as.matrix(read.csv('https://github.com/jmmorales/brady/raw/master/C.csv'))
Ys = read.csv('https://github.com/jmmorales/brady/raw/master/Y.csv')
TT = as.matrix(read.csv('https://github.com/jmmorales/brady/raw/master/TT.csv'))
X = as.matrix(read.csv('https://github.com/jmmorales/brady/raw/master/X.csv'))

dat <- list(N = 294,
                Y = Ys$Y,
                total = Ys$total,
                X = X,
                K = 4,
                J = 36,
                jj = Ys$J_2,
                L = 2,
                TT = TT,
                C = C,
                ones = numeric(36) + 1)

fit <- stan(file='zibin_pp_o.stan', data = dat, iter = 100000,thin=100, chains = 3)
save(fit, file="fitb.RData")

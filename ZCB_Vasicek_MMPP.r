rm(list=ls())
library(dplyr)
library(ggplot2)

# parameter
theta = 0.5          # long-run mean
k = 0.5              # mean-reversion speed
sigma = 0.1          # volatility
N = 5                # number of state variable set = 2N+1
h = sigma / sqrt(N*k) # step size
T = 2                 # maturity
num_states = 2*N + 1

##########################################################################################################
# Given r0 and maturity T, find the ZCB price at time 0 (3.6.2節MMPP Approach)
# generator matrix A
A = matrix(0, nrow = 2*N+1, ncol = 2*N+1)
A[1, 1] = -1* (sigma^2 + h*k*(N+1- 1)*h)/(2*h^2)                # -u(1)
A[2, 1] = (sigma^2 - h*k*(N+1- 2)*h)/(2*h^2)                    # d(2)
A[2*N, 2*N+1] = (sigma^2 + h*k*(N+1- 2*N)*h)/(2*h^2)            # u(2N)
A[2*N+1, 2*N+1] = -1 * (sigma^2 - h*k*(N+1- (2*N+1))*h)/(2*h^2) # -d(2N+1)

num = 2*N
for(i in c(2:num)){
    A[i-1, i] = (sigma^2 + h*k*(N+1- (i-1))*h)/(2*h^2)                                            # u(i-1)
    A[i, i] = -1 * ((sigma^2 + h*k*(N+1- (i))*h)/(2*h^2) + (sigma^2 - h*k*(N+1- (i))*h)/(2*h^2))  # -(u(i)+d(i))
    A[i+1, i] = (sigma^2 - h*k*(N+1- (i+1))*h)/(2*h^2)                                            # d(i+1)
}

V = matrix(0, nrow = 2*N+1, ncol = 2*N+1)
for(i in c(1:num_states)){
    V[i,i] = theta - (N+1-i) * h # r(i)
}

Lambda = matrix(0, nrow = 2*N+1, ncol = 1)
for(i in c(1:num_states)){
    Lambda[i,1] = theta - (N+1-i) * h # r(i)
}

P_0 = matrix(0, nrow = 1, ncol = 2*N+1)
P_0[N+1] = 1 # assume r0 = theta, which is the (n+1)th states
q_0 = P_0

q_T = q_0 %*% expm((A-V)*T)
ZCB_0 = sum(q_T)

##########################################################################################################
# after 1 unit time, the probability distribution of r(1)
P_1 = P_0 %*% expm(A*1)
ZCB_1 = NULL
for (i in c(1:num_states)){
    q_0 = matrix(0, nrow = 1, ncol = 2*N+1)
    q_0[i] = 1
    q_T = q_0 %*% expm((A-V)*(T-1))
    ZCB_1 = c(ZCB_1, sum(q_T))
}
ZCB_1

# plot ZCB dist
df1 = rbind(P_1, matrix(ZCB_1, nrow = 1, ncol = num_states)) %>% t(.)
colnames(df1) <- c("prob", "ZCB_price")
df1 %>% as.data.frame(.) %>% ggplot() + geom_col(aes(ZCB_price, prob), col = "blue")

# VaR
cumsum(P_1)

rm(list=ls())
library(dplyr)
library(ggplot2)
library(Matrix)
library("magrittr")

# parameter
theta = 0.5          # long-run mean
k = 0.5              # mean-reversion speed
sigma = 0.1          # volatility
N = 10                # number of state variable set = 2N+1
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

# ZCB price under discrete Vasicek model (discrete OU process)
ZCB_0 = sum(q_T)

# check with the closed-form solution of Vasicek model (continuous)
B = (1-exp(-k*(T-0)))/k
a = (B - (T-0))*(k*(k*theta) - 0.5*sigma^2)/(k^2) - (sigma^2 * B^2)/(4*k)
ZCB_closed_form = exp(a - B*theta) # r(0) = theta


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
as.matrix(df1) %>% as.data.frame(.) %>% ggplot() + geom_col(aes(ZCB_price, prob), col = "blue")

# Since state of r(0) if from small to big, the state of ZCB is from big to small 
# (we have to reverse the sequence, or we will get wrong VaR, like the code below)
cdf = cumsum(P_1) 
plot(ZCB_1, cdf, type = "l") # wrong picture for VaR

# 1-day 95% VaR for loss
df2 = rbind(cumsum(P_1[num_states:1]), matrix(ZCB_1[num_states:1] - ZCB_0, nrow = 1, ncol = num_states)) %>% t(.) %>% as.data.frame(.)
colnames(df2) <- c("CDF", "ZCB_loss")
plot(df2$ZCB_loss, df2$CDF, type = "p") 
lines(df2$ZCB_loss, df2$CDF, type = "l") 
abline(h=0.05, col = "blue", lty = 2)

df2 %>% filter(CDF<=0.05) %>% slice_tail(n=1)
VaR = df2 %>% filter(CDF<=0.05) %>% slice_tail(n=1) %$% ZCB_loss

rm(list=ls())
library(dplyr)
library(ggplot2)
library(Matrix)
library("magrittr")

################################### Vasicek one-factor model #############################################
# parameter
theta = 0.5          # long-run mean
k = 0.5              # mean-reversion speed
sigma = 0.1          # volatility
N = 1                # number of state variable set = 2N+1
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

################################### Vasicek two-factor model #############################################
# Parameter
theta1 = 0.5
theta2 = 0.5
k1 = 0.5
k2 = 0.5
sigma1 = 0.1
sigma2 = 0.1
rho = 0
x10 = theta1
x20 = theta2
N1 = 10
N2 = 10
h1 = sigma1 / sqrt(N1 * k1)
h2 = sigma2 / sqrt(N2 * k2)
T = 2

# r1, r2, r3, r4
{if (rho > 0){
    r1 = rho * sigma1 * sigma2 / (2*h1*h2)
    r3 = rho * sigma1 * sigma2 / (2*h1*h2)
    r2 = 0
    r4 = 0
}else{
    r2 = rho * sigma1 * sigma2 / (2*h1*h2)
    r4 = rho * sigma1 * sigma2 / (2*h1*h2)
    r1 = 0
    r3 = 0
}
}

# lambda function coefficient
a0 = 0
a1 = 1
a2 = 1

##########################################################################################################
# Given r0 and maturity T, find the ZCB price at time 0 (4.3節MMPP Approach)
# generator matrix A (PROBLEM: N1, N2 變大時，ZCB沒有穩定收斂到一個值的感覺~)
A = matrix(0, nrow = (2*N2+1)*(2*N1+1), ncol =  (2*N2+1)*(2*N1+1))
col_num = (2*N2+1)*(2*N1+1)
col = c(1:col_num)

# TYPE1
k = 1 # x1 state
j = 1 # x2 state
A[1,1] =  -1*( ((sigma1^2+h1*k1*h1*(N1+1-k))/(2*h1^2) -r1-r4) + 
                ((sigma2^2+h2*k2*h2*(N2+1-j))/(2*h2^2) -r1-r2) + r1)# -(u1+u2+r1)
A[2,1] = (sigma2^2-h2*k2*h2*(N2+1-(j+1)))/(2*h2^2) -r3-r4 # d2
A[2*N2+2,1] = (sigma1^2-h1*k1*h1*(N1+1-(k+1)))/(2*h1^2) -r2-r3 # d1
A[2*N2+3] = r3 # r3

type1_col = 1
col = col[col %in% type1_col == FALSE]

# TYPE2
m = 2*N2
for (i in c(2:m)){
    k = ifelse(i %% (2*N2+1) == 0, i %/% (2*N2+1), i %/% (2*N2+1) +1) # x1 state
    j = ifelse(i %% (2*N2+1) == 0, 2*N2+1, i %% (2*N2+1)) # x2 state
    
    A[i,i] = -1*(((sigma1^2+h1*k1*h1*(N1+1-k))/(2*h1^2) -r1-r4) + 
                     ((sigma2^2+h2*k2*h2*(N2+1-j))/(2*h2^2) -r1-r2) +
                     ((sigma2^2-h2*k2*h2*(N2+1-j))/(2*h2^2) -r3-r4) + r1 + r4) # -(u1+u2+d2+r1+r4)
    A[i-1, i] = ((sigma2^2+h2*k2*h2*(N2+1-(j-1)))/(2*h2^2) -r1-r2) # u2
    A[i+1, i] = (sigma2^2-h2*k2*h2*(N2+1-(j+1)))/(2*h2^2) -r3-r4 # d2
    A[i+2*N2,i] = r2 # r2
    A[i+2*N2+1,i] = (sigma1^2-h1*k1*h1*(N1+1-(k+1)))/(2*h1^2) -r2-r3 # d1
    A[i+2*N2+2,i] = r3 # r3
    col = col[col %in% i == FALSE]
}
# TYPE3
k = 1
j = 2*N2 + 1
A[2*N2+1, 2*N2+1] = -1*(((sigma1^2+h1*k1*h1*(N1+1-k))/(2*h1^2) -r1-r4) + 
                            ((sigma2^2-h2*k2*h2*(N2+1-j))/(2*h2^2) -r3-r4) + r4) # -(u1+d2+r4)
A[2*N2, 2*N2+1] = ((sigma2^2+h2*k2*h2*(N2+1-(j-1)))/(2*h2^2) -r1-r2)# u2
A[2*N2 + 1 + 2*N2, 2*N2+1] = r2 # r2
A[2*N2 + 1 + 2*N2 + 1, 2*N2+1] =  (sigma1^2-h1*k1*h1*(N1+1-(k+1)))/(2*h1^2) -r2-r3 # d1

type3_col = 2*N2+1
col = col[col %in% type3_col == FALSE]

# TYPE4
aa = 2*N2+1
bb = 2*N1-1
m = 1 + aa * (1:bb)
for (i in m){
    k = ifelse(i %% (2*N2+1) == 0, i %/% (2*N2+1), i %/% (2*N2+1) +1) # x1 state
    j = ifelse(i %% (2*N2+1) == 0, 2*N2+1, i %% (2*N2+1)) # x2 state
    
    A[i,i] = -1*(((sigma1^2+h1*k1*h1*(N1+1-k))/(2*h1^2) -r1-r4) +
                ((sigma2^2+h2*k2*h2*(N2+1-j))/(2*h2^2) -r1-r2) +
                ((sigma1^2-h1*k1*h1*(N1+1-k))/(2*h1^2) -r2-r3) + r1 +r2) # -(u1+u2+d1+r1+r2)
    A[i - 2*N2 - 1, i] =  ((sigma1^2+h1*k1*h1*(N1+1-(k-1)))/(2*h1^2) -r1-r4) # u1
    A[i - 2*N2, i] = r4 # r4
    A[i+1, i] = (sigma2^2-h2*k2*h2*(N2+1-(j+1)))/(2*h2^2) -r3-r4  # d2
    A[i + 2*N2 +1, i] =  ((sigma1^2-h1*k1*h1*(N1+1-(k+1)))/(2*h1^2) -r2-r3)# d1
    A[i + 2*N2 + 2, i] = r3 # r3
    col = col[col %in% i == FALSE]
}

# TYPE6
aa = 2*N2+1
bb = 2*N1
m = aa * (2:bb)
for (i in m){
    k = ifelse(i %% (2*N2+1) == 0, i %/% (2*N2+1), i %/% (2*N2+1) +1) # x1 state
    j = ifelse(i %% (2*N2+1) == 0, 2*N2+1, i %% (2*N2+1)) # x2 state
    
    A[i,i] = -1 * ( ((sigma1^2+h1*k1*h1*(N1+1-k))/(2*h1^2) -r1-r4) + 
                        ((sigma1^2-h1*k1*h1*(N1+1-k))/(2*h1^2) -r2-r3) +
                        ((sigma2^2-h2*k2*h2*(N2+1-j))/(2*h2^2) -r3-r4) + r3 + r4 ) # -(u1+d1+d2+r3+r4)
    A[i - 2*N2 - 2, i] = r1 # r1
    A[i - 2*N2 - 1, i] = ((sigma1^2+h1*k1*h1*(N1+1-(k-1)))/(2*h1^2) -r1-r4) # u1
    A[i-1, i] = ((sigma2^2+h2*k2*h2*(N2+1-(j-1)))/(2*h2^2) -r1-r2) # u2
    A[i+2*N2,i] = r2 # r2
    A[i+2*N2+1,i] =  ((sigma1^2-h1*k1*h1*(N1+1-(k+1)))/(2*h1^2) -r2-r3) # d1
    col = col[col %in% i == FALSE]
}

# TYPE7
k = 2*N2+1
j = 1
type7_col = (2*N2+1)*2*N1 + 1
A[type7_col, type7_col] = -1*( ((sigma2^2+h2*k2*h2*(N2+1-j))/(2*h2^2) -r1-r2) +  
                                   ((sigma1^2-h1*k1*h1*(N1+1-k))/(2*h1^2) -r2-r3) + r2) # -(u2+d1+r2)
A[type7_col - 2*N2-1, type7_col] = ((sigma1^2+h1*k1*h1*(N1+1-(k-1)))/(2*h1^2) -r1-r4) # u1
A[type7_col - 2*N2, type7_col] = r4 #r4
A[type7_col+1, type7_col] =  ((sigma2^2-h2*k2*h2*(N2+1-(j+1)))/(2*h2^2) -r3-r4) #d2

col = col[col %in% type7_col == FALSE]

# TYPE8
m1 = (2*N2+1) * 2*N1 +2
m2 = (2*N2+1)*(2*N1+1) -1
for (i in c(m1:m2)){
    k = ifelse(i %% (2*N2+1) == 0, i %/% (2*N2+1), i %/% (2*N2+1) +1) # x1 state
    j = ifelse(i %% (2*N2+1) == 0, 2*N2+1, i %% (2*N2+1)) # x2 state
    
    A[i,i] = -1*(((sigma2^2+h2*k2*h2*(N2+1-j))/(2*h2^2) -r1-r2) +  
                     ((sigma1^2-h1*k1*h1*(N1+1-k))/(2*h1^2) -r2-r3) +
                     ((sigma2^2-h2*k2*h2*(N2+1-j))/(2*h2^2) -r3-r4) +r2 +r3 ) # -(u2+d1+d2+r2+r3)
    A[i - 2*N2 - 2, i] = r1 # r1
    A[i - 2*N2 - 1, i] =  ((sigma1^2+h1*k1*h1*(N1+1-(k-1)))/(2*h1^2) -r1-r4)# u1
    A[i - 2*N2, i] = r4 # r4
    A[i-1, i] = ((sigma2^2+h2*k2*h2*(N2+1-(j-1)))/(2*h2^2) -r1-r2) # u2
    A[i+1, i] = ((sigma2^2-h2*k2*h2*(N2+1-(j+1)))/(2*h2^2) -r3-r4)# d2
    col = col[col %in% i == FALSE]
}

# TYPE9
k = 2*N2 + 1
j = 2*N2 + 1
type9_col = (2*N2+1)*(2*N1+1)
A[type9_col, type9_col] = -1*(((sigma1^2-h1*k1*h1*(N1+1-k))/(2*h1^2) -r2-r3) +  
                            ((sigma2^2-h2*k2*h2*(N2+1-j))/(2*h2^2) -r3-r4) + r3) # -(d1+d2+r3)
A[type9_col - 2*N2 -2, type9_col] = r1 #r1
A[type9_col - 2*N2 -1, type9_col] = ((sigma1^2+h1*k1*h1*(N1+1-(k-1)))/(2*h1^2) -r1-r4) #u1
A[type9_col-1, type9_col] = ((sigma2^2+h2*k2*h2*(N2+1-(j-1)))/(2*h2^2) -r1-r2)  #u2

col = col[col %in% type9_col == FALSE]

# TYPE5
for (i in col){
    k = ifelse(i %% (2*N2+1) == 0, i %/% (2*N2+1), i %/% (2*N2+1) +1) # x1 state
    j = ifelse(i %% (2*N2+1) == 0, 2*N2+1, i %% (2*N2+1)) # x2 state
    
    A[i,i] = -1 * ( ((sigma1^2+h1*k1*h1*(N1+1-k))/(2*h1^2) -r1-r4) + 
                        ((sigma2^2+h2*k2*h2*(N2+1-j))/(2*h2^2) -r1-r2) +
                        ((sigma1^2-h1*k1*h1*(N1+1-k))/(2*h1^2) -r2-r3) + 
                        ((sigma2^2-h2*k2*h2*(N2+1-j))/(2*h2^2) -r3-r4) +
                        r1 +r2 +r3 +r4) # -(u1+u2+d1+d2+r1+r2+r3+r4)
    A[i - 2*N2 - 2, i] = r1 # r1
    A[i - 2*N2 - 1, i] = ((sigma1^2+h1*k1*h1*(N1+1-(k-1)))/(2*h1^2) -r1-r4) # u1
    A[i - 2*N2, i] = r4 # r4
    A[i-1, i] = ((sigma2^2+h2*k2*h2*(N2+1-(j-1)))/(2*h2^2) -r1-r2) # u2
    A[i+1, i] = ((sigma2^2-h2*k2*h2*(N2+1-(j+1)))/(2*h2^2) -r3-r4) # d2
    A[i+2*N2,i] = r2 # r2
    A[i+2*N2+1,i] = ((sigma1^2-h1*k1*h1*(N1+1-(k+1)))/(2*h1^2) -r2-r3) # d1
    A[i + 2*N2 + 2, i] = r3 # r3
}

# V matrix (Lambda = a0 + a1*x1 + a2*x2)
V = matrix(0, nrow = (2*N2+1)*(2*N1+1), ncol =  (2*N2+1)*(2*N1+1))
m1 = 2*N1+1
m2 = 2*N2+1
num = 1
for (k in c(1:m1)){
    for (j in c(1:m2)){
        V[num,num] = a0 + a1 * (x10 + (k-(N1+1))*h1) + a2 * (x20 + (j-(N2+1))*h2)
        num = num + 1
    }
}

# prob
P_0 = matrix(0, nrow = 1, ncol = (2*N2+1)*(2*N1+1))
N_middle = 0.5*((2*N2+1)*(2*N1+1) + 1)
P_0[N_middle] = 1 # assume r0 = theta, which is the (n+1)th states
q_0 = P_0
q_T = q_0 %*% expm((A-V)*T)

# ZCB price under discrete Vasicek two factor model (discrete OU process)
ZCB_0 = sum(q_T)
ZCB_0

# check with the closed-form solution of Vasicek model (continuous)
B1 = (1 - exp(-k1*(T-0)))/k1
B2 = (1 - exp(-k2*(T-0)))/k2
AA =(B1 - (T-0))*(k1*(k1*theta1) - 0.5*sigma1^2)/(k1^2) - (sigma1^2 * B1^2)/(4*k1) +
    (B2 - (T-0))*(k2*(k2*theta2) - 0.5*sigma2^2)/(k2^2) - (sigma2^2 * B2^2)/(4*k2) +
    rho*sigma1*sigma2*(T-0-B1-B2+ (1-exp(-(k1+k2)*(T-0)))/(k1+k2)) / (k1*k2)

ZCB2_closed_form = exp(AA - B1*x10 - B2*x20)
ZCB2_closed_form
##############################
# after 1 unit time, the probability distribution of r(1)
P_1 = P_0 %*% expm(A*1)
num_states = (2*N2+1)*(2*N1+1)
ZCB_1 = NULL
for (i in c(1:num_states)){
    q_0 = matrix(0, nrow = 1, ncol = (2*N2+1)*(2*N1+1))
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

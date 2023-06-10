rm(list=ls())
library(dplyr)
library(ggplot2)
library(Matrix)
library("magrittr")
library(patchwork)
############################## Parameter ######################################
# OU process
theta1 = 0.02 # long-run mean
theta2 = 0.0261
k1 = 1.1829
k2 = 2.567
sigma1 = 0.01
sigma2 = 0.2392
rho = -0.12547

# short rate function
L = 0.01
U = 0.03

# lambda function
K = 0.025
a = 0.246
b = 4.67
c = -0.253
a_ = 0.039
b_ = 12.835
c_ = -0.243

# test linear
a_linear = 0.222373
b_linear = 7.049177
c_linear = -0.291903

# ita function
alpha = 0.6
beta = 0.4

# initial state
x10 = theta1
x20 = theta2

# number of states (discretization)
N1 = 3
N2 = N1
h1 = sigma1 / sqrt(N1 * k1)
h2 = sigma2 / sqrt(N2 * k2)
T = 1

# jump rate (斜的jump) r1, r2, r3, r4
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

##########################################################################################################
# Given r0 and maturity T, find the ZCB price at time 0 (4.3節MMPP Approach)
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

########################################################
# V matrix (Lambda = a0 + a1*x1 + a2*x2)
#L = c(-0.0025, -0.005, 0, 0.0025, 0.005, 0.0075, 0.01)
#U = c(0.1, 0.1,0.1,0.1,0.1,0.1,0.1)
#len = length(L)

#L = c(-0.0025, -0.0025, -0.0025, -0.0025, -0.0025, -0.0025, -0.0025)
#U = c(0.0025, 0.005,0.01,0.015,0.02,0.04,0.06)
#len = length(L)
# flatten and steepen (relationship btw interest rate and default rate)
#b = c(6,5,4,3,2,1)
#b_ = c(8,9,10,11,12,13)
#len = length(b)

# RR
#beta = seq(0, 0.7, 0.1)
#len = length(beta)
zcb = NULL
risk = NULL
len = 1
for (o in c(1:len)){
    V = matrix(0, nrow = (2*N2+1)*(2*N1+1), ncol =  (2*N2+1)*(2*N1+1))
    m1 = 2*N1+1
    m2 = 2*N2+1
    num = 1
    
    for (k in c(1:m1)){
        for (j in c(1:m2)){
            x1k = x10 + (k-(N1+1))*h1
            x2j = x20 + (j-(N2+1))*h2
            r_t = L * as.numeric(x1k <= L) + x1k * as.numeric(x1k >= L & x1k <= U) + U * as.numeric(x1k >= U)
            #r_t =  L * as.numeric(x1k <= L) + x1k * as.numeric(x1k >= L)
            #r_t = x1k
            lambda_t = (a + b*r_t + c*x2j) * as.numeric(r_t <= K) + (a_ + b_*r_t + c_*x2j) * as.numeric(r_t >= K)
            #lambda_t = a_linear + b_linear * r_t + c_linear * x2j
            ita_t = alpha + beta * exp(-1*lambda_t)
            
            V[num,num] = r_t + lambda_t * (1-ita_t)
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
    zcb = c(zcb, ZCB_0)
    
    ##############################
    # after 1 unit time, the probability distribution of r(1)
    P_1 = P_0 %*% expm(A*1)
    num_states = (2*N2+1)*(2*N1+1)
    ZCB_1 = NULL
    for (i in c(1:num_states)){
        q_0 = matrix(0, nrow = 1, ncol = (2*N2+1)*(2*N1+1))
        q_0[i] = 1
        q_T = q_0 %*% expm((A-V)*(T-1/12))
        ZCB_1 = c(ZCB_1, sum(q_T))
    }
    #ZCB_1
    
    # plot ZCB dist
    df1 = rbind(P_1, matrix(ZCB_1, nrow = 1, ncol = num_states)) %>% t(.)
    colnames(df1) <- c("prob", "ZCB_price")
    #as.matrix(df1) %>% as.data.frame(.) %>% ggplot() + geom_col(aes(ZCB_price, prob), col = "blue")
    
    # Since state of r(0) if from small to big, the state of ZCB is from big to small 
    # (we have to reverse the sequence, or we will get wrong VaR, like the code below)
    #cdf = cumsum(P_1) 
    #plot(ZCB_1, cdf, type = "l") # wrong picture for VaR
    
    # 1-day 95% VaR for loss
    df2 = rbind(cumsum(P_1[num_states:1]), matrix(ZCB_1[num_states:1] - ZCB_0, nrow = 1, ncol = num_states)) %>% t(.) %>% as.data.frame(.)
    colnames(df2) <- c("CDF", "ZCB_loss")
    #plot(df2$ZCB_loss, df2$CDF, type = "p") 
    #lines(df2$ZCB_loss, df2$CDF, type = "l") 
    #abline(h=0.05, col = "blue", lty = 2)
    
    #df2 %>% filter(CDF<=0.05) %>% slice_tail(n=1)
    VaR = df2 %>% filter(CDF<=0.01) %>% slice_tail(n=1) %$% ZCB_loss
    risk = c(risk, VaR)
}

plot_df = bind_cols(as.matrix(L),as.matrix(U), as.matrix(zcb), as.matrix(risk))
names(plot_df) = c("L", "U","ZCB", "VaR")
plot_df$IR = as.factor(paste(100*as.numeric(plot_df$L), 100*as.numeric(plot_df$U), sep = " ~ "))

# different L (TEST DIFFERENT INTERVAL)
# blue is Vasicek; red is CIR
p1 = plot_df %>% 
    ggplot()+geom_point(aes(IR, ZCB), col = "blue") + labs(x = "Interest Rate Interval (%)", y = "ZCB price at time 0 ($)") + theme_bw()

p2 = plot_df %>% 
    ggplot()+geom_point(aes(IR, VaR), col = "blue") + labs(x = "Interest Rate Interval (%)", y = "1day 95% VaR ($)") + theme_bw()
p1 | p2

###################
plot_df= bind_cols(as.matrix(b),as.matrix(b_), as.matrix(zcb), as.matrix(risk))
names(plot_df) = c("b", "b_","ZCB", "VaR")
plot_df

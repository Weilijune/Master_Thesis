## Principal Component Analysis
# https://rpubs.com/skydome20/R-Note7-PCA
# prcomp()：主成份分析的基本函式；plot()：繪製陡坡圖(screet plot)，選擇多少個主成份
# dotchart()：繪製主成份負荷圖(PCA loadings plot)；biplot()：繪製主成份負荷圖(PCA loadings plot)

# data("USArrests")
# head(USArrests)

# CDS Spread
CDS_df <- read.csv("D:/作業/台大課程/論文/data/Sovereign_CDS_quote.csv")
CDS_df$Year = CDS_df$Date %>% as.Date(.) %>% year(.)
CDS_df$Month = CDS_df$Date %>% as.Date(.) %>% month(.)
CDS_before_covid19 = CDS_df %>% filter(Year %in% c(2018,2019))
CDS_covid19 = CDS_df %>% filter(Year %in% c(2020))
CDS_post_covid19 = CDS_df %>% filter(Year %in% c(2021,2022,2023))

str_c(names(CDS_before_covid19 %>% select(-Date, -Year, -Month)), collapse = "+")   #選擇變數 (直接把結果貼到formula那邊)
# before的資料drop Japan (2018 too many na)
pca <- prcomp(formula = ~ Brazil+Colombia+Mexico+Chile+Peru+Panama+Turkey+South_Africa+Poland+Hungary+Qatar+Romania+Bulgaria+Slovakia+Israel+Croatia+South_Korea+China+Malaysia+Philippines+Thailand+Pakistan+Japan,
              data =CDS_df ,                          # 資料
              scale = TRUE)                              # 正規化資料
pca <- prcomp(formula = ~ Brazil+Colombia+Mexico+Chile+Peru+Panama+Turkey+South_Africa+Poland+Hungary+Qatar+Romania+Bulgaria+Slovakia+Israel+Croatia+South_Korea+China+Malaysia+Philippines+Thailand+Pakistan+Japan,
              data = CDS_post_covid19   ,                          # 資料
              scale = TRUE)
pca  

str_c(names(CDS_df %>% select(-Date, -Year, -Month)), collapse = ",")
CDS_america = CDS_df %>% select(Brazil,Colombia,Mexico,Chile,Peru,Panama)
CDS_EMEA = CDS_df %>% select(Turkey,South_Africa,Poland,Hungary,Qatar,Romania,Bulgaria,Slovakia,Israel,Croatia)
CDS_Asia = CDS_df %>% select(South_Korea,China,Malaysia,Philippines,Thailand,Pakistan,Japan)
CDS_G7 = CDS_df %>% select(Japan,USA,Canada,UK,France,Germany,Italy)

str_c(names(CDS_G7), collapse = "+")
pca <- prcomp(formula = ~ Turkey+South_Africa+Poland+Hungary+Qatar+Romania+Bulgaria+Slovakia+Israel+Croatia,
              data =CDS_EMEA ,                          # 資料
              scale = TRUE)     


###############################
# 累積解釋圖
vars <- (pca$sdev)^2  # 從pca中取出標準差(pca$sdev)後再平方，計算variance(特徵值)
vars
props <- vars / sum(vars) # 計算每個主成分的解釋比例 = 各個主成分的特徵值/總特徵值    
props
cumulative.props <- cumsum(props)
cumulative.props
plot(cumulative.props)

################################

# US stock market
us_stock_df <- read.csv("D:/作業/台大課程/論文/data/us_stock_market.csv") %>% na.omit()
names(us_stock_df)
us_stock_df$VIX_changeratio = as.numeric(us_stock_df$VIX_changeratio)
pca <- prcomp(formula = ~ S.P500_return+VIX_changeratio,
              data =us_stock_df ,                          # 資料
              scale = TRUE) 

pca$rotation
us_stock_df$pc1 = 0.7071068*us_stock_df$S.P500_return - 0.7071068*us_stock_df$VIX_changeratio
adf.test(us_stock_df$pc1)

us_stock_df$year = us_stock_df$Date %>% as.Date(.) %>% year(.)
us_stock_df$month = us_stock_df$Date %>% as.Date(.) %>% month(.)
us_stock_df = us_stock_df %>% mutate(quarter = case_when(month %in% c(1,2,3) ~ 1,
                                           month %in% c(4,5,6) ~ 2,
                                           month %in% c(7,8,9) ~ 3,
                                           TRUE ~ 4))
# change frequency
us_stock_df = us_stock_df %>% group_by(year, quarter) %>% summarise(return_mean = mean(S.P500_return),
                                                      VIX_Q3 = quantile(VIX_changeratio, 0.75)) %>% ungroup(.)

macro_df <- read.csv("D:/作業/台大課程/論文/data/brazil_test.csv")
macro_df$change_xrate = macro_df$Brazil_exchange_rate/lag(macro_df$Brazil_exchange_rate) - 1
macro_df$change_fr = macro_df$Brazil_foreign_reserves/lag(macro_df$Brazil_foreign_reserves) - 1
macro_df$GDPgrowth = macro_df$Brazil_realGDP/lag(macro_df$Brazil_realGDP) - 1
macro_df$currencycrisis_index = macro_df$Brazil_foreign_reserves/ macro_df$Brazil_M2

us_stock_df$date = paste(us_stock_df$year, us_stock_df$quarter, sep = "Q")
merged_df = merge(us_stock_df, macro_df, by="date")
names(merged_df)
pca_test <- prcomp(formula = ~ return_mean+VIX_Q3+change_xrate+change_fr+GDPgrowth+currencycrisis_index,
              data =merged_df, scale = TRUE) 

merged_df$PC1 = 0.39937494*merged_df$return_mean -0.48517372*merged_df$VIX_Q3-0.49488493*merged_df$change_xrate+
    0.45737554*merged_df$change_fr+0.38856307*merged_df$GDPgrowth +0.00462474*merged_df$currencycrisis_index

adf.test(merged_df$PC1)

###############################
# how many PC should we choose?
#  Kaiser eigenvalue-greater-than-one rule
plot(pca,type="line", main="Kaiser eigenvalue-greater-than-one rule") 
abline(h=1, col="blue")

###############################
# 累積解釋圖
vars <- (pca$sdev)^2  # 從pca中取出標準差(pca$sdev)後再平方，計算variance(特徵值)
vars
props <- vars / sum(vars) # 計算每個主成分的解釋比例 = 各個主成分的特徵值/總特徵值    
props
cumulative.props <- cumsum(props)
cumulative.props
plot(cumulative.props)

################################
# eigen vector(the linear combination of original variables)，loading factors
pca$rotation

# pick three PC to finish PCA loadings plot
top3.pca.eigenvector <- pca$rotation[, 1:3]
top3.pca.eigenvector
first.pca <- top3.pca.eigenvector[, 1]   #  First PC
second.pca <- top3.pca.eigenvector[, 2]  #  Second Pc
third.pca <- top3.pca.eigenvector[, 3]   #  Third PC

# First PC：Order (from small to large) the loading factor (coefficient)
first.pca[order(first.pca, decreasing=FALSE)] 
# Plot the PCA loading factor
dotchart(first.pca[order(first.pca, decreasing=FALSE)] , main="Loading Plot for PC1", xlab="Variable Loadings", col="red")
# Choose PC1 and PC2 to plot PCA loading factor
biplot(pca, choices=1:2)  

CDS_df %$% sd(Brazil, na.rm = T)
###############################################
# monthly variable_macro v.s. CDS spread's PC1
monthly_df <- read.csv("D:/作業/台大課程/論文/data/monthly_credit_spread.csv")
monthly_df %$% adf.test(HY_spread)
monthly_df$ym = paste(monthly_df$Year, monthly_df$Month, sep = "m")

aa = first.pca %>% as.matrix(., ncol = 1)
aa2 = second.pca %>% as.matrix(., ncol = 1)
aa3 = third.pca %>% as.matrix(., ncol = 1)
bb = CDS_df %>% select(.,Brazil:Japan) %>% as.matrix(.)
CDS_df$PC1 = bb %*% aa
CDS_df$PC2 = bb %*% aa2
CDS_df$PC3 = bb %*% aa3

#CDS_df = CDS_df %>% na.omit(.)

#test = CDS_df %>% group_by(Year, Month) %>% summarize(PC1_month = mean(PC1), PC2_month = mean(PC2), PC3_month = mean(PC3)) %>% ungroup(.)
#test$ym = paste(test$Year, test$Month, sep = "m")
#test = merge(test, monthly_df, by = "ym")
#test %$% cor(PC1_month, VIX)

#TEST = merge(us_stock_df, CDS_df, by = "Date")
#TEST %$% cor(PC1, VIX_changeratio)


# PLOT
loading_plot = matrix(c(names(CDS_df %>% select(.,Brazil:Japan)), aa,aa2,aa3), ncol = 4) %>% as.data.frame(.)
names(loading_plot) = c("country", "PC1", "PC2", "PC3")
loading_plot = loading_plot %>% pivot_longer(!country) %>% summarize(country, name, value = as.numeric(value)) 
loading_plot %>% ggplot(.) + geom_col(aes(country, value), fill = "purple") + facet_wrap(~name, nrow = 3, scale="free_y") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + labs(y = "factor loadings")

#############################################
# Richardson extrapolation
m1 = N1
m2 = N2
VaR_CONT = (m2*VaR_m2 - m1*VaR_m1) / (m2 - m1)

############################################
# Cross correlation function (CCF) AND Granger Causality Test
rm(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(tseries)
library(zoo)
library(lmtest)
library(patchwork)
library("magrittr")

# Read files
df = read.csv("D:/作業/台大課程/論文/data/default_interest_rate.csv")
df$month = df$date %>% as.Date(.) %>% month(.)
# The lag k value returned by ccf(x, y) estimates the correlation between x[t+k] and y[t].
# 1year default rate and the lagged value (0-12 period) 3 month interest rate 
ccf(df$X1y_dr_NorthAmerica, df$TB3MS, lag.max = 12, plot = F) %>% 
plot(., main = "Cross correlation function",xlim = c(0,12), lwd = 1)

df$date %>% first()
# quarterly data
quarter_month = c(1,4,7,10)
df_q = df %>% filter(month %in% quarter_month)
ccf(df_q$X1y_dr_NorthAmerica, df_q$TB3MS, lag.max = 12, plot = F) %>% 
    plot(., main = "Cross correlation function",xlim = c(0,12), lwd = 1)

df$rmean_lambda = rollapply(df$X1y_dr_NorthAmerica, 3, mean, align = "right", fill = 0)
df_mean_q = df %>% filter(month %in% quarter_month)
ccf(df_mean_q$rmean_lambda, df_mean_q$TB3MS, lag.max = 12, plot = F) %>% 
    plot(., main = "Cross correlation function",xlim = c(0,12), lwd = 1)
ccf(df_mean_q$TB3MS, df_mean_q$rmean_lambda, lag.max = 12, plot = F) %>% 
    plot(., main = "Cross correlation function",xlim = c(0,12), lwd = 1)

grangertest(DTB6~rmean_lambda, order = 6, data = df_mean_q)
grangertest(rmean_lambda~DTB6, order = 6, data = df_mean_q)

lm(rmean_lambda ~ lag(rmean_lambda, 1) +lag(rmean_lambda, 2) +lag(rmean_lambda, 3)+lag(rmean_lambda, 4) +lag(rmean_lambda, 5)+lag(rmean_lambda, 6)   + lag(TB3MS, 0)+lag(TB3MS, 1)+lag(TB3MS, 2)+lag(TB3MS, 3)+lag(TB3MS, 4)+lag(TB3MS, 5)+lag(TB3MS, 6), data = df_mean_q) %>% summary(.)

# 3 month rate and lagged value (0-12 preiod) 1 year default rate
# 過去的違約率與未來的利率的關係 (當過去的違約率越高，政府傾向調整為寬鬆政策，使企業壓力減低)
ccf(df_q$TB3MS, df_q$X1y_dr_NorthAmerica, lag.max = 12, plot = F) %>% 
    plot(., main = "Cross correlation function",xlim = c(0,12), lwd = 1)

# stationary test: ADF test
adf.test(df_q$TB3MS)
adf.test(df_q$X1y_dr_NorthAmerica)
adf.test(df_mean_q$rmean_lambda)

# Granger Causality Test
grangertest(TB3MS~X1y_dr_NorthAmerica, order = 6, data = df_q)
grangertest(X1y_dr_NorthAmerica~TB3MS, order = 6, data = df_q)

##################################################
# scatter plot
df$year = df$date %>% as.Date(.) %>% year(.)
df_dropna = df %>% filter(year != 1983)
df_before2000 = df_dropna %>% filter(year < 2000)
df_after2000 = df_dropna %>% filter(year >= 2000)

df_before2000 = df_before2000 %>% arrange(DTB6)
p1 = df_before2000 %>% ggplot()+geom_point(aes(DTB6, X1y_dr_NorthAmerica), col = "orange")+
    theme_bw() + labs(x = "DTB6 (%)", y = "default rate", title = "default rate v.s 6M interest rate (before 2000)")

df_before2000 = df_before2000 %>% arrange(TB3MS)
p2 = df_before2000 %>% ggplot()+geom_point(aes(TB3MS, X1y_dr_NorthAmerica), col = "blue")+
    theme_bw() + labs(x = "TB3MS (%)", y = "default rate", title = "default rate v.s 3M interest rate (before 2000)")

df_after2000 = df_after2000 %>% arrange(DTB6)
p3 = df_after2000 %>% ggplot()+geom_point(aes(DTB6, X1y_dr_NorthAmerica), col = "orange")+
    theme_bw() + labs(x = "DTB6 (%)", y = "default rate", title = "default rate v.s 6M interest short rate (after 2000)")

df_after2000 = df_after2000 %>% arrange(TB3MS)
p4 = df_after2000 %>% ggplot()+geom_point(aes(TB3MS, X1y_dr_NorthAmerica), col = "blue")+
    theme_bw() + labs(x = "TB3MS (%)", y = "default rate", title = "default rate v.s 3M interest rate (after 2000)")

df_dropna = df_dropna %>% arrange(DTB6)
p5 = df_dropna%>% ggplot()+geom_point(aes(DTB6, X1y_dr_NorthAmerica), col = "orange")+
    theme_bw() + labs(x = "DTB6 (%)", y = "default rate", title = "default rate v.s 6M interest short rate (full sample)")

df_dropna = df_dropna %>% arrange(TB3MS)
p6 = df_dropna %>% ggplot()+geom_point(aes(TB3MS, X1y_dr_NorthAmerica), col = "blue")+
    theme_bw() + labs(x = "TB3MS (%)", y = "default rate", title = "default rate v.s 3M interest rate (full sample)")


(p5|p1|p3)/(p6|p2|p4)

#####################################
# find different coefficient estimate or correlation (higher v.s lower rate)
df_dropna = df_dropna %>% filter(!year %in% c(1984, 1985, 2008, 2009))
df_highrate = df_dropna %>% filter(lag1y_TB3MS >= quantile(df_dropna$lag1y_TB3MS, 0.5))
df_lowrate = df_dropna %>% filter(lag1y_TB3MS < quantile(df_dropna$lag1y_TB3MS, 0.5))
lm(X1y_dr_NorthAmerica~lag1y_TB3MS, data = df_highrate) %>% summary(.)
lm(X1y_dr_NorthAmerica~lag1y_TB3MS, data = df_lowrate) %>% summary(.)
cor(df_highrate$lag1y_TB3MS, df_highrate$X1y_dr_NorthAmerica)
cor(df_lowrate$lag1y_TB3MS, df_lowrate$X1y_dr_NorthAmerica)

df_highrate = df_dropna %>% filter(TB3MS >= quantile(df_dropna$TB3MS, 0.5))
df_lowrate = df_dropna %>% filter(TB3MS < quantile(df_dropna$TB3MS, 0.5))
lm(X1y_dr_NorthAmerica~TB3MS, data = df_highrate) %>% summary(.)
lm(X1y_dr_NorthAmerica~TB3MS, data = df_lowrate) %>% summary(.)
cor(df_highrate$TB3MS, df_highrate$X1y_dr_NorthAmerica)
cor(df_lowrate$TB3MS, df_lowrate$X1y_dr_NorthAmerica)

#########################################################
# PD and RR regression outcome
# Read files
df1 = read.csv("D:/作業/台大課程/論文/data/df_rr_pd.csv") %>% na.omit()

df1 %>% names()
df1$exp_default_1y = exp(-1*df1$X1yPD)
df1$ln_RR = log(df1$RR)

# df1$default_rate = ifelse(df1$X1yPD == 1, 1, -1*log(1-df1$X1yPD))
df1$exp_lambda_1y = exp(-1*df1$X1y_lambda)
df1$exp_lambda_3y = exp(-3*df1$X3y_lambda)
df1$exp_lambda_5y = exp(-5*df1$X5y_lambda)
df1$exp_lambda_10y = exp(-10*df1$X10y_lambda)

lm(RR~X1yPD, data = df1) %>% summary(.)
lm(RR~exp_default_1y, data = df1) %>% summary(.)
lm(ln_RR~X1yPD, data = df1) %>% summary(.)


lm(RR ~ poly(X1y_lambda, 3), data = df1) %>% summary(.)


lm(RR~X1y_lambda, data = df1) %>% summary(.)
lm(RR~exp_lambda_1y, data = df1) %>% summary(.)
lm(ln_RR~X1y_lambda, data = df1) %>% summary(.)

lm(RR~X3y_lambda, data = df1) %>% summary(.)
lm(RR~exp_lambda_3y, data = df1) %>% summary(.)
lm(ln_RR~X3y_lambda, data = df1) %>% summary(.)

lm(RR~X5y_lambda, data = df1) %>% summary(.)
lm(RR~exp_lambda_5y, data = df1) %>% summary(.)
lm(ln_RR~X5y_lambda, data = df1) %>% summary(.)

lm(RR~X10y_lambda, data = df1) %>% summary(.)
lm(RR~exp_lambda_10y, data = df1) %>% summary(.)
lm(ln_RR~X10y_lambda, data = df1) %>% summary(.)

cor(df1$RR, df1$X1yPD)
cor(df1$RR, df1$exp_default)

plot(df1$X1y_lambda, df1$RR, ylim = c(0,1), lwd = 2, 
     main = "1 year Sovereign Default rate (lambda) v.s Recovery rate",xlab = "1y Default rate", ylab = "Recovery rate")
L = seq(0,1,0.01)
lines(L, 0.54898 - 0.27014*L, lwd = 2)
lines(L, 0.1333 + 0.4231*exp(-1*L), col = "red" , lwd = 2)
lines(L, exp(-0.69811 - 0.62011 * L), col = "blue" , lwd = 2)
lines(L, 0.4961 - 0.55325 * L + 0.05661*L^2 + 0.26443*L^3, col = "green" , lwd = 2)

###############################
test = read.csv("D:/作業/台大課程/論文/data/test.csv")
test %>% group_by(Year, Month) %>% summarize(DTB6 = ifelse(first(DTB6) == "#N/A", nth(DTB6, 2), first(DTB6))) %>%
    mutate(Date = paste(Year, Month, 1, sep = "/")) %>% 
    write.csv(., file = "D:/作業/台大課程/論文/data/CDS/6m.csv")

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
CDS_USA <- read.csv("D:/作業/台大課程/論文/data/CDS/usa_180501_230501.csv")
df1 = read.csv("D:/作業/台大課程/論文/data/CDS/usa_cycleindex.csv")
df2 = read.csv("D:/作業/台大課程/論文/data/CDS/usa_EFFR.csv")
df3 = read.csv("D:/作業/台大課程/論文/data/CDS/usa_sp100.csv")
df4 = read.csv("D:/作業/台大課程/論文/data/CDS/usa_stock.csv")
df5 = read.csv("D:/作業/台大課程/論文/data/CDS/usa_ted_spread.csv")

df1$Date = as.Date(df1$Date)
df2$Date = as.Date(df2$Date)
df3$Date = as.Date(df3$Date)
df4$Date = as.Date(df4$Date)
df5$Date = as.Date(df5$Date)

df = merge(df1, df2, by = "Date") %>% merge(., df3, by = "Date") %>% merge(., df4, by = "Date") %>% merge(., df5, by = "Date")
df %>% View(.)

str_c(names(CDS_USA %>% select(-Date)), collapse = "+")
pca <- prcomp(formula = ~ CDS_1Y+CDS_2Y+CDS_3Y+CDS_4Y+CDS_5Y+CDS_7Y+CDS_10Y,
              data =CDS_USA ,                          # 資料
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
# eigen vector(the linear combination of original variables)，loading factors
pca$rotation

# pick three PC to finish PCA loadings plot
top3.pca.eigenvector <- pca$rotation[, 1:3]
top3.pca.eigenvector
first.pca <- top3.pca.eigenvector[, 1]   #  First PC
second.pca <- top3.pca.eigenvector[, 2]  #  Second Pc
third.pca <- top3.pca.eigenvector[, 3]   #  Third PC

aa = first.pca %>% as.matrix(., ncol = 1)
aa2 = second.pca %>% as.matrix(., ncol = 1)
aa3 = third.pca %>% as.matrix(., ncol = 1)
BB = CDS_USA %>% select(.,-Date) %>% as.matrix(.)
CDS_USA$PC1 = BB %*% aa
CDS_USA$PC2 = BB %*% aa2
CDS_USA$PC3 = BB %*% aa3

# PLOT
loading_plot = matrix(c(names(CDS_USA %>% select(.,CDS_1Y:CDS_10Y)), aa,aa2,aa3), ncol = 4) %>% as.data.frame(.)
names(loading_plot) = c("tenor", "PC1", "PC2", "PC3")
loading_plot = loading_plot %>% pivot_longer(!tenor) %>% summarize(tenor, name, value = as.numeric(value)) 
loading_plot$tenor = loading_plot$tenor %>% as.factor(.) %>% fct_relevel(., "CDS_10Y", after = Inf)
loading_plot %>% ggplot(.) + geom_col(aes(tenor, value), fill = "purple") + facet_wrap(~name, nrow = 3, scale="free_y") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + labs(x = "Tenor", y = "Factor Loadings of PC")

CDS_USA %>% View()

#######################################################################
# regression output
CDS_USA$Date = as.Date(CDS_USA$Date)
df = merge(CDS_USA, df, by = "Date")
df %>% names(.)

# S.P500_return + VIX_change + USDX_change
lm(PC1 ~ EFFR + TED_Spread + S.P100PE + S.P500_index + VIX + USDX + BBG_fin_index + cycle_index, data = df) %>% summary(.)
lm(PC2 ~ EFFR + TED_Spread + S.P100PE + S.P500_index + VIX + USDX + BBG_fin_index + cycle_index, data = df) %>% summary(.)
lm(PC3 ~ EFFR + TED_Spread + S.P100PE + S.P500_index + VIX + USDX + BBG_fin_index + cycle_index, data = df) %>% summary(.)

# best subset approach
library(leaps)
regsubsets(PC1 ~ EFFR + TED_Spread + S.P100PE + S.P500_index + VIX + USDX + BBG_fin_index + cycle_index, data = df) %>% summary(.)
regsubsets(PC2 ~ EFFR + TED_Spread + S.P100PE + S.P500_index + VIX + USDX + BBG_fin_index + cycle_index, data = df) %>% summary(.)
regsubsets(PC3 ~ EFFR + TED_Spread + S.P100PE + S.P500_index + VIX + USDX + BBG_fin_index + cycle_index, data = df) %>% summary(.)

lm(PC1 ~ EFFR + cycle_index, data = df) %>% summary(.)
lm(PC2 ~ EFFR + S.P500_index+ USDX, data = df) %>% summary(.)
lm(PC3 ~ EFFR + S.P500_index+ cycle_index, data = df) %>% summary(.)

cor(df$PC1, df$EFFR) # -0.7
cor(df$PC1, df$cycle_index) #0.24
cor(df$PC1, df$TED_Spread) #-0.23
cor(df$PC1, df$USDX) #-0.47
cor(df$PC1, df$S.P100PE) #0.45
cor(df$PC1, df$BBG_fin_index) #0.28
cor(df$PC1, df$S.P500_index) # 0.03
cor(df$PC1, df$VIX) #-0.07

cor(df$PC2, df$EFFR) # 0.5
cor(df$PC2, df$cycle_index) #0.13
cor(df$PC2, df$TED_Spread) #-0.09
cor(df$PC2, df$USDX) #0.2
cor(df$PC2, df$S.P100PE) #-0.17
cor(df$PC2, df$BBG_fin_index) #0.1
cor(df$PC2, df$S.P500_index) # 0.32
cor(df$PC2, df$VIX) #-0.17

cor(df$PC3, df$EFFR) # 0.47
cor(df$PC3, df$cycle_index) #-0.08
cor(df$PC3, df$TED_Spread) #-0.01
cor(df$PC3, df$USDX) #0.22
cor(df$PC3, df$S.P100PE) #-0.13
cor(df$PC3, df$BBG_fin_index) #-0.05
cor(df$PC3, df$S.P500_index) # 0.27
cor(df$PC3, df$VIX) #-0.002


###################################################################################
# read file (with more data from 2013-2023)
test1 = read.csv("D:/作業/台大課程/論文/data/CDS/test1.csv") %>% na.omit()
test2 = read.csv("D:/作業/台大課程/論文/data/CDS/test2.csv") %>% na.omit()
test3 = read.csv("D:/作業/台大課程/論文/data/CDS/test3.csv")%>% na.omit()

# adf test
adf.test(test1$經濟景氣指標)
adf.test(test2$Ted_spread)
adf.test(test2$BBG_fin_index)
adf.test(test3$EFFR)

adf.test(df$cycle_index)
adf.test(df$S.P500_return)
adf.test(df$BBG_fin_index)

######################################################################################
# EFFR v.s 1Y CDS Spread
# http://www.worldgovernmentbonds.com/sovereign-cds/
df$year = df$Date %>% year(.)
df$implied_default_rate = df$CDS_5Y/(1-0.4)/100 # on the hypothesis of 40% recovery rate
df$EFFR = df$EFFR / 100
df %>% arrange (EFFR) %>% ggplot()+geom_point(aes(EFFR, implied_default_rate, col = year), lwd = 2) + theme_bw() + 
    scale_color_continuous(direction = -1, type = 'viridis') +
    labs(x = "Effective Fed Fund Rate", y = "Implied default rate from 5yr USA CDS spread", title = "Relation between default rates and short-term interest rates")
lm(implied_default_rate ~ EFFR + cycle_index, data = df) %>% summary(.)

df_highEFFR = df %>% filter(EFFR >= 0.024)
df_lowEFFR = df %>% filter(EFFR < 0.024)
lm(implied_default_rate ~ EFFR + cycle_index, data = df_highEFFR) %>% summary(.)
lm(implied_default_rate ~ EFFR + cycle_index, data = df_lowEFFR) %>% summary(.)

lm(implied_default_rate ~ EFFR + cycle_index, data = df) %>% summary(.)

cor(df$EFFR, df$cycle_index)

ccf(df$implied_default_rate, df$EFFR, lag.max = 12, plot = F) %>% 
    plot(., main = "Cross correlation function", lwd = 1)

#############################################
# ESTIMATE OU process: est.ou(data, method = "Hessian", days = 360, significanceLevel = 0.95)
library("SMFI5")
est.ou(df$cycle_index, method = "Hessian", days = 360, significanceLevel = 0.95)
est.ou(df$EFFR, method = "Hessian", days = 360, significanceLevel = 0.95)
cor(df$EFFR, df$cycle_index)
sd(df$EFFR)

test3$Date = as.Date(test3$Date, format = "%m/%d/%Y")
test1$Date = as.Date(test1$Date)
test = merge(test1, test3, by = "Date")
cor(test$EFFR, test$經濟景氣指標)

test3$EFFR = test3$EFFR /100
est.ou(test1$經濟景氣指標, method = "Hessian", days = 360, significanceLevel = 0.95)
est.ou(test3$EFFR, method = "Hessian", days = 360, significanceLevel = 0.95)

############################################
# summary statistic
df %>% names()
df %>% select(starts_with("CDS")) %>% pivot_longer(starts_with("CDS")) %>% group_by(name) %>%
    summarize(mean = mean(value, na.rm = T),
              std = sd(value, na.rm = T),
              min = min(value),
              max = max(value),
              median = median(value),
              n = n()
              ) %>% View()

df %>% names()
df %>% select(Date, EFFR, TED_Spread, USDX, S.P500_index, VIX, S.P100PE, cycle_index, BBG_fin_index) %>% pivot_longer(-Date) %>% group_by(name) %>%
    summarize(mean = mean(value, na.rm = T),
              std = sd(value, na.rm = T),
              min = min(value),
              max = max(value),
              median = median(value),
              n = n()
    ) %>% View()
df$Date

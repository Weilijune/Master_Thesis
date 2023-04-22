## Principal Component Analysis
# https://rpubs.com/skydome20/R-Note7-PCA
# prcomp()：主成份分析的基本函式；plot()：繪製陡坡圖(screet plot)，選擇多少個主成份
# dotchart()：繪製主成份負荷圖(PCA loadings plot)；biplot()：繪製主成份負荷圖(PCA loadings plot)

# data("USArrests")
# head(USArrests)

# CDS Spread
CDS_df <- read.csv("D:/作業/台大課程/論文/data/Sovereign_CDS_quote.csv")
CDS_before_covid19 = CDS_df %>% filter(Year %in% c(2018,2019))
CDS_covid19 = CDS_df %>% filter(Year %in% c(2020))
CDS_post_covid19 = CDS_df %>% filter(Year %in% c(2021,2022,2023))

str_c(names(CDS_before_covid19 %>% select(-Date, -Year, -Month)), collapse = "+")   #選擇變數 (直接把結果貼到formula那邊)
# before的資料drop Japan (2018-2019 too many na)
pca <- prcomp(formula = ~ Brazil+Chile+Argentina+Peru+Mexico+Panama+USA+Turkey+South_Africa+Israel+Hungary+Romania+Bulgaria+Poland+Qatar+Slovakia+South_Korea+China+Malaysia+Thailand+Philippines+Pakistan,
              data =CDS_before_covid19 ,                          # 資料
              scale = TRUE)                              # 正規化資料
pca <- prcomp(formula = ~ Brazil+Chile+Argentina+Peru+Mexico+Panama+USA+Turkey+South_Africa+Israel+Hungary+Romania+Bulgaria+Poland+Qatar+Slovakia+Japan+South_Korea+China+Malaysia+Thailand+Philippines+Pakistan,
              data = CDS_post_covid19   ,                          # 資料
              scale = TRUE)
pca  

str_c(names(CDS_df %>% select(-Date, -Year, -Month)), collapse = ",")
CDS_america = CDS_df %>% select(Brazil,Chile,Argentina,Peru,Mexico,Panama,USA)
CDS_EMEA = CDS_df %>% select(Turkey,South_Africa,Israel,Hungary,Romania,Bulgaria,Poland,Qatar,Slovakia)
CDS_Asia = CDS_df %>% select(Japan,South_Korea,China,Malaysia,Thailand,Philippines,Pakistan)

pca <- prcomp(formula = ~ Japan+South_Korea+China+Malaysia+Thailand+Philippines+Pakistan,
              data =CDS_Asia ,                          # 資料
              scale = TRUE)     

# US stock market
us_stock_df <- read.csv("D:/作業/台大課程/論文/data/us_stock_market.csv") %>% na.omit()
names(us_stock_df)
us_stock_df$VIX_changeratio = as.numeric(us_stock_df$VIX_changeratio)
pca <- prcomp(formula = ~ S.P500_return+VIX_changeratio,
              data =us_stock_df ,                          # 資料
              scale = TRUE) 

pca$rotation
us_stock_df$pc1 = 0.7071068*us_stock_df$S.P500_return - 0.7071068*us_stock_df$VIX_changeratio
adf.test(us_stock_df$pc1 )

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

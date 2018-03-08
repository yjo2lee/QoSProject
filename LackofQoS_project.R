Data <- read.table(file ="/Users/sjlee/Desktop/Graduate/intern_etri/171116_1800.csv", header = TRUE, sep = ",")
Final_Data <- read.csv(file = '/Users/sjlee/Desktop/Graduate/intern_etri/ML_Analysis/Final_Data.csv', header = TRUE, sep = ",")
Data <- Data[Data$RB>Data$RNTI,]
Predicted_Data = 0
Predicted_Data <- read.csv(file = '/Users/sjlee/Desktop/Graduate/intern_etri/ML_Analysis/1116_1800_predicted.csv', header = TRUE, sep = ",")

# eigendecomposition
Trans_Data_sum <- Trans_Data_sum[, 1:6]
Trans_Data_sum[["eigenvalues_11"]] <- rep(0,10)
Trans_Data_sum[["eigenvalues_12"]]<-rep(0,10)
Trans_Data_sum[["eigenvalues_22"]]<-rep(0,10)

for(i in 1:25){
        temp_mat <- matrix(c(Data_Final[i,4], Data_Final[i,6], Data_Final[i,6], Data_Final[i,5]), nrow = 2, ncol = 2)
        temp_var <- temp_mat%*%t(temp_mat)
        Data_Final[i,14]<-temp_var[1,1]
        Data_Final[i,15] <- temp_var[1,2]
        Data_Final[i,16]<-temp_var[2,2]
}


for(i in 1:10){
        temp_mat <- matrix(c(Trans_Data_sum[i,4], Trans_Data_sum[i,6], Trans_Data_sum[i,6], Trans_Data_sum[i,5]), nrow = 2, ncol = 2)
        temp_ev<- eigen(temp_mat)
        temp_B <- temp_ev$vectors%*%matrix(c(sqrt(temp_ev$values[1]), 0, 0, sqrt(temp_ev$values[2])), nrow = 2, ncol =2)%*%t(temp_ev$vectors)
        Trans_Data_sum[i,7]<-temp_B[1,1]
        Trans_Data_sum[i, 8] <- temp_B[1,2]
        Trans_Data_sum[i,9]<-temp_B[2,2]
}

Trans_Data_selected <- Trans_Data_sum[c(3,4,5,7,8,9,10),]


#Predicted_Data <- Predicted_Data[1:15,1:7]
Data_Final[["mu_marginal_50"]] = rep(0,25)
Data_Final[["var_marginal_50"]] = rep(0,25)
Data_Final[["mu_marginal_60"]] = rep(0,25)
Data_Final[["var_marginal_60"]] = rep(0,25)
Data_Final[["mu_marginal_70"]] = rep(0,25)
Data_Final[["var_marginal_70"]] = rep(0,25)
Data_Final[["mu_marginal_80"]] = rep(0,25)
Data_Final[["var_marginal_80"]] = rep(0,25)
Data_Final[["mu_marginal_90"]] = rep(0,25)
Data_Final[["var_marginal_90"]] = rep(0,25)
Data_Final["prob_50"] = rep(0,25)
Data_Final["prob_60"] = rep(0,25)
Data_Final["prob_70"] = rep(0,25)
Data_Final["prob_80"] = rep(0,25)
Data_Final["prob_90"] = rep(0,25)


RB_list = c(50,60,70,80,90)



for (j in 1:5){
	RB_temp = RB_list[j]
	RB_temp_trans = log(RB_temp/(100-RB_temp))
	for (i in 1:25){
        temp_matrix = matrix(c(Data_Final[i,4], Data_Final[i, 5], Data_Final[i,5], Data_Final[i,6]), nrow=2, ncol = 2)
        temp_matrix = temp_matrix %*% t(temp_matrix)
        inverse_matrix = solve(temp_matrix)
        Data_Final[i, 13+2*j-1] = Data_Final[i,3] - (inverse_matrix[2,1]/inverse_matrix[2,2])*(RB_temp_trans - Data_Final[i,2])
        Data_Final[i, 14+2*j-1] = 1/inverse_matrix[2,2]
        #Data_Final[i, 8+j]<- pnorm((log(Data_Final[i,1]*2000)-Data_Final[i,7])/sqrt(Data_Final[i,8]))
}
}


for(i in 1:20){
	
}
	
}
RB_temp_trans = log(RB_temp/(selected100-RB_temp))
for (i in 1:20){
        temp_matrix = matrix(c(Data_Final[i,4], Data_Final[i, 5], Data_Final[i,5], Data_Final[i,6]), nrow=2, ncol = 2)
        temp_matrix = temp_matrix %*% t(temp_matrix)
        inverse_matrix = solve(temp_matrix)
        Data_Final[i, 7] = Data_Final[i,3] - (inverse_matrix[2,1]/inverse_matrix[2,2])*(RB_temp_trans - Data_Final[i,2])
        Data_Final[i, 8] = 1/inverse_matrix[2,2]
}

###########################

Data_RB1090[["mu_marginal_10"]] = rep(0,25)
Data_RB1090[["var_marginal_10"]] = rep(0,25)
Data_RB1090[["mu_marginal_20"]] = rep(0,25)
Data_RB1090[["var_marginal_20"]] = rep(0,25)
Data_RB1090[["mu_marginal_30"]] = rep(0,25)
Data_RB1090[["var_marginal_30"]] = rep(0,25)
Data_RB1090[["mu_marginal_40"]] = rep(0,25)
Data_RB1090[["var_marginal_40"]] = rep(0,25)
Data_RB1090[["mu_marginal_50"]] = rep(0,25)
Data_RB1090[["var_marginal_50"]] = rep(0,25)
Data_RB1090[["mu_marginal_60"]] = rep(0,25)
Data_RB1090[["var_marginal_60"]] = rep(0,25)
Data_RB1090[["mu_marginal_70"]] = rep(0,25)
Data_RB1090[["var_marginal_70"]] = rep(0,25)
Data_RB1090[["mu_marginal_80"]] = rep(0,25)
Data_RB1090[["var_marginal_80"]] = rep(0,25)
Data_RB1090[["mu_marginal_90"]] = rep(0,25)
Data_RB1090[["var_marginal_90"]] = rep(0,25)


RB_list = c(10,20,30,40,50,60,70,80,90)

for (j in 1:9){
	RB_temp = RB_list[j]
	RB_temp_trans = log(RB_temp/(100-RB_temp))
	for (i in 1:25){
        temp_matrix = matrix(c(Data_RB1090[i,4], Data_RB1090[i, 5], Data_RB1090[i,5], Data_RB1090[i,6]), nrow=2, ncol = 2)
        temp_matrix = temp_matrix %*% t(temp_matrix)
        inverse_matrix = solve(temp_matrix)
        Data_RB1090[i, 2*j+5] = Data_RB1090[i,3] - (inverse_matrix[2,1]/inverse_matrix[2,2])*(RB_temp_trans - Data_RB1090[i,2])
        Data_RB1090[i, 2*j+6] = 1/inverse_matrix[2,2]
        #Data_Final[i, 8+j]<- pnorm((log(Data_Final[i,1]*2000)-Data_Final[i,7])/sqrt(Data_Final[i,8]))
}
}

for(i in 1:25){
        temp_mat <- matrix(c(Data_Final[i,4], Data_Final[i,6], Data_Final[i,6], Data_Final[i,5]), nrow = 2, ncol = 2)
        temp_var <- temp_mat%*%t(temp_mat)
        Data_Final[i,14]<-temp_var[1,1]
        Data_Final[i,15] <- temp_var[1,2]
        Data_Final[i,16]<-temp_var[2,2]
}


for(i in 1:20){
	
}
	
}
RB_temp_trans = log(RB_temp/(selected100-RB_temp))
for (i in 1:20){
        temp_matrix = matrix(c(Data_Final[i,4], Data_Final[i, 5], Data_Final[i,5], Data_Final[i,6]), nrow=2, ncol = 2)
        temp_matrix = temp_matrix %*% t(temp_matrix)
        inverse_matrix = solve(temp_matrix)
        Data_Final[i, 7] = Data_Final[i,3] - (inverse_matrix[2,1]/inverse_matrix[2,2])*(RB_temp_trans - Data_Final[i,2])
        Data_Final[i, 8] = 1/inverse_matrix[2,2]
}


Data_Final["prob_90"] = rep(0,20)
for(i in 1:20){
	Data_Final[i, 13]<- pnorm((log(Data_Final[i,1]*2000)-Data_Final[i,7])/sqrt(Data_Final[i,8]))
}

When_3 <- data.frame('RB' = Data$RB[Data$RNTI == 3], 'Bit' = Data$Bit_Cell[Data$RNTI == 3])
When_4 <- data.frame('RB' = Data$RB[Data$RNTI == 4], 'Bit' = Data$Bit_Cell[Data$RNTI == 4])
When_5 <- data.frame('RB' = Data$RB[Data$RNTI == 5], 'Bit' = Data$Bit_Cell[Data$RNTI == 5])
When_6 <- data.frame('RB' = Data$RB[Data$RNTI == 6], 'Bit' = Data$Bit_Cell[Data$RNTI == 6])
When_7 <- data.frame('RB' = Data$RB[Data$RNTI == 7], 'Bit' = Data$Bit_Cell[Data$RNTI == 7])
When_8 <- data.frame('RB' = Data$RB[Data$RNTI == 8], 'Bit' = Data$Bit_Cell[Data$RNTI == 8])
When_9 <- data.frame('RB' = Data$RB[Data$RNTI == 9], 'Bit' = Data$Bit_Cell[Data$RNTI == 9])
When_10 <- data.frame('RB' = Data$RB[Data$RNTI == 10], 'Bit' = Data$Bit_Cell[Data$RNTI == 10])
When_11 <- data.frame('RB' = Data$RB[Data$RNTI == 11], 'Bit' = Data$Bit_Cell[Data$RNTI == 11])
When_12 <- data.frame('RB' = Data$RB[Data$RNTI == 12], 'Bit' = Data$Bit_Cell[Data$RNTI == 12])

Data_list <- list(When_3, When_4, When_5, When_6, When_7, When_8, When_9, When_10, When_11, When_12)

# Shapiro_test_RB&TP
Shapiro_test_pvalue<-data.frame('Num_RNTI' = rep(0,10), 'RB_pvalue' = rep(0, 10), 'Bit_pvalue' = rep(0,10), 'result' = rep(0,10))
for(i in 1:10) {
        Shapiro_test_pvalue[i,1] <- i+2
        Shapiro_test_pvalue[i,2] <- round(shapiro.test(Data_list[[i]][,1])$p.value, digits = 4)
        Shapiro_test_pvalue[i,3] <- round(shapiro.test(Data_list[[i]][,2])$p.value, digits = 4)
        Shapiro_test_pvalue[i,4] <- (Shapiro_test_pvalue[i,2]>0.05)&&(Shapiro_test_pvalue[i,3]>0.05)
}

# QQPlot_RB_normality test
QQPlot_RB<-par(mfrow = c(2,5))
for(i in 1:10) {
        qqnorm(Data_list[[i]][,1])
}

# QQPlot_TP_normality test
QQPlot_Bit<-par(mfrow=c(2,5))
for(i in 1:10) {
        qqnorm(Data_list[[i]][,2])
}

# Shapiro_test_RB_Trans
Shapiro_test_RB_Trans <- data.frame('Num_RNTI' = rep(0,10), 'RB_Trans_pvalue' = rep(0,10), 'result' = rep(0,10))
for (i in 1:10) {
        Shapiro_test_RB_Trans[i,1] <- i+2
        Shapiro_test_RB_Trans[i,2] <- round(shapiro.test(Trans_Data_list[[i]][,1])$p.value, digits = 4)
        Shapiro_test_RB_Trans[i,3] <- (Shapiro_test_RB_Trans[i,2]>0.05)
        if (Shapiro_test_RB_Trans[i,3] == 1){
                Shapiro_test_RB_Trans[i,3] <- "TRUE"
        }
}

# QQPlot_RBTrans_normality test
QQPlot_RB_Trans <- par(mfrow = c(2,5))
for (i in 1:10){
        qqnorm(Trans_Data_list[[i]][,1])
}

# Shapiro_test_TP_Trans
Shapiro_test_TP_Trans <- data.frame('Num_RNTI' = rep(0,10), 'TP_Trans_pvalue' = rep(0,10), 'result' = rep(0,10))
for (i in 1:10) {
        Shapiro_test_TP_Trans[i,1] <- i+2
        Shapiro_test_TP_Trans[i,2] <- round(shapiro.test(Trans_Data_list[[i]][,2])$p.value, digits = 4)
        Shapiro_test_TP_Trans[i,3] <- (Shapiro_test_TP_Trans[i,2] > 0.05)
        if (Shapiro_test_TP_Trans[i,3] == 1){
                Shapiro_test_TP_Trans[i,3] <- "TRUE"
        } else if (Shapiro_test_TP_Trans[i,3] == 0) {
                Shapiro_test_TP_Trans[i,3]<- "FALSE"
        }
}

# QQPlot_BTTrans_normality test
QQPlot_RB_Trans <- par(mfrow = c(2,5))
for (i in 1:10) {
        qqnorm(Trans_Data_list[[i]][,2])
}
Max_Bit <- 75236

# RB/TP Transformation
When_3_Trans <- data.frame('RB_Trans' = log(When_3$RB/(100-When_3$RB)), 'Bit_Trans' = log(When_3$Bit))
When_4_Trans <- data.frame('RB_Trans' = log(When_4$RB/(100-When_4$RB)), 'Bit_Trans' = log(When_4$Bit))
When_5_Trans <- data.frame('RB_Trans' = log(When_5$RB/(100-When_5$RB)), 'Bit_Trans' = log(When_5$Bit))
When_6_Trans <- data.frame('RB_Trans' = log(When_6$RB/(100-When_6$RB)), 'Bit_Trans' = log(When_6$Bit))
When_7_Trans <- data.frame('RB_Trans' = log(When_7$RB/(100-When_7$RB)), 'Bit_Trans' = log(When_7$Bit))
When_8_Trans <- data.frame('RB_Trans' = log(When_8$RB/(100-When_8$RB)), 'Bit_Trans' = log(When_8$Bit))
When_9_Trans <- data.frame('RB_Trans' = log(When_9$RB/(100-When_9$RB)), 'Bit_Trans' = log(When_9$Bit))
When_10_Trans <- data.frame('RB_Trans' = log(When_10$RB/(100-When_10$RB)), 'Bit_Trans' = log(When_10$Bit))
When_11_Trans <- data.frame('RB_Trans' = log(When_11$RB/(100-When_11$RB)), 'Bit_Trans' = log(When_11$Bit))
When_12_Trans <- data.frame('RB_Trans' = log(When_12$RB/(100-When_12$RB)), 'Bit_Trans' = log(When_12$Bit))
Trans_Data_list <- list(When_3_Trans, When_4_Trans, When_5_Trans, When_6_Trans, When_7_Trans, When_8_Trans, When_9_Trans, When_10_Trans, When_11_Trans, When_12_Trans)

Trans_Data_sum <- data.frame('RNTI' = 3:12, 'Avg_RB' = rep(0,10),
                             'Avg_Bit' = rep(0,10), 'Var_RB' = rep(0,10), 'Var_Bit' = rep(0,10), 'Covar' = rep(0,10))
for (i in 1:10){
        Trans_Data_sum[i, 1] = i+2
        Trans_Data_sum[i, 2] = mean(Trans_Data_list[[i]][,1])
        Trans_Data_sum[i, 3] = mean(Trans_Data_list[[i]][,2])
        Trans_Data_sum[i, 4] = var(Trans_Data_list[[i]][,1])
        Trans_Data_sum[i, 5] = var(Trans_Data_list[[i]][,2])
        Trans_Data_sum[i, 6] = cov(Trans_Data_list[[i]][,1], Trans_Data_list[[i]][,2])
}

installed.packages('MVN')
library(MVN)

# MVN test
MVN_test <- data.frame('Mardia' = rep(0,10), 'hz' = rep(0,10), 'royston' = rep(0,10))
mardiaTest(Trans_Data_list[[i]], qqplot = FALSE)
hzTest(Trans_Data_list[[i]], qqplot = FALSE)
roystonTest(Trans_Data_list[[i]], qqplot = FALSE)

mardiaTest(When_12_Trans, qqplot = FALSE)
hzTest(When_12_Trans, qqplot = FALSE)
roystonTest(When_12_Trans, qqplot = FALSE)



means <- data.frame('RNTI' = 3:13, 'RB' = rep(0,11), 'Bit' = rep(0,11))

means$RB[1] = mean(When_3_Trans$RB_Trans)
means$Bit[1] = mean(When_3_Trans$Bit_Trans)
means$RB[2] = mean(When_4_Trans$RB_Trans)
means$Bit[2] = mean(When_4_Trans$Bit_Trans)
means$RB[3] = mean(When_5_Trans$RB_Trans)
means$Bit[3] = mean(When_5_Trans$Bit_Trans)
means$RB[4] = mean(When_6_Trans$RB_Trans)
means$Bit[4] = mean(When_6_Trans$Bit_Trans)
means$RB[5] = mean(When_7_Trans$RB_Trans)
means$Bit[5] = mean(When_7_Trans$Bit_Trans)
means$RB[6] = mean(When_8_Trans$RB_Trans)
means$Bit[6] = mean(When_8_Trans$Bit_Trans)
means$RB[7] = mean(When_9_Trans$RB_Trans)
means$Bit[7] = mean(When_9_Trans$Bit_Trans)
means$RB[8] = mean(When_10_Trans$RB_Trans)
means$Bit[8] = mean(When_10_Trans$Bit_Trans)
means$RB[9] = mean(When_11_Trans$RB_Trans)
means$Bit[9] = mean(When_11_Trans$Bit_Trans)
means$RB[10] = mean(When_12_Trans$RB_Trans)
means$Bit[10] = mean(When_12_Trans$Bit_Trans)
means$RB[11] = mean(When_13_Trans$RB_Trans)
means$Bit[11] = mean(When_13_Trans$Bit_Trans)
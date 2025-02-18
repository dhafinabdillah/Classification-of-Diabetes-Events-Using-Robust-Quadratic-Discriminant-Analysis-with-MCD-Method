# Library
library(ggplot2) # Membuat chart
library(RColorBrewer) # Memberi warna pada chart
library(dplyr) # Membuat data frame di uji normal multivariat
library(biotools) # Menghitung nilai chi di uji kesamaan matriks kovarian
library(MASS) # Membuat aturan klasifikasi
library(rrcov) # Memuat MCD estimator
library(mvoutlier) # Mendeteksi outlier

# Import Dataset
Diabetes <- read.csv("D:/Dhafin's/Kuliah/Skripsi/R/diabetes.csv")

# Data Pre processing
Diabetes$Outcome[Diabetes$Outcome == "1"] <- "2"
Diabetes$Outcome[Diabetes$Outcome == "0"] <- "1" # Data Transformation

Diabetes <- Diabetes %>%
  filter(
    !if_any(-Pregnancies, ~ . == 0 | is.na(.))
  )
rownames(Diabetes) <- NULL # Data Cleaning

# Exploratory Data Analysis
colors <- brewer.pal(2, "Set2")

ggplot(Diabetes, aes(x = Outcome, fill = Outcome)) +
  geom_bar(position = "dodge", width = 0.5, color = "black") +
  scale_fill_manual(
    values = brewer.pal(2, "Set2"),
    labels = c("1" = "Negatif Diabetes", "2" = "Positif Diabetes")
  ) +
  labs(
    x = "Kondisi Diabetes", 
    y = "Jumlah Penduduk", 
    fill = "Kondisi",
    title = "Distribusi Kondisi Diabetes Berdasarkan Hasil Pemeriksaan"
  ) +
  geom_text(
    stat = "count", 
    aes(label = ..count..), 
    position = position_dodge(width = 0.5), 
    vjust = -0.5, 
    color = "black", size = 3
  ) +
  theme_minimal() # Histogram Variabel Respon

descriptive_stats <- Diabetes %>%
  group_by(Outcome) %>%                     
  summarise(across(everything(),            
                   list(mean = ~mean(.),    
                        sd = ~sd(.)),       
                   .names = "{col}_{fn}"))

print(descriptive_stats) # Statistika deskriptif variabel penjelas

# Uji Perbedaan Signifikan Antar Kelompok
group1 <- Diabetes[Diabetes$Outcome == "1", c("Pregnancies", "Glucose","BloodPressure","SkinThickness","Insulin","BMI","DiabetesPedigreeFunction","Age")]
group2 <- Diabetes[Diabetes$Outcome == "2", c("Pregnancies", "Glucose","BloodPressure","SkinThickness","Insulin","BMI","DiabetesPedigreeFunction","Age")]

mean_group1 <- colMeans(group1)
mean_group2 <- colMeans(group2) # Vektor rata-rata sampel

cov_group1 <- cov(group1)
cov_group2 <- cov(group2) # Matriks kovarian sampel

n1 <- nrow(group1)
n2 <- nrow(group2)

spooled <- ((n1 - 1) * cov_group1 + (n2 - 1) * cov_group2) / (n1 + n2 - 2) # Calculate the pooled covariance matrix

mean_diff <- as.matrix(mean_group1 - mean_group2)

t2 <- (n1 * n2) / (n1 + n2) * t(mean_diff) %*% solve(spooled) %*% mean_diff
t2 # T2 Hotelling

F_critical <- qf(0.95, 8, 383)

T2_critical <- (8 * (392 - 2)) / (392 - 8 - 1) * F_critical
T2_critical

# Uji Asumsi
md_Diabetes<-mahalanobis(Diabetes[1:8],colMeans(Diabetes[1:8]),cov(Diabetes[1:8]))

chikuadrat_Diabetes <- qchisq((nrow(Diabetes[1:8]) - seq_len(nrow(Diabetes[1:8])) + 0.5) / nrow(Diabetes[1:8]), 
                              df = ncol(Diabetes[1:8]))

qqplot(chikuadrat_Diabetes,md_Diabetes,main="Q-Q Plot",ylab="Jarak Kuadrat",xlab="Chi-square")
abline(a=0,b=1) # Plot normal multivariat

normalitas_Diabetes <- data.frame(md_Diabetes)
normalitas_Diabetes$chikuadrat_Diabetes <- chikuadrat_Diabetes
normalitas_Diabetes$results<-ifelse(normalitas_Diabetes$md_Diabetes <= normalitas_Diabetes$chikuadrat_Diabetes, 'True', 
                                    ifelse(normalitas_Diabetes$chikuadrat_Diabetes > normalitas_Diabetes$md_Diabetes, 'No', 'None'))

sum_md_less_than_chi_Diabetes<-length(which(normalitas_Diabetes$results=="True"))
n_Diabetes<-nrow(normalitas_Diabetes)
(sum_md_less_than_chi_Diabetes/n_Diabetes)*100 # Normal multivariat

boxM <- boxM(as.matrix(Diabetes[1:8]), Diabetes$Outcome)
boxM

v<-0.5*ncol(Diabetes[1:8])*(ncol(Diabetes[1:8])+1)*(2-1)
chisq<-qchisq(c(0.05),df=v,lower.tail=FALSE)
chisq # Kesamaan matriks kovarian

# Detecting Outliers
outlier <- dd.plot(Diabetes[, 1:8], quan = 0.5114796)
outlier_detect <- outlier$outliers

Diabetes$outlier <- outlier_detect
mahalanobis_outlier<-outlier$md.cla
robust_distances<-outlier$md.rob

outlierr <- data.frame(mahalanobis_outlier)
outlierr$robust_distances <- robust_distances

Outliers <- subset(Diabetes, outlier == 'TRUE')
No_Outliers <- subset(Diabetes, outlier == 'FALSE')
str(Outliers)
str(No_Outliers)

# Analisis Diskriminan Kuadratik
det_Sigma1 <- det(cov_group1)
det_Sigma2 <- det(cov_group2) # Determinan matriks kovarian

inv_Sigma1 <- solve(cov_group1) 
inv_Sigma2 <- solve(cov_group2) # Inverse matriks kovarian

term1 <- 0.5 * log(det_Sigma1 / det_Sigma2)
term2 <- 0.5 * (t(mean_group1) %*% inv_Sigma1 %*% mean_group1 - t(mean_group2) %*% inv_Sigma2 %*% mean_group2)
k <- term1 + term2

compute_score <- function(row) {
  lhs <- -0.5 * t(row) %*% (inv_Sigma1 - inv_Sigma2) %*% row +
    (t(mean_group1) %*% inv_Sigma1 - t(mean_group2) %*% inv_Sigma2) %*% row - 
    k
  return(lhs)
} # Aturan klasifikasi analisis diskriminan kuadratik

data_matrix <- as.matrix(Diabetes[1:8])
scores <- apply(data_matrix, 1, compute_score)

predicted_classes <- ifelse(scores >= 0, "1", "2")

adk <- data.frame(
  score = as.numeric(scores),
  decision_boundary = 0,
  predicted_class = predicted_classes,
  actual_class = Diabetes$Outcome
)

head(adk)

# Recall
confusion_matrix <- table(Predicted = adk$predicted_class, Actual = adk$actual_class)
print(confusion_matrix) # Confusion matrix

FN <- confusion_matrix["1", "2"] 
TP <- confusion_matrix["2", "2"] 

Recall_adk <- (TP / (TP + FN))* 100
print(paste("Recall:", Recall_adk,"%")) # Recall

# Analisis diskriminan kuadratik robust
mcd1 <- CovMcd(group1, alpha = 0.5171756, nsamp = 500)
mcd2 <- CovMcd(group2, alpha = 0.5346154, nsamp = 500) # alpha = h1/n

mean_mcd_group1 <- mcd1@center
mean_mcd_group2 <- mcd2@center

cov_mcd_group1 <- mcd1@cov
cov_mcd_group2 <- mcd2@cov

det_mcd_Sigma1 <- det(cov_mcd_group1)
det_mcd_Sigma2 <- det(cov_mcd_group2)

inv_mcd_Sigma1 <- solve(cov_mcd_group1)
inv_mcd_Sigma2 <- solve(cov_mcd_group2)

term1 <- 0.5 * log(det_mcd_Sigma1 / det_mcd_Sigma2)
term2 <- 0.5 * (t(mean_mcd_group1) %*% inv_mcd_Sigma1 %*% mean_mcd_group1 - t(mean_mcd_group2) %*% inv_mcd_Sigma2 %*% mean_mcd_group2)
k <- term1 + term2

compute_score_mcd <- function(row) {
  lhs <- -0.5 * t(row) %*% (inv_mcd_Sigma1 - inv_mcd_Sigma2) %*% row +
    (t(mean_mcd_group1) %*% inv_mcd_Sigma1 - t(mean_mcd_group2) %*% inv_mcd_Sigma2) %*% row -
    k
  return(lhs)
}

scores_mcd <- apply(data_matrix, 1, compute_score_mcd)

predicted_classes_mcd <- ifelse(scores_mcd >= 0, "1", "2")

adkr <- data.frame(
  score = as.numeric(scores_mcd),
  decision_boundary = 0,
  predicted_class = predicted_classes_mcd,
  actual_class = Diabetes$Outcome
)

head(adkr)

# Recall
confusion_matrix_mcd <- table(Predicted = adkr$predicted_class, Actual = adkr$actual_class)
print(confusion_matrix_mcd) # Confusion matrix

FN <- confusion_matrix_mcd["1", "2"] 
TP <- confusion_matrix_mcd["2", "2"] 

Recall_adkr <- (TP / (TP + FN))* 100
print(paste("Recall:", Recall_adkr,"%")) # Recall
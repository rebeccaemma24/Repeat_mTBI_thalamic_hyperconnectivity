# Code accompanying paper "Repeat traumatic brain injury exacerbates acute thalamic hyperconnectivity in humans"
# Average global thalamocortical FC 

library(R.matlab)
library(neuroCombat)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(wesanderson)
library(reshape2)
library(MASS)
library(rstatix)
library(car)
library(factoextra)
library(corrplot)
library(MatchIt)
library(clinfun)


# Data Load and Setup -----------------------------------------------------


# Load in data for HC and single-event TBI
regressors <- read.csv('/data/all_cov.csv')
beta1 <- read.table("/data/avBETA_GM_cortex_001.txt", quote="\"", comment.char="")
beta2 <- read.table("/data/avBETA_GM_cortex_002.txt", quote="\"", comment.char="")
beta3 <- read.table("/data/avBETA_GM_cortex_003.txt", quote="\"", comment.char="")
beta4 <- read.table("/data/avBETA_GM_cortex_004.txt", quote="\"", comment.char="")
beta5 <- read.table("/data/avBETA_GM_cortex_005.txt", quote="\"", comment.char="")
beta6 <- read.table("/data/avBETA_GM_cortex_006.txt", quote="\"", comment.char="")
beta7 <- read.table("/data/avBETA_GM_cortex_007.txt", quote="\"", comment.char="")
beta8 <- read.table("/data/avBETA_GM_cortex_008.txt", quote="\"", comment.char="")
beta9 <- read.table("/data/avBETA_GM_cortex_009.txt", quote="\"", comment.char="")
beta10 <- read.table("/data/avBETA_GM_cortex_010.txt", quote="\"", comment.char="")
beta11 <- read.table("/data/avBETA_GM_cortex_011.txt", quote="\"", comment.char="")
beta12 <- read.table("/data/avBETA_GM_cortex_012.txt", quote="\"", comment.char="")
beta13 <- read.table("/data/avBETA_GM_cortex_013.txt", quote="\"", comment.char="")
beta14 <- read.table("/data/avBETA_GM_cortex_014.txt", quote="\"", comment.char="")
beta15 <- read.table("/data/avBETA_GM_cortex_015.txt", quote="\"", comment.char="")
beta16 <- read.table("/data/avBETA_GM_cortex_016.txt", quote="\"", comment.char="")

beta<-matrix(0, nrow = 16, ncol = 184)
beta <- t(cbind( beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9, beta10, beta11, beta12, beta13, beta14, beta15, beta16))
rm(beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9, beta10, beta11, beta12, beta13, beta14, beta15, beta16)

regressors_repeat <- read.csv("/data/ALL_concussion_cohort.csv")

beta1 <- read.table("/repeat_data/avBETA_GM_cortex_001.txt", quote="\"", comment.char="")
beta2 <- read.table("/repeat_data/avBETA_GM_cortex_002.txt", quote="\"", comment.char="")
beta3 <- read.table("/repeat_data/avBETA_GM_cortex_003.txt", quote="\"", comment.char="")
beta4 <- read.table("/repeat_data/avBETA_GM_cortex_004.txt", quote="\"", comment.char="")
beta5 <- read.table("/repeat_data/avBETA_GM_cortex_005.txt", quote="\"", comment.char="")
beta6 <- read.table("/repeat_data/avBETA_GM_cortex_006.txt", quote="\"", comment.char="")
beta7 <- read.table("/repeat_data/avBETA_GM_cortex_007.txt", quote="\"", comment.char="")
beta8 <- read.table("/repeat_data/avBETA_GM_cortex_008.txt", quote="\"", comment.char="")
beta9 <- read.table("/repeat_data/avBETA_GM_cortex_009.txt", quote="\"", comment.char="")
beta10 <- read.table("/repeat_data/avBETA_GM_cortex_010.txt", quote="\"", comment.char="")
beta11 <- read.table("/repeat_data/avBETA_GM_cortex_011.txt", quote="\"", comment.char="")
beta12 <- read.table("/repeat_data/avBETA_GM_cortex_012.txt", quote="\"", comment.char="")
beta13 <- read.table("/repeat_data/avBETA_GM_cortex_013.txt", quote="\"", comment.char="")
beta14 <- read.table("/repeat_data/avBETA_GM_cortex_014.txt", quote="\"", comment.char="")
beta15 <- read.table("/repeat_data/avBETA_GM_cortex_015.txt", quote="\"", comment.char="")
beta16 <- read.table("/repeat_data/avBETA_GM_cortex_016.txt", quote="\"", comment.char="")

beta_repeat <-matrix(0, nrow = 16, ncol = 23)
beta_repeat <- t(cbind( beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9, beta10, beta11, beta12, beta13, beta14, beta15, beta16))
rm(beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9, beta10, beta11, beta12, beta13, beta14, beta15, beta16)


# Setup regressors
regressors_repeat$Group <- "repeat"
regressors_repeat %>% dplyr::select("subjectId", "Group","Subject.SiteCode", "Subject.Age",  "Subject.Sex", "Subject.EduYrCt") -> reg_r2
colnames(reg_r2) <- names(regressors)[-1]
cohort1 <- cbind(regressors[,-1], t(beta))
cohort2 <- cbind(reg_r2, t(beta_repeat))
ALL_DATA <- rbind(cohort1, cohort2)
ALL_DATA$Sex[ALL_DATA$Sex == "M"] <- "0"
ALL_DATA$Sex[ALL_DATA$Sex == "F"] <- "1"

# Add site codes manufacturer
ALL_DATA$SiteCode[ALL_DATA$SiteCode=="0ade7c"] <- "0ade7c_b"
ALL_DATA$SiteCode[ALL_DATA$SiteCode=="54ceb9"] <- "54ceb9_b"
ALL_DATA$SiteCode[ALL_DATA$SiteCode=="59129a"] <- "59129a_c1"
ALL_DATA$SiteCode[ALL_DATA$SiteCode=="8effee"] <- "8effee_b"
ALL_DATA$SiteCode[ALL_DATA$SiteCode=="9109c8"] <- "9109c8_c"
ALL_DATA$SiteCode[ALL_DATA$SiteCode=="9e6a55"] <- "9e6a55_a"
ALL_DATA$SiteCode[ALL_DATA$SiteCode=="a72b20"] <- "a72b20_b1"
ALL_DATA$SiteCode[ALL_DATA$SiteCode=="ac3478"] <- "ac3478_b"


# Multisite Harmonisation -------------------------------------------------

# SiteCode
ALL_DATA$SiteCode %>%
  as.factor(.) %>%
  as.numeric(.) -> SiteCode
# Model for age, sex, yrEdu, and Group
ALL_DATA2 <- ALL_DATA
ALL_DATA2$X <- 1
ALL_DATA2$Group[ALL_DATA2$Group == "HC"]<-"0"
ALL_DATA2$Group[ALL_DATA2$Group == "mTBI"]<-"1"
ALL_DATA2$Group[ALL_DATA2$Group == "repeat"]<-"2"
reg <- ALL_DATA2[,c("X", "Group", "Age", "Sex")]
mod <- apply(as.matrix.noquote(reg), 2, as.numeric)

beta_all <- cbind(beta, beta_repeat)
beta.harmonized <- neuroCombat(dat=beta_all, batch=SiteCode, mod=mod)
beta.combat <- as.data.frame(t(beta.harmonized[["dat.combat"]]))
colnames(beta.combat) <- paste0("beta", seq(1:16))

# Add harmonised data to existing dataframe and cleanup 
ALL_DATA2[,c(7:22)] <- beta.combat
ALL_DATA2$X <- NULL
names <- c('L-Thalamus', 'R-Thalamus', 'L-Pulvinar', 'L-Anterior', 'L-mDorsal', 'L-vlDorsal', 'L-Central', 'L-vAnterior', 'L-vlVentral', 'R-Pulvinar', 'R-Anterior', 'R-mDorsal', 'R-vlDorsal', 'R-Central', 'R-vAnterior', 'R-vlVentral')
colnames(ALL_DATA2)[c(7:22)] <- names


# Repeat injury = history of 2 or more ------------------------------------

num_repeat <- merge(x=ALL_DATA2, y=regressors_repeat[,c(2,18)], by.x="GUPI", by.y="subjectId", all.x=TRUE)
num_repeat %>% mutate(combined_concussion = dplyr::recode(combined_concussion, "1-2" = "1", "2-3 times" = "2", "4-5"="4", "Several times, no official diagnosis"="3")) -> num_repeat
num_repeat$combined_concussion <- as.numeric(num_repeat$combined_concussion)

# Edit factor structure: HC, mTBI (original), mTBI (2+)
num_repeat$combined_concussion[num_repeat$combined_concussion == 1] <- NA
num_repeat$combined_concussion[num_repeat$Group == "0"] <- 0 
num_repeat$combined_concussion[num_repeat$Group == "1"] <- 1
# Combine 2,3,4,and 5
num_repeat$combined_concussion[num_repeat$combined_concussion > 1] <- 2

# Remove remaining NAs
num_repeat %>% filter(!is.na(.$combined_concussion)) -> num_repeat

# Match for single and repeat injury on age and sex

num_repeat$Group <- as.numeric(num_repeat$Group)
num_repeat %>% group_by(.$GUPI) %>%
  arrange(desc(.$Group)) %>%
  slice(1L) %>%
  ungroup() -> num_repeat2
num_repeat2$`.$GUPI` <- NULL

# Run between single and repeat concussion
num_repeat2 %>% filter(Group > 0) -> num_repeat3
num_repeat3$Group[num_repeat3$Group == 1] <- 0
num_repeat3$Group[num_repeat3$Group == 2] <- 1

# Assess balance pre-matching 
m.out0 <- matchit(Group ~ Age + Sex, data = num_repeat3,
                  method = NULL, distance = "glm")

m.out1 <- matchit(Group ~ Age + Sex, data = num_repeat3,
                  method = "nearest", distance = "glm")
# Checking balance after matching
summary(m.out1)
plot(summary(m.out1))
plot(m.out1, type = "jitter", interactive = FALSE)
plot(m.out1, type = "density", interactive = FALSE,
     which.xs = ~ Age + Sex )

# Extract out these data
m.data <- match.data(m.out1)

# Get vector of included/excluded single event group
v <- as.data.frame(cbind(num_repeat3$GUPI, m.out1$treat, m.out1$weights))
colnames(v) <- c("GUPI", "Group", "Weight")
v %>% filter(Group == "0") -> v
v2 <- merge(x=cohort1[c(77:184),c(1:2)], y=v, by.x="GUPI", by.y="GUPI", all.x=TRUE) 
v2$Weight[is.na(v2$Weight)] <- 0
v2$Group.x <- NULL
v2$Group.y <- NULL

# Redefine factor structure of single and repeat
m.data$Group[m.data$Group == 1] <- 2
m.data$Group[m.data$Group == 0] <- 1

# Original n=76 controls
m.controls <- num_repeat[(num_repeat$Group == 0),]
m.controls$distance <- NA
m.controls$weights <- NA
m.controls$subclass <- NA

matched_data <- rbind(m.controls, m.data)
matched_data <- as.data.frame(matched_data)



# Analysis & Plot ----------------------------------------------------------------


# Statistical tests between groups- Jonckheere-Terpstra test
p2_jt <- replicate(16, 0)
JT <- replicate(16,0)
for (i in 1:16){
  test <- jonckheere.test(matched_data[,(i+6)], matched_data$combined_concussion, alternative = "increasing", nperm=1000)
  p2_jt[i] <- test$p.value
  JT[i] <- test$statistic
}
p2_jt_fdr <- p.adjust(p2_jt, method="fdr")

# Plot
matched_data$combined_concussion <- as.factor(matched_data$combined_concussion)
plot_data <- matched_data[,c(23,7:22)]
plot_data.m <- melt(plot_data)
names(plot_data.m) <- c("Group", "Source", "Value")
plot_data.m$Group <- as.character(plot_data.m$Group)
plot_data.m$Group[plot_data.m$Group == "0"] <- "HC"
plot_data.m$Group[plot_data.m$Group == "1"] <- "mTBI"
plot_data.m$Group[plot_data.m$Group == "2"] <- "Repeat"
plot_data.m$Group <- factor(plot_data.m$Group, levels=c("HC", "mTBI", "Repeat"))

ggplot(plot_data.m, aes(x=Source, y=Value, color=Group)) + 
  theme_classic(base_size = 20) +
  geom_boxplot() +
  #geom_point(size=0.3, position=position_jitterdodge()) +
  labs(x ="", y ="Thalamocortical Functional Connectivity") +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.text= element_text(color="black")) +
  scale_color_manual(values=c("coral1", "cyan3","mediumpurple"))

# Are there explicit differences between the mTBI groups?- linear model
matched_data %>% filter(Group != "0") -> matched_data2
p3 <- replicate(16, 0)
F3 <- replicate(16, 0)
for (i in 1:16){ 
  lm_thal <- lm(matched_data2[,(i+6)] ~ matched_data2$combined_concussion, data = matched_data2)
  aov <- Anova(lm_thal, type=3)
  p3[i] <- aov[2,4]
  F3[i] <- aov[2,3]
}
p3b <- p3[c(6,8,13,15)]
p3_fdr<-p.adjust(p3b, method = "fdr")

# Statistical tests between groups- Jonckheere-Terpstra test
p3_jt <- replicate(16, 0)
JT3<- replicate(16, 0)
for (i in 1:16){
  test <- jonckheere.test(matched_data2[,(i+6)], as.numeric(as.character(matched_data2$combined_concussion)), alternative = "increasing")
  p3_jt[i] <- test$p.value
  JT3[i] <- test$statistic
}
p3b <- p3_jt[c(1,6,8,13,15)]
p3_jt_fdr <- p.adjust(p3b, method="fdr")

# Plot out signifiant variables for manuscript
plot_data.m %>% filter(Source == "L-Thalamus" | Source == "L-vlDorsal" | Source == "L-vAnterior" | Source == "R-vlDorsal" | Source == "R-vAnterior" ) -> plot_data2.m
labels <- data.frame(Source=names[c(1,6,8,13,15)], label=round(p2_jt_fdr[c(1,6,8,13,15)], digits=3))

labels$label <- paste("p= ", labels$label)
myplot <- ggplot(plot_data2.m, aes(x=Group, y=Value, color=Group)) + 
  geom_boxplot(outlier.shape=NA) +
  theme_classic(base_size = 15) +
  geom_point(size=0.5, position=position_jitterdodge()) +
  xlab("") +
  ylab("Thalamocortical Functional Connectivity") +
  scale_color_manual(values=c("coral1", "cyan3","mediumpurple"))+
  theme(axis.text= element_text(color="black"))+
  facet_wrap(~ Source, nrow=1) 
myplot + 
  geom_label(data = labels, aes(label=label), x = Inf, y = -Inf, hjust=1, vjust=-2, inherit.aes = FALSE, label.size = NA)

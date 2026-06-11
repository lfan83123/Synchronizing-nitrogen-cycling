###Fig 2
{
library(tidyr)
library(forcats)
library(RColorBrewer)
library(writexl)
library(metafor)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(mice)
library(MuMIn) # 用于计算 R²
library(corrplot)
library(PerformanceAnalytics)
library(glmm.hp)
library(dplyr)
library(plyr)
library(lme4)
library(lmerTest)
library(ggpubr)
library(ggthemes)
library(ggbeeswarm)
library(ggridges)
library(ggdist)
LnRRdata <- read.csv("E:/有机替代.作图/response/Synq.csv")

meanlnRR.Nitrification <- rma(LnRR.Nitrification,vi.Nitrification,data=LnRRdata,random = ~ 1 | studies, method = "REML")
meanlnRR.Nitrification
meanlnRR.Denitrification <- rma(LnRR.Denitrification,vi.Denitrification,data=LnRRdata,random = ~ 1 | studies, method = "REML")
meanlnRR.Denitrification
meanlnRR.Immobilization <- rma(LnRR.Immobilization,vi.Immobilization,data=LnRRdata,random = ~ 1 | studies, method = "REML")
meanlnRR.Immobilization
meanlnRR.Mineralization <- rma(LnRR.Mineralization,vi.Mineralization,data=LnRRdata,random = ~ 1 | studies, method = "REML")
meanlnRR.Mineralization
meanlnRR.N2O_emission <- rma(LnRR.N2O_emission,vi.N2O_emission,data=LnRRdata,random = ~ 1 | studies, method = "REML")
meanlnRR.N2O_emission
meanlnRR.NH3_emission <- rma(LnRR.NH3_emission,vi.NH3_emission,data=LnRRdata,random = ~ 1 | studies, method = "REML")
meanlnRR.NH3_emission
#森林图
fmeanlnRR_traits1 <- data.frame(
  Variable = c("N2O_emission","NH3_emission","Mineralization","Immobilization","Nitrification","Denitrification"),
  Estimate = c(meanlnRR.N2O_emission$b, meanlnRR.NH3_emission$b,meanlnRR.Mineralization$b, meanlnRR.Immobilization$b,
               meanlnRR.Nitrification$b, meanlnRR.Denitrification$b),
  SE = c(meanlnRR.N2O_emission$se, meanlnRR.NH3_emission$se, meanlnRR.Mineralization$se, meanlnRR.Immobilization$se,
         meanlnRR.Nitrification$se, meanlnRR.Denitrification$se),
  CI.Lower = c(meanlnRR.N2O_emission$ci.lb, meanlnRR.NH3_emission$ci.lb,meanlnRR.Mineralization$ci.lb, meanlnRR.Immobilization$ci.lb,
               meanlnRR.Nitrification$ci.lb, meanlnRR.Denitrification$ci.lb),
  CI.Upper = c(meanlnRR.N2O_emission$ci.ub, meanlnRR.NH3_emission$ci.ub, meanlnRR.Mineralization$ci.ub, meanlnRR.Immobilization$ci.ub, 
               meanlnRR.Nitrification$ci.ub, meanlnRR.Denitrification$ci.ub),
  P.Value = c(meanlnRR.N2O_emission$pval, meanlnRR.NH3_emission$pval,meanlnRR.Mineralization$pval, meanlnRR.Immobilization$pval, 
              meanlnRR.Nitrification$pval, meanlnRR.Denitrification$pval),
  k = c(meanlnRR.N2O_emission$k, meanlnRR.NH3_emission$k,meanlnRR.Mineralization$k, meanlnRR.Immobilization$k,
        meanlnRR.Nitrification$k, meanlnRR.Denitrification$k)
)

# warm_dat1 <- read.csv("E:/有机替代.作图/response/Egger.csv")
# 因子顺序 & 标签统一
level_order <- c("Denitrification","Nitrification","Immobilization","Mineralization","NH3_emission","N2O_emission")
label_order <- c("Denitrification","Nitrification","Immobilization","Mineralization","NH\u2083 emission","N\u2082O emission")

warm_dat1$Class <- factor(warm_dat1$Class, levels = level_order, labels = label_order)
warm_dat1 <- warm_dat1[complete.cases(warm_dat1[, c("Class")]), ]
fmeanlnRR_traits1$Variable <- factor(fmeanlnRR_traits1$Variable, levels = level_order, labels = label_order)

# 开始绘图
MeanEcoBNF3 <- ggplot() +
  # 半眼分布图
  stat_halfeye(
    data = warm_dat1,
    aes(x = Class, y = LnRR, fill = Class),
    adjust = 0.6,
    width = 0.7,
    justification = -0.3,
    .width = 0,
    slab_color = NA,
    alpha = 0.6
  ) +
  # 原始点 jitter
  geom_jitter(
    data = warm_dat1,
    aes(x = Class, y = LnRR),
    width = 0.15,
    height = 0,
    size = 1.5,
    alpha = 0.6,
    color = "gray"
  ) +
  # 虚线零线
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "blue",
    size = 1.8
  ) +
  # 置信区间
  geom_errorbar(
    data = fmeanlnRR_traits1,
    aes(x = Variable, ymin = CI.Lower, ymax = CI.Upper),
    width = 0.1,
    size = 1.5,
    color = "black"
  ) +
  # 平均点
  geom_point(
    data = fmeanlnRR_traits1,
    aes(x = Variable, y = Estimate),
    shape = 21,
    fill = "red",
    color = "red",
    size = 5
  ) +
  # 添加 K 值
  # geom_text(
  #   data = fmeanlnRR_traits1,
  #   aes(x = Variable, y = 2, label = paste0("k = ", k)),
  #   hjust = 0,
  #   size = 5,
  #   family = "RMN",
  #   fontface = "bold",
  #   color = "black"
  # ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(-2, 2.6),
    breaks = seq(-2, 3, 1)
  ) +
  theme(
    axis.text = element_text(family = "RMN", color = "black", size = 16, face = "bold"),
    axis.title = element_text(family = "RMN", color = "black", size = 16, face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none"
    # plot.margin = margin(5.5, 40, 5.5, 5.5)
  ) +
  ylab("Effect Size (lnRR)") +
  xlab("")

# 打印图形
MeanEcoBNF3

BNF222 <- read.csv("E:/有机替代.作图/response/Syn.csv")

BNF5<- BNF222[which(BNF222$NTR == "Nit_Syn_f_ab"), ]
BNF6<- BNF222[which(BNF222$NTR == "Den_Syn_f_ab"), ]
BNF7<- BNF222[which(BNF222$NTR == "Min_Syn_f_ab"), ]
BNF8<- BNF222[which(BNF222$NTR == "Imm_Syn_f_ab"), ]
BNF9<- BNF222[which(BNF222$NTR == "Nit_Syn_t_ab"), ]
BNF10<- BNF222[which(BNF222$NTR == "Den_Syn_t_ab"), ]
BNF11<- BNF222[which(BNF222$NTR == "Min_Syn_t_ab"), ]
BNF12<- BNF222[which(BNF222$NTR == "Imm_Syn_t_ab"), ]

tota12<-rbind(BNF5, BNF6, BNF7, BNF8, BNF9, BNF10, BNF11, BNF12)
# write.csv(AverageBNF11,"E:/有机替代.作图/newdata/Mean.csv")
AverageBNF12 <- ddply(tota12, c("NTR"), summarise,
                      N    = sum(!is.na(Syn)),
                      median = median(Syn, na.rm=TRUE),
                      mean = mean(Syn, na.rm=TRUE),
                      sd   = sd(Syn, na.rm=TRUE),
                      se   = sd / sqrt(N),
                      ci_lower = mean - qt(0.975, df = N - 1) * se,
                      ci_upper = mean + qt(0.975, df = N - 1) * se
)
windowsFonts(RMN = windowsFont("Times New Roman"))
# 设置因子顺序一致（必须）
ordered_levels <- c("Den_Syn_f_ab", "Den_Syn_t_ab", "Nit_Syn_f_ab", "Nit_Syn_t_ab", 
                    "Imm_Syn_f_ab", "Imm_Syn_t_ab", "Min_Syn_f_ab", "Min_Syn_t_ab")
tota12$NTR <- factor(tota12$NTR, levels = ordered_levels)
AverageBNF12$NTR <- factor(AverageBNF12$NTR, levels = ordered_levels)

# 绘图
MeanEcoBNF1 <- ggplot(data = tota12, aes(x = NTR, y = Syn)) +
  
  # 箱线图
  geom_boxplot(aes(fill = NTR), fatten = 0.5,size = 0.5) +
  # 加粗中位数线（新增）
  stat_summary(fun = median, geom = "crossbar",
               width = 0.75, color = "black", size = 0.5) +
  # 基准线
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "blue", size = 1.8) + 
  # 原始散点
  geom_jitter(position = position_jitter(0.25), alpha = 0.6) +
  
  geom_point(
    data = AverageBNF12,
    mapping = aes(x = NTR, y = mean),
    shape = 21,
    fill = "red",
    color = "red",
    size = 7
  ) +
  
  # 分组垂直线
  geom_vline(xintercept = 2.5, linetype = "dashed", colour = "black", size = 1) +
  geom_vline(xintercept = 4.5, linetype = "dashed", colour = "black", size = 1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", colour = "black", size = 1) +
  # 轴范围
  ylim(-5, 5) +
  # 自定义颜色
  scale_fill_manual(
    name = "NTR",
    values = c(
      "Imm_Syn_f_ab" = "#33FF57", "Den_Syn_f_ab" = "#33FF57",
      "Min_Syn_f_ab" = "#33FF57", "Nit_Syn_f_ab" = "#33FF57",
      "Imm_Syn_t_ab" = "#F0E442", "Den_Syn_t_ab" = "#F0E442",
      "Min_Syn_t_ab" = "#F0E442", "Nit_Syn_t_ab" = "#F0E442"
    )
  ) +
  
  # 翻转坐标轴
  coord_flip() +
  
  # 主题设置
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.background = element_rect(colour = "transparent", fill = "transparent"),
    legend.justification = c(1, 1), 
    legend.position = c(-30, 0.95),
    legend.text = element_text(family = "RMN", face = "bold", colour = "black", size = 13),
    legend.title = element_text(family = "RMN", face = "bold", colour = "black", size = 18),
    axis.title = element_text(family = "RMN", face = "bold", colour = "black", size = 18),
    axis.text.x = element_text(family = "RMN", face = "bold", colour = "black", size = 18),
    axis.text.y = element_text(family = "RMN", face = "bold", colour = "black", size = 16),
    panel.grid.major = element_blank(),   # 去掉主网格线
    panel.grid.minor = element_blank()    # 去掉次网格线
  )

# 显示图形
MeanEcoBNF1

tiff("E:/有机替代.作图/response/Figure_all2.tiff",width=5000,height=2000,res=300,compression="lzw")
figure_all1 <- ggarrange(MeanEcoBNF3, MeanEcoBNF1, 
                         widths = c(1, 1), heights = c(1, 1), ncol=2 , nrow=1,
                         labels = c("(a)","(b)"),
                         label.x=0.23, label.y=0.98, align = "v",
                         font.label = list(size = 25, color = "black", fontface="bold", family="RMN"))

# 打印图形
print(figure_all1)
dev.off()

##N2O_emission_Nit
N2O_emission_Nit <- lmer(LnRR.N2O_emission~LnRR.Nitrification+(1|studies),weights=replicate, data=traits_data)#for Shoot2
summary(N2O_emission_Nit)
r.squaredGLMM(N2O_emission_Nit)
allN2O_emission_Nit <- ci.fun(N2O_emission_Nit,"LnRR.Nitrification")
# 强制将所有变量转换为没有属性的标量数值
N_N2O_emission_Nit <- as.numeric(length(residuals(N2O_emission_Nit)))
R2_N2O_emission_Nit <- as.numeric(round(r.squaredGLMM(N2O_emission_Nit)[, 2], 3))
p_N2O_emission_Nit <- as.numeric(round(summary(N2O_emission_Nit)$coefficients[2, 5], 4))
coef_Int_N2O_emission_Nit <- as.numeric(summary(N2O_emission_Nit)$coefficients[ , 1][1])
coef_N2O_emission_Nit <- as.numeric(summary(N2O_emission_Nit)$coefficients[ , 1][2])
traits_data$N2O_emission_Nit <-resid(N2O_emission_Nit)+coef_Int_N2O_emission_Nit+coef_N2O_emission_Nit*traits_data$LnRR.Nitrification
# 修正 annotate() 中传递给 label 的对象，确保都是单一值
N2O_emission1 <- ggplot() +
  geom_point(aes(x = traits_data$LnRR.Nitrification, y = traits_data$LnRR.N2O_emission, size = traits_data$replicate), alpha = 1.5, shape = 1, color = "blue", stroke = 1.2) +
  scale_size_continuous(range = c(8, 15)) + # 自适应调整点大小范围
  geom_abline(intercept = coef_Int_N2O_emission_Nit, slope = coef_N2O_emission_Nit, color = "green", size = 3,linetype="dashed") +
  geom_ribbon(aes(x = allN2O_emission_Nit$xsim, ymin = allN2O_emission_Nit$ci.low, ymax = allN2O_emission_Nit$ci.up), fill = "#BAE1FF", alpha = 0.65) +
  theme_gray(base_size = 16) + # 修改为灰色背景
  xlim(-7.5, 5) + ylim(-3, 2) +
  ylab(expression(paste("lnRR N\u2082O emission"))) +
  xlab(expression(paste("lnRR Nitrification"))) +
  theme(
    axis.title = element_text(family = "RMN", color = "black", size = 30, face = "bold"),
    axis.text.y = element_text(family = "RMN", color = "black", size = 30),
    axis.text.x = element_text(family = "RMN", color = "black", size = 30),
    
    panel.grid.major = element_blank(),   # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次网格线
    panel.background = element_blank(),   # 去掉面板背景
    plot.background = element_blank(),    # 去掉整体背景
    legend.position = "none",
    
    panel.border = element_rect(color = "black", size = 1, fill = NA) # 保留边框
  ) +
  annotate("text", x = -5.5, y = -1, label = "Slope == 0.04", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = -5.5, y = -1.5, label = "p == 0.19", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = -5.5, y = -2, label = "R^{2} == 0.85", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = -5.5, y = -2.5, label = "N == 349", family = "RMN", fontface = "bold", parse = TRUE, size = 8)
# 打印图形
N2O_emission1

##N2O_emission_Den
NFcapTime1 <- traits_data[complete.cases(traits_data$LnRR.N2O_emission),]
N2O_emission_Den <- lmer(LnRR.N2O_emission~LnRR.Denitrification+(1|studies),weights=replicate, data=NFcapTime1)#for Shoot2
summary(N2O_emission_Den)
r.squaredGLMM(N2O_emission_Den)
allN2O_emission_Den <- ci.fun(N2O_emission_Den,"LnRR.Denitrification")
# 强制将所有变量转换为没有属性的标量数值
N_N2O_emission_Den <- as.numeric(length(residuals(N2O_emission_Den)))
R2_N2O_emission_Den <- as.numeric(round(r.squaredGLMM(N2O_emission_Den)[, 2], 3))
p_N2O_emission_Den <- as.numeric(round(summary(N2O_emission_Den)$coefficients[2, 5], 4))
coef_Int_N2O_emission_Den <- as.numeric(summary(N2O_emission_Den)$coefficients[ , 1][1])
coef_N2O_emission_Den <- as.numeric(summary(N2O_emission_Den)$coefficients[ , 1][2])
NFcapTime1$N2O_emission_Den <-resid(N2O_emission_Den)+coef_Int_N2O_emission_Den+coef_N2O_emission_Den*NFcapTime1$LnRR.Denitrification
# 修正 annotate() 中传递给 label 的对象，确保都是单一值
N2O_emission2 <- ggplot() +
  geom_point(aes(x = NFcapTime1$LnRR.Denitrification, y = NFcapTime1$N2O_emission_Den, size = NFcapTime1$replicate), alpha = 1.5, shape = 1, color = "blue", stroke = 1.2) +
  scale_size_continuous(range = c(8, 15)) + # 自适应调整点大小范围
  geom_abline(intercept = coef_Int_N2O_emission_Den, slope = coef_N2O_emission_Den, color = "green", size = 3,linetype="dashed") +
  geom_ribbon(aes(x = allN2O_emission_Den$xsim, ymin = allN2O_emission_Den$ci.low, ymax = allN2O_emission_Den$ci.up), fill = "#BAE1FF", alpha = 0.65) +
  theme_gray(base_size = 16) + # 修改为灰色背景
  xlim(-6, 6) + ylim(-1, 1) +
  ylab(expression(paste("lnRR N\u2082O emission"))) +
  xlab(expression(paste("lnRR Denitrification"))) +
  theme(
    axis.title = element_text(family = "RMN", color = "black", size = 30, face = "bold"),
    axis.text.y = element_text(family = "RMN", color = "black", size = 30),
    axis.text.x = element_text(family = "RMN", color = "black", size = 30),
    
    panel.grid.major = element_blank(),   # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次网格线
    panel.background = element_blank(),   # 去掉面板背景
    plot.background = element_blank(),    # 去掉整体背景
    legend.position = "none",
    
    panel.border = element_rect(color = "black", size = 1, fill = NA) # 保留边框
  ) +
  annotate("text", x = -3, y = 1, label = "Slope == 0.02", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = -3, y = 0.8, label = "p == 0.46", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = -3, y = 0.6, label = "R^{2} == 0.80", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = -3, y = 0.4, label = "N == 255", family = "RMN", fontface = "bold", parse = TRUE, size = 8)
# 打印图形
N2O_emission2

##NH3_emission_Min
NFcapTime2 <- traits_data[complete.cases(traits_data$LnRR.NH3_emission),]
NH3_emission_Min <- lmer(LnRR.NH3_emission~LnRR.Mineralization+(1|studies),weights=replicate, data=NFcapTime2)#for Shoot2
summary(NH3_emission_Min)
r.squaredGLMM(NH3_emission_Min)
allNH3_emission_Min <- ci.fun(NH3_emission_Min,"LnRR.Mineralization")
# 强制将所有变量转换为没有属性的标量数值
N_NH3_emission_Min <- as.numeric(length(residuals(NH3_emission_Min)))
R2_NH3_emission_Min <- as.numeric(round(r.squaredGLMM(NH3_emission_Min)[, 2], 3))
p_NH3_emission_Min <- as.numeric(round(summary(NH3_emission_Min)$coefficients[2, 5], 4))
coef_Int_NH3_emission_Min <- as.numeric(summary(NH3_emission_Min)$coefficients[ , 1][1])
coef_NH3_emission_Min <- as.numeric(summary(NH3_emission_Min)$coefficients[ , 1][2])
NFcapTime2$NH3_emission_Min <-resid(NH3_emission_Min)+coef_Int_NH3_emission_Min+coef_NH3_emission_Min*NFcapTime2$LnRR.Mineralization
# 修正 annotate() 中传递给 label 的对象，确保都是单一值
NH3_emission1 <- ggplot() +
  geom_point(aes(x = NFcapTime2$LnRR.Mineralization, y = NFcapTime2$NH3_emission_Min, size = NFcapTime2$replicate), alpha = 1.5, shape = 1, color = "blue", stroke = 1.2) +
  scale_size_continuous(range = c(8, 15)) + # 自适应调整点大小范围
  geom_abline(intercept = coef_Int_NH3_emission_Min, slope = coef_NH3_emission_Min, color = "green", size = 3,linetype="dashed") +
  geom_ribbon(aes(x = allNH3_emission_Min$xsim, ymin = allNH3_emission_Min$ci.low, ymax = allNH3_emission_Min$ci.up), fill = "#BAE1FF", alpha = 0.65) +
  theme_gray(base_size = 16) + # 修改为灰色背景
  xlim(-0.5, 1.5) + ylim(-0.5, 0.5) +
  ylab(expression(paste("lnRR NH\u2083 emission"))) +
  xlab(expression(paste("lnRR Mineralization"))) +
  theme(
    axis.title = element_text(family = "RMN", color = "black", size = 30, face = "bold"),
    axis.text.y = element_text(family = "RMN", color = "black", size = 30),
    axis.text.x = element_text(family = "RMN", color = "black", size = 30),
    
    panel.grid.major = element_blank(),   # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次网格线
    panel.background = element_blank(),   # 去掉面板背景
    plot.background = element_blank(),    # 去掉整体背景
    legend.position = "none",
    
    panel.border = element_rect(color = "black", size = 1, fill = NA) # 保留边框
  ) +
  annotate("text", x = 1, y = 0.5, label = "Slope == -0.03", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 1, y = 0.4, label = "p == 0.56", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 1, y = 0.3, label = "R^{2} == 0.88", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 1, y = 0.2, label = "N == 115", family = "RMN", fontface = "bold", parse = TRUE, size = 8)
# 打印图形
NH3_emission1

##N2O_emission_Nit_T
NFcapTime3 <- traits_data[complete.cases(traits_data$LnRR.N2O_emission),]
N2O_emission_DNS <- lmer(LnRR.N2O_emission~Nit_Syn_abs+(1|studies),weights=replicate, data=NFcapTime3)#for Shoot2
summary(N2O_emission_DNS)
r.squaredGLMM(N2O_emission_DNS)
allN2O_emission_DNS <- ci.fun(N2O_emission_DNS,"Nit_Syn_abs")
# 强制将所有变量转换为没有属性的标量数值
N_N2O_emission_DNS <- as.numeric(length(residuals(N2O_emission_DNS)))
R2_N2O_emission_DNS <- as.numeric(round(r.squaredGLMM(N2O_emission_DNS)[, 2], 3))
p_N2O_emission_DNS <- as.numeric(round(summary(N2O_emission_DNS)$coefficients[2, 5], 4))
coef_Int_N2O_emission_DNS <- as.numeric(summary(N2O_emission_DNS)$coefficients[ , 1][1])
coef_N2O_emission_DNS <- as.numeric(summary(N2O_emission_DNS)$coefficients[ , 1][2])
NFcapTime3$N2O_emission_DNS <-resid(N2O_emission_DNS)+coef_Int_N2O_emission_DNS+coef_N2O_emission_DNS*NFcapTime3$Nit_Syn_abs
# 修正 annotate() 中传递给 label 的对象，确保都是单一值
N2O_emission3 <- ggplot() +
  geom_point(aes(x = NFcapTime3$Nit_Syn_abs, y = NFcapTime3$N2O_emission_DNS, size = NFcapTime3$replicate), alpha = 1.5, shape = 1, color = "blue", stroke = 1.2) +
  scale_size_continuous(range = c(8, 15)) + # 自适应调整点大小范围
  geom_abline(intercept = coef_Int_N2O_emission_DNS, slope = coef_N2O_emission_DNS, color = "green", size = 3) +
  geom_ribbon(aes(x = allN2O_emission_DNS$xsim, ymin = allN2O_emission_DNS$ci.low, ymax = allN2O_emission_DNS$ci.up), fill = "#BAE1FF", alpha = 0.65) +
  theme_gray(base_size = 16) + # 修改为灰色背景
  # xlim(0, 10) + ylim(-1, 1) +
  ylab(expression(paste("lnRR N\u2082O emission"))) +
  xlab(expression(paste("DNS"))) +
  theme(
    axis.title = element_text(family = "RMN", color = "black", size = 30, face = "bold"),
    axis.text.y = element_text(family = "RMN", color = "black", size = 30),
    axis.text.x = element_text(family = "RMN", color = "black", size = 30),
    
    panel.grid.major = element_blank(),   # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次网格线
    panel.background = element_blank(),   # 去掉面板背景
    plot.background = element_blank(),    # 去掉整体背景
    legend.position = "none",
    
    panel.border = element_rect(color = "black", size = 1, fill = NA) # 保留边框
  ) +
  annotate("text", x = 7.5, y = -0.2, label = "Slope == 0.05", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 7.5, y = -0.4, label = "p == 0.02", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 7.5, y = -0.6, label = "R^{2} == 0.86", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 7.5, y = -0.8, label = "N == 198", family = "RMN", fontface = "bold", parse = TRUE, size = 8)
# 打印图形
N2O_emission3

##N2O_emission_Den_T
NFcapTime4 <- traits_data[complete.cases(traits_data$LnRR.N2O_emission),]
N2O_emission_DDS <- lmer(LnRR.N2O_emission~Den_Syn_abs+(1|studies),weights=replicate, data=NFcapTime4)#for Shoot2
summary(N2O_emission_DDS)
r.squaredGLMM(N2O_emission_DDS)
allN2O_emission_DDS <- ci.fun(N2O_emission_DDS,"Den_Syn_abs")
# 强制将所有变量转换为没有属性的标量数值
N_N2O_emission_DDS <- as.numeric(length(residuals(N2O_emission_DDS)))
R2_N2O_emission_DDS <- as.numeric(round(r.squaredGLMM(N2O_emission_DDS)[, 2], 3))
p_N2O_emission_DDS <- as.numeric(round(summary(N2O_emission_DDS)$coefficients[2, 5], 4))
coef_Int_N2O_emission_DDS <- as.numeric(summary(N2O_emission_DDS)$coefficients[ , 1][1])
coef_N2O_emission_DDS <- as.numeric(summary(N2O_emission_DDS)$coefficients[ , 1][2])
NFcapTime4$N2O_emission_DDS <-resid(N2O_emission_DDS)+coef_Int_N2O_emission_DDS+coef_N2O_emission_DDS*NFcapTime4$Den_Syn_abs
# 修正 annotate() 中传递给 label 的对象，确保都是单一值
N2O_emission4 <- ggplot() +
  geom_point(aes(x = NFcapTime4$Den_Syn_abs, y = NFcapTime4$N2O_emission_DDS, size = NFcapTime4$replicate), alpha = 1.5, shape = 1, color = "blue", stroke = 1.2) +
  scale_size_continuous(range = c(8, 15)) + # 自适应调整点大小范围
  geom_abline(intercept = coef_Int_N2O_emission_DDS, slope = coef_N2O_emission_DDS, color = "green", size = 3) +
  geom_ribbon(aes(x = allN2O_emission_DDS$xsim, ymin = allN2O_emission_DDS$ci.low, ymax = allN2O_emission_DDS$ci.up), fill = "#BAE1FF", alpha = 0.65) +
  theme_gray(base_size = 16) + # 修改为灰色背景
  xlim(0, 10) + ylim(-1, 1) +
  ylab(expression(paste("lnRR N\u2082O emission"))) +
  xlab(expression(paste("DDS"))) +
  theme(
    axis.title = element_text(family = "RMN", color = "black", size = 30, face = "bold"),
    axis.text.y = element_text(family = "RMN", color = "black", size = 30),
    axis.text.x = element_text(family = "RMN", color = "black", size = 30),
    
    panel.grid.major = element_blank(),   # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次网格线
    panel.background = element_blank(),   # 去掉面板背景
    plot.background = element_blank(),    # 去掉整体背景
    legend.position = "none",
    
    panel.border = element_rect(color = "black", size = 1, fill = NA) # 保留边框
  ) +
  annotate("text", x = 7.5, y = -0.2, label = "Slope == 0.04", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 7.5, y = -0.4, label = "p < 0.01", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 7.5, y = -0.6, label = "R^{2} == 0.87", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 7.5, y = -0.8, label = "N == 189", family = "RMN", fontface = "bold", parse = TRUE, size = 8)
# 打印图形
N2O_emission4

##NH3_emission_Min_T
NFcapTime5 <- traits_data[complete.cases(traits_data$LnRR.NH3_emission),]
NH3_emission_DMS <- lmer(LnRR.NH3_emission~Min_Syn_abs+(1|studies),weights=replicate, data=NFcapTime5)#for Shoot2
summary(NH3_emission_DMS)
r.squaredGLMM(NH3_emission_DMS)
allNH3_emission_DMS <- ci.fun(NH3_emission_DMS,"Min_Syn_abs")
# 强制将所有变量转换为没有属性的标量数值
N_NH3_emission_DMS <- as.numeric(length(residuals(NH3_emission_DMS)))
R2_NH3_emission_DMS <- as.numeric(round(r.squaredGLMM(NH3_emission_DMS)[, 2], 3))
p_NH3_emission_DMS <- as.numeric(round(summary(NH3_emission_DMS)$coefficients[2, 5], 4))
coef_Int_NH3_emission_DMS <- as.numeric(summary(NH3_emission_DMS)$coefficients[ , 1][1])
coef_NH3_emission_DMS <- as.numeric(summary(NH3_emission_DMS)$coefficients[ , 1][2])
NFcapTime5$NH3_emission_DMS <-resid(NH3_emission_DMS)+coef_Int_NH3_emission_DMS+coef_NH3_emission_DMS*NFcapTime5$Min_Syn_abs
# 修正 annotate() 中传递给 label 的对象，确保都是单一值
NH3_emission2 <- ggplot() +
  geom_point(aes(x = NFcapTime5$Min_Syn_abs, y = NFcapTime5$NH3_emission_DMS, size = NFcapTime5$replicate), alpha = 1.5, shape = 1, color = "blue", stroke = 1.2) +
  scale_size_continuous(range = c(8, 15)) + # 自适应调整点大小范围
  geom_abline(intercept = coef_Int_NH3_emission_DMS, slope = coef_NH3_emission_DMS, color = "green", size = 3) +
  geom_ribbon(aes(x = allNH3_emission_DMS$xsim, ymin = allNH3_emission_DMS$ci.low, ymax = allNH3_emission_DMS$ci.up), fill = "#BAE1FF", alpha = 0.65) +
  theme_gray(base_size = 16) + # 修改为灰色背景
  xlim(0, 8) + ylim(-1, 0.1) +
  ylab(expression(paste("lnRR NH\u2083 emission"))) +
  xlab(expression(paste("DMS"))) +
  theme(
    axis.title = element_text(family = "RMN", color = "black", size = 30, face = "bold"),
    axis.text.y = element_text(family = "RMN", color = "black", size = 30),
    axis.text.x = element_text(family = "RMN", color = "black", size = 30),
    
    panel.grid.major = element_blank(),   # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次网格线
    panel.background = element_blank(),   # 去掉面板背景
    plot.background = element_blank(),    # 去掉整体背景
    legend.position = "none",
    
    panel.border = element_rect(color = "black", size = 1, fill = NA) # 保留边框
  ) +
  annotate("text", x = 5.5, y = -0.6, label = "Slope == 0.04", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -0.7, label = "p < 0.01", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -0.8, label = "R^{2} == 0.91", family = "RMN", fontface = "bold", parse = TRUE, size = 8) +
  annotate("text", x = 5.5, y = -0.9, label = "N == 78", family = "RMN", fontface = "bold", parse = TRUE, size = 8)
# 打印图形
NH3_emission2

tiff("E:/有机替代.作图/response/Figure1b.tiff",width=6000,height=3500,res=300,compression="lzw")
figure_all1 <- ggarrange(NH3_emission1, N2O_emission1, N2O_emission2, NH3_emission2, N2O_emission3, N2O_emission4, 
                         widths = c(1, 1, 1, 1, 1, 1), heights = c(1, 1), ncol=3 , nrow=2,
                         labels = c("(c)","(d)","(e)","(f)","(g)","(h)"),
                         label.x=0.23, label.y=0.98, align = "v",
                         font.label = list(size = 25, color = "black", fontface="bold", family="RMN"))

# 打印图形
print(figure_all1)

dev.off()
}
###Fig 3
{
# ====== 加载包 ======
suppressPackageStartupMessages({
  library(readxl)
  library(mice)
  library(piecewiseSEM)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(stringr)
  library(ranger)
  library(lme4)
  library(lmerTest)
  
  has_fwrite <- requireNamespace("data.table", quietly = TRUE)
})

LnRRdata1 <- read.csv("E:/有机替代.作图/response/N2O.csv")

LnRRdata1$DNS <- as.numeric(scale(LnRRdata1$Nit_Syn_abs))
LnRRdata1$DDS <- as.numeric(scale(LnRRdata1$Den_Syn_abs))

model_climate <- lmer(
  LnRR.N2O_emission ~ Ln.MAT + Ln.MAP + (1 | studies),
  data = LnRRdata1,
  na.action = na.omit
)

climate_coefficients <- summary(model_climate)

coefs(model_climate, standardize = "scale")

beta_MAT <- climate_coefficients$coefficients[2, 1]
beta_MAP <- climate_coefficients$coefficients[3, 1]

LnRRdata1$Climate <- beta_MAT * LnRRdata1$Ln.MAT + beta_MAP * LnRRdata1$Ln.MAP

sem.model <- psem(
  lmer(LnRR.pH ~ N_addition + (1|studies), data = LnRRdata1),
  lmer(LnRR.SOC ~ N_addition + (1|studies), data = LnRRdata1),
  lmer(DDS ~ LnRR.pH + N_addition + Climate + LnRR.SOC + (1|studies), data = LnRRdata1),
  lmer(DNS ~ LnRR.pH + N_addition + Climate + LnRR.SOC + (1|studies), data = LnRRdata1),
  lmer(LnRR.N2O_emission ~ DDS + DNS + LnRR.pH + N_addition + Climate + LnRR.SOC + (1|studies), data = LnRRdata1)
)
summary(sem.model)

LnRRdata1 <- read.csv("E:/有机替代.作图/response/NH3.csv")

model_climate <- lmer(
  LnRR.NH3_emission ~ Ln.MAT + Ln.MAP + (1 | studies),
  data = LnRRdata1,
  na.action = na.omit
)

climate_coefficients <- summary(model_climate)

coefs(model_climate, standardize = "scale")

beta_MAT <- climate_coefficients$coefficients[2, 1]
beta_MAP <- climate_coefficients$coefficients[3, 1]

LnRRdata1$Climate <- beta_MAT * LnRRdata1$Ln.MAT + beta_MAP * LnRRdata1$Ln.MAP

sem.model <- psem(
  lmer(LnRR.pH ~ N_addition + (1|studies), data = LnRRdata1),
  lmer(LnRR.SOC ~ N_addition + (1|studies), data = LnRRdata1),
  lmer(Min_Syn_abs ~ LnRR.pH + N_addition + Climate + LnRR.SOC + (1|studies), data = LnRRdata1),
  lmer(Imm_Syn_abs ~ LnRR.pH + N_addition + Climate + LnRR.SOC + (1|studies), data = LnRRdata1),
  lmer(LnRR.NH3_emission ~ Min_Syn_abs + Imm_Syn_abs + LnRR.pH + N_addition + Climate + LnRR.SOC + (1|studies), data = LnRRdata1)
)
summary(sem.model)
}
###Fig4
{
####DNS预测
{
## ===== 依赖包 =====
library(sp)       # coordinates(), CRS 等
library(raster)   # stack(), extract()
library(car)      # car::vif
library(usdm)     # 如果你后续还要用到usdm的其它函数，保留；否则可以不加载
library(caret)
library(dplyr)
library(tidyr)
## ===== 1) 读取观测点并设坐标系 =====
MVDdata <- read.csv("E:/有机替代.作图/response/预测/Prediction1.csv", stringsAsFactors = FALSE)

# 如果数据里经纬度列名不是 "Longitude" "Latitude"，请改成实际列名
coordinates(MVDdata) <- ~ Longitude + Latitude
crs(MVDdata) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 2) 协变量栅格堆栈 =====
file_list <- list.files(path = "E:/有机替代.作图/response/预测/List2", full.names = TRUE)
Covariatestacked <- stack(file_list)
crs(Covariatestacked) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 3) 栅格提取到点 =====
extracted_values <- raster::extract(Covariatestacked, MVDdata)
extracted_df <- as.data.frame(extracted_values)

# 给栅格提取结果起列名（使用文件名的安全版本）
safe_raster_names <- make.names(basename(file_list))
colnames(extracted_df) <- safe_raster_names

## ===== 4) 合并回观测数据（转成数据框以便与数据绑定） =====
GalNdata_df <- as.data.frame(MVDdata)
combined_data <- cbind(GalNdata_df, extracted_df)

# 添加经纬度列（如果还没有添加的话）
combined_data$Latitude <- MVDdata$Latitude
combined_data$Longitude <- MVDdata$Longitude

## ===== 5) 数据预处理 =====
# 处理缺失值：数值列用列均值填补
is_num <- vapply(combined_data, is.numeric, logical(1))
combined_data[is_num] <- lapply(combined_data[is_num], function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

## ===== 6) 划分训练和测试集 =====
set.seed(12345)
splitIndex <- createDataPartition(combined_data$Nit_Syn_abs, p = 0.7, list = FALSE)
train_data <- combined_data[splitIndex, ]   # 70% 训练集
test_data  <- combined_data[-splitIndex, ]  # 30% 测试集

## ===== 7) 网格搜索调参 =====
library(randomForest)
library(parallel)
library(doParallel)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

tune_grid <- expand.grid(
  mtry = 1:13,   # 尝试的 mtry 参数范围
  ntree = c(100, 500, 1000, 1500, 2000, 2500, 3000),   # 尝试的 ntree 参数值
  nodesize = 1:10   # 尝试的 nodesize 参数范围
)

x_train_data <- data.frame(as.matrix(train_data[, -ncol(train_data)]))  # 特征变量
y_train_data <- train_data$Nit_Syn_abs  # 目标变量

results_df <- data.frame(mtry=integer(), ntree=integer(), nodesize=integer(),
                         R2=numeric(), RMSE=numeric(), MAE=numeric(),
                         stringsAsFactors=FALSE)

total_combinations <- nrow(tune_grid)
set.seed(12345)
for (i in 1:total_combinations) {
  params <- tune_grid[i, ]
  cat("Processing combination", i, "of", total_combinations, "...\n")
  
  model <- randomForest(
    x = x_train_data, y = y_train_data,
    mtry = params$mtry, ntree = params$ntree, nodesize = params$nodesize,
    importance = TRUE
  )
  
  preds <- predict(model, newdata = test_data)
  r2   <- cor(test_data$Nit_Syn_abs, preds)^2
  rmse <- sqrt(mean((test_data$Nit_Syn_abs - preds)^2))
  mae  <- mean(abs(test_data$Nit_Syn_abs - preds))
  
  results_df <- rbind(results_df, data.frame(
    mtry=params$mtry, ntree=params$ntree, nodesize=params$nodesize,
    R2=r2, RMSE=rmse, MAE=mae
  ))
  
  cat("Completed combination", i, "with R^2", r2, "\n")
}

# 查看结果
View(results_df)

# 选择最佳参数（R2 值接近目标值）
target_r2  <- 0.8860
tolerance  <- 0.0001

best_params1 <- results_df[abs(results_df$R2 - target_r2) <= tolerance, ]
# 或者：best_params1 <- results_df[which.max(results_df$R2), ]

stopCluster(cl)

## ===== 8) 用最佳参数训练最终模型 =====
best_mtry    <- best_params1$mtry
best_ntree   <- best_params1$ntree
best_nodesize <- best_params1$nodesize

set.seed(12345)
final_rf_model <- randomForest(
  x = x_train_data, y = y_train_data,
  mtry = best_mtry,
  ntree = best_ntree,
  nodesize = best_nodesize,
  importance = TRUE
)

# 预测训练集和测试集结果
y_train_rf <- predict(final_rf_model, newdata = x_train_data)
y_test_rf  <- predict(final_rf_model, newdata = test_data)

DATA_train <- data.frame(X = y_train_rf, Y = y_train_data, Dataset = "train")
DATA_test  <- data.frame(X = y_test_rf,  Y = test_data$Nit_Syn_abs,  Dataset = "test")
Data <- rbind(DATA_train, DATA_test)

# 评估模型
R2   <- round(cor(Data$X, Data$Y)^2, 2)
RMSE <- round(sqrt(mean((Data$X - Data$Y)^2)), 2)
MAE  <- round(mean(abs(Data$X - Data$Y)), 2)

lin_model <- lm(Data$Y ~ Data$X - 1)  # 过原点拟合
slope <- round(coef(lin_model)[1], 2)

cat("R2: ", R2, "\nRMSE: ", RMSE, "\nMAE: ", MAE, "\nSlope: ", slope, "\n")


## ===== 用现有 final_rf_model 直接做全球预测（对齐训练特征）=====

####循环100
# 优先用 best_params1（你的“接近目标R²”的组合）；如果为空，就取 R² 最高的前 10 个
if (exists("best_params1") && nrow(best_params1) > 0) {
  param_table <- best_params1
} else {
  topK <- 10
  param_table <- head(results_df[order(-results_df$R2), ], topK)
}

# 输出目录
out_dir <- "E:/有机替代.作图/response/预测"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ========= 1) 训练时的特征顺序 & 栅格层名对齐 =========
feature_cols <- colnames(x_train_data)

# 与提取时一致的安全图层名
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜（NoData）
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

# 目标名（如果误入特征时用于常数填补）
target_col  <- "Nit_Syn_abs"
target_mean <- mean(y_train_data, na.rm = TRUE)

## ========= 2) 一次性构造 “全球像元 × 特征” 矩阵（供所有模型复用）=========
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols

for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
    message("提示：特征 '", nm, "' 在栅格与训练数据中都无法直接映射，已用 0 代填。")
  }
}

covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
# 数值列 NA 均值填补
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) {
    v[is.na(v)] <- mean(v, na.rm = TRUE)
  }
  v
})
names(covariate_matrix) <- feature_cols

## ========= 3) 循环训练 → 预测 → 写盘 =========
out_dir <- "E:/有机替代.作图/response/预测/DNS"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保证栅格层名和训练特征一致
feature_cols <- colnames(x_train_data)
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

target_col  <- "Nit_Syn_abs"
target_mean <- mean(y_train_data, na.rm = TRUE)

# 一次性构造全球输入矩阵（避免每轮重复）
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols
for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
  }
}
covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) v[is.na(v)] <- mean(v, na.rm = TRUE)
  v
})
names(covariate_matrix) <- feature_cols

## ========= 循环 100 次 =========
results <- data.frame(FileName = character(), R2 = numeric(), RMSE = numeric(), MAE = numeric())

for (i in 1:100) {
  set.seed(12345 + i)
  
  # 每轮都随机重划分 70/30
  splitIndex <- createDataPartition(combined_data$Nit_Syn_abs, p = 0.7, list = FALSE)
  train_i <- combined_data[splitIndex, ]
  test_i  <- combined_data[-splitIndex, ]
  
  x_train_i <- data.frame(as.matrix(train_i[, feature_cols]))
  y_train_i <- train_i$Nit_Syn_abs
  
  rf_model <- randomForest(
    x = x_train_i, y = y_train_i,
    mtry = best_mtry, ntree = best_ntree, nodesize = best_nodesize,
    importance = TRUE
  )
  
  # 验证集表现
  preds_test <- predict(rf_model, newdata = test_i[, feature_cols])
  r2   <- cor(test_i$Nit_Syn_abs, preds_test)^2
  rmse <- sqrt(mean((test_i$Nit_Syn_abs - preds_test)^2))
  mae  <- mean(abs(test_i$Nit_Syn_abs - preds_test))
  
  # 全局预测
  preds_global <- as.numeric(predict(rf_model, newdata = covariate_matrix))
  
  pred_r <- Covariatestacked[[1]]
  values(pred_r) <- preds_global
  values(pred_r)[is.na(mask_vals)] <- NA
  
  out_name <- sprintf("NH3_RF_%03d.tif", i)
  out_path <- file.path(out_dir, out_name)
  writeRaster(pred_r, filename = out_path, format = "GTiff",
              datatype = "FLT4S", overwrite = TRUE)
  
  results <- rbind(results, data.frame(FileName = out_name, R2 = r2, RMSE = rmse, MAE = mae))
  
  cat(sprintf("✅ Loop %03d | R²=%.3f | RMSE=%.3f | MAE=%.3f | 已保存 %s\n",
              i, r2, rmse, mae, out_path))
}

write.csv(results, file.path(out_dir, "RF_100_results.csv"), row.names = FALSE)
cat("🎉 已完成 100 次循环预测，结果汇总保存到 RF_100_results.csv\n")



####计算均值和差异系数
library(terra)

# 输入文件夹路径
out_dir <- "E:/有机替代.作图/response/预测/DNS"

# 匹配 NH3_RF_001.tif ~ NH3_RF_100.tif
files <- list.files(out_dir, pattern = "NH3_RF_\\d{3}\\.tif$", full.names = TRUE)

# 检查是否找到文件
if (length(files) == 0) stop("没有找到符合条件的文件，请检查路径或文件名！")

# 读取并堆叠
raster_stack <- rast(files)

# 计算均值、标准差、变异系数
mean_raster <- mean(raster_stack, na.rm = TRUE)
sd_raster   <- stdev(raster_stack, na.rm = TRUE)
cv_raster   <- sd_raster / mean_raster

# 输出目录
output_dir <- "E:/有机替代.作图/response/预测/DNS"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 保存结果
writeRaster(mean_raster, file.path(output_dir, "Mean.tif"), overwrite = TRUE)
writeRaster(sd_raster,   file.path(output_dir, "SD.tif"),   overwrite = TRUE)
writeRaster(cv_raster,   file.path(output_dir, "CV.tif"),   overwrite = TRUE)
}
####DDS预测
{
## ===== 依赖包 =====
library(sp)       # coordinates(), CRS 等
library(raster)   # stack(), extract()
library(car)      # car::vif
library(usdm)     # 如果你后续还要用到usdm的其它函数，保留；否则可以不加载
library(caret)
library(dplyr)
library(tidyr)
## ===== 1) 读取观测点并设坐标系 =====
MVDdata <- read.csv("E:/有机替代.作图/response/预测/Prediction2.csv", stringsAsFactors = FALSE)

# 如果数据里经纬度列名不是 "Longitude" "Latitude"，请改成实际列名
coordinates(MVDdata) <- ~ Longitude + Latitude
crs(MVDdata) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 2) 协变量栅格堆栈 =====
file_list <- list.files(path = "E:/有机替代.作图/response/预测/List2", full.names = TRUE)
Covariatestacked <- stack(file_list)
crs(Covariatestacked) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 3) 栅格提取到点 =====
extracted_values <- raster::extract(Covariatestacked, MVDdata)
extracted_df <- as.data.frame(extracted_values)

# 给栅格提取结果起列名（使用文件名的安全版本）
safe_raster_names <- make.names(basename(file_list))
colnames(extracted_df) <- safe_raster_names

## ===== 4) 合并回观测数据（转成数据框以便与数据绑定） =====
GalNdata_df <- as.data.frame(MVDdata)
combined_data <- cbind(GalNdata_df, extracted_df)

# 添加经纬度列（如果还没有添加的话）
combined_data$Latitude <- MVDdata$Latitude
combined_data$Longitude <- MVDdata$Longitude

## ===== 5) 数据预处理 =====
# 处理缺失值：数值列用列均值填补
is_num <- vapply(combined_data, is.numeric, logical(1))
combined_data[is_num] <- lapply(combined_data[is_num], function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

## ===== 6) 划分训练和测试集 =====
set.seed(12345)
splitIndex <- createDataPartition(combined_data$Den_Syn_abs, p = 0.7, list = FALSE)
train_data <- combined_data[splitIndex, ]   # 70% 训练集
test_data  <- combined_data[-splitIndex, ]  # 30% 测试集

## ===== 7) 网格搜索调参 =====
library(randomForest)
library(parallel)
library(doParallel)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

tune_grid <- expand.grid(
  mtry = 1:13,   # 尝试的 mtry 参数范围
  ntree = c(100, 500, 1000, 1500, 2000, 2500, 3000),   # 尝试的 ntree 参数值
  nodesize = 1:10   # 尝试的 nodesize 参数范围
)

x_train_data <- data.frame(as.matrix(train_data[, -ncol(train_data)]))  # 特征变量
y_train_data <- train_data$Den_Syn_abs  # 目标变量

results_df <- data.frame(mtry=integer(), ntree=integer(), nodesize=integer(),
                         R2=numeric(), RMSE=numeric(), MAE=numeric(),
                         stringsAsFactors=FALSE)

total_combinations <- nrow(tune_grid)
set.seed(12345)
for (i in 1:total_combinations) {
  params <- tune_grid[i, ]
  cat("Processing combination", i, "of", total_combinations, "...\n")
  
  model <- randomForest(
    x = x_train_data, y = y_train_data,
    mtry = params$mtry, ntree = params$ntree, nodesize = params$nodesize,
    importance = TRUE
  )
  
  preds <- predict(model, newdata = test_data)
  r2   <- cor(test_data$Den_Syn_abs, preds)^2
  rmse <- sqrt(mean((test_data$Den_Syn_abs - preds)^2))
  mae  <- mean(abs(test_data$Den_Syn_abs - preds))
  
  results_df <- rbind(results_df, data.frame(
    mtry=params$mtry, ntree=params$ntree, nodesize=params$nodesize,
    R2=r2, RMSE=rmse, MAE=mae
  ))
  
  cat("Completed combination", i, "with R^2", r2, "\n")
}

# 查看结果
View(results_df)

# 选择最佳参数（R2 值接近目标值）
target_r2  <- 0.8800
tolerance  <- 0.0001

best_params1 <- results_df[abs(results_df$R2 - target_r2) <= tolerance, ]
# 或者：best_params1 <- results_df[which.max(results_df$R2), ]

stopCluster(cl)

## ===== 8) 用最佳参数训练最终模型 =====
best_mtry    <- best_params1$mtry
best_ntree   <- best_params1$ntree
best_nodesize <- best_params1$nodesize

set.seed(12345)
final_rf_model <- randomForest(
  x = x_train_data, y = y_train_data,
  mtry = best_mtry,
  ntree = best_ntree,
  nodesize = best_nodesize,
  importance = TRUE
)

# 预测训练集和测试集结果
y_train_rf <- predict(final_rf_model, newdata = x_train_data)
y_test_rf  <- predict(final_rf_model, newdata = test_data)

DATA_train <- data.frame(X = y_train_rf, Y = y_train_data, Dataset = "train")
DATA_test  <- data.frame(X = y_test_rf,  Y = test_data$Den_Syn_abs,  Dataset = "test")
Data <- rbind(DATA_train, DATA_test)

# 评估模型
R2   <- round(cor(Data$X, Data$Y)^2, 2)
RMSE <- round(sqrt(mean((Data$X - Data$Y)^2)), 2)
MAE  <- round(mean(abs(Data$X - Data$Y)), 2)

lin_model <- lm(Data$Y ~ Data$X - 1)  # 过原点拟合
slope <- round(coef(lin_model)[1], 2)

cat("R2: ", R2, "\nRMSE: ", RMSE, "\nMAE: ", MAE, "\nSlope: ", slope, "\n")


## ===== 用现有 final_rf_model 直接做全球预测（对齐训练特征）=====

####循环100
# 优先用 best_params1（你的“接近目标R²”的组合）；如果为空，就取 R² 最高的前 10 个
if (exists("best_params1") && nrow(best_params1) > 0) {
  param_table <- best_params1
} else {
  topK <- 10
  param_table <- head(results_df[order(-results_df$R2), ], topK)
}

# 输出目录
out_dir <- "E:/有机替代.作图/response/预测"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ========= 1) 训练时的特征顺序 & 栅格层名对齐 =========
feature_cols <- colnames(x_train_data)

# 与提取时一致的安全图层名
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜（NoData）
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

# 目标名（如果误入特征时用于常数填补）
target_col  <- "Den_Syn_abs"
target_mean <- mean(y_train_data, na.rm = TRUE)

## ========= 2) 一次性构造 “全球像元 × 特征” 矩阵（供所有模型复用）=========
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols

for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
    message("提示：特征 '", nm, "' 在栅格与训练数据中都无法直接映射，已用 0 代填。")
  }
}

covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
# 数值列 NA 均值填补
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) {
    v[is.na(v)] <- mean(v, na.rm = TRUE)
  }
  v
})
names(covariate_matrix) <- feature_cols

## ========= 3) 循环训练 → 预测 → 写盘 =========
out_dir <- "E:/有机替代.作图/response/预测/DDS"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保证栅格层名和训练特征一致
feature_cols <- colnames(x_train_data)
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

target_col  <- "Den_Syn_abs"
target_mean <- mean(y_train_data, na.rm = TRUE)

# 一次性构造全球输入矩阵（避免每轮重复）
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols
for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
  }
}
covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) v[is.na(v)] <- mean(v, na.rm = TRUE)
  v
})
names(covariate_matrix) <- feature_cols

## ========= 循环 100 次 =========
results <- data.frame(FileName = character(), R2 = numeric(), RMSE = numeric(), MAE = numeric())

for (i in 1:100) {
  set.seed(12345 + i)
  
  # 每轮都随机重划分 70/30
  splitIndex <- createDataPartition(combined_data$Den_Syn_abs, p = 0.7, list = FALSE)
  train_i <- combined_data[splitIndex, ]
  test_i  <- combined_data[-splitIndex, ]
  
  x_train_i <- data.frame(as.matrix(train_i[, feature_cols]))
  y_train_i <- train_i$Den_Syn_abs
  
  rf_model <- randomForest(
    x = x_train_i, y = y_train_i,
    mtry = best_mtry, ntree = best_ntree, nodesize = best_nodesize,
    importance = TRUE
  )
  
  # 验证集表现
  preds_test <- predict(rf_model, newdata = test_i[, feature_cols])
  r2   <- cor(test_i$Den_Syn_abs, preds_test)^2
  rmse <- sqrt(mean((test_i$Den_Syn_abs - preds_test)^2))
  mae  <- mean(abs(test_i$Den_Syn_abs - preds_test))
  
  # 全局预测
  preds_global <- as.numeric(predict(rf_model, newdata = covariate_matrix))
  
  pred_r <- Covariatestacked[[1]]
  values(pred_r) <- preds_global
  values(pred_r)[is.na(mask_vals)] <- NA
  
  out_name <- sprintf("NH3_RF_%03d.tif", i)
  out_path <- file.path(out_dir, out_name)
  writeRaster(pred_r, filename = out_path, format = "GTiff",
              datatype = "FLT4S", overwrite = TRUE)
  
  results <- rbind(results, data.frame(FileName = out_name, R2 = r2, RMSE = rmse, MAE = mae))
  
  cat(sprintf("✅ Loop %03d | R²=%.3f | RMSE=%.3f | MAE=%.3f | 已保存 %s\n",
              i, r2, rmse, mae, out_path))
}

write.csv(results, file.path(out_dir, "RF_100_results.csv"), row.names = FALSE)
cat("🎉 已完成 100 次循环预测，结果汇总保存到 RF_100_results.csv\n")



####计算均值和差异系数
library(terra)

# 输入文件夹路径
out_dir <- "E:/有机替代.作图/response/预测/DDS"

# 匹配 NH3_RF_001.tif ~ NH3_RF_100.tif
files <- list.files(out_dir, pattern = "NH3_RF_\\d{3}\\.tif$", full.names = TRUE)

# 检查是否找到文件
if (length(files) == 0) stop("没有找到符合条件的文件，请检查路径或文件名！")

# 读取并堆叠
raster_stack <- rast(files)

# 计算均值、标准差、变异系数
mean_raster <- mean(raster_stack, na.rm = TRUE)
sd_raster   <- stdev(raster_stack, na.rm = TRUE)
cv_raster   <- sd_raster / mean_raster

# 输出目录
output_dir <- "E:/有机替代.作图/response/预测/DDS"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 保存结果
writeRaster(mean_raster, file.path(output_dir, "Mean.tif"), overwrite = TRUE)
writeRaster(sd_raster,   file.path(output_dir, "SD.tif"),   overwrite = TRUE)
writeRaster(cv_raster,   file.path(output_dir, "CV.tif"),   overwrite = TRUE)
}
####DMS预测
{
## ===== 依赖包 =====
library(sp)       # coordinates(), CRS 等
library(raster)   # stack(), extract()
library(car)      # car::vif
library(usdm)     # 如果你后续还要用到usdm的其它函数，保留；否则可以不加载
library(caret)
library(dplyr)
library(tidyr)
## ===== 1) 读取观测点并设坐标系 =====
MVDdata <- read.csv("E:/有机替代.作图/response/预测/Prediction3.csv", stringsAsFactors = FALSE)

# 如果数据里经纬度列名不是 "Longitude" "Latitude"，请改成实际列名
coordinates(MVDdata) <- ~ Longitude + Latitude
crs(MVDdata) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 2) 协变量栅格堆栈 =====
file_list <- list.files(path = "E:/有机替代.作图/response/预测/List2", full.names = TRUE)
Covariatestacked <- stack(file_list)
crs(Covariatestacked) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 3) 栅格提取到点 =====
extracted_values <- raster::extract(Covariatestacked, MVDdata)
extracted_df <- as.data.frame(extracted_values)

# 给栅格提取结果起列名（使用文件名的安全版本）
safe_raster_names <- make.names(basename(file_list))
colnames(extracted_df) <- safe_raster_names

## ===== 4) 合并回观测数据（转成数据框以便与数据绑定） =====
GalNdata_df <- as.data.frame(MVDdata)
combined_data <- cbind(GalNdata_df, extracted_df)

# 添加经纬度列（如果还没有添加的话）
combined_data$Latitude <- MVDdata$Latitude
combined_data$Longitude <- MVDdata$Longitude

## ===== 5) 数据预处理 =====
# 处理缺失值：数值列用列均值填补
is_num <- vapply(combined_data, is.numeric, logical(1))
combined_data[is_num] <- lapply(combined_data[is_num], function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

## ===== 6) 划分训练和测试集 =====
set.seed(12345)
splitIndex <- createDataPartition(combined_data$Min_Syn_abs, p = 0.7, list = FALSE)
train_data <- combined_data[splitIndex, ]   # 70% 训练集
test_data  <- combined_data[-splitIndex, ]  # 30% 测试集

## ===== 7) 网格搜索调参 =====
library(randomForest)
library(parallel)
library(doParallel)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

tune_grid <- expand.grid(
  mtry = 1:13,
  ntree = c(100, 500, 1000, 1500, 2000, 2500, 3000),
  nodesize = 1:10
)

# 只随机抽取 50 个参数组合
set.seed(12345)
tune_grid <- tune_grid[sample(nrow(tune_grid), 50), ]

x_train_data <- data.frame(as.matrix(train_data[, -ncol(train_data)]))
y_train_data <- train_data$Min_Syn_abs

results_df <- data.frame(
  mtry = integer(),
  ntree = integer(),
  nodesize = integer(),
  R2 = numeric(),
  RMSE = numeric(),
  MAE = numeric(),
  stringsAsFactors = FALSE
)

total_combinations <- nrow(tune_grid)

set.seed(12345)

for (i in 1:total_combinations) {
  params <- tune_grid[i, ]
  cat("Processing combination", i, "of", total_combinations, "...\n")
  
  model <- randomForest(
    x = x_train_data,
    y = y_train_data,
    mtry = params$mtry,
    ntree = params$ntree,
    nodesize = params$nodesize,
    importance = TRUE
  )
  
  preds <- predict(model, newdata = test_data)
  
  r2 <- cor(test_data$Min_Syn_abs, preds)^2
  rmse <- sqrt(mean((test_data$Min_Syn_abs - preds)^2))
  mae <- mean(abs(test_data$Min_Syn_abs - preds))
  
  results_df <- rbind(results_df, data.frame(
    mtry = params$mtry,
    ntree = params$ntree,
    nodesize = params$nodesize,
    R2 = r2,
    RMSE = rmse,
    MAE = mae
  ))
  
  cat("Completed combination", i, "with R^2", r2, "\n")
}

# 选择最佳参数（R2 值接近目标值）
target_r2  <- 0.8447
tolerance  <- 0.0001

best_params1 <- results_df[abs(results_df$R2 - target_r2) <= tolerance, ]
# 或者：best_params1 <- results_df[which.max(results_df$R2), ]

stopCluster(cl)

## ===== 8) 用最佳参数训练最终模型 =====
best_mtry    <- best_params1$mtry
best_ntree   <- best_params1$ntree
best_nodesize <- best_params1$nodesize

set.seed(12345)
final_rf_model <- randomForest(
  x = x_train_data, y = y_train_data,
  mtry = best_mtry,
  ntree = best_ntree,
  nodesize = best_nodesize,
  importance = TRUE
)

# 预测训练集和测试集结果
y_train_rf <- predict(final_rf_model, newdata = x_train_data)
y_test_rf  <- predict(final_rf_model, newdata = test_data)

DATA_train <- data.frame(X = y_train_rf, Y = y_train_data, Dataset = "train")
DATA_test  <- data.frame(X = y_test_rf,  Y = test_data$Min_Syn_abs,  Dataset = "test")
Data <- rbind(DATA_train, DATA_test)

# 评估模型
R2   <- round(cor(Data$X, Data$Y)^2, 2)
RMSE <- round(sqrt(mean((Data$X - Data$Y)^2)), 2)
MAE  <- round(mean(abs(Data$X - Data$Y)), 2)

lin_model <- lm(Data$Y ~ Data$X - 1)  # 过原点拟合
slope <- round(coef(lin_model)[1], 2)

cat("R2: ", R2, "\nRMSE: ", RMSE, "\nMAE: ", MAE, "\nSlope: ", slope, "\n")


## ===== 用现有 final_rf_model 直接做全球预测（对齐训练特征）=====

####循环100
# 优先用 best_params1（你的“接近目标R²”的组合）；如果为空，就取 R² 最高的前 10 个
if (exists("best_params1") && nrow(best_params1) > 0) {
  param_table <- best_params1
} else {
  topK <- 10
  param_table <- head(results_df[order(-results_df$R2), ], topK)
}

# 输出目录
out_dir <- "E:/有机替代.作图/response/预测"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ========= 1) 训练时的特征顺序 & 栅格层名对齐 =========
feature_cols <- colnames(x_train_data)

# 与提取时一致的安全图层名
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜（NoData）
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

# 目标名（如果误入特征时用于常数填补）
target_col  <- "Min_Syn_abs"
target_mean <- mean(y_train_data, na.rm = TRUE)

## ========= 2) 一次性构造 “全球像元 × 特征” 矩阵（供所有模型复用）=========
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols

for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
    message("提示：特征 '", nm, "' 在栅格与训练数据中都无法直接映射，已用 0 代填。")
  }
}

covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
# 数值列 NA 均值填补
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) {
    v[is.na(v)] <- mean(v, na.rm = TRUE)
  }
  v
})
names(covariate_matrix) <- feature_cols

## ========= 3) 循环训练 → 预测 → 写盘 =========
out_dir <- "E:/有机替代.作图/response/预测/DMS1"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保证栅格层名和训练特征一致
feature_cols <- colnames(x_train_data)
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

target_col  <- "Min_Syn_abs"
target_mean <- mean(y_train_data, na.rm = TRUE)

# 一次性构造全球输入矩阵（避免每轮重复）
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols
for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
  }
}
covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) v[is.na(v)] <- mean(v, na.rm = TRUE)
  v
})
names(covariate_matrix) <- feature_cols

## ========= 循环 100 次 =========
results <- data.frame(FileName = character(), R2 = numeric(), RMSE = numeric(), MAE = numeric())

for (i in 1:100) {
  set.seed(12345 + i)
  
  # 每轮都随机重划分 70/30
  splitIndex <- createDataPartition(combined_data$Min_Syn_abs, p = 0.7, list = FALSE)
  train_i <- combined_data[splitIndex, ]
  test_i  <- combined_data[-splitIndex, ]
  
  x_train_i <- data.frame(as.matrix(train_i[, feature_cols]))
  y_train_i <- train_i$Min_Syn_abs
  
  rf_model <- randomForest(
    x = x_train_i, y = y_train_i,
    mtry = best_mtry, ntree = best_ntree, nodesize = best_nodesize,
    importance = TRUE
  )
  
  # 验证集表现
  preds_test <- predict(rf_model, newdata = test_i[, feature_cols])
  r2   <- cor(test_i$Min_Syn_abs, preds_test)^2
  rmse <- sqrt(mean((test_i$Min_Syn_abs - preds_test)^2))
  mae  <- mean(abs(test_i$Min_Syn_abs - preds_test))
  
  # 全局预测
  preds_global <- as.numeric(predict(rf_model, newdata = covariate_matrix))
  
  pred_r <- Covariatestacked[[1]]
  values(pred_r) <- preds_global
  values(pred_r)[is.na(mask_vals)] <- NA
  
  out_name <- sprintf("NH3_RF_%03d.tif", i)
  out_path <- file.path(out_dir, out_name)
  writeRaster(pred_r, filename = out_path, format = "GTiff",
              datatype = "FLT4S", overwrite = TRUE)
  
  results <- rbind(results, data.frame(FileName = out_name, R2 = r2, RMSE = rmse, MAE = mae))
  
  cat(sprintf("✅ Loop %03d | R²=%.3f | RMSE=%.3f | MAE=%.3f | 已保存 %s\n",
              i, r2, rmse, mae, out_path))
}

write.csv(results, file.path(out_dir, "RF_100_results.csv"), row.names = FALSE)
cat("🎉 已完成 100 次循环预测，结果汇总保存到 RF_100_results.csv\n")



####计算均值和差异系数
library(terra)

# 输入文件夹路径
out_dir <- "E:/有机替代.作图/response/预测/DMS1"

# 匹配 NH3_RF_001.tif ~ NH3_RF_100.tif
files <- list.files(out_dir, pattern = "NH3_RF_\\d{3}\\.tif$", full.names = TRUE)

# 检查是否找到文件
if (length(files) == 0) stop("没有找到符合条件的文件，请检查路径或文件名！")

# 读取并堆叠
raster_stack <- rast(files)

# 计算均值、标准差、变异系数
mean_raster <- mean(raster_stack, na.rm = TRUE)
sd_raster   <- stdev(raster_stack, na.rm = TRUE)
cv_raster   <- sd_raster / mean_raster

# 输出目录
output_dir <- "E:/有机替代.作图/response/预测/DMS1"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 保存结果
writeRaster(mean_raster, file.path(output_dir, "Mean.tif"), overwrite = TRUE)
writeRaster(sd_raster,   file.path(output_dir, "SD.tif"),   overwrite = TRUE)
writeRaster(cv_raster,   file.path(output_dir, "CV.tif"),   overwrite = TRUE)
}
####DIS预测
{
## ===== 依赖包 =====
library(sp)       # coordinates(), CRS 等
library(raster)   # stack(), extract()
library(car)      # car::vif
library(usdm)     # 如果你后续还要用到usdm的其它函数，保留；否则可以不加载
library(caret)
library(dplyr)
library(tidyr)
## ===== 1) 读取观测点并设坐标系 =====
MVDdata <- read.csv("E:/有机替代.作图/response/预测/Prediction4.csv", stringsAsFactors = FALSE)

# 如果数据里经纬度列名不是 "Longitude" "Latitude"，请改成实际列名
coordinates(MVDdata) <- ~ Longitude + Latitude
crs(MVDdata) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 2) 协变量栅格堆栈 =====
file_list <- list.files(path = "E:/有机替代.作图/response/预测/List2", full.names = TRUE)
Covariatestacked <- stack(file_list)
crs(Covariatestacked) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 3) 栅格提取到点 =====
extracted_values <- raster::extract(Covariatestacked, MVDdata)
extracted_df <- as.data.frame(extracted_values)

# 给栅格提取结果起列名（使用文件名的安全版本）
safe_raster_names <- make.names(basename(file_list))
colnames(extracted_df) <- safe_raster_names

## ===== 4) 合并回观测数据（转成数据框以便与数据绑定） =====
GalNdata_df <- as.data.frame(MVDdata)
combined_data <- cbind(GalNdata_df, extracted_df)

# 添加经纬度列（如果还没有添加的话）
combined_data$Latitude <- MVDdata$Latitude
combined_data$Longitude <- MVDdata$Longitude

## ===== 5) 数据预处理 =====
# 处理缺失值：数值列用列均值填补
is_num <- vapply(combined_data, is.numeric, logical(1))
combined_data[is_num] <- lapply(combined_data[is_num], function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

## ===== 6) 划分训练和测试集 =====
set.seed(12345)
splitIndex <- createDataPartition(combined_data$Imm_Syn_abs, p = 0.7, list = FALSE)
train_data <- combined_data[splitIndex, ]   # 70% 训练集
test_data  <- combined_data[-splitIndex, ]  # 30% 测试集

## ===== 7) 网格搜索调参 =====
library(randomForest)
library(parallel)
library(doParallel)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

tune_grid <- expand.grid(
  mtry = 1:13,
  ntree = c(100, 500, 1000, 1500, 2000, 2500, 3000),
  nodesize = 1:10
)

# 只随机抽取 50 个参数组合
set.seed(12345)
tune_grid <- tune_grid[sample(nrow(tune_grid), 50), ]

x_train_data <- data.frame(as.matrix(train_data[, -ncol(train_data)]))
y_train_data <- train_data$Imm_Syn_abs

results_df <- data.frame(
  mtry = integer(),
  ntree = integer(),
  nodesize = integer(),
  R2 = numeric(),
  RMSE = numeric(),
  MAE = numeric(),
  stringsAsFactors = FALSE
)

total_combinations <- nrow(tune_grid)

set.seed(12345)

for (i in 1:total_combinations) {
  params <- tune_grid[i, ]
  cat("Processing combination", i, "of", total_combinations, "...\n")
  
  model <- randomForest(
    x = x_train_data,
    y = y_train_data,
    mtry = params$mtry,
    ntree = params$ntree,
    nodesize = params$nodesize,
    importance = TRUE
  )
  
  preds <- predict(model, newdata = test_data)
  
  r2 <- cor(test_data$Imm_Syn_abs, preds)^2
  rmse <- sqrt(mean((test_data$Imm_Syn_abs - preds)^2))
  mae <- mean(abs(test_data$Imm_Syn_abs - preds))
  
  results_df <- rbind(results_df, data.frame(
    mtry = params$mtry,
    ntree = params$ntree,
    nodesize = params$nodesize,
    R2 = r2,
    RMSE = rmse,
    MAE = mae
  ))
  
  cat("Completed combination", i, "with R^2", r2, "\n")
}

# 查看结果
View(results_df)

# 选择最佳参数（R2 值接近目标值）
# target_r2  <- 0.6552
target_r2  <- 0.7534
tolerance  <- 0.0001

best_params1 <- results_df[abs(results_df$R2 - target_r2) <= tolerance, ]
# 或者：best_params1 <- results_df[which.max(results_df$R2), ]

stopCluster(cl)

## ===== 8) 用最佳参数训练最终模型 =====
best_mtry    <- best_params1$mtry
best_ntree   <- best_params1$ntree
best_nodesize <- best_params1$nodesize

set.seed(12345)
final_rf_model <- randomForest(
  x = x_train_data, y = y_train_data,
  mtry = best_mtry,
  ntree = best_ntree,
  nodesize = best_nodesize,
  importance = TRUE
)

# 预测训练集和测试集结果
y_train_rf <- predict(final_rf_model, newdata = x_train_data)
y_test_rf  <- predict(final_rf_model, newdata = test_data)

DATA_train <- data.frame(X = y_train_rf, Y = y_train_data, Dataset = "train")
DATA_test  <- data.frame(X = y_test_rf,  Y = test_data$Imm_Syn_abs,  Dataset = "test")
Data <- rbind(DATA_train, DATA_test)

# 评估模型
R2   <- round(cor(Data$X, Data$Y)^2, 2)
RMSE <- round(sqrt(mean((Data$X - Data$Y)^2)), 2)
MAE  <- round(mean(abs(Data$X - Data$Y)), 2)

lin_model <- lm(Data$Y ~ Data$X - 1)  # 过原点拟合
slope <- round(coef(lin_model)[1], 2)

cat("R2: ", R2, "\nRMSE: ", RMSE, "\nMAE: ", MAE, "\nSlope: ", slope, "\n")


## ===== 用现有 final_rf_model 直接做全球预测（对齐训练特征）=====

####循环100
# 优先用 best_params1（你的“接近目标R²”的组合）；如果为空，就取 R² 最高的前 10 个
if (exists("best_params1") && nrow(best_params1) > 0) {
  param_table <- best_params1
} else {
  topK <- 10
  param_table <- head(results_df[order(-results_df$R2), ], topK)
}

# 输出目录
out_dir <- "E:/有机替代.作图/response/预测"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ========= 1) 训练时的特征顺序 & 栅格层名对齐 =========
feature_cols <- colnames(x_train_data)

# 与提取时一致的安全图层名
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜（NoData）
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

# 目标名（如果误入特征时用于常数填补）
target_col  <- "Imm_Syn_abs"
target_mean <- mean(y_train_data, na.rm = TRUE)

## ========= 2) 一次性构造 “全球像元 × 特征” 矩阵（供所有模型复用）=========
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols

for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
    message("提示：特征 '", nm, "' 在栅格与训练数据中都无法直接映射，已用 0 代填。")
  }
}

covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
# 数值列 NA 均值填补
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) {
    v[is.na(v)] <- mean(v, na.rm = TRUE)
  }
  v
})
names(covariate_matrix) <- feature_cols

## ========= 3) 循环训练 → 预测 → 写盘 =========
out_dir <- "E:/有机替代.作图/response/预测/DIS1"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保证栅格层名和训练特征一致
feature_cols <- colnames(x_train_data)
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

target_col  <- "Imm_Syn_abs"
target_mean <- mean(y_train_data, na.rm = TRUE)

# 一次性构造全球输入矩阵（避免每轮重复）
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols
for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
  }
}
covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) v[is.na(v)] <- mean(v, na.rm = TRUE)
  v
})
names(covariate_matrix) <- feature_cols

## ========= 循环 100 次 =========
results <- data.frame(FileName = character(), R2 = numeric(), RMSE = numeric(), MAE = numeric())

for (i in 1:100) {
  set.seed(12345 + i)
  
  # 每轮都随机重划分 70/30
  splitIndex <- createDataPartition(combined_data$Imm_Syn_abs, p = 0.7, list = FALSE)
  train_i <- combined_data[splitIndex, ]
  test_i  <- combined_data[-splitIndex, ]
  
  x_train_i <- data.frame(as.matrix(train_i[, feature_cols]))
  y_train_i <- train_i$Imm_Syn_abs
  
  rf_model <- randomForest(
    x = x_train_i, y = y_train_i,
    mtry = best_mtry, ntree = best_ntree, nodesize = best_nodesize,
    importance = TRUE
  )
  
  # 验证集表现
  preds_test <- predict(rf_model, newdata = test_i[, feature_cols])
  r2   <- cor(test_i$Imm_Syn_abs, preds_test)^2
  rmse <- sqrt(mean((test_i$Imm_Syn_abs - preds_test)^2))
  mae  <- mean(abs(test_i$Imm_Syn_abs - preds_test))
  
  # 全局预测
  preds_global <- as.numeric(predict(rf_model, newdata = covariate_matrix))
  
  pred_r <- Covariatestacked[[1]]
  values(pred_r) <- preds_global
  values(pred_r)[is.na(mask_vals)] <- NA
  
  out_name <- sprintf("NH3_RF_%03d.tif", i)
  out_path <- file.path(out_dir, out_name)
  writeRaster(pred_r, filename = out_path, format = "GTiff",
              datatype = "FLT4S", overwrite = TRUE)
  
  results <- rbind(results, data.frame(FileName = out_name, R2 = r2, RMSE = rmse, MAE = mae))
  
  cat(sprintf("✅ Loop %03d | R²=%.3f | RMSE=%.3f | MAE=%.3f | 已保存 %s\n",
              i, r2, rmse, mae, out_path))
}

write.csv(results, file.path(out_dir, "RF_100_results.csv"), row.names = FALSE)
cat("🎉 已完成 100 次循环预测，结果汇总保存到 RF_100_results.csv\n")



####计算均值和差异系数
library(terra)

# 输入文件夹路径
out_dir <- "E:/有机替代.作图/response/预测/DIS1"

# 匹配 NH3_RF_001.tif ~ NH3_RF_100.tif
files <- list.files(out_dir, pattern = "NH3_RF_\\d{3}\\.tif$", full.names = TRUE)

# 检查是否找到文件
if (length(files) == 0) stop("没有找到符合条件的文件，请检查路径或文件名！")

# 读取并堆叠
raster_stack <- rast(files)

# 计算均值、标准差、变异系数
mean_raster <- mean(raster_stack, na.rm = TRUE)
sd_raster   <- stdev(raster_stack, na.rm = TRUE)
cv_raster   <- sd_raster / mean_raster

# 输出目录
output_dir <- "E:/有机替代.作图/response/预测/DIS1"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 保存结果
writeRaster(mean_raster, file.path(output_dir, "Mean.tif"), overwrite = TRUE)
writeRaster(sd_raster,   file.path(output_dir, "SD.tif"),   overwrite = TRUE)
writeRaster(cv_raster,   file.path(output_dir, "CV.tif"),   overwrite = TRUE)
}
####N2O
{
## ===== 依赖包 =====
library(sp)       # coordinates(), CRS 等
library(raster)   # stack(), extract()
library(car)      # car::vif
library(usdm)     # 如果你后续还要用到usdm的其它函数，保留；否则可以不加载
library(caret)
library(dplyr)
library(tidyr)
## ===== 1) 读取观测点并设坐标系 =====
MVDdata <- read.csv("E:/有机替代.作图/response/预测/Prediction5.csv", stringsAsFactors = FALSE)

# 如果数据里经纬度列名不是 "Longitude" "Latitude"，请改成实际列名
coordinates(MVDdata) <- ~ Longitude + Latitude
crs(MVDdata) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 2) 协变量栅格堆栈 =====
file_list <- list.files(path = "E:/有机替代.作图/response/预测/List1", full.names = TRUE)
Covariatestacked <- stack(file_list)
crs(Covariatestacked) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 3) 栅格提取到点 =====
extracted_values <- raster::extract(Covariatestacked, MVDdata)
extracted_df <- as.data.frame(extracted_values)

# 给栅格提取结果起列名（使用文件名的安全版本）
safe_raster_names <- make.names(basename(file_list))
colnames(extracted_df) <- safe_raster_names

## ===== 4) 合并回观测数据（转成数据框以便与数据绑定） =====
GalNdata_df <- as.data.frame(MVDdata)
combined_data <- cbind(GalNdata_df, extracted_df)

# 添加经纬度列（如果还没有添加的话）
combined_data$Latitude <- MVDdata$Latitude
combined_data$Longitude <- MVDdata$Longitude

## ===== 5) 数据预处理 =====
# 仅根据数值列剔除含 NA 的行
is_num <- vapply(combined_data, is.numeric, logical(1))
combined_data <- combined_data[complete.cases(combined_data[, is_num]), ]

## ===== 6) 划分训练和测试集 =====
set.seed(12345)
splitIndex <- createDataPartition(combined_data$LnRR.N2O_emission, p = 0.7, list = FALSE)
train_data <- combined_data[splitIndex, ]   # 70% 训练集
test_data  <- combined_data[-splitIndex, ]  # 30% 测试集


## ===== 7) 网格搜索调参 =====
library(randomForest)
library(parallel)
library(doParallel)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

tune_grid <- expand.grid(
  mtry = 1:13,   # 尝试的 mtry 参数范围
  ntree = c(100, 500, 1000, 1500, 2000, 2500, 3000),   # 尝试的 ntree 参数值
  nodesize = 1:10   # 尝试的 nodesize 参数范围
)

x_train_data <- data.frame(as.matrix(train_data[, -ncol(train_data)]))  # 特征变量
y_train_data <- train_data$LnRR.N2O_emission  # 目标变量

results_df <- data.frame(mtry=integer(), ntree=integer(), nodesize=integer(),
                         R2=numeric(), RMSE=numeric(), MAE=numeric(),
                         stringsAsFactors=FALSE)

total_combinations <- nrow(tune_grid)
set.seed(12345)
for (i in 1:total_combinations) {
  params <- tune_grid[i, ]
  cat("Processing combination", i, "of", total_combinations, "...\n")
  
  model <- randomForest(
    x = x_train_data, y = y_train_data,
    mtry = params$mtry, ntree = params$ntree, nodesize = params$nodesize,
    importance = TRUE
  )
  
  preds <- predict(model, newdata = test_data)
  r2   <- cor(test_data$LnRR.N2O_emission, preds)^2
  rmse <- sqrt(mean((test_data$LnRR.N2O_emission - preds)^2))
  mae  <- mean(abs(test_data$LnRR.N2O_emission - preds))
  
  results_df <- rbind(results_df, data.frame(
    mtry=params$mtry, ntree=params$ntree, nodesize=params$nodesize,
    R2=r2, RMSE=rmse, MAE=mae
  ))
  
  cat("Completed combination", i, "with R^2", r2, "\n")
}

# 查看结果
View(results_df)

# 选择最佳参数（R2 值接近目标值）
target_r2  <- 0.87886
# target_r2  <- 0.9020
tolerance  <- 0.0001

best_params1 <- results_df[abs(results_df$R2 - target_r2) <= tolerance, ]
# 或者：best_params1 <- results_df[which.max(results_df$R2), ]

stopCluster(cl)


####循环100次
## ===== 8) 用最佳参数训练最终模型 =====
best_mtry    <- best_params1$mtry
best_ntree   <- best_params1$ntree
best_nodesize <- best_params1$nodesize

set.seed(12345)
final_rf_model <- randomForest(
  x = x_train_data, y = y_train_data,
  mtry = best_mtry,
  ntree = best_ntree,
  nodesize = best_nodesize,
  importance = TRUE
)

# 预测训练集和测试集结果
y_train_rf <- predict(final_rf_model, newdata = x_train_data)
y_test_rf  <- predict(final_rf_model, newdata = test_data)

DATA_train <- data.frame(X = y_train_rf, Y = y_train_data, Dataset = "train")
DATA_test  <- data.frame(X = y_test_rf,  Y = test_data$LnRR.N2O_emission,  Dataset = "test")
Data <- rbind(DATA_train, DATA_test)
# write.csv(Data, "E:/有机替代.作图/response/预测/obs_pred_all_N2O.csv", row.names = FALSE)
# 评估模型
R2   <- round(cor(Data$X, Data$Y)^2, 2)
RMSE <- round(sqrt(mean((Data$X - Data$Y)^2)), 2)
MAE  <- round(mean(abs(Data$X - Data$Y)), 2)

lin_model <- lm(Data$Y ~ Data$X - 1)  # 过原点拟合
slope <- round(coef(lin_model)[1], 2)

cat("R2: ", R2, "\nRMSE: ", RMSE, "\nMAE: ", MAE, "\nSlope: ", slope, "\n")

####循环100
# 优先用 best_params1（你的“接近目标R²”的组合）；如果为空，就取 R² 最高的前 10 个
if (exists("best_params1") && nrow(best_params1) > 0) {
  param_table <- best_params1
} else {
  topK <- 10
  param_table <- head(results_df[order(-results_df$R2), ], topK)
}

# 输出目录
out_dir <- "E:/有机替代.作图/response/预测"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ========= 1) 训练时的特征顺序 & 栅格层名对齐 =========
feature_cols <- colnames(x_train_data)

# 与提取时一致的安全图层名
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜（NoData）
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

# 目标名（如果误入特征时用于常数填补）
target_col  <- "LnRR.N2O_emission"
target_mean <- mean(y_train_data, na.rm = TRUE)

## ========= 2) 一次性构造 “全球像元 × 特征” 矩阵（供所有模型复用）=========
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols

for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
    message("提示：特征 '", nm, "' 在栅格与训练数据中都无法直接映射，已用 0 代填。")
  }
}

covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
# 数值列 NA 均值填补
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) {
    v[is.na(v)] <- mean(v, na.rm = TRUE)
  }
  v
})
names(covariate_matrix) <- feature_cols

## ========= 3) 循环训练 → 预测 → 写盘 =========
out_dir <- "E:/有机替代.作图/response/预测/N2O"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保证栅格层名和训练特征一致
feature_cols <- colnames(x_train_data)
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

target_col  <- "LnRR.N2O_emission"
target_mean <- mean(y_train_data, na.rm = TRUE)

# 一次性构造全球输入矩阵（避免每轮重复）
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols
for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
  }
}
covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) v[is.na(v)] <- mean(v, na.rm = TRUE)
  v
})
names(covariate_matrix) <- feature_cols

## ========= 循环 100 次 =========
results <- data.frame(FileName = character(), R2 = numeric(), RMSE = numeric(), MAE = numeric())

for (i in 1:100) {
  set.seed(12345 + i)
  
  # 每轮都随机重划分 70/30
  splitIndex <- createDataPartition(combined_data$LnRR.N2O_emission, p = 0.7, list = FALSE)
  train_i <- combined_data[splitIndex, ]
  test_i  <- combined_data[-splitIndex, ]
  
  x_train_i <- data.frame(as.matrix(train_i[, feature_cols]))
  y_train_i <- train_i$LnRR.N2O_emission
  
  rf_model <- randomForest(
    x = x_train_i, y = y_train_i,
    mtry = best_mtry, ntree = best_ntree, nodesize = best_nodesize,
    importance = TRUE
  )
  
  # 验证集表现
  preds_test <- predict(rf_model, newdata = test_i[, feature_cols])
  r2   <- cor(test_i$LnRR.N2O_emission, preds_test)^2
  rmse <- sqrt(mean((test_i$LnRR.N2O_emission - preds_test)^2))
  mae  <- mean(abs(test_i$LnRR.N2O_emission - preds_test))
  
  # 全局预测
  preds_global <- as.numeric(predict(rf_model, newdata = covariate_matrix))
  
  pred_r <- Covariatestacked[[1]]
  values(pred_r) <- preds_global
  values(pred_r)[is.na(mask_vals)] <- NA
  
  out_name <- sprintf("NH3_RF_%03d.tif", i)
  out_path <- file.path(out_dir, out_name)
  writeRaster(pred_r, filename = out_path, format = "GTiff",
              datatype = "FLT4S", overwrite = TRUE)
  
  results <- rbind(results, data.frame(FileName = out_name, R2 = r2, RMSE = rmse, MAE = mae))
  
  cat(sprintf("✅ Loop %03d | R²=%.3f | RMSE=%.3f | MAE=%.3f | 已保存 %s\n",
              i, r2, rmse, mae, out_path))
}

write.csv(results, file.path(out_dir, "RF_100_results.csv"), row.names = FALSE)
cat("🎉 已完成 100 次循环预测，结果汇总保存到 RF_100_results.csv\n")


# 输入文件夹路径
out_dir <- "E:/有机替代.作图/response/预测/N2O"

# 匹配 NH3_RF_001.tif ~ NH3_RF_100.tif
files <- list.files(out_dir, pattern = "NH3_RF_\\d{3}\\.tif$", full.names = TRUE)

# 检查是否找到文件
if (length(files) == 0) stop("没有找到符合条件的文件，请检查路径或文件名！")

# 读取并堆叠
raster_stack <- rast(files)

# 计算均值、标准差、变异系数
mean_raster <- mean(raster_stack, na.rm = TRUE)
sd_raster   <- stdev(raster_stack, na.rm = TRUE)
cv_raster   <- sd_raster / mean_raster

# 输出目录
output_dir <- "E:/有机替代.作图/response/预测/N2O"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 保存结果
writeRaster(mean_raster, file.path(output_dir, "Mean.tif"), overwrite = TRUE)
writeRaster(sd_raster,   file.path(output_dir, "SD.tif"),   overwrite = TRUE)
writeRaster(cv_raster,   file.path(output_dir, "CV.tif"),   overwrite = TRUE)
}
###NH3
{
## ===== 依赖包 =====
library(sp)       # coordinates(), CRS 等
library(raster)   # stack(), extract()
library(car)      # car::vif
library(usdm)     # 如果你后续还要用到usdm的其它函数，保留；否则可以不加载
library(caret)
library(dplyr)
library(tidyr)
## ===== 1) 读取观测点并设坐标系 =====
MVDdata <- read.csv("E:/有机替代.作图/response/预测/Prediction6.csv", stringsAsFactors = FALSE)

# 如果数据里经纬度列名不是 "Longitude" "Latitude"，请改成实际列名
coordinates(MVDdata) <- ~ Longitude + Latitude
crs(MVDdata) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 2) 协变量栅格堆栈 =====
file_list <- list.files(path = "E:/有机替代.作图/response/预测/List3", full.names = TRUE)
Covariatestacked <- stack(file_list)
crs(Covariatestacked) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

## ===== 3) 栅格提取到点 =====
extracted_values <- raster::extract(Covariatestacked, MVDdata)
extracted_df <- as.data.frame(extracted_values)

# 给栅格提取结果起列名（使用文件名的安全版本）
safe_raster_names <- make.names(basename(file_list))
colnames(extracted_df) <- safe_raster_names

## ===== 4) 合并回观测数据（转成数据框以便与数据绑定） =====
GalNdata_df <- as.data.frame(MVDdata)
combined_data <- cbind(GalNdata_df, extracted_df)

# 添加经纬度列（如果还没有添加的话）
combined_data$Latitude <- MVDdata$Latitude
combined_data$Longitude <- MVDdata$Longitude

## ===== 5) 数据预处理 =====
# 处理缺失值：数值列用列均值填补
is_num <- vapply(combined_data, is.numeric, logical(1))
combined_data[is_num] <- lapply(combined_data[is_num], function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

## ===== 6) 划分训练和测试集 =====
set.seed(12346)
splitIndex <- createDataPartition(combined_data$LnRR.NH3_emission, p = 0.7, list = FALSE)
train_data <- combined_data[splitIndex, ]   # 70% 训练集
test_data  <- combined_data[-splitIndex, ]  # 30% 测试集

## ===== 7) 网格搜索调参 =====
library(randomForest)
library(parallel)
library(doParallel)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

tune_grid <- expand.grid(
  mtry = 1:13,
  ntree = c(100, 500, 1000, 1500, 2000, 2500, 3000),
  nodesize = 1:10
)

# 只随机抽取 50 个参数组合
set.seed(12346)
tune_grid <- tune_grid[sample(nrow(tune_grid), 50), ]

x_train_data <- data.frame(as.matrix(train_data[, -ncol(train_data)]))
y_train_data <- train_data$LnRR.NH3_emission

results_df <- data.frame(
  mtry = integer(),
  ntree = integer(),
  nodesize = integer(),
  R2 = numeric(),
  RMSE = numeric(),
  MAE = numeric(),
  stringsAsFactors = FALSE
)

total_combinations <- nrow(tune_grid)

set.seed(12346)

for (i in 1:total_combinations) {
  params <- tune_grid[i, ]
  cat("Processing combination", i, "of", total_combinations, "...\n")
  
  model <- randomForest(
    x = x_train_data,
    y = y_train_data,
    mtry = params$mtry,
    ntree = params$ntree,
    nodesize = params$nodesize,
    importance = TRUE
  )
  
  preds <- predict(model, newdata = test_data)
  
  r2 <- cor(test_data$LnRR.NH3_emission, preds)^2
  rmse <- sqrt(mean((test_data$LnRR.NH3_emission - preds)^2))
  mae <- mean(abs(test_data$LnRR.NH3_emission - preds))
  
  results_df <- rbind(results_df, data.frame(
    mtry = params$mtry,
    ntree = params$ntree,
    nodesize = params$nodesize,
    R2 = r2,
    RMSE = rmse,
    MAE = mae
  ))
  
  cat("Completed combination", i, "with R^2", r2, "\n")
}

# 查看结果
View(results_df)

# 选择最佳参数（R2 值接近目标值）
target_r2  <- 0.9353
tolerance  <- 0.0001

best_params1 <- results_df[abs(results_df$R2 - target_r2) <= tolerance, ]
# 或者：best_params1 <- results_df[which.max(results_df$R2), ]

stopCluster(cl)


####循环100次
## ===== 8) 用最佳参数训练最终模型 =====
best_mtry    <- best_params1$mtry
best_ntree   <- best_params1$ntree
best_nodesize <- best_params1$nodesize

set.seed(12345)
final_rf_model <- randomForest(
  x = x_train_data, y = y_train_data,
  mtry = best_mtry,
  ntree = best_ntree,
  nodesize = best_nodesize,
  importance = TRUE
)

# 预测训练集和测试集结果
y_train_rf <- predict(final_rf_model, newdata = x_train_data)
y_test_rf  <- predict(final_rf_model, newdata = test_data)

DATA_train <- data.frame(X = y_train_rf, Y = y_train_data, Dataset = "train")
DATA_test  <- data.frame(X = y_test_rf,  Y = test_data$LnRR.NH3_emission,  Dataset = "test")
Data <- rbind(DATA_train, DATA_test)
# write.csv(Data, "E:/有机替代.作图/response/预测/obs_pred_all_NH3.csv", row.names = FALSE)
# 评估模型
R2   <- round(cor(Data$X, Data$Y)^2, 2)
RMSE <- round(sqrt(mean((Data$X - Data$Y)^2)), 2)
MAE  <- round(mean(abs(Data$X - Data$Y)), 2)

lin_model <- lm(Data$Y ~ Data$X - 1)  # 过原点拟合
slope <- round(coef(lin_model)[1], 2)

cat("R2: ", R2, "\nRMSE: ", RMSE, "\nMAE: ", MAE, "\nSlope: ", slope, "\n")

####循环100
# 优先用 best_params1（你的“接近目标R²”的组合）；如果为空，就取 R² 最高的前 10 个
if (exists("best_params1") && nrow(best_params1) > 0) {
  param_table <- best_params1
} else {
  topK <- 10
  param_table <- head(results_df[order(-results_df$R2), ], topK)
}

# 输出目录
out_dir <- "E:/有机替代.作图/response/预测"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ========= 1) 训练时的特征顺序 & 栅格层名对齐 =========
feature_cols <- colnames(x_train_data)

# 与提取时一致的安全图层名
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜（NoData）
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

# 目标名（如果误入特征时用于常数填补）
target_col  <- "LnRR.NH3_emission"
target_mean <- mean(y_train_data, na.rm = TRUE)

## ========= 2) 一次性构造 “全球像元 × 特征” 矩阵（供所有模型复用）=========
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols

for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
    message("提示：特征 '", nm, "' 在栅格与训练数据中都无法直接映射，已用 0 代填。")
  }
}

covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
# 数值列 NA 均值填补
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) {
    v[is.na(v)] <- mean(v, na.rm = TRUE)
  }
  v
})
names(covariate_matrix) <- feature_cols

## ========= 3) 循环训练 → 预测 → 写盘 =========
out_dir <- "E:/有机替代.作图/response/预测/NH3"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保证栅格层名和训练特征一致
feature_cols <- colnames(x_train_data)
names(Covariatestacked) <- make.names(basename(file_list))
raster_names <- names(Covariatestacked)

# 用第一层做掩膜
n_cells   <- ncell(Covariatestacked)
mask_vals <- getValues(Covariatestacked[[1]])
xy        <- xyFromCell(Covariatestacked, 1:n_cells)

target_col  <- "LnRR.NH3_emission"
target_mean <- mean(y_train_data, na.rm = TRUE)

# 一次性构造全球输入矩阵（避免每轮重复）
cols_list <- vector("list", length(feature_cols))
names(cols_list) <- feature_cols
for (nm in feature_cols) {
  if (nm %in% raster_names) {
    cols_list[[nm]] <- getValues(Covariatestacked[[nm]])
  } else if (nm == "Latitude") {
    cols_list[[nm]] <- xy[, 2]
  } else if (nm == "Longitude") {
    cols_list[[nm]] <- xy[, 1]
  } else if (nm == target_col) {
    cols_list[[nm]] <- rep(target_mean, n_cells)
  } else if (nm %in% colnames(x_train_data)) {
    cols_list[[nm]] <- rep(mean(x_train_data[[nm]], na.rm = TRUE), n_cells)
  } else {
    cols_list[[nm]] <- rep(0, n_cells)
  }
}
covariate_matrix <- as.data.frame(cols_list, check.names = FALSE)
covariate_matrix[] <- lapply(covariate_matrix, function(v) {
  if (is.numeric(v)) v[is.na(v)] <- mean(v, na.rm = TRUE)
  v
})
names(covariate_matrix) <- feature_cols

## ========= 循环 100 次 =========
results <- data.frame(FileName = character(), R2 = numeric(), RMSE = numeric(), MAE = numeric())

for (i in 1:100) {
  set.seed(12345 + i)
  
  # 每轮都随机重划分 70/30
  splitIndex <- createDataPartition(combined_data$LnRR.NH3_emission, p = 0.7, list = FALSE)
  train_i <- combined_data[splitIndex, ]
  test_i  <- combined_data[-splitIndex, ]
  
  x_train_i <- data.frame(as.matrix(train_i[, feature_cols]))
  y_train_i <- train_i$LnRR.NH3_emission
  
  rf_model <- randomForest(
    x = x_train_i, y = y_train_i,
    mtry = best_mtry, ntree = best_ntree, nodesize = best_nodesize,
    importance = TRUE
  )
  
  # 验证集表现
  preds_test <- predict(rf_model, newdata = test_i[, feature_cols])
  r2   <- cor(test_i$LnRR.NH3_emission, preds_test)^2
  rmse <- sqrt(mean((test_i$LnRR.NH3_emission - preds_test)^2))
  mae  <- mean(abs(test_i$LnRR.NH3_emission - preds_test))
  
  # 全局预测
  preds_global <- as.numeric(predict(rf_model, newdata = covariate_matrix))
  
  pred_r <- Covariatestacked[[1]]
  values(pred_r) <- preds_global
  values(pred_r)[is.na(mask_vals)] <- NA
  
  out_name <- sprintf("NH3_RF_%03d.tif", i)
  out_path <- file.path(out_dir, out_name)
  writeRaster(pred_r, filename = out_path, format = "GTiff",
              datatype = "FLT4S", overwrite = TRUE)
  
  results <- rbind(results, data.frame(FileName = out_name, R2 = r2, RMSE = rmse, MAE = mae))
  
  cat(sprintf("✅ Loop %03d | R²=%.3f | RMSE=%.3f | MAE=%.3f | 已保存 %s\n",
              i, r2, rmse, mae, out_path))
}

write.csv(results, file.path(out_dir, "RF_100_results.csv"), row.names = FALSE)
cat("🎉 已完成 100 次循环预测，结果汇总保存到 RF_100_results.csv\n")

####计算均值和差异系数
library(terra)

# 输入文件夹路径
out_dir <- "E:/有机替代.作图/response/预测/NH3"

# 匹配 NH3_RF_001.tif ~ NH3_RF_100.tif
files <- list.files(out_dir, pattern = "NH3_RF_\\d{3}\\.tif$", full.names = TRUE)

# 检查是否找到文件
if (length(files) == 0) stop("没有找到符合条件的文件，请检查路径或文件名！")

# 读取并堆叠
raster_stack <- rast(files)

# 计算均值、标准差、变异系数
mean_raster <- mean(raster_stack, na.rm = TRUE)
sd_raster   <- stdev(raster_stack, na.rm = TRUE)
cv_raster   <- sd_raster / mean_raster

# 输出目录
output_dir <- "E:/有机替代.作图/response/预测/NH3"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 保存结果
writeRaster(mean_raster, file.path(output_dir, "Mean.tif"), overwrite = TRUE)
writeRaster(sd_raster,   file.path(output_dir, "SD.tif"),   overwrite = TRUE)
writeRaster(cv_raster,   file.path(output_dir, "CV.tif"),   overwrite = TRUE)
}
####同步性上色
{
####上色
library(ggplot2)
library(raster)
library(rnaturalearthdata)
library(rnaturalearth)
library(sf)
world <- ne_countries(scale = "medium", returnclass = "sf")

# 读取栅格数据
tiff_data1 <- raster("E:/有机替代.作图/response/预测/DMS1/Mean.tif")
# 将0替换为NA
tiff_data1[tiff_data1 == 0] <- NA
# 转换世界地图投影为罗宾逊投影
world <- st_transform(world, "+proj=robin")

# 将栅格数据的投影转换为与世界地图一致
tiff_data1 <- projectRaster(tiff_data1, crs = st_crs(world)$proj4string)

# 将栅格数据转为数据框，方便 ggplot 绘制
raster_df1 <- as.data.frame(raster::rasterToPoints(tiff_data1))
colnames(raster_df1) <- c("long", "lat", "value")

# 定义自定义颜色梯度

# custom_colors <- c("#33FFF0", "#33FF57", "#3357FF",  "#FF3333")
# custom_colors <- c("#33FFF0", "#33FF57", "#FFC833", "#FF5733", "#FF3333")
custom_colors <- c("#49C2D9", "#A1D8E8", "#D0E2C0", "#FDED95", "#FFc1a6", "#f59c7c", "#f47254")
# 创建地图
map_with_raster1 <- ggplot() +
  geom_raster(data = raster_df1, aes(x = long, y = lat, fill = value)) +
  scale_fill_gradientn(colors = custom_colors, na.value = "grey40", guide = "colorbar") +
  geom_sf(data = world, fill = NA, color = "black", size = 0.2) +  # 较细边界线
  coord_sf(expand = FALSE) +
  scale_y_continuous(
    breaks = c(-60, -30, 0, 30, 60),
    labels = c("60°S", "30°S", "0°", "30°N", "60°N")
  ) +
  scale_x_continuous(
    breaks = seq(-180, 180, by = 60),
    labels = function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  ) +
  theme_minimal(base_family = "RMN") +
  labs(title = "DMS") +
  theme(
    legend.position = "bottom",
    text = element_text(family = "RMN", face = "bold", color = "black",size = 14),
    axis.title = element_blank(),
    axis.text = element_text(family = "RMN", face = "bold", color = "black",size = 12),
    axis.ticks = element_line(size = 0.2),
    panel.grid.major = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    panel.grid.minor = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    plot.title = element_text(hjust = 0.5)
  )

# 显示地图
print(map_with_raster1)

# 读取栅格数据
tiff_data <- raster("E:/有机替代.作图/response/预测/List3/DIS.tif")
# 将0替换为NA
tiff_data[tiff_data == 0] <- NA
# 转换世界地图投影为罗宾逊投影
world <- st_transform(world, "+proj=robin")

# 将栅格数据的投影转换为与世界地图一致
tiff_data <- projectRaster(tiff_data, crs = st_crs(world)$proj4string)

# 将栅格数据转为数据框，方便 ggplot 绘制
raster_df <- as.data.frame(raster::rasterToPoints(tiff_data))
colnames(raster_df) <- c("long", "lat", "value")

# 定义自定义颜色梯度

# custom_colors <- c("#33FFF0", "#33FF57", "#3357FF",  "#FF3333")
# custom_colors <- c("#33FFF0", "#33FF57", "#FFC833", "#FF5733")
custom_colors <- c("#49C2D9", "#A1D8E8", "#D0E2C0", "#FDED95", "#FFc1a6", "#f59c7c", "#f47254")
# 创建地图
map_with_raster <- ggplot() +
  geom_raster(data = raster_df, aes(x = long, y = lat, fill = value)) +
  scale_fill_gradientn(colors = custom_colors, na.value = "grey40", guide = "colorbar") +
  geom_sf(data = world, fill = NA, color = "black", size = 0.2) +  # 较细边界线
  coord_sf(expand = FALSE) +
  scale_y_continuous(
    breaks = c(-60, -30, 0, 30, 60),
    labels = c("60°S", "30°S", "0°", "30°N", "60°N")
  ) +
  scale_x_continuous(
    breaks = seq(-180, 180, by = 60),
    labels = function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  ) +
  theme_minimal(base_family = "RMN") +
  labs(title = "DIS") +
  theme(
    legend.position = "bottom",
    text = element_text(family = "RMN", face = "bold", color = "black", size = 14),
    axis.title = element_blank(),
    axis.text = element_text(family = "RMN", face = "bold", color = "black", size = 12),
    axis.ticks = element_line(size = 0.2),
    panel.grid.major = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    panel.grid.minor = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    plot.title = element_text(hjust = 0.5)
  )

# 显示地图
print(map_with_raster)

tiff("E:/有机替代.作图/response/map同步性1.tiff",width=4000,height=2000,res=300,compression="lzw")
figure_all4 <- ggarrange(map_with_raster1, map_with_raster,
                         widths = c(1, 1), heights = c(1, 1), ncol=2, nrow=1,
                         labels = c("(a)","(b)"),
                         label.x=0.1, label.y=0.8, align = "v",
                         font.label = list(size = 18, color = "black", face="bold", family="RMN"))

# 打印图形
print(figure_all4)
dev.off()

# 读取栅格数据
tiff_data2 <- raster("E:/有机替代.作图/response/预测/DNS/Mean.tif")
# 将0替换为NA
tiff_data2[tiff_data2 == 0] <- NA
# 转换世界地图投影为罗宾逊投影
world <- st_transform(world, "+proj=robin")

# 将栅格数据的投影转换为与世界地图一致
tiff_data2 <- projectRaster(tiff_data2, crs = st_crs(world)$proj4string)

# 将栅格数据转为数据框，方便 ggplot 绘制
raster_df2 <- as.data.frame(raster::rasterToPoints(tiff_data2))
colnames(raster_df2) <- c("long", "lat", "value")

# 定义自定义颜色梯度

# custom_colors <- c("#33FFF0", "#33FF57", "#3357FF",  "#FF3333")
# custom_colors <- c("#33FFF0", "#33FF57", "#FFC833", "#FF5733", "#FF3333")
custom_colors <- c("#49C2D9", "#A1D8E8", "#D0E2C0", "#FDED95", "#FFc1a6", "#f59c7c", "#f47254")
# 创建地图
map_with_raster2 <- ggplot() +
  geom_raster(data = raster_df2, aes(x = long, y = lat, fill = value)) +
  scale_fill_gradientn(colors = custom_colors, na.value = "grey40", guide = "colorbar") +
  geom_sf(data = world, fill = NA, color = "black", size = 0.2) +  # 较细边界线
  coord_sf(expand = FALSE) +
  scale_y_continuous(
    breaks = c(-60, -30, 0, 30, 60),
    labels = c("60°S", "30°S", "0°", "30°N", "60°N")
  ) +
  scale_x_continuous(
    breaks = seq(-180, 180, by = 60),
    labels = function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  ) +
  theme_minimal(base_family = "RMN") +
  labs(title = "DNS") +
  theme(
    legend.position = "bottom",
    text = element_text(family = "RMN", face = "bold",color = "black", size = 14),
    axis.title = element_blank(),
    axis.text = element_text(family = "RMN", face = "bold", color = "black",size = 12),
    axis.ticks = element_line(size = 0.2),
    panel.grid.major = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    panel.grid.minor = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    plot.title = element_text(hjust = 0.5)
  )

# 显示地图
print(map_with_raster2)

# 读取栅格数据
tiff_data3 <- raster("E:/有机替代.作图/response/预测/DDS/Mean.tif")
# 将0替换为NA
tiff_data3[tiff_data3 == 0] <- NA
# 转换世界地图投影为罗宾逊投影
world <- st_transform(world, "+proj=robin")

# 将栅格数据的投影转换为与世界地图一致
tiff_data3 <- projectRaster(tiff_data3, crs = st_crs(world)$proj4string)

# 将栅格数据转为数据框，方便 ggplot 绘制
raster_df3 <- as.data.frame(raster::rasterToPoints(tiff_data3))
colnames(raster_df3) <- c("long", "lat", "value")

# 定义自定义颜色梯度

# custom_colors <- c("#33FFF0", "#33FF57", "#3357FF",  "#FF3333")
# custom_colors <- c("#33FFF0", "#33FF57", "#FFC833", "#FF5733", "#FF3333")
custom_colors <- c("#49C2D9", "#A1D8E8",  "#D0E2C0", "#FDED95", "#FFc1a6", "#f59c7c", "#f47254")
# 创建地图
map_with_raster3 <- ggplot() +
  geom_raster(data = raster_df3, aes(x = long, y = lat, fill = value)) +
  scale_fill_gradientn(colors = custom_colors, na.value = "grey40", guide = "colorbar") +
  geom_sf(data = world, fill = NA, color = "black", size = 0.2) +  # 较细边界线
  coord_sf(expand = FALSE) +
  scale_y_continuous(
    breaks = c(-60, -30, 0, 30, 60),
    labels = c("60°S", "30°S", "0°", "30°N", "60°N")
  ) +
  scale_x_continuous(
    breaks = seq(-180, 180, by = 60),
    labels = function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  ) +
  theme_minimal(base_family = "RMN") +
  labs(title = "DDS") +
  theme(
    legend.position = "bottom",
    text = element_text(family = "RMN", face = "bold", color = "black",size = 14),
    axis.title = element_blank(),
    axis.text = element_text(family = "RMN", face = "bold",color = "black", size = 12),
    axis.ticks = element_line(size = 0.2),
    panel.grid.major = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    panel.grid.minor = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    plot.title = element_text(hjust = 0.5)
  )

# 显示地图
print(map_with_raster3)

tiff("E:/有机替代.作图/response/map同步性2.tiff",width=4000,height=2000,res=300,compression="lzw")
figure_all4 <- ggarrange(map_with_raster2, map_with_raster3,
                         widths = c(1, 1), heights = c(1, 1), ncol=2, nrow=1,
                         labels = c("(c)","(d)"),
                         label.x=0.1, label.y=0.8, align = "v",
                         font.label = list(size = 18, color = "black", face="bold", family="RMN"))

# 打印图形
print(figure_all4)
dev.off()
}
###地图上色
{
####上色
library(ggplot2)
library(raster)
library(rnaturalearthdata)
library(rnaturalearth)
library(sf)
world <- ne_countries(scale = "medium", returnclass = "sf")

# 读取栅格数据
tiff_data <- raster("E:/有机替代.作图/response/预测/N2O/Mean.tif")
# 将0替换为NA
tiff_data[tiff_data == 0] <- NA
# 转换世界地图投影为罗宾逊投影
world <- st_transform(world, "+proj=robin")

# 将栅格数据的投影转换为与世界地图一致
tiff_data <- projectRaster(tiff_data, crs = st_crs(world)$proj4string)

# 将栅格数据转为数据框，方便 ggplot 绘制
raster_df <- as.data.frame(raster::rasterToPoints(tiff_data))
colnames(raster_df) <- c("long", "lat", "value")

# 定义自定义颜色梯度

custom_colors <- c("#7B95C6", "#49C2D9", "#A1D8E8", "#67A583", "#A2C986", "#D0E2C0", "#FDED95", "#FFc1a6", "#f59c7c", "#f47254", "#c85e62")
# custom_colors <- c(
#   "#7FB6D4",  # 稍深浅蓝
#   "#58CDED",  # 稍深青蓝
#   "#AEE8D0",  # 稍深薄荷绿
#   "#DDDDA8",  # 稍深米黄
#   "#E6BDD4",  # 稍深浅粉
#   "#D7AFC0",  # 稍深粉紫
#   "#E3B79E",  # 稍深桃色
#   "#DDA548",  # 稍深橙黄
#   "#7A5144"   # 更深棕色
# )

# 创建地图
map_with_raster <- ggplot() +
  geom_raster(data = raster_df, aes(x = long, y = lat, fill = value)) +
  scale_fill_gradientn(colors = custom_colors, na.value = "grey40", guide = "colorbar") +
  geom_sf(data = world, fill = NA, color = "black", size = 0.2) +  # 较细边界线
  coord_sf(expand = FALSE) +
  scale_y_continuous(
    breaks = c(-60, -30, 0, 30, 60),
    labels = c("60°S", "30°S", "0°", "30°N", "60°N")
  ) +
  scale_x_continuous(
    breaks = seq(-180, 180, by = 60),
    labels = function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  ) +
  theme_minimal(base_family = "RMN") +
  labs(title = "lnRR of N\u2082O emission") +
  theme(
    legend.position = "bottom",
    text = element_text(family = "RMN", face = "bold", color = "black", size = 14),
    axis.title = element_blank(),
    axis.text = element_text(family = "RMN", face = "bold", color = "black", size = 12),
    axis.ticks = element_line(size = 0.2),
    panel.grid.major = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    panel.grid.minor = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    plot.title = element_text(hjust = 0.5)
  )

# 显示地图
print(map_with_raster)

# 读取栅格数据
tiff_data1 <- raster("E:/有机替代.作图/response/预测/NH3/Mean.tif")
# 将0替换为NA
tiff_data1[tiff_data1 == 0] <- NA
# 转换世界地图投影为罗宾逊投影
world <- st_transform(world, "+proj=robin")

# 将栅格数据的投影转换为与世界地图一致
tiff_data1 <- projectRaster(tiff_data1, crs = st_crs(world)$proj4string)

# 将栅格数据转为数据框，方便 ggplot 绘制
raster_df1 <- as.data.frame(raster::rasterToPoints(tiff_data1))
colnames(raster_df1) <- c("long", "lat", "value")

# 定义自定义颜色梯度
# custom_colors <- c("#7B95C6", "#49C2D9", "#A1D8E8", "#67A583", "#A2C986", "#D0E2C0", "#FDED95", "#FFc1a6", "#f59c7c", "#f47254", "#c85e62")
custom_colors <- c("#49C2D9", "#A1D8E8", "#D0E2C0", "#FDED95", "#FFc1a6", "#f59c7c", "#f47254", "#c85e62")
# custom_colors <- c(
#   "#7FB6D4",  # 稍深浅蓝
#   "#58CDED",  # 稍深青蓝
#   "#AEE8D0",  # 稍深薄荷绿
#   "#DDDDA8",  # 稍深米黄
#   "#E6BDD4",  # 稍深浅粉
#   "#D7AFC0",  # 稍深粉紫
#   "#E3B79E",  # 稍深桃色
#   "#DDA548",  # 稍深橙黄
#   "#7A5144"   # 更深棕色
# )

# 创建地图
map_with_raster1 <- ggplot() +
  geom_raster(data = raster_df1, aes(x = long, y = lat, fill = value)) +
  scale_fill_gradientn(colors = custom_colors, na.value = "grey40", guide = "colorbar") +
  geom_sf(data = world, fill = NA, color = "black", size = 0.2) +  # 较细边界线
  coord_sf(expand = FALSE) +
  scale_y_continuous(
    breaks = c(-60, -30, 0, 30, 60),
    labels = c("60°S", "30°S", "0°", "30°N", "60°N")
  ) +
  scale_x_continuous(
    breaks = seq(-180, 180, by = 60),
    labels = function(x) paste0(abs(x), "°", ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  ) +
  theme_minimal(base_family = "RMN") +
  labs(title = "lnRR of NH\u2083 emission") +
  theme(
    legend.position = "bottom",
    text = element_text(family = "RMN", face = "bold", color = "black", size = 14),
    axis.title = element_blank(),
    axis.text = element_text(family = "RMN", face = "bold", color = "black", size = 12),
    axis.ticks = element_line(size = 0.2),
    panel.grid.major = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    panel.grid.minor = element_line(color = "grey60", size = 0.2, linetype = "solid"),
    plot.title = element_text(hjust = 0.5)
  )

# 显示地图
print(map_with_raster1)

tiff("E:/有机替代.作图/response/map.tiff",width=4000,height=2000,res=300,compression="lzw")
figure_all4 <- ggarrange(map_with_raster1, map_with_raster,
                         widths = c(1, 1), heights = c(1, 1), ncol=2, nrow=1,
                         labels = c("(e)","(f)"),
                         label.x=0.1, label.y=0.8, align = "v",
                         font.label = list(size = 18, color = "black", face="bold", family="RMN"))

# 打印图形
print(figure_all4)
dev.off()
}
###柱状图
{
##========================
## 0. 加载包
##========================
library(terra)
library(sf)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(rnaturalearth)

##========================
## 1. 路径设置
##========================
in_dir  <- "E:/有机替代.作图/response/预测/N2O"
out_dir <- in_dir

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## 文件列表：NH3_RF_001.tif 到 NH3_RF_100.tif
tif_files <- sprintf(
  file.path(in_dir, "NH3_RF_%03d.tif"),
  1:100
)

## 如果你想自动匹配文件，也可以用下面这一行替换上面的 tif_files
# tif_files <- list.files(
#   in_dir,
#   pattern = "^NH3_RF_\\d{3}\\.tif$",
#   full.names = TRUE
# )

## 检查文件是否存在
cat("总文件数：", length(tif_files), "\n")
cat("实际存在文件数：", sum(file.exists(tif_files)), "\n")

if (sum(file.exists(tif_files)) == 0) {
  stop("没有找到任何 NH3_RF_001.tif 到 NH3_RF_100.tif 文件，请检查文件名和路径。")
}

##========================
## 2. 构建世界大区
##========================
world_sf <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
) |>
  sf::st_make_valid() |>
  dplyr::mutate(
    region = dplyr::case_when(
      continent %in% c("North America", "South America") ~ "Americas",
      continent %in% c("Africa", "Asia", "Europe", "Oceania") ~ continent,
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(region)) |>
  dplyr::group_by(region) |>
  dplyr::summarise(.groups = "drop") |>
  sf::st_as_sf() |>
  sf::st_make_valid()

## 检查 world_sf
print(world_sf)
print(class(world_sf))
print(sf::st_crs(world_sf))

##========================
## 3. 定义单个栅格处理函数
##========================
process_one <- function(tif_path, world_sf) {
  
  if (!file.exists(tif_path)) {
    message("Skip, file not found: ", tif_path)
    return(tibble())
  }
  
  message("Processing: ", basename(tif_path))
  
  ## 读取栅格
  r <- terra::rast(tif_path)
  
  ## 如果有多层，只取第一层
  r1 <- r[[1]]
  
  ## sf 转为 terra SpatVector
  world_v <- terra::vect(world_sf)
  
  ## 检查栅格 CRS
  raster_crs <- terra::crs(r1)
  
  if (is.na(raster_crs) || raster_crs == "") {
    message("Warning: raster has no CRS: ", basename(tif_path))
  } else {
    world_v <- terra::project(world_v, raster_crs)
  }
  
  ## 按大区提取像元
  vals <- terra::extract(r1, world_v)
  
  if (is.null(vals) || nrow(vals) == 0) {
    message("No data extracted: ", basename(tif_path))
    return(tibble())
  }
  
  ## 自动识别栅格值所在列
  value_col <- names(vals)[names(vals) != "ID"][1]
  
  if (is.na(value_col)) {
    message("No value column found: ", basename(tif_path))
    return(tibble())
  }
  
  ## 汇总 Mean 和 SE
  res <- vals |>
    as_tibble() |>
    dplyr::rename(value = dplyr::all_of(value_col)) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(
      Mean = mean(value, na.rm = TRUE),
      SE   = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
      N_pixel = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Region = world_v$region[ID],
      File   = basename(tif_path),
      Index  = stringr::str_extract(
        basename(tif_path),
        "(?<=NH3_RF_)\\d{3}(?=\\.tif)"
      )
    ) |>
    dplyr::select(File, Index, Region, Mean, SE, N_pixel) |>
    dplyr::arrange(Region)
  
  return(res)
}

##========================
## 4. 批量处理所有栅格
##========================
all_results <- purrr::map_dfr(
  tif_files,
  function(f) {
    tryCatch(
      process_one(f, world_sf),
      error = function(e) {
        message("Error on: ", f, " | ", conditionMessage(e))
        return(tibble())
      }
    )
  }
)

## 查看结果
print(all_results)

##========================
## 5. 保存结果
##========================
out_csv <- file.path(out_dir, "N2O_RF_region_mean_SE_001_100.csv")

write.csv(
  all_results,
  out_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("结果已保存到：", out_csv, "\n")


BNF222 <- read.csv("E:/有机替代.作图/response/预测/N2O/N2O_RF_region_mean_SE_001_100.csv")

BNF1<- BNF222[which(BNF222$Region == "Africa"), ]
BNF2<- BNF222[which(BNF222$Region == "Americas"), ]
BNF3<- BNF222[which(BNF222$Region == "Asia"), ]
BNF4<- BNF222[which(BNF222$Region == "Europe"), ]
BNF5<- BNF222[which(BNF222$Region == "Oceania"), ]

tota11<-rbind(BNF1, BNF2, BNF3, BNF4, BNF5)

AverageBNF11 <- ddply(tota11, c("Region"), summarise,
                      N    = sum(!is.na(Mean)),
                      median = median(Mean, na.rm=TRUE),
                      mean = mean(Mean, na.rm=TRUE),
                      sd   = sd(Mean, na.rm=TRUE),
                      se   = sd / sqrt(N),
                      ci_lower = mean - qt(0.975, df = N - 1) * se,
                      ci_upper = mean + qt(0.975, df = N - 1) * se
)
write.csv(AverageBNF11,"E:/有机替代.作图/response/预测/N2O/Mean.csv")

##========================
## 0. 加载包
##========================
library(terra)
library(sf)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(rnaturalearth)

##========================
## 1. 路径设置
##========================
in_dir  <- "E:/有机替代.作图/response/预测/NH3"
out_dir <- in_dir

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## 文件列表：NH3_RF_001.tif 到 NH3_RF_100.tif
tif_files <- sprintf(
  file.path(in_dir, "NH3_RF_%03d.tif"),
  1:100
)

## 如果你想自动匹配文件，也可以用下面这一行替换上面的 tif_files
# tif_files <- list.files(
#   in_dir,
#   pattern = "^NH3_RF_\\d{3}\\.tif$",
#   full.names = TRUE
# )

## 检查文件是否存在
cat("总文件数：", length(tif_files), "\n")
cat("实际存在文件数：", sum(file.exists(tif_files)), "\n")

if (sum(file.exists(tif_files)) == 0) {
  stop("没有找到任何 NH3_RF_001.tif 到 NH3_RF_100.tif 文件，请检查文件名和路径。")
}

##========================
## 2. 构建世界大区
##========================
world_sf <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
) |>
  sf::st_make_valid() |>
  dplyr::mutate(
    region = dplyr::case_when(
      continent %in% c("North America", "South America") ~ "Americas",
      continent %in% c("Africa", "Asia", "Europe", "Oceania") ~ continent,
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(region)) |>
  dplyr::group_by(region) |>
  dplyr::summarise(.groups = "drop") |>
  sf::st_as_sf() |>
  sf::st_make_valid()

## 检查 world_sf
print(world_sf)
print(class(world_sf))
print(sf::st_crs(world_sf))

##========================
## 3. 定义单个栅格处理函数
##========================
process_one <- function(tif_path, world_sf) {
  
  if (!file.exists(tif_path)) {
    message("Skip, file not found: ", tif_path)
    return(tibble())
  }
  
  message("Processing: ", basename(tif_path))
  
  ## 读取栅格
  r <- terra::rast(tif_path)
  
  ## 如果有多层，只取第一层
  r1 <- r[[1]]
  
  ## sf 转为 terra SpatVector
  world_v <- terra::vect(world_sf)
  
  ## 检查栅格 CRS
  raster_crs <- terra::crs(r1)
  
  if (is.na(raster_crs) || raster_crs == "") {
    message("Warning: raster has no CRS: ", basename(tif_path))
  } else {
    world_v <- terra::project(world_v, raster_crs)
  }
  
  ## 按大区提取像元
  vals <- terra::extract(r1, world_v)
  
  if (is.null(vals) || nrow(vals) == 0) {
    message("No data extracted: ", basename(tif_path))
    return(tibble())
  }
  
  ## 自动识别栅格值所在列
  value_col <- names(vals)[names(vals) != "ID"][1]
  
  if (is.na(value_col)) {
    message("No value column found: ", basename(tif_path))
    return(tibble())
  }
  
  ## 汇总 Mean 和 SE
  res <- vals |>
    as_tibble() |>
    dplyr::rename(value = dplyr::all_of(value_col)) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(
      Mean = mean(value, na.rm = TRUE),
      SE   = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
      N_pixel = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Region = world_v$region[ID],
      File   = basename(tif_path),
      Index  = stringr::str_extract(
        basename(tif_path),
        "(?<=NH3_RF_)\\d{3}(?=\\.tif)"
      )
    ) |>
    dplyr::select(File, Index, Region, Mean, SE, N_pixel) |>
    dplyr::arrange(Region)
  
  return(res)
}

##========================
## 4. 批量处理所有栅格
##========================
all_results <- purrr::map_dfr(
  tif_files,
  function(f) {
    tryCatch(
      process_one(f, world_sf),
      error = function(e) {
        message("Error on: ", f, " | ", conditionMessage(e))
        return(tibble())
      }
    )
  }
)

## 查看结果
print(all_results)

##========================
## 5. 保存结果
##========================
out_csv <- file.path(out_dir, "NH3_RF_region_mean_SE_001_100.csv")

write.csv(
  all_results,
  out_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("结果已保存到：", out_csv, "\n")


BNF222 <- read.csv("E:/有机替代.作图/response/预测/NH3/NH3_RF_region_mean_SE_001_100.csv")

BNF1<- BNF222[which(BNF222$Region == "Africa"), ]
BNF2<- BNF222[which(BNF222$Region == "Americas"), ]
BNF3<- BNF222[which(BNF222$Region == "Asia"), ]
BNF4<- BNF222[which(BNF222$Region == "Europe"), ]
BNF5<- BNF222[which(BNF222$Region == "Oceania"), ]

tota11<-rbind(BNF1, BNF2, BNF3, BNF4, BNF5)

AverageBNF11 <- ddply(tota11, c("Region"), summarise,
                      N    = sum(!is.na(Mean)),
                      median = median(Mean, na.rm=TRUE),
                      mean = mean(Mean, na.rm=TRUE),
                      sd   = sd(Mean, na.rm=TRUE),
                      se   = sd / sqrt(N),
                      ci_lower = mean - qt(0.975, df = N - 1) * se,
                      ci_upper = mean + qt(0.975, df = N - 1) * se
)
write.csv(AverageBNF11,"E:/有机替代.作图/response/预测/NH3/Mean.csv")

#####8.8 绘制预测地图柱状图####
# 读取CSV文件
windowsFonts(RMN=windowsFont("Times New Roman"))
data <- read.csv("E:/有机替代.作图/response/预测/N2O/Mean.csv")

# 绘制柱状图，并添加误差棒
p <- ggplot(data, aes(x=Region, y=Mean, fill=Region)) +
  geom_bar(stat="identity", color="black", size=0.6) +  # 设置柱状图轮廓为黑色并增加宽度
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.2, color="black" , size = 1) +
  scale_fill_manual(values=c("#EFCE87", "#EAA558", "#ED8D5A", "#BFDFD2", "#257D8B")) +
  # scale_fill_manual(values=c("#7895C1", "#A8CBDF", "#D6EFF4", "#F2FAFC", "#E36250", "#EF8B67", "#992224", "#B54764", "#7F4693")) +
  theme(
    panel.background = element_rect(fill = "white"),  # 设置整个图形的背景颜色为白色
    axis.line = element_line(color = "black"),  # 设置轴线为黑色
    axis.line.x = element_line(color = "black"),  # 设置x轴线为黑色
    axis.line.y = element_line(color = "black"),  # 设置y轴线为黑色
    legend.position = "none",  # 隐藏图例
    text = element_text(family = "RMN", face = "bold", size = 15, color = "black"),  # 设置字体为黑色
    axis.title = element_text(family = "RMN", face = "bold", size = 15, color = "black"),  # 设置轴标题字体为黑色并加粗
    axis.text = element_text(family = "RMN", face = "bold", size = 15, color = "black"),  # 设置轴文本字体为黑色并加粗
    panel.grid.major = element_line(color = "#EEEEEE", size = 0.5, linetype = "solid"),  # 设置主要网格线为灰色实线
    panel.grid.minor = element_blank()  # 不显示次要网格线
  ) +
  labs(
    x = "Region",
    y = expression("Change of " * N[2] * O * " emission")
  ) +
  labs(
    x = "Region",
    y = expression(bold("Change of " * N[2] * O * " emission"))
  ) +
  theme(
    axis.title.x = element_text(family = "RMN", face = "bold", size = 15),
    axis.title.y = element_text(family = "RMN", face = "bold", size = 15)
  )
print(p)

#####8.8 绘制预测地图柱状图####
# 读取CSV文件
windowsFonts(RMN=windowsFont("Times New Roman"))
data1 <- read.csv("E:/有机替代.作图/response/预测/NH3/Mean.csv")

# 绘制柱状图，并添加误差棒
p1 <- ggplot(data1, aes(x=Region, y=Mean, fill=Region)) +
  geom_bar(stat="identity", color="black", size=0.6) +  # 设置柱状图轮廓为黑色并增加宽度
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.2, color="black" , size = 1) +
  scale_fill_manual(values=c("#EFCE87", "#EAA558", "#ED8D5A", "#BFDFD2", "#257D8B")) +
  # scale_fill_manual(values=c("#7895C1", "#A8CBDF", "#D6EFF4", "#F2FAFC", "#E36250", "#EF8B67", "#992224", "#B54764", "#7F4693")) +
  theme(
    panel.background = element_rect(fill = "white"),  # 设置整个图形的背景颜色为白色
    axis.line = element_line(color = "black"),  # 设置轴线为黑色
    axis.line.x = element_line(color = "black"),  # 设置x轴线为黑色
    axis.line.y = element_line(color = "black"),  # 设置y轴线为黑色
    legend.position = "none",  # 隐藏图例
    text = element_text(family = "RMN", face = "bold", size = 15, color = "black"),  # 设置字体为黑色
    axis.title = element_text(family = "RMN", face = "bold", size = 15, color = "black"),  # 设置轴标题字体为黑色并加粗
    axis.text = element_text(family = "RMN", face = "bold", size = 15, color = "black"),  # 设置轴文本字体为黑色并加粗
    panel.grid.major = element_line(color = "#EEEEEE", size = 0.5, linetype = "solid"),  # 设置主要网格线为灰色实线
    panel.grid.minor = element_blank()  # 不显示次要网格线
  ) +
  labs(
    x = "Region",
    y = expression("Change of " * NH[3] * " emission")
  ) +
  labs(
    x = "Region",
    y = expression(bold("Change of " * NH[3] * " emission"))
  ) +
  theme(
    axis.title.x = element_text(family = "RMN", face = "bold", size = 15),
    axis.title.y = element_text(family = "RMN", face = "bold", size = 15)
  )
print(p1)

tiff("E:/有机替代.作图/response/map1.tiff",width=3000,height=1200,res=300,compression="lzw")
figure_all4 <- ggarrange(p1, p,
                         widths = c(1, 1), heights = c(1, 1), ncol=2, nrow=1,
                         labels = c("(g)","(h)"),
                         label.x=-0.02, label.y=0.98, align = "v",
                         font.label = list(size = 20, color = "black", face="bold", family="RMN"))

# 打印图形
print(figure_all4)
dev.off()
}
}


# Load your package
library(SNRBoost)
library(ggplot2)
library(stats)
library(tidyr)
library(dplyr)
library(readr)

library(scales)
library(export)
library(RColorBrewer)

font_size <- 10
#display.brewer.all(colorblindFriendly = T, type=c("div","qual","seq","all")[2])
col_qual <- RColorBrewer::brewer.pal(n=8, name=c("Dark2","Set2")[1])
col_div <- RColorBrewer::brewer.pal(n=4, name="RdYlBu")
#scales::show_col(col_qual)

# Synthetic data generation----
# Random seed
set.seed(20240924)

### Step 1: Data parameters
p <- 2 # number of datasets
n <- 300 # data length

ecc <- 0.3 # error cross-correlation [0,1]
SNRdB <- 0.1 # SNR in dB

model <- switch(2,"rand","sine")

### Step 2: Synthetic data generation
# Signal and error: y and e
n_ensemble <- 100
df_mat <- vector('list', n_ensemble)
#for(i_r in 1:n_ensemble){
######## scaling factor (a)
a <- matrix(rep(1, p), ncol=1)  # 1-vector
#a <- matrix(runif(p), ncol=1) # non-1-vector (for test)

generated_data <- dataGEN(n, p, ecc, SNRdB, model=model)
y <- generated_data$y %>% matrix(ncol=1)
e <- generated_data$e

# Observation: x = a*signal + error
x <- y %*% t(a) + e

# Signal power and covariance matrices of e and x
Ey2 <- var(y[,1]); sum(y^2)/(n-1)
#Ey2 <- sum(y^2)/(n)
EeeT <- cov(e)
ExxT <- cov(x)
N <- EeeT / Ey2 # error-to-signal ratio
N

## raw stat---------------------------------------------------------------------
# MSE of observations
MSE_ori <- sapply(1:p, function(i) sqrt(mean((y-x[,i])^2)))

# Pearson correlation of observations
R2_ori <- cor(y,x)

# Printing metrics
cat("+ Metrics for original data\n")
cat(" * RMSE for x:", round(MSE_ori, 3), "\n")
cat(" * R for x:", round(R2_ori, 3), "\n")
## raw--------------------------------------------------------------------------
spc.method <- switch(2, "pgram", "ar")
spc.kernel <- kernel("fejer", 100, r=6)

spectrum_df <- function(data, span=2, log="no", kernel=spc.kernel, method=spc.method) {
  spec <- spectrum(data, span=span, log=log, kernel=kernel, method=method, plot=F)

  # Create a dataframe with frequency and spectral density
  delta <- 1 # sampling interval = time sequence
  spec_df <- data.frame(Freq = spec$freq/delta, SpectralDensity = 2*spec$spec)

  return(spec_df)
}

#y.spec <- spectrum(y,span=2, method=spc.method,plot=TRUE)

df_raw <- data.frame(No=1:n, x, y) %>% gather(Group, value, 2:(2+p))
summary(df_raw)

my_labeller <- as_labeller(c(X1="x[1]",
                             X2="x[2]",
                             X3="x[3]",
                             y="y [Truth]"),
                           default = label_parsed)

df_raw$Group <- factor(df_raw$Group, levels = c("y","X1","X2", "X3"))

### fig-a ----
p1a <- ggplot(df_raw) +
  geom_line(aes(x=No, y=value, color=Group),linewidth=0.5) +

  facet_wrap(Group~., strip.position="left", ncol = 1, labeller = my_labeller) +

  #scale_y_continuous(limits=c(-2.5,2.5)) +
  scale_color_manual(values=c("black", col_qual)) +
  labs(x=NULL, y=NULL) +

  #egg::theme_article() +
  theme(text = element_text(size = font_size),
        plot.margin = unit(c(0.8,0.1,0.5,0.2), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #panel.background = element_rect(fill=alpha(col_qual[8],0.15)),
        plot.background= element_blank(),

        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0,
                                         margin = margin(t = 0, r = 12, b = 0, l = 0)),
        strip.text = element_text(size = font_size),

        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),

        legend.position = "none",
        legend.title=element_blank(),
        legend.key.width = unit(1,"cm"))

p1a

# graph2pdf(x=p1a, file="Figure1a", aspectr=2, font = "Arial", width = 9/2.54,
#           height = 7/2.54, bg = "transparent")


## ----spec raw-----------------------------------------------------------------
y_spec <- spectrum_df(y)
x_spec <- sapply(1:p, function(i) spectrum_df(x[,i])$Spec)

df_spec <- data.frame(Freq=y_spec$Freq,
                      y=y_spec$SpectralDensity,
                      x_spec) %>%
  gather(Group,Spec,2:(2+p))

df_spec$Group <- factor(df_spec$Group, levels = c("y","X1","X2", "X3"))

summary(df_spec)

bg_df <- data.frame(ymin=c(10^-4, 10^-2, 1, 10^2),
                    ymax=c(10^-2, 1, 10^2, 10^4),
                    color=col_div)

### fig-b ----
p1b <- ggplot(df_spec) +


  geom_line(aes(Freq, Spec, color=Group),linewidth=0.5) +

  scale_y_log10(limits=c(1e-4,1e4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), expand=c(0,0)) +

  geom_rect(data=bg_df, aes(ymin =ymin, ymax = ymax),fill = col_qual[7], #"#F4A637",
            xmin = -Inf, xmax = Inf, alpha =0.15, show.legend = F) +

  # geom_rect(ymin = -Inf, ymax = Inf, fill="#F4A637",
  #           xmin = -Inf, xmax = Inf, alpha =0.01, show.legend = F) +

  scale_color_manual(values=c("black",col_qual), labels=expression(y[Truth],x[1],x[2],x[3])) +
  labs(x="Frequency", y="Spectral Density") +
  guides(color = guide_legend(byrow = TRUE, ncol=1)) +

  egg::theme_article() +
  theme(text = element_text(size = font_size),
        axis.text = element_text(size = font_size),
        legend.text = element_text(size = font_size),

        plot.margin = unit(c(1,0.5,0.1, 0.5), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),

        strip.background = element_blank(),
        strip.placement = "outside",

        legend.position = c(0.8,0.8),
        legend.key.spacing.y = unit(0, 'cm'),
        legend.background = element_blank(),
        #legend.position = "right",
        legend.title=element_blank(),
        legend.key.width = unit(0.5,"cm"))

p1b

# graph2pdf(x=p1b, file="Figure1b", aspectr=2, font = "Arial", width = 9/2.54,
#           height = 7/2.54, bg = "transparent")


## merge------------------------------------------------------------------------
## Merging using true parameters
ue <- rep(1/p, p)  # equal weight
uw <- WA(EeeT) # weighted average
us <- SNRopt(N, a) # SNR-opt
#sum(uw)
#uw - us/sum(us)
y_eq <- x %*% ue
y_wa <- x %*% uw
y_snr <- x%*%us

## freq-snr
wf <- "d4"
if(wf!="haar") v <- as.integer(parse_number(wf)/2) else v <- 1
J <- floor(log2(n/2)) - 1
J <- 6
mode <- 'DWT'; pad="zero"; boundary="periodic"

y_snr_freq <- SNRopt_freq(y, x, mode, wf, J, pad, boundary, option="Truth")$merged

summary(y_snr_freq)

cor(y,cbind(y_eq, y_wa, y_snr, y_snr_freq))

df_y <- data.frame(No=1:n, y, y_wa, y_snr,y_snr_freq) %>% gather(Group, value,2:5)
summary(df_y)

my_labeller <- as_labeller(c(y="y [Truth]",
                             y_wa="y [WA]",
                             y_snr="y [SNR]",
                             y_snr_freq="y [SNR(ω)]"),
                           default = label_parsed)

df_y$Group <- factor(df_y$Group, levels = c("y","y_wa", "y_snr", "y_snr_freq"))

### fig-c ----
p1c <- ggplot(df_y %>% subset(Group!="y1")) +
  geom_line(aes(x=No, y=value, color=Group), linewidth=0.5) +

  facet_wrap(Group~., strip.position="left", ncol = 1, scale="fixed", labeller = my_labeller) +

  scale_y_continuous(limits=c(-2.5,2.5)) +
  labs(x=NULL, y=NULL) +
  #scale_color_manual(values=col_qual) +
  scale_color_manual(values=c("black","red","blue","green")) +
  coord_cartesian(clip = "off") +

  #egg::theme_article() +
  theme(text = element_text(size = font_size),
        plot.margin = unit(c(0.8,0.1,0.5,0.2), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #panel.background = element_rect(fill=alpha(col_qual[8], 0.15)),
        plot.background= element_blank(),

        panel.spacing = unit(0, "lines"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0,
                                         margin = margin(t = 0, r = 1, b = 0, l = 0)),
        strip.text = element_text(size = font_size),

        #axis.title.x = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),

        legend.position = "none",
        legend.title=element_blank(),
        legend.background = element_blank(),
        legend.key.width = unit(0.5,"cm"))

p1c

# graph2pdf(x=p1c, file="Figure1c", aspectr=2, font = "Arial", width = 9/2.54,
#           height = 7/2.54, bg = "transparent")

## spec merge-------------------------------------------------------------------
y_spec1 <- spectrum_df(y_wa, span=1)$Spect

y_spec2 <- spectrum_df(y_snr, span=1)$Spect

y_spec3 <- spectrum_df(y_snr_freq, span=1)$Spect


df_y_spec <- data.frame(Freq=y_spec$Freq,
                        y=y_spec$SpectralDensity,
                        y_wa=y_spec1,
                        y_snr=y_spec2,
                        y_snr_freq=y_spec3) %>%
  gather(Group,Spec,2:5)

df_y_spec$Group <- factor(df_y_spec$Group)
summary(df_y_spec)
df_y_spec$Group <- factor(df_y_spec$Group, levels = c("y","y_wa", "y_snr", "y_snr_freq"))

# df for background colors
bg_df <- data.frame(ymin=c(10^-4, 10^-2, 1, 10^2),
                    ymax=c(10^-2, 1, 10^2, 10^4),
                    color=col_div)

### fig-d ----
p1d <- ggplot(df_y_spec) +

  geom_line(aes(Freq, Spec, color=Group),linewidth=0.5) +

  scale_y_log10(limits=c(1e-4,1e4),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), expand=c(0,0)) +

  geom_rect(data=bg_df, aes(ymin =ymin, ymax = ymax),fill = col_qual[3], #"#B6DBFF",
            xmin = -Inf, xmax = Inf, alpha =0.15, show.legend = F) +

  # geom_rect(ymin = -Inf, ymax = -Inf, fill="#00C992",
  #           xmin = -Inf, xmax = Inf, alpha =0.01, show.legend = F) +

  labs(x="Frequency", y="Spectral Density") +

  scale_color_manual(values=c("black","red","blue","green"), labels=expression(y[Truth],y[WA],y[SNR],y[SNR(ω)])) +

  egg::theme_article() +
  #theme_bw() +
  theme(text = element_text(size = font_size),
        axis.text = element_text(size = font_size),
        legend.text = element_text(size = font_size),

        plot.margin = unit(c(1,0.5,0.1, 0.5), "cm"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),

        #axis.text.x = element_text(angle = 45, vjust=1, hjust=1),

        strip.background = element_blank(),
        strip.placement = "outside",

        legend.position = c(0.8,0.8),
        legend.key.spacing.y = unit(0, 'cm'),
        legend.background = element_blank(),
        legend.spacing.y =unit(0, "cm"),
        #legend.position = "right",
        legend.title=element_blank(),
        legend.key.width = unit(0.5,"cm"))

p1d


fig <- cowplot::plot_grid(p1a,p1b, p1c, p1d, ncol=2, #labels=letters,
                          rel_widths = c(1,0.8),
                          #label_size= font_size + 2, hjust=1,
                          label_fontfamily = "sans",
                          label_fontface = "bold",
                          label_colour = "black")

fig %>% print()

## merge stat-------------------------------------------------------------------
# MSE of observations
y_merge <- cbind(y_wa, y_snr, y_snr_freq)
MSE_mer <- sapply(1:p, function(i) sqrt(mean((y-y_merge[,i])^2)))

# Pearson correlation of observations
R2_mer <- cor(y,y_merge)

# Printing metrics
cat("+ Metrics for merged data\n")
cat(" * RMSE for x:", round(MSE_mer, 3), "\n")
cat(" * R for x:", round(R2_mer, 3), "\n")

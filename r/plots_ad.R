## ## ######################### SOME BASIC PLOTS

## ### Plots for one specific model
here::i_am("r/plots_ad.R")
library(here)

## ## Do we want to save or/as well as display plotting output?
save.flag <- TRUE

### Crucial Specifications
# mod.to.plot <- 6
mod.to.plot <- 4 ## Describes the type of spline used to model incidence
data.to.plot <- 6 ## Describes the exclusion criteria used in compiling the data.
data.id <- 1 ### data.rm index. data.id = 1 considers all data.
yrs.to.plot <- 2012:2022 ## Which years do we want to feature in the plots

### Loading results and diagnosis data
load(here("data/postproc_ad.RData"))
load(here("data/test_data_ad.RData"))

# models.ind <- sapply(out, function(x) as.numeric(x$model))
# data.ind <- sapply(out, function(x) as.numeric(x$data))

models.ind <- 4
data.ind <- 1

ind.to.plot <- which(data.ind == data.id & models.ind == mod.to.plot)

res <- out
rm(out)

library(tidyverse)
library(ggplot2)
library(RColorBrewer)

yrs <- 1995:(1994 + ncol(res$infs))
qts <- seq(1995, (1994 + length(yrs) + 0.75), .25)

min.y <- min(yrs.to.plot)
ind.y <- which(yrs %in% yrs.to.plot)
ind.q <- ((4 * ind.y[1]) - 3):(4 * ind.y[length(ind.y)])
qts.to.plot <- qts[ind.q]
nq <- length(qts.to.plot)
plt.y <- substr(yrs.to.plot, 3, 4)
m <- models.ind[ind.to.plot]

# add named values
data.95$num.quarters <- data.95$nquar # PK
data.95$num.years <- data.95$num.quarters / 4
data.95$AIDS.diagnoses <- data.95$AIDS
data.95$HIV.diagnoses <- data.95$HIV
data.95$CD4.cell.proportions <- data.95$CD4

### Data manipultation
yr.dx <- tapply(rowSums(data.95$HIV.diagnoses + data.95$AIDS.diagnoses), rep(1:length(yrs), each = 4), sum)
yr.15.dx <- tapply(rowSums((data.95$HIV.diagnoses + data.95$AIDS.diagnoses)[, 1:10]), rep(1:length(yrs), each = 4), sum)
yr.25.dx <- tapply(rowSums((data.95$HIV.diagnoses + data.95$AIDS.diagnoses)[, 11:20]), rep(1:length(yrs), each = 4), sum)
yr.35.dx <- tapply(rowSums((data.95$HIV.diagnoses + data.95$AIDS.diagnoses)[, 21:30]), rep(1:length(yrs), each = 4), sum)
yr.45.dx <- tapply(rowSums((data.95$HIV.diagnoses + data.95$AIDS.diagnoses)[, 31:52]), rep(1:length(yrs), each = 4), sum)

## ### Grepping number of MSM
## x <- which(msm.number.data$Age=="Total" & msm.number.data$Region=="Total")
## n.MSM <- msm.number.data[x,"Median"]
## names(n.MSM) <-  msm.number.data[x,"Year"]

##### INFECTIONS ALL AGES

### Data frame infections
df1 <- data.frame(
  "yrs" = factor(rep(plt.y, 3), levels = plt.y, ordered = TRUE),
  "infs" = c(res$infs[3, ind.y], res$infs[1, ind.y], res$infs[4, ind.y]),
  "grp" = rep(c("m", "l", "u"), each = length(yrs.to.plot)),
  "grp1" = rep(c("m", "c", "c"), each = length(yrs.to.plot))
)

### Data frame dx
df2 <- data.frame("yrs" = plt.y, "dx" = yr.dx[ind.y], "type" = rep(1, length(yrs.to.plot)))

### Data frame
p <- df1 |>
  select(-grp1) |>
  pivot_wider(names_from = grp, values_from = infs) |>
  ggplot() +
  aes(x = factor(yrs), y = m, fill = "Estimated posterior mean number of new infections") +
  # geom_line(data=df1, aes(x=factor(yrs), y=infs, group=grp, linetype=grp1), col="cornflowerblue") +
  geom_col() +
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.5) +
  labs(x = "Year", y = "Expected number of new HIV infections") +
  ylim(0, 4000) +
  theme_minimal() +
  # scale_linetype_manual(name="", breaks=c("m","c"), values=c("dashed", "solid"), labels=c("Estimated posterior mean number of new infections", "Posterior 95% credible intervals")) +
  geom_point(data = df2, aes(x = factor(yrs), y = dx, shape = as.factor(type)), size = 2, color = "black") +
  scale_shape_manual(name = "", values = c(4), labels = c("Number of new observed diagnoses")) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_fill_discrete(name = "") +
  ggtitle("Expected infections over time") +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE))

if (save.flag) ggsave(filename = here(paste0("figures/Infs_mod", m, "_dat_", data.to.plot, "_", min.y, ".png")), plot = p, width = 12) else print(p)

###### INFECTIONS AGE CLASS PLOTS

### Data frame infections
age.strings <- c("15-24", "25-34", "35-44", "45+")
grp.t <- rep(c("m", "l", "u"), each = length(yrs.to.plot))
grp1.t <- rep(c("m", "c", "c"), each = length(yrs.to.plot))
ac.t <- rep(age.strings, each = length(yrs.to.plot) * 3)

df1 <- data.frame(
  "yrs" = factor(rep(plt.y, 12), levels = plt.y, ordered = TRUE),
  "infs" = c(
    res$infs.1524[3, ind.y], res$infs.1524[1, ind.y], res$infs.1524[4, ind.y],
    res$infs.2534[3, ind.y], res$infs.2534[1, ind.y], res$infs.2534[4, ind.y],
    res$infs.3544[3, ind.y], res$infs.3544[1, ind.y], res$infs.3544[4, ind.y],
    res$infs.45[3, ind.y], res$infs.45[1, ind.y], res$infs.45[4, ind.y]
  ),
  "grp" = rep(grp.t, 4),
  "grp1" = rep(grp1.t, 4),
  "ac" = ac.t
)

df2 <- data.frame(
  "yrss" = substr(yrs.to.plot, 3, 4),
  "dx" = c(yr.15.dx[ind.y], yr.25.dx[ind.y], yr.35.dx[ind.y], yr.45.dx[ind.y]),
  "ac" = rep(age.strings, each = length(yrs.to.plot))
)


p <- df1 |>
  select(-grp1) |>
  pivot_wider(names_from = grp, values_from = infs) |>
  ggplot() +
  aes(x = factor(yrs), y = m, fill = "Estimated posterior mean number of new infections") +
  # geom_line(data=df1, aes(x=factor(yrs), y=infs, group=grp, linetype=grp1), col="cornflowerblue") +
  geom_col() +
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.5) +
  # p <- ggplot() +
  #  geom_line(data=df1, aes(x=factor(yrs), y=infs, group=grp, linetype=grp1), col="cornflowerblue") +
  # facet_wrap( ~ ac) +
  labs(x = "Year", y = "Expected number of new HIV infections") +
  ylim(0, 1500) +
  theme_minimal() +
  #  scale_linetype_manual(name="", breaks=c("m","c"), values=c("dashed", "solid"), labels=c("Estimated posterior mean number of new infections", "Posterior 95% credible intervals")) +
  geom_point(data = df2, aes(x = factor(yrss), y = dx, shape = as.factor(1)), size = 2, color = "black") +
  facet_wrap(~ac) +
  scale_shape_manual(name = "", values = c(4), labels = c("Number of new observed diagnoses")) +
  scale_fill_discrete(name = "") +
  ggtitle("Expected infections over time, by age group") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold"))
# guides(linetype=guide_legend(nrow=2,byrow=TRUE))

if (save.flag) ggsave(filename = here(paste0("figures/InfsAge_mod", m, "_dat_", data.to.plot, "_", min.y, ".png")), plot = p, width = 12) else print(p)

## ### Peak distribution
## df1 <- data.frame("pk"=res$peak+2006 )

## p <- ggplot() + geom_bar(data=df1, aes(x=pk, y = ..prop.., fill = factor(..x..), stats = "count")) +
##   labs(x = "Year", y= "Empirical frequency of peaks in the posterior - Data 2016") +
##   guides(fill=guide_legend(title="Year"))+
##   scale_y_continuous(labels=scales::percent) +
##   scale_x_continuous(breaks=2007:2016,
##                      labels=substr(2007:2016,3,4),
##                      limits=c(2007,2016.5)) +
##   theme(plot.title = element_text(hjust = 0.5, face="bold"))+
##   ggtitle("Frequency of peak in new infections")

## if(save.flag) ggsave(filename=paste0(save.dir,"PeakInfs.png"), plot=p) else print(p)

## ### by age
## df1 <- data.frame("pk"=c(res$peak.1524+2006,res$peak.2534+2006,
##                          res$peak.3544+2006,res$peak.45+2006),
##                   "ac"=rep(c("15-24","25-34","35-44","45+"), each=length(res$peak.1524)))

## p <- ggplot() + geom_bar(data=df1, aes(x=pk, y = ..prop.., fill = factor(..x..), stats = "count")) +
##   facet_wrap( ~ ac) +
##   scale_x_continuous(breaks=2007:2016,
##                      labels=substr(2007:2016,3,4)) +
##   labs(x = "Year", y= "Empirical distribution of peaks in the posterior") +
##   guides(fill=guide_legend(title=""))+
##   scale_y_continuous(labels=scales::percent) +
##   ggtitle("Frequency of peak in incidence in last 10 years by age class") +
##   theme(plot.title = element_text(hjust = 0.5, face="bold")) +
##   guides(fill=FALSE)
## if(save.flag) ggsave(filename=paste0(save.dir,"PeakInfsAge.png"), plot=p) else print(p)

#### Probability of decrease in the last 3 years
df1 <- data.frame("pk" = as.factor(res$trend.2))

p <- ggplot() +
  geom_bar(data = df1, aes(x = pk, y = (..count..) / sum(..count..), fill = factor(..x..))) +
  labs(x = "Trend", y = "Frequency of trends in the last 2 years in posterior distribution") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Posterior distribution trend in last 2 years") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(fill = "none")
if (save.flag) ggsave(filename = here("figures/TrendInfs.png"), plot = p, width = 12) else print(p)

### By age group
df1 <- data.frame(
  "pk" = c(
    res$trend.2.1524, res$trend.2.2534,
    res$trend.2.3544, res$trend.2.45
  ),
  "ac" = rep(c("15-24", "25-34", "35-44", "45+"), each = length(res$trend.2.45))
)

p <- ggplot() +
  geom_bar(data = df1, aes(x = pk, y = ..prop.., fill = factor(..x..), group = ac)) +
  facet_wrap(~ac) +
  labs(x = "Trend", y = "Frequency of trends in the last 2 years in posterior distribution") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = "none") +
  ggtitle("Posterior distribution trend in last 2 years by age group") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
if (save.flag) ggsave(filename = here("figures/TrendInfsAgeClass.png"), plot = p, width = 12) else print(p)

## ### Incidence plots
## df1 <- data.frame("yrs"= substr(rep(2012:2016, 3),3,4), "infs"=c(res$inc[3,],res$inc[1,],res$inc[4,]),
##                   "grp"=rep(c("m","l","u"), each=5),"grp1"=rep(c("m","c","c"), each=5))

## p <- ggplot() +
##   geom_line(data=df1, aes(x=factor(yrs), y=infs, group=grp, linetype=grp1), col="cornflowerblue") +
##   labs(x = "Year", y= "Expected HIV incidence") +
##   scale_linetype_manual(name="", breaks=c("m","c"), values=c("dashed", "solid"), labels=c("Estimated incidence of new infections", "Posterior 95% credible intervals")) +
##   theme(legend.position="bottom", legend.box = "vertical",legend.box.just = "left", plot.title = element_text(hjust = 0.5, face="bold")) +
##   guides(linetype=guide_legend(nrow=2,byrow=TRUE)) +
##   ggtitle("Expected Incidence")

## if(save.flag) ggsave(filename=paste0(save.dir,"Incidence.png"), plot=p) else print(p)

## ### By age group
## yrs.to.plot <- 2012:2016  ### years of incidence
## grp.t <- rep(c("m","l","u"), each=length(yrs.to.plot))
## grp1.t <- rep(c("m","c","c"), each=length(yrs.to.plot))
## ac.t <- rep(c("15-34","35-44","45+"), each=length(yrs.to.plot)*3)

## df1 <- data.frame("yrs"=rep(substr(yrs.to.plot,3,4),9),
##                   "infs"=c(res$inc.1534[3,],res$inc.1534[1,],res$inc.1534[4,],
##                            res$inc.3544[3,],res$inc.3544[1,],res$inc.3544[4,],
##                            res$inc.45[3,],res$inc.45[1,],res$inc.45[4,]),
##                   "grp"=rep(grp.t,3),
##                   "grp1"=rep(grp1.t,3),
##                   "ac"=ac.t)


## p <- ggplot() +
##   geom_line(data=df1, aes(x=factor(yrs), y=infs, group=grp, linetype=grp1), col="cornflowerblue") +
##   facet_wrap( ~ ac) +
##   labs(x = "Year", y= "Expected HIV incidence") +
##   scale_linetype_manual(name="", breaks=c("m","c"), values=c("dashed", "solid"), labels=c("Estimated expected HIV incidence", "Posterior 95% credible intervals")) +
##   theme(legend.position="bottom", legend.box = "vertical",legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"),plot.title = element_text(hjust = 0.5, face="bold")) +
##   guides(linetype=guide_legend(nrow=2,byrow=TRUE))+
##   ggtitle("Expected Incidence by Age")

## if(save.flag) ggsave(filename=paste0(save.dir,"IncAge.png"), plot=p) else print(p)

#### TABLE REQUESTED
yr.ind <- seq(3, length(yrs) * 4, 4)

df <- data.frame(
  yrs,
  yr.dx,
  res$infs[3, ],
  res$prev.st1[3, yr.ind] + res$prev.st2[3, yr.ind],
  res$prev.st3[3, yr.ind] + res$prev.st4[3, yr.ind],
  yr.15.dx, res$infs.1524[3, ], res$prev.st1.1524[3, yr.ind] + res$prev.st2.1524[3, yr.ind], res$prev.st3.1524[3, yr.ind] + res$prev.st4.1524[3, yr.ind],
  yr.25.dx, res$infs.2534[3, ], res$prev.st1.2534[3, yr.ind] + res$prev.st2.2534[3, yr.ind], res$prev.st3.2534[3, yr.ind] + res$prev.st4.2534[3, yr.ind],
  yr.35.dx, res$infs.3544[3, ], res$prev.st1.3544[3, yr.ind] + res$prev.st2.3544[3, yr.ind], res$prev.st3.3544[3, yr.ind] + res$prev.st4.3544[3, yr.ind],
  yr.45.dx, res$infs.45[3, ], res$prev.st1.45[3, yr.ind] + res$prev.st2.45[3, yr.ind], res$prev.st3.45[3, yr.ind] + res$prev.st4.45[3, yr.ind]
)

df <- data.frame(
  yrs,
  yr.dx,
  res$infs[3, ], res$infs[1, ], res$infs[4, ],
  yr.15.dx, res$infs.1524[3, ], res$infs.1524[1, ], res$infs.1524[4, ],
  yr.25.dx, res$infs.2534[3, ], res$infs.2534[1, ], res$infs.2534[4, ],
  yr.35.dx, res$infs.3544[3, ], res$infs.3544[1, ], res$infs.3544[4, ],
  yr.45.dx, res$infs.45[3, ], res$infs.45[1, ], res$infs.45[4, ]
)

colnames(df) <- c(
  "Year", "Diagnoses", "Infections", "Infections (2.5% CrI)", "Infections (97.5% CrI)",
  "Diagnoses_15-24", "Infections_15-24", "Infections_15-24 (2.5% CrI)", "Infections_15-24 (97.5% CrI)",
  "Diagnoses_25-34", "Infections_25-34", "Infections_25-34 (2.5% CrI)", "Infections_25-34 (97.5% CrI)",
  "Diagnoses_35-44", "Infections_35-44", "Infections_35-44 (2.5% CrI)", "Infections_35-44 (97.5% CrI)",
  "Diagnoses_45+", "Infections_45+", "Infections_45+ (2.5% CrI)", "Infections_45+ (97.5% CrI)"
)

df <- df[df$Year %in% yrs.to.plot, ]
df <- round(df)

library(gridExtra)
png(file = here(paste0("figures/Table_mod", m, "_dat_", data.to.plot, "_", min.y, ".png")), height = 10, width = 20)
p <- tableGrob(t(df))
grid.arrange(p)
dev.off()

library(tidyverse)
df <- tribble(
  ~"Year", ~"Age", ~"Diagnoses", ~"Infections", ~"Infections (2.5% CrI)", ~"Infections (97.5% CrI)",
  yrs, "all", c(yr.dx), res$infs[3, ], res$infs[1, ], res$infs[4, ]
) |>
  bind_rows(tribble(
    ~"Year", ~"Age", ~"Diagnoses", ~"Infections", ~"Infections (2.5% CrI)", ~"Infections (97.5% CrI)",
    yrs, "15-24", c(yr.15.dx), res$infs.1524[3, ], res$infs.1524[1, ], res$infs.1524[4, ]
  )) |>
  bind_rows(tribble(
    ~"Year", ~"Age", ~"Diagnoses", ~"Infections", ~"Infections (2.5% CrI)", ~"Infections (97.5% CrI)",
    yrs, "25-34", c(yr.25.dx), res$infs.2534[3, ], res$infs.2534[1, ], res$infs.2534[4, ]
  )) |>
  bind_rows(tribble(
    ~"Year", ~"Age", ~"Diagnoses", ~"Infections", ~"Infections (2.5% CrI)", ~"Infections (97.5% CrI)",
    yrs, "35-44", c(yr.35.dx), res$infs.3544[3, ], res$infs.3544[1, ], res$infs.3544[4, ]
  )) |>
  bind_rows(tribble(
    ~"Year", ~"Age", ~"Diagnoses", ~"Infections", ~"Infections (2.5% CrI)", ~"Infections (97.5% CrI)",
    yrs, "45+", c(yr.45.dx), res$infs.45[3, ], res$infs.45[1, ], res$infs.45[4, ]
  ))

library(readr)
df <- tibble(df)
write_csv(df, file = here(paste0("figures/Table_mod", m, "_dat_", data.to.plot, "_", min.y, ".csv")))


#### UNDIAGNOSED PREVALENCE
## Data frame with prevalence
df1 <- data.frame(
  "qts" = rep(qts.to.plot, 3),
  "infs" = c(res$prev[3, ind.q], res$prev[1, ind.q], res$prev[4, ind.q]),
  "grp" = rep(c("m", "l", "u"), each = length(ind.q)),
  "grp1" = rep(c("m", "c", "c"), each = length(ind.q))
)

p <- ggplot() +
  geom_line(data = df1, aes(x = qts, y = infs, group = grp, linetype = grp1), col = "cornflowerblue") +
  labs(x = "Year", y = "Expected number of undiagnosed individuals living with HIV") +
  ylim(0, max(df1$infs)) +
  theme_minimal() +
  scale_linetype_manual(name = "", breaks = c("m", "c"), values = c("dashed", "solid"), labels = c("Estimated expected number of individuals undiagnosed with HIV", "Posterior 95% credible intervals")) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm")) +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Undiagnosed Prevalence") +
  ## scale_x_discrete() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
    # axis.line = element_line(colour = "darkblue",size = 1.5, linetype = "solid")
  )

if (save.flag) ggsave(filename = here(paste0("figures/UDiag_", min.y, ".png")), plot = p, width = 12) else print(p)


### BY AGE
grp.t <- rep(c("m", "l", "u"), each = length(qts.to.plot))
grp1.t <- rep(c("m", "c", "c"), each = length(qts.to.plot))
ac.t <- rep(age.strings, each = length(qts.to.plot) * 3)
na <- length(age.strings)
df1 <- data.frame(
  "qts" = rep(qts.to.plot, 3 * na),
  "infs" = c(
    res$prev.1524[3, ind.q], res$prev.1524[1, ind.q], res$prev.1524[4, ind.q],
    res$prev.2534[3, ind.q], res$prev.2534[1, ind.q], res$prev.2534[4, ind.q],
    res$prev.3544[3, ind.q], res$prev.3544[1, ind.q], res$prev.3544[4, ind.q],
    res$prev.45[3, ind.q], res$prev.45[1, ind.q], res$prev.45[4, ind.q]
  ),
  "grp" = rep(grp.t, na),
  "grp1" = rep(grp1.t, na),
  "ac" = ac.t
)

p <- ggplot() +
  geom_line(data = df1, aes(x = qts, y = infs, group = grp, linetype = grp1), col = "cornflowerblue") +
  facet_wrap(~ac) +
  labs(x = "Year", y = "Expected Number of undiagnosed individuals living with HIV") +
  ylim(0, max(df1$infs)) +
  theme_minimal() +
  scale_linetype_manual(name = "", breaks = c("m", "c"), values = c("dashed", "solid"), labels = c("Estimated expected number of individuals undiagnosed with HIV", "Posterior 95% credible intervals")) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold")) + # axis.line = element_line(colour = "darkblue",size = 1.5, linetype = "solid")) +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Undiagnosed Prevalence, stratified by age group")

if (save.flag) ggsave(filename = here(paste0("figures/UdiagAge_", min.y, ".png")), plot = p, width = 12) else print(p)

### In percentage
df1 <- data.frame(
  "qts" = rep(qts.to.plot, 3 * na),
  "infs" = c(
    res$prev.1524[3, ind.q],
    res$prev.2534[3, ind.q],
    res$prev.3544[3, ind.q],
    res$prev.45[3, ind.q]
  ),
  "grp" = rep(grp.t, na),
  "grp1" = rep(grp1.t, na),
  "ac" = ac.t
)

p <- ggplot() +
  geom_bar(data = df1, aes(x = qts, weight = infs, fill = ac), position = "fill") +
  labs(x = "Year", y = "Proportion of Undiagnosed") +
  scale_fill_manual(name = "", values = c("skyblue", "royalblue", "blue", "navy")) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Proportion of undiagnosed individuals, stratified by age group, over time")

if (save.flag) ggsave(filename = here(paste0("figures/UdiagAgePercent_", min.y, ".png")), plot = p, width = 12) else print(p)


### BY STATE
state.strings <- c("CD4 > 500", "350 < CD4 < 500", "350 < CD4 < 200", "CD4 < 200")
ns <- length(state.strings)
ac.t <- rep(state.strings, each = length(qts.to.plot) * 3)

df1 <- data.frame(
  "qts" = rep(qts.to.plot, 3 * ns),
  "infs" = c(
    res$prev.st1[3, ind.q], res$prev.st1[1, ind.q], res$prev.st1[4, ind.q],
    res$prev.st2[3, ind.q], res$prev.st2[1, ind.q], res$prev.st2[4, ind.q],
    res$prev.st3[3, ind.q], res$prev.st3[1, ind.q], res$prev.st3[4, ind.q],
    res$prev.st4[3, ind.q], res$prev.st4[1, ind.q], res$prev.st4[4, ind.q]
  ),
  "grp" = rep(grp.t, ns),
  "grp1" = rep(grp1.t, ns),
  "ac" = ac.t
)

### Trick to sort order of facets
df1$ac <- factor(df1$ac, levels = state.strings)

p <- ggplot() +
  geom_line(data = df1, aes(x = qts, y = infs, group = grp, linetype = grp1), col = "cornflowerblue") +
  facet_wrap(~ac) +
  labs(x = "Year", y = "Expected Number of undiagnosed individual living with HIV") +
  ylim(0, 5000) +
  scale_linetype_manual(name = "", breaks = c("m", "c"), values = c("dashed", "solid"), labels = c("Estimated posterior mean number of HIV undiagnosed individuals", "Posterior 95% credible intervals")) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Undiagnosed Prevalence by State")

if (save.flag) ggsave(filename = here(paste0("figures/UdiagState_", min.y, ".png")), plot = p, width = 12) else print(p)

### As a percentage
df1 <- data.frame(
  "qts" = rep(qts.to.plot, ns),
  "infs" = c(
    res$prev.st1[3, ind.q],
    res$prev.st2[3, ind.q],
    res$prev.st3[3, ind.q],
    res$prev.st4[3, ind.q]
  ),
  "ac" = rep(state.strings, each = length(qts.to.plot))
)

### Trick to sort order of groups
df1$ac <- factor(df1$ac, levels = state.strings)

p <- ggplot() +
  geom_bar(data = df1, aes(x = qts, weight = infs, fill = ac), position = "fill") +
  labs(x = "Year", y = "Proportion of Undiagnosed") +
  scale_fill_manual(name = "", values = c("skyblue", "royalblue", "blue", "navy")) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Proportion of undiagnosed individuals, stratified by state, over time")

if (save.flag) ggsave(filename = here(paste0("figures/UdiagStatePercent_", min.y, ".png")), plot = p, width = 12) else print(p)

## By state and age - NEW
df1 <- data.frame(
  "qts" = rep(qts.to.plot, 3 * ns * na),
  "infs" = c(
    res$prev.st1.1524[3, ind.q], res$prev.st1.1524[1, ind.q], res$prev.st1.1524[4, ind.q],
    res$prev.st1.2534[3, ind.q], res$prev.st1.2534[1, ind.q], res$prev.st1.2534[4, ind.q],
    res$prev.st1.3544[3, ind.q], res$prev.st1.3544[1, ind.q], res$prev.st1.3544[4, ind.q],
    res$prev.st1.45[3, ind.q], res$prev.st1.45[1, ind.q], res$prev.st1.45[4, ind.q],
    res$prev.st2.1524[3, ind.q], res$prev.st2.1524[1, ind.q], res$prev.st2.1524[4, ind.q],
    res$prev.st2.2534[3, ind.q], res$prev.st2.2534[1, ind.q], res$prev.st2.2534[4, ind.q],
    res$prev.st2.3544[3, ind.q], res$prev.st2.3544[1, ind.q], res$prev.st2.3544[4, ind.q],
    res$prev.st2.45[3, ind.q], res$prev.st2.45[1, ind.q], res$prev.st2.45[4, ind.q],
    res$prev.st3.1524[3, ind.q], res$prev.st3.1524[1, ind.q], res$prev.st3.1524[4, ind.q],
    res$prev.st3.2534[3, ind.q], res$prev.st3.2534[1, ind.q], res$prev.st3.2534[4, ind.q],
    res$prev.st3.3544[3, ind.q], res$prev.st3.3544[1, ind.q], res$prev.st3.3544[4, ind.q],
    res$prev.st3.45[3, ind.q], res$prev.st3.45[1, ind.q], res$prev.st3.45[4, ind.q],
    res$prev.st4.1524[3, ind.q], res$prev.st4.1524[1, ind.q], res$prev.st4.1524[4, ind.q],
    res$prev.st4.2534[3, ind.q], res$prev.st4.2534[1, ind.q], res$prev.st4.2534[4, ind.q],
    res$prev.st4.3544[3, ind.q], res$prev.st4.3544[1, ind.q], res$prev.st4.3544[4, ind.q],
    res$prev.st4.45[3, ind.q], res$prev.st4.45[1, ind.q], res$prev.st4.45[4, ind.q]
  ),
  "ac" = rep(rep(age.strings, each = 3 * nq), ns),
  "st" = rep(state.strings, each = nq * 3 * na),
  "col" = rep(1:ns, each = nq * 3 * na),
  "lty" = rep(rep(c("m", "l", "u"), each = nq), ns * na)
)

### Trick to sort order of groups
df1$ac <- factor(df1$ac, levels = age.strings)
df1$st <- factor(df1$st, levels = state.strings)

p <- ggplot() +
  geom_line(data = df1, aes(x = qts, y = infs, colour = ac, group = interaction(ac, lty), linetype = (lty != "m"))) +
  facet_wrap(~st) +
  labs(x = "Year", y = "Expected Number of undiagnosed individual living with HIV") +
  ylim(0, 2000) +
  scale_linetype_manual(name = "", breaks = c(FALSE, TRUE), values = c("solid", "dashed"), labels = c("Estimated posterior mean number of HIV undiagnosed individuals", "Posterior 95% credible intervals")) +
  scale_color_discrete(name = "", breaks = unique(df1$ac), labels = age.strings) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Undiagnosed Prevalence by Age and State")

if (save.flag) ggsave(filename = here(paste0("figures/UdiagStateAge_", min.y, ".png")), plot = p, width = 12) else print(p)

## Plot the other way-around
p <- ggplot() +
  geom_line(data = df1, aes(x = qts, y = infs, colour = st, group = interaction(st, lty), linetype = (lty != "m"))) +
  facet_wrap(~ac) +
  labs(x = "Year", y = "Expected Number of undiagnosed individual living with HIV") +
  scale_linetype_manual(name = "", breaks = c(FALSE, TRUE), values = c("solid", "dashed"), labels = c("Estimated posterior mean number of HIV undiagnosed individuals", "Posterior 95% credible intervals")) +
  scale_color_discrete(name = "", breaks = unique(df1$st), labels = state.strings) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Undiagnosed Prevalence by Age and State")

if (save.flag) ggsave(filename = here(paste0("figures/UdiagStateAge2_", min.y, ".png")), plot = p, width = 12) else print(p)

### As a percentage
df.sub <- df1[df1$lty == "m", ]

p <- ggplot() +
  geom_bar(data = df.sub, aes(x = qts, weight = infs, fill = st), position = "fill") +
  facet_wrap(~ac) +
  labs(x = "Year", y = "Proportion of Undiagnosed") +
  scale_fill_manual(name = "", values = c("skyblue", "royalblue", "blue", "navy")) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Proportion of undiagnosed individuals, stratified by state and age group, over time")

if (save.flag) ggsave(filename = here(paste0("figures/UdiagStateAgePercent_", min.y, ".png")), plot = p, width = 12) else print(p)

## And the other way around
p <- ggplot() +
  geom_bar(data = df.sub, aes(x = qts, weight = infs, fill = ac), position = "fill") +
  facet_wrap(~st) +
  labs(x = "Year", y = "Proportion of Undiagnosed") +
  scale_fill_manual(name = "", values = c("skyblue", "royalblue", "blue", "navy")) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(-0.25, "cm"), plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Proportion of undiagnosed individuals, stratified by state and age group, over time")

if (save.flag) ggsave(filename = here(paste0("figures/UdiagAgeStatePercent_", min.y, ".png")), plot = p, width = 12) else print(p)


require(abind)
diag.plot <- function(k, d1, d2, d3, d4, plt.idx = ind.q, qnames = qts.to.plot) {
  rcols <- rainbow(4)
  d <- abind(d1, d2, d3, d4, along = 0)
  y.max <- max(d)
  plot(qnames, d[1, 2, plt.idx], type = "l", col = rcols[1], main = paste("Diagnosis probabilities, ", k), ylim = c(0, y.max), xlab = "Year", ylab = "Quarterly probability")
  lines(qnames, d[1, 1, plt.idx], lty = 3, col = rcols)
  lines(qnames, d[1, 4, plt.idx], lty = 3, col = rcols)
  for (a in 2:4) {
    lines(qnames, d[a, 2, plt.idx], col = rcols[a])
    lines(qnames, d[a, 1, plt.idx], lty = 3, col = rcols[a])
    lines(qnames, d[a, 4, plt.idx], lty = 3, col = rcols[a])
  }
}

png(here("figures/DiagProbs.png"), width = 800)
par(mfrow = c(2, 2))
attach(res)
if ("d1.q.y" %in% names(res)) { ## Diagnosis probabilities governed by logistic regression with age as a four-level factor
  diag.plot("CD4>500", d1.q.y, d1.q.my, d1.q.mo, d1.q.o)
  diag.plot("CD4 350-500", d2.q.y, d2.q.my, d2.q.mo, d2.q.o)
  diag.plot("CD4 200-350", d3.q.y, d3.q.my, d3.q.mo, d3.q.o)
  diag.plot("CD4<200", d4.q.y, d4.q.my, d4.q.mo, d4.q.o)
  legend("bottomright", legend = c("15-24 years", "25-34 years", "35-44 years", ">44 years"), lty = 1, col = rainbow(4), cex = 0.75, bty = "n")
} else { ## Diagnosis probabilities modelled as a spline, with continuous response over age
  pda <- c(20, 30, 40, 53) - 14
  diag.plot("CD4>500", d1.q[, , pda[1]], d1.q[, , pda[2]], d1.q[, , pda[3]], d1.q[, , pda[4]])
  diag.plot("CD4 350-500", d2.q[, , pda[1]], d2.q[, , pda[2]], d2.q[, , pda[3]], d2.q[, , pda[4]])
  diag.plot("CD4 200-350", d3.q[, , pda[1]], d3.q[, , pda[2]], d3.q[, , pda[3]], d3.q[, , pda[4]])
  diag.plot("CD4<200", d4.q[, , pda[1]], d4.q[, , pda[2]], d4.q[, , pda[3]], d4.q[, , pda[4]])
  legend("bottomright", legend = paste(pda + 14, "years"), lty = 1, col = rainbow(4), cex = 0.75, bty = "n")
}
detach(res)
dev.off()

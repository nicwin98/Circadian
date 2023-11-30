library(tidyverse)
library(cosinoRmixedeffects)
library(lme4)
library(emmeans)
library(limma)
library(cosinor)
library(circacompare)
library(patchwork)
library(cosinor2)

df <- read.csv2("xxx")

cpep <- df %>% filter(type %in% c("c-pep"))

glu <- df %>% filter(type %in% c("glu"))

gcg <- df %>% filter(type %in% c("gcg"))

glp1 <- df %>% filter(type %in% c("glp1"))

gip <- df %>% filter(type %in% c("gip"))


#Cosinor analyses:

#C-peptid#

#Create model
cpepmod <- cosinor.lm(value ~ time(tid), data=cpep, period=24)
summary(cpepmod)

#Peak time
(cscpep <- circa_single(x=cpep, col_time = "tid", col_outcome = "value", period=24)$summary)

#test the overall significance of the cosinormodel. rythm detection test/zero amplitude test (cornelissen,2014):
cosinor.detect(cpepmod)


timeintervals <- seq(0, 24, 0.1)
newdata <- data.frame(time = timeintervals, rrr = cos(2 * pi * timeintervals / 24),
                      sss = sin(2 * pi * timeintervals / 24))


cpeppred <- as.data.frame(predict(cpepmod$fit, newdata, interval = "confidence"))
cpeppred$time <- timeintervals

#Plot
(cpepplot <- ggplot(cpeppred, aes(x=time, y=fit)) +
    annotate("rect", xmin=14, xmax=23, ymin=0, ymax=max(cpep$value, na.rm=TRUE), fill = "Grey", alpha = 0.4) +
    geom_segment(aes(x = 0.5, y = 0, xend = 0.5, yend = -(max(cpep$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 4, y = 0, xend = 4, yend = -(max(cpep$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 10, y = 0, xend = 10, yend = -(max(cpep$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_point(data=cpep, aes(x=tid, y=value), shape = 21, fill = "Grey", color = "Black", alpha = 0.9) +
    geom_boxplot(data=cpep, aes(x=tid, y=value, group=tid), width = 0.4, position= position_nudge(x=-1), outlier.shape = NA) +
    geom_line(color = "Red", linewidth = 1) +
    geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill = "Red", alpha = 0.2) +
    theme_classic() +
    theme(axis.text.x=element_text(size=rel(2))) +
    theme(axis.text.y=element_text(size=rel(2))) +
    theme(axis.title=element_text(size=16,face="bold")) +
    theme(title=element_text(size=18, face='bold')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(-1.5,24.5), breaks = seq(0,24,3), labels = c("09", "12", "15", "18", "21", "00", "03", "06", "09")) +
    scale_y_continuous(expand=c(0.01,0)) +
    labs(x="Time (Hours)", y="C-peptide conc. (pmol/L)", title = "C-Peptide\n p-value <0.0001"))


#Glucose#

#Create model
glumod <- cosinor.lm(value ~ time(tid), data=glu, period=24)
summary(glumod)

#test the overall significance of the cosinormodel. rythm detection test/zero amplitude test (cornelissen,2014):
cosinor.detect(glumod)

#Peak time
(csglu <- circa_single(x=glu, col_time = "tid", col_outcome = "value", period=24)$summary)

timeintervals <- seq(0, 24, 0.1)
newdata <- data.frame(time = timeintervals, rrr = cos(2 * pi * timeintervals / 24),
                      sss = sin(2 * pi * timeintervals / 24))

glupred <- as.data.frame(predict(glumod$fit, newdata, interval = "confidence"))
glupred$time <- timeintervals

#Plot
(gluplot <- ggplot(glupred, aes(x=time, y=fit)) +
    annotate("rect", xmin=14, xmax=23, ymin=2, ymax=max(glu$value, na.rm=TRUE), fill = "Grey", alpha = 0.4) +
    geom_segment(aes(x = 0.5, y = 2, xend = 0.5, yend = 2-(max(glu$value, na.rm=TRUE)*0.031)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 4, y = 2, xend = 4, yend = 2-(max(glu$value, na.rm=TRUE)*0.031)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 10, y = 2, xend = 10, yend = 2-(max(glu$value, na.rm=TRUE)*0.031)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_point(data=glu, aes(x=tid, y=value),  shape = 21, fill = "Grey", color = "Black", alpha = 0.9) +
    geom_boxplot(data=glu, aes(x=tid, y=value, group=tid), width = 0.4, position= position_nudge(x=-1), outlier.shape = NA) +
    geom_line(color = "Red", linewidth = 1) +
    geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill = "Red", alpha = 0.2) +
    theme_classic() +
    theme(axis.text.x=element_text(size=rel(2))) +
    theme(axis.text.y=element_text(size=rel(2))) +
    theme(axis.title=element_text(size=16,face="bold")) +
    theme(title=element_text(size=18, face='bold')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(-1.5,24.5), breaks = seq(0,24,3), labels = c("09", "12", "15", "18", "21", "00", "03", "06", "09")) +
    scale_y_continuous(expand=c(0.01,0)) +
    labs(x="Time (Hours)", y="Glucose conc. (mmol/L)", title = "Glucose\n p-value <0.0001"))


#Glucagon#

#Create model
gcgmod <- cosinor.lm(value ~ time(tid), data=gcg, period=24)
summary(gcgmod)

#test the overall significance of the cosinormodel. rythm detection test/zero amplitude test (cornelissen,2014):
cosinor.detect(gcgmod)

#Peak time
(csgcg <- circa_single(x=gcg, col_time = "tid", col_outcome = "value", period=24)$summary)

timeintervals <- seq(0, 24, 0.1)
newdata <- data.frame(time = timeintervals, rrr = cos(2 * pi * timeintervals / 24),
                      sss = sin(2 * pi * timeintervals / 24))

gcgpred <- as.data.frame(predict(gcgmod$fit, newdata, interval = "confidence"))
gcgpred$time <- timeintervals

#Plot
(gcgplot <- ggplot(gcgpred, aes(x=time, y=fit)) +
    annotate("rect", xmin=14, xmax=23, ymin=0, ymax=max(gcg$value, na.rm=TRUE), fill = "Grey", alpha = 0.4) +
    geom_segment(aes(x = 0.5, y = 0, xend = 0.5, yend = -(max(gcg$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 4, y = 0, xend = 4, yend = -(max(gcg$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 10, y = 0, xend = 10, yend = -(max(gcg$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_point(data=gcg, aes(x=tid, y=value), shape = 21, fill = "Grey", color = "Black", alpha = 0.9) +
    geom_boxplot(data=gcg, aes(x=tid, y=value, group=tid), width = 0.4, position= position_nudge(x=-1), outlier.shape = NA) +
    geom_line(color = "Red", linewidth = 1) +
    geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill = "Red", alpha = 0.2) +
    theme_classic() +
    theme(axis.text.x=element_text(size=rel(2))) +
    theme(axis.text.y=element_text(size=rel(2))) +
    theme(axis.title=element_text(size=16,face="bold")) +
    theme(title=element_text(size=18, face='bold')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(-1.5,24.5), breaks = seq(0,24,3), labels = c("09", "12", "15", "18", "21", "00", "03", "06", "09")) +
    scale_y_continuous(expand=c(0.01,0)) +
    labs(x="Time (Hours)", y="Glucagon conc. (pmol/L)", title = "Glucagon\n p-value <0.0001"))


#GLP-1#

#Create model
glp1mod <- cosinor.lm(value ~ time(tid), data=glp1, period=24)
summary(glp1mod)

#test the overall significance of the cosinormodel. rythm detection test/zero amplitude test (cornelissen,2014):
cosinor.detect(glp1mod)

#Peak time
(glp1cg <- circa_single(x=glp1, col_time = "tid", col_outcome = "value", period=24)$summary)

timeintervals <- seq(0, 24, 0.1)
newdata <- data.frame(time = timeintervals, rrr = cos(2 * pi * timeintervals / 24),
                      sss = sin(2 * pi * timeintervals / 24))

glp1pred <- as.data.frame(predict(glp1mod$fit, newdata, interval = "confidence"))
glp1pred$time <- timeintervals

#Plot
(glp1plot <- ggplot(glp1pred, aes(x=time, y=fit)) +
    annotate("rect", xmin=14, xmax=23, ymin=0, ymax=max(glp1$value, na.rm=TRUE), fill = "Grey", alpha = 0.4) +
    geom_segment(aes(x = 0.5, y = 0, xend = 0.5, yend = -(max(glp1$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 4, y = 0, xend = 4, yend = -(max(glp1$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 10, y = 0, xend = 10, yend = -(max(glp1$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_point(data=glp1, aes(x=tid, y=value), shape = 21, fill = "Grey", color = "Black", alpha = 0.9) +
    geom_boxplot(data=glp1, aes(x=tid, y=value, group=tid), width = 0.4, position= position_nudge(x=-1), outlier.shape = NA) +
    geom_line(color = "Red", linewidth = 1) +
    geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill = "Red", alpha = 0.2) +
    theme_classic() +
    theme(axis.text.x=element_text(size=rel(2))) +
    theme(axis.text.y=element_text(size=rel(2))) +
    theme(axis.title=element_text(size=16,face="bold")) +
    theme(title=element_text(size=18, face='bold')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(-1.5,24.5), breaks = seq(0,24,3), labels = c("09", "12", "15", "18", "21", "00", "03", "06", "09")) +
    scale_y_continuous(expand=c(0.01,0)) +
    labs(x="Time (Hours)", y="GLP-1 conc. (pmol/L)", title = "GLP-1\n p-value <0.0001"))

#GIP#

#Create model
gipmod <- cosinor.lm(value ~ time(tid), data=gip, period=24)
summary(gipmod)

#test the overall significance of the cosinormodel. rythm detection test/zero amplitude test (cornelissen,2014):
cosinor.detect(gipmod)

#Peak time
(gipcg <- circa_single(x=gip, col_time = "tid", col_outcome = "value", period=24)$summary)

timeintervals <- seq(0, 24, 0.1)
newdata <- data.frame(time = timeintervals, rrr = cos(2 * pi * timeintervals / 24),
                      sss = sin(2 * pi * timeintervals / 24))

gippred <- as.data.frame(predict(gipmod$fit, newdata, interval = "confidence"))
gippred$time <- timeintervals

#Plot
(gipplot <- ggplot(gippred, aes(x=time, y=fit)) +
    annotate("rect", xmin=14, xmax=23, ymin=0, ymax=max(gip$value, na.rm=TRUE), fill = "Grey", alpha = 0.4) +
    geom_segment(aes(x = 0.5, y = 0, xend = 0.5, yend = -(max(gip$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 4, y = 0, xend = 4, yend = -(max(gip$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = 10, y = 0, xend = 10, yend = -(max(gip$value, na.rm=TRUE)*0.041)),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_point(data=gip, aes(x=tid, y=value), shape = 21, fill = "Grey", color = "Black", alpha = 0.9) +
    geom_boxplot(data=gip, aes(x=tid, y=value, group=tid), width = 0.4, position= position_nudge(x=-1), outlier.shape = NA) +
    geom_line(color = "Red", linewidth = 1) +
    geom_ribbon(aes(x=time, ymin=lwr, ymax=upr), fill = "Red", alpha = 0.2) +
    theme_classic() +
    theme(axis.text.x=element_text(size=rel(2))) +
    theme(axis.text.y=element_text(size=rel(2))) +
    theme(axis.title=element_text(size=16,face="bold")) +
    theme(title=element_text(size=18, face='bold')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(-1.5,24.5), breaks = seq(0,24,3), labels = c("09", "12", "15", "18", "21", "00", "03", "06", "09")) +
    scale_y_continuous(expand=c(0.01,0)) +
    labs(x="Time (Hours)", y="GIP conc. (pmol/L)", title = "GIP\n p-value <0.0001"))




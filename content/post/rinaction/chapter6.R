
# chapter6. Basic graphs


# This chapter covers
#   Bar, box, and dot plots
#   Pie and fan charts
#   Histograms and kernel density plots


# Remove most objects from the working environment
rm(list = ls()) 
options(stringsAsFactors = F)


# 6.1. Bar plots
library(vcd)
head(Arthritis)
dim(Arthritis)

# 6.1.1. Simple bar plots
counts <- table(Arthritis$Improved)
counts

# code listing 6.1
barplot(counts,
        main = "Simple Bar Plot",
        xlab = "Improvement",
        ylab = "Frequency")

barplot(counts,
        main = "Horizontal Bar Plot",
        ylab = "Improvement",
        xlab = "Frequency",
        horiz = T)

# 6.1.2. Stacked and grouped bar plots
counts <- table(Arthritis$Improved, Arthritis$Treatment)
counts

# code listing 6.2
barplot(counts,
        main = "Stacked Bar Plot",
        xlab = "Treatment",
        ylab = "Frequency",
        col = c("red", "yellow", "green"),
        legend=rownames(counts))

barplot(counts,
        main = "Grouped Bar Plot",
        xlab = "Treatment",
        ylab = "Frequency",
        col = c("red", "yellow", "green"),
        legend=rownames(counts), beside = T)

# 6.1.3. Mean bar plots
# code listing 6.3
states <- data.frame(state.region, state.x77)
head(states)

means <- aggregate(states$Illiteracy, by=list(state.region), FUN=mean)
means

means <- means[order(means$x),]

barplot(means$x, names.arg = means$Group.1, title("Mean Illiteracy Rate"))

# 6.1.4. Tweaking bar plots
# code listing 6.4
par(mar=c(5,8,4,2))
par(las=2)
counts <- table(Arthritis$Improved)

barplot(counts,
        main="Treatment Outcome",
        horiz=TRUE, cex.names=0.8,
        names.arg=c("No Improvement", "Some Improvement",
                    "Marked Improvement"))

# 6.1.5. Spinograms spine()

library(vcd)
attach(Arthritis)
counts <- table(Treatment, Improved)
spine(counts, main="Spinogram Example")
detach(Arthritis)

# 6.2. Pie charts












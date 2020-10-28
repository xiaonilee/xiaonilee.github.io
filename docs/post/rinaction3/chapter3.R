
# Chapter3. Getting Started with Graphs


# This chapter covers
#   Creating and saving graphs
#   Customizing symbols, lines, colors, and axes
#   Annotating with text and titles
#   Controlling a graph’s dimensions
#   Combining multiple graphs into one

# Remove most objects from the working environment
rm(list = ls()) 
options(stringsAsFactors = F)



# Section 3.3 Graphical Parameters
?par 
# equal
help("par")

## Section 3.3.2 Colors
library(RColorBrewer)
n <- 7
mycolors <- brewer.pal(n, "Set1")
barplot(rep(1,n), col = mycolors)

n <- 10
mycolors <- rainbow(n)
pie(rep(1,n),labels = mycolors,col = mycolors)
mygrays <- gray(0:n/n)
pie(rep(1,n), labels = mygrays, col = mygrays)

# context
quartzFonts()
names(postscriptFonts())

## Section 3.3.4. Graph and margin dimensions
## code listing 3.1
dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
drugB <- c(15, 18, 25, 31, 40)

plot(dose, drugA, type = "b")

opar <- par(no.readonly = TRUE)

# par(): control global graphs
par(pin=c(2,3))
par(lwd=2, cex=1.5)
par(cex.axis=.75, font.axis=3)

# plot(): control current settings
plot(dose, drugA, type = "b", pch=19, lty=2, col="red")
plot(dose, drugB, type = "b", pch=23, lty=6, col="blue", bg="green")

# At end of plotting, reset to previous settings:
par(opar)


# code listing 3.3
dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
drugB <- c(15, 18, 25, 31, 40)

opar <- par(no.readonly = TRUE)
opar

par(lwd=2, cex=1.5, font.lab=2)

plot(dose, drugA, type = "b", 
     pch=15, lty=1, col="red",ylim = c(0,60),
     main = "Drug A vs. Drug B",
     xlab = "Drug Dosage", ylab = "Drug Response")

lines(dose, drugB, type = "b", 
      pch=17, lty=2, col="blue")
# Add Reference lines
abline(h=c(30), lwd=1.5, lty=2, col="gray")

# Add Minor Tick Marks
library(Hmisc)
minor.tick(nx=3, ny=3, tick.ratio=0.5)

legend("topleft", inset=0.05, title = "Drug Type", c("A", "B"),
       lty=c(1, 2), pch=c(15, 17), col=c("red", "blue"))

par(opar)

# Math Annotations
demo(plotmath)

# 3.4.5. Text annotations
text(location, "text to place", pos, ...)
mtext("text to place", side, line=n, ...)

# 3.5. Combining graphs
help("layout")

# 3.5.1. Creating a figure arrangement with fine control
# code listing 3.4
opar <- par(no.readonly = T)
par(fig=c(0, 0.8, 0, 0.8))
plot(mtcars$wt, mtcars$mpg,
     xlab = "Miles Per Gallon",
     ylab = "Car Weight")

par(fig=c(0, 0.8, 0.55, 1), new=TRUE)
boxplot(mtcars$wt, horizontal = T, axes=FALSE)

par(fig=c(0.65, 1, 0, 0.8), new=TRUE)
boxplot(mtcars$mpg, axes=F)

mtext("Enhanced Scatterplot", side = 3, outer = T, line = -3)
par(opar)

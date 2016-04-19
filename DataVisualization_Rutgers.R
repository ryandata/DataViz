# Data Visualization and R
# Ryan Womack, Data Librarian, rwomack@rutgers.edu
# R script to accompany 
# Rutgers Data Visualization Workshop
# Spring 2016 (2016-03-31 version)
# For all files, see
# http://libguides.rutgers.edu/data_R

#######################
# SECTION 1
# INSTALLING PACKAGES
#######################

install.packages("lattice", dependencies=TRUE)
install.packages("ggplot2", dependencies=TRUE)

# these come with ggplot2, but might need to be installed separately if 
#     dependencies fail
# install.packages("RColorBrewer", dependencies=TRUE)
# install.packages("maps", dependencies=TRUE)
# install.packages("latticeExtra", dependencies=TRUE)
# install.packages("hexbin", dependencies=TRUE)
# install.packages("plyr", dependencies=TRUE)
# install.packages("reshape2", dependencies=TRUE)
# install.packages("scales", dependencies=TRUE)
# install.packages("mapproj", dependencies=TRUE)

install.packages("party", dependencies=TRUE)
# install.packages("vcd", dependencies=TRUE)

install.packages("klaR", dependencies=TRUE)
# install.packages("som", dependencies=TRUE)

install.packages("kohonen", dependencies=TRUE)

install.packages("timeline", dependencies=TRUE)
# install.packages("shiny", dependencies=TRUE)

install.packages("vcdExtra", dependencies=TRUE)
# install.packages("rgl", dependencies=TRUE)
# install.packages("grid", dependencies=TRUE)

install.packages("animation", dependencies=TRUE)
install.packages("tabplot", dependencies=TRUE)
install.packages("googleVis", dependencies=TRUE)
install.packages("rainfreq", dependencies=TRUE)

#playwith requires GTK+2.0
#install via your OS package manager, or download separately
install.packages("playwith", dependencies=TRUE)

#rggobi requires ggobi
#install via your OS package manager, or download separately
#this appears to be a real pain on the Mac - not actively maintained
install.packages("rggobi", dependencies=TRUE)

# rattle has many dependencies to be fully functional
# don't worry too much for demo purposes
# you may need a version JDK (Java Development Kit) and ODBC installed separately too
install.packages("rattle", dependencies=TRUE)
# install.packages("latticist", dependencies=TRUE)

install.packages("aplpack", dependencies=TRUE)
install.packages("memoise", dependencies=TRUE)
install.packages("devtools", dependencies=TRUE)
install.packages("dplyr", dependencies=TRUE)


#devtools allows the installation of unofficial packages from locations like github
devtools::install_github("ramnathv/rCharts")
  # see http://ramnathv.github.io/rCharts/
  # for more info

#on Windows systems ---
    #for bigvis, you will need to install Rtools separately 
    #see cran.r-project.org/bin/windows/Rtools
    #requires restart of R too
devtools::install_github("hadley/bigvis")
devtools::install_github(c("hadley/testthat", "rstudio/ggvis"))



#######################
# SECTION 2
# lattice
#######################

# 2.1 load data

library(lattice)
library(ggplot2)
data(diamonds)
?diamonds
attach(diamonds)

# 2.2 scatterplot

xyplot(price~carat)
xyplot(x*y*z~carat)
xyplot(price~carat | clarity)

# 2.3 barchart

barchart(table(cut))
barchart(table(clarity), horizontal=FALSE)
barchart(table(clarity,cut),horizontal=FALSE, stack=FALSE)

# 2.4 groups

xyplot(price~carat, groups=cut)
xyplot(price~carat | cut + clarity)
xyplot(price~carat | cut , groups=clarity, auto.key=list(space="right"))

# 2.5  store and modify lattice object

diamondgraph<-xyplot(price~carat | cut)
update(diamondgraph, aspect="fill", layout=c(1,5))
update(diamondgraph, panel=panel.barchart)
update(diamondgraph, col="tomato") 

# 2.6 shingles

diamonds$shingle <- equal.count(price, number=5, overlap=0)
xyplot(carat~depth | diamonds$shingle)
xyplot(carat~depth | diamonds$shingle, strip=strip.custom(strip.levels=TRUE, strip.names=FALSE))

# 2.7 tweaking, regression

xyplot(price~carat, col="steelblue", pch=3, main="Diamond Data", xlab="weight of diamond in carats", ylab="price of diamond in dollars", xlim=c(0,3), scales=list(tick.number=10))
xyplot(price~carat, type=c("p","r"))

# this command will take a long time to run, especially in Rstudio
# you may want to save for later
splom(diamonds, groups=clarity)
parallelplot(diamonds, groups=clarity, horizontal.axis=FALSE)

# 2.8 modifying settings in lattice

show.settings()
colors()
mytheme<-trellis.par.get("plot.line")
mytheme$col="tomato"
trellis.par.set("plot.line",mytheme) 
print(xyplot(price~carat, type=c("p","r")))


#######################
# SECTION 3
# ggplot2
#######################

# 3.1 scatterplot

ggplot(diamonds, aes(carat, price)) + geom_point()
ggplot(diamonds, aes(carat, x*y*z)) + geom_point()
ggplot(diamonds, aes(carat, price)) + geom_point() + facet_grid(.~clarity)
ggplot(diamonds, aes(carat, price)) + geom_point(aes(color=cut))
ggplot(diamonds, aes(carat,price)) + xlim(0,3) + geom_point(colour="steelblue", pch=3) + labs (x="weight of diamond in carats", y="price of diamond in dollars", title="Diamond Data") 

# 3.2 histogram

ggplot(diamonds, aes(depth))+geom_histogram()
ggplot(diamonds, aes(depth))+geom_histogram(aes(fill = ..count..))

# 3.3 barchart

ggplot(diamonds, aes(cut) )+ geom_bar(position="stack") 
ggplot(diamonds, aes(clarity) )+ geom_bar(position="stack") 
ggplot(diamonds, aes(clarity)) +facet_grid(.~cut) + geom_bar(position="dodge")

# 3.4 theme tweaks

library(RColorBrewer)
ggplot(diamonds, aes(clarity)) +facet_grid(.~cut) + geom_bar(position="dodge", fill="purple")+theme(panel.background = element_rect(fill='pink', colour='green'))
ggplot(diamonds, aes(clarity)) +facet_grid(.~cut) + geom_bar(position="dodge")+theme(panel.background = element_rect(fill='white', colour='black'))
ggplot(diamonds, aes(x=clarity, fill=factor(clarity))) +facet_grid(.~cut) + geom_bar(position="dodge")+ scale_fill_brewer(palette="Reds")

# this option can be used in some contexts
# scale_color_manual(values = c("yellow","orange","pink","red","purple"))


# 3.5 regression

ggplot(diamonds, aes(carat, price)) + geom_point() + geom_smooth(method=lm)
ggplot(diamonds, aes(carat, price)) + geom_point() + stat_smooth()
ggplot(mtcars, aes(mpg, disp)) + geom_point() + stat_smooth()
# we can also plot customized confidence interval bands, but this requires computing them separately [see ggplot2 help]

# 3.6 using the power

mydata <- ggplot(diamonds, aes(clarity)) +facet_grid(.~cut) 
mytheme <- theme(panel.background = element_rect(fill='lightblue', colour='darkgrey'))
mychart <- geom_bar(position="dodge", fill="thistle", color="black")

mydata+mytheme+mychart

# 3.7 Minard in ggplot

troops <- read.table("http://ryanwomack.com/IASSIST/DataViz/Data/troops.txt", header=T)
cities <- read.table("http://ryanwomack.com/IASSIST/DataViz/Data/cities.txt", header=T)
temps <- read.table("http://ryanwomack.com/IASSIST/DataViz/Data/temps.txt", header=T)
temps$date <- as.Date(strptime(temps$date,"%d%b%Y"))

library(maps)
borders <- data.frame(map("world", xlim=c(10,50), ylim=c(40, 80), plot=F)[c("x","y")])

xlim <- scale_x_continuous(limits = c(24, 39))

ggplot(cities, aes(x = long, y = lat)) + 
  geom_path(
    aes(size = survivors, colour = direction, group = group), 
    data=troops
  ) + 
  geom_point() + 
  geom_text(aes(label = city), hjust=0, vjust=1, size=4) + 
  scale_size(range = c(1, 10)) + 
  scale_colour_manual(values = c("grey50","red")) +
  xlim


ggsave(file = "march.pdf", width=16, height=4)

qplot(long, temp, data=temps, geom="line") + 
  geom_text(aes(label = paste(day, month)), vjust=1) + xlim

ggsave(file = "temps.pdf", width=16, height=4)


##########################
# SECTION 4
# A Miscellany of Methods
##########################

# 4.1 exporting images

pdf(file="output.pdf")
xyplot(carat~depth | diamonds$shingle, strip=strip.custom(strip.levels=TRUE, strip.names=FALSE))
mydata+mytheme+mychart
dev.off()

jpeg(file="output.jpg", width = 800, height = 600, quality=100)
mydata+mytheme+mychart
dev.off()

# 4.2 dot plots in lattice

attach(mtcars)

  # create factors with value labels

gear.f<-factor(gear,levels=c(3,4,5),
               labels=c("3gears","4gears","5gears"))
cyl.f <-factor(cyl,levels=c(4,6,8),
               labels=c("4cyl","6cyl","8cyl"))

  # dotplot for each combination of two factors

dotplot(cyl.f~mpg|gear.f,
        main="Dotplot Plot by Number of Gears and Cylinders",
        xlab="Miles Per Gallon",
        layout=(c(1,3)))

  # Dotplot: Grouped Sorted and Colored
  # Sort by mpg, group and color by cylinder

dotplot(reorder(row.names(mtcars),mpg)~mpg,
        pch=21, gpch=21, cex=.9,groups= cyl,
        main="Gas Mileage for Car Models\ngrouped by cylinder",
        xlab="Miles Per Gallon", auto.key=list(space="right")) 

  # multiple variable dotplot

combined <- cbind(mpg,disp/10)
dotplot(reorder(row.names(mtcars), disp) ~ mpg + disp/10,
        pch=c(1,2), gpch=c(2,3), cex=.9,
        main="Gas Mileage for Car Models\ngrouped by cylinder",
        xlab="Miles Per Gallon and Displacement/10") 


# 4.3 kernel density plots in lattice

  # kernel density plot

densityplot(~mpg,
            main="Density Plot",
            xlab="Miles per Gallon")

  # kernel density plots by factor level

densityplot(~mpg|cyl.f,
            main="Density Plot by Number of Cylinders",
            xlab="Miles per Gallon")

  # kernel density plots by factor level (alternate layout)

densityplot(~mpg|cyl.f,
            main="Density Plot by Number of Cylinders",
            xlab="Miles per Gallon",
            layout=c(1,3))

# 4.4 scatterplot matrix in lattice

splom(mtcars[c(1,3,4,5,6)],
      main="MTCARS Data")

splom(mtcars[c(1,3,4,5,6)], groups=cyl.f,
      main="MTCARS Data")

# 4.5 box and whiskers plots in lattice

bwplot(cyl.f~mpg|gear.f,
       ylab="Cylinders", xlab="Miles per Gallon",
       main="Mileage by Cylinders and Gears",
       layout=(c(1,3)))

# 4.6 stem and leaf

stem(cyl)
stem(mpg)
stem(mpg, scale=3)
stem(carat)

# 4.7 violin plot and "dot plot" in ggplot2
      # note there is also a vioplot package
# also heatmap

ggplot(diamonds, aes(cut, price)) + geom_violin()
ggplot(diamonds, aes(clarity, price)) + geom_violin()
ggplot(diamonds, aes(color, price)) + geom_violin()

ggplot(mtcars, aes(x = factor(vs), fill = factor(cyl), y = mpg)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge")

# heatmap
x  <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
hv <- heatmap(x, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              xlab = "specification variables", ylab =  "Car Models",
              main = "heatmap(<Mtcars data>, ..., scale = \"column\")")
utils::str(hv) # the two re-ordering index vectors

## no column dendrogram (nor reordering) at all:
heatmap(x, Colv = NA, col = cm.colors(256), scale = "column",
        RowSideColors = rc, margins = c(5,10),
        xlab = "specification variables", ylab =  "Car Models",
        main = "heatmap(<Mtcars data>, ..., scale = \"column\")")

## "no nothing"
heatmap(x, Rowv = NA, Colv = NA, scale = "column", col = cm.colors(256),
        main = "heatmap(*, NA, NA) ~= image(t(x))")



# 4.8 tree plots (party package)

library(party)

data(airquality)
airq <- subset(airquality, !is.na(Ozone))
airct <- ctree(Ozone ~ ., data = airq,
               controls = ctree_control(maxsurrogate = 3))
plot(airct)

# 4.9 self-organizing maps (som and kohonen packages)

  # som package

library(klaR)
library(som)
data(countries)
logcount <- log(countries[,2:7])
sdlogcount <- apply(logcount, 2, sd)
logstand <- t((t(logcount) / sdlogcount) * c(1,2,6,5,5,3))
cclasses <- cutree(hclust(dist(logstand)), k = 6)
countryEDAM <- EDAM(logstand, classes = cclasses, sa = FALSE, 
                    iter.max = 10, random = FALSE)
plot(countryEDAM, vertices = FALSE, label = TRUE, stck = FALSE)

  # kohonen package

library(kohonen)
data(wines)
set.seed(7)

kohmap <- xyf(scale(wines), classvec2classmat(wine.classes),
              grid = somgrid(5, 5, "hexagonal"), rlen=100)
plot(kohmap, type="changes")
plot(kohmap, type="codes", main = c("Codes X", "Codes Y"))
plot(kohmap, type="counts")


# 4.10 clustering

  # reduce size of data to speed processing

mydata <- diamonds[sample(1:nrow(diamonds),50),]

  # Ward Hierarchical Clustering

d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
  
  # draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=5, border="red") 

# 4.11 mosaic plots
  #note there is a mosaicplot command in base r, but we use the
  #slightly more polished vcd package here

library(vcd)
data(Titanic)
ftable(Titanic, row.vars = 3:1)
mosaic(~ Survived, data = Titanic,
       main = "Survival on the Titanic", shade = TRUE, legend = TRUE)
mosaic(~ Sex + Survived, data = Titanic,
       main = "Survival on the Titanic", shade = TRUE, legend = TRUE)
mosaic(~ Sex + Age + Survived, data = Titanic,
       main = "Survival on the Titanic", shade = TRUE, legend = TRUE)
mosaic(~ Sex + Age + Class + Survived, data = Titanic,
       main = "Survival on the Titanic", shade = TRUE, legend = TRUE)
mosaic(~ Sex + Class + Survived, data = Titanic,
       main = "Survival on the Titanic", shade = TRUE, legend = TRUE)
mosaic(Survived ~ ., data = Titanic,
       main = "Survival on the Titanic")
spineplot(cut~clarity)

#see also other mosaic implementations
mosaicplot(Titanic)
strucplot(Titanic)

# 4.12 timeline

  #uses ggplot2 as a basis

library(timeline)
data(ww2)
timeline(ww2, ww2.events, event.spots=2, event.label="", event.above=FALSE)

# 4.13 choropleth maps

detach(package:kohonen)
library(lattice)
library(maps)
county.map <- map("county", plot = FALSE, fill = TRUE)
data(ancestry, package = "latticeExtra")
ancestry <- subset(ancestry, !duplicated(county))
rownames(ancestry) <- ancestry$county 
freq <- table(ancestry$top)
keep <- names(freq)[freq > 10] 
ancestry$mode <- with(ancestry, factor(ifelse(top %in% keep, top, "Other"))) 
modal.ancestry <- ancestry[county.map$names, "mode"] 
library("RColorBrewer") 
colors <- brewer.pal(n = nlevels(ancestry$mode), name = "Pastel1")
xyplot(y ~ x, county.map, aspect = "iso", scales = list(draw = FALSE), xlab = "", ylab = "", par.settings = list(axis.line = list(col = "transparent")), col = colors[modal.ancestry], border = NA, panel = panel.polygon, key = list(text = list(levels(modal.ancestry), adj = 1), rectangles = list(col = colors), x = 1, y = 0, corner = c(1, 0))) 

# 4.14 choropleth maps w/projection correction

rad <- function(x) { pi * x / 180 } 
county.map$xx <- with(county.map, cos(rad(x)) * cos(rad(y)))
county.map$yy <- with(county.map, sin(rad(x)) * cos(rad(y))) 
county.map$zz <- with(county.map, sin(rad(y))) 
panel.3dpoly <- function (x, y, z, rot.mat = diag(4), distance, ...)
{
  m <- ltransform3dto3d(rbind(x, y, z), rot.mat, distance)
  panel.polygon(x = m[1, ], y = m[2, ], ...)
}
aspect <- with(county.map, c(diff(range(yy, na.rm = TRUE)), diff(range(zz, na.rm = TRUE))) / diff(range(xx, na.rm = TRUE)))
cloud(zz ~ xx * yy, county.map, par.box = list(col = "grey"), aspect = aspect, panel.aspect = 0.6, lwd = 0.01, panel.3d.cloud = panel.3dpoly, col = colors[modal.ancestry], screen = list(z = 10, x = -30), key = list(text = list(levels(modal.ancestry), adj = 1), rectangles = list(col = colors), space = "top", columns = 4), scales = list(draw = FALSE), zoom = 1.1, xlab = "", ylab = "", zlab = "")


# 4.15 choropleth maps with mapproj

library("latticeExtra")
library("mapproj")
data(USCancerRates)
rng <- with(USCancerRates, range(rate.male, rate.female, finite = TRUE)) 
nbreaks <- 50 
breaks <- exp(do.breaks(log(rng), nbreaks)) 
mapplot(rownames(USCancerRates) ~ rate.male + rate.female, data = USCancerRates, breaks = breaks, map = map("county", plot = FALSE, fill = TRUE, projection = "tetra"), scales = list(draw = FALSE), xlab = "", main = "Average yearly deaths due to cancer per 100000") 

#4.16 rainfreq
library(rainfreq)
library(maps)
x_se <- extract_freq()
print(x_se)

x_mw <- extract_freq(region_name = "mw", storm_RP = 1000, storm_duration = "48h")
print(x_mw)

x_hi <- extract_freq(region_name = "hi", storm_RP = 10, storm_duration = "6h")
print(x_hi)

data(rain_max_usa)
head(rain_max_usa)

data(rain_max_world)
head(rain_max_world)

x_se <- x_se * 0.001
x_mw <- x_mw * 0.001
x_hi <- x_hi * 0.001

# Here is a plot of the three rainfall estimates obtained so far.
# State boundaries are added for spatial reference.
# southeast
plot(x_se, breaks = c(6, 9, 12, 15, 18),
     col = c("red", "yellow", "green", "blue"),
     main = "100-year 24-hour Rainfall for Southeast USA (inches)")
map('state', region = c('florida', 'arkansas', 'louisiana', 'mississippi',
                        'alabama', 'georgia'), add = TRUE)
# midwest
plot(x_mw, breaks = c(2, 5, 10, 15, 20),
     col = c("red", "yellow", "green", "blue"),
     main = "1000-year 48-hour Rainfall for Midwest USA (inches)")
map('state', region = c('colorado', 'north dakota', 'south dakota', 'nebraska',
                        'oklahoma', 'minnesota', 'iowa', 'missouri',
                        'wisconsin', 'michigan'), add = TRUE)
# hawaii
plot(x_hi, breaks = c(1, 3, 6, 9, 12),
     col = c("red", "yellow", "green", "blue"),
     main = "10-year 6-hour Rainfall for Hawaii (inches)")


# 4.17 Chernoff faces 

library(aplpack)
set.seed(17)
faces(matrix(sample(1:1000,128,),16,8),main="random faces")


##########################
# SECTION 5
# 3-D VIZUALIZATION
##########################

# 5.1 cloud 3D scatterplot in lattice 
  #note there is a scatterplot3d package too

cloud(price~carat*table)
cloud(price~carat*table | clarity)
cloud(price~carat*table | clarity, screen=list(z=-60,x=-60))
cloud(mpg~wt*qsec|cyl,
      main="3D Scatterplot by Cylinders")

    # cloud.table

cloud(prop.table(Titanic, margin = 1:3),
      type = c("p", "h"), strip = strip.custom(strip.names = TRUE),
      scales = list(arrows = FALSE, distance = 2), panel.aspect = 0.7,
      zlab = "Proportion")[, 1]

# 5.2 levelplot and contourplot in lattice 

  # setup some smooth data
x <- seq(pi/4, 5 * pi, length.out = 100)
y <- seq(pi/4, 5 * pi, length.out = 100)
r <- as.vector(sqrt(outer(x^2, y^2, "+")))
grid <- expand.grid(x=x, y=y)
grid$z <- cos(r^2) * exp(-r/(pi^3))

  # and plot

levelplot(z~x*y, grid, cuts = 100, scales=list(log="e"), xlab="",
          ylab="", main="Weird Function", sub="with log scales",
          colorkey = FALSE, region = TRUE)

  # contourplot setup data

library(stats)
attach(environmental)
ozo.m <- loess((ozone^(1/3)) ~ wind * temperature * radiation,
               parametric = c("radiation", "wind"), span = 1, degree = 2)
w.marginal <- seq(min(wind), max(wind), length.out = 50)
t.marginal <- seq(min(temperature), max(temperature), length.out = 50)
r.marginal <- seq(min(radiation), max(radiation), length.out = 4)
wtr.marginal <- list(wind = w.marginal, temperature = t.marginal,
                     radiation = r.marginal)
grid <- expand.grid(wtr.marginal)
grid[, "fit"] <- c(predict(ozo.m, grid))

  # contourplot command

contourplot(fit ~ wind * temperature | radiation, data = grid,
            cuts = 10, region = TRUE,
            xlab = "Wind Speed (mph)",
            ylab = "Temperature (F)",
            main = "Cube Root Ozone (cube root ppb)")

  #ggplot2 uses geom_tile() and geom_contour() for these functions
  #ggplot2 does not perform true 3-D visualization

# 5.3 wireframe in lattice

wireframe(volcano, shade = TRUE,
          aspect = c(61/87, 0.4),
          light.source = c(10,0,10))

g <- expand.grid(x = 1:10, y = 5:15, gr = 1:2)
g$z <- log((g$x^g$gr + g$y^2) * g$gr)
wireframe(z ~ x * y, data = g, groups = gr,
          scales = list(arrows = FALSE),
          drape = TRUE, colorkey = TRUE,
          screen = list(z = 30, x = -60))

# 5.4 rgl

library(rgl)

example(surface3d)
example(plot3d)

plot3d(price~carat*color, col=rainbow(1000))

  # mosaic3d

library(vcdExtra)
mosaic3d(structable(Sex + Class ~ Survived + Age, data = Titanic))

  # persp3d and play3d

lat <- matrix(seq(90,-90, len=50)*pi/180, 50, 50, byrow=TRUE)
long <- matrix(seq(-180, 180, len=50)*pi/180, 50, 50)

r <- 6378.1 # radius of Earth in km
x <- r*cos(lat)*cos(long)
y <- r*cos(lat)*sin(long)
z <- r*sin(lat)

open3d()
persp3d(x, y, z, col="white", 
        texture=system.file("textures/worldsmall.png",package="rgl"), 
        specular="black", axes=FALSE, box=FALSE, xlab="", ylab="", zlab="",
        normal_x=x, normal_y=y, normal_z=z)
play3d(spin3d(axis=c(0,0,1), rpm=2), duration=60)

  # another persp3d example

set.seed(1)
x <- seq(1,10)
y <- seq(1,10)
w <- runif(100)
z <- runif(100)

wcolors <- rainbow(length(w))[rank(w)]
zmat <- matrix(z, length(x),length(y))

persp3d(x=x, y=y, z=zmat, col = wcolors)


######################
# SECTION 6
# ANIMATION
######################

library(animation)

# basic example
# gif requires ImageMagick
saveGIF(
  {brownian.motion()
  }
, movie.name="mygraphs.gif", img.name="brownian")

# more complex example
# illustrating a probability distribution

#create a range to plot with
values<-seq(0,10,0.01)
j<-values

#we can look at individual ranges for pdf, holding one variable constant and letting the other fluctuate
# the label on the x axis is the value of the median
for (i in 1:10) plot (x=values, y=dweibull(values,i,1), xlab=qweibull(0.5,i,j), ylab="proportion", main=c("Weibull pdf, scale and shape",i,1),pch=20)
for (i in 1:10) plot (x=values, y=dweibull(values,1,i), xlab=qweibull(0.5,i,j), ylab="proportion", main=c("Weibull pdf, scale and shape",1,i),pch=20)

#and do the same thing for CDF
for (i in 1:10) plot (x=values, y=pweibull(values,i,1), xlab=qweibull(0.5,i,j), ylab="proportion", main=c("Weibull pdf, scale and shape",i,1),pch=20)
for (i in 1:10) plot (x=values, y=pweibull(values,1,i), xlab=qweibull(0.5,i,j), ylab="proportion", main=c("Weibull pdf, scale and shape",1,i),pch=20)

#if we want to just output the median values of weibull
for (i in 1:10) 
  for (j in 1:10)
    print(qweibull(0.5,i,j))

#now plot all 100 combinations, shape from 1 to 10, scale from 1 to 10 
#use "animation" package to generate HTML version


saveHTML({
  for (i in 1:10) {
    for (j in 1:10) {
      plot(x=values, y=dweibull(values,i,j),xlab=qweibull(0.5,i,j), ylab="proportion", main=c("Weibull pdf, scale and shape",i,j),pch=20)
      abline(v=qweibull(0.5,i,j))
    }}}
  , img.name="Weibullpdf.html")

saveHTML({
  for (i in 1:10) {
    for (j in 1:10) {    plot(x=values, y=pweibull(values,i,j),xlab=qweibull(0.5,i,j), ylab="proportion", main=c("Weibull CDF, scale and shape",i,j),pch=20)
                         abline(v=qweibull(0.5,i,j))
    }}}
  , img.name="WeibullCDF.html")



######################
# SECTION 7
# INTERACTIVE DATA
######################

  # linked data panels are used in several contexts

# 7.1 time series

data(AirPassengers)
AP<-AirPassengers
decompose(AP)
plot(decompose(AP))

# 7.2 tableplot

library(tabplot)
diamonds2<-diamonds[,1:7]
tableplot(diamonds2)


# 7.3 googleVis

library(googleVis)

# World Population
# May not work due to Flash security restrictions

WorldPopulation=data.frame(Country=Population$Country, 
                           Population.in.millions=round(Population$Population/1e6,0),
                           Rank=paste(Population$Country, "Rank:", Population$Rank))

G5 <- gvisGeoMap(WorldPopulation, "Country", "Population.in.millions", "Rank", 
                 options=list(dataMode="regions", width=1800, height=900))
plot(G5)

# motion charts

M <- gvisMotionChart(Fruits, idvar="Fruit", timevar="Year", options=list(width=2000, height=1000)) 
plot(M)

#be Hans Rosling

demo(WorldBank)


# 7.4 playwith

library(playwith)

  # you can try playwith with the diamonds data
  # but it is often sluggish with a dataset this size
  # playwith(cloud(price~carat*table))
  
  # or use this smaller dataset for less memory demands and faster response
playwith(cloud(mpg~wt*qsec))

  # works for any function
playwith(xyplot(mpg~wt))

# 7.5 rggobi

  # this package is no longer updated and getting harder to install all requirements
library(rggobi)
ggobi(mtcars)

# 7.6 rattle 

library(rattle)
rattle()

# 7.7 latticist

library(latticist)
latticist(mtcars)

#7.8  rCharts

  #see http://ramnathv.github.io/rCharts/ for more info

library(rCharts)

# first, mouseovers, rPlot
## Example 1 Facetted Scatterplot
names(iris) = gsub("\\.", "", names(iris))
rPlot(SepalLength ~ SepalWidth | Species, data = iris, color = 'Species', type = 'point')

## Example 2 Facetted Barplot
hair_eye = as.data.frame(HairEyeColor)
rPlot(Freq ~ Hair | Eye, color = 'Eye', data = hair_eye, type = 'bar')

#slightly more interactivity, mPlot
data(economics, package = "ggplot2")
econ <- transform(economics, date = as.character(date))
mPlot(x = "date", y = c("psavert", "uempmed"), type = "Line", data = econ)

#see the web site for more examples

#7.9 shiny
library(shiny)
system.file("examples", package="shiny")

#try any of the following built in examples
runExample("01_hello") # a histogram
runExample("02_text") # tables and data frames
runExample("03_reactivity") # a reactive expression
runExample("04_mpg") # global variables
runExample("05_sliders") # slider bars
runExample("06_tabsets") # tabbed panels
runExample("07_widgets") # help text and submit buttons
runExample("08_html") # shiny app built from HTML
runExample("09_upload") # file upload wizard
runExample("10_download") # file download wizard
runExample("11_timer") # an automated timer

#7.10 ggvis
library(dplyr)
library(ggvis)

mtcars %>%
  ggvis(~wt, ~mpg) %>%
  layer_smooths(span = input_slider(0.5, 1, value = 1)) %>%
  layer_points(size := input_slider(100, 1000, value = 100))

mtcars %>% ggvis(x = ~wt) %>%
  layer_densities(
    adjust = input_slider(.1, 2, value = 1, step = .1, label = "Bandwidth adjustment"),
    kernel = input_select(
      c("Gaussian" = "gaussian",
        "Epanechnikov" = "epanechnikov",
        "Rectangular" = "rectangular",
        "Triangular" = "triangular",
        "Biweight" = "biweight",
        "Cosine" = "cosine",
        "Optcosine" = "optcosine"),
      label = "Kernel")
  )


##############################
# SECTION 8
# BIG DATA
##############################

airline<-read.csv("http://rci.rutgers.edu/~rwomack/IASSIST/DataViz/airlineJan.csv")
head(airline)
nrow(airline)
dim(airline)

library(bigvis)
library(ggplot2)
library(plyr)
library(grid)
library(reshape2)
library(scales)
library(memoise)
library(hexbin)

delay <- airline$ArrDelay
dist <- airline$Distance
time <- airline$AirTime

speed <- dist / time * 60

xyplot(time~dist)
hexbinplot(time~dist)

# Determine number of missing values ---------------------

bad_speed <- speed < 0 | speed > 761
bad_time <- time < 0
bad_dist <- dist > 2724

tabulate(bad_speed + bad_time + bad_time + 1, 4)

# Compare distributions of missings to non-missings ---------------

dist[bad_dist] <- NA
speed[bad_speed] <- NA

sd <- condense(bin(dist, 10), bin(speed, 10))

sd2 <- transform(sd, dist = is.na(dist))
sd2 <- ddply(sd2, "dist", mutate, .count = .count / sum(.count))

qplot(speed, .count, data = sd2, geom = "line", colour = factor(dist)) +
  scale_colour_hue("Distance", breaks = c(0, 1),
                   labels = c("Present", "Missing")) +
  theme(
    legend.position = c(0, 1),
    legend.justification = c("left", "top"),
    plot.margin = unit(c(0, 0.5, 0, 0.5), "lines")) +
  ylab("Proportion")
ggsave("speed-distance.pdf", width = 6, height = 3)
ggsave("speed-distance.jpg", width = 6, height = 3)

#Condense

tweak <- list(
  scale_x_continuous("Distance (miles)", breaks = c(500, 1000, 1500, 2000, 2500)),
  ylab(NULL),
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "lines"),
    text = element_text(size = 18),
    panel.margin = unit(0.25, "cm")
  )
)

dist_sum <- condense(bin(dist, 10), z = speed, summary = "sd")
dist_sum$.count <- dist_sum$.count / 1e6

dsm <- data.frame(
  dist = rep(dist_sum$dist, 3),
  summary = rep(c("count (millions)", "mean", "sd"), each = nrow(dist_sum)),
  value = c(dist_sum$.count, dist_sum$.mean, dist_sum$.sd)
)

qplot(dist, value, data = dsm, geom = "line") +
  facet_grid(summary ~ ., scales = "free_y") +
  tweak
ggsave("condense.pdf", width = 8, height = 7)
ggsave("condense.jpg", width = 8, height = 7)


dist_sum <- condense(bin(dist, 10), z = speed, summary = "sd")

# Figure out optimal smoothing characteristics --------------------

best_h(dist_sum, var = ".count")
best_h(dist_sum, var = ".mean")
best_h(dist_sum, var = ".sd")

dist_sum2 <- transform(dist_sum, 
   .count = scale(.count), 
   .mean = scale(.count), 
   .sd = scale(.sd))
 
 best <- list(
   count = best_h(dist_sum2, var = ".count"),
   mean = best_h(dist_sum2, var = ".mean"),
   sd = best_h(dist_sum2, var = ".sd"))
 dist_summary<-unlist(best)
 bestdf <- data.frame(
   summary = names(best), 
   dist_summary
 )


grid <- h_grid(dist_sum2, max = 10, n = 50)
rmses <- list(
  count = rmse_cvs(dist_sum2, grid, var = ".count"),
  mean = rmse_cvs(dist_sum2, grid, var = ".mean"),
  sd = rmse_cvs(dist_sum2, grid, var = ".sd")  
)
rmses <- Map(function(x, n) {
  x$summary <- n
  x
}, rmses, names(rmses))
rmse <- do.call(rbind, rmses)

qplot(dist, pmin(err, 2), data = rmse, colour  = summary, geom = "line") +
  ylab("rmse") +
  scale_x_continuous("Bandwidth (miles)") + 
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "lines")
  )
ggsave("smooth-rmse.pdf", width = 6, height = 3)
ggsave("smooth-rmse.jpg", width = 6, height = 3)

# Display smoothed data in one plot -----------------------------------

smoothes <- list(
  count = smooth(dist_sum, 50, var = ".count", type = "robust"),
  mean = smooth(dist_sum, 50, var = ".mean", type = "robust"),
  sd = smooth(dist_sum, 50, var = ".sd", type = "robust"))
smoothes <- Map(function(x, n) {
  names(x)[2] <- "value"
  x$summary <- n
  x
}, smoothes, names(smoothes))
smooth <- do.call(rbind, smoothes)
smooth$summary <- factor(smooth$summary, levels = c("count", "mean", "sd"))
levels(smooth$summary)[1] <- "count (millions)"

qplot(dist, value, data = smooth, geom = "line") +
  facet_grid(summary ~ ., scales = "free_y") +
  tweak
ggsave("smooth.pdf", width = 8, height = 7)
ggsave("smooth.jpg", width = 8, height = 7)

# 2d ----------------

ds <- condense(bin(dist, 50), bin(speed, 20))

ggplot(ds, aes(dist, speed, fill = .count)) +
  geom_raster() +
  scale_fill_gradient("Count\n(x 1000)", low = "grey90", high = "black",
                      breaks = c(1, 2, 3, 5, 10) * 1e5,
                      labels = c("100", "200", "300", "500", "1000"),
                      trans = mt_trans(0.5), guide = guide_colorbar(
                        title.vjust = 0.75, barwidth = unit(5, "inches")
                      )
  ) +
  scale_x_continuous("Distance (miles)",
                     breaks = c(500, 1000, 1500, 2000, 2500)) +
  ylab("Speed (mph)") +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "lines"),
    text = element_text(size = 18),
    legend.position = "bottom"
  )
ggsave("condense-2d.pdf", width = 8, height = 6)
ggsave("condense-2d.jpg", width = 8, height = 6)

# Explore 1d distributions -----------------------------------------------------
speed_sum <- condense(bin(speed, 5))
autoplot(speed_sum)
autoplot(smooth(speed_sum, 26))

dist_sum <- condense(bin(dist, 10))
autoplot(dist_sum)
autoplot(smooth(dist_sum, 55))

ds <- condense(bin(dist, 50), bin(speed, 25))
autoplot(ds)
autoplot(peel(ds))

teaser <- list(
  theme(
    legend.position = "bottom",
    plot.margin = unit(c(0, 0.5, 0, 0.5), "lines"),
    legend.key.width = unit(1.45, "inches"),
    text = element_text(size = 24)
  ),
  labs(x = NULL, y = NULL, fill = NULL))

# Make teaser images -----------------------------------------------------------
dsd <- condense(bin(dist, 20), bin(speed, 20), z = delay)
dsd <- subset(dsd, dist < 2700)
autoplot(dsd) + teaser
ggsave("teaser-1.pdf", width = 8, height = 6)
ggsave("teaser-1.jpg", width = 8, height = 6)

dsd2 <- peel(dsd, .995)
autoplot(dsd2) + teaser
ggsave("teaser-2.pdf", width = 8, height = 6)
ggsave("teaser-2.jpg", width = 8, height = 6)

autoplot(dsd2) + scale_fill_gradient2(trans = mt_trans(0.25),
                                      breaks = c(-50, -10, 0, 10, 40, 100, 200, 400)) + teaser
ggsave("teaser-3.pdf", width = 8, height = 6)
ggsave("teaser-3.jpg", width = 8, height = 6)


# Other exploration ------------------------------------------------------------

dsdp <- peel(dsd)
ggplot(dsdp, aes(dist, speed)) +
  geom_point(aes(colour = .mean, size = .count)) +
  scale_size_area() +
  scale_colour_gradient2()

ggplot(dsdp, aes(dist, speed)) +
  geom_raster(aes(fill = mt(.mean, 0.25))) +
  scale_fill_gradient2()


ggplot(dsdp, aes(dist, speed)) +
  geom_raster(aes(fill = .mean)) +
  geom_contour(aes(z = .mean), colour = "grey50") +
  scale_fill_gradient2()

dsd_speed <- transform(dsd, speed = is.na(speed))
qplot(dist, .mean, data = subset(dsd_speed, .count > 100), geom = "line", colour = factor(speed))


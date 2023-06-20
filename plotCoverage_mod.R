rm(list=ls())

source("/home/lisa/Scripts/skewPlot_mod.R")

require("lattice")
require("signal")
require("RColorBrewer")

Args <- commandArgs();
argument <- Args[4];

# data acquisition
coverage <- read.table(argument, h=F)
names(coverage) <- c("pos", "cov")

# data packaging
analyses <- list(filters = list("Savitsky-Golay" = sg_filt),
                  # "Moving average" = av_filt),
                  params = list(c(51, 301, 501)),# c(51, 101, 201)),
                 colors = c("Blues")
                 )
data <- list(coverage = list(
               x = coverage$pos,
               y = coverage$cov),
               main = paste("Coverage:", argument) 
             )

# Control zone
#pdf <- FALSE
pdf <- TRUE
#file <- "test.jpeg"
file <- paste(argument,".pdf", sep = "")
xlim <- c(0, max(data[[1]]$x))
#xlim <- c(0, 500000)

# Plot. Uncomment pdf to print to pdf
#if (pdf) jpeg(file, w=850, h=600)#, paper="a4r")
if (pdf) pdf(file,paper="a4r")
plotCoverage(data, analyses, xlim=xlim)
if (pdf) dev.off()



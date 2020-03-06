# especially for plotting #

library(ggplot2)
library(RColorBrewer)
library("gplots")

set1 = brewer.pal(9,"Set1")
set2 = brewer.pal(9,"Set2")

theme_grid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            # panel.grid.major = element_blank(),
                             # panel.grid.minor = element_blank(),
                             panel.background = element_blank())

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.background = element_blank())


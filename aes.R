# especially for plotting #
packages = c(
  "ggplot2",
  "RColorBrewer",
  "Polychrome"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

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


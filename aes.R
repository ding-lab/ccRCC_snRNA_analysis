# especially for plotting #

# source ------------------------------------------------------------------
packages = c(
  "ggplot2",
  "RColorBrewer",
  "Polychrome"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

theme_grid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            # panel.grid.major = element_blank(),
                             # panel.grid.minor = element_blank(),
                             panel.background = element_blank())

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.background = element_blank())


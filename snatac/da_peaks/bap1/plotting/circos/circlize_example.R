library(circlize)

load(paste(system.file(package = "circlize"), "/extdata/DMR.RData", sep=""))

# rainfall
circos.initializeWithIdeogram(plotType = c("axis", "labels"))

bed_list = list(DMR_hyper, DMR_hypo)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))

circos.genomicDensity(bed_list[[1]], col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(bed_list[[2]], col = c("#0000FF80"), track.height = 0.1)

circos.clear()

bed = generateRandomBed(500)
bed = generateRandomBed(nr = 300, nc = 4)

## add y -axis
op = par(no.readonly = TRUE)

sectors = letters[1:8]
circos.par(points.overflow.warning = FALSE)
circos.par(gap.degree = 8)
circos.initialize(sectors, xlim = c(0, 10))
circos.trackPlotRegion(sectors, ylim = c(0, 10), track.height = 0.5)
par(cex = 0.8)
for(a in letters[2:4]) {
  circos.yaxis(side = "left", sector.index = a)
}
for(a in letters[5:7]) {
  circos.yaxis(side = "right", sector.index = a)
}
circos.clear()

par(op)


rand_col = function(k) {
  return(rgb(runif(k), runif(k), runif(k)))
}

library(circlize)

par(lwd = 0.5)
circos.par("cell.padding" = c(0, 0, 0, 0))
circos.initializeWithIdeogram(chromosome.index = "chr1", plotType = NULL)

posTransform.fun = function(region) {
  return(region)
}
bed = generateRandomBed(nr = 400, fun = function(k) rep("gene", k))
circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), panel.fun = function(region, value, ...) {
  circos.genomicText(region, value, y = 0, labels.column = 1, facing = "clockwise", adj = c(0, 0.5), cex = 0.5, posTransform = posTransform.fun)
}, track.height = 0.05, bg.border = NA)
circos.genomicPosTransformLines(bed, posTransform = posTransform.fun, track.height = 0.04, col = "orange")
cytoband = read.cytoband()$df
circos.genomicTrackPlotRegion(cytoband, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = cytoband.col(value[, 2]), border = NA, ...)
  cell.xlim = get.cell.meta.data("cell.xlim")
  cell.ylim = get.cell.meta.data("cell.ylim")
  circos.rect(cell.xlim[1], cell.ylim[1], cell.xlim[2], cell.ylim[2], border = "black")
  major.at = seq(0, cell.xlim[2], by = 5000000)
  major.labels = major.at/1000000
  l = major.at %% 50000000 == 0
  major.labels[l] = ""
  circos.axis("top", major.at = major.at, labels = major.labels, labels.facing = "clockwise", labels.cex = 0.4)
  circos.text(major.at[l], rep(1.7, sum(l)), paste0(major.at[l]/1000000, "MB"), cex = 0.8, facing = "clockwise", adj = c(0, 0.5), niceFacing = TRUE)
}, bg.border = NA, track.height = 0.02)

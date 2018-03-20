library(tidyverse)
library(rgl)
setwd("~/Dropbox/PhysiCell/output/")

dat <- read.delim("cells_119.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dat$genotype <- as.factor(dat$genotype)

allcolors = grDevices::colors()[grep('(gr(a|e)y)|(white)', grDevices::colors(), invert = T)]

# overall picture of tumor
#spheres3d(dat$x, dat$y, dat$z, radius = dat$radius, col = as.numeric(dat$genotype))

# movie of cross sections
dat$layers <- cut(dat$z, 100)

par3d(windowRect = c(0, 0, 800, 800))
xlim <- range(dat$x)
ylim <- range(dat$y)
zlim <- range(dat$z)
for(i in 1:length(levels(dat$layers)))
{
  rgl.clear()
  tempdat <- dat %>% filter(layers == levels(dat$layers)[i], genotype != "0")
  rgl.lines(xlim, c(0, 0), c(0, 0), color = "black")
  rgl.lines(c(0, 0), ylim, c(0, 0), color = "black")
  rgl.lines(c(0, 0), c(0, 0), zlim, color = "black")
  spheres3d(tempdat$x, tempdat$y, tempdat$z, radius = tempdat$radius, col = allcolors[as.numeric(tempdat$genotype) %% length(allcolors)], alpha = ifelse(as.numeric(tempdat$genotype) == 1, 0.2, 1))
  rgl.snapshot(paste("img_", sprintf("%04d", i), ".png", sep = ""))
}

location <- "~/Dropbox/PhysiCell/output"
system(paste("mogrify -format jpg ", location, "/*.png", sep = ""))
system(paste("ffmpeg -r 10 -y -i ", location, "/img_%04d.jpg -f mp4 -vcodec libx264 -pix_fmt yuv420p ", location, "/zaxis_noslice.mp4", sep = ""))
system("rm img*")


boundarydat <- dat %>% mutate(distance = sqrt(x^2 + y^2 + z^2)) %>% filter(distance > 150)
# movie of subclones
par3d(windowRect = c(0, 0, 800, 800))
xlim <- range(dat$x)
ylim <- range(dat$y)
zlim <- range(dat$z)
for(i in 2:length(table(dat$genotype)))
{
  rgl.clear()
  tempdat <- dat %>% filter(genotype == names(sort(table(dat$genotype), decreasing = TRUE))[i])
  if(nrow(tempdat) <= 10) break
  rgl.lines(xlim, c(0, 0), c(0, 0), color = "black")
  rgl.lines(c(0, 0), ylim, c(0, 0), color = "black")
  rgl.lines(c(0, 0), c(0, 0), zlim, color = "black")
  #spheres3d(boundarydat$x, boundarydat$y, boundarydat$z, radius = 1)
  spheres3d(tempdat$x, tempdat$y, tempdat$z, radius = tempdat$radius, col = allcolors[as.numeric(tempdat$genotype) %% length(allcolors)])
  rgl.snapshot(paste("img_", sprintf("%04d", i), ".png", sep = ""))
}

location <- "~/Dropbox/PhysiCell/output"
system(paste("mogrify -format jpg ", location, "/*.png", sep = ""))
system(paste("ffmpeg -r 10 -y -i ", location, "/img_%04d.jpg -f mp4 -vcodec libx264 -pix_fmt yuv420p ", location, "/subclone.mp4", sep = ""))
system("rm img*")


# split cells up into alleles
no_ancestor <- dat %>% filter(genotype != "0")
alleles <- strsplit(as.character(no_ancestor$genotype), ">")
alleles <- lapply(alleles, function(x) as.numeric(x)[-1])

# to make adjacency matrix
# d1 <- stack(setNames(alleles, seq_along(alleles)))
# Un1 <- unlist(alleles)
# m1 <- matrix(0, nrow=length(alleles), ncol=max(Un1))
# m1[cbind(as.numeric(d1$ind), d1$values)] <- 1

# reduce to top 20 alleles
d1 <- stack(setNames(alleles, no_ancestor$ID))
d1$ind <- as.numeric(as.character(d1$ind))
d2 <- d1[d1$values %in% as.numeric(names(table(d1$values))[table(d1$values) >= 30]),]
names(d2) <- c("gene", "ID")

variants <- no_ancestor %>% right_join(d2, by = "ID")
variants$gene <- as.factor(variants$gene)

par3d(windowRect = c(0, 0, 800, 800))
for(i in 1:length(unique(variants$gene)))
{
  rgl.clear()
  tempdat <- variants %>% filter(gene == unique(variants$gene)[i])
  rgl.lines(xlim, c(0, 0), c(0, 0), color = "black")
  rgl.lines(c(0, 0), ylim, c(0, 0), color = "black")
  rgl.lines(c(0, 0), c(0, 0), zlim, color = "black")
  #spheres3d(boundarydat$x, boundarydat$y, boundarydat$z, radius = 1)
  spheres3d(tempdat$x, tempdat$y, tempdat$z, radius = 10, col = as.numeric(tempdat$gene)+1)
  rgl.snapshot(paste("img_", sprintf("%04d", i), ".png", sep = ""))
}

location <- "~/Dropbox/PhysiCell/output"
system(paste("mogrify -format jpg ", location, "/*.png", sep = ""))
system(paste("ffmpeg -r 10 -y -i ", location, "/img_%04d.jpg -f mp4 -vcodec libx264 -pix_fmt yuv420p ", location, "/alleles.mp4", sep = ""))
system("rm img*")


head(variants)


# movie of cross sections
dat<- dat %>% mutate(distance = sqrt(x^2 + y^2 + z^2))
dat$layers <- cut(dat$distance, 50)

par3d(windowRect = c(0, 0, 800, 800))
xlim <- range(dat$x)
ylim <- range(dat$y)
zlim <- range(dat$z)
for(i in 1:length(levels(dat$layers)))
{
  rgl.clear()
  tempdat <- dat %>% filter(layers == levels(dat$layers)[i])
  rgl.lines(xlim, c(0, 0), c(0, 0), color = "black")
  rgl.lines(c(0, 0), ylim, c(0, 0), color = "black")
  rgl.lines(c(0, 0), c(0, 0), zlim, color = "black")
  spheres3d(tempdat$x, tempdat$y, tempdat$z, radius = tempdat$radius, col = as.numeric(tempdat$genotype), alpha = ifelse(as.numeric(tempdat$genotype) == 1, 0.2, 1))
  rgl.snapshot(paste("img_", sprintf("%04d", i), ".png", sep = ""))
}

location <- "~/Dropbox/PhysiCell/output"
system(paste("mogrify -format jpg ", location, "/*.png", sep = ""))
system(paste("ffmpeg -r 10 -y -i ", location, "/img_%04d.jpg -f mp4 -vcodec libx264 -pix_fmt yuv420p ", location, "/center_to_out.mp4", sep = ""))
system("rm img*")

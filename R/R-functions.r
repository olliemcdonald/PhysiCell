setwd("~/Dropbox/PhysiCell/output/")
library(tidyverse)
trajectory <- read.delim("../simulation_report.txt")
plot(trajectory$simulated.time, log(trajectory$num.cells))


cellfiles <- list.files()
cellfiles <- cellfiles[grepl(".txt", cellfiles)]

cells <- data.frame(curr_time = numeric(),
                    X = integer(),
                    ID = integer(),
                    x = numeric(),
                    y = numeric(), 
                    z = numeric(),
                    radius = numeric(),
                    phenotype = integer())

for(i in 1:length(cellfiles))
{
  curr_time <- as.numeric(strsplit(cellfiles[i], "_|[.]")[[1]][2])
  temp_cells <- read.delim(cellfiles[i])
  temp_cells <- temp_cells %>%
    cbind(curr_time) %>%
    select(curr_time, X, ID, x, y, z, radius, phenotype)
  cells <- rbind(cells, temp_cells)
}

system("rm ./plots/*.png")
j <- 1
for(i in sort(unique(cells$curr_time)))
{
  g <- cells %>% filter(curr_time == i) %>%
    ggplot(aes(x = x, y = y, size = z - min(cells$z))) +
    geom_point(alpha = 0.2) + scale_radius() + xlim(-200, 200) + ylim(-200,200) +
    theme_minimal() + theme(legend.position = "none")
  ggsave(g, filename = paste("./plots/img", sprintf("%04d", j), ".png", sep = ""),
         width = 3, height = 3)
  graphics.off()
  j <- j + 1
}

system("cd plots && ffmpeg -r 20 -f image2 -i img%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test.mp4 -y")

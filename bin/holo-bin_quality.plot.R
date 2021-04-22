library("argparse")
library("ggplot2")
library("tidyverse")

# Parse inputs
parser <-  ArgumentParser(description='Runs Holoflow.')
parser$add_argument('-cov_data', dest='cov_data', help='coverage data', required=TRUE)
parser$add_argument('-qual_data', dest='qual', help='quality data', required=TRUE)
parser$add_argument('-ID', dest='ID', help='ID', required=TRUE)
parser$add_argument('-out_path', dest='out_path', help='directory to redirect output', required=TRUE)
args <- parser$parse_args()

# Define variables
cov_file <- args$cov_data
qual <- args$qual
ID <- args$ID
out_path <- args$out_path


# Run
cov_data <- read.table(file=cov_file,header = T,stringsAsFactors = F) # fields 1,3
cov_data <- cov_data[,c(1,3)]
colnames(cov_data) <- c('MAGName','totalAvgDepth')

qual_data <- read.delim(file = qual,header = T, stringsAsFactors = F)
qual_data <- qual_data[,c(2,13,14)]
colnames(qual_data) <- c("ID","Completeness","Contamination")

# Generate df to plot: MAGid, completeness, contamination, avg coverage
# Ensure total avg depth correspond to given contamination/completeness
qual_data$avg_depth <- cov_data$totalAvgDepth[match(qual_data$ID,cov_data$MAGName)]



qual <- ggplot() + geom_point(data=qual_data, aes(x=Completeness, y=Contamination, size=avg_depth, col=avg_depth), alpha=0.5) +
labs(colour= "Total Average Depth", size="Total Average Depth")



dpi <- 96
ggsave(plot = qual,filename = paste0(out_path,'/',ID,'_quality.coverage_Plot.pdf'), width = 1800 / dpi, height = 900 / dpi,dpi = dpi)

library("argparse")
library("ggplot2")
library("tidyverse")

# Parse inputs
parser <-  ArgumentParser(description='Runs Holoflow.')
parser$add_argument('-cov_data', dest='cov', help='coverage data', required=TRUE)
parser$add_argument('-qual_data', dest='qual', help='quality data', required=TRUE)
parser$add_argument('-ID', dest='ID', help='ID', required=TRUE)
parser$add_argument('-out_path', dest='out_path', help='directory to redirect output', required=TRUE)
args <- parser$parse_args()

# Define variables
cov <- args$cov_data
qual <- args$qual_data
ID <- args$ID
out_path <- args$out_path



# Run
cov_data <- read.table(file=cov,header = T,quote = F,stringsAsFactors = F) # fields 1,3
qual_data <- read.delim(file = qual,header = T, stringsAsFactors = F)
qual_data <- as.data.frame(cbind(qual_data$Bin.Id,qual_data$Completeness,qual_data$Contamination))
colnames(qual_data) <- c("ID","Completeness","Contamination")

# Generate df to plot: MAGid, completeness, contamination, avg coverage
# Ensure total avg depth correspond to given contamination/completeness
qual_data$avg_depth <- cov_data$totalAvgDepth[match(qual_data$ID,cov_data$MAGName)]


qual <- ggplot()+geom_point(data=qual_data, aes(x=Completeness, y=Contamination, colour=avg_depth), size = 2)+
scale_colour_gradient(low="#566643", high="#eb1c1c", "Total Average Depth")

ggsave(plot = qual,filename = paste0(out_path,'/',ID,'_quality.coverage_Plot.pdf'))

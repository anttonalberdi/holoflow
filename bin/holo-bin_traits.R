library("argparse")
library("tidyverse")
################################### NOT IN USE NOW ##################################

# Parse inputs
parser <-  ArgumentParser(description='Runs Holoflow.')
parser$add_argument('-ar_summ', dest='gtdbtk_ar', help='archaeal gtdbtk', required=TRUE)
parser$add_argument('-bac_summ', dest='gtdbtk_bac', help='bacterial gtdbtk', required=TRUE)
parser$add_argument('-ID', dest='ID', help='ID', required=TRUE)
parser$add_argument('-out_file', dest='out_file', help='file to redirect output', required=TRUE)
args <- parser$parse_args()

# Define variables
gtdbtk_ar <- args$gtdbtk_ar
gtdbtk_bac <- args$gtdbtk_bac
ID <- args$ID
out_path <- args$out_path

# Run

# Read data
traits <- read.csv("/home/projects/ku-cbd/data/HoloFood/bacteria-archaea-traits/output/condensed_traits_GTDB.csv",stringsAsFactors = F)

gtdbtk_summary_ar <- read.delim(gtdbtk_ar,stringsAsFactors = F)
gtdbtk_summary_bac <- read.delim(gtdbtk_bac,stringsAsFactors = F)


# Initialize data for matching
ar_data <- as.data.frame(cbind(gtdbtk_summary_ar[,1],str_split_fixed(gtdbtk_summary_ar$classification,";",7)))
ar_data <- as.data.frame(sapply(ar_data,sub,pattern = "[a-z]{1}__",replacement=""))

bac_data <- as.data.frame(cbind(gtdbtk_summary_bac[,1],str_split_fixed(gtdbtk_summary_bac$classification,";",7)))
bac_data <- as.data.frame(sapply(bac_data,sub,pattern = "[a-z]{1}__",replacement=""))


mag_data <- as.data.frame(rbind(bac_data,ar_data))
colnames(mag_data) <- c("MAG_ID","superkingdom","phylum","class","order","family","genus","species")
mag_data[mag_data == ""] <- NA

# Split mag_data data frame into small df
by_species <- subset(mag_data, !is.na(mag_data$species))
by_genus <- subset(mag_data, is.na(mag_data$species) & !is.na(mag_data$genus))
by_family <- subset(mag_data, is.na(mag_data$genus))


# Find traits for MAGs given taxonomy
# traits columns 10-27
traits_byspecies<- traits[,c(5:6,10:27)][match(by_species$species, traits$species),]
by_species <- cbind(by_species,traits_byspecies)

traits_bygenus <- traits[,c(5:6,10:27)][match(by_genus$genus, traits$genus),]
by_genus <- cbind(by_genus,traits_bygenus)

traits_byfamily <- traits[,c(5:6,10:27)][match(by_family$family, traits$family),]
by_family <- cbind(by_family,traits_byfamily)



output <- rbind(by_species,by_genus,by_family)
write.csv(x = output,file = out_file,quote = F,sep = "\t",col.names = F,row.names = F)

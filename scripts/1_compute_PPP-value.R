################################################################################
#
# R-script: Plotting the MAP trees from the partitioned UCE analyses
#
#
# authors: Sebastian Hoehna
#
################################################################################

library(RColorBrewer)

# Get arguments
args = commandArgs(trailingOnly=TRUE)

BASE_DIR = args[1]

OUTPUT_DIR = paste0(BASE_DIR, "/results")
FIGS_DIR   = paste0(BASE_DIR, "/figures")

MISSING_SEQ_TYPE = "without"
MISSING_SEQ_TYPE = "with"

data <- read.table(file=paste0(OUTPUT_DIR,"/data_summary_",MISSING_SEQ_TYPE,"_missing.csv"),sep=",",header=TRUE)

col.names <- c("multinomialProfileLikelihood",
               "numInvSitesWithAmbiguous",
               "numInvSitesWithoutAmbiguous",
               "numVarSitesWithAmbiguous",
               "numVarSitesWithoutAmbiguous",
               "maxGcWithAmbiguous",
               "maxGcWithoutAmbiguous",
               "maxInvBLWithAmbiguous",
               "maxInvBLWithoutAmbiguous",
               "maxPDWithAmbiguous",
               "maxPDWithoutAmbiguous",
               "maxVarBLWithAmbiguous",
               "maxVarBLWithoutAmbiguous",
               "meanGcWithAmbiguous",
               "meanGcWithoutAmbiguous",
               "minGcWithAmbiguous",
               "minGcWithoutAmbiguous",
               "minPDWithAmbiguous",
               "minPDWithoutAmbiguous",
               "numInvBlocksWithAmbiguous",
               "numInvBlocksWithoutAmbiguous",
               "varGcWithAmbiguous",
               "varGcWithoutAmbiguous")

PP_p_values = matrix(NA,nrow=436,ncol=length(col.names)+1)
colnames(PP_p_values) <- c("Locus",col.names)

for ( ds in 1:436 ) {

    fn = paste0(OUTPUT_DIR,"/PPS_data_summary_Lampyridae_",ds,"_",MISSING_SEQ_TYPE,"_missing.csv")
    if ( file.exists(fn) == FALSE ) next

    PPS <- read.table(file=fn,sep=",",header=TRUE)

    for ( i in 1:length(col.names) ) {
        var.name = col.names[i]
        sim_values <- PPS[[var.name]]
        emp_value <- as.numeric(data[ds,][[var.name]])
        p_value <- mean(sim_values > emp_value) + 0.5 * mean(sim_values == emp_value)
        PP_p_values[ds,i+1] <- p_value
    }
    PP_p_values[ds,1] <- ds

}
write.table(PP_p_values,file=paste0(OUTPUT_DIR,"/PPS_",MISSING_SEQ_TYPE,"_missing.csv"),sep=",",row.names=FALSE)

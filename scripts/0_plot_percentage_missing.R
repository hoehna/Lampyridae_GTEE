################################################################################
#
# R Script: Plotting the sequence coverage
#
# authors: Sebastian Hoehna
#
################################################################################

CSV_FILE = "data_AHE.csv"

library("plot.matrix")
library("gdata")

FIGS_DIR = "0_data_exploration/figures"

## Create figures directory
dir.create(FIGS_DIR, showWarnings=FALSE)

data <- read.csv(CSV_FILE)


m <- matrix(NA, ncol=(ncol(data)-1), nrow=nrow(data))
for ( i in 1:nrow(data) ) {
  for ( j in 1:(ncol(data)-1) ) {
    m[i,j] <- 1-data[i,j+1]
  }
}
colnames(m) <- names(data)[-1]
disp_names <- c()
for (i in 1:length(colnames(m))) {
    disp_names[i] <- ifelse( startsWith(colnames(m)[i], "Photinus", trim=FALSE, ignore.case=FALSE), "P", "" )
    disp_names[i] <- ifelse( startsWith(colnames(m)[i], "Ellychnia", trim=FALSE, ignore.case=FALSE), "E", disp_names[i] )
}
pdf(paste0(FIGS_DIR,"/data_AHE_missing.pdf"),width=8,height=8)
par(mar=c(3.1, 4.1, 4.1, 4.1)) # adapt margins
plot(m, digits=NA, cex=2.5, ylab="", xlab="", main="Percentage of sequence missing for UCE dataset",col=gray.colors(n=10, start = 0.0, end = 1.0, gamma = 2.2, rev = TRUE), border=NA, breaks=10, fmt.key="%.1f", axis.col=NULL, axis.row=NULL)
axis(2, lwd.tick=1, at=436-seq(50,436,by=50), labels=seq(50,436,by=50), lwd=0)
axis(1, lwd.tick=1, at=which(disp_names != ""), labels=disp_names[disp_names != ""], lwd=0)
mtext("Taxon", side=1, line=1.0, cex=1.25)
mtext("UCE", side=2, line=2.5, cex=1.25)
dev.off()

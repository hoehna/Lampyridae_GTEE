################################################################################
#
# Function to check for convergence for continuous parameters.
#
#
# author: Sebastian Hoehna
#
################################################################################

library(convenience)

source("scripts/convergence_helper.R")

# Get arguments
args = commandArgs(trailingOnly=TRUE)

cat(args,"\n")

# Analysis settings
N_REPS       = as.numeric(args[1])
OUTPUT_DIR   = args[2]
FIGS_DIR     = args[3]
DS           = args[4]

## Create figures directory
dir.create(FIGS_DIR, showWarnings=FALSE)

treefile_names = paste0(OUTPUT_DIR,"/Lampyridae_unrooted_",DS,"_run_",1:N_REPS,".trees")
logfile_names  = paste0(OUTPUT_DIR,"/Lampyridae_unrooted_",DS,"_run_",1:N_REPS,".log")
format         = "revbayes"

conv.control = makeControl( tracer = NULL,
                            burnin = 0.2,
                            precision = NULL,
                            namesToExclude = NULL
                           )

check_conv <- checkConvergence( list_files = c(logfile_names, treefile_names ),
                                control = conv.control,
                                format=format )

pdf( paste0(FIGS_DIR,"/convergence_",DS,".pdf"), height=8, width=8 )
  par(mfrow=c(2,2))

  plotEssContinuous(check_conv)
  plotKS(check_conv)
  plotEssSplits(check_conv)
  plotDiffSplits(check_conv)
dev.off()

save(check_conv, file = paste0(FIGS_DIR,"/convergence_",DS,".Rdata") )

SOFTWARE       <- "RevBayes"
chain_names    <- paste0("Run ",1:N_REPS)
plot_name      <- paste0(FIGS_DIR,"/split_frequencies_",DS,".pdf")
plot.split.frequencies( treefile_names, chain_names, plot_name, tolower(SOFTWARE) )

################################################################################
#
# Function to plot the split frequencies of different runs.
#
#
# author: Sebastian Hoehna
#
################################################################################

library(convenience)


plot.split.frequencies <- function( filenames, chain_names, plot.filename, format="revbayes" ) {

  clade_names <- c()
  splits <- list()
  index <- 1

  for (f in filenames) {
    my_trees <- list()
    my_trees[[1]] <- loadTrees(file=f, format = format)
    tmp_splits <- splitFreq(my_trees, windows=FALSE)

    for (j in 1:length(tmp_splits[1,])) {
      this_clade_name = as.character(tmp_splits[1,j])
      if ( sum( clade_names == this_clade_name ) == 0 ) {
        clade_names <- c(clade_names, this_clade_name)
      }
    }

    splits[[index]] <- tmp_splits

    index <- index+1
  }

  num_traces <- length( splits )
  split_freqs <- matrix(0,length(clade_names),num_traces)
  for (i in 1:num_traces) {
    tmp_splits <- splits[[i]]
    for (j in 1:length(tmp_splits[1,])) {
      this_clade_name = as.character(tmp_splits[1,j])
      index = which( clade_names == this_clade_name )
      pp  = as.numeric(tmp_splits[2,j])
      split_freqs[index,i] = pp
    }
  }

  layout_mat = t(matrix(1:(num_traces*num_traces), nrow=num_traces))
#filename <- paste0(RESULTS_DIR,"/",GENE_NAME,".pdf")
#filename <- paste0("0_mcmc_tree_moves/Splits.pdf")
  pdf(plot.filename, height=2.4*num_traces, width=2.4*num_traces)
  layout(layout_mat)
  par(mar=c(0,0.0,0.0,0.0), oma=c(4,4.0,5,0.1))

  for (i in 1:num_traces) {
    for (j in 1:num_traces) {
      plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", las=1, bty="n", xaxt="n", yaxt="n")

      box()
      if ( i==j ) {
        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
      } else {
        abline(a=0,b=1)
        points( x=split_freqs[,i], y=split_freqs[,j], pch=16, cex=2.5 )

        max_diff <- max( abs(split_freqs[,i] - split_freqs[,j]) )
        text( 0, 0.9, sprintf("M = %.2f", max_diff) )
      }
      if ( i == 1 ) {
        mtext( chain_names[j], side=3, line=1, cex=1.5 )
      }
      if ( j == 1 ) {
        mtext( chain_names[i], side=2, line=1, cex=1.5 )
      }

    }
  }

  dev.off()
}


plot.diff.splits <- function(output, plot.y.axis = FALSE, plot.x.axis = FALSE, plot.legend = FALSE) {

  minimumESS <- 625
  per_run = FALSE
  col_threshold <- "gray69"
  fill_color <- "seagreen4"
  fdir <- system.file("thresholds/expectedDiff_625.rds", package = "convenience")
  exp_diff_runs <- readRDS(fdir)

  ## Calculate the minimum ESS between runs for each split
  ess_min_between_runs <- matrix(ncol=length(output$tree_parameters$frequencies), nrow=nrow(output$tree_parameters$ess))
  n_runs <- length(output$tree_parameters$ess)
  count <-1
  for (i in 1:(n_runs-1)) {
    for (j in (i+1):n_runs) {
      for (z in 1:nrow(output$tree_parameters$ess)) {
        ess_min_between_runs[z,count] <- min(output$tree_parameters$ess[z,i], output$tree_parameters$ess[z,j], na.rm = T)
      }
      count <- count + 1
    }
  }
  row.names(ess_min_between_runs) <- row.names(output$tree_parameters$ess)

  list_freq <- list()
  list_diff <- list()
  list_freq_low_ess <- list()
  list_diff_low_ess <- list()
  for (i in 1:length(output$tree_parameters$frequencies)) {
    frequencies <- vector()
    differences <- vector()
    freq_low_ess <- vector()
    diff_low_ess <- vector()
    for (j in names(output$tree_parameters$frequencies[[i]])) {
      if ( j %in% row.names(ess_min_between_runs)){
        if ( ess_min_between_runs[j, i] < minimumESS ){
          freq_low_ess <- c( freq_low_ess, as.vector(unlist(output$tree_parameters$frequencies[[i]][j])) )
          diff_low_ess <- c( diff_low_ess, as.vector(unlist(output$tree_parameters$compare_runs[[i]][j])) )
        }
        else{
          frequencies <- c( frequencies, as.vector(unlist(output$tree_parameters$frequencies[[i]][j])) )
          differences <- c( differences, as.vector(unlist(output$tree_parameters$compare_runs[[i]][j])) )
        }
      }
      else{
        frequencies <- c( frequencies, as.vector(unlist(output$tree_parameters$frequencies[[i]][j])) )
        differences <- c( differences, as.vector(unlist(output$tree_parameters$compare_runs[[i]][j])) )
      }
    }
    list_freq[[i]] <- frequencies
    list_diff[[i]] <- differences
    list_freq_low_ess[[i]] <- freq_low_ess
    list_diff_low_ess[[i]] <- diff_low_ess

  }

  x_axis <- exp_diff_runs[1,]
  y_axis <- exp_diff_runs[2,]

  y_lim <- max(unlist(list_diff), y_axis, unlist(list_diff_low_ess), na.rm = T)
  y_lim <- y_lim + y_lim*0.5
  y_lim <- 1

  plot <- plot(NA,
               xlab = NA,
               ylab = NA,
               xaxt="n",
               yaxt="n",
               cex.main = 0.9,
               xlim = c(0.0,1.0),
               ylim = c(0.0, y_lim),
               las = 1)

  plot <- polygon(x_axis, y_axis, col = "gray89", border = NA)
  x_extra <- c(0.0, 0.01, 0.99, 1.0, 0.0)
  y_extra <- c(0.0, min(y_axis), min(y_axis), 0.0, 0.0)

  plot <- polygon(x_extra, y_extra, border = NA, col = "gray89")
  plot <- lines(x_axis, y_axis, col = col_threshold, lwd=2)

  if ( plot.x.axis ) axis(1, lwd.tick=1, lwd=0, cex.axis=1.1)
  if ( plot.y.axis ) axis(2, lwd.tick=1, lwd=0, cex.axis=1.1)

#  if( is.null(xlab)) xlab <- "Split frequency" else xlab <- xlab
#  if( is.null(ylab)) ylab <- "Difference between splits" else ylab <- ylab
#  title(xlab = xlab, outer = F, line = 2.5)
#  title(ylab = ylab, outer = F, line = 3.5)


  if( length(output$tree_parameters$frequencies) > 1 & per_run == TRUE ){

    for (i in 1:length(output$tree_parameters$frequencies)) {
      plot <- points(unlist(list_freq[i]), unlist(list_diff[i]), pch = i, col = fill_color)
      if (length(list_freq_low_ess[i]) > 0){
        plot <- points(unlist(list_freq_low_ess[i]), unlist(list_diff_low_ess[i]), pch = i, col = "tan1")
        legend("topleft",
               legend = c(paste("ESS",expression("<"),minimumESS),paste("ESS",expression(">"),minimumESS)),
               pch = c(16,16),
               col = c("tan1",fill_color),
               box.lty = 2,
               box.col="gray7",
               cex = 0.7,
               inset = 0.01)
      }
      legend("topright",
             legend = names(output$tree_parameters$frequencies),
             pch = c(1:length(output$tree_parameters$frequencies)),
             box.lty = 2,
             box.col = "gray7",
             cex = 0.7,
             inset = 0.01)
    }

  } else{
    frequencies <- unlist(list_freq)
    differences <- unlist(list_diff)
    freq_low_ess <- unlist(list_freq_low_ess)
    diff_low_ess <- unlist(list_diff_low_ess)
    plot <- points(frequencies, differences, pch = 16, col = fill_color)
    if (length(freq_low_ess) > 0){
      plot <- points(freq_low_ess, diff_low_ess, pch = 16, col = "tan1")
      if ( plot.legend == TRUE ) {
        legend("topright",
               legend = c(paste("ESS",expression("<"),minimumESS),paste("ESS",expression(">"),minimumESS)),
               pch = c(16,16),
               col = c("tan1",fill_color),
               box.lty = 2,
               box.col="gray7",
               cex = 0.7,
               inset = 0.01)
      }

    }

  }

}

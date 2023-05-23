#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



# test if there is at least one argument: if not, return an error
if (length(args) !=2) {
  stop("Please add first the steem root of the treemix output as first argument and the pop order as second.", call.=FALSE)
  }
#add the plotting functions used by treemix
source("/ebio/abt6_projects8/Tarvense_TE/code/genetic_analysis/plotting_funcs.R")

steem=args[1]
pop_order=args[2]

#name the pdf to be written the plot
pdf(file=paste0("treemix.tree.",steem,".pdf"),12,8)
  #plot tree
  plot_tree(steem)
  dev.off()

#name the pdf to be written the plot
pdf(file=paste0("treemix.residuals.",steem,".pdf"),12,8)
   #plot residuals
  plot_resid(steem, pop_order = pop_order)
  dev.off()

#plot_tree("/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/treemix_results_nomig_k100", plotnames = )
#plot_resid("/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/",
#           pop_order="/ebio/abt6_projects8/Tarvense_TE/data/Tarvense_ancestry/treemix/pop_order.txt")

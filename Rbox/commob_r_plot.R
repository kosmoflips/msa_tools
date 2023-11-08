# ////////////////////////////////////////////////////
# ----- commonly used graphic libs -----
# ////////////////////////////////////////////////////
library(ggplot2)
library(cowplot) # for place multiple plots in one file
library(grid) # for draw legend
#library(scales)
#library("RColorBrewer")


# draw a plot with title only
mktitle <- function(titlestring='') {
	plottitle <- ggdraw() + 
		draw_label( titlestring,
								fontface = 'bold',
								x = 0,
								hjust = 0
		) +
		theme(
			# add margin on the left of the drawing canvas,
			# so title is aligned with left edge of first plot
			plot.margin = margin(0, 0, 0, 7)
		)
	return (plottitle)
}


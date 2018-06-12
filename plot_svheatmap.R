#!/usr/bin/env Rscript

###
# Plot structure variants heatmap in multiple genomes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The input is a CSV file in which columns corresponds to each genome and rows
# are the location of variants.
#
# copyright: (c) 2018 by Ali Ghaffaari.
# license: MIT, see LICENSE for more details.
###

library('tibble')
library('scales')
library('plyr')
library('reshape2')
library('ggplot2')
library('optparse')

# Default values
in_file <- 'all_sv_loci.csv'
p_title <- 'Distribution of structural variations throught out the genome'
p_subtitle  <- ''
nofbins <- 1000
base_size <- 9
option_list <- list(make_option(c('-f', '--file'), type='character', default=in_file, help='Input CSV file\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-t', '--title'), type='character', default=p_title, help='Plot title\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-s', '--subtitle'), type='character', default=p_subtitle, help='Plot subtitle pattern\n\t\t[default=\'%default\']', metavar='character'),
                    make_option(c('-n', '--binsno'), type='integer', default=nofbins, help='Number of bins\n\t\t[default=\'%default\']', metavar='integer'),
                    make_option(c('-b', '--basesize'), type='integer', default=base_size, help='Font base size\n\t\t[default=\'%default\']', metavar='integer'))
opt_parser <- OptionParser(option_list=option_list)
# Parsing command-line arguments
opt <- parse_args(opt_parser)

# Reading input file
d <- read.table(opt$file, sep=',', header=TRUE)
minloc <- min(d, na.rm=TRUE)
maxloc <- max(d, na.rm=TRUE)
bins <- seq(0, maxloc, length=nofbins+1)
h <- tibble(bin=1:nofbins)
for (g in names(d)) {
    h[,g] <- as.data.frame(table(cut(d[,g], bins, include.lowest=TRUE)))$Freq
}
h <- melt(h, id=c('bin'))
h <- ddply(h, .(variable), transform, rescale=rescale(value))
g <- ggplot(h, aes(variable, bin, fill=rescale)) +
     geom_tile() +
     labs(title=opt$title,
          subtitle=opt$subtitle,
          x='',
          y='') +
     scale_fill_gradient(low='white', high='maroon') +
     scale_x_discrete(expand=c(0, 0)) +
     scale_y_discrete(expand=c(0, 0)) +
     theme_gray(base_size=base_size) +
     theme(legend.position="none", axis.ticks=element_blank(),
           axis.text.x=element_text(size=rel(0.8), angle=330, hjust=0, colour="grey50"))
ggsave('svheatmap.pdf')

save_pheatmap_pdf_size <- function(x, filename, width, height) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

save_pheatmap_png_size <- function(x, filename, width, height) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   png(filename, width=width, height=height, units = 'in', res = 300)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

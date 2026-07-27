#!/usr/bin/env Rscript
## combine_genome_wide_plots.R
##
## Combine ichorCNA Pass 1 and Pass 2 genomeWide plots into a single
## multi-page PDF for side-by-side QC review.
##
## Layout per A4 portrait page:
##   3 rows × 2 columns = 6 plots per page
##   Left  column : Pass 1 plot
##   Right column : Pass 2 plot
##   Each row     : one sample
##
## Usage:
##   Rscript combine_genome_wide_plots.R <output.pdf> \
##       <p1_file_1> ... <p1_file_N>  \
##       <p2_file_1> ... <p2_file_N>
##
## Positional arguments:
##   args[1]          : output PDF path
##   args[2 .. N+1]   : Pass 1 genomeWide PDF files  (N files)
##   args[N+2 .. 2N+1]: Pass 2 genomeWide PDF files  (N files, same sample order)

suppressPackageStartupMessages({
    library(pdftools)
    library(png)
    library(grid)
    library(gridExtra)
})

# ---- Parse arguments --------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3L || (length(args) - 1L) %% 2L != 0L) {
    stop(
        "Usage: combine_genome_wide_plots.R output.pdf ",
        "p1_file_1 ... p1_file_N  p2_file_1 ... p2_file_N\n",
        "(number of pass1 and pass2 files must be equal)"
    )
}

output_file <- args[1L]
n           <- (length(args) - 1L) / 2L
p1_files    <- args[seq(2L,       n + 1L)]
p2_files    <- args[seq(n + 2L, 2L * n + 1L)]

# ---- Helper: extract sample ID from filename --------------------------------

sample_id_from_p1 <- function(f) sub("_pass1_genomeWide\\.pdf$", "", basename(f))
sample_id_from_p2 <- function(f) sub("_genomeWide\\.pdf$",       "", basename(f))

# ---- Helper: render one PDF page to a rasterGrob ---------------------------
#
# pdf_render_page(..., numeric = TRUE) returns an array of dimensions
# [channels(4), width, height] with values 0–255.
# We transpose to [height, width, 3] (RGB only) and scale to [0, 1] before
# passing to grDevices::as.raster().

render_pdf_grob <- function(pdf_path, label = NULL, dpi = 100L) {
    # Convert the PDF page to a PNG in a temp file, then read it back with
    # png::readPNG.  This avoids all class/scaling ambiguity in pdftools'
    # bitmap/nativeRaster return types across versions.
    tmp <- tempfile(fileext = ".png")
    on.exit(unlink(tmp), add = TRUE)
    pdftools::pdf_convert(pdf_path, format = "png", pages = 1L,
                          filenames = tmp, dpi = dpi, verbose = FALSE)
    img <- png::readPNG(tmp)          # [height, width, channels], 0-1

    # Use absolute inch dimensions so that aspect ratio is preserved
    # regardless of viewport nesting (e.g. the arrangeGrob label wrapper
    # creates an inner viewport whose npc != the outer cell's npc, which
    # caused horizontal stretching when npc units were used for both axes).
    #
    # Available area for the image after subtracting the page-level title bar,
    # the per-plot label line, and a small inter-cell margin.
    cell_w_in <- A4_W / 2            - CELL_MARGIN
    cell_h_in <- (A4_H - PAGE_HEADER_H) / ROWS_PER_PAGE - PLOT_LABEL_H - CELL_MARGIN

    img_asp <- ncol(img) / nrow(img)  # image width / height
    w_in    <- cell_w_in
    h_in    <- w_in / img_asp
    if (h_in > cell_h_in) {           # taller than cell: constrain height instead
        h_in <- cell_h_in
        w_in <- h_in * img_asp
    }

    g <- grid::rasterGrob(img,
                           width       = grid::unit(w_in, "inches"),
                           height      = grid::unit(h_in, "inches"),
                           just        = "centre",
                           interpolate = TRUE)
    if (!is.null(label)) {
        g <- gridExtra::arrangeGrob(
            g,
            top = grid::textGrob(
                label,
                gp = grid::gpar(fontsize = 7, col = "grey25")
            )
        )
    }
    g
}

# ---- Write multi-page PDF ---------------------------------------------------

ROWS_PER_PAGE  <- 3L
A4_W           <- 8.27    # page width,  inches
A4_H           <- 11.69   # page height, inches
PAGE_HEADER_H  <- 0.40    # height of the page-level title bar (approx), inches
PLOT_LABEL_H   <- 0.20    # height of the per-plot sample-id label (approx), inches
CELL_MARGIN    <- 0.05    # padding between adjacent cells, inches

message("Combining ", n, " sample pair(s) into ", output_file, " ...")

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
grDevices::pdf(output_file, width = A4_W, height = A4_H)

for (page_start in seq(1L, n, by = ROWS_PER_PAGE)) {

    idx   <- seq(page_start, min(page_start + ROWS_PER_PAGE - 1L, n))
    grobs <- vector("list", 6L)

    for (k in seq_along(idx)) {
        j       <- idx[k]
        row_pos <- 2L * k - 1L   # positions: 1, 3, 5 (left = pass1)
        grobs[[row_pos]]     <- render_pdf_grob(
            p1_files[j],
            label = paste0(sample_id_from_p1(p1_files[j]), "  [Pass 1]")
        )
        grobs[[row_pos + 1L]] <- render_pdf_grob(
            p2_files[j],
            label = paste0(sample_id_from_p2(p2_files[j]), "  [Pass 2]")
        )
    }

    # Pad incomplete last page with blank grobs
    for (i in seq_along(grobs)) {
        if (is.null(grobs[[i]])) grobs[[i]] <- grid::nullGrob()
    }

    gridExtra::grid.arrange(
        grobs  = grobs,
        ncol   = 2L,
        nrow   = 3L,
        top    = grid::textGrob(
            paste0(
                "ichorCNA 2-pass genomeWide plots",
                "     \u2502     left: Pass 1   \u2022   right: Pass 2"
            ),
            gp = grid::gpar(fontsize = 10, fontface = "bold")
        )
    )

    # Release memory between pages
    rm(grobs)
    invisible(gc(verbose = FALSE))
}

grDevices::dev.off()
message("Done.")

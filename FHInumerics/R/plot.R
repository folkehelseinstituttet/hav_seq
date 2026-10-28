# Tree plot with metadata labels and colouring.

plot_tree <- function(res, colour_by = "none", label_cols = character(0), show_support = TRUE) {
  meta <- res$tip_meta
  meta$label_text <- meta$id
  for (col in label_cols) {
    val <- ifelse(meta$is_query, "", ifelse(is.na(meta[[col]]), "?", meta[[col]]))
    meta$label_text <- ifelse(meta$is_query, meta$label_text, paste(meta$label_text, val, sep = " | "))
  }
  meta$label_text <- ifelse(meta$is_query, paste0(meta$id, "  (QUERY)"), meta$label_text)

  label_offset <- 0.015 * max(ape::node.depth.edgelength(res$tree))
  p <- ggtree::ggtree(res$tree) %<+% meta
  if (colour_by != "none") {
    p <- p + ggtree::geom_tippoint(ggplot2::aes(colour = .data[[colour_by]]), size = 2.5,
                                   data = function(d) d[d$isTip & !d$is_query, ]) +
      ggplot2::labs(colour = colour_by)
  }
  p <- p +
    ggtree::geom_tippoint(data = function(d) d[d$isTip & d$is_query, ], shape = 18, size = 4, colour = "#d7301f") +
    ggtree::geom_tiplab(ggplot2::aes(label = label_text), data = function(d) d[d$isTip & !d$is_query, ],
                        size = 3.2, offset = label_offset) +
    ggtree::geom_tiplab(ggplot2::aes(label = label_text), data = function(d) d[d$isTip & d$is_query, ],
                        size = 3.4, offset = label_offset, colour = "#d7301f", fontface = "bold") +
    ggtree::geom_treescale(fontsize = 3) +
    ggtree::hexpand(0.6) +
    ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 10))
  if (show_support && !is.null(res$tree$node.label)) {
    p <- p + ggtree::geom_text2(ggplot2::aes(subset = !isTip & suppressWarnings(as.numeric(label)) >= 70,
                                             label = label), hjust = 1.2, vjust = -0.4, size = 2.6, colour = "grey40")
  }
  # ggtree's print method only drops this class and prints; without it, Shiny sees a
  # ggplot and reports mouse positions in the tree's own coordinates (needed for hover).
  class(p) <- setdiff(class(p), "ggtree")
  p
}

`%<+%` <- ggtree::`%<+%`

# ---- Hover --------------------------------------------------------------------------
# Shiny's hover input gives the pointer position in data units and in CSS pixels, plus
# the panel's data range (domain) and pixel range. These convert plot points to CSS
# pixels, so the window can find the point under the pointer and draw a ring around it.
data_to_css <- function(x, y, hover) {
  dom <- hover$domain
  rng <- hover$range
  ratio <- hover$img_css_ratio %||% list(x = 1, y = 1)
  list(x = (rng$left + (x - dom$left) / (dom$right - dom$left) * (rng$right - rng$left)) / ratio$x,
       y = (rng$bottom + (y - dom$bottom) / (dom$top - dom$bottom) * (rng$top - rng$bottom)) / ratio$y)
}

# The point under the pointer: the nearest point whose centre is within `radius` CSS
# pixels (one value, or one per point). Returns list(i, x, y, r) in CSS pixels, or NULL.
point_at <- function(x, y, hover, radius = 8) {
  if (is.null(hover$coords_css) || is.null(hover$domain) || !length(x)) return(NULL)
  pos <- data_to_css(x, y, hover)
  dist <- sqrt((pos$x - hover$coords_css$x)^2 + (pos$y - hover$coords_css$y)^2)
  radius <- rep_len(radius, length(x))
  hit <- which(dist <= radius)
  if (!length(hit)) return(NULL)
  i <- hit[which.min(dist[hit])]
  list(i = i, x = pos$x[i], y = pos$y[i], r = radius[i])
}

# Radius in CSS pixels of the circles drawn by layer `layer` of plot `p` (geom_point,
# solid circles). Shiny draws plots at 72 dpi per CSS pixel, so 1 pt = 1 CSS pixel.
point_radius_css <- function(p, layer) {
  d <- ggplot2::layer_data(p, layer)
  0.375 * (d$size * ggplot2::.pt + (d$stroke %||% 0.5) * ggplot2::.stroke / 2)
}

# The tree tip whose point is under the pointer. `tree_data` is the plot's data (p$data).
# Returns list(label, x, y, r) or NULL.
tip_at <- function(tree_data, hover) {
  tips <- tree_data[tree_data$isTip, ]
  hit <- point_at(tips$x, tips$y, hover, radius = 8)
  if (is.null(hit)) return(NULL)
  c(list(label = tips$label[hit$i]), hit[c("x", "y", "r")])
}

# What the hover box shows for one tip: a named character vector.
# Database samples: all their metadata. Queries: the closest sequence(s) in the tree.
tip_details <- function(res, label) {
  m <- res$tip_meta[res$tip_meta$id == label, , drop = FALSE]
  if (!nrow(m)) return(c(Sample = label))
  if (isTRUE(m$is_query[1])) {
    snp <- res$snp[label, setdiff(colnames(res$snp), res$queries), drop = TRUE]
    if (!length(snp)) return(c(Query = label))
    closest <- names(snp)[snp == min(snp)]
    more <- if (length(closest) > 3) sprintf(" and %d more", length(closest) - 3) else ""
    return(c(Query = label, `Closest in tree` = paste0(paste(utils::head(closest, 3), collapse = ", "), more),
             SNPs = as.character(min(snp))))
  }
  vals <- unlist(as.list(m[1, setdiff(names(m), c("id", "is_query")), drop = FALSE]))
  c(Sample = label, ifelse(is.na(vals), "?", vals))
}

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
  p
}

`%<+%` <- ggtree::`%<+%`

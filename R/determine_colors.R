
determine_colors = function(factors,nclass){
  library(RColorBrewer)
  # Set base color palette
  base_palette <- colorRampPalette(brewer.pal(8, "Set2"))(10)
  
  row.colors <- base_palette[seq_along(factors)]
  
  
  # Dynamically extend color palette if needed
  if (nclass > length(row.colors)) {
    extras_needed <- nclass - length(row.colors)
    remaining_colors <- base_palette[!base_palette %in% row.colors]
    
    additional_colors <- if (length(remaining_colors) < extras_needed) {
      colorRampPalette(brewer.pal(8, "Set3"))(extras_needed)
    } else {
      remaining_colors[1:extras_needed]
    }
    
    col.colors <- c(row.colors, additional_colors)
  } else {
    col.colors <- row.colors[1:nclass]
  }
  
  return(list(row.colors=row.colors,col.colors=col.colors))
}


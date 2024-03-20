plot_contribution_relative_bar <- function (contribution, signatures = NA, colorsig = NA, index = NA) {
  if (!is.na(index)) {
    contribution <- contribution[, index, drop = FALSE]
  }
  if (is.na(colorsig)) {
    colorsig <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                  "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                  "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                  "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                  "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                  "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
  }
  if (length(colorsig) < dim(contribution)[1]) {
    stop(paste("Provide colors vector with length ", as.character(dim(mut_mat)[2])), call. = FALSE)
  }
  
  tb <- contribution %>% as.data.frame() %>% tibble::rownames_to_column("Signature") %>% 
    tidyr::pivot_longer(-Signature, names_to = "Sample", 
                        values_to = "Contribution") %>% dplyr::mutate(Sample = factor(Sample, 
                                                                                      levels = unique(Sample)), Signature = factor(Signature, 
                                                                                                                                   levels = unique(Signature)))
  
  bar_geom <- geom_bar(position = "fill", stat = "identity", 
                       colour = "black")
  y_lab <- "Relative contribution"
  
  present_sigs <- tb %>% dplyr::filter(Contribution != 0) %>% 
    dplyr::pull(Signature) %>% unique()
  plot <- ggplot(tb, aes(x = Sample, y = Contribution, fill = Signature)) + 
    bar_geom + 
    labs(x = "", y = y_lab) +
    scale_fill_manual(values = colorsig)+
    theme_classic() + 
    theme(axis.ticks.x = element_line(size=1),
          axis.ticks.y =  element_line(size=1),
          axis.line.x = element_line(size=1),
          axis.line.y = element_line(size=1),
          strip.background.x = element_blank(),
          strip.background.y = element_blank(),
          axis.text.x = element_text(angle = 90, size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15))+
    scale_y_continuous(expand = c(0,0))
  
  return(plot)
}
plot_96_profile_M <- function (mut_matrix, colors = NA, colorsig = NA, size, condensed = FALSE) 
{
  # size
  # 1. signature?????Ŀ??ȣ?
  # 2. signature??????С
  # 3. ??ͷType?????·?ͼ?εľ???
  # 4. ??ͷType?Ĵ?С
  # 5. ??ͷ?Ĵ?С
  # 6. x????ǩ?��Ĵ?С
  freq <- full_context <- substitution <- context <- NULL
  if (is.na(colors)) {
    colors <- c("#5ABD85", "#36637A", "#FAC076", "#3B3D20", "#886DE0", "#DF8574")
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }
  if (is.na(colorsig)) {
    colorsig <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
                  "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                  "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
                  "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
                  "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
                  "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")
  }
  if (length(colors) < dim(mut_matrix)[2]) {
    stop(paste("Provide colors vector with length ", as.character(dim(mut_mat)[2])), call. = FALSE)
  }
  norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
  tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>% 
    dplyr::mutate(substitution = stringr::str_replace(full_context, "\\w\\[(.*)\\]\\w", "\\1"), 
                  context = stringr::str_replace(full_context,"\\[.*\\]", "\\.")) %>%
    dplyr::select(-full_context) %>% 
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample",values_to = "freq") %>%
    dplyr::mutate(sample = factor(sample,levels = unique(sample)))
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }  else {
    width <- 0.6
    spacing <- 0.5
  }
  
  SIGNATURE <- unique(tb$sample)
  maxfreq <- max(tb$freq)
  maxcont <- length(unique(tb$context))
  p <- list()
  for (i in 1:length(SIGNATURE)) {
    s <- SIGNATURE[i]
    tbs <- tb[tb$sample == s, ]
    tbs$label <- factor(rownames(mut_matrix), levels =  rownames(mut_matrix))
    maxnumber <- dim(tbs)[1]
    p[[i]] <- ggplot(data = tbs, aes(x = label, y = freq, fill = substitution)) + 
      geom_bar(stat = "identity", size = 1) +
      scale_fill_manual(values = colors) +
      xlab("")+ylab("")+
      coord_cartesian(ylim = c(0, maxfreq))+
      theme_classic()+ # 坐标轴的粗细
      theme(panel.spacing.x = unit(0, "cm"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 15),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(size=1),
            axis.ticks.length = unit(0.3, "cm"),
            panel.grid.minor.y = element_line( color="black", linetype  = 2 ),
            axis.line = element_line( color="black"),
            legend.position = 'none',
            legend.box.background = element_rect(fill = NA, color= "black", linetype = 1),
            panel.background = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(size=1),
            strip.background.x = element_blank(),
            strip.background.y = element_blank(),
            strip.text.x = element_blank())+
      scale_x_discrete(expand = c(0.04,0))+#?ı?X??
      geom_rect(aes(xmin = maxnumber-20, xmax = maxnumber,
                    ymin = maxfreq-size[1], ymax = maxfreq),
                fill = colorsig[i])+
      annotate("text", x = maxnumber-10, y= maxfreq-size[1]/2, label = s, size = size[2])
  }
 
  
  
  tbs$freq <- 0
  pt <- ggplot(data = tbs, aes(x = label, y = freq, fill = substitution)) + 
    geom_bar(stat = "identity", size = 1) +
    scale_fill_manual(values = colors) +
    xlab("")+ylab("")+
    coord_cartesian(ylim = c(0, maxfreq/5))+
    theme_classic()+ # 坐标轴的粗细
    theme(panel.spacing.x = unit(0, "cm"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 15,color="white"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(size=1,color="white"),
          axis.ticks.length = unit(0.3, "cm"),
          # panel.grid.minor.y = element_line( color="black", linetype  = 2 ),
          axis.line = element_line( color="black"),
          legend.position = 'none',
          legend.box.background = element_rect(fill = NA, color= "black", linetype = 1),
          panel.background = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(size=1,color="white"),
          strip.background.x = element_blank(),
          strip.background.y = element_blank(),
          strip.text.x = element_blank())+
    scale_x_discrete(expand = c(0.04,0))+#?ı?X??
    geom_rect(aes(xmin = maxcont*0+1, xmax = maxcont-0.5,
                  ymin = 0, ymax = maxfreq/15),
              fill = colors[1])+
    geom_rect(aes(xmin = maxcont*1+1, xmax = maxcont*2-0.5,
                  ymin = 0, ymax = maxfreq/15),
              fill = colors[2])+
    geom_rect(aes(xmin = maxcont*2+1, xmax = maxcont*3-0.5,
                  ymin = 0, ymax = maxfreq/15),
              fill = colors[3])+
    geom_rect(aes(xmin = maxcont*3+1, xmax = maxcont*4-0.5,
                  ymin = 0, ymax = maxfreq/15),
              fill = colors[4])+
    geom_rect(aes(xmin = maxcont*4+1, xmax = maxcont*5-0.5,
                  ymin = 0, ymax = maxfreq/15),
              fill = colors[5])+
    geom_rect(aes(xmin = maxcont*5+1, xmax = maxcont*6-0.5,
                  ymin = 0, ymax = maxfreq/15),
              fill = colors[6])+
    annotate("text", x = maxcont/2+0.25, y=maxfreq/15+size[3] , label = TYPE[1], size = size[4])+
    annotate("text", x = maxcont+maxcont/2+0.25, y=maxfreq/15+size[3] , label = TYPE[1], size = size[4])+
    annotate("text", x = maxcont*2+maxcont/2+0.25, y=maxfreq/15+size[3] , label = TYPE[2], size = size[4])+
    annotate("text", x = maxcont*3+maxcont/2+0.25, y=maxfreq/15+size[3] , label = TYPE[3], size = size[4])+
    annotate("text", x = maxcont*4+maxcont/2+0.25, y=maxfreq/15+size[3] , label = TYPE[4], size = size[4])+
    annotate("text", x = maxcont*5+maxcont/2+0.25, y=maxfreq/15+size[3] , label = TYPE[5], size = size[4])
  
  
  tbs$freq <- 0
  px <- ggplot(data = tbs, aes(x = label, y = freq, fill = substitution)) + 
    geom_bar(stat = "identity", size = 1) +
    scale_fill_manual(values = colors) +
    xlab("")+ylab("")+
    coord_cartesian(ylim = c(0,maxfreq))+
    theme_classic()+ # 坐标轴的粗细
    theme(panel.spacing.x = unit(0, "cm"),
          axis.text.y = element_text(size = 15, co, color="white"  
          axis.text.x = element_text(angle = 90,hjust = 0.5, vjust = 0.5, margin = margin(-2,0,0,0,'cm')),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(size=1, color="white"),
          axis.ticks.length = unit(0.3, "cm"),
          legend.position = 'none',
          legend.box.background = element_rect(fill = NA, color= "black", linetype = 1),
          panel.background = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(color="white"),
          strip.background.x = element_blank(),
          strip.background.y = element_blank(),
          strip.text.x = element_blank())+
    scale_x_discrete(expand = c(0.04,0))+
    scale_y_continuous(expand = c(0,0)   expp <- c()
  hexpp <- c(size[5],rep(3,length(SIGNATURE)),size[6])
  for (i in 1:length(SIGNATURE)) {
    expp <- paste(expp,"p[[", as.character(i),"]],", sep = "", collapse = "")
  }
  
  p1 <- eval(parse(text = paste("grid.arrange(pt,", expp,"px,heights = hexpp)", sep = "", collapse = "") ))
  return(p1)
  
  
}
  
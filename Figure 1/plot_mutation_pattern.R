plot_mutation_pattern <- function (mut_mat, colors = NA, size, condensed = FALSE){
  
  # size: ???·????͵???????С
  freq <- full_context <- substitution <- context <- NULL
  if (is.na(colors)) {
    colors <- c("#5ABD85", "#36637A", "#FAC076", "#3B3D20", "#886DE0", "#DF8574")#c("#085A9C","#EF0808","#526373","#FFFFE7","#FF9418","#219431") #"#9C52AD"
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6", call. = FALSE)
  }
  
  norm_mut_matrix <- mut_mat #apply(mut_matrix, 2, function(x) x/sum(x))
  tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>% 
    dplyr::mutate(substitution = stringr::str_replace(full_context,"\\w\\[(.*)\\]\\w", " \\1 "), 
                  context = stringr::str_replace(full_context,"\\[.*\\]", "\\.")) %>% 
    dplyr::select(-full_context) %>% 
    tidyr::pivot_longer(c(-substitution, -context), names_to = "sample",values_to = "freq") %>%
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))
  
  group <- c()
  sample_names <- as.character(tb$sample)
  for (i in 1:length(sample_names)){
    s <- sample_names[i]
    group <- c(group, substr(s, nchar(s),nchar(s)))
  }
  group[group == "C"] <- "Cancer"
  group[group == "H"] <- "High grade"
  group[group == "L"] <- "Low grade"
  
  tb <- mutate(tb, group)
  
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }else {
    width <- 0.6
    spacing <- 0.5
  }
  
  TYPE <- unique(tb$substitution)
  p <- list()
  my_g0 <- grobTree(rectGrob(gp=gpar(col = "White",fill="white")))
  ymax <- 2000
  
  tbbar <- data.frame(c("0", "0"), c("Cancer", "High grade"), c(0,0))
  colnames(tbbar) <- c("context", "group", "freq")
  for (i in 1:length(TYPE)) {
    
    tbtype <- tb[tb$substitution == TYPE[i], ]
    tbtype <- aggregate(x = tbtype$freq, by=list(tbtype$context,tbtype$group),sum)
    colnames(tbtype) <- c("context", "group", "freq")
    tbtype$context <- paste(TYPE[i], tbtype$context)
    tbtype <- tbtype[order(tbtype$context), ]
    
    tbbar <- rbind(tbbar, tbtype)
    tbbar0 <- data.frame(c(as.character(i), as.character(i)), c("Cancer", "High grade"), c(0,0))
    colnames(tbbar0) <- c("context", "group", "freq")
    tbbar <- rbind(tbbar, tbbar0)
  }
  tbbar$number <- c(1:dim(tbbar)[1])
  tbbar$context = factor(tbbar$context, levels=unique(tbbar$context))
  
  p <- ggplot(data = tbbar, aes(x = number, y = freq, fill = group)) + 
    geom_bar(stat = "identity", size = 1) +
    scale_fill_manual(values=c("#085A9C","#943126"))+
    xlab("")+ylab("")+
    theme_prism(palette = "colors",
                base_fontface = "plain", # 字体样式，可??? bold, plain, italic
                base_family = "serif", # 字体格式，可??? serif, sans, mono, Arial???
                base_size = 15,  # 图形的字体大???
                base_line_size = 0.8)+ # 坐标轴的粗细
    theme(panel.spacing.x = unit(0, "cm"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.ticks.y = element_blank(),
          legend.position = c(0.9,0.7),
          legend.box.background = element_rect(fill = NA, color= "black", linetype = 1),
          panel.background = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          strip.background.x = element_blank(),
          strip.background.y = element_blank(),
          strip.text.x = element_blank())+
    scale_x_discrete(labels = tbbar$context)+#?ı?X??
    geom_rect(aes(xmin = 2.5, xmax = 34.5,
                  ymin = -30, ymax = -10),
              fill = "#5F6A6A")+
    geom_rect(aes(xmin = 36.5, xmax = 68.5,
                  ymin = -30, ymax = -10),
              fill = "#5F6A6A")+
    geom_rect(aes(xmin = 70.5, xmax = 102.5,
                  ymin = -30, ymax = -10),
              fill = "#5F6A6A")+
    geom_rect(aes(xmin = 104.5, xmax = 136.5,
                  ymin = -30, ymax = -10),
              fill = "#5F6A6A")+
    geom_rect(aes(xmin = 138.5, xmax = 170.5,
                  ymin = -30, ymax = -10),
              fill = "#5F6A6A")+
    geom_rect(aes(xmin = 172.5, xmax = 204.5,
                  ymin = -30, ymax = -10),
              fill = "#5F6A6A")+
    geom_rect(aes(xmin = 70.5, xmax = 136.5,
                  ymin = -190, ymax = -170),
              fill = "#5F6A6A")+
    annotate("text", x = 103.5, y= -250, label = paste0(unique(tb$context),sep="",collapse="\n"), 
             angle =90, size = size) +
    geom_segment(aes(x = 18.5, y = -30, xend = 103.5, yend = -170), colour = "#5F6A6A", lty=2) +
    geom_segment(aes(x = 52.5, y = -30, xend = 103.5, yend = -170), colour = "#5F6A6A", lty=2) +
    geom_segment(aes(x = 86.5, y = -30, xend = 103.5, yend = -170), colour = "#5F6A6A", lty=2) +
    geom_segment(aes(x = 120.5, y = -30, xend = 103.5, yend = -170), colour = "#5F6A6A", lty=2)+
    geom_segment(aes(x = 154.5, y = -30, xend = 103.5, yend = -170), colour = "#5F6A6A", lty=2)+
    geom_segment(aes(x = 186.5, y = -30, xend = 103.5, yend = -170), colour = "#5F6A6A", lty=2)+
    geom_rect(aes(xmin = 2.5, xmax = 34.5,
                  ymin = 2000, ymax = 2150),
              fill = colors[1])+
    geom_rect(aes(xmin = 36.5, xmax = 68.5,
                  ymin = 2000, ymax = 2150),
              fill = colors[2])+
    geom_rect(aes(xmin = 70.5, xmax = 102.5,
                  ymin = 2000, ymax = 2150),
              fill = colors[3])+
    geom_rect(aes(xmin = 104.5, xmax = 136.5,
                  ymin = 2000, ymax = 2150),
              fill = colors[4])+
    geom_rect(aes(xmin = 138.5, xmax = 170.5,
                  ymin = 2000, ymax = 2150),
              fill = colors[5])+
    geom_rect(aes(xmin = 172.5, xmax = 204.5,
                  ymin = 2000, ymax = 2150),
              fill = colors[6])+
    annotate("text", x = 18.5, y= 2075, label = TYPE[1], color= "white", size = 5)+
    annotate("text", x = 52.5, y= 2075, label = TYPE[2], color= "white", size = 5)+
    annotate("text", x = 86.5, y= 2075, label = TYPE[3], color= "white", size = 5)+
    annotate("text", x = 120.5, y= 2075, label = TYPE[4], color= "white", size = 5)+
    annotate("text", x = 154.5, y= 2075, label = TYPE[5], color= "white", size = 5)+
    annotate("text", x = 188.5, y= 2075, label = TYPE[6], color= "white", size = 5)+
    geom_rect(aes(xmin = -0.5, xmax = 0.4,ymin = 0, ymax = 2000),fill = "black")

  
  return(p)
}
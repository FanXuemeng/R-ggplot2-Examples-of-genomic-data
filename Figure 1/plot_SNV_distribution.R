plot_SNV_distribution <- function(type_occurrences ){
  
  ### Figure 1C
  group <- c()
  type <- c()
  value <- c()
  for (i in 1:length(sample_names)){
    s <- sample_names[i]
    group <- c(group, substr(s, nchar(s),nchar(s)))
  }
  group[group == "C"] <- "Cancer"
  group[group == "H"] <- "High grade"
  group[group == "L"] <- "Low grade"
  
  
  
  Valmean <- matrix(0, length(GROUP), length(TYPE))
  for(i in 1:length(GROUP)){
    al <- rowSums(type_occurrences[group == GROUP[i],])
    for (j in 1:length(TYPE)) {
      Valmean[i,j] <- mean(type_occurrences[group == GROUP[i], j]/al)
    }
  }
  rownames(Valmean) <- GROUP
  colnames(Valmean) <- TYPE
  
  SE <- matrix(0, length(GROUP), length(TYPE))
  for(i in 1:length(GROUP)){
    al <- rowSums(type_occurrences[group == GROUP[i],])
    for (j in 1:length(TYPE)) {
      SE[i,j] <- sd(type_occurrences[group == GROUP[i], j]/al)
    }
  }
  rownames(SE) <- GROUP
  colnames(SE) <- TYPE
  
  type <- rep(TYPE, length(GROUP))
  group <- c()
  for (i in 1:length(GROUP)) {
    group <- c(group, rep(GROUP[i], length(TYPE)))
  }
  se <- as.vector(SE)
  valmean <- as.vector(Valmean)
  data <- data.frame(group, type, valmean,se)
  
  
  p <- ggplot(data,aes(group,valmean, fill=group))+
    facet_grid(~type,scales = 'free',switch ="x")+
    geom_bar(stat="identity",position="dodge")+ #绘制柱状???
    geom_errorbar(aes(ymin = valmean - se, ymax = valmean + se), width = 0.2)+
    scale_fill_manual(values=c("#085A9C","#943126","#B9770E"))+
    labs(x="Samples",y=NULL)+#标题
    theme_prism(palette = "colors",
                base_fontface = "plain", # 字体样式，可??? bold, plain, italic
                base_family = "serif", # 字体格式，可??? serif, sans, mono, Arial???
                base_size = 20,  # 图形的字体大???
                base_line_size = 0.8)+ # 坐标轴的粗细
    xlab("")+ylab("Percent of mutation distributions")+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = c(0.9,0.9),
          legend.box.background = element_rect(fill = NA, color= "black", linetype = 1),
          panel.background = element_blank(),
          axis.line.x = element_blank(),
          strip.background.x = element_blank())+
    geom_hline(aes(yintercept = 0), size = 0.8)
  p
  return(p)
}

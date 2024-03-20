#!/usr/bin/env Rscript
library(readr)
library(ggplot2)


ChrLength <- c(
248956422,242193529,
198295559,190214555,181538259,170805979,159345973,145138636,138394717,
133797422,135086622,133275309,114364328,107043718,101991189,
90338345,83257441,80373285,58617616,64444167,46709983,50818468,
156040895,57227415)
ChrLengthCS <- C(0, cumsum(ChrLength))
chrname <- 1:24
names(chrname) <- c(1:22,'X','Y')

png(filename = paste("total_bam_ratio_", sample,".png",sep = ""), width = 1180, height = 1180,
    units = "px", pointsize = 20, bg = "white", res = NA)
plot(1:10)

ratio <- as.data.frame(read.csv("total_bam_ratio.csv", header=TRUE))
ratio <- ratio[ratio$X.ID == sample, ]
ploidy <- type.convert(c("2"))
sample_names <- ratio$X.ID
sample <- c()
group <- c()
unisam <- unique(sample_names)
for (i in 1:length(unisam)){
  print(i)
  s <- unisam[i]
  sample <- c(sample, rep(substr(s, 1,nchar(s)-1), sum(sample_names == s)))
  group <- c(group, rep(substr(s, nchar(s),nchar(s)), sum(sample_names == s)))
}
group[group == "C"] <- "Cancer"
group[group == "H"] <- "High grade"
group[group == "L"] <- "Low grade"


ratio <- cbind(ratio, data.frame(sample, group))
s <- "C14"

ratiosample <- ratio[ratio$sample == s, ]
color <- rep(0, dim(ratiosample)[1])
shape <- rep(0, dim(ratiosample)[1])
size <- rep(0, dim(ratiosample)[1])
group <- ratiosample$group
unigroup <- unique(group)
start <- ratiosample$Start
CN <- log2(ratiosample$Ratio)

color[ratiosample$CopyNumber>ploidy] <- 1
color[ratiosample$CopyNumber<ploidy & ratiosample$CopyNumber!= -1] = 2
for (i in c(1:22,'X','Y')) {
  tt <- which(ratiosample$Chromosome==i)
  if (length(tt)>0) {
    if(chrname[i]>1) {start[tt] <- start[tt] + ChrLengthCS[chrname[i]]}
  }
}

# 添加指示染色体起始位置的点
color <- as.factor(c(color, rep(rep(3, 24), length(unigroup))))
shape <- as.factor(c(shape, rep(rep(1, 24), length(unigroup))))
size <- as.factor(c(size, rep(rep(1, 24), length(unigroup))))
start <- c(start, rep(ChrLengthCS[2:25], length(unigroup)))
for (i in 1:length(unigroup)) {
  group <- c(group, rep(unigroup[i], 24))
}
CN <- c(CN, rep(rep(ploidy, 24), length(unigroup)))

RatioPlot <- data.frame(start, CN, color, shape, size, group)

ggplot(RatioPlot, aes(x = start, y = CN, colour = color, shape = shape, size = size)) +
  geom_point() +
  facet_grid(group ~ ., scales = "free", space = "free") +
  geom_hline(aes(yintercept=ploidy)) +
  ylab("Copy Number (log2)") + xlab("") +
  scale_shape_manual(values = c(16, 15)) +
  scale_color_manual(values = c(colors()[c(88, 136, 461)], "black")) +
  scale_size_manual(values = c(0.5, 1)) + theme_classic()+ # 坐标轴的粗细
  theme(panel.spacing.x = unit(0, "cm"),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size=1),
        axis.line = element_line( color="black"),
        legend.position = 'none',
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(size=1),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, size = 15),
        strip.background = element_rect(color = "white", fill = "white"))
  

dev.off()





png(filename = paste("b",".png",sep = ""), width = 1180, height = 1180,
    units = "px", pointsize = 20, bg = "white", res = NA)

BAF <-as.data.frame(read.csv("total_bam_BAF.csv", header=TRUE))

sample_names <- BAF$ID
sample <- c()
group <- c()
unisam <- unique(sample_names)
for (i in 1:length(unisam)){
  print(i)
  s <- unisam[i]
  sample <- c(sample, rep(substr(s, 1,nchar(s)-1), sum(sample_names == s)))
  group <- c(group, rep(substr(s, nchar(s),nchar(s)), sum(sample_names == s)))
}
group[group == "C"] <- "Cancer"
group[group == "H"] <- "High grade"
group[group == "L"] <- "Low grade"


BAF <- cbind(BAF, data.frame(sample, group))
s <- "C14"

Position <- c()
value <- c()
color <- c()
size <- c()
group <- c()
for (i in c(1:22,'X','Y')) {
  tt <- which(BAF$Chromosome==i)
  if (length(tt)>0){
    lBAF <-BAF[tt,]
    # 
    # plot(lBAF$Position,lBAF$BAF,ylim = c(-0.1,1.1),xlab = paste ("position, chr",i),
    #      ylab = "BAF",pch = ".",col = colors()[1])
    tt <- which(lBAF$A==0.5)		
    Position <- c(Position, lBAF$Position[tt]+ChrLengthCS[chrname[i]])
    value <- c(value, lBAF$BAF[tt])
    color <- c(color, rep(1, length(tt)))
    size <- c(size, rep(1, length(tt)))
    group <- c(group, lBAF$group[tt])
    # points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[92])
    tt <- which(lBAF$A!=0.5 & lBAF$A>=0)
    Position <- c(Position, lBAF$Position[tt]+ChrLengthCS[chrname[i]])
    value <- c(value, lBAF$BAF[tt])
    color <- c(color, rep(2, length(tt)))
    size <- c(size, rep(1, length(tt)))
    group <- c(group, lBAF$group[tt])
    # points(lBAF$Position[tt],lBAF$BAF[tt],pch = ".",col = colors()[62])
    tt <- 1
    pres <- 1
    
    if (length(lBAF$A)>4) {
      for (j in c(2:(length(lBAF$A)-pres-1))) {
        if (lBAF$A[j]==lBAF$A[j+pres]) {	
          tt[length(tt)+1] <- j 
        }
      }
      Position <- c(Position, lBAF$Position[tt]+ChrLengthCS[chrname[i]])
      value <- c(value, lBAF$A[tt])
      color <- c(color, rep(3, length(tt)))
      size <- c(size, rep(2, length(tt)))
      group <- c(group, lBAF$group[tt])
      # points(lBAF$Position[tt],lBAF$A[tt],pch = ".",col = colors()[24],cex=4)
      Position <- c(Position, lBAF$Position[tt]+ChrLengthCS[chrname[i]])
      value <- c(value, lBAF$B[tt])
      color <- c(color, rep(3, length(tt)))
      size <- c(size, rep(2, length(tt)))
      group <- c(group, lBAF$group[tt])
      # points(lBAF$Position[tt],lBAF$B[tt],pch = ".",col = colors()[24],cex=4)	
    }
    
    tt <- 1
    pres <- 1
    if (length(lBAF$FittedA)>4) {
      for (j in c(2:(length(lBAF$FittedA)-pres-1))) {
        if (lBAF$FittedA[j]==lBAF$FittedA[j+pres]) {	
          tt[length(tt)+1] <- j 
        }
      }
      Position <- c(Position, lBAF$Position[tt]+ChrLengthCS[chrname[i]])
      value <- c(value, lBAF$FittedA[tt])
      color <- c(color, rep(4, length(tt)))
      size <- c(size, rep(2, length(tt)))
      group <- c(group, lBAF$group[tt])
      
      # points(lBAF$Position[tt],lBAF$FittedA[tt],pch = ".",col = colors()[463],cex=4)
      Position <- c(Position, lBAF$Position[tt]+ChrLengthCS[chrname[i]])
      value <- c(value, lBAF$FittedB[tt])
      color <- c(color, rep(4, length(tt)))
      size <- c(size, rep(2, length(tt)))
      group <- c(group, lBAF$group[tt])
      # points(lBAF$Position[tt],lBAF$FittedB[tt],pch = ".",col = colors()[463],cex=4)	
    }
  }
}



BAFPlot <- data.frame(Position, value, color, size, group)

ggplot(BAFPlot, aes(x = Position, y = value, colour = factor(color), size = factor(size))) +
  geom_point() +
  facet_grid(group ~ ., scales = "free", space = "free") +
  ylab("BAF") + xlab("") +
  ylim(c(-0.1, 1.1)) +
  scale_color_manual(values = colors()[c(92, 62, 24, 463)]) +
  scale_size_manual(values = c(0.5, 1)) + theme_classic()+ # 坐标轴的粗细
  theme(panel.spacing.x = unit(0, "cm"),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size=1),
        axis.line = element_line( color="black"),
        legend.position = 'none',
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(size=1),
        strip.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, size = 15),
        strip.background = element_rect(color = "white", fill = "white"))


dev.off()
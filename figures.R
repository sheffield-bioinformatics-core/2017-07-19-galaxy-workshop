library(ggplot2)
library(dplyr)
library(gridExtra)

A <- rnorm(15,mean = 5, sd = 1)
B <- rnorm(15, mean=3, sd=1)
my_pal <- c(rgb(0,159,218,maxColorValue = 255),
            rgb(31,20,93,maxColorValue = 255))
            
p1 <- data.frame(A,B) %>% 
  tidyr::gather() %>% 
  ggplot(aes(x=key,y=value,fill=key)) + geom_boxplot(alpha=0.7,outlier.shape = NA) + geom_jitter(width=0.1) + scale_fill_manual(name="Condition",values=my_pal) + xlab("Condition") + ylab("Expression level of gene X")
ggsave(p1, file="media/de_boxplot.png")
p2 <- data.frame(A,B) %>% 
  tidyr::gather() %>% 
  ggplot(aes(x=value,fill=key)) + geom_density(alpha=0.7) + scale_fill_manual(name="Condition",values=my_pal) + xlab("Condition") + ylab("Expression level of gene X")
grid.arrange(p1,p2)

ggsave(grid.arrange(p1,p2), file="media/de_explained.png")

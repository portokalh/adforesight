library(ggplot2)
myrmse<-read.csv('/Users/alex/AlexBadea_MyPapers/MEMRI/xlsshets/GOLDEN_RMSE_4R.csv')
output.file<-'/Users/alex/AlexBadea_MyPapers/MEMRI/xlsshets/RMSE_4folda'
myrmse<-myrmse[which(myrmse$full == "0"),]  


# p <- ggplot(myrmse, aes(x=model, y=RMSE,fill=as.factor(full), alpha=0.5),facet_grid(~contrast)) + 
#   #scale_color_manual(values=c("red", "green")+
#   geom_boxplot()+
#   facet_wrap(~ contrast) +
#   theme_bw() +
#   theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
#         axis.text.y = element_text(face="bold", size=14, angle=0),
#         axis.line.x = element_line(colour = 'black', size=0.5),
#         axis.line.y = element_line(colour = 'black', size=0.5),
#         # panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank())

#p <- ggplot(myrmse, aes(x=contrast, y=RMSE,fill=as.factor(~steps), alpha=0.5)) + 
  #scale_color_manual(values=c("red", "green")+
  p <- ggplot(myrmse, aes(x=contrast, y=RMSE,fill=as.factor(steps), alpha=0.75)) + 
  #scale_color_manual(values=c("red", "green")+
  geom_boxplot()+
  theme_bw() +
  theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
        axis.text.y = element_text(face="bold", size=14, angle=0),
        axis.line.x = element_line(colour = 'black', size=0.5),
        axis.line.y = element_line(colour = 'black', size=0.5),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.background = element_blank())

p
# Tourner le box plot
p + coord_flip()


ggsave(paste(output.file,'allinone.pdf',sep=''), plot = last_plot(), device = 'pdf',
       scale = 1, width = 10, height = 5, units = c("in"),dpi = 300)


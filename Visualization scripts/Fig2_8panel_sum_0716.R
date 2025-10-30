# 2025-0716
# GUO
# Fig2 a-h
fig2_df <- read.table("/Users/xr/Desktop/00_PhD/05_Q40/data/figures/Fig2_8panel_summary_0716.csv")
head(fig2_df)
unique(fig2_df$Metric)
rownames(fig2_df) <- 1:542
which.min(fig2_df[fig2_df$Group=="DNA"&fig2_df$Metric=="Mapping quality","Value"])
fig2_df <- fig2_df[-446,]
fig2_df$Metric <-factor(levels = c("Phred quality score","Duplication (%)",
                                   "Aligned reads","Mapping quality","bias"),fig2_df$Metric)

p <- ggplot(data = fig2_df ,       
            aes(x=Quality, y=Value, fill=Quality)) +
  #geom_half_violin(side="r")+
  #geom_half_boxplot(side = "r") +    
  # geom_half_point_panel(side = "l")+
  geom_half_violin(side = "r", color=NA, alpha=0.35) +    
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, 
                    outlier.size = 1,
                    outlier.stroke = 0.2,
                    outlier.shape = NA, 
                    width=0.2, linewidth=0.5) +    
  geom_half_point_panel(aes(fill = Quality),side = "l", 
                        range_scale = .85,
                        shape=21, size=1.5, color="white")+
  scale_fill_manual(values = c(pal_simpsons()(16)[c(2,8)])) +
  #scale_fill_simpsons()+
  theme_bw()+
  theme(
    axis.text.x = element_text(face = "bold",color = "black"),
    axis.text.y = element_text(face = "bold",color = "black"),
    axis.title  = element_text(face = "bold",color = "black"),
    panel.grid.minor  =  element_blank(),
    legend.background =  element_blank(),
    legend.title = element_text(face = "bold",color = "black"),
    legend.text = element_text(color = "black")
  )+
  scale_y_continuous(expand = expansion(mult=c(0.1,0.1)))+
  stat_compare_means(aes(group = Quality),
                     label="p.signif",
                     method="t.test",
                     show.legend = F)+
  facet_wrap(Group~Metric,scales = "free",ncol = 4);p
ggsave("/Users/xr/Desktop/Quartet_as/new/Fig6_rain_box_0711.png",p,width = 9,height = 5.5)
topptx(p,"/Users/xr/Desktop/Quartet_as/new/Fig6_rain_box_0711.pptx",width = 7.2,height = 4.4)
topptx(p,"/Users/xr/Desktop/00_PhD/05_Q40/plot/Fig2a-h_box_0716.pptx",width=7,height = 6)
??stat_compare_means

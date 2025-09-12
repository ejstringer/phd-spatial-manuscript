ind %>% 
  arrange(npp) %>% 
  mutate(species = ifelse(species == 'Pseudomys hermannsburgensis',
                          'P. hermannsburgensis',
                          'S. youngsoni'),
         phase = factor(phase, levels = c('low', 'increase', 'decrease'))) %>% 
  ggplot(aes(km, euclidean, colour = species,
             fill = npp, group = species))+
  #facet_wrap(~species, scale = 'free')+
  facet_grid(species~phaseNo, scale = 'free_y')+
  geom_point(pch = 21, colour = 'grey98', stroke = 0.25)+
  scale_shape_manual(values = c(21:24))+ 
  scale_fill_gradientn(colours = nppCol,
                       na.value = 'white',
                       name = "NPP")+
  # geom_smooth(method = 'loess', se = F, show.legend = F)+
  scale_color_manual(values = c('black', 'grey20'),
                     name = 'Comparisons')+
  theme_bw()+
  scale_x_continuous(breaks = c(0, 40, 80))+
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(1, units = 'cm'),
        legend.key.height = unit(0.35, units = 'cm'),
        strip.text.x = element_text(face = 'bold'),
        strip.text.y = element_blank(),# element_text(face = 'italic'),
        strip.background = element_blank())+  
  ylab('genetic distance')+
  xlab('geographical distance (km)')+
  # scale_x_log10()+
  theme(panel.grid = element_blank()) -> con

l1 <- g_legend(con)

lay <- rbind(
  #c(rep(1,10),2,2),
  #c(rep(1,10),2,2),
  c(rep(1,10),2,2),
  c(rep(1,10),2,2),
  c(rep(1,10),NA,NA),
  #c(rep(1,10),3,3),
  #c(rep(1,10),3,3),
  c(rep(1,10),3,3),
  c(rep(1,10),3,3),
  c(rep(c(1,1),each = 5),NA,NA),
  c(rep(c(1,1),each = 5),NA,NA))
lay <- rbind(
  #c(rep(1,8),2,2),
  #c(rep(1,8),2,2),
  c(rep(1,8),2,2),
  c(rep(1,8),2,2),
  c(rep(1,8),3,3),
  #c(rep(1,8),3,3),
  #c(rep(1,8),3,3),
  c(rep(1,8),3,3),
  c(rep(1,8),3,3),
  c(rep(c(1),each = 8),NA,NA),
  c(rep(c(1),each = 8),NA,NA))

fig1 <- grid.arrange(con + theme(strip.text.y = element_blank()),
                     pic, pic2,
                     layout_matrix = lay)

fig1 <- con
ggsave('./figures/IBD_raw.png', plot = fig1,
       height =10, width = 16, units = 'cm', dpi = 600)



# load --------------------------------------------------------------------

source('./code/libraries.R')
source('./code/functions_data-tables.R')
source('./code/functions_genetic-data-creation.R')

phaseCol <- c(
  increase = terrain.colors(10)[1],
  decrease = terrain.colors(10)[4],
  low = terrain.colors(10)[7],
  stable = 'grey70')

nppCol <- rev(terrain.colors(225)[c(seq(1,80,2),
                                    seq(81,225,10))])[-1]

fstdata <- read.csv('./output/fst_analysis.csv') %>% 
  mutate(phaseNo = factor(phaseNo, levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                   rep(1:3, each = 3))),
         phase = ifelse(grepl('young', species), 'stable', phase),
         ph = grepl('herm', species),
         phase = factor(phase, levels = c('increase', 'decrease','low', 'stable')))

ind <- read.csv('./output/individual_analysis.csv') 

ibdcorr <- read.csv('./output/ibd_analysis.csv') %>% 
  mutate(phase = case_when(
    phase == 'L' ~ 'low',
    phase == 'D' ~ 'decrease',
    phase == 'I' ~ 'increase'),
  phase = ifelse(grepl('young', species), 'stable', phase)) %>% 
  left_join(unique(ind[,c('phaseNo', 'species', 'captures')]))

# npp ~ captures ----------------------------------------------------------

## data -----------------
cn <- fstdata |>
  mutate(ph = (grepl('herm', species))) |> 
  dplyr::select(species,ph, phaseNo, npp.log,npp, phase, captures) |>
  unique() |>
  glimpse() %>% 
  arrange(species)
## correlation ----------
cn %>% 
  group_by(species) %>% 
  summarise(cor = cor(captures, npp.log)) 

## models ---------------

m_parameters_ph <- (lm(captures ~ npp.log, data = filter(cn, species == "Pseudomys hermannsburgensis")))
m_parameters_sy <-(lm(captures ~ npp.log, data = filter(cn, species == "Sminthopsis youngsoni")))

summary(m_parameters_ph)
summary(m_parameters_sy)

## plots ----------------
figcor <- ggplot(cn, aes(y = captures, x = npp, group = species)) + 
  stat_smooth(data = filter(cn, ph),
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  geom_point(aes(fill = phase), size = 5, pch = 21) +
  stat_smooth(method = lm, se = F, colour = 'grey30', lty = 1) +
  facet_wrap(~species, ncol = 1, scale = "free_x") +
  theme_bw()+
  scale_x_log10()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.78,0.3),
        legend.key.size = unit(0.5, units = 'cm'),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 12))+
  xlab(expression('Net primary productivity '[log-scale]))+
  scale_fill_manual(values = phaseCol,
                    labels = str_to_title(names(phaseCol)),
                    name = 'Population phase')+
  ylab('Mean captures'); figcor


# cor plot ----------------------------------------------------------------


ggsave('./figures/caps_npp_correlation.png', plot = figcor,
       units = "cm", width = 10, height = 14, dpi = 600)

# fst npp captures --------------------------------------------------------

# fst by captures
# mean fst
mfc <- fstdata |>
  group_by(species, phaseNo, phase, captures) |>
  summarise(mean.fst = mean(fst)) |>
  glimpse()

figCaps <- ggplot(fstdata, aes(y = fst, x = captures, fill = phase, group = species)) +
  stat_smooth(data = filter(fstdata, ph),
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  geom_point(aes(colour = phase, fill = phase), size = 3, alpha = 0.2) +
  geom_point(data = mfc, aes(y = mean.fst), colour = 'black', 
             shape = 23, size = 4) +
  stat_smooth(method = lm, se = F, colour = 'grey30') +
  facet_wrap(~ species, ncol = 1, scale = "free_x") +
  theme_bw()+
  xlab('Mean captures')+
  ylab(expression(italic('F')[ST]))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',#c(0.85,0.73),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 12))+
  scale_fill_manual(values = phaseCol,
                    name = 'Population phase',
                    labels = str_to_title(names(phaseCol)))+
  scale_colour_manual(values = phaseCol,
                      labels = str_to_title(names(phaseCol)),
                      name = 'Population phase');figCaps


# fst by npp
# mean fst
mfn <- fstdata |>
  group_by(species, phaseNo, phase, npp.log, npp) |>
  summarise(mean.fst = mean(fst),
            x = npp.log) |>
  glimpse()

figNPP <- ggplot(dat, aes(y = fst, x = npp, fill = phase)) +
  stat_smooth(data = filter(dat, ph), show.legend = F,
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  geom_point(aes(colour = phase, fill = phase), size = 3, show.legend = F,
             alpha = 0.2) +
  geom_point(data = mfn, aes(y = mean.fst), 
             colour = 'black', shape = 23, size = 4) +
  stat_smooth(method = lm, se = F, colour = 'grey30', show.legend = F,
              aes(group = species)) +
  facet_wrap(~ species, ncol = 1, scale = "free_x") +
  scale_colour_discrete(na.translate = F)+
  theme_bw()+
  scale_x_log10()+
  xlab(expression('Net primary productivity '[log-scale]))+
  ylab(expression(italic('F')[ST]))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.83,0.28),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 12))+
  scale_fill_manual(values = phaseCol,
                    name = 'Population phase', na.translate = FALSE,
                    labels = str_to_title(names(phaseCol)))+
  scale_colour_manual(values = phaseCol,
                      labels = str_to_title(names(phaseCol)),
                      name = 'Population phase');figNPP


# FST arrange plot --------

figCapsNPP <- grid.arrange(figCaps, 
                           figNPP + theme(legend.position = c(0.7,0.3),
                                          legend.key.size = unit(0.5, 
                                                                 units = 'cm')),
                           ncol = 2)

ggsave(filename = './figures/caps_npp_fst_spp.png', figCapsNPP,
       units = "cm", width = 16, height = 14, dpi = 600)

# models -----------------

#https://stats.stackexchange.com/questions/570148/is-there-anything-like-a-two-way-anova-but-for-continuous-independent-variables
m_fst_anova_ph <- anova(lm(fst ~ captures + npp.log,
                           data = filter(datScaled,
                                         species == "Pseudomys hermannsburgensis")))

m_fst_anova_sy <- anova(lm(fst ~ captures + npp.log,
                           data = filter(datScaled,
                                         species == "Sminthopsis youngsoni")))

m_fst_anova_ph
m_fst_anova_sy

m_fst_ph <- (lmer(fst ~ captures + npp.log + (1|phaseNo),
            data = filter(datScaled,
                          species == "Pseudomys hermannsburgensis")))

m_fst_sy <- (lmer(fst ~ captures + npp.log + (1|phaseNo),
                           data = filter(datScaled,
                                         species == "Sminthopsis youngsoni")))

summary(m_fst_ph)
summary(m_fst_sy)


# individual --------------------------------------------------------------

plot( unique(ind[ind$species2 == 'ph',c('npp.log', 'captures', 'phaseNo')])$npp.log,
      predict(m_fst_ph, newdata = unique(ind[ind$species2 == 'ph',c('npp.log', 'captures', 'phaseNo')])))

# correlograms ------------------------------------------------------------

mantel.results <- ibdcorr %>% 
  mutate(
    sig = Pr.corrected. < 0.05 & Mantel.cor > 0,
    d.class = exp(breaksto)
    #sig = Pr.Mantel. < 0.05
  ) %>% 
  filter(complete.cases(Pr.corrected.))  

xibdcorr <- ibdcorr %>% 
  filter(class.index < 0) %>% 
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0, 
                             Mantel.cor)) 



## plot --------------------------------------------------------------------

figCorrelogram <- mantel.results %>%  
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0,
                             Mantel.cor)) %>%
  ggplot(aes(exp(breaksto), Mantel.cor, 
             fill= sig, colour = npp,
             group = phaseNo))+
  geom_hline(yintercept = 0, #size = 1,
             colour = 'grey', 
             linetype = 'dashed')+
  geom_line(linewidth = 0.75)+
  geom_point( shape = 21, size = 2, colour = 'black')+
  
  scale_x_log10(breaks = unique(exp(mantel.results$breaksto)))+
   scale_fill_manual(values = c('white', 'pink'),
                    labels = c('No IBD', 'IBD'),
                    name = NULL)+
  scale_colour_gradientn(colours = nppCol,
                         na.value = 'white',
                         name = "NPP")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = 'italic', size = 10),
        strip.text.y = element_text(size = 12))+
  xlab('Distance classes (km)')+
  ylab(expression(italic('r')))+
  facet_grid(~species, scale = 'free_x');figCorrelogram


fig_r_captures <- xibdcorr %>% 
  ggplot(aes(captures, Mantel.cor, fill = phase))+
  stat_smooth(data = filter(xibdcorr, grepl('herm', species)), show.legend = F,
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  scale_shape_manual(values = c(24,22))+
  geom_point(pch = 21, #aes(pch = sex_pairs), 
             size = 3,
             colour = 'black',
             #alpha = 0.5
  )+
  theme_bw()+
  guides(shape = 'none')+
  geom_smooth(method = 'lm', se = F, show.legend = F, aes(colour = species,
                                                          group = species),
              colour = 'grey30')+
  facet_grid(~species, scale = 'free') +
  ylab('mantel correlation (<1km)')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.89,0.65),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 11),
        legend.key.size = unit(0.5, units = 'cm'),
        strip.text.x = element_text(face = 'italic', size = 11),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank())+
  scale_fill_manual(values = phaseCol,
                    name = 'Population phase', na.translate = FALSE,
                    labels = str_to_title(names(phaseCol)))+
  scale_colour_manual(values = phaseCol,
                      labels = str_to_title(names(phaseCol)),
                      name = 'Population phase')+
  ylab(expression(italic("r")*' (distance class 0 - 1km)'))+
  xlab(expression('Mean captures')); fig_r_captures

fig_r_npp <- xibdcorr %>% 
  ggplot(aes(npp, Mantel.cor, fill = phase))+
  stat_smooth(data = filter(xibdcorr, grepl('herm', species)), show.legend = F,
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  scale_shape_manual(values = c(24,22))+
  geom_point(pch = 21, #aes(pch = sex_pairs), 
             size = 3,
             colour = 'black',
             #alpha = 0.5
  )+
  theme_bw()+
  guides(shape = 'none')+
  geom_smooth(method = 'lm', se = F, show.legend = F, aes(colour = species,
                                                          group = species),
              colour = 'grey30')+
  scale_x_continuous(trans = log_trans(), 
                     breaks = 8 * 2^seq(0,5,1))+
  facet_grid(~species) +
  ylab('mantel correlation (<1km)')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.885,0.63),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 11),
        legend.key.size = unit(0.5, units = 'cm'),
        strip.text.x = element_text(face = 'italic', size = 11),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank())+
  scale_fill_manual(values = phaseCol,
                    name = 'Population phase', na.translate = FALSE,
                    labels = str_to_title(names(phaseCol)))+
  scale_colour_manual(values = phaseCol,
                      labels = str_to_title(names(phaseCol)),
                      name = 'Population phase')+
  ylab(expression(italic("r")*' (distance class 0 - 1km)'))+
  xlab(expression('Net primary productivity '[log-scale])); fig_r_npp


# r arrange plot ----------------------------------------------------------

r_plot <- grid.arrange(figCorrelogram,
             fig_r_captures+theme(legend.position = 'none'),
             fig_r_npp + theme(legend.title = element_text(size = 10), 
                               legend.text = element_text(size = 9)))

ggsave(filename = './figures/caps_npp_r_spp.png',r_plot,
       units = "cm", width = 16, height = 19, dpi = 600)


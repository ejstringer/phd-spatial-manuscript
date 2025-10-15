em.Nm <- function(fst) (1/fst-1)/4

d<-fstdata %>% 
  group_by(species, phaseNo, phase, captures, npp.log, npp) %>% 
  summarise(fst = mean(fst),
            Nm = em.Nm(fst)) %>%
  mutate(mSrgt = Nm/(captures)) %>% 
 # filter(phase != 'stable') %>% 
  arrange(phase)
d

d %>% 
  ggplot(npp, Nm)+
  geom_point()+
  facet_grid(~species)+
  geom_smooth(method= 'lm')+
  theme_classic()

fstdata2<- fstdata %>% left_join(data.frame(species = d$species,
                                            phaseNo = d$phaseNo,
                                            mSrgt = (d$mSrgt))) %>% 
  filter(!is.infinite(mSrgt))

meandata <- fstdata2 |>
  group_by(species, phaseNo, phase, npp.log, captures, npp, mSrgt) |>
  summarise(mean.fst = mean(fst))

fstdata2 %>% 
  mutate(logmSrgt = log(mSrgt)) %>% 
  pivot_longer(cols = c(logmSrgt, captures)) %>% 
  ggplot(aes(value, fst,  group = phase))+
  geom_point(aes(colour = phase))+
  facet_wrap(species~name, scale = 'free')+
#scale_x_log10()+
  theme_bw()+
  geom_smooth(method = 'lm',aes(group = species), colour = 'black', se = F)+
  geom_smooth(method = 'lm', colour = 'grey', se = F)



# big fig Nm N m---------------
freee <- 'free'

nmplot <- ggplot(fstdata, aes(y = fst, x = npp, fill = phase)) +
  stat_smooth(data = filter(fstdata, ph), show.legend = F,
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  geom_point(aes(colour = phase, fill = phase), size = 3, show.legend = F,
             alpha = 0.2) +
  geom_point(data = mfn, aes(y = mean.fst), 
             colour = 'black', shape = 23, size = 4) +
  stat_smooth(method = lm, se = F, colour = 'grey30', show.legend = F,
              aes(group = species)) +
  facet_wrap(~ species, ncol = 2, scale = freee) +
  scale_colour_discrete(na.translate = F)+
  theme_bw()+
  scale_x_log10()+
  xlab(expression('Net primary productivity '[log-scale]))+
  ylab(expression(italic('F')[ST]))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position ='none',
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 12))+
  scale_fill_manual(values = phaseCol,
                    name = 'Population phase', na.translate = FALSE,
                    labels = str_to_title(names(phaseCol)))+
  scale_colour_manual(values = phaseCol,
                      labels = str_to_title(names(phaseCol)),
                      name = 'Population phase')
nmplot
mplot <-fstdata2 %>% 
  ggplot(aes(mSrgt, fst, colour = phase, fill = phase, group = phase))+
  scale_x_log10()+
  theme_bw()+
  facet_wrap(~species,scale = freee, ncol = 2)+
  stat_smooth(data = filter(fstdata2, grepl('hermann', species)), 
              show.legend = F,
              method = lm, se = F, aes(group = phase),
              lwd = 2)+
  geom_point(alpha = 0.2)+
  geom_point(data = meandata, aes(y = mean.fst), 
             colour = 'black', shape = 23, size = 3)+
  stat_smooth(method = lm, se = F, colour = 'black', show.legend = F,
              aes(group = species)) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.87,0.73),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 12))+
  scale_colour_manual(values = phaseCol,
                      labels = str_to_title(names(phaseCol)),
                      name = 'Population phase')+
  scale_fill_manual(values = phaseCol,
                    name = 'Population phase', na.translate = FALSE,
                    labels = str_to_title(names(phaseCol)))+
  ylab(expression(italic('F')[ST]))+
  xlab(expression('Migration rate '[(Nm / captures )~ log-scale]))
mplot
nplot <- fstdata2 %>% 
  ggplot(aes(captures, fst, colour = phase, fill = phase, group = phase))+
  facet_wrap(~species,scale = freee, ncol = 2)+
  theme_bw()+
  geom_smooth(data = filter(fstdata2, grepl('hermann', species)),
              method = 'lm', se = F, linewidth = 2)+
  geom_point(alpha = 0.2)+
  geom_point(data = meandata, aes(y = mean.fst), 
             colour = 'black', shape = 23, size = 3)+
  geom_smooth(method = 'lm', se = F, aes(group = species), colour = 'black')+
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 12))+
  scale_colour_manual(values = phaseCol,
                      labels = str_to_title(names(phaseCol)),
                      name = 'Population phase')+
  scale_fill_manual(values = phaseCol,
                    name = 'Population phase', na.translate = FALSE,
                    labels = str_to_title(names(phaseCol)))+
  ylab(expression(italic('F')[ST]))+
  xlab(expression('Mean captures '[100 /trap-nights]))

mplot
fstmn <- grid.arrange(nplot, 
                      mplot+theme(legend.position = c(0.7,0.3),
                                  legend.key.size = unit(0.5, units = 'cm')), 
                                  ncol = 2)

textNm <- grobTree(textGrob('Nm',
                                 x=0.01,  
                                 y= 0.5,
                                 hjust=0, rot = 0,     #was fontsize = 12
                                 gp=gpar(col="black", fontsize=18, fontface="italic")))
ggplot(fstdata)+
  textNm


textN <- grobTree(textGrob('N',
                            x=0.01,  
                            y= 0.5,
                            hjust=0, rot = 0,     #was fontsize = 12
                            gp=gpar(col="black", fontsize=18, fontface="italic")))

# 
# library(extrafont) 
# #font_import()
# loadfonts(device = "win")

textm <- grobTree(textGrob('m',
                            x=0.01,  
                            y= 0.5,
                            hjust=0, rot = 0,     #was fontsize = 12
                            gp=gpar(col="black", fontsize=18, fontface="italic", fontfamily = "Zapfino")))


repx <- 12
layout <- rbind(c(rep(1,repx), 4),
                c(rep(2,repx), 5),
                c(rep(3,repx), 6))
fstmn <- grid.arrange(nmplot,nplot, 
                      mplot+theme(legend.position = c(0.86,0.67),
                                  legend.key.size = unit(0.45, units = 'cm')),
                      textNm,textN, textm,
                      layout_matrix = layout)


ggsave(filename = './figures/caps_m_fst_spp.png', plot = fstmn,
       units = "cm", width = 16, height = 21, dpi = 600)

fstdata2$x <- fstdata2$captures
mn <-  lmer(fst ~ x * phase + (1|phaseNo), data = filter(fstdata2,
                                                         species2 == 'ph')) 
fstdata2$x <- log(fstdata2$mSrgt)
mm <- lmer(fst ~ x * phase + (1|phaseNo), data = filter(fstdata2,
                                                           species2 == 'ph')) 
fstdata2$x <- fstdata2$npp.log
mnm <- lmer(fst ~ x * phase + (1|phaseNo), data = filter(fstdata2,
                                                               species2 == 'ph')) 
                                                                
sjPlot::tab_model(mm, mn,mnm, dv.labels = c('m', 'N', 'Nm'))

fstdata2$x <- fstdata2$captures
mn <-  lmer(fst ~ x  + (1|phaseNo), data = filter(fstdata2,
                                                         species2 == 'ph')) 
fstdata2$x <- log(fstdata2$mSrgt)
mm <- lmer(fst ~ x  + (1|phaseNo), data = filter(fstdata2,
                                                        species2 == 'ph')) 
fstdata2$x <- fstdata2$npp.log
mnm <- lmer(fst ~ x  + (1|phaseNo), data = filter(fstdata2,
                                                         species2 == 'ph')) 


sjPlot::tab_model(mm, mn,mnm, dv.labels = c('m', 'N', 'Nm'))

MuMIn::r.squaredGLMM(mn)
MuMIn::r.squaredGLMM(mm)

summary(mn)
summary(mm)


mn <-  lmer(fst ~ captures +  (1|phaseNo), data = filter(fstdata2,
                                                                species2 == 'sy')) 

mm <- lmer(fst ~ log(mSrgt) + (1|phaseNo), data = filter(fstdata2,
                                                                 species2 == 'sy')) 


MuMIn::r.squaredGLMM(mn)
MuMIn::r.squaredGLMM(mm)

summary(mn)
summary(mm)




d$mscale <- scale(d$mSrgt, scale = T, center = F)
d$cscale <- scale(d$captures, scale = T, center = F) 

d %>% 
  pivot_longer(cols = c(mSrgt, captures)) %>% 
  ggplot(aes(npp, value, colour = phase, group = species))+
  geom_point()+
  facet_wrap(~name, scale = 'free')+
  scale_x_log10()+
  theme_classic()+
  geom_smooth(method = 'lm', se = T)


lm(npp ~ cscale*mscale, data = d) %>% 
  summary

lm(mSrgt ~ npp.log, data = d) %>% 
  summary

lm(captures ~ npp, data = d) %>% 
  summary




d %>% 
  ggplot(aes(captures, Nm, colour = phase, group = species))+
  geom_point()+
  #scale_y_log10()+
  theme_classic()+
  geom_smooth(method = 'lm', se = T)

plot(d$npp.log, d$fst)
plot(d$captures, d$fst)
plot(log(d$Nm), log(d$fst))
plot(d$npp.log,log(d$Nm))
plot(d$captures,d$nmlog)

npp.model <- lm(Nm~npp, d)
cap.model <- lm(Nm ~ captures, d)

m <- d$Nm/(d$captures*100)
npp.model %>% summary
summary(cap.model)

dp <- d %>%ungroup() %>%  
  mutate(nm.pred = predict(cap.model, newdata = .),
         residuals = Nm-nm.pred)
dp %>% ggplot(aes(captures, nm.pred))+
  geom_smooth(colour = 'grey90')+
  geom_point()+
  geom_point(aes(y = Nm), colour = 'grey50',size = 2)+
  theme_bw()+
  geom_segment(aes(xend = captures, yend = Nm, color = "resid"), colour = 'darkred')


lm(Nm ~ npp, d) %>% summary








# random variables --------------------------------------------------------
dp$residuals

reps <- 1000
library(lhs)

lower <- 0.001  
upper <- 0.999


sample_data <- randomLHS(reps, 9)
colnames(sample_data) <- d$phaseNo
migration <- (lower + sample_data * (upper - lower))

par(mfrow = c(3,3))
apply(migration, 2, hist)
par(mfrow = c(1,1))

population <- sapply(1:9, function(x) dp$Nm[x]/migration[,x])

x <- 500
dp$Nm/migration[x,]==population[x,]

best <- as.numeric(reps)
for (i in 1:reps) {
  
  N <- population[i,]
  rep.model <- lm(dp$nmlog~N)
  nm.pred <- predict(rep.model, newdata = data.frame(N))
  residx <- dp$nmlog-nm.pred
  
  best[i] <- sum(abs(dp$residuals -residx))
  
}

b<- which.min(best)
best[b]
population[b,]*migration[b,]
dp$Nm

plot(population[b,], dp$Nm)
dp %>% left_join(data.frame(phaseNo = dp$phaseNo, Nsim = population[b,],
                            m = migration[b,])) %>% 
  ggplot(aes(Nsim, Nm))+
  geom_point()


# ibd ---------------------------------------------------------------------


# raw ibd 5km ----------------
triangle <- data.frame(x = c(1,1,5),
                       y = c(88,60, 88, 42,28,42), sex_pairs = 'm-m',
                       npp = 100,species = rep(c('Pseudomys hermannsburgensis',
                                                 'Sminthopsis youngsoni'), each = 3))
trianglec <- data.frame(x = c(1,5),
                        y = c(88, 88,42,42), sex_pairs = 'm-m',
                        npp = 100,species = rep(c('Pseudomys hermannsburgensis',
                                                  'Sminthopsis youngsoni'), each = 2))
ggplot(ind %>% filter(km < 5.1) %>% 
         arrange(npp),
       aes(km, euclidean,
           fill = npp, group = species))+
  facet_wrap(~species,ncol = 2, scale = 'free')+
  geom_point(aes(pch = sex_pairs), colour = 'grey90',
             size = 3)+
  geom_line(data = triangle, aes(x,y), colour = 'grey55',
            linewidth = 0.75, lty = 2, alpha = 0.8)+
  scale_shape_manual(values = c(24:21),
                     name = 'Comparison')+ 
  geom_line(data = trianglec, aes(x,y), colour = 'grey55',
            linewidth = 0.75, lty = 2, alpha = 0.8)+
  guides(shape = guide_legend( 
    override.aes=list(colour = 'black')))+
  scale_fill_gradientn(colours = nppCol,
                       na.value = 'white',
                       name = "NPP")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = 'italic', size = 10.5))+  
  ylab('Genetic distance')+
  xlab('Geographical distance (km)')+
  theme(panel.grid = element_blank()) -> rawslope

ggplot(ind ,aes(km, euclidean, fill = npp))+
  facet_wrap(~species,ncol = 2, scale = 'free')+
  geom_point(aes(pch = sex_pairs), colour = 'grey90',
             size = 1.5, alpha = 0.9, pch = 21)+
  guides(shape = guide_legend( 
    override.aes=list(colour = 'black')))+
  scale_fill_gradientn(colours = nppCol,
                       na.value = 'white',
                       name = "NPP")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        strip.text.x = element_text(face = 'italic', size = 10.5))+  
  ylab('Genetic distance')+
  xlab('Geographical distance (km)')+
  theme(panel.grid = element_blank()) -> rawslope80


#ggsave(units = 'cm', './figures/IBD/fig_ibd_80km', height = 8, width = 16, dpi = 300)

tiff(filename = './figures/fig_ibd_80km.tiff', 
     units = 'cm', height = 16, width = 16, res = 600)
grid.arrange(rawslope80, rawslope, ncol = 1)
dev.off()



# correlograms ------------------------------------------------------------
ibdcorr %>%  
  mutate(
    sig = Pr.corrected. < 0.05 & Mantel.cor > 0,
    d.class = exp(breaksto)
    #sig = Pr.Mantel. < 0.05
  ) %>% 
  filter(complete.cases(Pr.corrected.)) -> mantel.results2 

 #ibdcorrsex%>% 
ibdcorrsexes_join %>% 
mutate(
    sig = Pr.corrected. < 0.05 & Mantel.cor > 0,
    #sig = Pr.Mantel. < 0.05
  ) %>% 
  filter(complete.cases(Pr.corrected.)) -> mantel.results 

# figure
xibdcorrexes <- do.call('rbind', ibdcorrsexes) %>% 
  filter(class.index < 0) %>% 
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0, 
                             Mantel.cor),
         sex_pairs = ifelse(sex_pairs == 'f-f', 'Female', 'Male')) 

phaseCol2 <- phaseCol
names(phaseCol2) <- toupper(str_sub(names(phaseCol2), 1,1))

mantel.results %>% 
  bind_rows(mantel.results2) %>%
  mutate(sig = ifelse(sig, 1, 0),
         sig = ifelse(n.dist < 30, 2, sig),
         sig = factor(sig)) %>% 
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0, Mantel.cor),
         sex_pairs = case_when(
           sex_pairs == 'f-f' ~ 'Females',
           sex_pairs == 'm-m' ~ 'Males',
           sex_pairs == 'all' ~ 'All'
         )
  ) %>%
  ggplot(aes(exp(breaksto), Mantel.cor, 
             fill= sig, colour = npp,
             group = phaseNo))+
  geom_hline(yintercept = 0, #size = 1,
             colour = 'grey', 
             linetype = 'dashed')+
  geom_line(linewidth = 0.75)+
  geom_point( shape = 21, size = 2, colour = 'black')+
  
  scale_x_log10(breaks = unique(exp(mantel.results$breaksto)))+
  #guides(size = 'none')+
  # scale_colour_manual(values = unname(phaseCol),
  #                     labels = c('low', 'increase', 'decrease'))+
  scale_fill_manual(values = c('white', 'pink', 'grey'),
                    labels = c('No IBD', 'IBD', 'n < 10'),
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
  facet_grid(sex_pairs~species, scale = 'free_x') -> mantel.plot.sex
mantel.plot.sex



# mantel ------------------------------------------------------------------

ibdcorrsexes_join %>% 
  bind_rows(ibdcorr_join) %>% 
  filter(dist < 1) %>% 
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0,Mantel.cor)) %>% 
  mutate(sex_pairs = ifelse(grepl('all', sex_pairs), 'All', sex_pairs)) %>% 
  ggplot(aes(npp, Mantel.cor, colour = sex_pairs,
             fill = sex_pairs, group = sex_pairs))+
  scale_shape_manual(values = c(24,22))+
  geom_point(pch = 21, #aes(pch = sex_pairs), 
             size = 3,
             colour = 'black',
             alpha = 0.5
  )+
  theme_bw()+
  guides(shape = 'none')+
  geom_smooth(method = 'lm', se = F, show.legend = F)+
  # scale_x_continuous(trans = log_trans(), 
  #                    breaks = 8 * 2^seq(0,5,1))+
  scale_x_log10()+
  facet_grid(~species) +
  ylab('mantel correlation (<1km)')+
  scale_fill_manual(values = c('grey', 'red', 'black'), #unname(phaseCol),
                    name = 'Comparison')+
  scale_colour_manual(values = c('grey', 'red', 'black'),
                      name = 'Comparison',
                      #labels = c('females', 'males')
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.89,0.65),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        legend.key.size = unit(0.5, units = 'cm'),
        strip.text.x = element_text(face = 'italic', size = 11),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank())+
  ylab(expression(italic("r")*' (distance class 0 - 1km)'))+
  xlab(expression('Net primary productivity '[log-scale])) -> fig_r_npp_sexes2;fig_r_npp_sexes2


ggsave('./figures/fig_r_npp_sexes.png', plot = fig_r_npp_sexes2,
       units = "cm", width = 14, height = 9, dpi = 600)



# load --------------------------------------------------------------------

source('./code/libraries.R')
library(emmeans)
source('./code/functions_data-tables.R')
source('./code/functions_genetic-data-creation.R')

em.coef.models <-function(m, spp.ph = T, modelname = 'Interaction'){
  
  summary(m)$coefficients %>% 
    as.data.frame() %>% 
    mutate(Species = ifelse(spp.ph, 'P. hermanns',
                            'S. youngsoni'),
           t = qt(0.975, df),
           me = t*`Std. Error`,
           lower = round(Estimate - me, 3),
           upper = round(Estimate + me, 4),
           confidence = paste0(round(Estimate,3), ' (',
                               lower, ' - ', upper, ')'),
           Parameter = rownames(.),
           #Species = ifelse(duplicated(Species), NA, Species),
           sig = `Pr(>|t|)` < 0.05,
           R2 = round(MuMIn::r.squaredGLMM(m)[1], 3),
         #  R2 = ifelse(duplicated(R2), NA, R2),
           model = modelname,
           Model = model#ifelse(duplicated(model), NA, model)
           ) %>% 
    dplyr::select(Species, Model, Parameter, confidence, R2, sig)%>% 
    rename(`Estimate (95% CI)` = confidence)
}

em.coef.models2 <-function(m, spp.ph = T, modelname = 'Interaction'){
  
  df1 <- summary(m)$coefficients %>% 
    as.data.frame() 
  df1$df <- ifelse('df' %in% names(df1), df1$df, m$df.residual)
  df1 %>% 
    mutate(Species = ifelse(spp.ph, 'P. hermanns',
                            'S. youngsoni'),
           t = qt(0.975, df),
           me = t*`Std. Error`,
           lower = round(Estimate - me, 3),
           upper = round(Estimate + me, 3),
           confidence = paste0(round(Estimate,3), ' (',
                               lower, ' - ', upper, ')'),
           Parameter = rownames(.),
           #Species = ifelse(duplicated(Species), NA, Species),
           sig = `Pr(>|t|)` < 0.05,
           R2 = round(MuMIn::r.squaredGLMM(m)[2], 2),
           R2 = ifelse(duplicated(R2), NA, R2),
           model = modelname,
           Model = ifelse(duplicated(model), NA, model)) %>% 
    dplyr::select(Species, Model, Parameter, confidence, R2, sig)%>% 
    rename(`Estimate (95% CI)` = confidence)
}


phaseCol <- c(
  increase = terrain.colors(10)[1],
  decrease = terrain.colors(10)[4],
  low = terrain.colors(10)[7],
  stable = 'grey70')

nppCol <- rev(terrain.colors(225)[c(seq(1,80,2),
                                    seq(81,225,10))])[-1]

## fst ---------------
fstdata_NNm <- read.csv('./output/fst_analysis.csv') %>% 
  mutate(phaseNo = factor(phaseNo, levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                   rep(1:3, each = 3))),
         phase = ifelse(grepl('young', species), 'stable', phase),
         ph = grepl('herm', species),
         phase = factor(phase, levels = c('increase', 'decrease','low', 'stable')))
capadjust <- 25/1.69
em.Nm <- function(fst) (1/fst-1)/4
mfm <- fstdata_NNm |>
  group_by(species, phaseNo, phase, captures) |>
  summarise(mean.fst = mean(fst)) |>
  mutate(#mean.fst = ifelse(mean.fst == 0, 0.0001, mean.fst),
         Nm = em.Nm(mean.fst),
         mSrgt = Nm/(captures))|>
  glimpse()

mfn <- fstdata_NNm |>
  group_by(species, phaseNo, phase, npp.log, npp) |>
  summarise(mean.fst = mean(fst))

fstdata <- left_join(fstdata_NNm, mfm) %>%  filter(!is.infinite(mSrgt))

meandata <- fstdata |>
  group_by(species, phaseNo, phase, npp.log, captures, npp, mSrgt) |>
  summarise(mean.fst = mean(fst))

## ind ---------------
ind <- read.csv('./output/individual_analysis.csv') 

ibdcorr <- read.csv('./output/ibd_analysis_change_corrlog_function.csv') %>% 
  mutate(phase = case_when(
    phase == 'L' ~ 'low',
    phase == 'D' ~ 'decrease',
    phase == 'I' ~ 'increase'),
  phase = ifelse(grepl('young', species), 'stable', phase)) %>% 
  left_join(unique(ind[,c('phaseNo', 'species', 'captures')]))


ibdcorrsex <- read.csv('./output/ibd_sex_analysis.csv') %>% 
  mutate(phase = case_when(
    phase == 'L' ~ 'low',
    phase == 'D' ~ 'decrease',
    phase == 'I' ~ 'increase'),
    phase = ifelse(grepl('young', species), 'stable', phase)) %>% 
  left_join(unique(ind[,c('phaseNo', 'species', 'captures')]))

# FST ---------------------------------------------------------------------

## fst npp captures Nm --------------------------------------------------------

freee <- 'free'



### nm plot -----------------------------------------------------------------

nmplot <- ggplot(fstdata, aes(y = fst, x = npp, fill = phase)) +
  stat_smooth(data = filter(fstdata, ph), show.legend = F,
              method = lm, se = F, colour = 'grey75', aes(group = phase),
              lwd = 2) +
  geom_point(aes(colour = phase, fill = phase), size = 3, show.legend = F,
             alpha = 0.2) +
  geom_point(data = meandata, aes(y = mean.fst), 
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


### m plot ------------------------------------------------------------------


mplot <-fstdata %>% 
  ggplot(aes(mSrgt, fst, colour = phase, fill = phase, group = phase))+
  scale_x_log10()+
  theme_bw()+
  facet_wrap(~species,scale = freee, ncol = 2)+
  stat_smooth(data = filter(fstdata, grepl('hermann', species)), 
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
  xlab(expression(italic(Nm)~" / mean captures  "[ log-scale]))
mplot


### n plot ------------------------------------------------------------------


nplot <- fstdata %>% 
  ggplot(aes(captures, fst, colour = phase, fill = phase, group = phase))+
  facet_wrap(~species,scale = freee, ncol = 2)+
  theme_bw()+
  geom_smooth(data = filter(fstdata, grepl('hermann', species)),
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

nplot

### arrange plot --------

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

## models ---------

### sy ---------
fstdata$x <- fstdata$captures
mns <-  lmer(fst ~ x  + (1|phaseNo), data = filter(fstdata,
                                                  species2 == 'sy'),REML = F) 
fstdata$x <- log(fstdata$mSrgt)
mms <- lmer(fst ~ x  + (1|phaseNo), data = filter(fstdata, !is.infinite(mSrgt),
                                                 species2 == 'sy'),REML = F) 
fstdata$x <- fstdata$npp.log
mnms <- lmer(fst ~ x  + (1|phaseNo), data = filter(fstdata, 
                                                  species2 == 'sy'),REML = F) 

sjPlot::tab_model(mms, mns,mnms, dv.labels = c('m', 'N', 'Nm'), digits = 4)

aicsy<-AIC(mnms, mns, mms) %>% 
  as.data.frame %>% 
  mutate(delta = AIC - min(AIC),
         rel_lik = exp(-0.5 * delta),
         weights = (rel_lik / sum(rel_lik)),
         Model = c('Nm','N', 'm'),
         Species = 'S. youngsoni') %>% 
  arrange(AIC)
### no phase ---------
fstdata$x <- fstdata$captures
mn <-  lmer(fst ~ x  + (1|phaseNo), data = filter(fstdata,
                                                  species2 == 'ph'),REML = F) 
fstdata$x <- log(fstdata$mSrgt)
mm <- lmer(fst ~ x  + (1|phaseNo), data = filter(fstdata,
                                                 species2 == 'ph'),REML = F) 
fstdata$x <- fstdata$npp.log
mnm <- lmer(fst ~ x  + (1|phaseNo), data = filter(fstdata,
                                                  species2 == 'ph'),REML = F) 

sjPlot::tab_model(mm, mn,mnm, dv.labels = c('m', 'N', 'Nm'), digits = 4)




### phase -------
fstdata$x <- fstdata$captures
mnp <-  lmer(fst ~ x * phase + (1|phaseNo), data = filter(fstdata,
                                                         species2 == 'ph'),REML = F) 
fstdata$x <- log(fstdata$mSrgt)
mmp <- lmer(fst ~ x * phase + (1|phaseNo), data = filter(fstdata,
                                                        species2 == 'ph'),REML = F) 
fstdata$x <- fstdata$npp.log
mnmp <- lmer(fst ~ x * phase + (1|phaseNo), data = filter(fstdata,
                                                         species2 == 'ph'),REML = F) 

sjPlot::tab_model(mmp, mnp,mnmp, dv.labels = c('m', 'N', 'Nm'), digits = 3)
AIC(mnmp, mnp, mmp)

aicph<-AIC(mnm, mn, mm) %>% 
  as.data.frame %>% 
  mutate(delta = AIC - min(AIC),
         rel_lik = exp(-0.5 * delta),
         weights = (rel_lik / sum(rel_lik)),
         Model = c('Nm','N', 'm'),
         Species = 'P. hermanns') %>% 
  arrange(AIC)

aicphp<-AIC(mnmp, mnp, mmp) %>% 
  as.data.frame %>% 
  mutate(delta = AIC - min(AIC),
         rel_lik = exp(-0.5 * delta),
         weights = (rel_lik / sum(rel_lik)),
         Model = c('Nm_phase','N_phase', 'm_phase'),
         Species = 'P. hermanns') %>% 
  arrange(AIC)

aicall <- AIC(mnm, mn, mm,mnmp, mnp, mmp) %>% 
  as.data.frame %>% 
  mutate(delta = AIC - min(AIC),
         rel_lik = exp(-0.5 * delta),
         weights = (rel_lik / sum(rel_lik)),
         Model = c('Nm','N', 'm',
                   'Nm_phase','N_phase', 'm_phase'),
         Species = 'P. hermanns') %>% 
  arrange(AIC) %>% 
  rbind(aicsy) %>% 
  mutate(delta = round(delta, 1),
         weights = round(weights,2)) %>% 
  dplyr::select(-rel_lik)

# 
# aicall <- rbind(aicph, aicphp) %>% 
#   rbind(aicsy) %>% 
#   mutate(delta = round(delta, 1),
#          weights = round(weights,3)) %>% 
#   dplyr::select(-rel_lik)
aicall
MuMIn::r.squaredGLMM(mnp)
MuMIn::r.squaredGLMM(mmp)
### flextable ---------------------------------------------------------------
modeltbsx<-em.coef.models(mnm, T, 'Nm') %>% 
  rbind(em.coef.models(mn, T, 'N')) %>% 
  rbind(em.coef.models(mm, T, 'm')) %>%
  rbind(em.coef.models(mnmp, T, 'Nm_phase')) %>% 
  rbind(em.coef.models(mnp, T, 'N_phase')) %>% 
  rbind(em.coef.models(mmp, T, 'm_phase')) %>%
  rbind(em.coef.models(mnms, F, 'Nm')) %>% 
  rbind(em.coef.models(mns, F, 'N')) %>% 
  rbind(em.coef.models(mms, F, 'm')) %>%
  left_join(aicall[,3:6]) %>% 
  separate(Model, into = c('Model', 'interaction')) %>% 
  # mutate(Parameter = case_when(
  #   Model == 'Nm' & Parameter == 'x' ~ 'npp log',
  #   Model == 'N' & Parameter == 'x'~ 'Captures',
  #   Model == 'm' & Parameter == 'x'~ 'Nm/captures',
  #   .default = Parameter
  # )) %>% 
  mutate(R2 = ifelse(duplicated(R2), NA, R2),
         delta = ifelse(duplicated(paste(Species, Model,
                                         interaction)), NA, delta),
         weights = ifelse(duplicated(paste(Species, Model,
                                         interaction)), NA, weights),
         Model = ifelse(duplicated(paste(Species, Model,
                                         interaction)), NA, Model),
         
         Species = ifelse(duplicated(Species), NA, Species),
         Parameter = gsub('phase', '', Parameter),
         Parameter = gsub('\\(Intercept\\)', 'Intercept', Parameter)) %>% 
  dplyr::select(-interaction)
  


modeltbsx
mtb <- modeltbsx
modeltbsx %>% 
  dplyr::select(-sig) %>% 
  flextable() %>% 
  bold(part = 'header') %>%
  italic(j = 1) %>% 
  # italic(part = 'header', j = 2) %>% 
  border_remove() %>%  
  bold(j = 4, i = which(mtb$sig)) %>% 
  autofit() %>%
  italic(j = 2) %>% 
  #width(j = 3, width = 5, unit = 'cm') %>% 
  hline(i = c(nrow(mtb)),
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header',
        border = fp_border(color = "grey40", width = 2.5)) %>% 
  hline(i = which(!is.na(mtb$Model))[-1]-1,
        j = -1,
        border = fp_border(color = "grey80",
                           width = 1.5, 
                           style = 'solid')) %>% 
hline(i = which(mtb$Species == 'S. youngsoni')-1,
      border = fp_border(color = "grey60",
                         width = 2)) %>% 
  hline(i = which(!is.na(mtb$Model))[4]-1,
        j = -1,
        border = fp_border(color = "grey60",
                           width = 2)) %>% 
  compose(
    part = "header", j = 5,
    value = as_paragraph(('R'), as_sup('2'))) %>% 
  compose(
    part = "header", j = 6,
    value = as_paragraph(('AIC'), as_sub('delta'))) %>% 
  compose(
    part = "header", j = 7,
    value = as_paragraph(('AIC'), as_sub('weights'))) %>% 
  compose(
    part = "body", j = 3, i = grep('log', mtb$Parameter),
    value = as_paragraph(('log'), as_sub('e'), ('NPP'))) -> fxtb_fst;fxtb_fst


# individual --------------------------------------------------------------

## correlograms ------------------------------------------------------------
ibdcorr %>%  
  mutate(
    sig = Pr.corrected. < 0.05 & Mantel.cor > 0,
    d.class = exp(breaksto)
    #sig = Pr.Mantel. < 0.05
  ) %>% 
  filter(complete.cases(Pr.corrected.)) -> mantel.results2 

#ibdcorrsex%>% 
ibdcorrsex %>% 
  mutate(
    sig = Pr.corrected. < 0.05 & Mantel.cor > 0,
    #sig = Pr.Mantel. < 0.05
  ) %>% 
  filter(complete.cases(Pr.corrected.)) -> mantel.results 

# figure
xibdcorrexes <- ibdcorrsex %>% 
  filter(class.index < 0) %>% 
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0, 
                             Mantel.cor),
         sex_pairs = ifelse(sex_pairs == 'f-f', 'Female', 'Male')) 

phaseCol2 <- phaseCol
names(phaseCol2) <- toupper(str_sub(names(phaseCol2), 1,1))

mantel.results %>% 
  bind_rows(mantel.results2) %>%
  mutate(sig = ifelse(sig, 1, 0),
       #    sig = ifelse(n.dist < 10, 2, sig),
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
        strip.text.x = element_text(face = 'italic', size = 9),
        strip.text.y = element_text(size = 12))+
  xlab('Distance classes (km)')+
  ylab(expression(italic('r')))+
  facet_grid(sex_pairs~species, scale = 'free_x') -> mantel.plot.sex
mantel.plot.sex

# tiff('./figures/IBD/fig_correlogram_sexes.tiff', res = 600,
#      units = 'cm', width = 16, height = 18)
# mantel.plot.sex
# dev.off()

ggsave('./figures/fig_correlogram_sexes.png', plot = mantel.plot.sex,
       units = "cm", width = 14, height = 14, dpi = 300)

## mantel r ------------------------------------------------------------------

xibdcorr <- ibdcorrsex %>% 
  bind_rows(ibdcorr) %>% 
  filter(dist < 1) %>% 
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0,Mantel.cor)) %>% 
  mutate(sex_pairs = ifelse(grepl('all', sex_pairs), 'All', sex_pairs)) 

xibdcorr %>% 
  mutate(Mantel.cor = ifelse(Mantel.cor < 0, 0,Mantel.cor)) %>% 
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


### models ------------------------------------------------------------------

mod_r_ph <- lm(Mantel.cor ~ log(npp) , weights = (n.dist),
               data = filter(xibdcorr, 
                             grepl('herm', species),
                             sex_pairs == 'All')) 
mod_r_sy <- lm(Mantel.cor ~ log(npp), 
               data = filter(xibdcorr, 
                             grepl('young', species),
                             sex_pairs == 'All'))

mod_r_ph <- lm(Mantel.cor ~ log(npp) + sex_pairs,
               data = filter(xibdcorr, 
                             grepl('herm', species),
                             sex_pairs != 'All')) 
mod_r_sy <- lm(Mantel.cor ~ log(npp) + sex_pairs, 
               data = filter(xibdcorr, 
                             grepl('young', species),
                             sex_pairs != 'All'))

sjPlot::tab_model(mod_r_ph, mod_r_sy, digits = 3)

### flextables --------------------------------------------------------------

modeltbsx <- rbind(em.coef.models2(mod_r_ph, T,'P. hermann'),
                   em.coef.models2(mod_r_sy, F, 'S. youngsoni')) %>% 
  dplyr::select(-Model) %>% 
  mutate(Species = ifelse(duplicated(Species), NA, Species))

mtb <- modeltbsx
modeltbsx %>% 
  dplyr::select(-sig) %>% 
  flextable() %>% 
  bold(part = 'header') %>%
  italic(j = 1) %>% 
  # italic(part = 'header', j = 2) %>% 
  border_remove() %>%  
  bold(j = 3, i = which(mtb$sig)) %>% 
  autofit() %>%
  width(j = 3, width = 5, unit = 'cm') %>% 
  hline(i = c(nrow(mtb)),
        border = fp_border(color = "grey40", width = 3)) %>% 
  hline(part = 'header',
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline(i = which(mtb$Species == 'S. youngsoni')-1,
        border = fp_border(color = "grey60",
                           width = 2, 
                           style = 'dashed')) %>% 
  compose(
    part = "body", j = 2, i = grep('Intercept', mtb$Parameter),
    value = as_paragraph(('Intercept'), as_sub('mantel'))) %>% 
  compose(
    part = "body", j = 2, i = grep('cap', mtb$Parameter),
    value = as_paragraph(('Captures'))) %>% 
  compose(
    part = "header", j = 4,
    value = as_paragraph(('R'), as_sup('2'))) %>% 
  compose(
    part = "body", j = 2, i = grep('log', mtb$Parameter),
    value = as_paragraph(('log'), as_sub('e'), ('NPP'))) %>% 
  compose(
    part = "body", j = 2, i = grep('m-m', mtb$Parameter),
    value = as_paragraph(('Sex pair'), as_sub('m-m'), ('NPP')))-> fxtb_r_npp;fxtb_r_npp


## phase comparisons -------------------------------------------------------


### plot --------------------------------------------------------------------

ibdsexes <- ind %>% #rbind(ibdDataAll[[1]],ibdDataAll[[2]]) %>% 
  mutate(sites = factor(ifelse(km == 0, 'within', 
                               ifelse(km>1, 'beyond', 'between')),
                        levels = c('within', 'between', 'beyond')),
         ibd = factor(ifelse(sites == 'beyond', '> 1km', '< 1km'),
                      levels = c('< 1km', '> 1km')),
         phase = factor(phase, levels = c('low', 'increase', 'decrease'))) %>% 
  filter(sex_pairs %in% c('m-m',
                          # 'f-m', 
                          'f-f')) %>% filter(km <= 1)

ibdsexes %>% group_by(species) %>% 
  summarise(dist = mean(euclidean))

ibdplotting <- ibdsexes %>% 
  arrange((npp)) %>% 
  mutate(phase = factor(phase,
                        levels = c('low', 'increase', 'decrease', 'stable')),
         hline = ifelse(grepl('herm', species), 91, 39),
         sex_pairs = ifelse(sex_pairs == 'f-f', 'Females', 'Males'),
         species = ifelse(species == 'Pseudomys hermannsburgensis',
                          'P. hermannsburgensis',
                          'S. youngsoni'))  
(ibdplotting$phase %>% table %>% names)== names(phaseCol)
ggplot(ibdplotting,aes(km, euclidean, colour = phase,
                       fill = npp, group = phase))+
  facet_grid(species~sex_pairs, scale = 'free')+
  geom_point(pch = 21, colour = 'white', alpha = 0.7, fill = 'grey70')+
  scale_shape_manual(values = c(21:24))+ 
  scale_fill_gradientn(colours = nppCol,
                       na.value = 'white',
                       name = "NPP")+
  geom_hline(aes(yintercept = hline), alpha = 0.8,
             lty = 2, colour = 'grey70', lwd = 0.3)+
  geom_smooth(method = 'lm', se = F, lwd = 1)+
  scale_color_manual(values = phaseCol,
                     labels = c('Low','Increase', 'Decrease', 'Stable'),
                     name = 'Population phase')+
  #scale_color_manual(values = phaseCol)+
  theme_bw()+
  scale_x_continuous(limits = c(0,1), labels = c(seq(0,0.8, 0.25),'1'))+
  theme(panel.grid = element_blank(),
        #legend.position = 'right',
        legend.direction="vertical",
        #legend.box = "horizontal",
        legend.position = c(0.85, 0.14),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "grey"),
        legend.key.size = unit(0.35, units = 'cm'),
        legend.title = element_text(size = 10),
        strip.background = element_blank(),
        axis.title = element_text(size = 13),
        strip.text.y = element_text(face = 'italic', size = 11),
        strip.text.x = element_text(face =   'bold',size = 12))+  
  ylab('Genetic distance')+
  xlab('Geographical distance (km)') -> IBD_sexes_phase; IBD_sexes_phase

# 
# tiff(filename = './figures/IBD/fig_slope.tiff', 
#      units = 'cm', height = 12, width = 14, res = 600)
# IBD_sexes_phase
# dev.off()

ggsave('./figures/phasesexes_slope.png', plot = IBD_sexes_phase,
       units = 'cm', height = 12, width = 14, dpi = 300)



### models ------------------------------------------------------------------

phsex <- ind %>% filter(sex_pairs %in% c('f-f', 'm-m'),
                        species2== 'ph', km <= 1) %>% 
  mutate(phase = factor(phase, levels = c('low', 'increase', 'decrease')))

sysex <- ind %>% filter(sex_pairs %in% c('f-f', 'm-m'),
                        species2== 'sy', km <= 1)

names(phsex)

m <-lmer(euclidean ~ km * sex_pairs *  phase + (1|phaseNo), data = phsex)

m2 <-lmer(euclidean ~ km * sex_pairs + (1|phaseNo), data = sysex)

anovafx<-anova(m) %>% as.data.frame() %>% mutate_all(round, 3) %>% 
  mutate(parameter= rownames(.),
         species = 'P. hermann',
         `Pr(>F)` = ifelse(`Pr(>F)`<0.001, '<0.001', `Pr(>F)`)) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  relocate(species, parameter) %>% flextable() %>% 
  autofit() %>% 
  italic(j =1)

anovafxsy<-anova(m2) %>% as.data.frame() %>% mutate_all(round, 3) %>% 
  mutate(parameter= rownames(.),
         species = 'S. young',
         `Pr(>F)` = ifelse(`Pr(>F)`<0.001, '<0.001', `Pr(>F)`)) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  relocate(species, parameter) %>% flextable() %>% 
  autofit() %>% 
  italic(j =1)


#https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/


sydiff <- emmeans(m2, specs = pairwise ~ km:sex_pairs, at = list(km = 0))

differences <- emmeans(m, specs = pairwise ~ km:sex_pairs:phase, at = list(km = 0))

#### flextable ---------------------------------------------------------------


diffcontrasts <- as.data.frame(differences$contrasts) %>% 
  mutate(km.mean = 0,
         contrast = gsub('km0', '', contrast),
         contrast = gsub('\\(', '', contrast),
         contrast = gsub('\\)', '', contrast)) %>% 
  separate(contrast, into = c('con1', 'con2'),sep = ' - ') %>% 
  mutate(p.value = ifelse(p.value < 0.001, '<0.001',
                          sprintf("%.3f",round(p.value, 3))),
         df = '>3000') %>% 
  mutate_if(is.numeric, round, 2) %>% 
  dplyr::select(-df) %>% 
  mutate_all(trimws) %>% 
  separate(con1, into = c('sex1', 'phase1'), sep = ' ')%>% 
  separate(con2, into = c('sex2', 'phase2'), sep = ' ') %>% 
  mutate(sex = ifelse(sex1 == sex2, sex1, 'xdiff'),
         contrast = paste0('(', sex1, ' ', phase1, ')', ' - ',
                           '(', sex2, ' ', phase2, ')')) %>% 
  relocate(sex2, .after = sex1) %>% 
  arrange(sex, sex1) %>% 
  mutate(sex = ifelse(sex == 'xdiff', 'f-f:m-m', sex),
         estimate = as.numeric(estimate),
         z.ratio = as.numeric(z.ratio),
         contrast = ifelse(sex == 'f-f', paste0('(', phase1, ' - ', phase2,')'),
                           contrast),
         contrast = ifelse(sex == 'm-m', paste0('(', phase1, ' - ', phase2,')'),
                           contrast),
         contrast = ifelse(sex1 == 'm-m' & sex2 == 'f-f',
                           paste0('(', sex2, ' ', phase2, ')', ' - ',
                                  '(', sex1, ' ', phase1, ')'), contrast),
         estimate = ifelse(sex1 == 'm-m' & sex2 == 'f-f', -estimate, estimate),
         z.ratio = ifelse(sex1 == 'm-m' & sex2 == 'f-f', -z.ratio, z.ratio))%>% 
  dplyr::select(km.mean,sex,contrast, estimate, SE, z.ratio, p.value) %>% 
  rename(sex.comparison = sex)


fxtb_compair<-diffcontrasts %>% 
  mutate(km.mean = ifelse(duplicated(km.mean), NA, km.mean),
         sex.comparison = ifelse(duplicated(sex.comparison),
                                 NA, sex.comparison)) %>% 
  flextable() %>% 
  autofit() %>% 
  theme_alafoli() %>% 
  bold(j = 3:7, i = which(abs(diffcontrasts$z.ratio) >2.9)) %>% 
  hline(j = -1, i = c(3,6), 
        border= fp_border_default(width = 1.5, color = 'grey', style = 'dashed'))  



fxtb_compairsy<-(as.data.frame(sydiff$contrasts)) %>% 
  mutate(km.mean = 0,
         sex.comparison = '',
         contrast = gsub('km0', '', contrast),
         contrast = gsub('\\(', '', contrast),
         contrast = gsub('\\)', '', contrast),
         contrast = '(f-f) - (m-m)',
         p.value = ifelse(p.value < 0.001, '<0.001',
                          sprintf("%.3f",round(p.value, 3))),) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  relocate(km.mean, sex.comparison) %>% 
  dplyr::select(-df) %>% 
  flextable() %>% 
  autofit() %>% 
  theme_alafoli() %>% 
  bold(j = 3:7, i = which(abs(as.data.frame(sydiff$contrast)$t.ratio) > 2))
fxtb_compairsy


#### diff km = 1 -------------------------------------------------------------

sydiff <- emmeans(m2, specs = pairwise ~ km:sex_pairs, at = list(km = 1))

differences <- emmeans(m, specs = pairwise ~ km:sex_pairs:phase, at = list(km = 1))

#### flextable ---------------------------------------------------------------


diffcontrasts <- as.data.frame(differences$contrasts) %>% 
  mutate(km.mean = 1,
         contrast = gsub('km1', '', contrast),
         contrast = gsub('\\(', '', contrast),
         contrast = gsub('\\)', '', contrast)) %>% 
  separate(contrast, into = c('con1', 'con2'),sep = ' - ') %>% 
  mutate(p.value = ifelse(p.value < 0.001, '<0.001',
                          sprintf("%.3f",round(p.value, 3))),
         df = '>3000') %>% 
  mutate_if(is.numeric, round, 2) %>% 
  dplyr::select(-df) %>% 
  mutate_all(trimws) %>% 
  separate(con1, into = c('sex1', 'phase1'), sep = ' ')%>% 
  separate(con2, into = c('sex2', 'phase2'), sep = ' ') %>% 
  mutate(sex = ifelse(sex1 == sex2, sex1, 'xdiff'),
         contrast = paste0('(', sex1, ' ', phase1, ')', ' - ',
                           '(', sex2, ' ', phase2, ')')) %>% 
  relocate(sex2, .after = sex1) %>% 
  arrange(sex, sex1) %>% 
  mutate(sex = ifelse(sex == 'xdiff', 'f-f:m-m', sex),
         estimate = as.numeric(estimate),
         z.ratio = as.numeric(z.ratio),
         contrast = ifelse(sex == 'f-f', paste0('(', phase1, ' - ', phase2,')'),
                           contrast),
         contrast = ifelse(sex == 'm-m', paste0('(', phase1, ' - ', phase2,')'),
                           contrast),
         contrast = ifelse(sex1 == 'm-m' & sex2 == 'f-f',
                           paste0('(', sex2, ' ', phase2, ')', ' - ',
                                  '(', sex1, ' ', phase1, ')'), contrast),
         estimate = ifelse(sex1 == 'm-m' & sex2 == 'f-f', -estimate, estimate),
         z.ratio = ifelse(sex1 == 'm-m' & sex2 == 'f-f', -z.ratio, z.ratio))%>% 
  dplyr::select(km.mean,sex,contrast, estimate, SE, z.ratio, p.value) %>% 
  rename(sex.comparison = sex)


fxtb_compair1<-diffcontrasts %>% 
  mutate(km.mean = ifelse(duplicated(km.mean), NA, km.mean),
         sex.comparison = ifelse(duplicated(sex.comparison),
                                 NA, sex.comparison)) %>% 
  flextable() %>% 
  autofit() %>% 
  theme_alafoli() %>% 
  bold(j = 3:7, i = which(abs(diffcontrasts$z.ratio) >2.9)) %>% 
  hline(j = -1, i = c(3,6), 
        border= fp_border_default(width = 1.5, color = 'grey', style = 'dashed'))  



fxtb_compairsy1<-(as.data.frame(sydiff$contrasts)) %>% 
  mutate(km.mean = 0,
         sex.comparison = '',
         contrast = gsub('km0', '', contrast),
         contrast = gsub('\\(', '', contrast),
         contrast = gsub('\\)', '', contrast),
         contrast = '(f-f) - (m-m)',
         p.value = ifelse(p.value < 0.001, '<0.001',
                          sprintf("%.3f",round(p.value, 3))),) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  relocate(km.mean, sex.comparison) %>% 
  dplyr::select(-df) %>% 
  flextable() %>% 
  autofit() %>% 
  theme_alafoli() %>% 
  bold(j = 3:7, i = which(abs(as.data.frame(sydiff$contrast)$t.ratio) >2))
 
fxtb_compairsy1
fxtb_compair1





# time to equalibrium -----------------------------------------------------
capadjust <- 36/1.69mean(c(1.69,5.88))  #mean(c(37/7.3, 25/mean(c(1.69,5.88))))
51.5

Nest <- meandata$captures[grepl('herm', meandata$species)]*capadjust

em.equal <- function(m, ne) log(1/2)/(log(((1-m)^2)*(1-(1/(2*ne)))))



mi <- 0.4
ne <- 30
te <- log(1/2)/(log(((1-mi)^2)*(1-(1/(2*ne)))))

em.equalm <- function(te, ne) 1-sqrt(exp(log(1/2)/te)/(1-(1/(2*ne))))

em.equalm(1, nseq)

em.equal(0.1, Nest)
em.equal(seq(0.2,1,0.1), 100)

em.fst <- function(Nm) 1/(4*Nm+1)

time2equalibrium <- meandata %>% 
  filter(grepl('herm', species)) %>% 
  ungroup() %>% 
  group_by(phase) %>% 
  summarise(captures = ceiling(mean(captures)),
            mean.fst = round(mean(mean.fst),4)) %>%
  mutate(adjustment = round(capadjust),
         Nest = round(captures*capadjust)) %>%
  mutate(equ0.01 = em.equal(0.01, Nest),
         equ0.1 = em.equal(0.1, Nest),
         equ0.2 = em.equal(0.2, Nest),
         equ0.3 = em.equal(0.3, Nest),
         equ0.5 = em.equal(0.5, Nest)) %>% 
  pivot_longer(cols = equ0.01:equ0.5, names_to = 'm',
               values_to = 't') %>% 
  mutate(m = sub('equ', '', m),
         m = as.numeric(m),
         t = round(t, 1),
         phase = ifelse(duplicated(phase), NA, as.character(phase)),
         captures = ifelse(duplicated(captures), NA, captures),
         fst.est = round(em.fst(Nest*m),3),
         Nest = ifelse(duplicated(Nest), NA, Nest),
         adjustment = ifelse(duplicated(mean.fst), NA, adjustment),
         mean.fst = ifelse(duplicated(mean.fst), NA, mean.fst)) %>%
  relocate(adjustment, Nest, .after = captures)


time2equalibrium_fx <-  time2equalibrium %>% 
  flextable() %>% 
  autofit() %>% 
  theme_alafoli() %>% 
  #border_outer() %>% 
  bold(i = c(2,7,14)) %>% 
  hline(i = seq(5,14,5),
        border = fp_border(color = "grey70", width = 2)) %>% 
  vline(j = 5,
        border = fp_border(color = "grey90", width = 1.5)) %>% 
  hline(part = 'header',
        border = fp_border(color = "grey40", width = 2)) %>% 
  hline_bottom(border = fp_border(color = "grey40", width = 2))

time2equalibrium_fx

# time ------------

caps <- em.rain_caps_phase() %>% 
  mutate(phaseNo = factor(phaseNo, levels =paste0(rep(c('L', 'I', 'D'), 3),
                                                  rep(1:3, each = 3))))

captures <- caps %>% 
  group_by(phaseNo) %>% 
  summarise(captures = mean(captures, na.rm = T),
            capturesSy = mean(capturesSy, na.rm = T)) %>% 
  filter(complete.cases(phaseNo)) %>% 
  pivot_longer(cols = captures:capturesSy, 
               names_to = 'species', values_to = 'captures') %>% 
  mutate(species = ifelse(species == 'captures', 'Pseudomys hermannsburgensis', 
                           'Sminthopsis youngsoni')) %>% 
  arrange(species, phaseNo)

captures
mseq <- seq(0.01, 2, 0.00001)
nseq <- c(222, 89, 30) 
nseq <- rep(c(25,100,75),3) # 
nseq[5] <- 200

capadjust <- 25/1.69

time2fst <- fstdata_NNm |>
  group_by(species, phaseNo, phase) |>
  summarise(mean.fst = mean(fst)) %>% 
  #filter(grepl('herm', species)) %>% 
  ungroup() %>% 
  right_join(captures) %>% 
  group_by(species, phaseNo,phase) %>% 
  summarise(captures = round(mean(captures),2),
            mean.fst = round(mean(mean.fst),4)) %>%
  ungroup() %>% 
  mutate(adjustment = round(capadjust,2),
         Nest = round(captures*capadjust),
         Nest = ifelse(grepl('young', species), NA, Nest))
time2fst$phase[is.na(time2fst$phase)] <- 'stable'
time2fst$m <- NA
time2fst$t <- NA
#time2fst$Nest <- c(100,50,25)
for (i in 1:nrow(time2fst)) {
  d <- time2fst[i,]
  
  ft <- em.fst(d$Nest*mseq) 
  
  fttest <- round(ft, 4) == d$mean.fst
  table(fttest)
  d$m <- mean(mseq[fttest])
   
  time2fst[i,]$m <- d$m
  time2fst[i,]$t <- em.equal(d$m, d$Nest)
}

time2fst %>% 
  mutate(Nm = Nest*m,
         Nmf = em.Nm(mean.fst))
time2fstfx<- time2fst %>% 
  mutate(Nm = round(em.Nm(mean.fst),1),
         m = round(m,2),
         t= round(t,1),
         species = ifelse(grepl('herm', species),
                          'P. hermann', 'S. young'),
         species = ifelse(duplicated(species), NA, species)) %>% 
  dplyr::select(-adjustment) %>% 
  relocate(mean.fst, Nm, .after = Nest) %>%
  flextable() %>% 
  autofit() %>% 
  theme_alafoli() %>% 
  italic(j = 8, i = which(time2fst$m > 1)) %>% 
  bold(i = which(time2fst$t < 1), j= -1) %>%
  hline(i = 9,border = fp_border(color = "grey70", width = 1.5)) %>%  
  vline(j = 6,
        border = fp_border(color = "grey90", width = 1.5)) %>% 
  hline_bottom(border = fp_border(color = "grey40", width = 2)) %>% 
  hline_top(border = fp_border(color = "grey40", width = 2)) %>% 
  italic(j = 6:9, part = 'header') %>% 
  italic(j = 1)

time2fstfx
real <- data.frame(phase = time2equalibrium$phase[complete.cases(time2equalibrium$mean.fst)],
                   fst = time2equalibrium$mean.fst[complete.cases(time2equalibrium$mean.fst)],
                   N = nseq) %>% 
  mutate(m = em.Nm(fst)/N,
         t = em.equal(m, N))

# figure ------------------------------------------------------------------
mseq <- seq(0.01, 0.99, 0.01)
nseq <- c(222, 89, 30) 
nseq <- c(100,50,25) # 

real <- data.frame(phase = time2equalibrium$phase[complete.cases(time2equalibrium$mean.fst)],
                   fst = time2equalibrium$mean.fst[complete.cases(time2equalibrium$mean.fst)],
           N = nseq) %>% 
  mutate(m = em.Nm(fst)/N,
    t = em.equal(m, N))

time2equalibrium
time2df <-data.frame(m = mseq,
           N = rep(nseq, each = length(mseq))) %>% 
  mutate(t = em.equal(m, N),
         fst = em.fst(m*N)) 

time2df %>% 
  ggplot(aes(fst, t,colour = m, group = m))+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(size = 2)+
  facet_grid(~N)+
  geom_vline(data = real, aes(xintercept = fst),
             linetype = 1, linewidth = 1, alpha = 0.5,
             colour = 'red')+
  theme_bw()

time2df %>% 
  ggplot(aes(m,t,colour = factor(N), group = N))+
  #scale_x_log10()+
  #scale_y_log10()+
  scale_y_continuous(breaks = seq(0,6,1))+
  geom_hline(yintercept = c(1), colour = 'grey', linetype = 1,
             linewidth= 0.75)+
  geom_line(linewidth = 1.5)+
  facet_grid(~N)+
  scale_color_manual(values = rev(unname(phaseCol[1:3])))+
  geom_point(data = real, size = 3, colour = 'black',
             aes(shape = factor(fst)))+
  theme_bw()+
  theme(panel.grid = element_blank())

time2df %>% 
  mutate(t1 = t < 1) %>% 
  ggplot(aes(m,fst,colour = t1, group = N))+
  #scale_x_log10()+
  #scale_y_log10()+
  geom_point()+
  facet_grid(~N)+
  #scale_color_manual(values = rev(unname(phaseCol[1:3])))+
  geom_point(data = real, size = 3, colour = 'black',
             aes(shape = factor(fst)))+
  theme_bw()+
  theme(panel.grid = element_blank())


# save tables -------------------------------------------------------------

save_as_docx(fxtb_fst, fxtb_r_npp, anovafx, anovafxsy,
             fxtb_compair,fxtb_compairsy,
             fxtb_compair1, fxtb_compairsy1, time2equalibrium_fx,
             time2fstfx,
             path = './figures/model_tables.docx')



# load --------------------------------------------------------------------

source('./code/libraries.R')
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
aicph <- AIC(mnm, mn, mm,mnmp, mnp, mmp) %>% 
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
aicph
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
  left_join(aicph[,3:6]) %>% 
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
    value = as_paragraph(('log'), as_sub('e'), ('NPP'))) -> fxtb_r_npp;fxtb_r_npp

save_as_docx(fxtb_r_npp, path = './figures/fst_models.docx')

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

ibdcorrsex %>% 
  bind_rows(ibdcorr) %>% 
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


### models ------------------------------------------------------------------




### flextables --------------------------------------------------------------



# ibd models --------------------------------------------------------------



mod_r_ph <- lm(Mantel.cor ~ captures + log(npp) , #weights = (n.dist),
               data = filter(xibdcorr, grepl('herm', species))) 
mod_r_sy <- lm(Mantel.cor ~ captures+ log(npp), 
               data = filter(xibdcorr, grepl('young', species)))

sjPlot::tab_model(mod_r_ph, mod_r_sy, digits = 3,
                  df.method = 'satterthwaite')



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
    value = as_paragraph(('log'), as_sub('e'), ('NPP'))) -> fxtb_r_npp;fxtb_r_npp


### anovas -------------
mod_r_ph <- lm(Mantel.cor ~ captures + log(npp) , #weights = (n.dist),
               data = filter(xibdcorr, grepl('herm', species))) 
mod_r_sy <- lm(Mantel.cor ~ captures+ log(npp), 
               data = filter(xibdcorr, grepl('young', species)))

m5sy <- lm(Mantel.cor ~ log(npp) +captures, #weights = (n.dist),
               data = filter(xibdcorr, grepl('herm', species))) 
m5bsy <- lm(Mantel.cor ~ log(npp) + captures, 
               data = filter(xibdcorr, grepl('young', species)))



rbind(anova(mod_r_ph), anova(m5sy), anova(mod_r_sy),
      anova(m5bsy)) %>% as.data.frame %>% 
  mutate(pvalue = ifelse(`Pr(>F)` < 0.001,'<0.001', 
                         round(`Pr(>F)`, 4)),
         anova = rep(c(1:2,1:2), each = 3),
         parameter = rownames(.),
         parameter = sub('[[:digit:]]', '', parameter),
         species = rep(c('P. hermans', 'S. young'),each = 6),
         `F value` = round(`F value`, 1)) %>% 
  dplyr::select(-`Pr(>F)`) %>% 
  relocate(species, anova, parameter) %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(anova = ifelse(duplicated(paste(anova, species)), NA, anova),
         species = ifelse(duplicated(species), NA, species)) %>% 
  flextable() %>% 
  autofit() %>% 
  # align(j = 1, align = 'right') %>%
  #  align(j = 1, align = 'right', part = 'header') %>% 
  border_remove() %>% 
  bold(part = 'header') %>%
  border_outer() %>% 
  bg(bg = 'grey95', j = 1) %>%
  align(j = 3, align = 'right') %>% 
  # bg(bg = 'grey99', j = 2:3) %>%
  hline(i = seq(3, 9,3), j = -1,
        border = fp_border_default(width = 1, col = 'grey70')) %>% 
  hline(part = 'header', border = fp_border_default(width = 3, col = 'grey70')) %>% 
  hline(i = 6, border = fp_border_default(width = 2, col = 'grey70'))



# sex-density-dependence --------------------------------------------------

phsex <- ind %>% filter(sex_pairs %in% c('f-f', 'm-m'),
               species2== 'ph', km <= 1) %>% 
  mutate(phase = factor(phase, levels = c('low', 'increase', 'decrease')))

sysex <- ind %>% filter(sex_pairs %in% c('f-f', 'm-m'),
                        species2== 'sy', km <= 1)

names(phsex)

m<-lmer(euclidean ~ km * sex_pairs *  phase + (1|phaseNo), data = phsex)

m2 <-lmer(euclidean ~ km * sex_pairs + (1|phaseNo), data = sysex)

anova(m)
anova(m2)

summary(m2)

anova(lm(iris$Petal.Length ~ iris$Sepal.Length + iris$Species))
#https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/
library(emmeans)

emmeans(m2, specs = pairwise ~ km:sex_pairs)

differences <- emmeans(m, specs = pairwise ~ km:sex_pairs:phase)

emmeans(m, specs = pairwise ~ sex_pairs:phase)
emmeans(m, specs = pairwise ~ km:sex_pairs)

emmeans(m, specs = pairwise ~ km:phase)

ggplot(filter(ind, km <= 5), aes(km, kinship, colour = npp))+
  geom_point()+
  facet_wrap(~species, scale = 'free_x')


diffcontrasts <- as.data.frame(differences$contrasts) %>% 
  mutate(km.mean = 0.25,
         contrast = gsub('km0.253353317056273', '', contrast),
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


diffcontrasts %>% 
  mutate(km.mean = ifelse(duplicated(km.mean), NA, km.mean),
         sex.comparison = ifelse(duplicated(sex.comparison), NA, sex.comparison)) %>% 
  flextable() %>% 
  autofit() %>% 
  theme_alafoli() %>% 
  bold(j = 3:7, i = which(abs(diffcontrasts$z.ratio) >2.9)) %>% 
  hline(j = -1, i = c(3,6), border= fp_border_default(width = 1.5, color = 'grey', style = 'dashed'))



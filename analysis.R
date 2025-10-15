
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


ibdcorrsex <- read.csv('./output/ibd_sex_analysis.csv') %>% 
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
  summarise(cor = cor(captures, npp.log)) %>% 
  mutate(r = cor^2)

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

figNPP <- ggplot(fstdata, aes(y = fst, x = npp, fill = phase)) +
  stat_smooth(data = filter(fstdata, ph), show.legend = F,
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
                           data = filter(fstdata,
                                         species == "Pseudomys hermannsburgensis")))

m_fst_anova_sy <- anova(lm(fst ~ captures + npp.log,
                           data = filter(fstdata,
                                         species == "Sminthopsis youngsoni")))

m_fst_anova_ph
m_fst_anova_sy

m_fst_ph <- (lmer(fst ~ captures + npp.log + (1|phaseNo),
            data = filter(fstdata,
                          species == "Pseudomys hermannsburgensis")))

m_fst_sy <- (lmer(fst ~ captures + npp.log + (1|phaseNo),
                           data = filter(fstdata,
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
  scale_x_continuous(trans = scales::log_trans(), 
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

# flextables --------------------------------------------------------------

em.coef.models <-function(m, spp.ph = T, modelname = 'Interaction'){
  
  summary(m)$coefficients %>% 
    as.data.frame() %>% 
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


# multiple regression ----
m7 <- (lmer(fst ~ captures + npp.log + (1|phaseNo), 
            data = filter(fstdata, species == "Pseudomys hermannsburgensis")))
m8 <- (lmer(fst ~ captures + npp.log + (1|phaseNo), data = filter(fstdata, species == "Sminthopsis youngsoni")))
anova(m7)
sjPlot::tab_model(m7, m8, dv.labels = c('P. hermannsburgensis', 'S. youngsoni'),
                  title = 'Fst', df.method = "satterthwaite",
                  digits = 3,collapse.ci = T)

m5 <- (lm(fst ~ captures + npp.log, 
            data = filter(fstdata, species == "Pseudomys hermannsburgensis")))
m5b <- (lm(fst ~ npp.log+captures, 
          data = filter(fstdata, species == "Pseudomys hermannsburgensis")))

m5sy <- (lm(fst ~ captures + npp.log, 
          data = filter(fstdata, species == "Sminthopsis youngsoni")))
m5bsy <- (lm(fst ~ npp.log+captures, 
           data = filter(fstdata, species == "Sminthopsis youngsoni")))


m6 <- (lm(fst ~ captures + npp.log, data = filter(fstdata, species == "Sminthopsis youngsoni")))

am <- anova(m5)
anova(m6)
am2 <- anova(m5b)

rbind(am, am2) %>% as.data.frame %>% 
  mutate(Pr = ifelse(`Pr(>F)` < 0.001,'<0.001', 
                     round(`Pr(>F)`, 4)),
         anova = rep(paste0('captures', 1:2), c(nrow(am), nrow(am2))),
         parameter = rownames(.),
         parameter = sub('1', '', parameter)) %>% 
  rename(sumsq = `Sum Sq`) %>% 
  dplyr::select(anova, parameter, sumsq) %>% 
  filter(parameter != 'Residuals') %>% 
  pivot_wider(names_from = anova, values_from = sumsq)

am$`Pr(>F)` <- ifelse(am$`Pr(>F)` < 0.001,'<0.001', 
                      round(am$`Pr(>F)`, 4))

am %>% mutate_if(is.numeric, round, 4) %>% 
  mutate(Parameter = rownames(.)) %>% 
#  mutate(Parameter = c('Captures', 'NPP (log-scale)', 'Residuals')) %>% 
  relocate(Parameter) %>% 
  flextable() %>% 
  autofit() %>% 
  # align(j = 1, align = 'right') %>%
  #  align(j = 1, align = 'right', part = 'header') %>% 
  border_remove() %>% 
  bold(part = 'header') %>%
  border_outer() %>% 
  bg(bg = 'grey95', j = 1) %>% 
  hline(part = 'header', border = fp_border_default(width = 3, col = 'grey70'))


  rbind(am, am2, anova(m5sy), anova(m5bsy)) %>% as.data.frame %>% 
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
  

em.mRate(0.1, 20)


em.mRate <- function(fst, ne) (1/fst-1)/(4*ne)
em.fst <- function(Nm) 1/(4*Nm+1)
em.Nm <- function(fst) (1/fst-1)/4
em.mRate(0.1,100)
em.Nm(0.1)
em.fst(100*0.0225)
em.Nm(0.18)
em.Nm(0.08)


# lmer table --------------------------------------------------------------

modeltbsx <- rbind(em.coef.models(m7, T,'P. hermann'),
                   em.coef.models(m8, F, 'S. youngsoni')) %>% 
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
    value = as_paragraph(('Intercept'))) %>% 
  compose(
    part = "body", j = 2, i = grep('cap', mtb$Parameter),
    value = as_paragraph(('Captures'))) %>% 
  compose(
    part = "header", j = 4,
    value = as_paragraph(('R'), as_sup('2'))) %>% 
  compose(
    part = "body", j = 2, i = grep('log', mtb$Parameter),
    value = as_paragraph(('log'), as_sub('e'), ('NPP'))) -> fxtb_r_npp;fxtb_r_npp


# ibd models --------------------------------------------------------------
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



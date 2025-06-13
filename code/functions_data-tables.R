
em.period_table <- function(fstdata){ 
period_summaries <- fstdata %>%
  group_by(species, phaseNo, phase, npp, n, ngrids) %>% 
  summarise(nPairs = n(),
            sd = sd(fst),
            fst = mean(fst),
            se = sd/sqrt(nPairs),
            lowerf = round(fst - se*1.96,3),
            upperf = round(fst + se*1.96,3),
            fstsum = paste0(round(fst, 3), " (", lowerf, ' - ' ,upperf, ")" ),
            fstsum = ifelse(is.na(se), as.character(round(fst,3)), fstsum)) %>%
  dplyr::select(-lowerf, -upperf) %>% 
  ungroup() %>% 
  mutate(npp.log = log(npp)) %>% 
  mutate(#mRate = (1/ifelse(fst > 0, fst, NA)-1)/(4*(ne)),
         nMigrants = (1/ifelse(fst > 0, fst, NA)-1)/4,
         #mRate = round(mRate, 3),
         nMigrants = round(nMigrants),
         npp = round(npp,1),
         nPairs = ifelse(is.na(ngrids), NA, nPairs)) %>% 
  arrange(species, phaseNo)  %>% 
  mutate(fst = round(fst, 3)) %>% 
  dplyr::select(species,phaseNo,  phase, npp, nPairs,
                #fst,
                fstsum,
                n,
                nMigrants) ; period_summaries


periodft <- period_summaries %>% 
  mutate(species = ifelse(grepl('hermann', species),
                          'P. hermanns',
                          'S. youngsoni'))

colnames(periodft) <- c('species', 'period','phase',  'NPP', 
                        'ngrid', 'FST', 'n', 'migrants')

return(periodft)
}


em.coef.models <-function(m, spp.ph = T, modelname = 'Interaction'){
  
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
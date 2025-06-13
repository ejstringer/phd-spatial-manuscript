# analysis 2 - structure and dispersal

# create data

# filter data
source('./code/libraries.R')
source('./code/functions_genetic-data-creation.R')
source('./code/functions-npp-data-creation.R')
source('./code/functions_correlograms.R')


# filtering ---------------------------------------------------------------

ph <- em.filtering('ph')
sy <- em.filtering('sy')


# sample size -------------------------------------------------------------

neprep <- lapply(list(pherm = ph, syoung = sy), em.ne_prep)
n <- sapply(neprep, function(x) sapply(x, nInd)) %>% 
  as.data.frame() %>% 
  mutate(period = rownames(.),
         phaseNo = paste0(rep(c('L', 'I', 'D'), 3), rep(1:3, each = 3))) %>% 
  pivot_longer(pherm:syoung, names_to = 'species', values_to = 'n') %>% 
  mutate(species = ifelse(species == 'pherm', 'Pseudomys hermannsburgensis',
                          'Sminthopsis youngsoni')) %>% 
  arrange(species)
write.csv(n, './output/sampe_sizes.csv', row.names = FALSE)

# npp ---------------------------------------------------------------------


npp <- em.npp_layers(ph)

nppmean <- read.csv('./output/npp_means.csv') 


# period ------------------------------------------------------------------
periodID <- data.frame(phaseNo = paste0(rep(c('L', 'I', 'D'), 3),
                                        rep(1:3, each = 3)),
                       period = paste0('period ', 1:9))


# captures ----------------------------------------------------------------

caps <- em.rain_caps_phase() %>% 
  mutate(phaseNo = factor(phaseNo, levels =paste0(rep(c('L', 'I', 'D'), 3),
                                                  rep(1:3, each = 3))))

captures <- caps %>% 
  group_by(phaseNo, period, phase) %>% 
  summarise(captures = mean(captures, na.rm = T),
            capturesSy = mean(capturesSy, na.rm = T)) %>% 
  filter(complete.cases(phaseNo)) %>% 
  pivot_longer(cols = captures:capturesSy, 
               names_to = 'species2', values_to = 'captures') %>% 
  mutate(species2 = ifelse(species2 == 'captures', 'ph', 'sy')) %>% 
  arrange(period)


# genetics ----------------------------------------------------------------

system.time(gdist <- lapply(list(pherm = ph, syoung = sy), 
                            em.individual_relatedness)) # ~12mins

system.time(fst <- lapply(list(pherm = ph, syoung = sy), 
                          em.phase_fst)) # ~22mins


## individual --------------------------------------------------------------

system.time(individual <- gdist$pherm %>% 
  bind_rows(gdist$syoung) %>% 
    filter(within) %>% 
    left_join(nppmean, by = c('period.x' = 'period')) %>% 
    rename(phase = phase.x, period = period.y) %>% 
    left_join(periodID) %>% 
  dplyr::select(species, period, phaseNo, phase, kinship, ibdecent, euclidean, 
                metres, km, npp, mc, pairs, idpairs, trip.x, trip.y, 
                sex.x, sex.y, id.x, id.y)) 
names(individual)
head(individual)

individual$npp.log = log(individual$npp)
individual$phase <- factor(individual$phase,
                             levels = c('low', 'increase', 'decrease'))

individual$phaseNo <- factor(individual$phaseNo,
                           levels =  periodID$phaseNo)

individual$rep <- as.numeric(str_sub(individual$phaseNo, 2,2))

individual$daydist <- abs(ymd(individual$trip.x) - ymd(individual$trip.y))
individual$sex_pairs <- paste0(individual$sex.x, '-', individual$sex.y)

individual$sex_pairs <- ifelse(individual$sex_pairs %in% c('f-m', 'm-f'), 'f-m',
                               individual$sex_pairs)

individual$sex_pairs <- ifelse(grepl('NA', individual$sex_pairs), 'unknown', 
                               individual$sex_pairs)

individual$sex_pairs %>% table

individual_final <- individual %>% dplyr::select(-c(sex.x, sex.y, trip.x, trip.y)) %>% 
  mutate(species2 = ifelse(grepl('herman', species), 'ph', 'sy')) %>% 
  left_join(captures) %>% 
  mutate(km.log = ifelse(km == 0, log(km + 0.1), log(km)),
         phase = ifelse(species == 'Sminthopsis youngsoni', 'stable', phase),
         phase = factor(phase, levels = c('low', 'increase',
                                          'decrease','stable')))
glimpse(individual_final)

## fst ---------------------------------------------------------------------

ngrids <- c(
  sapply(split(fst$pherm, fst$pherm$period), function(x)
    length(unique(c(x$pop1, x$pop2)))),
  sapply(split(fst$syoung, fst$syoung$period), function(x)
    length(unique(c(filter(x, grepl('you', species))$pop1,
                    filter(x, grepl('you', species))$pop2))))) %>% 
  data.frame(ngrids = ., period = names(.), 
             species = rep(c('Pseudomys hermannsburgensis',
                             'Sminthopsis youngsoni'), each = 9)[-18])

fst_final <- fst$syoung %>% 
  bind_rows(fst$pherm) %>% 
  left_join(ngrids) %>%
  left_join(nppmean) %>%  
  left_join(periodID)  %>%   
  mutate(npp.log = log(npp),
         fst = ifelse(fst < 0, 0, fst),
         fst.log = log(fst + 0.00001),
         phase = factor(phase, levels = c("low", 'increase','decrease')),
         phaseNo =  factor(phaseNo, levels = periodID$phaseNo),
         rep = as.numeric(str_sub(phaseNo, 2,2)),
         species2 = ifelse(grepl('Pse', species), 'ph', 'sy'),
         # mRate = NA, #(1/ifelse(fst > 0, fst, NA)-1)/(4*(ne)),
         # mRate.log = NA, # log(mRate),
         nMigrants = (1/ifelse(fst > 0, fst, NA)-1)/4,
         nMigrants.log = log(nMigrants)
  ) %>% 
  left_join(captures)

head(as.data.frame(fst_final))


# correlograms ------------------------------------------------------------
ibdDataAll <- lapply(c('herm', 'young'), 
                     function(x) filter(individ, grepl(x, species),
                                        complete.cases(metres)))

## mantel correlog --------
system.time(ibdcorr <- lapply(ibdDataAll, 
                              em.ibd.period.correlog, 
                              test = 'pearson')) # 24 mins
ibdcorr_join <- do.call('bind_rows', ibdcorr)

# save --------------------------------------------------------------------

write.csv(individual_final, './output/individual_analysis.csv', row.names = F)
write.csv(fst_final, './output/fst_analysis.csv', row.names = FALSE)

write.csv(ibdcorr_join, './output/ibd_analysis.csv', row.names = FALSE)


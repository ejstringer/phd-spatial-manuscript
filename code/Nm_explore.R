
em.mRate <- function(fst, ne) (1/fst-1)/(4*ne)
em.fst <- function(Nm) 1/(4*Nm+1)
em.Nm <- function(fst) (1/fst-1)/4

em.N <- function(fst, m) (1/fst-1)/(4*m)

fst <- c(0.01, 0.02, 0.1, 0.2)
N <- c(25:100)

m <- sapply(fst, function(x) em.mRate(x, N))
colnames(m) <- paste0('fst', fst)
#m <- 0.05
em.N(0.2, 0.05)

Ne <- sapply(1:4, function(x) em.N(fst[x], m[,x]))

fst
plot(Ne,m, col = rep(1:4, each = nrow(m)), pch = 16)

cbind(m,N) %>% as.data.frame() %>% 
  pivot_longer(-N) %>% 
  mutate(name = (sub('fst','', name))) %>% 
  rename(fst =name, m = value) %>% 
  ggplot(aes(N, m, colour = fst, group = fst)) +
  geom_line(size = 1)+
  theme_classic()

fstmin <- min(fstdata$fst[fstdata$fst>0])
fstdata$Nm <- em.Nm(fstdata$fst+fstmin)

ggplot(filter(fstdata, fst > fstmin), aes(npp.log,Nm , color = phase))+
  geom_point()+
  #scale_y_log10()+
  #scale_x_log10()+
  geom_smooth(method = 'lm')+
  facet_grid(~species)

fstph<-fstdata %>% 
  group_by(phase, phaseNo, npp.log, captures, species) %>% 
  summarise(fst = mean(fst)) %>% 
  mutate(Nm = em.Nm(fst),
         N = captures*15,
         m = em.mRate(fst, N)) %>% 
  filter(!is.infinite(Nm), grepl('herm', species)) 

ggplot(fstph,aes(captures*15, Nm, colour = phase))+ 
  geom_point()+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3),
              se = F, aes(group = species))
  
ggplot(fstph,aes(npp.log, Nm, colour = phase))+ 
  geom_point()+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 5),
              se = F, aes(group = species))


ggplot(fstph,aes(N, m, colour = phase))+ 
  #geom_line(aes(group = species), linewidth = 1, colour = 'black')+
  geom_point(size = 3)+
  geom_smooth(span = 0.9, aes(group = species), se = F) +
  theme_classic()



AIC(lm(Nm ~ captures, data = fstph), lm(log(Nm) ~ captures,data = fstph))


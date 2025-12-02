
em.mRate <- function(fst, ne) (1/fst-1)/(4*ne)
em.fst <- function(Nm) 1/(4*Nm+1)
em.Nm <- function(fst) (1/fst-1)/4

em.N <- function(fst, m) (1/fst-1)/(4*m)
em.m <- function(fst, N) (1/fst-1)/(4*N)



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


# t explore ---------------------------------------------------------------


em.mRate <- function(fst, ne) (1/fst-1)/(4*ne)
em.fst <- function(Nm) 1/(4*Nm+1)
em.Nm <- function(fst) (1/fst-1)/4

em.N <- function(fst, m) (1/fst-1)/(4*m)

em.fstiso <- function(ne, t = 2) 1 - (1 - 1 / (2 * ne))^t


d<-fstdata %>% 
  group_by(species, phaseNo, phase, captures, npp.log, npp) %>% 
  summarise(fst = mean(fst),
            Nm = em.Nm(fst)) %>%
  mutate(mSrgt = Nm/(captures)) %>% 
  # filter(phase != 'stable') %>% 
  arrange(phase)

e <- d %>% filter(grepl('herm', species))

est.fst <- apply(data.frame(n = e$captures[4:9],
                            t = c(rep(1,3), 2, 2, 2)), 1,
                 function(x) em.fstiso(x[1], x[2]))


est.fst2 <- est.fst^2
lm(e$fst[4:9]~ est.fst) %>% summary

plot(e$fst[4:9]~ est.fst, col = e$phase[4:9], pch = 16)



em.equal <- function(m, ne) log(1/2)/(log(((1-m)^2)*(1-(1/(2*ne)))))
em.equal(seq(0,1,0.1), 25)
em.equal(seq(0.2,1,0.1), 100)

# Define parameters
m <- 0.01    # Migration rate
N <- 1000    # Population size

em.equal2 <- function(m, N){
  
  # Calculate the expression
numerator <- log(1 / 2)
denominator <- log((1 - m)^2 * (1 - 1 / (2 * N)))
result <- numerator / denominator

result
}
m.mice <- seq(0.02,0.5,0.02)
data.frame(ne = rep(c(10,25,100,250), each = length(m.mice)),
           m = m.mice) %>% 
  mutate(tte = em.equal2(m, ne)) %>% 
  ggplot(aes(m, tte, colour = ne, group = ne))+
  geom_point()+
  geom_smooth(se = F)+
  theme_classic()+
  geom_hline(yintercept = c(1,4), lty = 2, colour = 'grey')


em.equal(0.01, 250)

em.equal2(0.0001, 100)
em.equal2(0.0001, 10000)

em.equal2(0.1, 10000)
em.equal2(0.1, 100)

em.equal2(0.1, 10000)
em.equal2(0.1, 100)

em.equal2(0.8, 100)
em.equal2(0.1, 100)

em.equal2(0.1, 1000)
em.equal2(0.3, 100)




e %>% 
  mutate(tte = em.equal())





e %>% 
  mutate(m = em.mRate(fst, captures),
         m2 = (1-fst)/(4*captures*fst))
em.mRate()





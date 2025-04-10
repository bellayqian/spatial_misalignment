library(sf)
library(ggplot2)
library(USAboundaries)
library(RColorBrewer)
library(cowplot)

load('../Dashboard/data/final_misaligned_data.RData')

stbounds<-USAboundaries::us_states()
stbounds<-subset(stbounds,!(stusps %in% c("AK","HI","PR")))

cd<-us_congressional()
cd<-subset(cd,!(state_abbr%in% c('AK','HI','PR')))

aa<-brewer.pal(n = 8, name = "Dark2")

## make map of counties intersected by atoms ##
countybds<-ggplot() +
  geom_sf(data = county,fill='white',lwd=.2,color=aa[2])+
  geom_sf(data = stbounds,lwd=.2,color='black',fill=NA)+
  theme_void()

cdbds<-ggplot() +
  geom_sf(data = cd,fill='white',lwd=.2,color=aa[3])+
  geom_sf(data = stbounds,lwd=.2,color='black',fill=NA)+
  theme_void()

county$num_atoms_cat<-cut(county$num_atoms,breaks=c(0,1.5,2.5,5.5,10.5,20),labels=c('1','2','3-5','6-10','10+'))

numat_fill<-ggplot() +
  geom_sf(data = county,aes(fill=num_atoms_cat),lwd=.2,color=aa[2])+
  geom_sf(data = stbounds,lwd=.2,color='black',fill=NA)+
  scale_fill_brewer(palette = "YlGnBu")+
  theme_void()+
  labs(fill='Number of\nAtoms')

top_row <- plot_grid(countybds, cdbds, labels = c('A', 'B'), label_size = 12)
# pdf('/Users/rachel/OneDrive - Harvard University/Grant_Apps/Other_PI/krieger_feb23/working_paper/figures/maps_together.pdf',
#     width=10)
# print(cowplot::plot_grid(top_row,numat_fill,labels = c('', 'C'), label_size = 12, ncol = 1))
# dev.off()

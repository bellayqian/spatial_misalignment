
## compute adusted poverty counts using the county level denominators from the census
atom_noshp<-atom
st_geometry(atom_noshp)<-NULL
dc_county_pop<-aggregate(pop_total~county_geoid,data=atom,FUN=sum,na.rm=T)
names(dc_county_pop)<-c('geoid','pop_total')
county<-merge(county,dc_county_pop,by='geoid')
county$poverty_count_adjusted<-round((county$poverty_pct/100)*county$pop_total)

## create age-standardized denominators at the atom level using the county age distributions ##

## get age distribution in proportions (and remove counts and denominator) ##
atom_new<-atom
atom_new<-atom_new%>%mutate(
  age_denom=age0_24+age25_44+age45_64+age65_74+age75plus
)%>%mutate(
  age0_24_p=age0_24/age_denom,
  age25_44_p=age25_44/age_denom,  
  age45_64_p=age45_64/age_denom,
  age65_74_p=age65_74/age_denom,
  age75plus_p=age75plus/age_denom
) %>%select(-c(age0_24,age25_44,age45_64,age65_74,age75plus,age_denom))

## get age group specific counts for the census data by multiplying total atom total population counts by age proportions
atom_new<-atom_new%>%mutate(
  age0_24=pop_total*age0_24_p,
  age25_44=pop_total*age25_44_p,
  age45_64=pop_total*age45_64_p,
  age65_74=pop_total*age65_74_p,
  age75plus=pop_total*age75plus_p
)

## christian's function for getting pop standards ##
source('pull_pop_std_rates.R')
std<-load_std_population_sizes()

## now use these to get age-standardized denominators ##
atom_new<-atom_new%>%mutate(
  expected_count=age0_24*std$rate[1]+age25_44*std$rate[2]+age45_64*std$rate[3]+age65_74*std$rate[4]+age75plus*std$rate[5]
)

atom<-atom_new

## remove atoms with zero population ##
atom<-subset(atom,!(pop_total==0))

## recompute the num_atoms variable now because we are throwing out atoms without population ##
countyXatom<-data.frame(table(atom$county_geoid))
names(countyXatom)<-c('geoid','num_atoms')

## merge the number of atoms for each county into the county level dataset ##
county$num_atoms<-NULL
county<-merge(county,countyXatom,by='geoid')
county<-county[order(county$num_atoms,county$geoid),]

## merge the number of atoms for each county into the atom level dataset ##
atom$num_atoms<-NULL
atom<-merge(atom,countyXatom,by.x='county_geoid',by.y='geoid')
atom<-atom[order(atom$num_atoms,atom$county_geoid),]
atom$out_ind<-as.numeric(factor(atom$county_geoid,levels=county$geoid))

## compute population density at the atom level
atom$atom_area<-as.numeric(st_area(atom))
atom$popdensity<-100000*atom$pop_total/atom$atom_area

## probably have to round death counts since we're considering them poisson ##
county$deaths_round<-round(county$deaths)

## save all the misaligned data ##
save(county,atom,file='final_misaligned_data.RData')

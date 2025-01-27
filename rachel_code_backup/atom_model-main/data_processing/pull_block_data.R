## script to get atom-level variables and insert them into the atom dataset ##

variables_dict <-
  tibble::tribble(
    ~var,          ~shortname,      ~desc,
    "P1_001N",     'pop_total',                      'total population estimates (race composition table)',
    "P1_003N",     'race_white_alone',               "white alone population estimate",
    "P1_004N",     'race_black_alone',               "black alone population estimate",
    "P1_005N",     'race_am_indian_ak_native_alone', "american indian alone population estimate",
    "P1_006N",     'race_asian_alone',               "asian alone population estimate",
    "P1_007N",     'race_native_hawaiian_opi_alone', "native hawaiian or other pacific islander alone",
    "P2_001N",     'hispanic_denom',                 "denominator population for hispanic estimate",
    "P2_002N",     'hispanic_count',                 "population estimate for hispanic",
  )

atom<-st_transform(atom, 3857)
atom$atomID<-1:nrow(atom)

## set up empty variables in the atom dataset ##
atom$pop_total<-NA
atom$pct_white<-NA
atom$pct_black<-NA
atom$pct_aian<-NA
atom$pct_asian<-NA
atom$pct_hispanic<-NA

newcols<-which(names(atom)=='pop_total'):which(names(atom)=='pct_hispanic')

## most counties are atoms, so pull the county-level dataset for the whole US ##
dc2020_county<- get_decennial(geography = "county",
                              variables = variables_dict$var,year = 2020,sumfile = "pl")

# pivot to a wide format for renaming
dc2020_county <- dc2020_county %>% 
                 pivot_wider(names_from = variable, values_from = value)

rename_vars <- setNames(variables_dict$var, variables_dict$shortname)
dc2020_county <- dc2020_county %>% rename(!!rename_vars)

# process racial composition data
dc2020_county <- dc2020_county %>%
  mutate(
    pct_white = 100*race_white_alone / pop_total,
    pct_black = 100*race_black_alone / pop_total,
    pct_aian = 100*race_am_indian_ak_native_alone / pop_total,
    pct_asian = 100*race_asian_alone / pop_total,
    pct_hispanic = 100*hispanic_count / hispanic_denom) %>%
  select(geoid=GEOID,pop_total,pct_white,pct_black,pct_aian,pct_asian,pct_hispanic)

dc2020_county<-data.frame(dc2020_county)


for (i in 1:nrow(atom)){
  ## for a given county, if it only has 1 atom (county=atom) then pull county level variables ##
  if (atom$num_atoms[i]==1){
    
    atom[i,newcols]<-unlist(c(dc2020_county[which(dc2020_county$geoid==atom$county_geoid[i]),2:ncol(dc2020_county)]))
  
  } else{ 
    
    ## download, unzip, and read block shapefile for the state (these can't be pulled via tidycensus yet, ugh) ##
    #https://www2.census.gov/geo/tiger/TIGER2020/TABBLOCK20/tl_2020_01_tabblock20.zip
    ## if you've already loaded in the state's shapefile, skip this ##
    if (substr(atom$county_geoid[i-1],1,2)==substr(atom$county_geoid[i],1,2)){
      
    }else{
      temp <- tempfile()
      temp2<-tempfile()
      download.file(paste0('https://www2.census.gov/geo/tiger/TIGER2020/TABBLOCK20/tl_2020_',substr(atom$county_geoid[i],1,2),'_tabblock20.zip'),temp)
      unzip(zipfile = temp, exdir = temp2)
      blockshp<-st_read(file.path(temp2, paste0("tl_2020_",substr(atom$county_geoid[i],1,2),"_tabblock20.shp")))
    }
    
    ## if you've already loaded in the county's block level data, don't reload/re-process
    if (atom$county_geoid[i-1]==atom$county_geoid[i]){
      
    }else{
      
      dc2020_block<-get_decennial(geography = "block",
                    state=substr(atom$county_geoid[i],1,2),
                    county=substr(atom$county_geoid[i],3,5),
                    variables = variables_dict$var,year = 2020,sumfile = "pl")
      
      # pivot to a wide format for renaming
      dc2020_block <- dc2020_block %>% 
        pivot_wider(names_from = variable, values_from = value)
      
      dc2020_block <- dc2020_block %>% rename(!!rename_vars)
      
      ## get a county level block shapefile ##
      countyblocks<-subset(blockshp,COUNTYFP20==substr(atom$county_geoid[i],3,5))
      ## merge with the census block data ##
      countyblocks<-merge(countyblocks,dc2020_block,by.x='GEOID20',by.y='GEOID')
      
      ## make block polygons into points at centroids ##
      countyblocks<-st_transform(countyblocks, 3857)
      countyblocks$geom2 = st_point_on_surface(st_geometry(countyblocks))
      st_geometry(countyblocks) <- "geom2"
    }
    
    ## sum values within the atom ##
    pointsID <- subset(st_join(countyblocks, atom[i,c('atomID')]),!is.na(atomID))
    st_geometry(pointsID)<-NULL
    p_ag1 = data.frame(lapply(pointsID[,c("pop_total","race_white_alone", "race_black_alone",         
                                      "race_am_indian_ak_native_alone", "race_asian_alone","race_native_hawaiian_opi_alone",
                                      "hispanic_denom","hispanic_count")], sum,na.rm=T))
    p_ag1 <- p_ag1 %>%
      mutate(
        pct_white = 100*race_white_alone / pop_total,
        pct_black = 100*race_black_alone / pop_total,
        pct_aian = 100*race_am_indian_ak_native_alone / pop_total,
        pct_asian = 100*race_asian_alone / pop_total,
        pct_hispanic = 100*hispanic_count / hispanic_denom
      ) %>% select(pop_total,pct_white,pct_black,pct_aian,pct_asian,pct_hispanic)
    
    atom[i,newcols]<-unlist(c(p_ag1))
  }
}

## saving this dataset to be safe since it takes forever to process ##
save(atom,file='block_data.RData')

# Load packages ----------------------------------------------------------------
sapply(list.files("./R", full.names = T), source)

ipak(c("tidyverse",
       "sf"))

# Load move dataframe and make sf object
move_df <- read.csv("./data/land_use_dependent_dispersal_data/tidy/move_df.csv", row.names = 1)


mollcrs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
move_sf <- move_df %>% st_as_sf(coords = c("location.long", "location.lat"),
                                crs = projcrs) %>%
  st_transform(crs = mollcrs)



# Pull GEE data on NDVI

library("rgee")
ee_Initialize(gcs = T, drive = T)

# Modified from https://gis$stackexchange$com/questions/290892/google-earth-enginendvi-by-month-per-year
ndvi <- ee$ImageCollection("MODIS/061/MOD13Q1")$select("NDVI")

months <- ee$List$sequence(1, 12);
years <- ee$List$sequence(2001, 2021);

byMonthYear <- ee$ImageCollection$fromImages(
  years$map(
    ee_utils_pyfunc(
      function(y) {
        months$map(
          ee_utils_pyfunc(
            function(m) {
              ndvi$
                filter(ee$Filter$calendarRange(y, y, 'year'))$
                filter(ee$Filter$calendarRange(m, m, 'month'))$
                mean()$set('month', m)$set('year', y)
            }
          )
        )
      })
  )$flatten())


ee_print(byMonthYear$
           filter(ee$Filter$eq("month", 2))$
           filter(ee$Filter$eq("year", 2012)))

# Want to take the monthly averages of NDVI from GEE
# First mutate to get a month and year column
move_sf <- move_sf %>%
  mutate(month = substr(timestamp, 6, 7) %>% as.numeric(),
         year = substr(timestamp, 1, 4) %>% as.numeric()) %>%
  mutate(year = ifelse(year < 2001, 2001, year)) %>%
  mutate(year = ifelse(year > 2021, 2021, year)) # This will take the nearest year's monthly average

mm <- 1:12
yy <- 2001:2021
tt <- expand.grid(mm, yy)

# Create an empty column for ndvi that will be filled in with the following code
move_sf <- move_sf %>%
  mutate(ndvi_terra = NA)

# Will start this process up again if we encounter errors
last_i <- 0

#
for(i in (last_i+1):nrow(tt)){ # last_i
  month_i <- tt[i,1]
  year_i <- tt[i,2]

  time_match_inds <- which(move_sf$month == month_i & move_sf$year == year_i)

  if(length(time_match_inds) == 0){
    next()
  }

  # GEE doesn't want you to make more than 5000 requests at a time, I think
  inds_list <- split(time_match_inds, ceiling(seq_along(1:length(time_match_inds))/5000))# %>% unlist() %>% as.vector() %>% tail()

  for(j in 1:length(inds_list)){
    focal_inds <- inds_list[[j]]
    ee_dat <- ee_extract(x = byMonthYear$
                           filter(ee$Filter$eq("month", month_i))$
                           filter(ee$Filter$eq("year", year_i)),
                         y = move_sf[focal_inds,"event.id"],
                         fun = ee$Reducer$mean(),
                         scale = 250,
                         sf = T)

    if(all(ee_dat$event.id == move_sf$event.id[focal_inds])){
      if(length(grep("NDVI", colnames(ee_dat))) > 0){
        move_sf$ndvi_terra[focal_inds] <- ee_dat %>% pull(2)
      } else{
        next() #These should be cases where there are only ocean coordinates
      }
    } else{
      stop() # This would be some sort of error
    }
  }





  last_i <- i
  print(paste("Year:", year_i, " Month:", month_i, " Total coordinates:", length(time_match_inds)))

}
Sys.time() # This takes several hours

# # Write this out so that we have a copy of this in case anything doesn't work
# # with the next steps
# all(move_sf$event.id == move_df$event.id)
# move_df$ndvi_terra <- move_sf$ndvi_terra
# write.csv(move_sf, file = "./data/land_use_dependent_dispersal_data/tidy/move_df_with_ndvi.csv", row.names = F)


# Same thing for ocean (won't end up using this for now, but potentially
# relevant for a future analysis using this framework)

chlor <- ee$ImageCollection("NASA/OCEANDATA/MODIS-Terra/L3SMI")$select("chlor_a")

months <- ee$List$sequence(1, 12);
years <- ee$List$sequence(2001, 2021);

byMonthYear <- ee$ImageCollection$fromImages(
  years$map(
    ee_utils_pyfunc(
      function(y) {
        months$map(
          ee_utils_pyfunc(
            function(m) {
              chlor$
                filter(ee$Filter$calendarRange(y, y, 'year'))$
                filter(ee$Filter$calendarRange(m, m, 'month'))$
                mean()$set('month', m)$set('year', y)
            }
          )
        )
      })
  )$flatten())

ee_print(byMonthYear$
           filter(ee$Filter$eq("month", 2))$
           filter(ee$Filter$eq("year", 2012)))

# Create an empty column for chlor that will be filled in with the following code
move_sf <- move_sf %>%
  mutate(chlor_terra = NA)

# Will start this process up again if we encounter errors
last_i <- 0

#
for(i in (last_i + 1):nrow(tt)){ # last_i
  month_i <- tt[i,1]
  year_i <- tt[i,2]

  time_match_inds <- which(move_sf$month == month_i & move_sf$year == year_i)

  if(length(time_match_inds) == 0){
    next()
  }

  # GEE doesn't want you to make more than 5000 requests at a time, I think
  inds_list <- split(time_match_inds, ceiling(seq_along(1:length(time_match_inds))/5000))# %>% unlist() %>% as.vector() %>% tail()

  # For some reason, data wasn't available for July 2021. Will
  # instead use July 2020 values

  if(month_i == 7 & year_i == 2021){
    year_i <- 2020
  }

  for(j in 1:length(inds_list)){
    focal_inds <- inds_list[[j]]
    ee_dat <- ee_extract(x = byMonthYear$
                           filter(ee$Filter$eq("month", month_i))$
                           filter(ee$Filter$eq("year", year_i)),
                         y = move_sf[focal_inds,"event.id"],
                         fun = ee$Reducer$mean(),
                         scale = 250,
                         sf = T)

    if(all(ee_dat$event.id == move_sf$event.id[focal_inds])){
      if(length(grep("chlor", colnames(ee_dat))) > 0){
        move_sf$chlor_terra[focal_inds] <- ee_dat %>% pull(2)
      } else{
        next() #These should be cases where there are only ocean coordinates
      }
    } else{
      stop() # This would be some sort of error
    }
  }





  last_i <- i
  print(paste("Year:", year_i, " Month:", month_i, " Total coordinates:", length(time_match_inds)))

}
Sys.time() # Started around 5:28 starting with index 48


# Add the chlor values to move_df and write out
all(move_sf$event.id == move_df$event.id)
move_df$chlor_terra <- move_sf$chlor_terra
write.csv(move_df, file = "./data/land_use_dependent_dispersal_data/tidy/move_df_with_gee_covs.csv", row.names = F)


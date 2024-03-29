library(httr)
library(raster)
library(rnpn)
library(lubridate)
#start <- now()
library(tidyverse)
library(leaflet)
library(leaflet.extras)
library(RColorBrewer)
library(jsonlite)
library(ggmap)
library(sf)
library(geodata)
library(shiny)
#setwd('C:/Users/twmat/USA_NPN/Buffelgrass_v2/')



# Functions for calculating forecast values and acquiring RainLog data
wBufferStations <- function(vals, na.rm=T) { # This function uses inches as per RCC-ACIS
  # Initialize event and buffer trackers
  events <- 0
  low_ppt_days <- 4
  # Examine each preivous day's ppt value at a station
  for (ppt in vals) {
    if (!is.na(ppt)) {
      # Treat missing or delayed records as no ppt
      if (ppt == 'M' | ppt == 'S') {
        ppt <- 0
        # Treat trace records as low ppt
      } else if (ppt == 'T') {
        ppt <- 0.01
        # Treat accumulated records as true records
      } else if (substr(ppt, nchar(ppt), nchar(ppt)) == 'A') {
        ppt <- substr(ppt, 1, nchar(ppt)-1)
      }
      # Ensure value is numeric
      ppt <- as.numeric(ppt)
      
      # Case if ppt and buffer exceed thresholds
      if (ppt >= 0.25 & low_ppt_days >= 3) {
        # Count new event and reset buffer
        events <- events + 1
        low_ppt_days <- 0
        # Case if no ppt
      } else if (ppt == 0) {
        # Increment buffer
        low_ppt_days <- low_ppt_days + 1
        # Case if ppt occurred but is below threshold AND at least one no ppt day occurred since last event
      } else if (ppt < 0.25 & low_ppt_days != 0) {
        # Increment buffer
        low_ppt_days <- low_ppt_days + 1
        # Case if ppt occurred (above or below threshold) AND ppt occurred continuously since last event
      } else {
        # Reset buffer
        low_ppt_days <- 0
      }
    }
  }
  # After all days for a station have been checked, normalize scale by resetting
  # event counters above 4 to 4.
  if (events > 4) {
    events <- 4
  }
  return (events)
}

getRainLogReadings <- function(sdate, edate, offset = 0) {
  # This function obtains 1000 entries from RainLog gauges
  # INPUTS:
  # sdate: date as char, Start date of observation window
  # edate: date as char, End date of observation window
  # offset: int as char, index of POST to start at (increment by limit each iteration)
  
  headers <- c(
    'Content-Type' = 'application/json',
    'Accept' = 'application/json')
  
  body = paste0('{
    "quality": ["Good"],
    "pagination": {
      "offset": ',offset,',
      "limit": 1000
    },
    "dateRangeStart": "',sdate,'",
    "dateRangeEnd": "',edate,'",
    "region": {
      "type": "Rectangle",
      "westLng": -114.8154,
      "eastLng": -109.0449,
      "northLat": 31.32917,
      "southLat": 37.00459
    }
  }')
  raw.reading <- POST('https://rainlog.org/api/1.0/Reading/getFiltered', body = body, add_headers(headers))
  reading_data <- fromJSON(rawToChar(raw.reading$content))
  return (reading_data)
}

getRainLogGauges <- function(sdate, edate, offset = 0) {
  # This function obtains 1000 entries from RainLog gauges
  # INPUTS:
  # sdate: date as char, Start date of observation window
  # edate: date as char, End date of observation window
  # offset: int as char, index of POST to start at (increment by limit each iteration)
  
  headers <- c(
    'Content-Type' = 'application/json',
    'Accept' = 'application/json')
  
  body = paste0('{
    "pagination": {
      "offset": ',offset,',
      "limit": 1000
    },
    "dateRangeStart": "',sdate,'",
    "dateRangeEnd": "',edate,'",
    "region": {
      "type": "Rectangle",
      "westLng": -114.8154,
      "eastLng": -109.0449,
      "northLat": 31.32917,
      "southLat": 37.00459
    }
  }')
  raw.gauges <-POST('https://rainlog.org/api/1.0/GaugeRevision/getFiltered', body = body, add_headers(headers))
  gauges_data <- fromJSON(rawToChar(raw.gauges$content))
  return (gauges_data)
}

rainlogPrep <- function(rain_list, dates_list) {
  # Takes as input: readingDate, rainAmount for sub-array (grouped by gaugeId)
  new <- data.frame(readingDate = lookback, rainAmount = 0)
  new[which(as.character(lookback) %in% dates_list),2] <- rain_list
  val <- wBufferStations(new$rainAmount)
  return (val)
}

# Create list of days for forecasting
today <- as_date(today()) 
lookback <- as_date(today - days(30:1))


### PRISM BASE RASTER


prism_raster <- raster('https://geoserver.usanpn.org/geoserver/wcs?service=WCS&version=2.0.1&request=GetCoverage&coverageId=precip:buffelgrass_rshiny&format=geotiff')
#prism_raster <- raster('nimbus.tif')
# Assign values as factors for mapping
prism_raster <- setValues(prism_raster, as.factor(getValues(prism_raster)))
# Assign appropriate extents (temp fix)
prism_raster@extent@xmin <- -114.9792
prism_raster@extent@xmax <- -108.9792
prism_raster@extent@ymin <- 31.3125
prism_raster@extent@ymax <- 37.02083
# Get a shapefile for AZ borders and crop prism_raster to AZ borders
AZ <- readRDS('./gadm36_USA_1_sp.rds') %>% subset(NAME_1=="Arizona")
prism_raster <- prism_raster %>% mask(AZ)


### RCC-ACIS (NOAA) STATION DATA


# Download JSON to working directory and open
url <- paste0('http://data.rcc-acis.org/MultiStnData?state=AZ&sdate=',
              as_date(lookback[1]),'&edate=',as_date(today),'&elems=pcpn')
dest <- paste0(getwd(),'/station_data.json') 
download.file(url, dest)
station_data <- fromJSON('station_data.json')

# Format df_stations data frame
df_stations <- data.frame( # Initialize an empty df
  longitude=replicate(length(station_data$data$meta$ll),NA),
  latitude=replicate(length(station_data$data$meta$ll),NA))
df_stations$name <- station_data[["data"]][["meta"]][["name"]]
df_stations <- mutate(df_stations, coordinates = station_data$data$meta$ll) # Populate df_stations$coordinates with vectors in form c(long, lat)
df_stations <- mutate(df_stations, longitude = lapply(df_stations$coordinates, FUN = function(coords) {return(coords[1])})) # Separate longitude
df_stations <- mutate(df_stations, latitude = lapply(df_stations$coordinates, FUN = function(coords) {return(coords[2])}))  # Separate latitude
station_values <- lapply(station_data[["data"]][["data"]], FUN = wBufferStations) # Calculate likelihood values
df_stations <- mutate(df_stations, values = station_values) # Assign these values to df
df_stations <- filter(df_stations, longitude != 'NULL') # Remove null objects -- they refer to areas such as "Greater Tucson Area"
df_stations$longitude <- unlist(df_stations$longitude) # Convert longitude from list of length 1 to numeric
df_stations$latitude <- unlist(df_stations$latitude) # Convert latitude from list of length 1 to numeric
df_stations <- subset(df_stations, select=-coordinates) # Remove coordinates in c(long, lat) form
df_stations$values <- unlist(df_stations$values) # Convert values from single-item list to numeric
stations <- df_stations
stations[["name_mod"]] <- paste(as.character(stations[["values"]]),stations[["name"]],sep = ' : ')


### RAINLOG


# gagues stores locational information for each gauge, readings stores precip data
# Combining them allows for forecasting and mapping of station by unique gaugeId
# Initialize Readings array w/ first 1000 (API only delivers 1000 readings at once)
readings <- getRainLogReadings(lookback[1], lookback[(length(lookback))], "0")
i <- 1000
# Loop and iterate, calling API function for each 1000 entries until no more available entries
# (because array length %% is not 1000)
while (length(readings[[1]]) %% 1000 == 0) {
  temp_readings <- getRainLogReadings(lookback[1], lookback[(length(lookback))], as.character(i))
  i = i + 1000
  readings <- rbind(readings, temp_readings) }
# Initialize gauges array and loop until complete
gauges <- getRainLogGauges(lookback[1], lookback[(length(lookback))], "0")
# Get vectorized columns for latitude and longitude and drop data frame column to avoid rbind ERROR in 208
gauges$latitude <- gauges$position$lat
gauges$longitude <- gauges$position$lng
gauges <- gauges[, -7]

i <- 1000
while (length(gauges[[1]]) %% 1000 == 0) {
  temp_gauges <- getRainLogGauges(lookback[1], lookback[(length(lookback))], as.character(i))
  temp_gauges$latitude <- temp_gauges$position$lat
  temp_gauges$longitude <- temp_gauges$position$lng
  temp_gauges <- temp_gauges[, -7]
  i = i + 1000
  gauges <- rbind(gauges, temp_gauges) }

# Arrange both gauges and readings arrays by gaugeId
readings <- arrange(readings, gaugeId)
gauges <- arrange(gauges, gaugeId)
# Prep gauges for use
gauges1 <- group_by(gauges, gaugeId) %>%
  summarize(latitude = mean(latitude), longitude = mean(longitude))

# Calculate forecast values for rainlog stations
rainlog <- group_by(readings, gaugeId) %>%
  summarize(forecast_value = rainlogPrep(rainAmount, readingDate)) %>%
  mutate(source = 'RainLog')

freelist <- replicate(length(rainlog$gaugeId), NA)
rllat <- freelist
rllng <- freelist
rllat[which(rainlog$gaugeId %in% gauges1$gaugeId)] <- gauges1$latitude[which(gauges1$gaugeId %in% rainlog$gaugeId)]
rllng[which(rainlog$gaugeId %in% gauges1$gaugeId)] <- gauges1$longitude[which(gauges1$gaugeId %in% rainlog$gaugeId)]
rainlog <- mutate(rainlog, lat=rllat, lng = rllng)
rainlog[["name_mod"]] <- paste(as.character(rainlog[["forecast_value"]]), "RainLog User", sep=" : ")



# Set up colors needed for mapping
map_colors <- colorFactor(palette=c('hotpink','#fbf3de','#bae4b3','#74c476','#238b41', '#105e1e'), domain=c(1,2,3,4,5,6), na.color=NA)
legend_colors <- colorFactor(palette=c('#fbf3de','#bae4b3','#74c476','#238b41', '#105e1e'), domain=c(0,1,2,3,4), na.color=NA)

# Create leaflet
l <- leaflet() %>%
  addProviderTiles('Esri') %>%
  addRasterImage(prism_raster, colors = map_colors, opacity=.8, project=T) %>% #Overplot the forecast
  addLegend(pal= legend_colors, values=c(0,1,2,3,4), opacity = .8, title='# of Events') %>% # Include legend 
  addCircles(data=stations, lng= ~longitude, lat= ~latitude, # Overplot station data
             stroke = T, color='black', weight = 1, radius = 250,
             fillColor = ~legend_colors(values), fillOpacity = .8, 
             popup = ~name_mod) %>% # Allow to click on points and bring up station names
  addCircles(data=rainlog, lng= ~lng, lat= ~lat, # Overplot station data
             stroke = T, color='black', weight = 1, radius = 250,
             fillColor = ~legend_colors(forecast_value), fillOpacity = .8, 
             popup = ~name_mod)
l
#end <- now()
#print(end-start)



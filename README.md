# Event-based Buffelgrass Forecast
This repository contains the R code used to generate the USA-NPN's prototyped updated buffelgrass green-up forecast.

## Overview
The R shiny app contained within this repo is designed to generate a forecast for the green-up of buffelgrass in Arizona. This prototype calculates the number of distinct precipitation events in a 30-day window. In contrast, the existing [NPN buffelgrass forecast](https://usanpn.org/data/forecasts/Buffelgrass) uses cumulative preciptation in a 24-day window. Precipitation data used is from [PRISM](http://prism.oregonstate.edu/) (for Rasters) and [RCC-ACIS.](https://www.rcc-acis.org/index.html) (for spatial point data). 

## Rules
This app counts rainfall events. Within it, an event is defined as a calendar day with 0.25+ inches of precipitation. The event extends until the first succeeding day with no precipitation. After an event, a buffer of at least 3 days with precipitation below 0.25 inches must occur in order for the next event to count. For clarification, a buffer will always start on a day with no precipitation, and will include days with no or little (less than 0.25in) precipitation up until the next day with sufficient precipitation.

## Access
An in-browser instance of the Shiny App can be run at [https://travismatlock.shinyapps.io/buffelgrass_v2/](https://travismatlock.shinyapps.io/buffelgrass_v2/). The app is limited to 25 active hours per month.

## Recent Developments
In March 2024, improvements were made to the first version event-based forecast. These are:
1. A nightly script downloads and preprocesses the PRISM precipitation raster data for improved efficiency.
2. When clicking on points, the name of the source and the event value are displayed.

## Future Developments
There are a number of possible improvements that are being considered for this app. They are as follows:
1. Slider tools to allow for user-controlled timeframes, precipitation threshold, and buffer length.
2. Integration of spatial point data into nightly scripts for further improved efficiency

## Author
* **Travis Matlock** - *Coding* - Student software developer at [USA-NPN](https://github.com/usa-npn).

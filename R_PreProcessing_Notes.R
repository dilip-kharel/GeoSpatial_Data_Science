# libraries used in PreProcessing field data


library(foreach) # Provides Foreach Looping Construct for R
library(gstat) # Spatial and Spatio-Temporal Geostatistical Modelling, Prediction and Simulation
library(rgeos) # Interface to Geometry Engine - Open Source ('GEOS')
library(sp) # Classes and Methods for Spatial Data
library(raster) # Geographic Data Analysis and Modeling
library(rgdal) # Bindings for the 'Geospatial' Data Abstraction Library
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(RColorBrewer) # Provides color schemes for maps (and other graphics) designed by Cynthia Brewer as described at http://colorbrewer2.org


# This example works on Farm data.

poly <- readOGR(dsn = 'path/to/Boundary_File',
                layer = 'Boundary_File_Name')

head(poly)


# standard column structure for AgLeader Advanced text format data set
# that is exported from SMS

agldr.col <- c("Longitude", "Latitude", "Flow", "GPStime", "Interval", "Distance",
               "Swath",  "Moisture", "Header", "PassNum", "SerialNum", "FieldID",
               "LoadID", "GrainType", "GPSstatus", "PDOP", "Altitude")

# Set working directory where files are Corn silage 2013

setwd("path/to/Working_Directory")

cs13 <- list.files(pattern = "_2013_.*\\.txt$", recursive = TRUE)


ssf13.silage <- foreach (i = cs13, .combine = 'rbind') %do% {
  
  df <- read.table(i, sep =",", fill=TRUE, header = FALSE) 
  colnames(df) <- agldr.col
  df$FieldName <- unlist(strsplit(i, "\\.")[[1]])[1]
  
  # Extract field ID using regexp ----------------
  begin <- regexpr("ssf_", df$FieldName, TRUE)
  end <- regexpr("_2013_", df$FieldName, TRUE)
  field <- substr(df$FieldName, begin+4 , end-1)
  df$Field <- field
  print(cbind(dim(df), unique(df$Field)))
  return(df)
}

unique(ssf13.silage$Field)
head(ssf13.silage)

# prepare dataset of each years as shown in above 
#----------------------------------------------------------------
#----------------------------------------------------------------
# This function reads both yield monitor dataset and boundary file,
# and extract correct field id for each dataset from boundary file

# Determine what columns to use from field boundary file
head(poly)
unique(poly$Field)




# function start
getFiledId <- function(dataset, fieldBoundary){
  # add Lat Long column to retain original column intact
  dataset$lat <-dataset$Latitude
  dataset$long <- dataset$Longitude
  poly <- fieldBoundary
  
  # Convert to same projection
  wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  nad <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  utm18 <- "+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
  
  data.sp <- dataset
  coordinates(data.sp)<- ~ long +lat
  proj4string(data.sp ) <- wgs
  data.sp <- spTransform(data.sp, CRS(utm18))
  
  proj4string(poly) <- wgs
  poly.utm <- spTransform(poly, CRS(utm18))
  
  test1 <- over(data.sp, poly.utm)
  data.sp$field2 <- test1$Field # Column for field info
  summary(data.sp$field2)
  
  # Remove dataset with NA in field2
  dataa.sp <- data.sp[complete.cases(data.sp$field2), ]
  datab.sp <- data.sp[!complete.cases(data.sp$field2), ]
  
  
  unique(dataa.sp$field2)
  unique(datab.sp$Field)
  
  datab.sp <- foreach (i = unique(datab.sp$Field), .combine = 'rbind') %do% {
    field.ss <- subset(datab.sp, Field == i)
    field.ss$nrow <- nrow(field.ss)
    plot(field.ss, main = paste(i, 'No Boundary field'))
    print(paste(unique(field.ss$nrow), i ))
    return(field.ss)
  }
  # Several fields are there with complete dataset,
  # but the boundary file is not available
  # work with these files later
  
  dataa.sp$field2 <- droplevels(dataa.sp$field2)
  
  # Calculate area based on points and polygon and find fraction of area 
  # covered by points, if fraction area covered by points are too small,
  # remove those field for further processing.
  
  # find the mode of the dataset
  (ModeDist <- which.max(tabulate(dataa.sp$Distance)))
  (ModeInterval <- which.max(tabulate(dataa.sp$Interval)))
  (ModeSwath <- which.max(tabulate(dataa.sp$Swath)))
  (area_pt <- (ModeDist * ModeSwath)/134)
  (acre_pt <- area_pt/43560)
  
  
  # use as.character to correct error in follwoing for loop
  # Error in { : task 1 failed - "level sets of factors are different"
  poly.utm$Field <- as.character(poly.utm$Field) # look for the field info columns
  dataa.sp2 <- foreach (i = unique(dataa.sp$field2), .combine = 'rbind') %do% {
    field.ss <- subset(dataa.sp, field2 == i)
    field.ss$nrow <- nrow(field.ss)
    a <- subset(poly.utm, Field == i)
    # sqm <- area(a)
    # 1 acre = 4047 sq m (43560 sq ft)
    acre <- area(a)/4047
    field.ss$acre_poly <- acre
    return(field.ss)
  }
  
  dataa.sp2$acre_point <- dataa.sp2$nrow * acre_pt
  dataa.sp2$area_frac <- dataa.sp2$acre_point/dataa.sp2$acre_poly
  unique(dataa.sp2$area_frac)
  
  # dataa.sp3 <- subset(dataa.sp2, area_frac > 0.70)
  # length(unique(dataa.sp3$field2))
  # unique(dataa.sp3$nrow)
  # 
  
  dataa.sp3 <- subset(dataa.sp2, area_frac > 0.50)
  #--------
  
  # plot and see each field first
  
  foreach (i = unique(dataa.sp3$field2)) %do% {
    field <- subset(poly.utm, Field == i)
    a <- subset(dataa.sp3, field2 == i)
    plot(field, main = paste(i, "With Boundary") )
    points(a, col = 'red')
    plot(field, lwd = 2, add = TRUE)
  }
  
  
  data.df <- as.data.frame(dataa.sp3)
  return(data.df)
}

# function end
#---------------------------------------------------
#---------------------------------------------------
# use function with each years dataset


ssf13silage.df <- getFiledId(ssf13.silage, poly)

# Fields with no boundary . barn side Indian Rd_33.

unique(ssf13silage.df$field2)
# keep only needed column info
colnames(ssf13silage.df)
ssf13silage.df <- droplevels(ssf13silage.df[, -c(21:26)])


#set working directory to a new folder where csv ouput can be saved
# These CSV file is in AgLeader Format, but contains Column informaiton
setwd("path/to/Output_Directory")
foreach (i = unique(ssf13silage.df$field2)) %do% {
  field <- subset(ssf13silage.df, field2 == i)
  write.csv(field, file = paste('SsF13silage',i, 'rfYE.csv', sep = '_'), 
            row.names = FALSE)
}



# Remove column name, save in a folder
# These files are ready to run in Yield Editor
setwd("path/to/Output_Directory")
foreach (i = unique(ssf13silage.df$field2)) %do% {
  field <- subset(ssf13silage.df, field2 == i)
  write.table(field, file = paste('SsF13silage_',i,'_rfYE.csv'),  
              row.names = FALSE, col.names = FALSE, sep = ",")
}

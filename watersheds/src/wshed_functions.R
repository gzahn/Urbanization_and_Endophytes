#
# Watershed functions
# Olson M.
#
library(sf);library(terra);library(whitebox);library(parallel);library(dplyr)

#=========================================================
# Data Preprocessing
#=========================================================
#
# clean & convert table into lat/lon vector features
# average GPS coordinates for each site
# Projection: NAD83 / UTM zone 12N
# prepare elevation data for study area (source: )
# prepare land cover data for study area (source: )


#=========================================================
# Watershed Functions
#=========================================================
#
# # # # # # # # # # # # # # # # # # # # # # # 
## Function. Flow accumulation
wbt_flow_accum <- function(dempath, outpath, breach_dist=5, fill_opt=TRUE, stream_thresh=6000){
  # 
  # Citation: Gannon, J. P. (2023), VT-Hydroinformatics, https://github.com/VT-Hydroinformatics
  # Wrapper tool to determine flow direction in DEM grid
  # files are saved to outpath folder
  #
  # 1. breach and fill holes in DEM
  whitebox::wbt_breach_depressions_least_cost(
    dem = dempath,
    output = file.path(outpath,"breached.tif"),
    dist = breach_dist, 
    fill = fill_opt)
  whitebox::wbt_fill_depressions_wang_and_liu(
    dem = file.path(outpath,"breached.tif"),
    output = file.path(outpath,"breached_filled.tif")
  )
  # 2. create flow accumulation and pointer
  whitebox::wbt_d8_flow_accumulation(input = file.path(outpath,"breached_filled.tif"),
                           output = file.path(outpath,"D8FA.tif"))
  whitebox::wbt_d8_pointer(dem = file.path(outpath,"breached_filled.tif"),
                 output = file.path(outpath,"D8pointer.tif"))
  # 3. generate streams
  whitebox::wbt_extract_streams(flow_accum = file.path(outpath,"D8FA.tif"),
                                output = file.path(outpath,"raster_streams.tif"),
                                threshold = stream_thresh)
  
  print(cat("\n Flow accumulation and streams saved in:",
            outpath, "\n") )
}

# # # # # # # # # # # # # # # # # # # # # # # 
## Function. Watershed delineate
wbt_delineate <- function(outpath, shpath, name.ext="", snap_distance=0.5){
  #
  # Wrapper tool for watershed delineation
  # files output to outpath
  # careful with units for snap_dist, based on DEM units
  #
  # check path
  if(file.exists(file.path(outpath,"raster_streams.tif"))){
    outpath_source <- outpath
  }else{ outpath_source <- dirname(outpath) } # go up one level}
  # 4. snap pour points to stream
  whitebox::wbt_jenson_snap_pour_points(pour_pts =  shpath,
                              streams = file.path(outpath_source,"raster_streams.tif"),
                              output = file.path(outpath,paste0("snappedpp_",name.ext,".shp")),
                              snap_dist = snap_distance) 
  # 5. delineate watersheds
  whitebox::wbt_watershed(d8_pntr = file.path(outpath_source,"D8pointer.tif"),
                pour_pts = file.path(outpath,paste0("snappedpp_",name.ext,".shp")),
                output = file.path(outpath,paste0("watershed_",name.ext,".tif")))
  
  print(cat("\n Delineation saved in:",
            file.path(outpath,"watershed.tif"), "\n") )
  
}

# # # # # # # # # # # # # # # # # # # # # # # 
## Function. Loop for watershed tools
wbt_watersheds <- function(dempath, shpath, outpath, huc8shed=NULL, longitudinal=TRUE, merge.all=TRUE){
  #
  # Runs watershed delineation for each point individually or all at once
  #
  strt = Sys.time()
  # generate flowpath grid
  wbt_flow_accum(dempath, outpath)
  
  # longitudinal
  if(!longitudinal){
    # run altogether - keep default settings
    wbt_delineate(outpath, shpath)
    stop("!subsequent steps need updating")
    # polygonize
  }else{
    # run point-by-point - tmp directories
    # save new idx pts to tmpath
    shpath_ix <- long_pts_parallel(shpath, outpath)
    # delineate individually
    wbt_delineate_long_paralell(shpath_ix, outpath)
    # force align
    wbt_watershed_forced(outpath)
    # polygonize
    watershed_poly_pts_long(outpath, huc8shed)
  } 
  # extract nlcd layers
  nlcd_watershed(outpath, nlcdpath)
  #
  print("Finished")
  print(Sys.time() - strt)
}

# # # # # # # # # # # # # # # # # # # # # # # 
## Function. separate points
long_pts_parallel <- function(shpath, outpath){
  #
  # save shp points to individual temporary file
  #   requires site_ID field
  #
  pts = st_read(shpath)
  tmpath = file.path(outpath,"tmppts")
  if(!dir.exists(tmpath)){dir.create(tmpath)}
  if(is.null(pts$site_ID)){
    stop("!shp needs site_ID field")
  }
  num_cores <- detectCores() - 1
  pt_indices <- splitIndices(nrow(pts), num_cores)
  
  export_pts <- function(indices) {
    for(i in indices){
      st_write(pts[i,,drop=F], file.path(tmpath, paste0("spoint_",i,".shp")))  
    }
  }
  parallel::mclapply(pt_indices, export_pts, mc.cores=num_cores)
  return(tmpath)
}

# # # # # # # # # # # # # # # # # # # # # # # 
## Function. Watershed delineate
wbt_delineate_long_paralell <- function(shpathx, outpath){
  #
  # Paralleled individual watershed delineation
  #
  # Set up cluster
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores) 
  # 
  # Function to run in parallel
  run_watershed <- function(i, shpathx, outpath) {
    # Unique in/out paths 
    outpath_i <- file.path(outpath, paste0("out_",i) )
    dir.create(outpath_i)
    shpath_i <- file.path(shpathx, paste0("spoint_",i,".shp") )
    # Run functions
    wbt_delineate(outpath_i, shpath_i, name.ext = i)
    
  }
  # Set up input indices 
  indices <- seq_len(10)
  # Export functions  
  clusterExport(cl, "wbt_delineate")
  # Run in parallel
  parLapply(cl, indices, run_watershed, shpathx=shpathx, outpath=outpath)
  # Stop cluster
  stopCluster(cl)
}

# # # # # # # # # # # # # # # # # # # # # # # 
## Function. Forced alignment - needs update
wbt_watershed_forced <- function(outpath){
  #
  # Forcing alignment for points that are too far off stream
  # UPDATE eventually
  #
  # read in stream and snapped point data
  ppoints.all <- list.files(outpath, pattern="snappedpp_.+.shp",recursive = TRUE, full.names = TRUE)
  # re-order
  pnum <- as.numeric(gsub("\\D", "", basename(ppoints.all)))
  ppoints.all <- ppoints.all[match(sort(pnum),pnum)]
  
  # streams
  streams <- rast(file.path(outpath,"raster_streams.tif"))
  
  # extract cell values
  ext.cells <- unlist(lapply(ppoints.all, function(x) extract(streams, vect(x), cells=T)$cell))
  pcol <- unlist(lapply(ext.cells, function(x) rowColFromCell(streams, x )[1,2])) # extract col
  prow <- unlist(lapply(ext.cells, function(x) rowColFromCell(streams, x )[1,1])) # extract row
  # run ns to find stream (only runs up/down column)
  # st.cell <- which(streams[,pcol[1]]==1); pclose.row <- which.min(abs(st.cell- prow[1]))
  # st.cell[pclose.row]
  new.row <- unlist(lapply(1:length(prow), function(x) which(streams[,pcol[x]]==1)[which.min(abs(which(streams[,pcol[x]]==1) - prow[x]))]) )
  new.row[1]=new.row[1]-1
  # which to fix
  to.fix <- is.na(unlist(lapply(ppoints.all, function(x) extract(streams, vect(x))$raster_streams)))
  ppoints.all[to.fix]
  
  # rerun each individually
  for (i in which(to.fix)){
    # remove files
    list.files(dirname(ppoints.all[i]),pattern = "snappedpp",full.names = T) %>% unlink()
    # rewrite snapped points
    cellFromRowCol(streams,new.row[i],pcol[i]) %>% xyFromCell(streams,.) %>% as.data.frame() %>%
      st_as_sf(.,coords=c("x","y"), crs = crs(streams)) %>% st_write(.,ppoints.all[i] )
    # rerun wbt watershed
    # 5. delineate watersheds
    whitebox::wbt_watershed(d8_pntr = file.path(dirname(dirname(ppoints.all[i])),"D8pointer.tif"),
                            pour_pts = ppoints.all[i],
                            output = file.path(dirname(ppoints.all[i]),paste0("watershed_",i,".tif")))
    
    print(cat("\n forcing points:",which(to.fix),"\n until match \n"))
  }
}

# # # # # # # # # # # # # # # # # # # # # # # 
## Function. Watershed polygonize
watershed_poly_pts_long <-  function(outpath, huc8shedpath){
  #
  # merge all longitudinal polygons
  #
  # define outpath
  wspaths <- list.files(outpath, pattern="watershed",recursive = TRUE, full.names = TRUE)
  # re-order
  wsnum <- as.numeric(gsub("\\D", "", basename(wspaths)))
  wspaths <- wspaths[match(sort(wsnum),wsnum)]
  
  # lower provo watershed  (static)
  pvws <- st_read(huc8shedpath)
  # Set up cluster
  num_cores <- detectCores()-1;cl <- makeCluster(num_cores) 
  # 
  # Function to run in parallel
  watershed_poly <- function(i, outpath, wspaths, pvws) {
    # polygonize
    wsshape <- rast(wspaths[i]) %>% stars::st_as_stars() %>% st_as_sf(merge = T)
    # clip
    interpv = st_intersection(pvws, wsshape)
    # write
    st_write(interpv, file.path(outpath,paste0("out_",i),paste0("watershed_",i,".shp")))
  }
  # Set up input indices 
  indices <- seq_len(10)
  # Export functions  
  # clusterExport(cl, "dplyr")
  clusterEvalQ(cl,  {library(dplyr);library(sf);library(terra)})
  # Run in parallel
  parLapply(cl, indices, watershed_poly, outpath=outpath, wspaths=wspaths, pvws=pvws)
  # Stop cluster
  stopCluster(cl)
  
  # list new watersheds
  wshps <- list.files(outpath, pattern="watershed_.+.shp",recursive = TRUE, full.names = TRUE)
  wsnum <- as.numeric(gsub("\\D", "", basename(wshps)))
  wshps <- wshps[match(sort(wsnum),wsnum)]
  # merge
  single_sf <- dplyr::bind_rows(lapply(wshps, function(x) st_read(x)))
  wsh <- single_sf %>% mutate(site_ID=1:10) %>% select(site_ID)
  # save
  st_write(wsh,file.path(outpath, "watersheds.shp"))
  
}

# # # # # # # # # # # # # # # # # # # # # # # 
## Function. NLCD extract
nlcd_watershed <- function(outpath, nlcdpath){
  #
  # calculate impervious surface area by water shed
  #
  #
  # read in watersheds
  # wsh <- st_read(file.path(outpath, "watersheds.shp"))
  wsh <- vect(file.path(outpath, "watersheds.shp"))
  # new watershed out
  new_out <- file.path(outpath,"land_sheds")
  if(!file.exists(new_out)){dir.create(new_out)}
  #
  
  ### START HERE
  # - save nlcd then make table - upload
  nlcdesc = rast(file.path(nlcdpath, "nlcd_2021_impervious_desc_proj.tif")) 
  nlcd = rast(file.path(nlcdpath, "nlcd_2021_landcover_proj.tif")) 
  
  for (i in 1:10){
    nlcd_ix = nlcdesc %>% crop(wsh[i,]) %>% mask(wsh[i,])
    nlcdl_ix = nlcd %>% crop(wsh[i,]) %>% mask(wsh[i,])
    # write
    writeRaster(nlcd_ix, file.path(new_out, paste0("nlcd_",i,"_desc.tif")) )
    writeRaster(nlcdl_ix, file.path(new_out, paste0("nlcd_",i,"_lc.tif")) )
  }
}

# # # # # # # # # # # # # # # # # # # # # # # 
## Function. Final table
extract_lc_table <- function(outpath){
  # 
  # extract variables for each site
  # ws_area, nlcd2021, imper2021
  #
  wsh <- st_read(file.path(outpath, "watersheds.shp"))
  wsh.tbl <- wsh %>% as.data.frame()
  # outpath
  new_out <- file.path(outpath,"land_sheds")
  
  # iterate through points
  for (i in 1:nrow(wsh.tbl)){
    print(paste("...extracting data from", i, "of", nrow(wsh.tbl)))
    # extract info from layers
    interpv <- wsh[i,] 
    v.area <- interpv %>% st_area() %>% as.numeric() / 1000^2 # in km^2
    if(i==1){wsh.tbl$area_km2 = NA}
    wsh.tbl$area_km2[i] = v.area
    
    # extract land cover
    nl21 <- rast( file.path(new_out, paste0("nlcd_",i,"_lc.tif")) )
    nl21i <- rast( file.path(new_out, paste0("nlcd_",i,"_desc.tif")) )
     
    # summarize NLCD values
    nltbl = nl21 %>% as.data.frame() %>% na.omit() %>%
      group_by(`NLCD Land Cover Class`) %>% count() %>% mutate(n = (n*res(nl21)[1]^2)/1000^2 )
    # summarize Impervious values
    imtbl = nl21i %>% as.data.frame() %>% na.omit() %>%
      group_by(Class_Names) %>% count() %>% mutate(n = (n*res(nl21i)[1]^2)/1000^2 )

    # save table
    if (i==1){
      df.nl = nltbl
      df.i = imtbl
      
    }else{
      # df.nl[i, 7:ncol(df.nl)] <- nltbl
      # df.i[i, 7:ncol(df.i)] <- imtbl
      # account for any new classes
      df.nl <- full_join(df.nl,nltbl,by="NLCD Land Cover Class")
      df.i <- full_join(df.i,imtbl,by= "Class_Names")
    }
  }
  # # rename and join tables
  colnames(df.nl)[2:11] <- paste0("site_",1:10)
  colnames(df.i)[2:11] <- paste0("site_",1:10)
  
  # pivot
  dfnl2 = df.nl %>% t() %>% as.data.frame()
  colnames(dfnl2) = dfnl2[1,];dfnl2 <- dfnl2[2:nrow(dfnl2),];dfnl2[is.na(dfnl2)] <- 0.0
  dfi2 = df.i %>% t() %>% as.data.frame()
  colnames(dfi2) = dfi2[1,];dfi2 <- dfi2[2:nrow(dfi2),];dfi2[is.na(dfi2)] <- 0.0
  
  # bind
  dfinal.nlcd = cbind(wsh.tbl, dfnl2) %>% select(-geometry)
  dfinal.imperv = cbind(wsh.tbl, dfi2)%>% select(-geometry)
  
  # save tables
  write.csv(dfinal.nlcd,file.path(outpath, "tables", "nlcd_landcover.csv"), row.names=FALSE)
  write.csv(dfinal.imperv,file.path(outpath, "tables", "nlcd_impervious_desc.csv"), row.names=FALSE)
  
  print("Complete")
  
}

# # # # # # # # # # # # # # # # # # # # # # # 
## Function. Clean and generate sf object from csv
geospat_csv <- function(csvpath, reproj.crs=CRS("EPSG:26912"), write.path=NULL){
  #
  # Create sf object and match projection - SPECIFIC to sample-sites.csv
  #
  # read in and clean data
  dfp = read.csv(csvpath)
  tail(dfp)
  dfp = dfp %>% na.omit %>% select(site_ID, lat, lon, collect_date, process_date, sra_description)
  dfp2 = dfp %>% select (site_ID, lat, lon)
  # take average of GPS points & merge
  new_idx = as.numeric(gsub("\\D", "", dfp2$site_ID))
  dfp2$site_ID = new_idx
  dfp3 = dfp2 %>% group_by(site_ID) %>% summarize(lat = mean(lat), lon = mean(lon))
  dfp3 = cbind(dfp3,dfp[match(unique(new_idx),new_idx),4:6])
  # create sf object
  shp <- st_as_sf(x = dfp3,                         
                  coords = c("lon", "lat"),
                  crs = CRS("EPSG:4326")) 
  if(!is.null(reproj)){shp <- st_transform(shp, reproj.crs)}
  if(!is.null(write.path)){ st_write(shp, write.path, overwrite = TRUE)}
  return(shp)
}
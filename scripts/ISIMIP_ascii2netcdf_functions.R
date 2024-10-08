###################################################################
# title:    4C-ISIMIP_ascii2netcdf_functions
# purpose:  provide functions and file naming conventions for the
#           main script "4C-ISIMIP_ascii2netcdf"
# note:     foreign function writeSim2netCDF() from Friedrich Bohn
#           adapted from ProfoundData package
# author:   Mats Mahnken (PIK Potsdam) - mahnken@pik-potsdam.de
# date:     10.10.2019
###################################################################

### procedures for netcdf creation

#' @title A function to write netCDF-files
#' @description This function transforms simulation results into netCDF files following the ISIMIP2 protocol
#'
#' @param df A data.frame containing in the first three columns longitude latitude
#'  and time. These columns are followed by columns containing the output variables.
#'   The columns have to be named with the output variable name as required by the 2B
#'   protocol. See table 21.
#' @param modelname The name of the used forest model
#' @param GCM The climate model which created the used climate time series
#' @param RCP The RCP scenario
#' @param ses The scenario describing forest management. UMsoc equals the "nat"
#' settings and histsoc and 2005soc equal the "man" settings in the ISIMIP2a
#' protocol. Default value: "nat".
#' @param ss  "co2" for all experiments other than the sensitivity experiments
#' for which 2005co2 is explicitly written.  Note: even models in which CO2 has
#' no effect should use the co2 identifier relevant to the experiment.  Default
#' value: "co2const".
#' @param start the start year of the simulation. Default value: 1980.
#' @param timestep the temporal resolution of the variables in df. 
#' Note: must be either "annual", "monthly" or "daily" and follow %Y, %Y-%m
#' or %Y-%m-%d respectively.Default value: "annual".
#' @param region the region or site of the simulation
#' @param folder The folder in which all netCDF files will be written
#' @param contact Your mail address
#' @param institution Your institution
#' @param comment1 Optional comment regarding your simulation
#' @param comment2 Optional comment regarding your simulation
#'
#' @details
#' The function transforms your simulation output data frame into several netCDF
#' -files and writes them into the indicated folder using the naming convention of
#'  the ISIMIP2(B)-protocol (https://www.isimip.org/protocol/). Units and long names
#'  of variables (table 21) will be created automatically.
#' @note To report errors in the package or the data, please use the issue tracker
#' in the GitHub repository of ProfoundData \url{https://github.com/COST-FP1304-PROFOUND/ProfoundData}
#' @example /inst/examples/writeSim2netCDFHelp.R
#' @export
#' @author Friedrich J. Bohn

writeSim2netCDF<-function(df,
                          comment1=NA,
                          comment2=NA,
                          institution= 'PIK',
                          contact= 'isi-mip@pik-potsdam.de',
                          modelname="4c",
                          GCM="hadgem",
                          bc="localbc",
                          RCP="rcp26",
                          ses ="nat",
                          ss="co2const",
                          region="peitz",
                          start='1980',
                          timestep='annual',
                          soildepths=NA,
                          folder="ISIMIP",
			  ISIMIP_round="2A"){
  
  if(length(comment1)==1) {
    ifelse(is.na(comment1), comment1<-rep("",length(colnames(df))), comment1<-rep(comment1,length(colnames(df))))
  }
  if(length(comment2)==1) {
    ifelse(is.na(comment2), comment2<-rep("",length(colnames(df))), comment2<-rep(comment2,length(colnames(df))))
  }
  
  for( i in 4:length(colnames(df))){
    variable<-colnames(df)[i]
    code<-strsplit(variable,"_")[[1]]
    # swith code part 1
    {switch(code[1],
            dbh={
              unit<-"cm"
              variable_long<-"Mean DBH"
            },
            dbhdomhei={
              unit<-"cm"
              variable_long<-"Mean DBH of 100 highest trees"
            },
            hei={
              unit<-"m"
              variable_long<-"Stand Height"
            },
            height={
              unit<-"m"
              variable_long<-"Stand Height"
            },
            domhei={
              unit<-"m"
              variable_long<-"Dominant Height"
            },
            domheight={
              unit<-"m"
              variable_long<-"Dominant Height"
            },
            density={
              unit<-"ha-1"
              variable_long<-"Stand Density"
            },
            ba={
              unit<-"m2 ha-1"
              variable_long<-"Basal Area"
            },
            mort={
              unit<-"m3 ha-1"
              variable_long<-"Volume of Dead Trees"
            },
            harv={
              unit<-"m3 ha-1"
              variable_long<-"Harvest by dbh-class"
            },
            stemno={
              unit<-"ha-1"
              variable_long<-"Remaining stem number after disturbance and management by dbh class"
            },
            vol={
              unit<-"m3 ha-1"
              variable_long<-"Stand Volume"
            },
            cveg={
              unit<-"kg m-2"
              variable_long<-"Carbon Mass in Vegetation biomass"
            },
            cvegag={
              unit<-"kg m-2"
              variable_long<-"Carbon Mass in aboveground vegetation biomass"
            },
            cvegbg={
              unit<-"kg m-2"
              variable_long<-"Carbon Mass in belowground vegetation biomass"
            },
            clitter={
              unit<-"kg m-2"
              variable_long<-"Carbon Mass in Litter Pool"
            },
            csoil={
              unit<-"kg m-2"
              variable_long<-"Carbon Mass in Soil Pool per soil layer"
            },
            age={
              unit<-"yr"
              variable_long<-"Tree age by dbh class"
            },
            gpp={
              unit<-"kg m-2 s-1"
              variable_long<-"Gross Primary Production"
            },
            npp={
              unit<-"kg m-2 s-1"
              variable_long<-"Net Primary Production"
            },
            ra={
              unit<-"kg m-2 s-1"
              variable_long<-"Autotrophic (Plant) Respiration"
            },
            rh={
              unit<-"kg m-2 s-1"
              variable_long<-"Heterotrophic Respiration"
            },
            nee={
              unit<-"kg m-2 s-1"
              variable_long<-"Net Ecosystem Exchange"
            },
            mai={
              unit<-"m3 ha-1"
              variable_long<-"Mean Annual Increment"
            },
            fapar={
              unit<-"%"
              variable_long<-"Fraction of absorbed photosynthetically active radiation"
            },
            lai={
              unit<-"m2 m-2"
              variable_long<-"Leaf Area Index"
            },
            species={
              unit<-"%"
              variable_long<-"Species composition"
            },
            evap={
              unit<-"kg m-2 s-1"
              variable_long<-"Total Evapotranspiration "
            },
            intercept={
              unit<-"kg m-2 s-1"
              variable_long<-"Evaporation from Canopy (interception) "
            },
            esoil={
              unit<-"kg m-2 s-1"
              variable_long<-"Water Evaporation from Soil"
            },
            trans={
              unit<-"kg m-2 s-1"
              variable_long<-"Transpiration"
            },
            soilmoist={
              unit<-"kg m-2"
              variable_long<-"Soil Moisture per soil layer"
            },
            mortstemno={
              unit<-"ha-1"
              variable_long<-"Removed stem numbers by size class by natural mortality"
            },
            harvstemno={
              unit<-"ha-1"
              variable_long<-"Removed stem numbers by size class by natural management"
            },
            dist={
              unit<-"m3 ha-1"
              variable_long<-"Volume of disturbance damage"
            },
            nlit={
              unit<-"g m-2 a-1"
              variable_long<-"Nitrogen of annual Litter"
            },
            nsoil={
              unit<-"g m-2 a-1"
              variable_long<-"Nitrogen in Soil"
            },
            nppleaf={
              unit<-"kg m-2 s-1"
              variable_long<-"Net Primary Production allocated to leaf biomass"
            },
            npproot={
              unit<-"kg m-2 s-1"
              variable_long<-"Net Primary Production allocated to fine root biomass"
            },
            nppagwood={
              unit<-"kg m-2 s-1"
              variable_long<-"Net Primary Production allocated to above ground wood biomass"
            },
            nppbgwood={
              unit<-"kg m-2 s-1"
              variable_long<-"Net Primary Production allocated to below ground wood biomass"
            },
            rr={
              unit<-"kg m-2 s-1"
              variable_long<-"Root autotrophic respiration"
            },
            cleaf={
              unit<-"kg m-2"
              variable_long<-"Carbon Mass in Leaves"
            },
            cwood={
              unit<-"kg m-2"
              variable_long<-"Carbon Mass in Wood"
            },
            croot={
              unit<-"kg m-2"
              variable_long<-"Carbon Mass in Roots"
            },
            tsl={
              unit<-"K"
              variable_long<-"Temperature of Soil per soil layer"
            },
            {
              message(paste("error:colname",variable,"does not correspond to the ISI-MIP convention of the ISI.MIP simulation protocol 2.1"))
              return(FALSE)
            }
    )
    }
    # switch code part 2
    if (length(code)>1){
      switch(code[2],
             total={variable_long<-paste(variable_long,"Total")},
             fasy={variable_long<-paste(variable_long,"Fagus sylvatica")},
             quro={variable_long<-paste(variable_long,"Quercus robur")},
             qupe={variable_long<-paste(variable_long,"Quercus petraea")},
             pisy={variable_long<-paste(variable_long,"Pinus sylvestris")},
             piab={variable_long<-paste(variable_long,"Picea abies")},
             pipi={variable_long<-paste(variable_long,"Pinus pinaster")},
             lade={variable_long<-paste(variable_long,"Larix decidua")},
             eugl={variable_long<-paste(variable_long,"Eucalyptus globulus")},
             acpl={variable_long<-paste(variable_long,"Acer platanoides ")},
             bepe={variable_long<-paste(variable_long,"Betula pendula")},
             frex={variable_long<-paste(variable_long,"Fraxinus excelsior")},
             poni={variable_long<-paste(variable_long,"Populus nigra")},
             rops={variable_long<-paste(variable_long,"Robinia pseudoacacia")},
             hawo={variable_long<-paste(variable_long,"Fagus sylvatica")},
             psme={variable_long<-paste(variable_long,"Pseudotsuda menziesii")},
             domhei={},
             height={},
             {
               message(paste("error:colname",variable,"does not correspond to the ISI-MIP convention of the ISI.MIP simulation protocol 2.1"))
               return(FALSE)
             }
      )
    }
    
    # write netCDF file
    write.netCDF(df=df
                 ,variable=variable
                 ,variable_long=variable_long
                 ,unit= unit
                 ,timestep=timestep
                 ,comment1=comment1[i]
                 ,comment2=comment2[i]
                 ,institution=institution
                 ,contact= contact
                 ,modelname=modelname
                 ,GCM=GCM
                 ,bc=bc
                 ,RCP=RCP
                 ,ses =ses
                 ,ss=ss
                 ,region=region
                 ,soildepths = soildepths
                 ,folder=folder)
  }
  return (TRUE)
}


write.netCDF<-function(df,
                       variable,
                       variable_long=  'precipitation',
                       unit= 'kg m-2 s-1',
                       timestep='annual',
                       comment1=  "",
                       comment2="",
                       institution= 'PIK',
                       contact= 'isi-mip@pik-potsdam.de',
                       modelname="4c",
                       GCM="hadgem",
                       bc="localbc",
                       RCP="rcp26",
                       ses ="nat",
                       ss="co2const",
                       region="peitz",
                       soildepths=NA,
                       folder="ISIMIP") {
  # cheques
  if (colnames(df)[1] != "lon" || colnames(df)[2] != "lat" || colnames(df)[3] != "time") {
    message('Error: First 3 columns must be named "lon", "lat", "time".')
    return(FALSE)
  }
  if (length(unique(df$lon))>1||length(unique(df$lat))>1){
    message("Error: more than one coordinate in data")
    return(FALSE)
  }
  if (!(timestep %in% c("annual", "monthly", "daily"))) {
    message("Error: timestep must be either one of \"annual\", \"monthly\" or \"daily\"")
    return(FALSE)
  }
  
  # timestep unit construction, start, end years and conversion into days/months/years since 1661-01-01 00:00:00
  if(timestep == "annual") {
    time_unit <- "years since 1661-01-01 00:00:00"
    start<-df$time[1]
    end<-df$time[nrow(df)]
    df$time <- as.numeric(as.character(df$time))-1661
  }
  if(timestep == "monthly") {
    time_unit <- "months since 1661-01-01 00:00:00"
    start <- strsplit(as.character(df$time), "-")[[1]][1]
    end <- strsplit(as.character(df$time), "-")[[length(df$time)]][1]
    df$time <- 12 * as.numeric(as.character(gsub("\\-.*","",df$time))) - 1661 + as.numeric(as.character(gsub(".*-","",df$time))) - 1}
  if(timestep == "daily") {
    time_unit <- "days since 1661-01-01 00:00:00"
    start<-format(as.Date(df$time), "%Y")[1]
    end<-format(as.Date(df$time), "%Y")[nrow(df)]
    df$time <- as.numeric(difftime(as.Date(df$time), "1661-01-01"))
  }
  
  # filename construction
  if(ISIMIP_round == "3A") {
    filename<-paste(modelname,GCM,RCP,ses,ss,paste(gsub("_", "-", variable), region, sep="-"),timestep,start,end,sep="_")
    title<-paste(modelname,GCM,RCP,ses,ss,region)
  } else {
    filename<-paste(modelname,GCM,bc,RCP,ses,ss,gsub("_", "-", variable),region,timestep,start,end,sep="_")
    title<-paste(modelname,bc,GCM,RCP,ses,ss,region)
  }
  filename<-paste0(folder,"/",filename,".nc4")

  
  #create file
  ncout <- RNetCDF::create.nc(filename, format="classic4")
  
  # determine if variable is dbh or depth classed
  is_dbhclassed <- F
  if(grepl("harv", variable) || grepl("stemno", variable) || grepl("age", variable)) {is_dbhclassed <- T}
  
  is_depthclassed <- F
  if(grepl("csoil", variable) || grepl("soilmoist", variable) || grepl("tsl", variable)) {is_depthclassed <- T}
  
  # dimension definitions
  if(is_dbhclassed) {RNetCDF::dim.def.nc(ncout, "dbh_class", 29)}
  if(is_depthclassed) {RNetCDF::dim.def.nc(ncout, "depth", length(soildepths))}
  RNetCDF::dim.def.nc(ncout, "lon", 1)
  RNetCDF::dim.def.nc(ncout, "lat", 1)
  RNetCDF::dim.def.nc(ncout, "time", unlim=TRUE)
  
  # variables
  if(is_dbhclassed) {
    RNetCDF::var.def.nc(ncout, "dbh_class", "NC_FLOAT", "dbh_class")
    RNetCDF::att.put.nc(ncout, "dbh_class", "long_name", "NC_CHAR", "dbh_class")
    RNetCDF::att.put.nc(ncout, "dbh_class", "standard_name", "NC_CHAR", "dbh_class")
    RNetCDF::att.put.nc(ncout, "dbh_class", "units", "NC_CHAR", "cm")
    RNetCDF::att.put.nc(ncout, "dbh_class", "axis", "NC_CHAR", "Z")
  }
  
  if(is_depthclassed) {
    RNetCDF::var.def.nc(ncout, "depth", "NC_FLOAT", "depth")
    RNetCDF::att.put.nc(ncout, "depth", "long_name", "NC_CHAR", "depth_below_land")
    RNetCDF::att.put.nc(ncout, "depth", "standard_name", "NC_CHAR", "depth")
    RNetCDF::att.put.nc(ncout, "depth", "positive", "NC_CHAR", "down")
    RNetCDF::att.put.nc(ncout, "depth", "units", "NC_CHAR", "m")
    RNetCDF::att.put.nc(ncout, "depth", "axis", "NC_CHAR", "Z")
  }
  
  RNetCDF::var.def.nc(ncout, "lon", "NC_FLOAT", "lon")
  RNetCDF::att.put.nc(ncout, "lon", "long_name", "NC_CHAR", "longitude")
  RNetCDF::att.put.nc(ncout, "lon", "standard_name", "NC_CHAR", "longitude")
  RNetCDF::att.put.nc(ncout, "lon", "units", "NC_CHAR", "degrees_east")
  RNetCDF::att.put.nc(ncout, "lon", "axis", "NC_CHAR", "X")
  
  RNetCDF::var.def.nc(ncout, "lat", "NC_INT", "lat")
  RNetCDF::att.put.nc(ncout, "lat", "long_name", "NC_CHAR", "latitude")
  RNetCDF::att.put.nc(ncout, "lat", "standard_name", "NC_CHAR", "latitude")
  RNetCDF::att.put.nc(ncout, "lat", "units", "NC_CHAR", "degrees_north")
  RNetCDF::att.put.nc(ncout, "lat", "axis", "NC_CHAR", "Y")
  
  RNetCDF::var.def.nc(ncout, "time", "NC_INT", "time")
  RNetCDF::att.put.nc(ncout, "time", "long_name", "NC_CHAR", "time")
  RNetCDF::att.put.nc(ncout, "time", "units", "NC_CHAR", time_unit)
  RNetCDF::att.put.nc(ncout, "time", "calendar", "NC_CHAR", "proleptic_gregorian")
  RNetCDF::att.put.nc(ncout, "time", "axis", "NC_CHAR", "T")
  
  RNetCDF::var.def.nc(ncout, gsub("_", "-", variable), "NC_FLOAT", if(is_dbhclassed || is_depthclassed) {c(0:3)} else {c(0:2)})
  RNetCDF::att.put.nc(ncout, gsub("_", "-", variable), "_FillValue", "NC_FLOAT", 1.e+20)
  RNetCDF::att.put.nc(ncout, gsub("_", "-", variable), "missing_value", "NC_FLOAT", 1.e+20)
  RNetCDF::att.put.nc(ncout, gsub("_", "-", variable), "short_field_name", "NC_CHAR", gsub("_", "-", variable))
  RNetCDF::att.put.nc(ncout, gsub("_", "-", variable), "long_field_name", "NC_CHAR", variable_long)
  RNetCDF::att.put.nc(ncout, gsub("_", "-", variable), "units", "NC_CHAR", unit)
  
  RNetCDF::att.put.nc(ncout, "NC_GLOBAL", "title", "NC_CHAR", title)
  RNetCDF::att.put.nc(ncout, "NC_GLOBAL", "comment1", "NC_CHAR", comment1)
  RNetCDF::att.put.nc(ncout, "NC_GLOBAL", "comment2", "NC_CHAR", comment2)
  RNetCDF::att.put.nc(ncout, "NC_GLOBAL", "institute", "NC_CHAR", institution)
  RNetCDF::att.put.nc(ncout, "NC_GLOBAL", "contact", "NC_CHAR", contact)
  
  # data
  if(is_dbhclassed) {
    RNetCDF::var.put.nc(ncout,"dbh_class",seq(0, 140, 5),1,length(seq(0, 140, 5)))
    RNetCDF::var.put.nc(ncout,"lon",df$lon[1],1,1)
    RNetCDF::var.put.nc(ncout,"lat",df$lat[1],1,1)
    RNetCDF::var.put.nc(ncout,"time",df$time,1,length(df$time))
    collumn<-which(colnames(df)==variable)
    RNetCDF::var.put.nc(ncout,gsub("_", "-", variable), unlist(df[,collumn]),c(1,1,1,1),c(length(seq(0, 140, 5)),1,1,nrow(df)))
  } 
  
  if(is_depthclassed) {
    RNetCDF::var.put.nc(ncout,"depth",soildepths,1,length(soildepths))
    RNetCDF::var.put.nc(ncout,"lon",df$lon[1],1,1)
    RNetCDF::var.put.nc(ncout,"lat",df$lat[1],1,1)
    RNetCDF::var.put.nc(ncout,"time",df$time,1,length(df$time))
    collumn<-which(colnames(df)==variable)
    RNetCDF::var.put.nc(ncout,gsub("_", "-", variable), unlist(df[,collumn]),c(1,1,1,1),c(length(soildepths),1,1,nrow(df)))
  } 
  
  if(!is_dbhclassed && !is_depthclassed) {
    RNetCDF::var.put.nc(ncout,"lon",df$lon[1],1,1)
    RNetCDF::var.put.nc(ncout,"lat",df$lat[1],1,1)
    RNetCDF::var.put.nc(ncout,"time",df$time,1,length(df$time))
    collumn<-which(colnames(df)==variable)
    RNetCDF::var.put.nc(ncout,gsub("_", "-", variable),df[,collumn],c(1,1,1),c(1,1,nrow(df)))
  }
  
  
  RNetCDF::close.nc(ncout)
}
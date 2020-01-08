###########################################################################
# title:    ISIMIP2B-FORMASAM_descriptive-plots
# purpose:  producing descriptive plots of simulation results 
#           from ISIMIP2B-FORMASAM alternative management scenario runs
# note:     model output must be provided in ISIMIP2B netCDF standard in
#           the "indir" directory; reading and aggregating data from 
#           the netCDFs might take some time (several minutes)
#           depending on number of experiments simulated; variable
#           specifications can be read in the current ISIMIP2B protocol
# author:   Mats Mahnken, Potsdam-Institute for Climate Impact Research
# date:     07.01.2020
###########################################################################

# load required R-packages
library('ggplot2')
library('gridExtra')
library('cowplot')
library('ncdf4')
library('RNetCDF')


########################## settings
##########################

# input directory with model simulation outputs as netCDF following ISIMIP2B conventions on file naming and structure
indir <- "C:/Users/..."
# output directory where the graphics should be stored
outdir <- "C:/Users/..."

# enable ("warn=0") or disable ("warn=-1") warnings globally if you get too many incorrect warnings
options(warn=0)

# enable ("png") or disable ("pdf") producing single .png plots instead of aggregated .pdf
# !!! "png" option produces large amounts of individual .png files in outdir
outformat <- "pdf"

# define variables to be read (do not change)
varlist <- c("ba", "npp", "cveg-", "nee", "vol", "csoil", "harv-", "gpp")


########################## read and aggregate model output; aggregated data is stored in one RData-file per site in indir
########################## including all variables in varlist and from all simulated ISIMIP2B experiments

# read all file names in input directory except *.RData files
allfiles <- grep(list.files(indir, recursive = TRUE), pattern='.RData', inv=T, value=T)

# produce a summary of the simulation experiments available
sumfiles <- do.call(rbind,sapply(1:length(allfiles),function(i) {strsplit(allfiles[i],"\\/|\\_|\\.")}))
sumfiles <- sumfiles[,-12]
sumfiles <- as.data.frame(sumfiles)
colnames(sumfiles) <- c("modelname",
                        "gcm",
                        "bias_correction",
                        "climate_scenario",
                        "management",
                        "CO2sens_scenario",
                        "variable",
                        "region",
                        "timestep",
                        "start",
                        "end")

# read available sites
sitelist <- unique(sumfiles$region)

experiments <- NULL

# read data from netCDFs
for(site in sitelist) { # iterate over sites
  data_array <- data.frame(gcm=factor(), 
                           climate_scenario=factor(), 
                           co2sens_scenario=factor(), 
                           management=factor(), 
                           year=integer())
  new_experiments <- unique(sumfiles[which(sumfiles$region==site),c(2,4,5,6,8)])
  experiments <- rbind(experiments, new_experiments)
  for(var in varlist) { # iterate over variables
    data_array_var <- NULL
    var_sp_list <- unique(sumfiles$variable[which(sumfiles$region==site & grepl(var, sumfiles$variable))])
    for(var_sp in var_sp_list) {  # iterate over species/total files of variable
      data_array_var_sp <- NULL
      site_var_files <- which(sumfiles$region==site & sumfiles$variable==var_sp)
      for(file in site_var_files) { # iterate over all experiment
        
        # establish connection to netCDF
        ncin <- nc_open(paste(indir, "/",allfiles[file], sep=""))
        
        # extract data array from netCDF
        var_array <- ncvar_get(ncin,names(ncin$var))
        # extract time array from netCDF
        time_array <- ncvar_get(ncin,"time")
        tunits <- ncatt_get(ncin,"time",attname="units")
        tustr <- strsplit(tunits$value, " ")
        date_array <- as.Date(time_array,origin=unlist(tustr)[3])
        year_array <- as.numeric(format.Date(strptime(date_array,"%Y-%m-%d"),"%Y"))
        # extract unit from netCDF
        var_unit_array <- ncatt_get(ncin,names(ncin$var),attname="units")
        var_unit <- var_unit_array$value
        
        # close connection to netCDF
        nc_close(ncin)
        
        # unclass depth-classes of csoil and dbh-classes of harv
        if(length(dim(var_array)) > 1) {
          var_array <- colSums(var_array)
        }
        
        # read start and end years of simulation
        start_year <- as.integer(as.character(sumfiles$start[file])); end_year <- as.integer(as.character(sumfiles$end[file]))
        if(start_year!=year_array[1] || end_year!=year_array[length(year_array)]) {warning(paste0("start/end year in time dimension of netCDF does not match file name for ", allfiles[file]))}
        
        # aggregate daily data to annual data as the mean (supposed daily values represent fluxes)
        if(sumfiles$timestep[file]=="daily") {
          var_array <- as.vector(aggregate(var_array, list(year_array), mean)[,2])
        }
        
        # append var_array + scenario_info + year to data_array
        new_data <- data.frame(gcm=as.character(sumfiles$gcm[file]),
                               climate_scenario=as.character(sumfiles$climate_scenario[file]),
                               co2sens_scenario=as.character(sumfiles$CO2sens_scenario[file]), 
                               management=as.character(sumfiles$management[file]),
                               year=start_year:end_year,
                               var=var_array)
        colnames(new_data) <- c(colnames(new_data)[-length(colnames(new_data))], paste0(as.character(var_sp), " [", var_unit, "]"))
        data_array_var_sp <- rbind(data_array_var_sp, new_data)
      }
      # merge all variables for all experiments
      data_array <- merge(data_array, data_array_var_sp, all=T, by=c("gcm", "climate_scenario", "co2sens_scenario", "management", "year"))
    }
  }
  
  # harmonize management scenario names across climate-scenarios
  data_array$management <- as.character(data_array$management)
  data_array$management[which(grepl("socbe", data_array$management))] <- "socbe"
  data_array$management[which(grepl("sochwp", data_array$management))] <- "sochwp"
  data_array$management[which(grepl("soca", data_array$management))] <- "soca"
  data_array$management[which(grepl("socam", data_array$management))] <- "socam"
  
  # store aggregated data (data_array) in one .RData-file per site in indir
  save(data_array, file=paste0(indir, "/data-array_", site, ".RData"))
  
  # print progress to console
  print(paste0("finished aggregating data from ", site, ": found ", length(which(sumfiles$region==site)), " files from ", length(new_experiments[,1]), " simulation experiments"))
}


########################## plot descriptive graphics to outdir to get overview of model output
##########################

# define own color palette for different simulation experiment specifiers
cols <- c("picontrol.co2" = "grey", "rcp26.co2" = "green", "rcp60.co2" = "firebrick1", "rcp26.2005co2" = "seagreen" ,"rcp60.2005co2" = "firebrick",
          "picontrol" = "grey", "rcp26" = "green4", "rcp60" = "firebrick4")

# define own linetypes and shapes for different species
lty_sp <- c("total"=1, "fasy"=2, "quro"=3, "qupe"=3, "pisy"=4, "piab"=5, "pipi"=4, "lade"=6, "acpl"=2, "eugl"=6, "bepe"=5, "bepu"=5, "rops"=6, "frex"=6, "poni"=6, "soau"=6, "c3gr"=6, "hawo"=6, "psme"=6)
shp_sp <- c("total"=15, "fasy"=16, "quro"=17, "qupe"=17, "pisy"=18, "piab"=1, "pipi"=2, "lade"=3, "acpl"=4, "eugl"=5, "bepe"=6, "bepu"=6, "rops"=8, "frex"=9, "poni"=10, "soau"=11, "c3gr"=12, "hawo"=13, "psme"=7)


########################## plotting full time series

for(site in sitelist) { # iterate over all sites
  
  # read aggregated data from RData-file and open pdf
  load(paste0(indir, "/data-array_", site, ".RData"))
  if(outformat=="pdf") {pdf(paste0(outdir, "/timeseries_full_", site, ".pdf"),paper = "a4r", width = 0, height = 0)}
  
  for(var in varlist) { # iterate over variables
    
    # select all species/total specific variables linked to var
    vars <- colnames(data_array[6:length(colnames(data_array))])
    vars <- vars[grepl(var, vars)]
    
    tmp_var_data_array <- NULL
    
    # plotting if multiple species shall be displayed
    if(length(vars) > 1) {
      for(var_sp in vars) {  # iterate over all variables
        var_id <- which(colnames(data_array)==var_sp)
        tmp_var_sp_data_array <- data_array[, c(1:5, var_id)]
        colnames(tmp_var_sp_data_array)[6] <- var
        tmp_var_sp_data_array$species <- strsplit(var_sp, "-| ")[[1]][2]
        tmp_var_data_array <- rbind(tmp_var_data_array, tmp_var_sp_data_array)
      }
      
      # read variable unit for plotting yaxis title
      ylab_unit <- paste0("[", strsplit(var_sp, "[", fixed=T)[[1]][2])
      
      # remove all rows with NAs as variable values
      tmp_var_data_array <- tmp_var_data_array[which(!is.na(tmp_var_data_array[,6])), ]
      
      # define all simulated management scenarios
      mans <- unique(tmp_var_data_array$management)
      
      for(man in mans)  {   # iterate over all management scenarios
        tmp_var_man_data_array <- tmp_var_data_array[which(tmp_var_data_array$management==man),]
        
        if(all(is.na(tmp_var_man_data_array[,6]))) next   # skip loop if data is empty
        
        # open png
        if(outformat=="png") png(paste0(outdir, "/timeseries_full_", site, "_", var, "_", man, ".png"),width = 3600, height = 2000, res=200)
        
        # plotting
        p <- ggplot(tmp_var_man_data_array, aes(year, tmp_var_man_data_array[,6], group=factor(interaction(gcm, climate_scenario, co2sens_scenario, species)), linetype=gcm, shape=species, color=factor(interaction(climate_scenario, co2sens_scenario)))) +
          ylim(min(tmp_var_man_data_array[,6]), max(tmp_var_man_data_array[,6])) +
          scale_fill_manual(values = cols) +
          scale_color_manual(values = cols) +
          scale_shape_manual(values = shp_sp) +
          geom_line() +
          geom_point(size=1.5) +
          labs(title=paste(man, "management,", site), x="year", y=paste0(gsub("-", "", var), " ", ylab_unit)) +
          theme(legend.title=element_blank(), plot.title=element_text(size=18), axis.title=element_text(size=18), axis.text=element_text(size=18), legend.text = element_text(size = 18), legend.position = "right")
        
        print(p)
        
        # close png
        if(outformat=="png") dev.off()
      }
    }
    
    # plotting if only total is to be displayed
    if(length(vars) == 1) {
      for(var in vars) {  # iterate over all variables
        var_id <- which(colnames(data_array)==var)
        tmp_var_data_array <- data_array[, c(1:5, var_id)]
        
        # define all simulated management scenarios
        mans <- unique(tmp_var_data_array$management)
        
        for(man in mans)  {   # iterate over all management scenarios
          tmp_var_man_data_array <- tmp_var_data_array[which(tmp_var_data_array$management==man),]
          
          if(all(is.na(tmp_var_man_data_array[,6]))) next   # skip loop if data is empty
          
          # open png
          if(outformat=="png") png(paste0(outdir, "/timeseries_full_", site, "_", var, "_", man, ".png"),width = 3600, height = 2000, res=200)
          
          # plotting
          p <- ggplot(tmp_var_man_data_array, aes(year, tmp_var_man_data_array[,6], group=factor(interaction(gcm, climate_scenario, co2sens_scenario)), linetype=gcm, color=factor(interaction(climate_scenario, co2sens_scenario)))) +
            ylim(min(tmp_var_man_data_array[,6]), max(tmp_var_man_data_array[,6])) +
            scale_fill_manual(values = cols) +
            scale_color_manual(values = cols) +
            scale_shape_manual(values = shp_sp) +
            geom_line() +
            geom_point(size=1.5) +
            labs(title=paste(man, "management,", site), x="year", y=paste(var)) +
            theme(legend.title=element_blank(), plot.title=element_text(size=18), axis.title=element_text(size=18), axis.text=element_text(size=18), legend.text = element_text(size = 18), legend.position = "right")
          
          print(p)
          
          # close png
          if(outformat=="png") dev.off()
        }
      }
    }
  }
  
  if(outformat=="pdf") dev.off()  # close pdf
}


########################## plotting aggregated time series (aggregation over GCMs)

for(site in sitelist) { # iterate over all sites
  
  # read aggregated data from RData-file and open pdf
  load(paste0(indir, "/data-array_", site, ".RData"))
  vars <- colnames(data_array[6:length(colnames(data_array))])
  
  # open pdf
  if(outformat=="pdf") {pdf(paste0(outdir, "/timeseries_aggr_", site, ".pdf"),paper = "a4r", width = 0, height = 0)}
  
  for(var in varlist) { # iterate over variables
    
    # select all species/total specific variables linked to var
    vars <- colnames(data_array[6:length(colnames(data_array))])
    vars <- vars[grepl(var, vars)]
    
    tmp_var_data_array <- NULL
    
    # plotting if multiple species shall be displayed
    if(length(vars) > 1) {
      for(var_sp in vars) {  # iterate over all variables
        var_id <- which(colnames(data_array)==var_sp)
        tmp_var_sp_data_array <- data_array[, c(1:5, var_id)]
        colnames(tmp_var_sp_data_array)[6] <- var
        tmp_var_sp_data_array$species <- strsplit(var_sp, "-| ")[[1]][2]
        tmp_var_data_array <- rbind(tmp_var_data_array, tmp_var_sp_data_array)
      }
      # read variable unit for plotting yaxis title
      ylab_unit <- paste0("[", strsplit(var_sp, "[", fixed=T)[[1]][2])
      
      # remove all rows with NAs as variable values
      tmp_var_data_array <- tmp_var_data_array[which(!is.na(tmp_var_data_array[,6])), ]
      
      # define all simulated management scenarios
      mans <- unique(tmp_var_data_array$management)
      
      for(man in mans)  {   # iterate over all management scenarios
        tmp_var_man_data_array <- tmp_var_data_array[which(tmp_var_data_array$management==man),]
        
        if(all(is.na(tmp_var_man_data_array[,6]))) next   # skip loop if data is empty
        
        # open png
        if(outformat=="png") png(paste0(outdir, "/timeseries_full_", site, "_", var, "_", man, ".png"),width = 3600, height = 2000, res=200)
        
        # plotting
        p <- ggplot(tmp_var_man_data_array, aes(year, tmp_var_man_data_array[,6], group=factor(interaction(climate_scenario, co2sens_scenario, management, species)), linetype=species, color=factor(interaction(climate_scenario, co2sens_scenario)))) +
          ylim(min(tmp_var_man_data_array[,6]), max(tmp_var_man_data_array[,6])) +
          scale_color_manual(values = cols) +
          scale_linetype_manual(values=lty_sp) +
          stat_summary(geom = "line", fun.y = mean, na.rm=T, show.legend=T, size=1.2, alpha=0.7) +
          labs(x="", y=paste0(gsub("-", "", var), " ", ylab_unit)) +
          theme(legend.title=element_blank(), legend.text = element_text(size = 6), legend.position = "right") +
          ggtitle(paste0(man))
        
        print(p)
        
        # close png
        if(outformat=="png") dev.off()
      }
    }
    
    # plotting if only total is to be displayed
    if(length(vars) == 1) {
      for(var in vars) {  # iterate over all variables
        var_id <- which(colnames(data_array)==var)
        tmp_var_data_array <- data_array[, c(1:5, var_id)]
        mans <- unique(tmp_var_data_array$management)
        for(man in mans)  {   # iterate over all management scenarios
          tmp_var_man_data_array <- tmp_var_data_array[which(tmp_var_data_array$management==man),]
          
          if(all(is.na(tmp_var_man_data_array[,6]))) next   # next loop if data is empty
          
          # define all simulated management scenarios
          if(outformat=="png") png(paste0(outdir, "/timeseries_full_", site, "_", var, "_", man, ".png"),width = 3600, height = 2000, res=200)
          
          # plotting
          p <- ggplot(tmp_var_man_data_array, aes(year, tmp_var_man_data_array[,6], group=factor(interaction(climate_scenario, co2sens_scenario, management)), color=factor(interaction(climate_scenario, co2sens_scenario)))) +
            ylim(min(tmp_var_man_data_array[,6]), max(tmp_var_man_data_array[,6])) +
            scale_color_manual(values = cols) +
            scale_linetype_manual(values=lty_sp) +
            stat_summary(geom = "line", fun.y = mean, na.rm=T, show.legend=T, size=1.2, alpha=0.7) +
            labs(x="", y=paste0(var)) +
            theme(legend.title=element_blank(), legend.text = element_text(size = 6), legend.position = "right") +
            ggtitle(paste0(man))
          
          print(p)
          
          # close png
          if(outformat=="png") dev.off()
        }
      }
    }
  }
  if(outformat=="pdf") dev.off()  # close pdf
}


########################## summary plotting of total harvest, soil c accumulation and nee

# open pdf
if(outformat=="pdf") {pdf(paste0(outdir, "/summary.pdf"),paper = "a4r", width = 0, height = 0)}

for(site in sitelist) { # iterate over all sites
  
  # read data from RData file
  load(paste0(indir, "/data-array_", site, ".RData"))
  
  # select all respective variables
  var_id_h <- NA
  var_id_h <- which(grepl("harv-total", colnames(data_array)))
  if(length(var_id_h) == 0) {var_id_h <- which(grepl("harv", colnames(data_array)))}
  var_id_c <- NA
  var_id_c <- which(grepl("csoil-total", colnames(data_array)))
  if(length(var_id_c) == 0) {var_id_c <- which(grepl("csoil", colnames(data_array)))}
  var_id_n <- NA
  var_id_n <- which(grepl("nee-total", colnames(data_array)))
  if(length(var_id_n) == 0) {var_id_n <- which(grepl("nee", colnames(data_array)))}
  var_id <- c(var_id_h, var_id_c, var_id_n)
  tmp_var_data_array <- data_array[, c(1:5, var_id)]
  
  # aggregate harvest as sum over all simulation years
  tmp_aggregated_var_data_array <- aggregate(tmp_var_data_array[,6], list(tmp_var_data_array$gcm, tmp_var_data_array$climate_scenario, tmp_var_data_array$co2sens_scenario, tmp_var_data_array$management), sum)
  # aggregate csoil as the difference between csoil at simulation start and mean annual csoil (carbon stored in soil)
  tmp_aggregated_var_data_array <-  merge(tmp_aggregated_var_data_array, aggregate(tmp_var_data_array[,7], list(tmp_var_data_array$gcm, tmp_var_data_array$climate_scenario, tmp_var_data_array$co2sens_scenario, tmp_var_data_array$management), function(x) ((sum(x)/length(x))-x[1])), by=names(tmp_aggregated_var_data_array)[-5])
  # aggregate nee as mean nee of all simulation years
  tmp_aggregated_var_data_array <-  merge(tmp_aggregated_var_data_array, aggregate(tmp_var_data_array[,8], list(tmp_var_data_array$gcm, tmp_var_data_array$climate_scenario, tmp_var_data_array$co2sens_scenario, tmp_var_data_array$management), mean), by=names(tmp_aggregated_var_data_array)[c(-5,-6)])
  
  colnames(tmp_aggregated_var_data_array) <- c("gcm", "climate_scenario", "co2sens_scenario", "management", "harv", "csoil_accu", "nee_sum")
  
  # open png
  if(outformat=="png") png(paste0(outdir, "/summary_", site, ".png", sep=""), width = 2000, height = 3000, res=300)
  
  # plot harv figure
  p1 <- ggplot(tmp_aggregated_var_data_array, aes(management, harv, color=factor(interaction(climate_scenario, co2sens_scenario)))) +
    geom_point(position = position_jitterdodge(jitter.width=0), aes(fill=factor(interaction(climate_scenario, co2sens_scenario))), size=2, show.legend = F) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    ylab(paste0("sum of ", colnames(data_array)[var_id_h])) +
    xlab("") +
    labs(fill='') +
    theme(legend.title = element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle = 90)) +
    guides(fill=F) +
    ggtitle(paste0(site))
  
  # plot csoil figure
  p2 <- ggplot(tmp_aggregated_var_data_array, aes(management, csoil_accu, color=factor(interaction(climate_scenario, co2sens_scenario)))) +
    geom_point(position = position_jitterdodge(jitter.width=0), aes(fill=factor(interaction(climate_scenario, co2sens_scenario))), size=2, show.legend = F) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    ylab(paste0("mean ann. accumulation of ", colnames(data_array)[var_id_c])) +
    xlab("") +
    labs(fill='') +
    theme(legend.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(angle = 90)) +
    guides(fill=F)
  
  # plot nee figure
  p3 <- ggplot(tmp_aggregated_var_data_array, aes(management, nee_sum, color=factor(interaction(climate_scenario, co2sens_scenario)))) +
    geom_point(position = position_jitterdodge(jitter.width=0), aes(fill=factor(interaction(climate_scenario, co2sens_scenario))), size=2, show.legend = F) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    ylab(paste0("mean ", colnames(data_array)[var_id_n])) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE), breaks = scales::pretty_breaks(n = 3)) +
    xlab("") +
    labs(fill='') +
    theme(legend.title = element_blank(), axis.text.y = element_text(angle = 90)) +
    guides(fill=F)
  
  grid.arrange(p1, p2, p3, nrow=3)
  
  # close png
  if(outformat=="png") dev.off()
}
if(outformat=="pdf") dev.off()  # close pdf


########################## plotting GCM variable output

# open pdf
if(outformat=="pdf") {pdf(paste0(outdir, "/gcm-range.pdf"),paper = "a4r", width = 0, height = 0)}

# define variables to compare (do not change)
compvars <- c("gpp", "npp", "nee")

for(var in compvars) {  # iterate over c-flux variables (in compvars)
  gcmrangeall <- NULL
  
  for(site in sitelist) { # iterate over all sites
    
    # load data from RData file
    load(paste0(indir, "/data-array_", site, ".RData"))
    
    # extract variable
    var_id <- NA
    var_id <- which(grepl(paste0(var, "-total"), colnames(data_array)))
    if(length(var_id) == 0) {var_id <- which(grepl(var, colnames(data_array)))}
    if(length(var_id) != 1) {warning("none or multiple variables named \"", var, "\" or \"", paste0(var, "-total"), "\" for at least one experiment")}
    tmp_var_data_array <- data_array[, c(1:5, var_id)]
    
    # compute the range from all gcm simulation for each individual experiment
    tmp_mean_annual_var_data_array <- aggregate(tmp_var_data_array[,6], list(tmp_var_data_array$gcm, tmp_var_data_array$climate_scenario, tmp_var_data_array$co2sens_scenario, tmp_var_data_array$management), mean)
    tmp_mean_annual_range_var_data_array <- aggregate(tmp_mean_annual_var_data_array[,5], list(tmp_mean_annual_var_data_array[,2], tmp_mean_annual_var_data_array[,3], tmp_mean_annual_var_data_array[,4]), range)
    tmp_mean_annual_range_var_data_array$gcmrange <- apply(tmp_mean_annual_range_var_data_array[,4], 1, max) - apply(tmp_mean_annual_range_var_data_array[,4], 1, min)
    
    tmp_mean_annual_range_var_data_array$site <- site
    
    gcmrangeall <- rbind(gcmrangeall, tmp_mean_annual_range_var_data_array)
  }
  
  gcmrangeall <- as.data.frame(gcmrangeall)
  colnames(gcmrangeall) <- c("climate_scenario", "co2sens_scenario", "management", "range", "gcmrange", "site")
  
  # open png
  if(outformat=="png") png(paste0(outdir, "/gcm-range_", var, ".png"), width = 3000, height = 1500, res=300)
  
  # plotting
  p <- ggplot(gcmrangeall, aes(site, gcmrange, group=factor(interaction(site, climate_scenario)), col=climate_scenario)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(jitter.width=0)) +
    scale_color_manual(values = cols) +
    ylab(paste0("range of output for linked runs with different GCMs in ", colnames(data_array)[var_id])) +
    xlab("") +
    theme(legend.title = element_blank()) +
    ggtitle("output variability introduced by GCM")
  
  print(p)
  
  # close png
  if(outformat=="png") dev.off()
}
if(outformat=="pdf") dev.off()  # close pdf


########################## plotting potential CO2 fertilization effect

# open pdf
if(outformat=="pdf") {pdf(paste0(outdir, "/co2diff.pdf"),paper = "a4r", width = 0, height = 0)}

for(var in compvars) {  # iterate over c-flux variables (in compvars)
  co2diffall <- NULL
  
  for(site in sitelist) { # iterate over all sites
    
    # load data from RData file
    load(paste0(indir, "/data-array_", site, ".RData"))
    
    # extract variable
    var_id <- NA
    var_id <- which(grepl(paste0(var, "-total"), colnames(data_array)))
    if(length(var_id) == 0) {var_id <- which(grepl(var, colnames(data_array)))}
    if(length(var_id) != 1) {warning("none or multiple variables named \"", var, "\" or \"", paste0(var, "-total"), "\" for at least one experiment")}
    tmp_var_data_array <- data_array[, c(1:5, var_id)]
    
    # compute the difference between co2 and 2005co2 simulation runs for each individual experiment
    tmp_mean_annual_var_data_array <- aggregate(tmp_var_data_array[,6], list(tmp_var_data_array$gcm, tmp_var_data_array$climate_scenario, tmp_var_data_array$co2sens_scenario, tmp_var_data_array$management), mean)
    tmp_mean_annual_var_data_array <- tmp_mean_annual_var_data_array[which(!grepl("pic", tmp_mean_annual_var_data_array[,2])),] # exclude all picontrol runs
    tmp_mean_annual_diff_var_data_array <- merge(tmp_mean_annual_var_data_array[which(tmp_mean_annual_var_data_array[,3]=="co2"),], tmp_mean_annual_var_data_array[which(tmp_mean_annual_var_data_array[,3]!="co2"),], by=names(tmp_mean_annual_var_data_array)[c(-3,-5)])
    tmp_mean_annual_diff_var_data_array$co2diff <- tmp_mean_annual_diff_var_data_array[,5] - tmp_mean_annual_diff_var_data_array[,7]
    
    tmp_mean_annual_diff_var_data_array$site <- site
    
    co2diffall <- rbind(co2diffall, tmp_mean_annual_diff_var_data_array)
  }
  
  co2diffall <- as.data.frame(co2diffall[,-c(4,5,6,7)])
  colnames(co2diffall) <- c("gcm", "climate_scenario", "management", "co2diff", "site")
  
  # open png
  if(outformat=="png") png(paste0(outdir, "/co2diff_", var, ".png", sep=""), width = 1200, height = 1500, res=300)
  
  # plotting
  p <- ggplot(co2diffall, aes(gcm, co2diff, col=climate_scenario)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(jitter.width=0)) +
    scale_color_manual(values = cols) +
    ylab(paste0("average difference between linked co2 and 2005co2 runs in ", colnames(data_array)[var_id])) +
    xlab("") +
    theme(legend.title = element_blank()) +
    ggtitle("potential CO2 fertilization effect")
  
  print(p)
  
  # close png
  if(outformat=="png") dev.off()
}
if(outformat=="pdf") dev.off()  # close pdf


########################## plotting potential climate change effect

# open pdf
if(outformat=="pdf") {pdf(paste0(outdir, "/cceffect.pdf"),paper = "a4r", width = 0, height = 0)}

for(var in compvars) {  # iterate over c-flux variables (in compvars)
  cceffectall <- NULL
  
  for(site in sitelist) { # iterate over all sites
    
    # load data from RData file
    load(paste0(indir, "/data-array_", site, ".RData"))
    
    # extract variable
    var_id <- NA
    var_id <- which(grepl(paste0(var, "-total"), colnames(data_array)))
    if(length(var_id) == 0) {var_id <- which(grepl(var, colnames(data_array)))}
    if(length(var_id) != 1) {warning("none or multiple variables named \"", var, "\" or \"", paste0(var, "-total"), "\" for at least one experiment")}
    tmp_var_data_array <- data_array[, c(1:5, var_id)]
    
    # compute the difference between picontrol and rcp simulation runs for each individual experiment
    tmp_mean_annual_var_data_array <- aggregate(tmp_var_data_array[,6], list(tmp_var_data_array$gcm, tmp_var_data_array$climate_scenario, tmp_var_data_array$co2sens_scenario, tmp_var_data_array$management), mean)
    tmp_mean_annual_diff_var_data_array <- merge(tmp_mean_annual_var_data_array[which(!grepl("pic", tmp_mean_annual_var_data_array[,2])),], tmp_mean_annual_var_data_array[which(grepl("pic", tmp_mean_annual_var_data_array[,2])),], by=names(tmp_mean_annual_var_data_array)[c(-2,-3,-5)], all.x=T)
    tmp_mean_annual_diff_var_data_array$cceffect <- tmp_mean_annual_diff_var_data_array[,5] - tmp_mean_annual_diff_var_data_array[,8]
    
    tmp_mean_annual_diff_var_data_array$site <- site
    
    cceffectall <- rbind(cceffectall, tmp_mean_annual_diff_var_data_array)
  }
  
  cceffectall <- as.data.frame(cceffectall[,-c(5,6,7,8)])
  colnames(cceffectall) <- c("gcm", "management","climate_scenario", "co2sens_scenario", "cceffect", "site")
  
  if(var=="nee") {ylim <- c(NA, 0)} else {ylim <- c(0, NA)}
  
  # open png
  if(outformat=="png") png(paste0(outdir, "/cceffect_", var, ".png", sep=""), width = 2000, height = 1500, res=300)
  
  # plotting
  p <- ggplot(cceffectall, aes(gcm, cceffect, group=factor(interaction(gcm, climate_scenario, co2sens_scenario)), col=factor(interaction(climate_scenario, co2sens_scenario)))) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(jitter.width=0)) +
    scale_color_manual(values = cols) +
    ylab(paste0("average difference between picontrol run and linked rcp run in ", colnames(data_array)[var_id])) +
    ylim(ylim) +
    xlab("") +
    theme(legend.title = element_blank()) +
    ggtitle("potential climate change effect")
  
  print(p)
  
  #close png
  if(outformat=="png") dev.off()
}
if(outformat=="pdf") dev.off()  # close pdf


###################################################################################
# title:    4C-ISIMIP_ascii2netcdf
# purpose:  4C output conversion from ascii to netcdf following 
#           ISIMIP conventions on file naming and file content
# author:   Mats Mahnken (PIK Potsdam) - mahnken@pik-potsdam.de
# date:     19.12.2019
###################################################################################

#install.packages('RNetCDF', lib="/p/tmp/mahnken/Rlibs")
library('RNetCDF', lib.loc="/home/mahnken/p-tmp-mahnken/Rlibs")

# set working directory
setwd("/home/mahnken/p-tmp-mahnken")

# run the functions script to include all necessary functions for the data conversion (at the moment only writeSim2netCDF)
source("ascii2netcdf/4C-ISIMIP_ascii2netcdf_functions.R")

# define ISIMIP round
ISIMIP_round <- "3A"

# input files need to be named and organized under 
# wd/01_ascii_input/<ISIMIP_round>/<bias_correction>/<region>_<gcm>_<soc-scenario>_<climate-scenario>_<co2sens-scenario>
# where experiment specifiers can differ from the ones used in ISIMIP conventions (but order in filename needs to be followed)
input <- (paste("ascii2netcdf/01_ascii_input/", ISIMIP_round, sep=""))
# output files will be written to "wd/02_netcdf_output/<ISIMIP_round>/"
output <- (paste("ascii2netcdf/02_netcdf_output/", ISIMIP_round, sep=""))


###################################################################################
### 01: ISIMIP + 4C file naming and variable definition
###################################################################################

## ISIMIP file naming convention follows:
## <modelname>_<gcm/obs>_<bias-correction>_<climate-scenario>_<soc�scenario>_<co2sens�scenarios>_<variable>_<region>_<timestep>_<start�year>_<end�year>.nc4

# scenario naming
gcm_2A_naming <- c("gswp3", "princeton", "watch", "wfdei", "localclim")
gcm_2B_naming <- c("hadgem2-es", "ipsl-cm5a-lr", "miroc5", "gfdl-esm2m")
gcm_3A_naming <- c("chelsa-w5e5")
bias_correction_2A_naming <- c("wfd", "nobc")
bias_correction_2B_naming <- c("nobc", "localbc", "ewembi", "ewembi-isimip3basd")
bias_correction_3A_naming <- c("nobc")
climate_scenario_2A_naming <- c("rcp2p6", "rcp4p5", "rcp6p0", "rcp8p5", "hist", "presclim", "noclim")
climate_scenario_2B_naming <- c("picontrol", "historical", "rcp26", "rcp60", "rcp85")
climate_scenario_3A_naming <- c("obsclim")
soc_scenario_2A_naming <- c("nosoc", "pressoc", "varsoc", "nat")
soc_scenario_2B_naming <- c("nosoc", "1860soc", "histsoc", "2005soc", "rcp26soc", "rcp60soc", "2100rcp26soc", "2005socsite", "socbe", "sochwp", "soca", "socam")
soc_scenario_3A_naming <- c("hist")
co2sens_scenario_2A_naming <- c("co2", "noco2", "co2const")
co2sens_scenario_2B_naming <- c("co2", "2005co2")
co2sens_scenario_3A_naming <- c("30arcsec", "90arcsec", "300arcsec", "1800arcsec")
region_2A_naming <- c("hyytiala", "peitz", "solling-beech", "solling-spruce", "soro", "kroof", "le-bray", "collelongo", "bily-kriz")  # actually should be "solling_beech", "solling_spruce", "le_bray" and "bily_kriz"
region_2B_naming <- c("hyytiala", "peitz", "solling-beech", "solling-spruce", "soro", "kroof", "le-bray", "collelongo", "bily-kriz")
region_3A_naming <- c("hyytiala", "peitz", "solling-beech", "solling-spruce", "soro", "kroof", "le-bray", "collelongo", "bily-kriz")


# variable specifier naming
species_2A_naming <- c("fasy", "quro", "qupe", "pisy", "piab", "pipi", "lade", "acpl", "eugl", "bepe", "bepu", "rops", "frex", "poni", "soau", "hawo")
species_2B_naming <- c("fasy", "quro", "qupe", "pisy", "piab", "pipi", "lade", "acpl", "eugl", "bepe", "bepu", "rops", "frex", "poni", "soau", "c3gr", "hawo", "psme")
species_3A_naming <- c("fasy", "quro", "qupe", "pisy", "piab", "pipi", "lade", "acpl", "eugl", "bepe", "bepu", "rops", "frex", "poni", "soau", "c3gr", "hawo", "psme")

dbhclass_naming <- c("dbh_c0", "dbh_c5", "dbh_c10", "dbh_c15", "dbh_c20", "dbh_c25", "dbh_c30",
                     "dbh_c35", "dbh_c40", "dbh_c45", "dbh_c50", "dbh_c55", "dbh_c60", "dbh_c65",
                     "dbh_c70", "dbh_c75", "dbh_c80", "dbh_c85", "dbh_c90", "dbh_c95", "dbh_c100",
                     "dbh_c105", "dbh_c110", "dbh_c115", "dbh_c120", "dbh_c125", "dbh_c130", "dbh_c135",
                     "dbh_c140")

timestep_2A <- c("monthly", "annual", "decadal", "daily")
timestep_2B <- c("3hr", "daily", "monthly", "annual")
timestep_3A <- c("3hr", "daily", "monthly", "annual")

dist_name <- c("fi", "wi", "ins", "dr", "graz", "dis")

# variable naming (for mandatory(mand) and optional(opt) variables)
variable_mand_2A_naming<- c("dbh-<species/total>", 
                            "dbhdomhei", # actually "dbh-domhei" is correct
                            "height-<species/total>", 
                            "domheight",  # actually "dom-height" is correct
                            "density-<species/total>", 
                            "ba-<species/total>", 
                            "mort-<species/total>", 
                            "harv-<species/total>-<dbhclass/total>", 
                            "stemno-<species/total>-<dbhclass/total>", 
                            "vol-<species/total>",
                            "cveg-<species/total>", 
                            "clitter-<species/total>", 
                            "csoil", # actually "csoil-<species/total>" is correct
                            "age-<species/total>-<dbhclass/total>", 
                            "gpp-<species/total>", 
                            "npp-<species/total>", 
                            "ra-<species/total>", 
                            "rh-<total>", 
                            "nee-<total>", 
                            "mai-<species/total>", 
                            "fapar-<species/total>", 
                            "lai-<species/total>", 
                            "species-<species>", 
                            "evap", # actually "evap-<total>" is correct
                            "intercept-<species/total>", 
                            "esoil", 
                            "trans-<species/total>", 
                            "soilmoist")
variable_mand_2B_naming <- c("dbh-<species/total>", 
                             "dbhdomhei", 
                             "hei-<species/total>", 
                             "domhei", 
                             "density-<species/total>", 
                             "ba-<species/total>", 
                             "mort-<species/total>", 
                             "harv-<species/total>-<dbhclass/total>", 
                             "stemno-<species/total>-<dbhclass/total>", 
                             "vol-<species/total>",
                             "cveg-<species/total>", 
                             "cvegag-<species/total>", 
                             "cvegbg-<species/total>", 
                             "clitter-<species/total>", 
                             "csoil", # actually "csoil-<species/total>" is correct
                             "age-<species/total>-<dbhclass/total>", 
                             "gpp-<species/total>", 
                             "npp-<species/total>", 
                             "ra-<species/total>", 
                             "rh-<total>", 
                             "nee-<total>", 
                             "mai-<species/total>", 
                             "fapar-<species/total>", 
                             "lai-<species/total>", 
                             "species-<species>", 
                             "evap", 
                             "intercept-<species/total>", 
                             "esoil", 
                             "trans-<species/total>", 
                             "soilmoist")

variable_mand_3A_naming <- c("dbh-<species/total>", 
                             "dbhdomhei", 
                             "hei-<species/total>", 
                             "domhei", 
                             "density-<species/total>", 
                             "ba-<species/total>", 
                             "mort-<species/total>", 
                             "harv-<species/total>-<dbhclass/total>", 
                             "stemno-<species/total>-<dbhclass/total>", 
                             "vol-<species/total>",
                             "cveg-<species/total>", 
                             "cvegag-<species/total>", 
                             "cvegbg-<species/total>", 
                             "clitter-<species/total>", 
                             "csoil", # actually "csoil-<species/total>" is correct
                             "age-<species/total>-<dbhclass/total>", 
                             "gpp-<species/total>", 
                             "npp-<species/total>", 
                             "ra-<species/total>", 
                             "rh-<total>", 
                             "nee-<total>", 
                             "mai-<species/total>", 
                             "fapar-<species/total>", 
                             "lai-<species/total>", 
                             "species-<species>", 
                             "evap", 
                             "intercept-<species/total>", 
                             "esoil", 
                             "trans-<species/total>", 
                             "soilmoist")


variable_opt_2A_naming <- c("mortstemno-<species/total>-<dbhclass/total>", 
                            "harvstemno-<species/total>-<dbhclass/total>", 
                            "dist-<dist_name>", 
                            "nlit-<species/total>", 
                            "nsoil-<total>", 
                            "nppleaf-<species>", # actually "npp-landleaf-<species>" is correct
                            "npproot-<species>", # actually "npp-landroot-<species>" is correct
                            "nppagwood-<species>", # actually "npp-abovegroundwood-<species>" is correct
                            "nppbgwood-<species>", # actually "npp-belowgroundwood-<species>" is correct
                            "rr-<species/total>", 
                            "cleaf-<species>", 
                            "cwood-<species>", 
                            "croot-<species>", 
                            "tsl")
variable_opt_2B_naming <- c("mortstemno-<species/total>-<dbhclass/total>", 
                            "harvstemno-<species/total>-<dbhclass/total>", 
                            "dist-<dist_name>", 
                            "nlit-<species/total>", 
                            "nsoil-<total>", 
                            "nppleaf-<species>", 
                            "npproot-<species>", 
                            "nppagwood-<species>", 
                            "nppbgwood-<species>", 
                            "rr-<species/total>", 
                            "cleaf-<species>", 
                            "cwood-<species>", 
                            "croot-<species>", 
                            "tsl")

variable_opt_3A_naming <- c("mortstemno-<species/total>-<dbhclass/total>", 
                            "harvstemno-<species/total>-<dbhclass/total>", 
                            "dist-<dist_name>", 
                            "nlit-<species/total>", 
                            "nsoil-<total>", 
                            "nppleaf-<species>", 
                            "npproot-<species>", 
                            "nppagwood-<species>", 
                            "nppbgwood-<species>", 
                            "rr-<species/total>", 
                            "cleaf-<species>", 
                            "cwood-<species>", 
                            "croot-<species>", 
                            "tsl")


## file naming
# coordinates: "hyytiala", "peitz", "solling-beech", "solling-spruce", "soro", "kroof", "le-bray", "collelongo", "bily-kriz" (taken from ISIMIP2B protocol)
site_lats <- c(61.848, 51.917, 51.77, 51.77, 55.486, 48.25, 44.717, 41.849, 49.3)
site_lons <- c(24.295, 14.35, 9.57, 9.57, 11.645, 11.4, -0.769, 13.588, 18.32)

# scenario naming according to simulation round
if(ISIMIP_round == "2A") {
  gcm_naming <- gcm_2A_naming
  bias_correction_naming <- bias_correction_2A_naming
  climate_scenario_naming <- climate_scenario_2A_naming
  soc_scenario_naming <- soc_scenario_2A_naming
  co2sens_scenario_naming <- co2sens_scenario_2A_naming
  region_naming <- region_2A_naming
  species_naming <- species_2A_naming
  variable_mand_naming <- variable_mand_2A_naming
  variable_opt_naming <- variable_opt_2A_naming
}

if(ISIMIP_round == "2B") {
  gcm_naming <- gcm_2B_naming
  bias_correction_naming <- bias_correction_2B_naming
  climate_scenario_naming <- climate_scenario_2B_naming
  soc_scenario_naming <- soc_scenario_2B_naming
  co2sens_scenario_naming <- co2sens_scenario_2B_naming
  region_naming <- region_2B_naming
  species_naming <- species_2B_naming
  variable_mand_naming <- variable_mand_2B_naming
  variable_opt_naming <- variable_opt_2B_naming
}

if(ISIMIP_round == "3A") {
  gcm_naming <- gcm_3A_naming
  bias_correction_naming <- bias_correction_3A_naming
  climate_scenario_naming <- climate_scenario_3A_naming
  soc_scenario_naming <- soc_scenario_3A_naming
  co2sens_scenario_naming <- co2sens_scenario_3A_naming
  region_naming <- region_3A_naming
  species_naming <- species_3A_naming
  variable_mand_naming <- variable_mand_3A_naming
  variable_opt_naming <- variable_opt_3A_naming
}


## 4c (fs) experiment and variable specifier naming related to ISIMIP naming + simulated species per site (and management)
# for ISIMIP2A
if(ISIMIP_round == "2A") {
  # gcm_naming: "gswp3" "princeton" "watch" "wfdei" "localclim"
  fs_gcm_naming <- list(c("GSWP3", "gswp3"), c("pgfv2", "princton", "princ", "princeton"), c("watch", "wa"), c("watch+wfdei", "wa+wf"), c("obs", "lieb"))
  # bias_correction_naming: "wfd" "nobc"
  fs_bias_correction_naming <- list(c(NA), c("nobc"))
  # climate_scenario_naming: "rcp2p6" "rcp4p5" "rcp6p0" "rcp8p5" "hist" "presclim" "noclim"
  fs_climate_scenario_naming <- list(c(NA), c(NA), c(NA), c(NA), c("hist"), c(NA), c(NA))
  # soc_scenario_naming: "nosoc" "pressoc" "varsoc" "nat"
  fs_soc_scenario_naming <- list(c("om", "oman"), c(NA), c("m", "man"), c(NA))
  # co2sens_scenario_naimng: "co2" "noco2" "co2const"
  fs_co2sens_scenario_naming <- list(c("co2hist", "co2hist52", "52", "co2"), c(NA), c(NA))
  # region_naming: "hyytiala" "peitz" "solling_beech" "solling_spruce" "soro" "kroof" "le_bray" "collelongo" "bily_kriz"
  fs_region_naming <- list(c("hy"), c("peitz", "p", "pz"), c("Soll304", "S304"), c("Soll305", "S305"), c("Sor"), c("kb"), c(NA), c("colle", "c"), c("bilyk"))
  # species_naming: "fasy", "quro", "qupe", "pisy", "piab", "pipi", "lade", "acpl", "eugl", "bepe", "bepu", "rops", "frex", "poni", "soau", "hawo"
  fs_species_naming <- list(c("be"), c("oa"), c(NA), c("pi"), c("sp"), c(NA), c(NA), c(NA), c(NA), c("bi"), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA))
  # 4C species ids fs_species_naming: c("be"), c("oa"), c(NA), c("pi"), c("sp"), c(NA), c(NA), c(NA), c(NA), c("bi"), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c("dg")
  fs_species_ids <- list(c(1), c(4), c(NA), c(3), c(2), c(NA), c(NA), c(NA), c(NA), c(5), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(10))
}

# for ISIMIP2B
if(ISIMIP_round == "2B") {
  # gcm_naming: "hadgem2-es" "ipsl-cm5a-lr" "miroc5" "gfdl-esm2m"
  fs_gcm_naming <- list(c(NA), c("IPSLCM5ALR"), c("MIROC5"), c("GFDLESM2M"))
  # bias_correction_naming: "nobc" "localbc" "ewembi" "ewembi-isimip3basd"
  fs_bias_correction_naming <- list(c("nobc"), c("localbc"), c(NA), c(NA))
  # climate_scenario_naming: "picontrol" "historical" "rcp26" "rcp60" "rcp85"
  fs_climate_scenario_naming <- list(c("pic"), c(NA), c("rcp2p6"), c("rcp6p0"), c(NA))
  # soc_scenario_naming: "nosoc" "1860soc" "histsoc" "2005soc" "rcp26soc" "rcp60soc" "2100rcp26soc", "2005socsite", "socbe", "sochwp", "soca", "socam"
  fs_soc_scenario_naming <- list(c("om"), c(NA), c(NA), c("mm"), c(NA), c(NA), c(NA), c("CSSm", "CSS-MFAm", "CSS-HWPm"), c("BEm"), c("HWPm"), c("MFAm", "MFA1m"), c("MFA2m"))
  # co2sens_scenario_naimng: "co2" "2005co2"
  fs_co2sens_scenario_naming <- list(c("co2", "286"), c("oco2"))
  # region_naming: "hyytiala" "peitz" "solling-beech" "solling-spruce" "soro" "kroof" "le-bray" "collelongo" "bily-kriz"
  fs_region_naming <- list(c("hy"), c("peitz"), c("soll304"), c("soll305"), c("sor"), c("kb"), c(NA), c("collelongo"), c("bilyk"))
  # species_naming: "fasy", "quro", "qupe", "pisy", "piab", "pipi", "lade", "acpl", "eugl", "bepe", "bepu", "rops", "frex", "poni", "soau", "c3gr", "hawo", "psme"
  fs_species_naming <- list(c("be"), c("oa"), c(NA), c("pi"), c("sp"), c(NA), c(NA), c(NA), c(NA), c("bi"), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c("dg"))
  # 4C species ids fs_species_naming: c("be"), c("oa"), c(NA), c("pi"), c("sp"), c(NA), c(NA), c(NA), c(NA), c("bi"), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c("dg")
  fs_species_ids <- list(c(1), c(4), c(NA), c(3), c(2), c(NA), c(NA), c(NA), c(NA), c(5), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(10))
}

if(ISIMIP_round == "3A") {
  # gcm_naming: "chelsa-w5e5"
  fs_gcm_naming <- list(c("chelsa-w5e5"))
  # bias_correction_naming: "nobc"
  fs_bias_correction_naming <- list(c("nobc"))
  # climate_scenario_naming: "obsclim"
  fs_climate_scenario_naming <- list(c("obsclim"))
  # soc_scenario_naming: "histsoc"
  fs_soc_scenario_naming <- list(c("histsoc"))
  # co2sens_scenario_naimng: "30arcsec" "90arcsec" "300arcsec" "1800arcsec"
  fs_co2sens_scenario_naming <- list(c("30arcsec"), c("90arcsec"), c("300arcsec"), c("1800arcsec"))
  # region_naming: "hyytiala", "peitz", "solling-beech", "solling-spruce", "soro", "kroof", "le-bray", "collelongo", "bily-kriz"
  fs_region_naming <- list(c(NA), c(NA), c("sb"), c("ss"), c(NA), c(NA), c(NA), c("collelongo"), c("bily-kriz"))
  # species_naming: "fasy", "quro", "qupe", "pisy", "piab", "pipi", "lade", "acpl", "eugl", "bepe", "bepu", "rops", "frex", "poni", "soau", "c3gr", "hawo", "psme"
  fs_species_naming <- list(c("be"), c("oa"), c(NA), c("pi"), c("sp"), c(NA), c(NA), c(NA), c(NA), c("bi"), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c("dg"))
  # 4C species ids fs_species_naming: c("be"), c("oa"), c(NA), c("pi"), c("sp"), c(NA), c(NA), c(NA), c(NA), c("bi"), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c("dg")
  fs_species_ids <- list(c(1), c(4), c(NA), c(3), c(2), c(NA), c(NA), c(NA), c(NA), c(5), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(10))
}


# generic site specific species simulated
site_species <- list(c("pi", "sp"), # hy
                     c("pi"), # peitz
                     c("be"), # soll304
                     c("sp"), # soll305
                     c("be"), # sor
                     c("be", "sp"), # kb
                     c(NA), # le-bray
                     c("be"), # collelongo
                     c("sp")) # bilyk
# exceptions of simulated species for ISIMIP2A managment scenarios
site_species_peitzMFAm <- c("pi", "oa", "bi") # peitz MFAm
site_species_sollbeMFAonem <- c("be", "dg") # soll304 MFA1m
site_species_sollbeMFAtwom <- c("be", "sp", "bi", "dg") # soll304 MFA2m
site_species_sollspMFAm <- c("sp", "be", "bi", "oa") # soll305 MFAm
site_species_sorMFAm <- c("be", "dg") # sor MFAm
site_species_collelongoMFAm <- c("be", "oa") # collelongo MFAm
site_species_bilykMFAm <- c("sp", "be") # bilyk MFAm

# alphac values per species for biomass compartment fractions
# 4C species ids fs_species_naming: c("be"), c("oa"), c(NA), c("pi"), c("sp"), c(NA), c(NA), c(NA), c(NA), c("bi"), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c(NA), c("dg")
alphac_list <- c(0.48, 0.56, NA, 0.46, 0.50, NA, NA, NA, NA, 0.50, NA, NA, NA, NA, NA, NA, NA, 0.54)

# patch sizes of the sites in m²
# region_naming: "hyytiala" "peitz" "solling-beech" "solling-spruce" "soro" "kroof" "le-bray" "collelongo" "bily-kriz"
patch_size_list <- c(10000, 1000, 10000, 10000, 10000, 5000, NA, 2000, 2500)

###################################################################################
### 02: 4C output extraction and netcdf file creation
###################################################################################

### compile all experiments and associated files as simulated by 4c and available in input directory
allfiles <- list.files(path = input, recursive = TRUE)

# summary of files related to the experiments available
sumfiles <- suppressWarnings(do.call(rbind,sapply(1:length(allfiles),function(i) {strsplit(allfiles[i],"\\/|\\_|\\.")}))) # suppress warning: number of columns of result is not a multiple of vector length (arg 1)

if(ISIMIP_round %in% c("2B", "3A"))) {
  sumfiles <- sumfiles[,1:6]
  sumfiles <- as.data.frame(sumfiles)
  colnames(sumfiles) <- c("bias_correction","region","gcm","soc_scenario","climate_scenario","co2sens_scenario")
}

# correct climate-scenario and co2sens-scenario for ISIMIP2A because the climate-scenario is not in the file name
if(ISIMIP_round == "2A") {
  sumfiles <- sumfiles[,1:5]
  sumfiles <- as.data.frame(sumfiles)
  colnames(sumfiles) <- c("bias_correction","region","gcm","soc_scenario","co2sens_scenario")
  
  sumfiles$climate_scenario <- "hist"
  sumfiles$co2sens_scenario <- "co2"
  sumfiles$climate_scenario <- "hist"
}

# these are the unique experiments found in the input directory:
experiments <- unique(sumfiles)

for (tmp_experiment in 1:length(experiments[,1])) {

  print(experiments[tmp_experiment,])

  ## define current experiment
  tmp_gcm <- as.character(experiments[tmp_experiment, "gcm"])
  tmp_bias_correction <- as.character(experiments[tmp_experiment, "bias_correction"])
  tmp_climate_scenario <- as.character(experiments[tmp_experiment, "climate_scenario"])
  tmp_soc_scenario <- as.character(experiments[tmp_experiment, "soc_scenario"])
  tmp_co2sens_scenario <- as.character(experiments[tmp_experiment, "co2sens_scenario"])
  tmp_region <- as.character(experiments[tmp_experiment, "region"])
  
  # skip loop for soll304 MFA2m because it throws an error
  if(tmp_soc_scenario == "MFA2m" & tmp_region == "soll304") {next}
  # skip loop for oco co2sens_scenario because it throws an error
  if(tmp_co2sens_scenario == "oco") {next}


  # define species simulated on current site-management combination
  # usual case
  curr_exp_specs <- site_species[mapply('%in%', tmp_region, fs_region_naming)][[1]]
  # exceptions for certain management scenarios (mostly MFA) at specific sites in ISIMIP2B
  if(ISIMIP_round == "2B") {
    if(tmp_region == fs_region_naming[[2]] && tmp_soc_scenario == fs_soc_scenario_naming[[11]][1]) {curr_exp_specs <- site_species_peitzMFAm} # peitz MFAm
    if(tmp_region == fs_region_naming[[3]] && tmp_soc_scenario == fs_soc_scenario_naming[[11]][2]) {curr_exp_specs <- site_species_sollbeMFAonem} # soll304 MFA1m
    if(tmp_region == fs_region_naming[[3]] && tmp_soc_scenario == fs_soc_scenario_naming[[12]]) {curr_exp_specs <- site_species_sollbeMFAtwom} # soll304 MFA2m
    if(tmp_region == fs_region_naming[[4]] && tmp_soc_scenario == fs_soc_scenario_naming[[11]][1]) {curr_exp_specs <- site_species_sollspMFAm} # soll305 MFAm
    if(tmp_region == fs_region_naming[[5]] && tmp_soc_scenario == fs_soc_scenario_naming[[11]][1]) {curr_exp_specs <- site_species_sorMFAm} # sor MFAm
    if(tmp_region == fs_region_naming[[8]] && tmp_soc_scenario == fs_soc_scenario_naming[[11]][1]) {curr_exp_specs <- site_species_collelongoMFAm} # colellongo MFAm
    if(tmp_region == fs_region_naming[[9]] && tmp_soc_scenario == fs_soc_scenario_naming[[11]][1]) {curr_exp_specs <- site_species_bilykMFAm} # bilyk MFAm
  }
  
  # select temporary experiment path of the 4C outputs
  tmp_InputPath <- paste(input, "/", tmp_bias_correction, sep="")
  tmp_InputFiles <- allfiles[which(sumfiles$gcm==tmp_gcm 
                                   & sumfiles$bias_correction==tmp_bias_correction 
                                   & sumfiles$climate_scenario==tmp_climate_scenario 
                                   & sumfiles$soc_scenario==tmp_soc_scenario 
                                   & sumfiles$co2sens_scenario==tmp_co2sens_scenario 
                                   & sumfiles$region==tmp_region)]
  fs_tmp_experiment_output <- paste(input, "/", tmp_InputFiles, sep="")  
  
  # the outcommented two lines may be obsolete
  #tmp_InputName <- strsplit(tmp_InputFiles[grepl("*\\.ctr*", tmp_InputFiles)], "\\/|\\.")[[1]][2]
  #fs_tmp_experiment_output <- dir(tmp_InputPath, full.names=T, pattern = tmp_InputName)
  
  patch_size <- patch_size_list[mapply('%in%', tmp_region, fs_region_naming)]
  
  ###################################################################################
  ### 02a: read 4C output files for given experiment
  ###################################################################################

  ## mandatory
  # annual 4C output
  veg_file <- read.table(fs_tmp_experiment_output[grepl("_veg\\.out" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  veg_spec_files <- c(NULL)
  veg_spec_filenames <- fs_tmp_experiment_output[grepl("_veg_" ,fs_tmp_experiment_output)]  #select all veg_spec files of current experiment
  for (veg_spec_filename in veg_spec_filenames) {
    spec <- strsplit(veg_spec_filename, "\\/|\\_|\\.")[[1]][length(strsplit(veg_spec_filename, "\\/|\\_|\\.")[[1]])-1]  # extract 4c species
    if(spec %in% curr_exp_specs) { # only if current species was simulated in the current experiment, include it (info which species was simulated in the current experiment is in curr_exp_specs)
      veg_spec_files <- c(veg_spec_files, paste("veg_", spec, "_file", sep=""))
      eval(parse(text=paste("veg_", spec, "_file", " <- read.table(\"", veg_spec_filename,"\"" , ", header=F, comment.char=\"#\")", sep="")))
    }
  }
  
  classmvol_file <- read.table(fs_tmp_experiment_output[grepl("_classmvol" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  classd_file <- read.table(fs_tmp_experiment_output[grepl("_classd.out" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  litter_file <- read.table(fs_tmp_experiment_output[grepl("_litter" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  soil_file <- read.table(fs_tmp_experiment_output[grepl("_soil.out" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  classage_file <- read.table(fs_tmp_experiment_output[grepl("_classage" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  # cohort output files
  coh_atr_file <- read.table(fs_tmp_experiment_output[grepl("_atr\\." ,fs_tmp_experiment_output)], col.names = c("year", paste0("c",seq_len(5000))), header=F, comment.char="#", fill=T)
  coh_length <- length(coh_atr_file[length(coh_atr_file[,1]),][!is.na(coh_atr_file[length(coh_atr_file[,1]),])]) # determine how long the last row is
  coh_atr_file <- coh_atr_file[, 1:coh_length]; coh_atr_file[coh_atr_file==-99.9] <- NA # subset file to only filled cohort cols and set all -99.9 values do NA
  coh_atr_file[,-1] <- coh_atr_file[,-1] * 10000/patch_size # apply expansion factor to derive trees per hectare from trees on the patch
  coh_hei_file <- read.table(fs_tmp_experiment_output[grepl("_hei\\." ,fs_tmp_experiment_output)], col.names = c("year", paste0("c",seq_len(5000))), fill=T, header=F, comment.char="#")
  coh_hei_file <- coh_hei_file[, 1:coh_length]; coh_hei_file[coh_hei_file==-99.9] <- NA # subset file to only filled cohort cols and set all -99.9 values do NA
  coh_diam_file <- read.table(fs_tmp_experiment_output[grepl("_diam\\." ,fs_tmp_experiment_output)], col.names = c("year", paste0("c",seq_len(5000))), fill=T, header=F, comment.char="#")
  coh_diam_file <- coh_diam_file[, 1:coh_length]; coh_diam_file[coh_diam_file==-99.9] <- NA # subset file to only filled cohort cols and set all -99.9 values do NA
  
  copmfract_spec_file <- read.table(fs_tmp_experiment_output[grepl("_Copmfract\\." ,fs_tmp_experiment_output)], col.names = c("year", rep(c("spec", "c_opm_fol", "c_opm_tb", "c_opm_frt", "c_opm_crt", "c_opm_stm"), 60)), header=F, comment.char="#", fill=T)
  copmfract_spec_file <- copmfract_spec_file[colSums(!is.na(copmfract_spec_file)) > 0]
  
  copm_soil_file <- read.table(fs_tmp_experiment_output[grepl("_Copm\\.out" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  chum_soil_file <- read.table(fs_tmp_experiment_output[grepl("_Chum\\.out" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  coh_spn_file <- read.table(fs_tmp_experiment_output[grepl("_spn\\." ,fs_tmp_experiment_output)], col.names = c("year", paste0("c",seq_len(5000))), header=F, comment.char="#", fill=T)
  coh_spn_file <- coh_spn_file[, 1:coh_length]; coh_spn_file[coh_spn_file==-99.9] <- NA # subset file to only filled cohort cols and set all -99.9 values do NA
  
  water_res_file <- read.table(fs_tmp_experiment_output[grepl("_water\\.res" ,fs_tmp_experiment_output)], fill=T, header=T, comment.char="#")
  
  coh_respaut_file <- read.table(fs_tmp_experiment_output[grepl("_respaut\\." ,fs_tmp_experiment_output)], col.names = c("day", "year", paste0("c",seq_len(5000))), header=F, comment.char="#", fill=T)
  coh_respaut_file <- coh_respaut_file[, 1:(coh_length+1)]; coh_respaut_file[coh_respaut_file==-99.9] <- NA # subset file to only filled cohort cols and set all -99.9 values do NA
  
  coh_totfpar_file <- read.table(fs_tmp_experiment_output[grepl("_totfpar\\." ,fs_tmp_experiment_output)], col.names = c("day", "year", paste0("c",seq_len(5000))), header=F, comment.char="#", fill=T)
  coh_totfpar_file <- coh_totfpar_file[, 1:(coh_length+1)]; coh_totfpar_file[coh_totfpar_file==-99.9] <- NA # subset file to only filled cohort cols and set all -99.9 values do NA
  
  coh_aevi_file <- read.table(fs_tmp_experiment_output[grepl("_aevi\\." ,fs_tmp_experiment_output)], col.names = c("day", "year", paste0("c",seq_len(5000))), header=F, comment.char="#", fill=T)
  coh_aevi_file <- coh_aevi_file[, 1:(coh_length+1)]; coh_aevi_file[coh_totfpar_file==-99.9] <- NA # subset file to only filled cohort cols and set all -99.9 values do NA
  
  coh_respfrt_file <- read.table(fs_tmp_experiment_output[grepl("_respfrt\\." ,fs_tmp_experiment_output)], col.names = c("day", "year", paste0("c",seq_len(5000))), header=F, comment.char="#", fill=T)
  coh_respfrt_file <- coh_respfrt_file[, 1:(coh_length+1)]; coh_respfrt_file[coh_respfrt_file==-99.9] <- NA # subset file to only filled cohort cols and set all -99.9 values do NA
  
  coh_spn_day_file <- coh_respaut_file
  coh_atr_day_file <- coh_respaut_file
  for(year in unique(coh_respaut_file$year)[-length(unique(coh_respaut_file$year))]) {
    coh_spn_day_file[which(coh_spn_day_file$year==year), c(-1,-2)] <- coh_spn_file[which(coh_spn_file$year==year-1), c(-1)]
    coh_atr_day_file[which(coh_atr_day_file$year==year), c(-1,-2)] <- coh_atr_file[which(coh_atr_file$year==year-1), c(-1)]
  }
  
  # daily 4C output
  Cday_file <- read.table(fs_tmp_experiment_output[grepl("_Cday" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  day_file <- read.table(fs_tmp_experiment_output[grepl("_day\\.out" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  water_file <- read.table(fs_tmp_experiment_output[grepl("_water\\.out" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  # one-time 4C output: soil_ini needed for calculating variables related to soil horizons
  soil_ini_file <-  read.csv(fs_tmp_experiment_output[grepl("_soil.ini" ,fs_tmp_experiment_output)], sep="", header=F, comment.char="!", skip=2)
  is_of <- F
  rownumber <- 1
  ini_soil_horizons_depths <- NULL
  while(!is_of) { # extract initial soil horizon depths from soil_ini_file
    if(rownumber > length(soil_ini_file[,2])) {break}
    if(soil_ini_file[rownumber,2] == "of") {is_of <- T}
    if(!is_of) {ini_soil_horizons_depths <- c(ini_soil_horizons_depths, as.numeric(as.character(soil_ini_file[rownumber,2])))}
    rownumber <- rownumber + 1
  }
  ini_soil_horizons_depths <- ini_soil_horizons_depths / 100  # convert cm to m
  
  ## optional
  # annual 4C output
  classt_file <- read.table(fs_tmp_experiment_output[grepl("_classt" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  classdm_file <- read.table(fs_tmp_experiment_output[grepl("_classdm\\.out" ,fs_tmp_experiment_output)], header=F, comment.char="#")
  
  # daily 4C output
  temp_file <-  read.table(fs_tmp_experiment_output[grepl("_temp" ,fs_tmp_experiment_output)][!grepl("clim_temp" ,fs_tmp_experiment_output[grepl("_temp" ,fs_tmp_experiment_output)])], header=F, comment.char="#")
  
  ###################################################################################
  ### 02b: extract data for all relevant variables and write them into a joint data.frame
  ###      called df_ann, df_mon or df_day according to their temporal resolution
  ###################################################################################
  
  ## create a dataframe df that is used to store all variables for one single experiment with first cols lon, lat and time (in years)
  # specify lon and lat columns for site
  lon_site <- site_lons[mapply('%in%', tmp_region, fs_region_naming)]
  lat_site <- site_lats[mapply('%in%', tmp_region, fs_region_naming)]
  # specify time column for three different temporal resolutions: daily, monthly, annual
  date_sequence_days <- seq(as.Date(ISOdate(day_file[1,2], 01, 01)), as.Date(ISOdate(veg_file[length(veg_file[,1]),1], 12, 31)), "days")
  date_sequence_months <- do.call(paste, expand.grid(as.character(seq(veg_file[2,1], veg_file[length(veg_file[,1]),1])), c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")))
  date_sequence_months <- gsub(" ", "-", date_sequence_months)[order(gsub(" ", "-", date_sequence_months))]
  date_sequence_years <- as.character(seq(veg_file[2,1], veg_file[length(veg_file[,1]),1]))
  # construct three different data.frames df_* for daily, monthly and annually resolved variables
  df_day <- data.frame(lon = rep(lon_site, length(date_sequence_days)), lat = rep(lat_site, length(date_sequence_days)), time = date_sequence_days)
  df_mon <- data.frame(lon = rep(lon_site, length(date_sequence_months)), lat = rep(lat_site, length(date_sequence_months)), time = date_sequence_months)
  df_ann <- data.frame(lon = rep(lon_site, length(date_sequence_years)), lat = rep(lat_site, length(date_sequence_years)), time = date_sequence_years)
  
  # tmp data.frames for classed variables (cbh-classes and soil layer classes)
  df_ann_dbhclassed <- data.frame(time=df_ann[, "time"])
  df_ann_depthclassed <- data.frame(time=df_ann[, "time"])
  df_day_depthclassed <- data.frame(time=df_day[, "time"])
  
  ## mandatory variables (2A + 2B)
  # dbh-<species/total> from *.veg_spec/*.veg : mean diam
  df_ann[gsub("-<species/total>", "_total", variable_mand_naming[1])] <- veg_file[-length(veg_file[,1]),24] # shift of one year
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      var_spec <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),25]  # shift of one year
      df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[1])] <- var_spec
    }
  }
  
  # dbbdomhei from cohort output
  coh_atr_file[is.na(coh_atr_file)] <- 0
  
  dbhdomhei_timeseries <- NULL
  for(curr_year in 2:length(coh_atr_file[,1])) {  # iterate over all years;  shift of one year
    hei_order <- order(coh_hei_file[curr_year-1,-1], decreasing = T)
    trees <- NULL
    for(curr_coh in 1:length(hei_order)) {
      trees <- append(trees, rep(coh_diam_file[curr_year-1,-1][,hei_order][,curr_coh], coh_atr_file[curr_year-1,-1][,hei_order][,curr_coh]))
    }
    if(length(trees) >= 100) {new_dbhdomhei <- mean(trees[1:100])} else {new_dbhdomhei <- mean(trees)}
    
    dbhdomhei_timeseries <- append(dbhdomhei_timeseries, new_dbhdomhei)
  }
  df_ann[variable_mand_naming[2]] <- dbhdomhei_timeseries
  
  # hei-<species/total> / height-<species/total> from *.veg_spec/*.veg : mean height
  df_ann[gsub("-<species/total>", "_total", variable_mand_naming[3])] <- veg_file[-length(veg_file[,1]),25] /100  # convert from cm to m; shift of one year
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      var_spec <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),26]/100  # convert from cm to m
      df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[3])] <- var_spec
    }
  }
  
  # domhei / domheight from *.veg : domhei *1/100
  df_ann[variable_mand_naming[4]] <- veg_file[-length(veg_file[,1]),9] /100  # convert from cm to m
  
  # density-<species/total> from *.veg_spec/*.veg : Tree
  df_ann[gsub("-<species/total>", "_total", variable_mand_naming[5])] <- veg_file[-length(veg_file[,1]),4]
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      var_spec <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),4]
      df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[5])] <- var_spec
    }
  }
  
  # ba-<species/total> from *.veg_spec/*.veg : basal_area
  df_ann[gsub("-<species/total>", "_total", variable_mand_naming[6])] <- veg_file[-length(veg_file[,1]),26]
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      var_spec <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),27]
      df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[6])] <- var_spec
    }
  }
  
  # mort-<species/total> from *.veg_spec/*.veg : dead_stem_m3
  df_ann[gsub("-<species/total>", "_total", variable_mand_naming[7])] <- veg_file[-length(veg_file[,1]),27]
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      var_spec <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),28]
      df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[7])] <- var_spec
    }
  }
  
  # harv-<species/total>-<dbhclass/total> from *.classmvol
  #df_ann[gsub("-<species/total>-<dbhclass/total>", "_total_total", variable_mand_naming[8])] <- rowSums(classmvol_file[,-1])
  # calculate total per dbhclass
  for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
    var_spec <- rowSums(classmvol_file[-1,-1][, seq(tmp_dbh_class, (29*14), 29)])
    df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
  }
  list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
  list <- as.data.frame(matrix(list))
  df_ann[gsub("-<species/total>-<dbhclass/total>", "_total", variable_mand_naming[8])] <- as.vector(list)
  
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) { # iterate over all species
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
        var_spec <- classmvol_file[-1,-1][,(((curr_spec_id_fs-1) * 29)+1) + tmp_dbh_class]
        df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
      }
      list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
      list <- as.data.frame(matrix(list))
      df_ann[gsub("-<species/total>-<dbhclass/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[8])] <- as.vector(list)
      
      # calculate total per species
      #df_ann[gsub("-<species/total>-<dbhclass/total>", paste("_", tmp_species_ISIMIP, "_total", sep=""), variable_mand_naming[8])] <- rowSums(classmvol_file[,-1][,(((curr_spec_id_fs[[1]]-1) * 29)+1):(((curr_spec_id_fs[[1]]-1) * 29)+29)])
    }
  }
  
  # stemno-<species/total>-<dbhclass/total> from *.classd
  # do not calculate total over all dbhclasses because it is already in density variable
  # calculate total per dbhclass
  for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
    var_spec <- rowSums(classd_file[-1,-1][, seq(tmp_dbh_class, (29*14), 29)])
    df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
  }
  list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
  list <- as.data.frame(matrix(list))
  df_ann[gsub("-<species/total>-<dbhclass/total>", "_total", variable_mand_naming[9])] <- as.vector(list)
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) { # iterate over all species
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
        var_spec <- classd_file[-1,-1][,(((curr_spec_id_fs-1) * 29)) + tmp_dbh_class]
        df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
      }
      list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
      list <- as.data.frame(matrix(list))
      df_ann[gsub("-<species/total>-<dbhclass/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[9])] <- as.vector(list)
      # do not calculate total per species because it is already in density variable
    }
  }
  
  # vol-<species/total> from *.veg_spec/*.veg : Stemvol
  df_ann[gsub("-<species/total>", "_total", variable_mand_naming[10])] <- veg_file[-length(veg_file[,1]),15]
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      var_spec <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),15]
      df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[10])] <- var_spec
    }
  }
  
  # cveg-<species/total> from *.veg_spec/*.veg : biomass *1/20000
  df_ann[gsub("-<species/total>", "_total", variable_mand_naming[11])] <- veg_file[-length(veg_file[,1]),6] /20000 # conversion from kg_DW/ha biomass into kgC/m²
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      var_spec <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),6] /20000  # conversion from kg_DW/ha biomass into kgC/m²
      df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[11])] <- var_spec
    }
  }
  
  # clitter-<species/total> from *litter : tot_litter *1/10000
  df_ann[gsub("-<species/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 12, 14)])] <- litter_file[-length(litter_file[,1]),15] /10000 # conversion from kgC/ha to kgC/m²
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated TOCHECK
    for(veg_spec_file in veg_spec_files) {
      clitter_timeseries <- NULL
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      for(curr_year in 2:length(copmfract_spec_file[,1])) {
        curr_data <- copmfract_spec_file[curr_year-1,]
        # determine which cols correspond to current species
        if(!any(curr_data[1,]==curr_spec_id_fs, na.rm=T)) {
          clitter_timeseries <- append(clitter_timeseries, NA)
          next
        }
        pot_curr_spec_cols <- which(curr_data[1,]==curr_spec_id_fs)
        spec_cols <- which(grepl("spec", colnames(curr_data)))
        pot_cols <- c(pot_curr_spec_cols, spec_cols)
        spec_col <- pot_cols[duplicated(pot_cols)]
        
        new_clitter <- sum(curr_data[, (spec_col+1):(spec_col+5)])
        new_clitter <- new_clitter /1000  # conversion from gC/m² biomass into kgC/m²
        clitter_timeseries <- append(clitter_timeseries, new_clitter)
      }
      
      df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[ifelse(ISIMIP_round=="2A", 12, 14)])] <- clitter_timeseries
    }
  }
  
  # csoil-<species/total> from *soil : C_tot *1/1000 , split up for each soil layer
  for (layer in 1:length(copm_soil_file[1,-1])) {
    df_ann_depthclassed[, paste("soil_layer_", layer, sep="")] <- (copm_soil_file[,-1][,layer]+chum_soil_file[,-1][,layer])[-length(copm_soil_file[,1])] /1000 # conversion from gC/m² to kgC/m²
  }
  list <- as.list(as.data.frame(t(df_ann_depthclassed[,-1])))
  list <- as.data.frame(matrix(list))
  df_ann[paste(variable_mand_naming[ifelse(ISIMIP_round=="2A", 13, 15)], sep="")] <- as.vector(list)
  
  # age-<species/total>-<dbhclass/total> from *classage
  classage_file[classage_file==0] <- NA  # replace all 0s with NAs
  
  # df_ann[gsub("-<species/total>-<dbhclass/total>", "_total_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 14, 16)])] <- rowMeans(classage_file[,-1], na.rm=T)
  # calculate total per dbhclass
  for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
    var_spec <- rowMeans(classage_file[-length(classage_file[,1]),-1][, seq(tmp_dbh_class, (29*14), 29)], na.rm=T)
    df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
  }
  list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
  list <- as.data.frame(matrix(list))
  df_ann[gsub("-<species/total>-<dbhclass/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 14, 16)])] <- as.vector(list)
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) { # iterate over all species
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
        var_spec <- classage_file[-length(classage_file[,1]),][,(((curr_spec_id_fs-1) * 29)+1) + tmp_dbh_class]
        df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
      }
      list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
      list <- as.data.frame(matrix(list))
      df_ann[gsub("-<species/total>-<dbhclass/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[ifelse(ISIMIP_round=="2A", 14, 16)])] <- as.vector(list)
      
      # calculate total per species
      # df_ann[gsub("-<species/total>-<dbhclass/total>", paste("_", tmp_species_ISIMIP, "_total", sep=""), variable_mand_naming[ifelse(ISIMIP_round=="2A", 14, 16)])] <- rowMeans(classage_file[,-1][,(((curr_spec_id_fs[[1]]-1) * 29)+1):(((curr_spec_id_fs[[1]]-1) * 29)+29)], na.rm=T)
    }
  }
  
  # gpp-<species/total> from : *_Cday : GPP *1/(86400*1000)
  df_day[gsub("-<species/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 15, 17)])] <- Cday_file[,9] /86400000 # conversion from gC/m² into kgC/m²s
  
  # npp-<species/total> from : *_Cday : NPP *1/(86400*1000)
  df_day[gsub("-<species/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 16, 18)])] <- Cday_file[,8] /86400000 # conversion from gC/m² into kgC/m²s
  
  # ra-<species/total> from : *_Cday : resp_aut *1/(86400*1000)
  df_day[gsub("-<species/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 17, 19)])] <- Cday_file[,12] /86400000 # conversion from gC/m² into kgC/m²s
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated TOCHECK
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      
      var_spec <- coh_respaut_file[, c(-1,-2)]*coh_atr_day_file[, c(-1,-2)] # multiply ra per cohort with number of alive trees in cohorts
      coh_spn_day_file_mask <- coh_spn_day_file[, c(-1,-2)]
      coh_spn_day_file_mask[coh_spn_day_file_mask != curr_spec_id_fs] <- 0
      coh_spn_day_file_mask[coh_spn_day_file_mask == curr_spec_id_fs] <- 1
      var_spec <- var_spec * coh_spn_day_file_mask # remove all values from other species by multiplying with 0
      
      ra_timeseries <- rowSums(var_spec, na.rm=T)
      ra_timeseries <- ra_timeseries /864000000000 # conversion from gC/ha d into kgC/m²s
      
      df_day[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[ifelse(ISIMIP_round=="2A", 17, 19)])] <- ra_timeseries
    }
  }
  
  # rh-<total> from : *_Cday : resp_het *1/(86400*1000)
  df_day[gsub("-<total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 18, 20)])] <- Cday_file[,14] /86400000 # conversion from gC/m² into kgC/m²s
  
  # nee-<total> from : *_Cday : NEE *1/(86400*1000)
  df_day[gsub("-<total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 19, 21)])] <- Cday_file[,10] /86400000 # conversion from gC/m² into kgC/m²s
  
  # mai-<species/total> from *.veg_spec/*.veg : stem_inc
  df_ann[gsub("-<species/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 20, 22)])] <- veg_file[-1,28]
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      var_spec <- eval(parse(text=veg_spec_file))[-1,29]
      df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[ifelse(ISIMIP_round=="2A", 20, 22)])] <- var_spec
    }
  }
  
  # fapar-<species/total> from : *_Cday : FaPar (%)
  df_day[gsub("-<species/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 21, 23)])] <- Cday_file[,16]
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated TOCHECK
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      
      var_spec <- coh_totfpar_file[, c(-1,-2)]*coh_atr_day_file[, c(-1,-2)] # multiply ra per cohort with number of alive trees in cohorts
      coh_spn_day_file_mask <- coh_spn_day_file[, c(-1,-2)]
      coh_spn_day_file_mask[coh_spn_day_file_mask != curr_spec_id_fs] <- 0
      coh_spn_day_file_mask[coh_spn_day_file_mask == curr_spec_id_fs] <- 1
      var_spec <- var_spec * coh_spn_day_file_mask # remove all values from other species by multiplying with 0
      
      fapar_timeseries <- rowSums(var_spec, na.rm=T)
      fapar_timeseries <- fapar_timeseries *100  # convert from proportion to %
      
      df_day[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[ifelse(ISIMIP_round=="2A", 21, 23)])] <- fapar_timeseries
    }
  }
  
  # lai-<species/total> from : *_day : LAI , monthly averaging
  lai_day <- day_file[,21]
  lai_day_help <- data.frame(month=format(date_sequence_days, "%m"), year=format(date_sequence_days, "%Y"), lai=lai_day)
  lai_mon <- aggregate(lai ~ month + year, data=lai_day_help, max) # TOCHECK is lai the monthly mean or max value of daily values
  
  df_mon[gsub("-<species/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 22, 24)])] <- lai_mon[,3]
  
  # species-<species> from *.veg_spec/*.veg : Tree , calculate percentages
  ba_total <- veg_file[-length(veg_file[,1]),26]
  
  for(veg_spec_file in veg_spec_files) {
    curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
    tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
    var_spec <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),27] / ba_total * 100 # calculate the species portion as percentage of total basal area
    
    var_spec[which(!(var_spec >= 0 & var_spec <= 100))] <- NA   # check for inf, >100, <0 values and set NA
    var_spec[which(is.na(var_spec))] <- NA  # check for NaN values and set NA
    
    df_ann[gsub("-<species>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[ifelse(ISIMIP_round=="2A", 23, 25)])] <- var_spec
  }
  
  # evap from *_day : AET[mm/day] *1/86400
  df_day[variable_mand_naming[ifelse(ISIMIP_round=="2A", 24, 26)]] <- day_file[,9] /86400 # conversion from mm/d to mm/s
  
  # intercept-<species/total> from *_day : intercep [mm/day) *1/86400
  df_day[gsub("-<species/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 25, 27)])] <- day_file[,6] /86400 # conversion from mm/d to mm/s
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated TOCHECK
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      
      var_spec <- coh_aevi_file[, c(-1,-2)]*coh_atr_day_file[, c(-1,-2)] # multiply ra per cohort with number of alive trees in cohorts
      coh_spn_day_file_mask <- coh_spn_day_file[, c(-1,-2)]
      coh_spn_day_file_mask[coh_spn_day_file_mask != curr_spec_id_fs] <- 0
      coh_spn_day_file_mask[coh_spn_day_file_mask == curr_spec_id_fs] <- 1
      var_spec <- var_spec * coh_spn_day_file_mask # remove all values from other species by multiplying with 0
      
      intercept_timeseries <- rowSums(var_spec, na.rm=T)
      intercept_timeseries <- intercept_timeseries /86400  # convert from mm/day to mm/s
      
      df_day[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[ifelse(ISIMIP_round=="2A", 25, 27)])] <- intercept_timeseries
    }
  }
  
  # esoil
  df_day[variable_mand_naming[ifelse(ISIMIP_round=="2A", 26, 28)]] <- water_res_file[,20] /86400 # conversion from mm/d to mm/s
  
  # trans-<species/total> from *_day : transtree( mm/day) *1/86400 (TODO: check if trans is really the variable asked for)
  df_day[gsub("-<species/total>", "_total", variable_mand_naming[ifelse(ISIMIP_round=="2A", 27, 29)])] <- day_file[,11] /86400 # conversion from mm/d to mm/s
  
  # soilmoist from *_water : wats_1, wats_2,... , for all depth layers : layers from soil.ini
  for (layer in 1:length(water_file[1,c(-1,-2)])) {
    df_day_depthclassed[, paste("soil_layer_", layer, sep="")] <- water_file[,layer+2]
  }
  list <- as.list(as.data.frame(t(df_day_depthclassed[,-1])))
  list <- as.data.frame(matrix(list))
  df_day[paste(variable_mand_naming[ifelse(ISIMIP_round=="2A", 28, 30)], sep="")] <- as.vector(list)
  
  if(ISIMIP_round %in% c("2B", "3A")) {
    # cvegag-<species/total> (2B only) from veg-file different biomass compartments
    cvegag_total <- NULL
    
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      alphac <- alphac_list[mapply('%in%', curr_spec_fs, fs_species_naming)]
      biomass <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),6]
      fol_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),10]
      sap_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),11]
      hrt_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),13]
      frt_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),12]
      
      tbc_bio <- biomass - (fol_bio + sap_bio + hrt_bio + frt_bio)
      tb_bio <- (1-alphac) * tbc_bio
      vegag <- sap_bio + hrt_bio + fol_bio + tb_bio
      
      cvegag <- vegag /20000  # conversion from kg_DW/ha biomass into kgC/m²
      if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
        df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[12])] <- cvegag
      }
      if(is.null(cvegag_total)) {cvegag_total <- cvegag} else {cvegag_total <- cvegag_total + cvegag} # add up species values to total
    }
    df_ann[gsub("-<species/total>", "_total", variable_mand_naming[12])] <- cvegag_total
    
    # cvegbg-<species/total> (2B only)
    cvegbg_total <- NULL
    
    for(veg_spec_file in veg_spec_files) {
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      alphac <- alphac_list[mapply('%in%', curr_spec_fs, fs_species_naming)]
      biomass <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),6]
      fol_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),10]
      sap_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),11]
      hrt_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),13]
      frt_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),12]
      
      tbc_bio <- biomass - (fol_bio + sap_bio + hrt_bio + frt_bio)
      crt_bio <- alphac * tbc_bio
      vegbg <- frt_bio + crt_bio
      
      cvegbg <- vegbg /20000  # conversion from kg_DW/ha biomass into kgC/m²
      if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
        df_ann[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_mand_naming[13])] <- cvegbg
      }
      if(is.null(cvegbg_total)) {cvegbg_total <- cvegbg} else {cvegbg_total <- cvegbg_total + cvegbg} # add up species values to total
    }
    df_ann[gsub("-<species/total>", "_total", variable_mand_naming[13])] <- cvegbg_total
  }
  
  ## optional variables (2A + 2B)
  # mortstemno-<species/total>-<dbhclass/total> from *_classt
  # df_ann[gsub("-<species/total>-<dbhclass/total>", "_total_total", variable_opt_naming[1])] <- rowSums(classt_file[,-1])
  # calculate total per dbhclass
  for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
    var_spec <- rowSums(classt_file[-length(classt_file[,1]),-1][, seq(tmp_dbh_class, (29*14), 29)])
    df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
  }
  list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
  list <- as.data.frame(matrix(list))
  df_ann[gsub("-<species/total>-<dbhclass/total>", "_total", variable_opt_naming[1])] <- as.vector(list)
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) { # iterate over all species
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
        var_spec <- classt_file[-length(classt_file[,1]),-1][,(((curr_spec_id_fs-1) * 29)) + tmp_dbh_class]
        df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
      }
      list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
      list <- as.data.frame(matrix(list))
      df_ann[gsub("-<species/total>-<dbhclass/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_opt_naming[1])] <- as.vector(list)
      
      # calculate total per species
      # df_ann[gsub("-<species/total>-<dbhclass/total>", paste("_", tmp_species_ISIMIP, "_total", sep=""), variable_opt_naming[1])] <- rowSums(classt_file[,-1][,(((curr_spec_id_fs[[1]]-1) * 29)+1):(((curr_spec_id_fs[[1]]-1) * 29)+29)])
    }
  }
  
  # harvstemno-<species/total>-<dbhclass/total> from *_classdm
  # df_ann[gsub("-<species/total>-<dbhclass/total>", "_total_total", variable_opt_naming[2])] <- rowSums(classdm_file[,-1])
  # calculate total per dbhclass
  for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
    var_spec <- rowSums(classdm_file[-length(classdm_file[,1]),-1][, seq(tmp_dbh_class, (29*14), 29)])
    df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
  }
  list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
  list <- as.data.frame(matrix(list))
  df_ann[gsub("-<species/total>-<dbhclass/total>", "_total", variable_opt_naming[2])] <- as.vector(list)
  
  if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
    for(veg_spec_file in veg_spec_files) { # iterate over all species
      curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
      curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
      tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
      for(tmp_dbh_class in 1:29) { # iterate over all diameter classes
        var_spec <- classdm_file[-length(classdm_file[,1]),-1][,(((curr_spec_id_fs-1) * 29)) + tmp_dbh_class]
        df_ann_dbhclassed[, paste(dbhclass_naming[tmp_dbh_class], sep="")] <- var_spec
      }
      list <- as.list(as.data.frame(t(df_ann_dbhclassed[,-1])))
      list <- as.data.frame(matrix(list))
      df_ann[gsub("-<species/total>-<dbhclass/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_opt_naming[2])] <- as.vector(list)
      
      # calculate total per species
      # df_ann[gsub("-<species/total>-<dbhclass/total>", paste("_", tmp_species_ISIMIP, "_total", sep=""), variable_opt_naming[2])] <- rowSums(classdm_file[,-1][,(((curr_spec_id_fs[[1]]-1) * 29)+1):(((curr_spec_id_fs[[1]]-1) * 29)+29)])
    }
  }
  
  # dist-<dist_name> NOT SIMULATED
  
  # nlit-<species/total> from *_litter : tot_litter ( gN/ha) *1/10
  df_ann[gsub("-<species/total>", "_total", variable_opt_naming[4])] <- litter_file[-length(litter_file[,1]),21] /10  # conversion from kgN/ha to gN/m²
  
  # nsoil-<total> from *_soil : N_tot (g N/m²)
  df_ann[gsub("-<total>", "_total", variable_opt_naming[5])] <- soil_file[-length(soil_file[,1]),15] /10  # conversion from kgN/ha to gN/m²
  
  # nppleaf-<species> NOT AVAILABLE IN 4C
  
  # npproot-<species> NOT AVAILABLE IN 4C
  
  # nppagwood-<species> NOT AVAILABLE IN 4C
  
  # nppbgwood-<species> NOT AVAILABLE IN 4C
  
  # rr-<species/total> from cohort respfrt TOCHECK
  for(veg_spec_file in veg_spec_files) {
    curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
    curr_spec_id_fs <- fs_species_ids[mapply('%in%', curr_spec_fs, fs_species_naming)][[1]]
    tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
    
    var_spec <- coh_respfrt_file[, c(-1,-2)]*coh_atr_day_file[, c(-1,-2)] # multiply ra per cohort with number of alive trees in cohorts
    coh_spn_day_file_mask <- coh_spn_day_file[, c(-1,-2)]
    coh_spn_day_file_mask[coh_spn_day_file_mask != curr_spec_id_fs] <- 0
    coh_spn_day_file_mask[coh_spn_day_file_mask == curr_spec_id_fs] <- 1
    var_spec <- var_spec * coh_spn_day_file_mask # remove all values from other species by multiplying with 0
    
    rr_timeseries <- rowSums(var_spec, na.rm=T)
    rr_timeseries <- rr_timeseries /864000000000 # conversion from gC/ha d into kgC/m²s
    
    if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
      df_day[gsub("-<species/total>", paste("_", tmp_species_ISIMIP, sep=""), variable_opt_naming[10])] <- rr_timeseries
    } else {
      df_day[gsub("-<species/total>", "_total", variable_opt_naming[10])] <- rr_timeseries
    }
  }
  
  # cleaf-<species> from *.veg_spec : fol_bio ( kg DM/ha) *1/20000
  for(veg_spec_file in veg_spec_files) {
    curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
    tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
    var_spec <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),10] / 20000 # convert from 
    df_ann[gsub("-<species>", paste("_", tmp_species_ISIMIP, sep=""), variable_opt_naming[11])] <- var_spec
  }
  
  # cwood-<species> from *.veg_spec : sap_bio + hrt_bio (kg DM/ha) *1/20000
  for(veg_spec_file in veg_spec_files) {
    curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
    tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
    var_spec <- (eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),11] + eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),13]) / 20000 # convert from 
    df_ann[gsub("-<species>", paste("_", tmp_species_ISIMIP, sep=""), variable_opt_naming[12])] <- var_spec
  }
  
  # croot-<species> from *.veg_spec : ??? *1/20000
  cvegbg_total <- NULL
  
  for(veg_spec_file in veg_spec_files) {
    curr_spec_fs <- strsplit(veg_spec_file, "\\/|\\_|\\.")[[1]][2]
    tmp_species_ISIMIP <- species_naming[mapply('%in%', curr_spec_fs, fs_species_naming)]
    alphac <- alphac_list[mapply('%in%', curr_spec_fs, fs_species_naming)]
    biomass <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),6]
    fol_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),10]
    sap_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),11]
    hrt_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),13]
    frt_bio <- eval(parse(text=veg_spec_file))[-length(eval(parse(text=veg_spec_file))[,1]),12]
    
    tbc_bio <- biomass - (fol_bio + sap_bio + hrt_bio + frt_bio)
    crt_bio <- alphac * tbc_bio
    vegbg <- frt_bio + crt_bio
    
    cvegbg <- vegbg /20000  # conversion from kg_DW/ha biomass into kgC/m²
    if(length(curr_exp_specs) > 1) {  # only do species wise extraction if multiple species were simulated
      df_ann[gsub("-<species>", paste("_", tmp_species_ISIMIP, sep=""), variable_opt_naming[13])] <- cvegbg
    }
    if(is.null(cvegbg_total)) {cvegbg_total <- cvegbg} else {cvegbg_total <- cvegbg_total + cvegbg} # add up species values to total
  }
  df_ann[gsub("-<species>", "_total", variable_opt_naming[13])] <- cvegbg_total
  
  # tsl from *_temp : temp_1, temp2,...  , for all depth layers : layers from soil.ini
  for (layer in 1:length(temp_file[1,c(-1,-2,-3)])) {
    df_day_depthclassed[, paste("soil_layer_", layer, sep="")] <- temp_file[,layer+3] + 273.15  # conversion from C to K
  }
  list <- as.list(as.data.frame(t(df_day_depthclassed[,-1])))
  list <- as.data.frame(matrix(list))
  df_day[paste(variable_opt_naming[14], sep="")] <- as.vector(list)
  
  
  ## replace all NaN with NA
  #df_day[df_day=="NaN"] <- NA    # TODO throws an error but probably not needed
  df_mon[df_mon=="NaN"] <- NA
  df_ann[df_ann=="NaN"] <- NA
  
  ###################################################################################
  ### 02c: write netCDF files for annual, monthly and daily variables
  ###################################################################################
  
  ## link 4c scenario-specifier naming to ISIMIP naming conventions
  tmp_gcm_ISIMIP <- gcm_naming[mapply('%in%', tmp_gcm, fs_gcm_naming)]
  tmp_bias_correction_ISIMIP <- bias_correction_naming[mapply('%in%', tmp_bias_correction, fs_bias_correction_naming)]
  tmp_climate_scenario_ISIMIP <- climate_scenario_naming[mapply('%in%', tmp_climate_scenario, fs_climate_scenario_naming)]
  tmp_soc_scenario_ISIMIP <- soc_scenario_naming[mapply('%in%', tmp_soc_scenario, fs_soc_scenario_naming)]
  tmp_co2sens_scenario_ISIMIP <- co2sens_scenario_naming[mapply('%in%', tmp_co2sens_scenario, fs_co2sens_scenario_naming)]
  tmp_region_ISIMIP <- region_naming[mapply('%in%', tmp_region, fs_region_naming)]
  
  # change soc-scenario name if socbe, sochwp, soca, socam to rcpXXsocXX (only in 2B)
  if(tmp_soc_scenario_ISIMIP %in% c("socbe", "sochwp", "soca", "socam") & tmp_climate_scenario_ISIMIP %in% c("picontrol", "rcp26", "rcp60", "rcp85")) {
    tmp_soc_scenario_ISIMIP <- paste(tmp_climate_scenario_ISIMIP, tmp_soc_scenario_ISIMIP, sep="")
  }

  # write netCDFs using an adapted version of Friedrich Bohns writeSim2netCDF() function as available in the ProfoundData R-package
  for (temp_res in list(list(df_ann, "annual"), list(df_mon, "monthly"), list(df_day, "daily"))) {
    writeSim2netCDF(df=temp_res[[1]],
                    comment1=paste("Data prepared for ISIMIP", ISIMIP_round, sep=""),
                    comment2=NA,
                    institution= 'Potsdam-Institute for Climate Impact Research (PIK)',
                    contact= 'Mats Mahnken (mahnken@pik-potsdam.de), Martin Gutsch (gutsch@pik-potsdam.de), Christopher Reyer (reyer@pik-potsdam.de)',
                    modelname="4c",
                    GCM=tmp_gcm_ISIMIP,
                    bc=tmp_bias_correction_ISIMIP,
                    RCP=tmp_climate_scenario_ISIMIP,
                    ses=tmp_soc_scenario_ISIMIP,
                    ss=tmp_co2sens_scenario_ISIMIP,
                    region=tmp_region_ISIMIP,
                    start=NA,
                    timestep=temp_res[[2]],
                    soildepths=ini_soil_horizons_depths,
                    folder=output,
		    ISIMIP_round=ISIMIP_round)
  }
  
  
  ## print progress message
  print(paste("finished", tmp_experiment, "of", length(experiments[,1]), "conversions"))
}



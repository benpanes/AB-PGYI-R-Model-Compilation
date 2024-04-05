#/***************************/
#Compile plot level attributes
#/***************************/
start_time <- Sys.time()

treelist <- fread("GYPSY data/intermediate/i_tree_list.csv")
regen <- fread("GYPSY data/intermediate/i_regeneration.csv")
plot <- fread("GYPSY data/intermediate/i_plot.csv")
plot_mm <- fread("GYPSY data/intermediate/i_plot_measurement.csv")

### Export regen sized trees for inclusion in  regen density calculations
# export regen sized trees for inclusion in regen density calculations;
tagged_regen <- treelist %>%
  mutate(regeneration_plot_name = tree_type) %>%
  filter(height < 1.3 & grepl("R", tree_type, fixed = TRUE)) # Output if height<1.3 and "R" in tree_type

# KJS additions for regen start here
tagged_regen_count <- tagged_regen %>%
  arrange(company, company_plot_number, measurement_number, species, regeneration_plot_name) %>%
  group_by(company, company_plot_number, measurement_number, species, regeneration_plot_name) %>%
  summarise(tagregcount = n(), .groups = "drop")

tally_regen <- regen %>%
  mutate(species = toupper(species)) %>%
  filter(!(species %in% c("MS"))) %>%
  mutate(species = case_when(
    species %in% c('AW', 'AX') ~ 'AW',
    species %in% c('LA', 'LT', 'LW') ~ 'LT',
    species %in% c('P', 'PX', 'PF', 'PL', 'PW') ~ 'PL',
    species %in% c('FA', 'FB') ~ 'FB',
    species %in% c('SE', 'SW', 'SX') ~ 'SW',
    TRUE ~ species
  ))

tagged_regen_count <- tagged_regen_count %>%
  arrange(company, company_plot_number, measurement_number, species, regeneration_plot_name)

tally_regen <- tally_regen %>%
  arrange(company, company_plot_number, measurement_number, species, regeneration_plot_name)

allregencounts <- merge(tagged_regen_count, tally_regen, by = c("company", "company_plot_number", "measurement_number", "species", "regeneration_plot_name"), all = TRUE)

# Replace missing values with 0
allregencounts$tagregcount[is.na(allregencounts$tagregcount)] <- 0
allregencounts$regeneration_count[is.na(allregencounts$regeneration_count)] <- 0

# Calculate totalregen
allregencounts$totalregen <- allregencounts$tagregcount + allregencounts$regeneration_count
allregencounts <- allregencounts %>%
  group_by(company, company_plot_number, measurement_number, species) %>%
  mutate(totalregencount = sum(totalregen))

totalregen <- allregencounts
regendensity <- merge(totalregen, plot_mm, by = c("company", "company_plot_number", "measurement_number"), all.x = TRUE)

# Calculate sphRegen
regendensity <- regendensity %>%
  mutate(
    SphRegen = totalregencount * 10000 / (regen_plot_area * number_regen_plots),
    spp_grp = case_when(
      species %in% c("AW", "PB", "BW", "AX") ~ "AW",
      species %in% c("PL", "PJ", "PW", "PF", "PX", "PA", "LT", "LA", "LW") ~ "PL",
      species %in% c("SW", "SE", "SX", "FB", "FA", "FD") ~ "SW",
      species == "SB" ~ "SB",
      TRUE ~ NA_character_
    ),
    condec = case_when(
      species %in% c("AW", "PB", "BW", "AX", "PL", "PJ", "PW", "PF", "PX", "PA", "LT", "LA", "LW") ~ "DEC",
      species %in% c("SW", "SE", "SX", "FB", "FA", "FD", "SB") ~ "CON",
      TRUE ~ NA_character_
    )
  ) %>%
  select(company, company_plot_number, measurement_number, species, spp_grp, condec, SphRegen)

# In compilation code, 3 df are created for sppgrp, condec, and sph
# Calculate SphRegen for each spp_grp
sppgrpregden <- regendensity %>%
  arrange(company, company_plot_number, measurement_number, spp_grp) %>%
  group_by(company, company_plot_number, measurement_number, spp_grp) %>%
  summarise(SphRegen = sum(SphRegen), .groups = "drop")

# Calculate SphRegen for each condec
condecregden <- regendensity %>%
  group_by(company, company_plot_number, measurement_number, condec) %>%
  summarise(SphRegen = sum(SphRegen), .groups = "drop")

# Calculate total SphRegen
totregden <- regendensity %>%
  group_by(company, company_plot_number, measurement_number) %>%
  summarise(SphRegen = sum(SphRegen), .groups = "drop")

################################################################################
# Sort sasin.tree_list1 data
# Create plot at tree level
plot <- treelist  %>%
  arrange(company, company_plot_number, measurement_number, species)%>%
  mutate(
    sph = if_else(is.na(sph), first(sph, na.rm = TRUE), sph),
    vol_0000ha = vol_0000 * sph,
    vol_1307ha = vol_1307 * sph,
    vol_1510ha = vol_1510 * sph,
    baha = ba * sph,
    biomassha = biomass * sph,
    carbonha = carbon * sph,
    sphBH = ifelse(height >= 1.3 | tree_type %in% c('S1', 'T', 'ET'), sph, NA_real_),
    sphD15 = ifelse(dbh >= 1.5, sph, NA_real_),
    sphD91 = ifelse(dbh >= 9.1, sph, NA_real_)
  ) %>%
  filter(!(height < 1.3 | grepl("R", tree_type)))  # Remove regen without the saplings falling back under 1.3m




################################################################################
# Create plot2 at spp level
plot2 <- plot %>%
  arrange(company, company_plot_number, measurement_number, species) %>%
  group_by(company, company_plot_number, measurement_number, species) %>%
  mutate(
    scale = "species"
  ) %>%
  summarise(
    sph = sum(sph, na.rm = TRUE),
    sphBH = sum(sphBH, na.rm = TRUE),
    sphD15 = sum(sphD15, na.rm = TRUE),
    sphD91 = sum(sphD91, na.rm = TRUE),
    ba = sum(ba, na.rm = TRUE),
    vol_0000 = sum(vol_0000, na.rm = TRUE),
    vol_1307 = sum(vol_1307, na.rm = TRUE),
    vol_1510 = sum(vol_1510, na.rm = TRUE),
    biomass = sum(biomass, na.rm = TRUE),
    carbon = sum(carbon, na.rm = TRUE),
    .groups = "drop"
  )

# Merge regendensity and plot2
plot2$scale <- "species"
plot2a <- plot2 %>%
  left_join(regendensity, by = c("company", "company_plot_number", "measurement_number", "species"))

# Replace missing values with 0
plot2a[is.na(plot2a$sph), "sph"] <- 0
plot2a[is.na(plot2a$sphBH), "sphBH"] <- 0
plot2a[is.na(plot2a$sphD15), "sphD15"] <- 0
plot2a[is.na(plot2a$sphD91), "sphD91"] <- 0
plot2a[is.na(plot2a$ba), "ba"] <- 0
plot2a[is.na(plot2a$vol_0000), "vol_0000"] <- 0
plot2a[is.na(plot2a$vol_1307), "vol_1307"] <- 0
plot2a[is.na(plot2a$vol_1510), "vol_1510"] <- 0
plot2a[is.na(plot2a$biomass), "biomass"] <- 0
plot2a[is.na(plot2a$carbon), "carbon"] <- 0
################################################################################
# Create plot3 at spp group level
plot3 <- plot %>%
  arrange(company, company_plot_number, measurement_number, spp_grp) %>%
  group_by(company, company_plot_number, measurement_number, spp_grp) %>%
  summarise(
    sph = sum(sph, na.rm = TRUE),
    sphBH = sum(sphBH, na.rm = TRUE),
    sphD15 = sum(sphD15, na.rm = TRUE),
    sphD91 = sum(sphD91, na.rm = TRUE),
    ba = sum(ba),
    vol_0000 = sum(vol_0000, na.rm = TRUE),
    vol_1307 = sum(vol_1307, na.rm = TRUE),
    vol_1510 = sum(vol_1510, na.rm = TRUE),
    biomass = sum(biomass, na.rm = TRUE),
    carbon = sum(carbon, na.rm = TRUE),
    .groups = "drop"
  )

plot3$scale <- 'spp_grp'
plot3$species <- plot3$spp_grp
plot3a <- plot3 %>%
  left_join(sppgrpregden, by = c("company", "company_plot_number", "measurement_number", "spp_grp"))
################################################################################
# Create plot4 at condec group level
plot4 <- plot %>%
  arrange(company, company_plot_number, measurement_number, condec) %>%
  group_by(company, company_plot_number, measurement_number, condec) %>%
  summarise(
    sph = sum(sph, na.rm = TRUE),
    sphBH = sum(sphBH, na.rm = TRUE),
    sphD15 = sum(sphD15, na.rm = TRUE),
    sphD91 = sum(sphD91, na.rm = TRUE),
    ba = sum(ba, na.rm = TRUE),
    vol_0000 = sum(vol_0000, na.rm = TRUE),
    vol_1307 = sum(vol_1307, na.rm = TRUE),
    vol_1510 = sum(vol_1510, na.rm = TRUE),
    biomass = sum(biomass, na.rm = TRUE),
    carbon = sum(carbon, na.rm = TRUE),
    .groups = "drop"
  )

plot4 <- plot4 %>%
  mutate(
    scale = "con_dec",
    species = if_else(condec == "DEC", "DE", "CO")
  )

plot4a <- plot4 %>%
  left_join(condecregden, by = c("company", "company_plot_number", "measurement_number", "condec"))
################################################################################
# Create plot5 at total level
plot5 <- plot %>%
  arrange(company, company_plot_number, measurement_number, condec) %>%
  group_by(company, company_plot_number, measurement_number) %>%
  summarise(
    sph = sum(sph, na.rm = TRUE),
    sphBH = sum(sphBH, na.rm = TRUE),
    sphD15 = sum(sphD15, na.rm = TRUE),
    sphD91 = sum(sphD91, na.rm = TRUE),
    ba = sum(ba, na.rm = TRUE),
    vol_0000 = sum(vol_0000, na.rm = TRUE),
    vol_1307 = sum(vol_1307, na.rm = TRUE),
    vol_1510 = sum(vol_1510, na.rm = TRUE),
    biomass = sum(biomass, na.rm = TRUE),
    carbon = sum(carbon, na.rm = TRUE),
    .groups = "drop"
  )

# Add species and scale columns
plot5$species <- 'TO'
plot5$scale <- 'con_dec'

plot5a <- plot5 %>%
  left_join(totalregen, by = c("company", "company_plot_number", "measurement_number", "species"))
################################################################################
# combine

plot2b <- bind_rows(plot2a, plot3a)
plot2b <- plot2b %>%
  arrange(company, company_plot_number, measurement_number, scale, species)

plot2c <- bind_rows(plot2b, plot4a)
plot2c <- plot2c %>%
  arrange(company, company_plot_number, measurement_number, scale, species)

plot2d <- bind_rows(plot2c, plot5a)
plot2d <- plot2d %>%
  arrange(company, company_plot_number, measurement_number, scale, species)

##############################################################################
# Sort plot2a data
cols_to_drop <- c("company_stand_number", "establishment_month", "establishment_day", "measurement_day", 
                  "contractor", "cruiser_1_name", "cruiser_2_name", "plot_comment", "plot_measurement_comment", 
                  "avi_field_call", "shrub_cover", "herb_forb_cover", "grass_cover", "moss_lichen_cover",
                  "tree_plot_shape", "sapling_plot_shape", "regen_plot_shape", "plot_status", 
                  "ats_township", "ats_range", "ats_meridian", "ats_section")

plot <- fread("H:/Shared drives/Model Comparison/Compilation_Code_R/PGYI Compiled/i_plot.csv")
plot_mm <- fread("H:/Shared drives/Model Comparison/Compilation_Code_R/PGYI Compiled/i_plot_measurement.csv")

mmt <- plot_mm %>%
  left_join(plot, by = c("company", "company_plot_number")) %>%
  select(-any_of(cols_to_drop))

################################################################################
# Create mmt2 dataframe
mmt2 <- plot2d %>%
  select(
    company, company_plot_number, measurement_number,species, scale
  ) %>%
  left_join(mmt, by = c("company", "company_plot_number", "measurement_number")) 


# Merge mmt2 and plot3 dataframes
plot44 <- merge(mmt2, plot2d, by = c("company", "company_plot_number", "measurement_number", "scale", "species"))

################################################################################
plot6 <- plot44 %>%
  mutate(regendone = ifelse(!is.na(SphRegen), 1, 0)) %>%
  filter(species == "TO") %>%
  select(company, company_plot_number, measurement_number, regendone)

plot_level1 <- plot44 %>%
  left_join(plot6, by = c("company", "company_plot_number", "measurement_number"))

# Replace missing values with 0 
plot_level1 <- plot_level1 %>%
  mutate(
    sph = ifelse(is.na(sph), 0, sph),
    sphBH = ifelse(is.na(sphBH), 0, sphBH),
    sphD15 = ifelse(is.na(sphD15), 0, sphD15),
    sphD91 = ifelse(is.na(sphD91), 0, sphD91),
    ba = ifelse(is.na(ba), 0, ba),
    vol_0000 = ifelse(is.na(vol_0000), 0, vol_0000),
    vol_1307 = ifelse(is.na(vol_1307), 0, vol_1307),
    vol_1510 = ifelse(is.na(vol_1510), 0, vol_1510),
    biomass = ifelse(is.na(biomass), 0, biomass),
    carbon = ifelse(is.na(carbon), 0, carbon),
    SphRegen = ifelse(regendone == 1 & is.na(SphRegen), 0, SphRegen)
  )

# Remove rows where species is "NO" or blank
plot_level1 <- plot_level1 %>%
  filter(!(species %in% c("NO", ""))) %>%
  select(
    company,	company_plot_number,	plot_program,	establishment_year,	fmu,	fma,	opening_number,
    sampling_unit_number,	topographic_position,	elevation,	slope,	aspect,	x_coord,	y_coord,	utm_zone,	
    datum,	latitude,	longitude,	natural_subregion,	ecosite_guide,	ecosite,	ecosite_phase,	shared,	measurement_number,	
    measurement_year,	measurement_month,	stand_origin,	plot_type,	stand_type,	tree_plot_area,	tree_tagging_limit,	
    sapling_plot_area,	sapling_tagging_limit_dbh,	sapling_tagging_limit_height,	number_sapling_plots,	regen_plot_area,	
    regen_tagging_limit_conifer,	regen_tagging_limit_decid,	number_regen_plots,	scale,	species,	sph,	sphBH,	sphD15,	
    sphD91,	ba,	vol_0000,	vol_1307,	vol_1510,	biomass,	carbon,	SphRegen,	regendone
  )

fwrite("GYPSY data/intermediate/i_plot_level1.csv")

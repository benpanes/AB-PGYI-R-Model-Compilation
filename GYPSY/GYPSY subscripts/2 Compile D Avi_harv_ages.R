#/***************************/
#Compile plot level attributes
#/***************************/

library(tidyr)
start_time <- Sys.time()

photo.avi <- fread("GYPSY data/intermediate/i_photo_avi.csv")
plot <- fread("GYPSY data/intermediate/i_plot.csv")
plot_mmt <- fread("GYPSY data/intermediate/i_plot_measurement.csv")
treatment <- fread("GYPSY data/intermediate/i_treatment.csv")
disturbance <- fread("GYPSY data/intermediate/i_disturbance.csv")
plot_level1 <- fread("GYPSY data/intermediate/i_plot_level1.csv")
  
################################################################################
### STEP 1. append avi-based origin year to plot summary table;
################################################################################

# take species info from layer 1
photo1 <- photo.avi %>%
  filter(layer_type == 1) %>%
  select(company, company_plot_number, avi_version, polygon_number, year_photography, year_photo_call, 
         moist_reg, density, height, sp1, sp1_per, sp2, sp2_per, sp3, sp3_per, sp4, sp4_per, 
         sp5, sp5_per, struc, struc_val, origin, tpr, initials, nfl, nfl_per, nat_non, anth_veg, 
         anth_non, mod1, mod1_ext, mod1_yr, mod2, mod2_ext, mod2_yr, data, data_yr) %>%
  arrange(company, company_plot_number, avi_version, polygon_number, year_photography, year_photo_call)


# take species info from layer 2
photo2 <- photo.avi %>%
  filter(layer_type == 2) %>%
  select(-layer_type) %>%
  arrange(company, company_plot_number, avi_version, polygon_number, year_photography, year_photo_call)

photo2a <- photo2 %>%
  rename(
    umoist_reg = moist_reg,
    udensity = density,
    uheight = height,
    usp1 = sp1,
    usp1_per = sp1_per,
    usp2 = sp2,
    usp2_per = sp2_per,
    usp3 = sp3,
    usp3_per = sp3_per,
    usp4 = sp4,
    usp4_per = sp4_per,
    usp5 = sp5,
    usp5_per = sp5_per,
    ustruc = struc,
    uorigin = origin,
    utpr = tpr,
    uinitials = initials,
    unfl = nfl,
    unfl_per = nfl_per,
    unat_non = nat_non,
    uanth_veg = anth_veg,
    uanth_non = anth_non,
    umod1 = mod1,
    umod_ext = mod1_ext,
    umod_yr = mod1_yr,
    umod2 = mod2,
    umod2_ext = mod2_ext,
    umod2_yr = mod2_yr,
    udata = data,
    udata_yr = data_yr
  ) %>%
  select(company, company_plot_number, avi_version, polygon_number, year_photography, year_photo_call, 
         umoist_reg, udensity, uheight, usp1, usp1_per, usp2, usp2_per, usp3, usp3_per, usp4, usp4_per, 
         usp5, usp5_per, ustruc, struc_val, uorigin, utpr, uinitials, unfl, unfl_per, unat_non, uanth_veg, 
         uanth_non, umod1, umod_ext, umod_yr, umod2, umod2_ext, umod2_yr, udata, udata_yr)  %>%
  arrange(company, company_plot_number, avi_version, polygon_number, year_photography, year_photo_call)

# take nonforested and modifer info from both layers and remove duplicates
photo3 <- photo.avi %>%
  filter(layer_type == 1) %>%  # Select rows where layer_type is 1
  select(company, company_plot_number, avi_version, polygon_number, year_photography, year_photo_call, photo_avi_layer_comment) %>%
  arrange(company, company_plot_number, avi_version, polygon_number, year_photography, year_photo_call)

# merge photo1 and photo2a to fill missing info
photo4 <- merge(photo1, photo2a, by = c("company", "company_plot_number", "avi_version", "polygon_number", "year_photography", "year_photo_call"), all.x = TRUE)

# Merge photo4 with photo3 dataset
photo5 <- merge(photo4, photo3, by = c("company", "company_plot_number", "avi_version", "polygon_number", "year_photography", "year_photo_call"))


photo6 <- photo5 %>%
  filter(!(year_photography == 0 | is.na(year_photography))) %>%
  filter(!(origin == 0 | is.na(origin)) | !(mod1_yr == 0 | is.na(mod1_yr)) | !(uorigin == 0 | is.na(uorigin)))%>%
  mutate(
    avi_origin = case_when(
      origin %in% c(0, NA) & mod1_yr > 0 ~ mod1_yr,
      origin > 0 & mod1_yr %in% c(0, NA) ~ origin,
      origin > 0 & mod1_yr > 0 & mod1_ext %in% c(NA, 0, 1, 2, 3) ~ origin,
      origin > 0 & mod1_yr > 0 & mod1_ext %in% c(4, 5) ~ mod1_yr,
      TRUE ~ NA_real_
    )
  ) %>%
  select(company, company_plot_number, year_photography, avi_origin, mod1, mod1_ext) %>%
  arrange(company, company_plot_number, year_photography)

photo7 <- photo6 %>%
  group_by(company, company_plot_number) %>%
  mutate(avi_origin1 = avi_origin[1],
         avi_origin2 = ifelse(n() == 2, avi_origin[2], NA_real_),
         avi_origin3 = ifelse(n() == 3, avi_origin[3], NA_real_)) %>%
  distinct(company, company_plot_number, .keep_all = TRUE) %>%
  select(company, company_plot_number, avi_origin1, avi_origin2, avi_origin3)

photo8 <- photo6 %>%
  group_by(company, company_plot_number) %>%
  mutate(photo_year1 = year_photography[1],
         photo_year2 = ifelse(n() == 2, year_photography[2], NA_real_),
         photo_year3 = ifelse(n() == 3, year_photography[3], NA_real_)) %>%
  distinct(company, company_plot_number, .keep_all = TRUE) %>%
  select(company, company_plot_number, photo_year1, photo_year2, photo_year3)
  
photo6$mod1 <- as.character(photo6$mod1)
photo9 <- photo6 %>%
  group_by(company, company_plot_number) %>%
  mutate(mod1_1 = mod1[1],
         mod1_2 = ifelse(n() == 2, as.character(mod1[2]), NA_character_),
         mod1_3 = ifelse(n() == 3, as.character(mod1[3]), NA_character_)) %>%
  distinct(company, company_plot_number, .keep_all = TRUE) %>%
  select(company, company_plot_number, mod1_1, mod1_2, mod1_3)


photo10 <- photo6 %>%
  group_by(company, company_plot_number) %>%
  mutate(mod1_ext_1 = mod1_ext[1],
         mod1_ext_2 = ifelse(n() == 2, mod1_ext[2], NA_real_),
         mod1_ext_3 = ifelse(n() == 3, mod1_ext[3], NA_real_)) %>%
  distinct(company, company_plot_number, .keep_all = TRUE) %>%
  select(company, company_plot_number, mod1_ext_1, mod1_ext_2, mod1_ext_3)

photo11 <- inner_join(photo7, photo8, by = c("company", "company_plot_number"))

photo12 <- inner_join(photo11, photo9, by = c("company", "company_plot_number"))

photo13 <- inner_join(photo12, photo10, by = c("company", "company_plot_number")) %>%
  arrange(company, company_plot_number)

# lost 3 rows comparing with SAS code
mmt <- plot_mmt %>%
  left_join(plot, by = c("company", "company_plot_number")) %>%
  select(company, company_plot_number, measurement_year, stand_origin) %>%
  distinct(company, company_plot_number, measurement_year, .keep_all = TRUE) %>%
  arrange(company, company_plot_number)

mmt2 <- left_join(mmt, photo13, by = c("company", "company_plot_number"))

# same rules applied for managed and natural stands, why test if avi_origin is 'N' or not?
mmt3 <- mmt2 %>%
  mutate(avi_origin = ifelse(!is.na(avi_origin2),
                             ifelse(abs(avi_origin1 - avi_origin2) < 11, avi_origin2,
                                    ifelse(stand_origin == 'N',
                                           ifelse(!(mod1_ext_2 %in% c(4, 5)), avi_origin2,
                                                  ifelse(measurement_year > photo_year2, avi_origin2, avi_origin1)),
                                           ifelse(!(mod1_ext_2 %in% c(4, 5)), avi_origin2,
                                                  ifelse(measurement_year > photo_year2, avi_origin2, avi_origin1)))),
                             ifelse(!is.na(avi_origin1), avi_origin1, NA_real_))) %>%
  mutate(age_avi = measurement_year - avi_origin) %>%
  select(-avi_origin1, -avi_origin2, -photo_year1, -photo_year2, -mod1_ext_1, -mod1_ext_2, -mod1_1, -mod1_2, -avi_origin, -stand_origin)

################################################################################
# STEP 2.calculate age based on year of harvest
################################################################################

harv <- treatment %>%
  filter(
    treatment_code == 'H'
  )

harv2 <- harv %>%
  mutate(harv_orig = treatment_year) %>%
  select(company, company_plot_number, harv_orig)

mmt3 <- arrange(mmt3, company, company_plot_number)
harv2 <- arrange(harv2, company, company_plot_number)

# Merge mmt3 and harv2 data frames
harv3 <- mmt3 %>%
  left_join(harv2, by = c("company", "company_plot_number")) %>%
  distinct(company, company_plot_number, measurement_year, .keep_all = TRUE)

# Calculate age_harv
harv4 <- harv3 %>%
  mutate(age_harv = measurement_year - harv_orig) %>%
  select(-harv_orig)

###
fire <- disturbance %>%
  mutate(fire_orig = disturbance_year) %>%
  filter(disturbance_code == "DF") %>%
  select(company, company_plot_number, fire_orig)%>%
  arrange(fire_orig, company, company_plot_number)

fireharv <- merge(harv4, fire, by = c("company", "company_plot_number"), all.x = TRUE)

fireharv2 <- fireharv %>%
  mutate(age_fire = measurement_year - fire_orig) %>%
  select(-fire_orig)

# Sort fireharv2 dataset
fireharv2 <- arrange(fireharv2, company, company_plot_number, measurement_year)

# Merge sasin.plot_level1 and fireharv2 datasets
plot_level2 <- merge(plot_level1, fireharv2, by = c("company", "company_plot_number", "measurement_year"))

fwrite(plot_level2, "GYPSY data/intermediate/i_plot_level2.csv")

pgyi_compiled <- fread("GYPSY data/intermediate/i_pgyi_compiled.csv")

################################################################################
# set up to determine compositional strata by BA - note species level ba for Aw, Sw, Sb left blank so it doesn't double-count;

ba <- pgyi_compiled %>%
  mutate(
    totba = ifelse(species == "TO", ba, NA),
    totregden = ifelse(species == "TO", sphRegenOnly, NA),
    conba = ifelse(species == "CO", ba, NA),
    conregden = ifelse(species == "CO", sphRegenOnly, NA),
    decba = ifelse(species == "DE", ba, NA),
    decregden = ifelse(species == "DE", sphRegenOnly, NA),
    pinba = ifelse(scale == "spp_grp" & species == "PL", ba, NA),
    pinregden = ifelse(scale == "spp_grp" & species == "PL", sphRegenOnly, NA),
    sfba = ifelse(scale == "spp_grp" & species == "SW", ba, NA),
    sfregden = ifelse(scale == "spp_grp" & species == "SW", sphRegenOnly, NA),
    sbba = ifelse(scale == "spp_grp" & species == "SB", ba, NA),
    sbregden = ifelse(scale == "spp_grp" & species == "SB", sphRegenOnly, NA),
    pjba = ifelse(species == "PJ", ba, NA),
    pjregden = ifelse(species == "PJ", sphRegenOnly, NA),
    firba = ifelse(species %in% c("FB", "FA", "FD"), ba, NA),
    firregden = ifelse(species %in% c("FB", "FA", "FD"), sphRegenOnly, NA),
    bwba = ifelse(species == "BW", ba, NA),
    bwregden = ifelse(species == "BW", sphRegenOnly, NA),
    pbba = ifelse(species == "PB", ba, NA),
    pbregden = ifelse(species == "PB", sphRegenOnly, NA),
    co_plot = paste(company, plot, sep = "")
  ) 

plotba <- ba %>%
  group_by(company, plot, mmt_num) %>%
  summarise(
    age = first(age, na.rm = TRUE),
    type = first(stand_type, na.rm = TRUE),                 # what is "type" here?
    totba = sum(totba, na.rm = TRUE),
    conba = sum(conba, na.rm = TRUE),
    decba = sum(decba, na.rm = TRUE),
    pinba = sum(pinba, na.rm = TRUE),
    sfba = sum(sfba, na.rm = TRUE),
    sbba = sum(sbba, na.rm = TRUE),
    pjba = sum(pjba, na.rm = TRUE),
    firba = sum(firba, na.rm = TRUE),
    bwba = sum(bwba, na.rm = TRUE),
    pbba = sum(pbba, na.rm = TRUE),
    totregden = sum(totregden, na.rm = TRUE),
    conregden = sum(conregden, na.rm = TRUE),
    decregden = sum(decregden, na.rm = TRUE),
    pinregden = sum(pinregden, na.rm = TRUE),
    sfregden = sum(sfregden, na.rm = TRUE),
    sbregden = sum(sbregden, na.rm = TRUE),
    pjregden = sum(pjregden, na.rm = TRUE),
    bwregden = sum(bwregden, na.rm = TRUE),
    firregden = sum(firregden, na.rm = TRUE),
    pbregden = sum(pbregden, na.rm = TRUE),
    co_plot = first(co_plot, na.rm = TRUE)
  ) %>%
  ungroup()%>%
  arrange(company, plot, mmt_num)


# Replace missing values or "." with "0" for columns ending with "ba" or "regen"
# but no need

compos <- plotba %>%
  mutate(
  strata1 = if_else(
    totba > 0,
    if_else(decba / totba > 0.5, "Hw",
            if_else(pinba > (sfba + sbba), "Pl",
                    if_else(sfba >= sbba, "Sw",
                            if_else(sfba < sbba, "Sb", " ")))), ""),
  strata2 = if_else(
    totba > 0,
    if_else(decba / totba > 0.5 & decba / totba < 0.8 & sfba >= sbba, "Sw",
            if_else(decba / totba > 0.5 & decba / totba < 0.8 & sfba < sbba, "Sb",
                    if_else(decba / totba > 0.5 & decba / totba < 0.8 & (sfba + sbba) > 0 & pinba>(sfba+sbba), "Pl",
                                    if_else(decba / totba <= 0.5 & conba / totba < 0.8, "Hw", "")))), ""),
  substrata = if_else(
    totba > 0,
    if_else(decba / totba > 0.5 & bwba / decba > 0.5, "Bw",
            if_else(decba / totba > 0.5 & pbba / decba > 0.5, "Pb",
                    if_else(decba / totba <= 0.5 & pinba > 0 & pjba / pinba > 0.5, "Pj", 
                            if_else(decba / totba <= 0.5 & firba / conba > 0.5, "Fb", "")))), ""))%>%
  mutate (
    strata = paste0(strata1, strata2)
  )

compos_regen <- compos %>%
  filter(
    strata == ""
  ) %>%
  mutate(
    strata1 = if_else(
     totregden > 0,
      if_else(decregden / totregden > 0.5, "Hw",
              if_else(pinregden > (sfregden + sbregden), "Pl",
                      if_else(sfregden >= sbregden, "Sw",
                              if_else(sfregden < sbregden, "Sb", " ")))), ""),
    strata2 = if_else(
      totregden > 0,
      if_else(decregden / totregden > 0.5 & decregden / totregden < 0.8 & sfregden >= sbregden, "Sw",
              if_else(decregden / totregden > 0.5 & decregden / totregden < 0.8 & sfregden < sbregden, "Sb",
                      if_else(decregden / totregden > 0.5 & decregden / totregden < 0.8 & (sfregden + sbregden) > 0 & pinregden>(sfregden+sbregden), "Pl",
                              if_else(decregden / totregden <= 0.5 & conregden / totregden < 0.8, "Hw", "")))), ""),
    substrata = if_else(
      totregden > 0,
      if_else(decregden / totregden > 0.5 & bwregden / decregden > 0.5, "Bw",
              if_else(decregden / totregden > 0.5 & pbregden / decregden > 0.5, "Pb",
                      if_else(decregden / totregden <= 0.5 & pinregden > 0 & pjregden / pinregden > 0.5, "Pj", 
                              if_else(decregden / totregden <= 0.5 & firregden / conregden > 0.5, "Fb", "")))), ""))%>%
  mutate (
    strata = paste0(strata1, strata2)
  )

compos_nor <- compos %>%
  filter(
    strata != ""
  )

compos_2 <- bind_rows(compos_nor, compos_regen)

strata2 <- compos_2 %>%
  mutate(yrfrom50 = abs(age - 50)) %>%
  arrange(co_plot, yrfrom50) %>%
  group_by(co_plot) %>%
  slice(1) %>%
  select(
    co_plot, strata, substrata
  )

pgyi.compiledstratified <- ba %>%
  left_join(strata2, by = "co_plot") %>%
  select(-c(co_plot, totba, totregden, conba, conregden, decba, decregden, pinba, pinregden, sfba, sfregden, sbba, sbregden, pjba, pjregden, firba, firregden, bwba, pbba, pbregden, bwregden)) 

fwrite(pgyi.compiledstratified, "GYPSY data/final/PGYIcompiledstratified.csv")

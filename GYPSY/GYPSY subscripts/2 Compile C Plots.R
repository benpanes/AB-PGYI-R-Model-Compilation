#/***************************/
#Compile plot level attributes
#/***************************/
start_time <- Sys.time()

treelist <- fread("GYPSY data/intermediate/i_tree_list.csv")
regen <- fread("GYPSY data/intermediate/i_regeneration.csv")
plot <- fread("GYPSY data/intermediate/i_plot.csv")
plot_mm <- fread("GYPSY data/intermediate/i_plot_measurement.csv")

### Export regen sized trees for inclusion in  regen density calculations

tagged_regen <- # height < 1.3 and tree_type = R1-R10









fwrite("GYPSY data/intermediate/i_plot_level1.csv")
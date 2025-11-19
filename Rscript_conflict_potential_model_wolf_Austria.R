# *************************************************************
# Mapping Conflict Potential using Socio-Economic and Ecological Factors to Display conflicts about wolves in Austria
# *************************************************************
# version: 18.11.25
# r code authored by Fabian Knufinke, Jennifer Hatlauf, Florian Kunz
# script corresponds to publication: XXX [will be added as soon as accepted]

# How to cite: [will be added as soon as accepted]
  
require("terra")
require("tidyverse")
require("ape")
require("BAMMtools")
require("ggplot2")
require("ggspatial")

# 1) Read in input data - 17 variables
# **************************
path <- "C:\\Users\\" #set your path

{
  # Bevoelkerungsdichte
  bevoelk <- terra::rast(paste0(path, "popdenskm2_austria_100m_3416.tif"))
  # Ausgleichszahlungen
  ausgleich <- terra::rast(paste0(path, "compensation_austria_100m_3416.tif")) 
  # Hunde_allgemein
  hunde_a <- terra::rast(paste0(path, "hundekm2_austria_100m_3416.tif"))
  # Hunde_Jagd
  hunde_j <- terra::rast(paste0(path, "jagdhunde_austria_100m_3416.tif"))
  # Übernachtungszahlen Tourismus Winter
  nacht_wi <- terra::rast(paste0(path, "winternightskm2_austria_100m_3416.tif"))
  # Übernachtungszahlen Tourismus Sommer
  nacht_so <- terra::rast(paste0(path, "summernightskm2_austria_100m_3416.tif"))
  # Schalenwild_Abschuss
  schalenw_abs <- terra::rast(paste0(path, "DATEN/final_data/Schalenwilddichte/Dichte_Schalenwild_pro_100ha.tif"))
  # Schalenwild_Anzahl
  schalenw_anz <- terra::rast(paste0(path, "DATEN/final_data/Schalenwildartenanzahl/ungspecies_austria_100m_3416.tif"))
  # Wildruhezonen
  ruhe <- terra::rast(paste0(path, "Wildruhezonen_Flaechen.tif"))
    #Rotwildfütterungen
  fuetterung <- terra::rast(paste0(path, "Rotwildfütterungen_Flaechen.tif"))
  # Wildökologische Raumplanung
  woerp <- terra::rast(paste0(path, "WOERP_Flaechen.tif"))
  # Rissanfälligkeit
  riss_mod <- terra::rast("final_DRM.tif") 
  # set crs for le
  crs(riss_mod) <- "epsg:3416"
  # Herdenschutz Förderungszahlungen
  herden_s <- terra::rast(paste0(path, "protection_austria_100m_3416.tif"))
  # Nutztierdichte
  nutztier_d <- terra::rast(paste0(path, "ATfinal_Schafe_Rinder.tif"))
  # Gefährdete Nutztierrasse
  nutztier_r <- terra::rast(paste0(path, "gefnutzbinary_austria_100m_3416.tif"))
  # Gatterwild
  gatterw <- terra::rast(paste0(path, "zuchtwild_austria_100m_3416.tif"))
  # Schutzwald
  wald <- terra::rast(paste0(path, "Schutzwald_Flaechen_bml.tif")) 
  # Lebensraumpotenzial
  le_wolf <- terra::rast("final_HSI_avg.tif")  
  crs(le_wolf) <- "epsg:3416"
}

# 2) Variable reclassification
# **************************
# read in classification table
clas_table_Behörde <- read.csv2(paste0(path, "LEKO_classification_table_nur_Behörden.csv"), dec=",", stringsAsFactors = FALSE)
clas_table <- clas_table_Behörde 

# 2.1) bevoelk: con -> con
# **************************
values_bevoelk<- clas_table  %>%
  filter(variable == "bevoelk") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} %>% 
  {.[4,3] <- 10; .}

fun.reclassify <- splinefun(values_bevoelk$point, values_bevoelk$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}

pd_bevoelk <- terra::app(bevoelk, fun)

rm(list= c("values_bevoelk", "fun", "fun.reclassify"))

# 2.2) ausgleich: con -> con
# **************************
values_ausgleich <- clas_table  %>%
  filter(variable == "ausgleich") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} 

fun.reclassify <- splinefun(values_ausgleich$point, values_ausgleich$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}
pd_ausgleich <- terra::app(ausgleich, fun)

rm(list= c("values_ausgleich", "fun", "fun.reclassify"))

# 2.3) hunde_a: con -> con
# **************************
values_hunde_a <- clas_table  %>%
  filter(variable == "hunde_a") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} 

fun.reclassify <- splinefun(values_hunde_a$point, values_hunde_a$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}

pd_hunde_a <- terra::app(hunde_a, fun)

rm(list= c("values_hunde_a", "fun", "fun.reclassify"))

  # 2.4) hunde_j: con -> con
# **************************
values_hunde_j <- clas_table  %>%
  filter(variable == "hunde_j") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} 

fun.reclassify <- splinefun(values_hunde_j$point, values_hunde_j$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}

pd_hunde_j <- terra::app(hunde_j, fun)

rm(list= c("values_hunde_j", "fun", "fun.reclassify"))

  # 2.5a) nacht_wi: con -> con
# **************************
values_nacht_wi <- clas_table  %>%
  filter(variable == "nacht_wi") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} 

fun.reclassify <- splinefun(values_nacht_wi$point, values_nacht_wi$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}

pd_nacht_wi <- terra::app(nacht_wi, fun)

rm(list= c("values_nacht_wi", "fun", "fun.reclassify"))

  # 2.5b) nacht_so: con -> con
# **************************
values_nacht_so <- clas_table  %>%
  filter(variable == "nacht_so") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} 

fun.reclassify <- splinefun(values_nacht_so$point, values_nacht_so$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}

pd_nacht_so <- terra::app(nacht_so, fun)

rm(list= c("values_nacht_so", "fun", "fun.reclassify"))

# 2.6) schalenw_abs: con -> con
# **************************
values_schalenw_abs <- clas_table  %>%
  filter(variable == "schalenw_abs") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} 

fun.reclassify <- splinefun(values_schalenw_abs$point, values_schalenw_abs$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}

pd_schalenw_abs <- terra::app(schalenw_abs, fun)

rm(list= c("values_schalenw_abs", "fun", "fun.reclassify"))

# 2.7) schalenw_anz: con -> con
# **************************
values_schalenw_anz <- clas_table  %>%
  filter(variable == "schalenw_anz") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} 

fun.reclassify <- splinefun(values_schalenw_anz$point, values_schalenw_anz$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}

pd_schalenw_anz <- terra::app(schalenw_anz, fun)

rm(list= c("schalenw_anz", "fun", "fun.reclassify"))

# 2.8) ruhe: cat -> cat
# **************************
values_ruhe <- clas_table %>% 
  filter(variable == "ruhe") %>% 
  select(c("from", "value")) %>% 
  rename("is" = "from", "becomes" = "value")

pd_ruhe <- terra::classify(ruhe, as.matrix(values_ruhe))

rm(list = c("values_ruhe", "ruhe"))

# 2.9) fuetterung: cat -> cat
# **************************
values_fuetterung <- clas_table %>% 
  filter(variable == "fuetterung") %>% 
  select(c("from", "value")) %>% 
  rename("is" = "from", "becomes" = "value")

pd_fuetterung <- terra::classify(fuetterung, as.matrix(values_fuetterung))

rm(list = c("values_fuetterung", "fuetterung"))

  # 2.10) woerp: cat -> cat
# **************************
values_woerp <- clas_table %>% 
  filter(variable == "woerp") %>% 
  select(c("from", "value")) %>% 
  rename("is" = "from", "becomes" = "value")

pd_woerp <- terra::classify(woerp, as.matrix(values_woerp))

rm(list = c("values_woerp", "woerp"))

  # 2.11) riss_mod: con -> con
# **************************
values_riss_mod <- clas_table  %>%
  filter(variable == "riss_mod") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} 

fun.reclassify <- splinefun(values_riss_mod$point, values_riss_mod$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}

pd_riss_mod <- terra::app(riss_mod, fun)

rm(list= c("riss_mod", "fun", "fun.reclassify"))

 # 2.13) herden_s: con -> con
# **************************
values_herden_s <- clas_table  %>%
  filter(variable == "herden_s") %>% 
  select(c("from", "to", "value")) %>% 
  mutate(point = (from+to)/2) %>% 
  relocate(from, to, point, value) %>% 
  add_row(slice_head(.), .before = 1) %>% 
  add_row(slice_tail(.), .after = nrow(.)) %>% 
  {.[1,3] <- .[1,3]/2; .} %>%
  {.[1,4] <- .[1,4]; .} %>%
  {.[nrow(.),3] <- .[nrow(.),3]*1.1; .} %>%
  {.[nrow(.),4] <- .[nrow(.),4]; .} 

fun.reclassify <- splinefun(values_herden_s$point, values_herden_s$value, method = "monoH.FC")

fun <- function(x) {
  x[!is.na(x)] <- fun.reclassify(x[!is.na(x)])
  return(x)
}

pd_herden_s <- terra::app(herden_s, fun)

rm(list= c("herden_s", "fun", "fun.reclassify"))

# 2.14) nutztier_d: con -> con
# **************************
values_nutztier_d <- clas_table %>% 
  filter(variable == "nutztier_d") %>% 
  select(c("from", "value")) %>% 
  rename("is" = "from", "becomes" = "value")

pd_nutztier_d <- terra::classify(nutztier_d, as.matrix(values_nutztier_d))

rm(list = c("values_nutztier_d"))

# 2.15) nutztier_r: cat -> cat
# **************************
values_nutztier_r <- clas_table %>% 
  filter(variable == "nutztier_r") %>% 
  select(c("from", "value")) %>% 
  rename("is" = "from", "becomes" = "value")

pd_nutztier_r <- terra::classify(nutztier_r, as.matrix(values_nutztier_r))

rm(list = c("values_nutztier_r"))

# 2.16) gatterw: con -> con
# **************************
values_gatterw <- clas_table %>% 
  filter(variable == "gatterw") %>% 
  select(c("from", "value")) %>% 
  rename("is" = "from", "becomes" = "value")

pd_gatterw <- terra::classify(gatterw, as.matrix(values_gatterw))

rm(list = c("values_gatterw"))

# 2.17) wald: cat -> cat
# **************************
values_wald <- clas_table %>% 
  filter(variable == "wald") %>% 
  select(c("from", "value")) %>% 
  rename("is" = "from", "becomes" = "value")

pd_wald <- terra::classify(wald, as.matrix(values_wald))

rm(list = c("values_wald"))

# 3) Calculation of Model and Scenario setup
# **************************

# read in weights
{w_bevoelk <- clas_table$weight[clas_table$variable=="bevoelk"][1]
  w_ausgleich <- clas_table$weight[clas_table$variable=="ausgleich"][1]
  w_hunde_a <- clas_table$weight[clas_table$variable=="hunde_a"][1]
  w_hunde_j <- clas_table$weight[clas_table$variable=="hunde_j"][1]
  w_nacht_wi <- clas_table$weight[clas_table$variable=="nacht_wi"][1]
  w_nacht_so <- clas_table$weight[clas_table$variable=="nacht_so"][1]
  w_schalenw_abs <- clas_table$weight[clas_table$variable=="schalenw_abs"][1]
  w_schalenw_anz <- clas_table$weight[clas_table$variable=="schalenw_anz"][1]
  w_ruhe <- clas_table$weight[clas_table$variable=="ruhe"][1]
  w_fuetterung <- clas_table$weight[clas_table$variable=="fuetterung"][1]
  w_woerp <- clas_table$weight[clas_table$variable=="woerp"][1]
  w_riss_mod <- clas_table$weight[clas_table$variable=="riss_mod"][1]
  w_herden_s <-clas_table$weight[clas_table$variable=="herden_s"][1]
  w_nutztier_d <-clas_table$weight[clas_table$variable=="nutztier_d"][1]
  w_nutztier_r <-clas_table$weight[clas_table$variable=="nutztier_r"][1]
  w_gatterw <-clas_table$weight[clas_table$variable=="gatterw"][1]
  w_wald <-clas_table$weight[clas_table$variable=="wald"][1]
}

# calculate the sums
# Basic scenario for all seasons without scenarios
sum_all <- w_bevoelk*max(clas_table$value[clas_table$variable=="bevoelk"]) +
  w_hunde_a*max(clas_table$value[clas_table$variable=="hunde_a"]) +
  w_nacht_wi*max(clas_table$value[clas_table$variable=="nacht_wi"])/2 + # to not overaccount for tourism
  w_nacht_so*max(clas_table$value[clas_table$variable=="nacht_so"])/2 + # to not overaccount for tourism
  w_schalenw_abs*max(clas_table$value[clas_table$variable=="schalenw_abs"]) +
  w_schalenw_anz*max(clas_table$value[clas_table$variable=="schalenw_anz"]) +
  w_woerp*max(clas_table$value[clas_table$variable=="woerp"]) +
  w_riss_mod*max(clas_table$value[clas_table$variable=="riss_mod"]) +
  w_nutztier_d*max(clas_table$value[clas_table$variable=="nutztier_d"]) +
  w_gatterw*max(clas_table$value[clas_table$variable=="gatterw"]) +
  w_wald*max(clas_table$value[clas_table$variable=="wald"]) 

conflict_all <- (pd_bevoelk*w_bevoelk + 
                   pd_hunde_a*w_hunde_a +
                   (pd_nacht_wi*w_nacht_wi)/2 + 
                   (pd_nacht_so*w_nacht_so)/2 +
                   pd_schalenw_abs*w_schalenw_abs + 
                   pd_schalenw_anz*w_schalenw_anz + 
                   pd_woerp*w_woerp + 
                   pd_riss_mod*w_riss_mod +
                   pd_nutztier_d*w_nutztier_d + 
                   pd_gatterw*w_gatterw +
                   pd_wald*w_wald) /
                    sum_all

sum_summer <- w_bevoelk*max(clas_table$value[clas_table$variable=="bevoelk"]) +
  w_hunde_a*max(clas_table$value[clas_table$variable=="hunde_a"]) +
  w_nacht_so*max(clas_table$value[clas_table$variable=="nacht_so"]) +
  w_schalenw_abs*max(clas_table$value[clas_table$variable=="schalenw_abs"]) +
  w_schalenw_anz*max(clas_table$value[clas_table$variable=="schalenw_anz"]) +
  w_woerp*max(clas_table$value[clas_table$variable=="woerp"]) +
  w_riss_mod*max(clas_table$value[clas_table$variable=="riss_mod"]) +
  w_nutztier_d*max(clas_table$value[clas_table$variable=="nutztier_d"]) +
  w_gatterw*max(clas_table$value[clas_table$variable=="gatterw"]) +
  w_wald*max(clas_table$value[clas_table$variable=="wald"]) 

conflict_summer <- (pd_bevoelk*w_bevoelk + 
                      pd_hunde_a*w_hunde_a +
                      pd_nacht_so*w_nacht_so +
                      pd_schalenw_abs*w_schalenw_abs + 
                      pd_schalenw_anz*w_schalenw_anz + 
                      pd_woerp*w_woerp + 
                      pd_riss_mod*w_riss_mod +
                      pd_nutztier_d*w_nutztier_d + 
                      pd_gatterw*w_gatterw +
                      pd_wald*w_wald) /
                        sum_summer

sum_winter <- w_bevoelk*max(clas_table$value[clas_table$variable=="bevoelk"]) +
  w_hunde_a*max(clas_table$value[clas_table$variable=="hunde_a"]) +
  w_nacht_wi*max(clas_table$value[clas_table$variable=="nacht_wi"]) +
  w_schalenw_abs*max(clas_table$value[clas_table$variable=="schalenw_abs"]) +
  w_schalenw_anz*max(clas_table$value[clas_table$variable=="schalenw_anz"]) +
  w_woerp*max(clas_table$value[clas_table$variable=="woerp"]) +
  w_gatterw*max(clas_table$value[clas_table$variable=="gatterw"]) +
  w_wald*max(clas_table$value[clas_table$variable=="wald"]) 

conflict_winter <- (pd_bevoelk*w_bevoelk + 
                      pd_hunde_a*w_hunde_a +
                      pd_nacht_wi*w_nacht_wi + 
                      pd_schalenw_abs*w_schalenw_abs + 
                      pd_schalenw_anz*w_schalenw_anz + 
                      pd_woerp*w_woerp + 
                      pd_gatterw*w_gatterw +
                      pd_wald*w_wald) /
                        sum_winter

# 4) Transform the models into classes
# **************************
library(classInt)
library(patchwork)
library(cowplot) 
require(ggspatial)

# classify the models into 5 classes 
breaks_le <- classIntervals(values(le_wolf), n=5, style="jenks")$brks
breaks_co <- classIntervals(values(conflict_all), n=5, style="jenks")$brks
breaks_riss <- classIntervals(values(riss_mod), n=5, style="jenks")$brks

le_wolf_class <- classify(le_wolf, breaks_le, include.lowest=TRUE) 
conflict_all_class <- classify(conflict_all, breaks_co, include.lowest=TRUE) 
riss_mod_class <- classify(riss_mod, breaks_riss, include.lowest=TRUE) 


# 5) Plot the maps 
# **************************
colors_ko <- c("#f7f7f7", "#e2d9f3", "#cdbdf1", "#b7a1ef", "#a285ed")
colors_le <- c("#f7f7f7", "#d9f2d9", "#b2e6b2", "#8cda8c", "#66ce66")
colors_riss <- c("#f7f7f7",'#fecc5c','#fd8d3c','#f03b20','#bd0026') # be aware those are only used for the stand alone riss model
legend_labels_eng <- c("very low", "low", "medium", "high", "very high")

# Convert raster to data frame for ggplot
conflict_df <- as.data.frame(conflict_all_class, xy=TRUE)
colnames(conflict_df)[3] <- "value"

#make the eng plot
ko_plot_wo_title_eng <- ggplot(data = conflict_df) + 
  geom_raster(aes(x = x, y = y, fill = factor(value))) + 
  scale_fill_manual(values = colors_ko, labels = legend_labels_eng, name = "Conflict potential") + 
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")  +
  theme(legend.position = c(x = 0.15, y = 0.8),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_sf(crs= "epsg:3416") + 
  annotation_scale(location = "bl", width_hint = 0.5,
                   pad_x = unit(0.75, "cm"),
                   pad_y = unit(0.25, "cm"),
                   text_cex = 1)  +
  annotation_north_arrow(location = "br", which_north = "true", 
                         height = unit(1.5, "cm"), width = unit(1.5, "cm"),
                         style = north_arrow_fancy_orienteering,
                         pad_x = unit(1, "cm"),
                         pad_y = unit(0.25, "cm"),) 

le_plot_wo_title_eng <- ggplot(data = le_df) + 
  geom_raster(aes(x = x, y = y, fill = factor(value))) + 
  scale_fill_manual(values = colors_le, labels = legend_labels_eng, name = "Habitat potential") + 
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")  +
  theme(legend.position = c(x = 0.15, y = 0.8),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_sf(crs= "epsg:3416") + 
  annotation_scale(location = "bl", width_hint = 0.5,  
                   pad_x = unit(0.75, "cm"),
                   pad_y = unit(0.25, "cm"),
                   text_cex = 1)  +
  annotation_north_arrow(location = "br", which_north = "true", 
                         height = unit(1.5, "cm"), width = unit(1.5, "cm"),
                         style = north_arrow_fancy_orienteering,
                         pad_x = unit(1, "cm"),
                         pad_y = unit(0.25, "cm"),) 


riss_plot_wo_title_eng <- ggplot(data = riss_mod_df) + 
  geom_raster(aes(x = x, y = y, fill = factor(value))) + 
  scale_fill_manual(values = colors_riss, labels = legend_labels_eng, name = "Depredation susceptibility") + 
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")  +
  theme(legend.position = c(x = 0.25, y = 0.8),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  coord_sf(crs= "epsg:3416") + 
  annotation_scale(location = "bl", width_hint = 0.5,  
                   pad_x = unit(0.75, "cm"),
                   pad_y = unit(0.25, "cm"),
                   text_cex = 1)  +
  annotation_north_arrow(location = "br", which_north = "true", 
                         height = unit(1.5, "cm"), width = unit(1.5, "cm"),
                         style = north_arrow_fancy_orienteering,
                         pad_x = unit(1, "cm"),
                         pad_y = unit(0.25, "cm"),) 


# 6) Hybrid map with Le and Ko with pretty and understandable colors and classes
# **************************
raster_df <- as.data.frame(combined_classes, xy=TRUE)
colnames(raster_df)[3] <- "value"

# Create the bivariate color scale
bivariate_colors <- c(
  "0"  = "#f7f7f7",  "1"  = "#e2d9f3",  "2"  = "#cdbdf1",  "3"  = "#b7a1ef",  "4"  = "#a285ed",
  "10" = "#d9f2d9", "11" = "#c4d6e6", "12" = "#afbbe2", "13" = "#9a9fde", "14" = "#8583da",
  "20" = "#b2e6b2", "21" = "#9ccadb", "22" = "#86add3", "23" = "#7091cb", "24" = "#5a74c3",
  "30" = "#8cda8c", "31" = "#76becf", "32" = "#60a2c3", "33" = "#4a86b7", "34" = "#3469ab",
  "40" = "#66ce66", "41" = "#50b2c3", "42" = "#3a96b3", "43" = "#247aa3", "44" = "#0e5e93"
)

# create the legend as an extra grid
legend_data <- expand.grid(
  Le_Class = 0:4,  # X-Axis = habitat potential
  Ko_Class = 0:4   # Y-Axis = conflict potential
)

legend_data$class <- factor(legend_data$Le_Class * 10 + legend_data$Ko_Class, 
                            levels = names(bivariate_colors)) 

# Create and legend plot
legend_plot_eng <- ggplot(legend_data, aes(x = Le_Class, y = Ko_Class, fill = class)) +
  geom_tile() +
  scale_fill_manual(values = bivariate_colors, na.value = "white") +  # Falls NA, auf Weiß setzen
  labs(x = "Habitat potential", y = "Conflict potential") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),  # Explicitly set x-axis title
        axis.title.y = element_text(size = 12),  # Explicitly set y-axis title
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", # set those two to make the box around the plot disappear
                                        colour = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")
plot(legend_plot_eng)

# Create the plot wihtout the legend
hybrid_plot_eng <- ggplot(data = raster_df) + # define raster
  geom_raster(aes(x = x, y = y, fill = factor(value))) + # plot values
  scale_fill_manual(values = bivariate_colors) + # define colours
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  guides(fill="none") +
  coord_sf(crs= "epsg:3416") + # Define crs
  annotation_scale(location = "bl", width_hint = 0.5,  
                   pad_x = unit(0.75, "cm"),
                   pad_y = unit(0.25, "cm"),
                   text_cex = 1)  +
  annotation_north_arrow(location = "br", which_north = "true", 
                         height = unit(1.5, "cm"), width = unit(1.5, "cm"),
                         style = north_arrow_fancy_orienteering,
                         pad_x = unit(1, "cm"),
                         pad_y = unit(0.25, "cm"),) 

plot(hybrid_plot_eng)

combined_plot_eng <- ggdraw() +
  draw_plot(hybrid_plot_eng)  +
  draw_plot(legend_plot_eng,  x = 0.12, y = 0.6, width = 0.3, height = 0.3) # Define legend position
plot(combined_plot_eng)

## END
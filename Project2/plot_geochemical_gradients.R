library(tidyverse)

raw_data <- read_csv("Saanich_TimeSeries_Chemical_DATA.csv")

data <- select(raw_data, Cruise, Depth, CTD_O2, NO3, Mean_H2S, Mean_CH4) %>%
        rename(O2_uM = CTD_O2, NO3_uM = NO3, H2S_uM = Mean_H2S, CH4_uM = Mean_CH4)

data %>%
  filter(Cruise == 72) %>% 
  ggplot(aes(x = NO3_uM, y = Depth)) +
  geom_point(colour = "red") +
  scale_y_reverse(limits = c(200, 0)) +
  xlab("Total Concentration of NO3") + ylab("Depth In Meters")

data %>%
  filter(Cruise == 72) %>% 
  ggplot(aes(x = O2_uM, y = Depth)) +
  geom_point(colour = "purple") +
  scale_y_reverse(limits = c(200, 0)) +
  xlab("Total Concentration of O2") + ylab("Depth In Meters")

data %>%
  filter(Cruise == 72) %>% 
  ggplot(aes(x = CH4_uM, y = Depth)) +
  geom_point(colour = "darkgreen") +
  scale_y_reverse(limits = c(200, 0)) +
  xlab("Total Concentration of CH4") + ylab("Depth In Meters")

data %>%
  filter(Cruise == 72) %>% 
  ggplot(aes(x = H2S_uM, y = Depth)) +
  geom_point(colour = "blue") +
  scale_y_reverse(limits = c(200, 0)) +
  xlab("Total Concentration of H2S") + ylab("Depth In Meters")

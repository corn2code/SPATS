setwd("/Users/vladimir/Downloads")

library(tidyverse)
library(SpATS)
library(data.table)

# Read in data 
NE2020 <- fread("tpj16801-sup-0002-datasets1.csv", data.table = F)

NE2020$Anthesis <- as.numeric(as.character(NE2020$Anthesis))
NE2020$Silking <- as.numeric(as.character(NE2020$Silking))

head(NE2020)

# Define function
getSpatialCorrections <- function(data, response)
{
  #data <- NE2020
  #response <- "Anthesis"
  # Declare empty df
  df.sp <- tibble(plot = NULL, '{response}':= NULL)
  loc.df <- filter(data, !is.na(Row) & !is.na(Column) & !is.na(.data[[response]])) %>%
    mutate(as.factor(plot))
  
  ColumnKnots <- floor(max(loc.df$Column, na.rm = TRUE)/2) + 1
  RowKnots <- floor(max(loc.df$Row, na.rm = TRUE)/2) + 1
  model <- SpATS(response, genotype = 'plot', genotype.as.random = TRUE,
                 spatial = ~ SAP(Column, Row, nseg = c(ColumnKnots, RowKnots)),
                 data = loc.df)
  
  # Extract BLUPS
  intercept <- model$coeff['Intercept']
  sp <- as_tibble(model$coeff, rownames = 'plot') %>%
    filter(!is.na(plot)) %>%
    rowwise() %>%
    mutate(value = value + intercept) %>%
    rename('{response}':= value)
  # Bind to df
  df.sp <- bind_rows(df.sp, sp) %>%
    filter(!str_detect(plot, "[a-zA-Z]")) %>%  # Keep only rows without letters in `plot`
    mutate(plot = as.numeric(plot))  # Convert `plot` to numeric
  # Return df
  return(df.sp)
}

colnames(NE2020)


for(j in c('Anthesis', 'Silking')) {
  NE2020 <- full_join(NE2020, getSpatialCorrections(NE2020, j), join_by(plot), suffix = c('', '.sp'), keep = FALSE)
}

fwrite(NE2020, 'tpj16801-sup-0002-datasets1_coeff.csv')

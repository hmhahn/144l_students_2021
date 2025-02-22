EEMB 144L: Bacterial Abundance 144L
================
Hope Hahn
10/18/2021

``` r
#load packages
library(tidyverse)
library(readxl)
library(lubridate)
```

## Import Data

I loaded data from the 2021 Bacterial Abundance excel sheet, which had
two sheets: metadata and data. Each sheet was loaded into separate
dataframes called data and metadata.

``` r
# Use function excel_sheets and the argument is the pathway to get to the excel sheet
#excel_sheets("~/Desktop/github/144l_students_2021/Input_Data/week4/144L_2021_BactAbund.xlsx")

# Create dataset "metadata" by loading the specific sheet called metadata from the excel sheet
metadata <- read_excel("~/Desktop/github/144l_students_2021/Input_Data/week4/144L_2021_BactAbund.xlsx", sheet = "Metadata")

# Create dataset "data" by loading the specific sheet called data from the excel sheet
data <- read_excel("~/Desktop/github/144l_students_2021/Input_Data/week4/144L_2021_BactAbund.xlsx", sheet = "Data")

# join datasets
joined <- left_join(metadata, data)
```

# Prepare Data for Plotting

In this section, a new dataset called cells was made. I split the date
and time from the Datetime column, and converted the values to numeric
rather than character. Then, the all_cells_uL column values were
converted to numeric values in liters rather than microliters. The
interval between each time point in hours and days were calculated as
well.

``` r
cells <- joined %>%
  mutate(Datetime = ymd_hm(Datetime), #splits datetime
  cells_L = as.numeric(all_cells_uL)*1000000) %>% #convert cells uL to cells L in new column
  group_by(Treatment, Bottle) %>% #group dataset to calculate time elapsed by treatment and by bottle
  mutate(interv = interval(first(Datetime), Datetime), #get interval between each time point
         s = as.numeric(interv), #turn into numeric
         hours = s/3600, # gethours
         days = hours/24) %>% #get days
  ungroup() %>%
  select(Experiment: DNA_Sample, cells_L, hours, days) %>% #clean up
  drop_na(cells_L) #drop na values
```

## Plot Data

To help the plotting process, colors and levels were pre-assigned to the
treatments. I then used ggPlot and my custom colors/levels to plot the
data.

``` r
# assign different colors to different treatments
custom.colors <- c("Control" = "lightpink", "Kelp Exudate" = "cornflowerblue", "Kelp Exudate_Nitrate_Phosphate" = "darkseagreen", "Glucose_Nitrate_Phosphate" = "mediumpurple1")

# assign levels to control order in legend
levels <- c("Control", "Kelp Exudate", "Kelp Exudate_Nitrate_Phosphate", "Glucose_Nitrate_Phosphate")
```

``` r
# Plot data
cells %>% 
  mutate(dna = ifelse(DNA_Sample == T, "*", NA)) %>% #add new column with * and NA instead 
  ggplot(aes(x=days, y=cells_L, group = interaction(Treatment, Bottle))) + 
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1) + #add lines
  geom_point(aes(fill = factor(Treatment, levels = levels), stroke = NA), size = 3, shape = 21) + #add points
  geom_text(aes(label = dna), size = 12, color = "goldenrod1") + #add where dna time points are
  labs(x = "Days", y= expression(paste("Cells, L"^-1)), fill = "") + #new graph labels
  guides(color = "none") + 
  scale_color_manual(values = custom.colors) + #customize colors
  scale_fill_manual(values = custom.colors) + #customize colors
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
panel.background = element_blank(), axis.line = element_line(colour = "black")) #make graph nice looking
```

![](2021_bacterial_abundance_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
  #facet_grid(rows = "Treatment", scales = "free") 
```

## More Plotting Data

To visualize exponential growth, different plots needed to be made. A
new dataset called ln_cells was made, which had the same data, but had
the natural log of cells_L as well as the change in the natural log of
cells between each point.

``` r
# Create column with natural log
ln_cells <- cells %>%
  group_by(Treatment, Bottle) %>%
  mutate(ln_cells = log(cells_L),
         diff_ln_cells = ln_cells - lag(ln_cells, default = first(ln_cells)))
```

There are two plots: the change in ln(cells per liter) vs time, and just
ln(cells per liter). Both can be used to visualize when exponential
growth is occurring.

``` r
# plot change in natural log
ln_cells %>% 
   mutate(dna = ifelse(DNA_Sample == T, "*", NA)) %>% 
  ggplot(aes(x=days, y=diff_ln_cells, group = interaction(Treatment, Bottle))) + 
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1) + 
  geom_point(aes(fill = factor(Treatment, levels = levels), stroke = NA), size = 3, color = "black", shape = 21) + 
  geom_text(aes(label = dna), size = 12, color = "goldenrod1") + 
  labs(x = "Days", y= expression(paste("∆ln cells, L"^-1)), fill = "") + 
  guides(color = "none") + 
  scale_color_manual(values = custom.colors) + 
  scale_fill_manual(values = custom.colors) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap("Bottle", ncol = 2) 
```

![](2021_bacterial_abundance_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# plot natural log
ln_cells %>% 
   mutate(dna = ifelse(DNA_Sample == T, "*", NA)) %>% 
  ggplot(aes(x=days, y=ln_cells, group = interaction(Treatment, Bottle)))+ 
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1) + 
  geom_point(aes(fill = factor(Treatment, levels = levels), stroke = NA), size = 3, color = "black", shape = 21) + 
  geom_text(aes(label = dna), size = 12, color = "goldenrod1") + 
  labs(x = "Days", y= expression(paste("ln cells, L"^-1)), fill = "") + 
  guides(color = "none") + 
  scale_color_manual(values = custom.colors) + 
  scale_fill_manual(values = custom.colors) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap("Bottle", ncol = 2) 
```

![](2021_bacterial_abundance_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

By looking at these graphs, exponential growth seems to be occurring
early on in the experiment. For bottles E and G, there is no growth
exhibited in the beginning of the graph, but by looking at the replicate
bottles, F and H respectively, we can assume that E and G also exhibit
growth during days 0-1.

## Summary

This assignment introduced us to the basics of R as well as helped us to
visualize the data we collected for our experiment. There were some
complications with labeling the samples, and some of the labels came
off, but we are still able to make out overall patterns per treatment.
All treatments appear to have exponential growth in the beginning of the
experiment, which was expected considering that is when the additions to
the seawater were the most abundant.

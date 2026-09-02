#plot for streamflow in the mill river for wet center. 
library(readr)
library(ggplot2)
library(dplyr)
Amherst__Discharge <- read_csv("WetCenterFlowFigureData/Wet Center Flow April 2025 to Aug 2026.csv")

Amherst__Discharge <- as_tibble(Amherst__Discharge)

# Convert Time to actual datetime
library(lubridate)

Amherst__Discharge <- Amherst__Discharge %>%
  mutate(Time = as.POSIXct(Time, format = "%m/%d/%Y %H:%M"))

PFAS_inventory_poster <- read_csv("PFAS_inventory_poster.csv")

#create datetime column

PFAS_inventory_poster_DateTime <- PFAS_inventory_poster %>%
  filter(!is.na(Date), !is.na(Time)) %>%
  mutate(
    DateTime = parse_date_time(
      paste(Date, Time),
      orders = c("mdy IMS", "mdy HM", "mdy IM", "mdy HMS"),
      tz = "America/New_York"
    )
  ) %>%
  arrange(DateTime)


#2) filter inventory to just Wet Center and Stream type

# Extract Stream samples
stream_samples <- PFAS_inventory_poster_DateTime %>%
  filter(Location == "Wet Center", Type == "Stream") %>%
  arrange(DateTime)

# Sample dates for wet center pfas as a vector 
sample_dates <- stream_samples$DateTime

# Pull the corresponding flow values from your data
sample_points <- Amherst__Discharge %>%
  filter(as.Date(Time, format = "%m/%d/%Y %H:%M") %in% as.Date(sample_dates))

ggplot(Amherst__Discharge, aes(x = Time, y = Flow)) +
  geom_line(color = "steelblue", linewidth = 1) +
  #geom_point(data = sample_points, aes(x = Time, y = Flow),
             #color = "red", size = 3) +
  scale_y_log10() +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title = "Mill River Discharge in Amherst",
    x = "Date",
    y = "Discharge (m3/s)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold")
  )
#the warning is fine (removed 183 rows)

#prep info for second y axis
flow_range <- range(Amherst__Discharge$Flow, na.rm = TRUE)

library(readr)
raw <- read_csv("PFAS Head Group Data.csv", col_names = TRUE)
library(dplyr)
library(stringr)

headgroup_lookup <- tibble(
  column    = names(raw),
  Headgroup = str_remove(names(raw), "\\.\\.\\.\\d+$") %>% str_to_upper(),
  Compound  = as.character(unlist(raw[1, ]))
) %>%
  filter(!is.na(Compound)) %>%
  select(Compound, Headgroup)

headgroup_lookup

library(tidyr)

# ============================================================
# SETUP
# ============================================================
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)

# ============================================================
# LOAD DISCHARGE DATA
# ============================================================
Amherst__Discharge <- read_csv("WetCenterFlowFigureData/Wet Center Flow April 2025 to Aug 2026.csv") %>%
  as_tibble() %>%
  mutate(Time = as.POSIXct(Time, format = "%m/%d/%Y %H:%M"))

flow_range <- range(Amherst__Discharge$Flow, na.rm = TRUE)

# ============================================================
# LOAD PFAS INVENTORY DATA + CREATE DATETIME
# ============================================================
PFAS_inventory_poster <- read_csv("PFAS_inventory_poster.csv")
# drop rows with missing date/time
PFAS_inventory_poster_DateTime <- PFAS_inventory_poster %>%
  filter(!is.na(Date), !is.na(Time), !is.na(`PFOA Results`)) %>%          
  mutate(
    DateTime = parse_date_time(
      paste(Date, Time),
      orders = c("mdy IMS", "mdy HM", "mdy IM", "mdy HMS"),
      tz = "America/New_York"
    )
  ) %>%
  arrange(DateTime)

# ============================================================
# BUILD HEADGROUP LOOKUP TABLE
# ============================================================
raw <- read_csv("PFAS Head Group Data.csv", col_names = TRUE)

headgroup_lookup <- tibble(
  column    = names(raw),
  Headgroup = str_remove(names(raw), "\\.\\.\\.\\d+$") %>% str_to_upper(),
  Compound  = as.character(unlist(raw[1, ]))
) %>%
  filter(!is.na(Compound)) %>%
  select(Compound, Headgroup)

# add compounds not covered by the reference file (incl. FTCA/FTS precursors)
extra_lookup <- tibble(
  Compound = c("PFUnA", "PFMBA", "PFMPA", "HFPO-DA", "9Cl-PF3ONS",
               "NEtFOSA", "NEtFOSE", "NMeFOSE",
               "3-3 FTCA", "5-3 FTCA", "7-3 FTCA", "4:2FTS", "6:2FTS", "8:2FTS"),
  Headgroup = c("PFCA", "PFCA", "PFCA", "PFECA", "PFESA",
                "SULFONAMIDE", "SULFONAMIDE", "SULFONAMIDE",
                "FTCA", "FTCA", "FTCA", "FTS", "FTS", "FTS")
)

headgroup_lookup <- bind_rows(headgroup_lookup, extra_lookup)

# ============================================================
# BUILD PFAS_BY_HEADGROUP (long format, summed by headgroup per sample)
# ============================================================
pfas_by_headgroup <- PFAS_inventory_poster_DateTime %>%
  filter(Location == "Wet Center") %>%
  select(`Sample #`, DateTime, ends_with("Results")) %>%
  pivot_longer(
    cols = ends_with("Results"),
    names_to = "Compound",
    values_to = "Concentration"
  ) %>%
  mutate(
    Compound = str_remove(Compound, " Results$"),
    Concentration = as.numeric(str_extract(Concentration, "^[0-9.]+"))
  ) %>%
  left_join(headgroup_lookup, by = "Compound") %>%
  group_by(`Sample #`, DateTime, Headgroup) %>%
  summarize(TotalConcentration = sum(Concentration, na.rm = TRUE), .groups = "drop")

# exclude zero values so log10 doesn't choke
pfas_range <- range(pfas_by_headgroup$TotalConcentration[pfas_by_headgroup$TotalConcentration > 0], na.rm = TRUE)
log_pfas_range <- log10(pfas_range)

scale_factor <- diff(flow_range) / diff(log_pfas_range)

pfas_by_headgroup <- pfas_by_headgroup %>%
  mutate(ScaledConcentration = ifelse(
    TotalConcentration > 0,
    (log10(TotalConcentration) - log_pfas_range[1]) * scale_factor + flow_range[1],
    NA
  ))

library(patchwork)

# make sure discharge NAs are filtered first
Amherst__Discharge <- Amherst__Discharge %>%
  filter(!is.na(Time))

# panel 1: 
p1 <- ggplot(pfas_by_headgroup, aes(x = DateTime, y = TotalConcentration, color = Headgroup)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
  labs(y = "PFAS Concentration (ng/L)", x = "Date") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 14)
  )

# panel 2: 
p2 <- ggplot(Amherst__Discharge, aes(x = Time, y = Flow)) +
  geom_line(color = "steelblue", linewidth = 1) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(title = "Mill River Discharge and PFAS Concentration in Amherst",
       y = "Discharge (m3/s)", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold")
  )

p1 / p2 + plot_layout(heights = c(2, 1.5))



#plot for mill river disharge and PFAS, needs to be assessed
ggplot() +
  geom_line(data = Amherst__Discharge, aes(x = Time, y = Flow),
            color = "steelblue", linewidth = 1) +
  geom_point(data = pfas_by_headgroup, aes(x = DateTime, y = ScaledConcentration, color = Headgroup),
             size = 3) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b %Y") +
  scale_y_continuous(
    name = "Discharge (m3/s)",
    sec.axis = sec_axis(
      transform = ~ (. - flow_range[1]) / scale_factor + pfas_range[1],
      name = "PFAS Concentration (ng/L)"
    )
  ) +
  labs(
    title = "Mill River Discharge and PFAS Concentration in Amherst",
    x = "Date"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold")
  )

#plot pfas and flow for 3 month increments with temporal flow

library(patchwork)
library(lubridate)

make_discharge_pfas_plot <- function(start_date, end_date, title_suffix) {
  
  discharge_sub <- Amherst__Discharge %>%
    filter(Time >= as.POSIXct(start_date, tz = "America/New_York"),
           Time <= as.POSIXct(end_date, tz = "America/New_York"))
  
  pfas_sub <- pfas_by_headgroup %>%
    filter(DateTime >= as.POSIXct(start_date, tz = "America/New_York"),
           DateTime <= as.POSIXct(end_date, tz = "America/New_York"))
  
  p1 <- ggplot(discharge_sub, aes(x = Time, y = Flow)) +
    geom_line(color = "steelblue", linewidth = 1) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b %Y") +
    labs(title = paste("Mill River Discharge and PFAS Concentration:", title_suffix),
         y = "Discharge (m3/s)", x = NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
  
  p2 <- ggplot(pfas_sub, aes(x = DateTime, y = TotalConcentration, color = Headgroup)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b %Y") +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
    labs(y = "PFAS Concentration (ng/L)", x = "Date") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 14)
    )
  
  p1 / p2 + plot_layout(heights = c(2, 1.5))
}

# Four windows covering Jul 2025 - Aug 2026
plot_q1 <- make_discharge_pfas_plot("2025-07-01", "2025-09-30", "Jul-Sep 2025")
plot_q2 <- make_discharge_pfas_plot("2025-10-01", "2025-12-31", "Oct-Dec 2025")
plot_q3 <- make_discharge_pfas_plot("2026-01-01", "2026-04-30", "Jan-Apr 2026")
plot_q4 <- make_discharge_pfas_plot("2026-05-01", "2026-08-31", "May-Aug 2026")

# view them one at a time
plot_q1
plot_q2
plot_q3
plot_q4

#same as above, but with the same axes values

library(patchwork)
library(lubridate)
library(ggbreak)

# compute fixed ranges from the FULL dataset once, outside the function
flow_range_full <- range(Amherst__Discharge$Flow, na.rm = TRUE)
pfas_range_full <- range(pfas_by_headgroup$TotalConcentration, na.rm = TRUE)

make_discharge_pfas_plot <- function(start_date, end_date, title_suffix) {
  
  discharge_sub <- Amherst__Discharge %>%
    filter(Time >= as.POSIXct(start_date, tz = "America/New_York"),
           Time <= as.POSIXct(end_date, tz = "America/New_York"))
  
  pfas_sub <- pfas_by_headgroup %>%
    filter(DateTime >= as.POSIXct(start_date, tz = "America/New_York"),
           DateTime <= as.POSIXct(end_date, tz = "America/New_York"))
  
  p1 <- ggplot(discharge_sub, aes(x = Time, y = Flow)) +
    geom_line(color = "steelblue", linewidth = 1) +
    scale_x_break(
      c(as.Date("2025-04-13"), as.Date("2025-07-09"))
    ) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b %Y") +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),limits = flow_range_full) +
    labs(title = paste("Mill River Discharge and PFAS Concentration:", title_suffix),
         y = "Discharge (m3/s)", x = NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
  
  p2 <- ggplot(pfas_sub, aes(x = DateTime, y = TotalConcentration, color = Headgroup)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%b %Y") +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), limits = pfas_range_full) +
    labs(y = "PFAS Concentration (ng/L)", x = "Date") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 14)
    )
  
  p1 / p2 + plot_layout(heights = c(2, 1.5))
}

plot_q1 <- make_discharge_pfas_plot("2025-07-09", "2025-09-30", "Jul-Sep 2025")
plot_q2 <- make_discharge_pfas_plot("2025-10-01", "2025-12-31", "Oct-Dec 2025")
plot_q3 <- make_discharge_pfas_plot("2026-01-01", "2026-04-30", "Jan-Apr 2026")
plot_q4 <- make_discharge_pfas_plot("2026-05-01", "2026-08-31", "May-Aug 2026")

plot_q1
plot_q2
plot_q3
plot_q4

#make zoomed in plots on the squares




#create isotope data using Date Time, filter inventory to just Wet Center and Stream type

isotope_data <- PFAS_inventory_poster_DateTime %>%
  select(`Location`, `Type`, `Delta 2H`, `Delta 18O`, `Delta 17O`, DateTime)%>%
  arrange(DateTime)


#2) filter inventory to just Wet Center and Stream type
isotope_stream_samples <- isotope_data %>%
  filter(Location == "Wet Center", Type == "Stream") %>%
  arrange(DateTime)

# Sample dates for wet center pfas as a vector 
#usied befoe: sample_dates <- stream_samples$DateTime
isotope_sample_dates <- isotope_stream_samples$DateTime


# Pull the corresponding flow values from your data
#sample_points <- Amherst__Discharge %>%
  #filter(as.Date(Time, format = "%m/%d/%Y %H:%M") %in% as.Date(sample_dates))
isotope_sample_points <- Amherst__Discharge %>%
  filter(as.Date(Time, format = "%m/%d/%Y %H:%M") %in% as.Date(sample_dates))

#test for plots

#dont like this one, I want to break up the timeline and have 3 axes
ggplot() +
  # flow as the base line, on its natural axis
  geom_line(data = Amherst__Discharge,
            aes(x = as.POSIXct(Time, format = "%m/%d/%Y %H:%M"), y = Flow),
            color = "grey40") +
  # PFAS headgroup totals (scaled)
  geom_point(data = pfas_by_headgroup,
             aes(x = DateTime, y = ScaledConcentration, color = Headgroup),
             size = 2) +
  # isotopes (scaled), shaped differently so they're visually distinct from PFAS
  geom_point(data = iso_scaled,
             aes(x = DateTime, y = ScaledValue, shape = Isotope),
             color = "black", size = 2) +
  scale_y_continuous(
    name = "Flow",
    sec.axis = sec_axis(~ ., name = "Scaled PFAS / Isotope values")
  ) +
  labs(x = "Date", title = "Flow, PFAS Headgroups, and Isotopes Over Time") +
  theme_minimal()

#3 panels
library(ggplot2)
library(patchwork)

shared_xlim <- range(c(Amherst__Discharge$Time, pfas_by_headgroup$DateTime, iso_long$DateTime))

# ---- Panel 1: Flow (linear) ----
p_flow <- ggplot(Amherst__Discharge,
                  aes(x = as.POSIXct(Time, format = "%m/%d/%Y %H:%M"), y = Flow)) +
  geom_line(color = "grey40") +
  labs(y = "Flow") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# ---- Panel 2: PFAS (log scale) ----
p_pfas <- ggplot(pfas_by_headgroup,
                  aes(x = DateTime, y = TotalConcentration, color = Headgroup)) +
  geom_point(size = 2) +
  scale_y_log10(labels = label_number()) +
  labs(y = "PFAS Conc. (log scale)") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# ---- Panel 3: Isotopes (linear) ----
p_iso <- ggplot(iso_long, aes(x = DateTime, y = Value, shape = Isotope)) +
  geom_point(size = 2, color = "black") +
  labs(y = "Isotope (‰)", x = "Date") +
  theme_minimal()

# ---- Stack them, aligned x-axis ----
(p_flow / p_pfas / p_iso) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(title = "Flow, PFAS (log), and Isotopes Over Time")
            
#test

  make_segment_plot <- function(start_date, end_date) {
    
    flow_seg <- Amherst__Discharge %>%
      mutate(Time = as.POSIXct(Time, format = "%m/%d/%Y %H:%M")) %>%
      filter(Time >= start_date, Time <= end_date)
    
    pfas_seg <- pfas_by_headgroup %>%
      filter(DateTime >= start_date, DateTime <= end_date)
    
    iso_seg <- iso_long %>%
      filter(DateTime >= start_date, DateTime <= end_date)
    
    p_flow <- ggplot(flow_seg, aes(x = Time, y = Flow)) +
      geom_line(color = "grey40") +
      labs(y = "Flow") +
      xlim(start_date, end_date) +
      theme_minimal() +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    p_pfas <- ggplot(pfas_seg, aes(x = DateTime, y = TotalConcentration, color = Headgroup)) +
      geom_point(size = 2) +
      scale_y_log10(labels = label_number()) +
      xlim(start_date, end_date) +
      labs(y = "PFAS Conc. (log)") +
      theme_minimal() +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    p_iso <- ggplot(iso_seg, aes(x = DateTime, y = Value, shape = Isotope)) +
      geom_point(size = 2, color = "black") +
      xlim(start_date, end_date) +
      labs(y = "Isotope (\u2030)", x = "Date") +
      theme_minimal()
    
    (p_flow / p_pfas / p_iso) +
      plot_annotation(
        title = paste0(format(start_date, "%b %Y"), " \u2013 ", format(end_date, "%b %Y"))
      )
  }

# ============================================================
# QUARTERLY PLOTS (Jun-Sep / Oct-Dec / Jan-Apr / May-Aug)
# ============================================================
plot_q1 <- make_segment_plot(as.POSIXct("2025-06-24"), as.POSIXct("2025-09-30"))  # Jun-Sep 2025
plot_q2 <- make_segment_plot(as.POSIXct("2025-10-01"), as.POSIXct("2025-12-31"))  # Oct-Dec 2025
plot_q3 <- make_segment_plot(as.POSIXct("2026-01-01"), as.POSIXct("2026-04-30"))  # Jan-Apr 2026
plot_q4 <- make_segment_plot(as.POSIXct("2026-05-01"), as.POSIXct("2026-08-31"))  # May-Aug 2026

# view individually
plot_q1
plot_q2
plot_q3
plot_q4

#isotope stuff
gmwl <- function(d18O) 8 * d18O + 10

ggplot(isotope_stream_samples, aes(x = `Delta 18O`, y = `Delta 2H`)) +
  geom_abline(intercept = 10, slope = 8, color = "blue", linetype = "dashed") +
  annotate("text", x = min(isotope_stream_samples$`Delta 18O`, na.rm=TRUE),
           y = gmwl(min(isotope_stream_samples$`Delta 18O`, na.rm=TRUE)) + 5,
           label = "GMWL", color = "blue", hjust = 0) +
  geom_point(aes(color = DateTime), size = 3) +
  scale_color_gradient(low = "orange", high = "darkblue") +
  labs(x = expression(delta^18*O~"(\u2030)"), y = expression(delta^2*H~"(\u2030)"),
       title = "Dual Isotope Plot: Wet Center Stream Samples") +
  theme_minimal()

isotope_stream_samples_excess <- isotope_stream_samples %>%
  mutate(d_excess = `Delta 2H` - 8 * `Delta 18O`)

ggplot(isotope_stream_samples_excess, aes(x = DateTime, y = d_excess)) +
  geom_line(color = "grey50") +
  geom_point(size = 2) +
  labs(y = "d-excess (\u2030)", x = "Date", title = "Deuterium Excess Over Time") +
  theme_minimal()
            

suppressPackageStartupMessages({
  library(sf); library(igraph); library(ggplot2); library(scales); library(dplyr)
})
source("R/build_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6L) stop(paste(
  "usage: Rscript scripts/build_bogota_segment_base.R",
  "ROADS_RDS POLYGON_RDATA CRASH_CSV START_DATE END_DATE OUTPUT_DIR"
))
start_date <- as.Date(args[[4]]); end_date <- as.Date(args[[5]])
if (is.na(start_date) || is.na(end_date) || start_date > end_date) stop("Invalid dates")
output_dir <- args[[6]]; dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

roads <- readRDS(args[[1]])
if (!inherits(roads, "sf")) stop("ROADS_RDS is not sf")
# EPSG:3116 easting cutoff that removes the rural eastern extensions toward
# La Calera and Choachi while retaining the central urban road network.
urban_easting_max <- 1005000
southeast_easting_min <- 1003500
southeast_northing_max <- 999000
urban_northing_min <- 998000
polygon <- read_study_polygon(
  args[[2]], urban_easting_max,
  southeast_easting_min, southeast_northing_max, urban_northing_min
)
crashes_all <- read_unique_crashes(args[[3]])
crashes <- crashes_all[!is.na(crashes_all$analysis_date) &
                         crashes_all$analysis_date >= start_date &
                         crashes_all$analysis_date <= end_date, ]
polygon_crs <- st_transform(polygon, st_crs(crashes))
crashes <- crashes[lengths(st_within(crashes, polygon_crs)) > 0, ]
if (!nrow(crashes)) stop("No unique crashes remain after filtering")

supports_all <- extract_network_supports(build_road_network(roads, polygon))
excluded_disconnected_segments <- sum(!supports_all$segments$analysis_keep)
connected_lines <- supports_all$segments[supports_all$segments$analysis_keep, ]
connected_lines <- connected_lines[, intersect(
  c("source_row", "osm_id", "name", "highway", "geometry"),
  names(connected_lines)
)]
# Reconstruct every downstream object from the single main connected component.
# Thus isolated rural fragments cannot remain in maps, adjacency or model data.
supports <- extract_network_supports(as_sfnetwork(connected_lines, directed = FALSE))
assigned <- assign_crashes_to_segments(crashes, supports$segments)
segment_count <- tabulate(assigned$segment_id, nbins = nrow(supports$segments))
supports$segments$count <- as.integer(segment_count)
supports$segments$exposure_length_m <- supports$segments$length_m
supports$segments$log_exposure <- log(supports$segments$length_m)

if (sum(segment_count) != nrow(assigned)) stop("Crashes were not assigned exactly once")
if (any(segment_count < 0) || any(segment_count != round(segment_count))) stop("Invalid counts")
if (any(!is.finite(supports$segments$log_exposure))) stop("Invalid segment exposure")
if (vcount(supports$segment_graph) != nrow(supports$segments)) stop("Support mismatch")

segment_csv <- st_drop_geometry(supports$segments)
intersection_xy <- st_coordinates(supports$intersections)
intersection_csv <- st_drop_geometry(supports$intersections)
intersection_csv$x <- intersection_xy[, 1]; intersection_csv$y <- intersection_xy[, 2]
write.csv(segment_csv, file.path(output_dir, "bogota_segments_model_base.csv"), row.names = FALSE)
write.csv(intersection_csv, file.path(output_dir, "bogota_intersections.csv"), row.names = FALSE)
write.csv(supports$segment_adjacency,
          file.path(output_dir, "bogota_segment_adjacency.csv"), row.names = FALSE)
far_xy <- st_coordinates(assigned)
far_assignments <- st_drop_geometry(assigned)
far_assignments$x <- far_xy[, 1]; far_assignments$y <- far_xy[, 2]
far_assignments <- far_assignments[far_assignments$assignment_distance_m > 25, ]
write.csv(far_assignments,
          file.path(output_dir, "bogota_assignments_over_25m.csv"), row.names = FALSE)

saveRDS(list(
  road_graph = supports$road_graph,
  segment_graph = supports$segment_graph,
  segments = supports$segments,
  intersections = supports$intersections,
  segment_adjacency = supports$segment_adjacency,
  assigned_crashes = assigned,
  start_date = start_date, end_date = end_date,
  urban_easting_max_3116 = urban_easting_max,
  southeast_easting_min_3116 = southeast_easting_min,
  southeast_northing_max_3116 = southeast_northing_max,
  urban_northing_min_3116 = urban_northing_min,
  response_support = "road segments (vertices of the line graph)",
  adjacency_definition = "two segments are adjacent when they share a road intersection",
  exposure_definition = "segment length in metres"
), file.path(output_dir, "bogota_segment_model_data.rds"))

gpkg <- file.path(output_dir, "bogota_segment_spatial_layers.gpkg")
st_write(supports$segments, gpkg, "segments", delete_dsn = TRUE, quiet = TRUE)
st_write(supports$intersections, gpkg, "intersections", append = TRUE, quiet = TRUE)
st_write(assigned, gpkg, "assigned_crashes", append = TRUE, quiet = TRUE)

audit <- data.frame(
  source_person_rows = length(count.fields(args[[3]], skip = 1L)),
  source_unique_crashes = nrow(crashes_all), period_unique_crashes = nrow(crashes),
  urban_easting_max_3116 = urban_easting_max,
  southeast_easting_min_3116 = southeast_easting_min,
  southeast_northing_max_3116 = southeast_northing_max,
  urban_northing_min_3116 = urban_northing_min,
  excluded_disconnected_segments = excluded_disconnected_segments,
  assigned_unique_crashes = sum(segment_count), road_intersections = nrow(supports$intersections),
  road_segments = nrow(supports$segments), line_graph_vertices = vcount(supports$segment_graph),
  line_graph_edges = ecount(supports$segment_graph), road_components = components(supports$road_graph)$no,
  largest_component_segments = sum(supports$segments$analysis_keep),
  zero_count_segments = sum(segment_count == 0), maximum_segment_count = max(segment_count),
  median_assignment_distance_m = median(assigned$assignment_distance_m),
  q95_assignment_distance_m = unname(quantile(assigned$assignment_distance_m, .95)),
  assignments_over_25m = sum(assigned$assignment_distance_m > 25),
  assignments_over_50m = sum(assigned$assignment_distance_m > 50),
  assignments_over_100m = sum(assigned$assignment_distance_m > 100),
  maximum_assignment_distance_m = max(assigned$assignment_distance_m),
  total_segment_length_km = sum(supports$segments$length_m) / 1000
)
write.csv(audit, file.path(output_dir, "bogota_segment_base_audit.csv"), row.names = FALSE)

raw_map <- ggplot() +
  geom_sf(data = supports$segments, color = "grey76", linewidth = .15) +
  geom_sf(data = assigned, color = "#d73027", size = .18, alpha = .35) +
  coord_sf(datum = NA) + theme_void() +
  labs(title = "Unique crash locations and final road-segment network",
       subtitle = paste(start_date, "to", end_date))
ggsave(file.path(output_dir, "01_crashes_and_network.png"), raw_map,
       width = 10, height = 8, dpi = 220, bg = "white")

count_map <- ggplot(supports$segments) +
  geom_sf(aes(color = log1p(count)), linewidth = .42) +
  scale_color_viridis_c(name = "log(1 + crashes)") + coord_sf(datum = NA) +
  theme_void() + labs(title = "Crashes assigned once to the nearest road segment")
ggsave(file.path(output_dir, "02_segment_crash_counts.png"), count_map,
       width = 10, height = 8, dpi = 220, bg = "white")

component_map <- ggplot(supports$segments) +
  geom_sf(aes(color = analysis_keep), linewidth = .32) +
  scale_color_manual(values = c(`TRUE` = "#2166ac", `FALSE` = "grey82"),
                     name = "Largest component") + coord_sf(datum = NA) +
  theme_void() + labs(title = "Road segments retained for the primary connected analysis")
ggsave(file.path(output_dir, "03_segment_components.png"), component_map,
       width = 10, height = 8, dpi = 220, bg = "white")

degree_value <- degree(supports$segment_graph)
supports$segments$line_graph_degree <- degree_value
degree_map <- ggplot(supports$segments) +
  geom_sf(aes(color = line_graph_degree), linewidth = .40) +
  scale_color_viridis_c(name = "Segment neighbours") + coord_sf(datum = NA) +
  theme_void() + labs(title = "Line-graph degree: neighbouring road segments")
ggsave(file.path(output_dir, "04_line_graph_degree.png"), degree_map,
       width = 10, height = 8, dpi = 220, bg = "white")

distance_plot <- ggplot(st_drop_geometry(assigned), aes(assignment_distance_m)) +
  geom_histogram(bins = 60, fill = "#E87722", color = "white") +
  scale_y_continuous(labels = comma) + theme_minimal() +
  labs(title = "Crash-to-segment assignment distance", x = "Distance (m)", y = "Crashes")
ggsave(file.path(output_dir, "05_assignment_distances.png"), distance_plot,
       width = 8, height = 5, dpi = 220, bg = "white")

writeLines(c(
  "# Bogota segment-base report", "",
  paste0("Period: ", start_date, " to ", end_date), "",
  paste0("Eastern urban-network cutoff (EPSG:3116 easting): ", urban_easting_max),
  "The cutoff excludes the rural road extensions toward La Calera and Choachi.", "",
  paste0("Southeast rural exclusion (EPSG:3116): easting > ",
         southeast_easting_min, " and northing < ", southeast_northing_max), "",
  paste0("Southern urban-network cutoff (EPSG:3116 northing): ",
         urban_northing_min), "",
  paste0("Disconnected road segments excluded before modeling: ",
         excluded_disconnected_segments),
  "Only the principal connected road component is present in the final graph.", "",
  paste0("Unique crashes retained and assigned once: ", nrow(assigned)),
  paste0("Road segments (model vertices): ", nrow(supports$segments)),
  paste0("Segment adjacencies (line-graph edges): ", ecount(supports$segment_graph)),
  paste0("Largest connected segment component: ", sum(supports$segments$analysis_keep)),
  paste0("Zero-count segments: ", sum(segment_count == 0)), "",
  paste0("Assignments farther than 25 m (retained, flagged for review): ",
         sum(assigned$assignment_distance_m > 25)),
  paste0("Assignments farther than 50 m (retained, flagged for review): ",
         sum(assigned$assignment_distance_m > 50)), "",
  "Response: crash count per road segment.",
  "Exposure: road-segment length in metres.",
  "Adjacency: two road segments are neighbours when they share an intersection.",
  "Distance flags are diagnostics only: no crash is silently discarded."
), file.path(output_dir, "REPORT.md"))

print(audit, row.names = FALSE)

suppressPackageStartupMessages({
  library(sf); library(sfnetworks); library(igraph); library(tidygraph)
  library(dplyr); library(units)
})

central_localities <- c(
  "es:Chapinero", "es:Barrios Unidos (Bogotá)", "es:Los Mártires",
  "es:Teusaquillo", "es:La Candelaria", "es:Santa Fe (Bogotá)"
)

read_study_polygon <- function(path, urban_easting_max = 1005000,
                               southeast_easting_min = 1003500,
                               southeast_northing_max = 999000,
                               urban_northing_min = 998000) {
  env <- new.env(parent = emptyenv()); load(path, envir = env)
  if (!exists("bog_poly", env, inherits = FALSE)) stop("bog_poly is missing")
  x <- env$bog_poly
  selected <- x[x$wikipedia %in% central_localities, ]
  if (nrow(selected) != length(central_localities)) stop("Not all localities were found")
  selected <- st_transform(st_make_valid(st_union(selected)), 3116)
  bb <- st_bbox(selected)
  if (!is.null(urban_northing_min)) {
    bb[["ymin"]] <- max(bb[["ymin"]], urban_northing_min)
  }
  if (!is.null(urban_easting_max)) {
    bb[["xmax"]] <- min(bb[["xmax"]], urban_easting_max)
  }
  urban_window <- st_as_sfc(bb)
  selected <- suppressWarnings(st_intersection(selected, urban_window))
  # Remove the connected rural southeast spur visible below the urban mesh.
  # It cannot be removed by connected-component filtering because it joins the
  # principal road network near the eastern edge.
  if (!is.null(southeast_easting_min) && !is.null(southeast_northing_max)) {
    cut_bb <- st_bbox(c(
      xmin = southeast_easting_min, ymin = st_bbox(selected)[["ymin"]],
      xmax = st_bbox(selected)[["xmax"]], ymax = southeast_northing_max
    ), crs = st_crs(selected))
    selected <- suppressWarnings(st_difference(selected, st_as_sfc(cut_bb)))
  }
  selected
}

parse_arcgis_date <- function(x) {
  if (inherits(x, c("Date", "POSIXct", "POSIXt"))) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(as.POSIXct(x / 1000, origin = "1970-01-01", tz = "UTC")))
  as.Date(x)
}

read_unique_crashes <- function(path) {
  x <- read.csv(path, check.names = FALSE, colClasses = c(FORMULARIO = "character"))
  required <- c("FORMULARIO", "lon", "lat", "FECHA_HORA_ACC")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Crash CSV is missing: ", paste(missing, collapse = ", "))
  x <- x[is.finite(x$lon) & is.finite(x$lat) & !is.na(x$FORMULARIO) &
           nzchar(x$FORMULARIO), ]
  key <- paste(x$lon, x$lat)
  conflict <- tapply(key, x$FORMULARIO, function(z) length(unique(z)) > 1L)
  if (any(conflict)) stop("Some FORMULARIO values have conflicting coordinates")
  x$analysis_date <- parse_arcgis_date(x$FECHA_HORA_ACC)
  x$persons_involved <- ave(rep.int(1L, nrow(x)), x$FORMULARIO, FUN = length)
  x <- x[!duplicated(x$FORMULARIO), ]
  st_as_sf(x, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}

build_road_network <- function(roads, polygon) {
  roads <- st_transform(roads, 3116)
  polygon <- st_transform(polygon, st_crs(roads))
  roads <- suppressWarnings(st_intersection(roads, polygon))
  roads <- roads[roads$highway %in%
                   c("motorway", "primary", "secondary", "tertiary", "residential"), ]
  roads <- st_make_valid(roads)
  roads <- roads[!st_is_empty(roads) & st_dimension(roads) == 1, ]
  roads <- suppressWarnings(st_collection_extract(roads, "LINESTRING"))
  roads <- suppressWarnings(st_cast(roads, "LINESTRING"))
  roads <- roads[st_geometry_type(roads) == "LINESTRING", ]
  roads$source_length_m <- as.numeric(st_length(roads))
  roads <- roads[is.finite(roads$source_length_m) & roads$source_length_m > 5, ]
  roads$source_row <- seq_len(nrow(roads))
  keep <- intersect(c("source_row", "osm_id", "name", "highway", "geometry"), names(roads))
  net <- as_sfnetwork(roads[, keep], directed = FALSE)
  suppressWarnings(net %>% convert(to_spatial_subdivision))
}

extract_network_supports <- function(net) {
  road_graph <- as.igraph(net)
  intersections <- st_as_sf(net, "nodes")
  segments <- st_as_sf(net, "edges")
  intersections$intersection_id <- seq_len(nrow(intersections))
  segments$segment_id <- seq_len(nrow(segments))
  endpoint <- ends(road_graph, E(road_graph), names = FALSE)
  segments$intersection_from <- endpoint[, 1]
  segments$intersection_to <- endpoint[, 2]
  segments$length_m <- as.numeric(st_length(segments))
  comp <- components(road_graph)
  segments$component <- comp$membership[segments$intersection_from]
  edge_sizes <- tabulate(segments$component, nbins = comp$no)
  largest <- which.max(edge_sizes)
  segments$analysis_keep <- segments$component == largest
  intersections$component <- comp$membership
  intersections$analysis_keep <- intersections$component == largest

  segment_graph <- make_line_graph(road_graph)
  if (vcount(segment_graph) != nrow(segments)) stop("Line graph is not aligned with segments")
  V(segment_graph)$segment_id <- seq_len(vcount(segment_graph))
  pair <- ends(segment_graph, E(segment_graph), names = FALSE)
  adjacency <- data.frame(adjacency_id = seq_len(nrow(pair)),
                          segment_1 = pair[, 1], segment_2 = pair[, 2])
  list(road_graph = road_graph, segment_graph = segment_graph,
       intersections = intersections, segments = segments,
       segment_adjacency = adjacency, largest_component = largest)
}

assign_crashes_to_segments <- function(crashes, segments) {
  crashes <- st_transform(crashes, st_crs(segments))
  nearest <- st_nearest_feature(crashes, segments)
  crashes$segment_id <- segments$segment_id[nearest]
  crashes$assignment_distance_m <- as.numeric(
    st_distance(crashes, segments[nearest, ], by_element = TRUE)
  )
  crashes$assignment_quality <- cut(
    crashes$assignment_distance_m,
    breaks = c(-Inf, 25, 50, Inf),
    labels = c("within_25m", "25_to_50m", "over_50m")
  )
  crashes
}

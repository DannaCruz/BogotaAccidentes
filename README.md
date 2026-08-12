# Bogota road-segment line-graph base

This repository constructs the observational base used by the road-segment
model. A road segment is a vertex of the model graph. Two model vertices are
connected when the corresponding road segments share an intersection.

Each unique crash, identified by `FORMULARIO`, is assigned exactly once to its
nearest final road segment. Rows in the source CSV represent involved persons,
so counting rows directly would count multi-person crashes more than once.

The study polygon is clipped at easting 1,005,000 metres in EPSG:3116 before
filtering either roads or crashes. This explicit and reproducible urban-domain
cutoff removes the rural eastern road extensions toward La Calera and Choachi.
Change `urban_easting_max` near the beginning of the main script if an official
urban-perimeter layer is supplied later.

The connected rural spur at the lower-right of the study area is removed with
an additional EPSG:3116 exclusion window (`easting > 1,003,500` and
`northing < 999,000`). This mask is applied to the observational domain before
filtering both roads and crashes.

A southern cutoff at northing 998,000 metres removes the final road tail below
the urban mesh. The 2023 crash data have no observation south of this cutoff.

After clipping, the code keeps only the principal connected road component and
reconstructs the complete model graph from it. Disconnected road fragments are
therefore absent from the maps, adjacency list and model data—not merely hidden
by the plotting code.

## Inputs

Place or reference:

- `accidentes_bogota.csv`;
- `bogota_roads_sf.rds`;
- `bog_poly.RData`.

## Run

```bash
Rscript install_packages.R

Rscript scripts/build_bogota_segment_base.R \
  data/bogota_roads_sf.rds \
  data/bog_poly.RData \
  data/accidentes_bogota.csv \
  2023-01-01 \
  2023-12-31 \
  outputs/bogota-2023
```

## Model supports

- `road_graph`: intersections are vertices and road segments are edges; used
  only to reconstruct physical topology.
- `segment_graph`: road segments are vertices; its edges mean that two
  segments share an intersection. This is the graph used by the model.
- Response: unique crash count per segment.
- Initial exposure: segment length in metres.

## Outputs

- `bogota_segments_model_base.csv`: response, exposure, covariates and component.
- `bogota_segment_adjacency.csv`: line-graph edge list.
- `bogota_assignments_over_25m.csv`: assignments retained but flagged for
  spatial review; the 25 m threshold is diagnostic, not an exclusion rule.
- `bogota_segment_model_data.rds`: all aligned objects.
- `bogota_segment_spatial_layers.gpkg`: segments, intersections and assigned crashes.
- `bogota_segment_base_audit.csv` and `REPORT.md`.
- Five diagnostic maps/plots numbered in execution order.

The folder `example-output-2023` contains the report, audit table and figures
from a completed 2023 run, without redistributing the source data.

The exported network is the principal connected component; consequently every
exported segment has `analysis_keep == TRUE`.

# Bogota segment-base report

Period: 2023-01-01 to 2023-12-31

Eastern urban-network cutoff (EPSG:3116 easting): 1005000
The cutoff excludes the rural road extensions toward La Calera and Choachi.

Southeast rural exclusion (EPSG:3116): easting > 1003500 and northing < 999000

Southern urban-network cutoff (EPSG:3116 northing): 998000

Disconnected road segments excluded before modeling: 71
Only the principal connected road component is present in the final graph.

Unique crashes retained and assigned once: 2419
Road segments (model vertices): 12674
Segment adjacencies (line-graph edges): 27811
Largest connected segment component: 12674
Zero-count segments: 11050

Assignments farther than 25 m (retained, flagged for review): 88
Assignments farther than 50 m (retained, flagged for review): 42

Response: crash count per road segment.
Exposure: road-segment length in metres.
Adjacency: two road segments are neighbours when they share an intersection.
Distance flags are diagnostics only: no crash is silently discarded.

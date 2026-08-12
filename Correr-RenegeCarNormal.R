

source("funciones.R")
library(Metrics)
library(scales)



#=========================
# 1. CARGA DE BIBLIOTECAS
#=========================
library(sf)
library(sfnetworks)
library(tidygraph)
library(ggraph)
library(dplyr)
library(viridis)
library(Matrix)
library(igraph)
library(MASS)
library(stringr)
library(ggplot2)
library(coda)

set.seed(123)


load("Bogota_network_ready_for_models.RData")

# ============================================================
# PASO 1. Tomar un pedacito pequeño de la red
# ============================================================

library(sf)
library(dplyr)
library(purrr)
library(tibble)

# Tomamos solo los primeros 50 segmentos para probar.
# La idea es NO correr todavía sobre toda la red.
# ============================================================
# PASO 1. Tomar un pedacito espacial continuo
# ============================================================

library(sf)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(tidyr)
library(igraph)

# Aseguramos identificador
roads_all <- roads_ls %>%
  mutate(seg_id = row_number()) %>%
  select(
    seg_id,
    name,
    Y,
    length,
    highway,
    lanes,
    maxspeed,
    oneway,
    lit,
    bridge,
    tunnel,
    cycleway,
    sidewalk,
    bus,
    geometry
  )

# Escogemos un segmento semilla.
# Puedes cambiar este número: 1, 100, 500, 1000, etc.
seed_id <- 1

seed_seg <- roads_all %>%
  filter(seg_id == seed_id)

# Creamos un buffer alrededor del segmento semilla.
# OJO: si el CRS está en metros, 500 significa 500 metros.
# Si está en lon/lat, primero conviene transformar.
st_crs(roads_all)

# Transformamos a un CRS métrico para trabajar con distancias.
# Para Bogotá suele servir MAGNA-SIRGAS / Colombia Bogotá zone EPSG:3116.
roads_m <- st_transform(roads_all, 3116)

seed_seg_m <- roads_m %>%
  filter(seg_id == seed_id)

# Buffer de 500 metros alrededor del segmento
buffer_seed <- st_buffer(seed_seg_m, dist = 1000)

# Seleccionamos los segmentos que intersectan ese buffer
roads_test <- roads_m %>%
  filter(lengths(st_intersects(geometry, buffer_seed)) > 0)

# Revisamos cuántos segmentos quedaron
dim(roads_test)
# Revisamos que sigue siendo un objeto sf
class(roads_test)

# Miramos los primeros registros
head(roads_test)

# ============================================================
# PASO 2. Dibujar pedacito espacial continuo
# ============================================================

ggplot() +
  geom_sf(data = roads_test, linewidth = 0.7, color = "gray35") +
  geom_sf(data = seed_seg_m, linewidth = 1.5, color = "red") +
  geom_sf(data = buffer_seed, fill = NA, color = "blue", linewidth = 0.8) +
  theme_minimal() +
  labs(
    title = "Pedacito espacial continuo de la red vial",
    subtitle = "Rojo: segmento semilla; azul: buffer de selección"
  )

# ============================================================
# PASO 2. Extraer extremos de cada segmento
# ============================================================
# ============================================================
# RECONSTRUIR TODO DESDE roads_test
# ============================================================

library(sf)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(igraph)

# Función para extraer extremos
get_endpoints <- function(geom) {
  
  coords <- st_coordinates(geom)
  
  first_pt <- coords[1, c("X", "Y")]
  last_pt  <- coords[nrow(coords), c("X", "Y")]
  
  tibble(
    x = c(first_pt[1], last_pt[1]),
    y = c(first_pt[2], last_pt[2])
  )
}

# Extraer extremos del roads_test ACTUAL
endpoints_test <- map2_dfr(
  st_geometry(roads_test),
  roads_test$seg_id,
  ~ get_endpoints(.x) %>% mutate(seg_id = .y)
)

# Crear identificador de intersección
endpoints_test <- endpoints_test %>%
  mutate(
    x_round = round(x, 1),
    y_round = round(y, 1),
    inter_id = paste(x_round, y_round, sep = "_")
  )

# Intersecciones compartidas
intersections_shared_test <- endpoints_test %>%
  distinct(inter_id, seg_id) %>%
  group_by(inter_id) %>%
  summarise(
    segmentos = list(unique(seg_id)),
    n_segmentos = length(unique(seg_id)),
    .groups = "drop"
  ) %>%
  filter(n_segmentos >= 2)

# Crear lista de aristas segmento-segmento
edge_list_test <- intersections_shared_test %>%
  mutate(
    pares = map(segmentos, ~ {
      cmb <- combn(.x, 2)
      tibble(
        from = cmb[1, ],
        to   = cmb[2, ]
      )
    })
  ) %>%
  select(inter_id, pares) %>%
  unnest(pares) %>%
  mutate(
    from = as.character(from),
    to   = as.character(to)
  ) %>%
  distinct(from, to, .keep_all = TRUE)

# Crear tabla de vértices
vertices_test <- roads_test %>%
  st_drop_geometry() %>%
  select(seg_id) %>%
  mutate(seg_id = as.character(seg_id))

# Crear grafo
g_test <- graph_from_data_frame(
  d = edge_list_test %>% select(from, to),
  directed = FALSE,
  vertices = vertices_test
)

g_test

# ============================================================
# DIBUJO 1. Ver el pedacito de calles en su geometría real
# ============================================================


# ============================================================
# VERSIÓN ROBUSTA PARA CREAR LAS ARISTAS DEL GRAFO
# ============================================================

nodes_sf_test <- roads_test %>%
  mutate(seg_id = as.character(seg_id)) %>%
  st_centroid()

nodes_xy_test <- nodes_sf_test %>%
  mutate(
    x = st_coordinates(geometry)[, 1],
    y = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry()

edges_xy_test <- edge_list_test %>%
  select(from, to) %>%
  mutate(
    from = as.character(from),
    to   = as.character(to)
  ) %>%
  left_join(
    nodes_xy_test %>%
      rename(from = seg_id, x_from = x, y_from = y),
    by = "from"
  ) %>%
  left_join(
    nodes_xy_test %>%
      rename(to = seg_id, x_to = x, y_to = y),
    by = "to"
  ) %>%
  filter(
    !is.na(x_from),
    !is.na(y_from),
    !is.na(x_to),
    !is.na(y_to)
  )

edges_sf_test <- lapply(seq_len(nrow(edges_xy_test)), function(i) {
  
  st_linestring(
    matrix(
      c(
        edges_xy_test$x_from[i], edges_xy_test$y_from[i],
        edges_xy_test$x_to[i],   edges_xy_test$y_to[i]
      ),
      ncol = 2,
      byrow = TRUE
    )
  )
}) %>%
  st_sfc(crs = st_crs(roads_test)) %>%
  st_sf(
    edges_xy_test,
    geometry = .
  )


ggplot() +
  geom_sf(
    data = roads_test,
    linewidth = 0.8,
    color = "gray75"
  ) +
  geom_sf(
    data = edges_sf_test,
    linewidth = 0.7,
    color = "red",
    alpha = 0.8
  ) +
  geom_sf(
    data = nodes_sf_test,
    size = 1.8,
    color = "blue"
  ) +
  theme_minimal() +
  labs(
    title = "Grafo segmento-segmento sobre el mapa vial",
    subtitle = "Gris: segmentos reales | Azul: nodos del grafo | Rojo: aristas del grafo"
  )




#En esta formulación, cada tramo vial se considera una unidad espacial de análisis.
#La dependencia entre tramos no se define por vecindad areal, 
#sino por conectividad topológica: dos segmentos son vecinos si comparten una intersección.
#De esta forma, las intersecciones actúan como mecanismos de transmisión del riesgo entre segmentos, 
#mientras que el efecto latente RENeGe se define directamente sobre los tramos de la red vial.


# ============================================================
# MATRICES PARA EL PEDACITO roads_test
# ============================================================

library(Matrix)
library(igraph)
library(dplyr)

# ------------------------------------------------------------
# 1. Matriz de adyacencia entre segmentos: W_s
# ------------------------------------------------------------
# W_s[i,j] = 1 si el segmento i y el segmento j comparten intersección
# W_s[i,j] = 0 en caso contrario

W_s_test <- as_adjacency_matrix(
  g_test,
  sparse = TRUE,
  attr = NULL
)

# Revisar dimensión
dim(W_s_test)

# Revisar cuántas conexiones tiene cada segmento
grado_s_test <- Matrix::rowSums(W_s_test)

# ------------------------------------------------------------
# 2. Matriz diagonal de grados: M_s
# ------------------------------------------------------------

M_s_test <- Diagonal(
  x = grado_s_test
)

# ------------------------------------------------------------
# 3. Extraer aristas del grafo de segmentos
# ------------------------------------------------------------

edges_s_test <- igraph::as_data_frame(g_test, what = "edges") %>%
  mutate(edge_id = row_number())

head(edges_s_test)

n_seg_test <- vcount(g_test)
m_edge_test <- nrow(edges_s_test)

cat("Número de segmentos:", n_seg_test, "\n")
cat("Número de aristas segmento-segmento:", m_edge_test, "\n")



# ------------------------------------------------------------
# 4. Matriz de incidencia C
# ------------------------------------------------------------
# C tiene dimensión:
# número de segmentos x número de aristas del grafo

# OJO:
# En g_test los nombres de los nodos son caracteres.
# Necesitamos mapearlos a posiciones 1,...,n.

vertex_names <- V(g_test)$name

id_map <- tibble(
  seg_id = vertex_names,
  row_id = seq_along(vertex_names)
)

edges_s_num_test <- edges_s_test %>%
  left_join(id_map, by = c("from" = "seg_id")) %>%
  rename(i_from = row_id) %>%
  left_join(id_map, by = c("to" = "seg_id")) %>%
  rename(i_to = row_id)

head(edges_s_num_test)

# Construcción de C
C_test <- sparseMatrix(
  i = c(edges_s_num_test$i_from, edges_s_num_test$i_to),
  j = c(edges_s_num_test$edge_id, edges_s_num_test$edge_id),
  x = 1,
  dims = c(n_seg_test, m_edge_test)
)

dim(C_test)

# ------------------------------------------------------------
# 6. Grafo de aristas: line graph
# ------------------------------------------------------------
# Dos aristas son vecinas si comparten un segmento.

g_edge_test <- make_line_graph(g_test)

g_edge_test



# ------------------------------------------------------------
# 7. Matriz A^e
# ------------------------------------------------------------

A_e_test <- as_adjacency_matrix(
  g_edge_test,
  sparse = TRUE,
  attr = NULL
)

dim(A_e_test)


# ------------------------------------------------------------
# 8. Matriz M^e
# ------------------------------------------------------------

grado_e_test <- Matrix::rowSums(A_e_test)

M_e_test <- Diagonal(
  x = grado_e_test
)

dim(M_e_test)

summary(grado_e_test)


# ------------------------------------------------------------
# 9. Objeto final con matrices del pedacito
# ------------------------------------------------------------

matrices_test <- list(
  W_s = W_s_test,     # vecindad entre segmentos
  M_s = M_s_test,     # diagonal de grados de segmentos
  C   = C_test,       # incidencia segmento-arista
  A_e = A_e_test,     # vecindad entre aristas del grafo
  M_e = M_e_test      # diagonal de grados de aristas
)

matrices_test
#=========================
# 3. FUNCIONES AUXILIARES
#=========================



# ============================================================
# DATOS PARA STAN - PEDACITO
# ============================================================

library(Matrix)
library(cmdstanr)

# Asegurar que Y sea entero no negativo
Y_test <- as.integer(roads_test$Y)

# Offset por longitud
# Le ponemos una cota mínima para evitar log(0)
offset_test <- log(pmax(as.numeric(roads_test$length), 1e-6))

# Dimensiones
N_test <- nrow(roads_test)
E_test <- ncol(C_test)

# Convertimos matrices sparse a matrices densas
# Esto está bien para el pedacito.
C_dense_test <- as.matrix(C_test)
A_e_dense_test <- as.matrix(A_e_test)

# Grado de las aristas del grafo de aristas
degree_e_test <- as.numeric(Matrix::rowSums(A_e_test))



# ============================================================
# Covariables
# ============================================================




covariates_test <- roads_test %>%
  st_drop_geometry() %>%
  mutate(
    highway = ifelse(is.na(highway), "unknown", highway),
    
    lanes_num = as.numeric(lanes),
    maxspeed_num = as.numeric(str_extract(maxspeed, "\\d+")),
    
    is_oneway = ifelse(oneway %in% c("yes", "true", "1"), 1, 0),
    is_lit = ifelse(lit %in% c("yes", "true", "1"), 1, 0),
    is_bridge = ifelse(bridge %in% c("yes", "true", "1"), 1, 0),
    is_tunnel = ifelse(tunnel %in% c("yes", "true", "1"), 1, 0),
    
    has_cycleway = ifelse(!is.na(cycleway), 1, 0),
    has_sidewalk = ifelse(!is.na(sidewalk), 1, 0),
    has_bus = ifelse(!is.na(bus), 1, 0),
    
    degree_segment = as.numeric(grado_s_test)
  ) %>%
  select(
    highway,
    lanes_num,
    maxspeed_num,
    is_oneway,
    is_lit,
    is_bridge,
    is_tunnel,
    has_cycleway,
    has_sidewalk,
    has_bus,
    degree_segment
  )


covariates_test <- covariates_test %>%
  mutate(
    lanes_num = ifelse(
      is.na(lanes_num),
      median(lanes_num, na.rm = TRUE),
      lanes_num
    ),
    maxspeed_num = ifelse(
      is.na(maxspeed_num),
      median(maxspeed_num, na.rm = TRUE),
      maxspeed_num
    )
  )

X_test <- model.matrix(
  ~ highway +
    lanes_num +
    maxspeed_num +
    is_oneway +
    is_lit +
    is_bridge +
    is_tunnel +
    has_cycleway +
    has_sidewalk +
    has_bus +
    degree_segment,
  data = covariates_test
)

# Quitar intercepto
X_test <- X_test[, -1, drop = FALSE]

# Quitar columnas constantes
sd_cols <- apply(X_test, 2, sd)
X_test <- X_test[, sd_cols > 0, drop = FALSE]

# Estandarizar
X_test <- scale(X_test)

# Convertir a matriz normal
X_test <- as.matrix(X_test)

# Actualizar K
K_test <- ncol(X_test)




stan_data_test <- list(
  N = N_test,
  E = E_test,
  K = K_test,
  Y = Y_test,
  log_offset = offset_test,
  X = X_test,
  C = C_dense_test,
  A_e = A_e_dense_test,
  degree_e = degree_e_test
)

mod_renege_seg <- cmdstan_model("renege_segmentos_nb.stan")

fit_renege_seg <- mod_renege_seg$sample(
  data = stan_data_test,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123
)

saveRDS(fit_renege_seg, file = "fit_renege_seg.rds")

#load("fit_renege_seg.rds")



# ============================================================
# DATOS PARA MODELO CAR SOBRE SEGMENTOS
# ============================================================

W_dense_test <- as.matrix(W_s_test)
degree_s_test <- as.numeric(Matrix::rowSums(W_s_test))

stan_data_car_test <- list(
  N = N_test,
  K = K_test,
  Y = Y_test,
  log_offset = offset_test,
  X = X_test,
  W = W_dense_test,
  degree = degree_s_test
)

mod_car_seg <- cmdstan_model("car_segmentos_nb.stan")

fit_car_seg <- mod_car_seg$sample(
  data = stan_data_car_test,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123
)

saveRDS(fit_car_seg, file = "fit_car_seg.rds")



fit_renege_seg$summary(c("alpha", "beta", "sigma_rho", "gamma", "phi"))

beta_summary <- fit_renege_seg$summary("beta")

beta_results <- beta_summary %>%
  mutate(covariable = colnames(X_test)) %>%
  select(covariable, mean, sd, q5,q95, rhat, ess_bulk)

beta_results

beta_results <- beta_results %>%
  mutate(
    exp_mean = exp(mean),
    exp_q5 = exp(q5),
    exp_q95 = exp(q95)
  )

beta_results

lambda_summary <- fit_renege_seg$summary("lambda")

comparacion_test <- data.frame(
  seg_id = roads_test$seg_id,
  Y_obs = Y_test,
  lambda_hat = lambda_summary$mean,
  length = roads_test$length
)

head(comparacion_test)
# 
# ggplot(comparacion_test, aes(x = Y_obs, y = lambda_hat)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, color = "red") +
#   theme_minimal() +
#   labs(
#     x = "Accidentes observados",
#     y = "Accidentes estimados",
#     title = "Accidentes observados vs estimados"
#   )



theta_summary <- fit_renege_seg$summary("theta")

theta_results <- data.frame(
  seg_id = roads_test$seg_id,
  theta_mean = theta_summary$mean,
  theta_q5 = theta_summary$q5,
  theta_q95 = theta_summary$q95
)

head(theta_results)

# 
# ggplot(roads_effect_test) +
#   geom_sf(aes(color = theta_mean), linewidth = 1.3) +
#   scale_color_gradient2(
#     low = "#2c7bb6",
#     mid = "gray90",
#     high = "#d7191c",
#     midpoint = 0
#   ) +
#   theme_minimal() +
#   labs(
#     title = "Efecto latente de red sobre segmentos viales",
#     subtitle = "Valores positivos indican mayor riesgo latente",
#     color = "theta"
#   )


fit_renege_seg$diagnostic_summary()


lambda_summary <- fit_renege_seg$summary("lambda")

lambda_results <- data.frame(
  seg_id = roads_test$seg_id,
  Y_obs = Y_test,
  lambda_mean = lambda_summary$mean,
  lambda_q5 = lambda_summary$q5,
  lambda_q95 = lambda_summary$q95,
  length = roads_test$length
)

head(lambda_results)
# 
# 
# ggplot(roads_lambda_test) +
#   geom_sf(aes(color = lambda_mean), linewidth = 1.3) +
#   scale_color_viridis_c(option = "inferno") +
#   theme_minimal() +
#   labs(
#     title = "Accidentes esperados por segmento",
#     subtitle = "Media posterior de lambda",
#     color = "lambda"
#   )

######Para mapear el efecto por segmento,


theta_summary <- fit_renege_seg$summary("theta")

roads_effect_test <- roads_test %>%
  mutate(
    theta_mean = theta_summary$mean,
    theta_q5 = theta_summary$q5,
    theta_q95 = theta_summary$q95
  )

# ggplot(roads_effect_test) +
#   geom_sf(aes(color = theta_mean), linewidth = 1.2) +
#   scale_color_viridis_c(option = "plasma") +
#   theme_minimal() +
#   labs(
#     title = "Efecto latente de red sobre segmentos",
#     subtitle = "theta = C * rho",
#     color = "theta"
#   )


#También puedes mapear lambda, que es el número esperado total de accidentes:

lambda_summary <- fit_renege_seg$summary("lambda")

roads_lambda_test <- roads_test %>%
  mutate(
    lambda_mean = lambda_summary$mean,
    lambda_q5 = lambda_summary$q5,
    lambda_q95 = lambda_summary$q95
  )

# ggplot(roads_lambda_test) +
#   geom_sf(aes(color = lambda_mean), linewidth = 1.2) +
#   scale_color_viridis_c(option = "inferno") +
#   theme_minimal() +
#   labs(
#     title = "Accidentes esperados por segmento",
#     subtitle = "lambda posterior media",
#     color = "lambda"
#   )

#Para el efecto sobre las conexiones/intersecciones del grafo, usa rho:


rho_summary <- fit_renege_seg$summary("rho")

edges_effect_test <- edges_sf_test %>%
  mutate(
    rho_mean = rho_summary$mean,
    rho_q5 = rho_summary$q5,
    rho_q95 = rho_summary$q95
  )


# ggplot() +
#   geom_sf(
#     data = roads_test,
#     linewidth = 0.7,
#     color = "gray80"
#   ) +
#   geom_sf(
#     data = edges_effect_test,
#     aes(color = rho_mean),
#     linewidth = 1.1,
#     alpha = 0.9
#   ) +
#   geom_sf(
#     data = nodes_sf_test,
#     size = 1.5,
#     color = "black"
#   ) +
#   scale_color_viridis_c(option = "magma") +
#   theme_minimal() +
#   labs(
#     title = "Efecto latente sobre conexiones del grafo",
#     subtitle = "rho vive sobre las aristas segmento-segmento",
#     color = "rho"
#   )
# 
# 
# 
# ggplot() +
#   geom_sf(
#     data = roads_effect_test,
#     aes(color = theta_mean),
#     linewidth = 1.2
#   ) +
#   geom_sf(
#     data = edges_effect_test,
#     aes(linewidth = abs(rho_mean)),
#     color = "gray20",
#     alpha = 0.12
#   ) +
#   scale_color_gradient2(
#     low = "#2c7bb6",
#     mid = "gray90",
#     high = "#d7191c",
#     midpoint = 0
#   ) +
#   scale_linewidth_continuous(range = c(0.05, 0.6)) +
#   theme_minimal() +
#   labs(
#     title = "Efecto de red sobre segmentos y conexiones",
#     subtitle = "Color: theta en segmentos | líneas grises: conexiones del grafo",
#     color = "theta",
#     linewidth = "|rho|"
#   )





# tasa esperada por metro 

roads_lambda_test <- roads_lambda_test %>%
  mutate(
    lambda_per_km = lambda_mean / as.numeric(length) * 1000
  )

# ggplot(roads_lambda_test) +
#   geom_sf(aes(color = lambda_per_km), linewidth = 1.3) +
#   scale_color_viridis_c(option = "inferno") +
#   theme_minimal() +
#   labs(
#     title = "Tasa esperada de accidentes por segmento",
#     subtitle = "Accidentes esperados por kilómetro",
#     color = "lambda/km"
#   )

comparacion_efectos <- roads_effect_test %>%
  st_drop_geometry() %>%
  select(seg_id, theta_mean) %>%
  left_join(
    roads_lambda_test %>%
      st_drop_geometry() %>%
      select(seg_id, lambda_mean, lambda_per_km),
    by = "seg_id"
  )

ggplot(comparacion_efectos, aes(x = theta_mean, y = lambda_per_km)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "theta",
    y = "lambda por km",
    title = "Relación entre efecto latente y tasa esperada"
  )

#####




####

calles_riesgo <- roads_test %>%
  st_drop_geometry() %>%
  mutate(
    lambda_mean = lambda_summary$mean,
    lambda_per_km = lambda_mean / as.numeric(length) * 1000,
    theta_mean = theta_summary$mean
  ) %>%
  select(
    seg_id,
    name,
    highway,
    Y,
    length,
    lambda_mean,
    lambda_per_km,
    theta_mean
  )


calles_riesgo %>%
  arrange(desc(lambda_per_km)) %>%
  head(20)


calles_riesgo_nombre <- calles_riesgo %>%
  filter(!is.na(name), name != "") %>%
  group_by(name) %>%
  summarise(
    n_segmentos = n(),
    accidentes_obs = sum(Y, na.rm = TRUE),
    longitud_total_m = sum(as.numeric(length), na.rm = TRUE),
    lambda_total = sum(lambda_mean, na.rm = TRUE),
    lambda_por_km = lambda_total / longitud_total_m * 1000,
    theta_promedio = mean(theta_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(lambda_por_km))



calles_riesgo_nombre %>%
  head(20)


#Si quieres la calle con más riesgo latente de red, usa:
calles_riesgo_nombre %>%
  arrange(desc(theta_promedio)) %>%
  head(20)

#Si quieres la calle con más accidentes esperados totales, usa:
calles_riesgo_nombre %>%
  arrange(desc(lambda_total)) %>%
  head(20)

# ============================================================
# GUARDAR MAPAS DEL MODELO EN PNG
# ============================================================

# Crear carpeta para guardar las figuras
dir.create("figures", showWarnings = FALSE)

# ============================================================
# 1. MAPA: TASA ESPERADA DE ACCIDENTES POR KM
# ============================================================





# ============================================================
# COMPONENTES DEL MODELO
# ============================================================

# Resumen posterior de lambda: accidentes esperados por segmento
lambda_summary <- fit_renege_seg$summary("lambda")

# Resumen posterior de theta: efecto latente de red por segmento
theta_summary <- fit_renege_seg$summary("theta")

# Media posterior de beta: coeficientes de covariables
beta_mean <- fit_renege_seg$summary("beta")$mean

# Media posterior de alpha: intercepto
alpha_mean <- fit_renege_seg$summary("alpha")$mean

# Componente lineal de covariables: X beta
xbeta_mean <- as.numeric(X_test %*% beta_mean)

# Tabla con componentes del modelo por segmento
componentes_test <- data.frame(
  seg_id = roads_test$seg_id,
  name = roads_test$name,
  Y_obs = Y_test,
  length = as.numeric(roads_test$length),
  theta_mean = theta_summary$mean,
  xbeta_mean = xbeta_mean,
  lambda_mean = lambda_summary$mean
) %>%
  mutate(
    lambda_per_km = lambda_mean / length * 1000,
    eta_sin_offset = alpha_mean + xbeta_mean + theta_mean
  )

head(componentes_test)

# Agregamos la tasa esperada por km al objeto espacial
roads_total_test <- roads_test %>%
  mutate(lambda_per_km = componentes_test$lambda_per_km)

# Creamos el mapa de tasa esperada de accidentes por km
p_lambda_km <- ggplot(roads_total_test) +
  geom_sf(aes(color = lambda_per_km), linewidth = 1.3) +
  scale_color_viridis_c(option = "inferno") +
  theme_minimal() +
  labs(
    title = "Expected Accident Rate",
    subtitle = "Posterior mean of expected accidents per kilometer",
    color = "Expected accidents/km"
  )

# Guardamos el mapa en PNG
ggsave(
  filename = "figures/map_expected_accident_rate.png",
  plot = p_lambda_km,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# 2. MAPA: EFECTO DE COVARIABLES
# ============================================================

# Agregamos el componente X beta al objeto espacial
roads_xbeta_test <- roads_test %>%
  mutate(xbeta_mean = componentes_test$xbeta_mean)

# Creamos el mapa del efecto de covariables
p_xbeta <- ggplot(roads_xbeta_test) +
  geom_sf(aes(color = xbeta_mean), linewidth = 1.3) +
  scale_color_gradient2(
    low = "#2c7bb6",
    mid = "gray90",
    high = "#d7191c",
    midpoint = 0
  ) +
  theme_minimal() +
  labs(
    title = "Covariate Effect on Accident Risk",
    subtitle = "Linear predictor component: X beta",
    color = "X beta"
  )

# Guardamos el mapa en PNG
ggsave(
  filename = "figures/map_covariate_effect.png",
  plot = p_xbeta,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# 3. MAPA: EFECTO LATENTE DE RED
# ============================================================

# Agregamos el efecto latente theta al objeto espacial
roads_theta_test <- roads_test %>%
  mutate(theta_mean = componentes_test$theta_mean)

# Creamos el mapa del efecto latente de red
p_theta <- ggplot(roads_theta_test) +
  geom_sf(aes(color = theta_mean), linewidth = 1.3) +
  scale_color_gradient2(
    low = "#2c7bb6",
    mid = "gray90",
    high = "#d7191c",
    midpoint = 0
  ) +
  theme_minimal() +
  labs(
    title = "Latent Network Effect",
    subtitle = "Segment-level latent component: theta",
    color = "theta"
  )

# Guardamos el mapa en PNG
ggsave(
  filename = "figures/map_latent_network_effect.png",
  plot = p_theta,
  width = 9,
  height = 7,
  dpi = 300
)




# ============================================================
# MAPA: ACCIDENTES OBSERVADOS
# ============================================================

# Agregamos los accidentes observados al objeto espacial
roads_observed_test <- roads_test %>%
  mutate(
    observed_accidents = Y_test,
    observed_accidents_per_km = Y_test / as.numeric(length) * 1000
  )

# Mapa de conteo observado de accidentes por segmento
p_observed_counts <- ggplot(roads_observed_test) +
  geom_sf(aes(color = observed_accidents), linewidth = 1.3) +
  scale_color_viridis_c(option = "inferno") +
  theme_minimal() +
  labs(
    title = "Observed Accidents",
    subtitle = "Observed accident counts by road segment",
    color = "Observed accidents"
  )

# Guardar mapa en PNG
ggsave(
  filename = "figures/map_observed_accidents.png",
  plot = p_observed_counts,
  width = 9,
  height = 7,
  dpi = 300
)



# ============================================================
# MAPA: TASA OBSERVADA DE ACCIDENTES POR KM
# ============================================================

p_observed_rate <- ggplot(roads_observed_test) +
  geom_sf(aes(color = observed_accidents_per_km), linewidth = 1.3) +
  scale_color_viridis_c(option = "inferno") +
  theme_minimal() +
  labs(
    title = "Observed Accident Rate",
    subtitle = "Observed accidents per kilometer by road segment",
    color = "Observed accidents/km"
  )

# Guardar mapa en PNG
ggsave(
  filename = "figures/map_observed_accident_rate.png",
  plot = p_observed_rate,
  width = 9,
  height = 7,
  dpi = 300
)

#####
#CAR
########

lambda_car_summary <- fit_car_seg$summary("lambda")

comparacion_car <- data.frame(
  seg_id = roads_test$seg_id,
  Y_obs = Y_test,
  lambda_car = lambda_car_summary$mean,
  length = as.numeric(roads_test$length)
) %>%
  mutate(lambda_car_per_km = lambda_car / length * 1000)

head(comparacion_car)


roads_car_test <- roads_test %>%
  mutate(lambda_car_per_km = comparacion_car$lambda_car_per_km)

ggplot(roads_car_test) +
  geom_sf(aes(color = lambda_car_per_km), linewidth = 1.3) +
  scale_color_viridis_c(option = "inferno") +
  theme_minimal() +
  labs(
    title = "CAR Expected Accident Rate",
    subtitle = "Posterior mean of expected accidents per kilometer",
    color = "Expected accidents/km"
  )





#######

# MEDIDAS
###############


library(loo)
library(posterior)
library(dplyr)


log_lik_renege <- fit_renege_seg$draws("log_lik", format = "matrix")
log_lik_car <- fit_car_seg$draws("log_lik", format = "matrix")


waic_renege <- waic(log_lik_renege)
waic_car <- waic(log_lik_car)

waic_renege
waic_car

loo_compare(waic_renege, waic_car)

loo_renege <- loo(log_lik_renege)
loo_car <- loo(log_lik_car)

loo_renege
loo_car

loo_compare(loo_renege, loo_car)


calc_dic_from_loglik <- function(log_lik_matrix) {
  # log_lik_matrix: draws x N
  
  # Log-verosimilitud total por iteración
  loglik_total <- rowSums(log_lik_matrix)
  
  # Deviance por iteración
  deviance <- -2 * loglik_total
  
  # D_bar
  D_bar <- mean(deviance)
  
  # Aproximación de D_hat usando la media posterior de log-lik por observación
  # Es una aproximación, no el DIC exacto clásico.
  D_hat <- -2 * sum(colMeans(log_lik_matrix))
  
  # Número efectivo de parámetros
  p_D <- D_bar - D_hat
  
  # DIC
  DIC <- D_bar + p_D
  
  data.frame(
    D_bar = D_bar,
    D_hat = D_hat,
    p_D = p_D,
    DIC = DIC
  )
}

dic_renege <- calc_dic_from_loglik(log_lik_renege)
dic_car <- calc_dic_from_loglik(log_lik_car)

dic_renege
dic_car

comparacion_modelos <- bind_rows(
  data.frame(
    modelo = "RENeGe",
    waic = waic_renege$estimates["waic", "Estimate"],
    p_waic = waic_renege$estimates["p_waic", "Estimate"],
    elpd_loo = loo_renege$estimates["elpd_loo", "Estimate"],
    p_loo = loo_renege$estimates["p_loo", "Estimate"],
    looic = loo_renege$estimates["looic", "Estimate"],
    DIC = dic_renege$DIC
  ),
  data.frame(
    modelo = "CAR",
    waic = waic_car$estimates["waic", "Estimate"],
    p_waic = waic_car$estimates["p_waic", "Estimate"],
    elpd_loo = loo_car$estimates["elpd_loo", "Estimate"],
    p_loo = loo_car$estimates["p_loo", "Estimate"],
    looic = loo_car$estimates["looic", "Estimate"],
    DIC = dic_car$DIC
  )
)

comparacion_modelos



# ============================================================
# MAPA CAR: TASA ESPERADA DE ACCIDENTES POR KM
# ============================================================

# Crear carpeta para guardar figuras
dir.create("figures", showWarnings = FALSE)

# Extraer lambda del modelo CAR
lambda_car_summary <- fit_car_seg$summary("lambda")

# Crear tabla con tasa esperada por km para el modelo CAR
comparacion_car <- data.frame(
  seg_id = roads_test$seg_id,
  Y_obs = Y_test,
  lambda_car = lambda_car_summary$mean,
  length = as.numeric(roads_test$length)
) %>%
  mutate(
    lambda_car_per_km = lambda_car / length * 1000
  )

# Agregar la tasa esperada CAR al objeto espacial
roads_car_test <- roads_test %>%
  mutate(
    lambda_car_per_km = comparacion_car$lambda_car_per_km
  )

# Crear mapa CAR
p_car_lambda_km <- ggplot(roads_car_test) +
  geom_sf(aes(color = lambda_car_per_km), linewidth = 1.3) +
  scale_color_viridis_c(option = "inferno") +
  theme_minimal() +
  labs(
    title = "CAR Expected Accident Rate",
    subtitle = "Posterior mean of expected accidents per kilometer",
    color = "Expected accidents/km"
  )

# Guardar mapa en PNG
ggsave(
  filename = "figures/map_car_expected_accident_rate.png",
  plot = p_car_lambda_km,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# MAPA RENeGe: EFECTO LATENTE SOBRE CONEXIONES / INTERSECCIONES
# ============================================================

# Crear carpeta para guardar figuras
dir.create("figures", showWarnings = FALSE)

# Extraer resumen posterior de rho
rho_summary <- fit_renege_seg$summary("rho")

# Agregar rho a las aristas del grafo segmento-segmento
edges_renege_test <- edges_sf_test %>%
  mutate(
    rho_mean = rho_summary$mean,
    rho_q5 = rho_summary$q5,
    rho_q95 = rho_summary$q95,
    rho_abs = abs(rho_mean)
  )

# Mapa del efecto RENeGe sobre conexiones
p_renege_connections <- ggplot() +
  geom_sf(
    data = roads_test,
    linewidth = 0.6,
    color = "gray85"
  ) +
  geom_sf(
    data = edges_renege_test,
    aes(color = rho_mean, linewidth = rho_abs),
    alpha = 0.85
  ) +
  scale_color_gradient2(
    low = "#2c7bb6",
    mid = "gray90",
    high = "#d7191c",
    midpoint = 0
  ) +
  scale_linewidth_continuous(
    range = c(0.1, 1.2)
  ) +
  theme_minimal() +
  labs(
    title = "RENeGe Latent Effect on Network Connections",
    subtitle = "Posterior mean of rho on segment-to-segment connections",
    color = "rho",
    linewidth = "|rho|"
  )

# Guardar mapa en PNG
ggsave(
  filename = "figures/map_renege_network_connections.png",
  plot = p_renege_connections,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# MAPA RENeGe: PUNTOS DE CONEXIÓN CON EFECTO RHO
# ============================================================

# Calcular punto medio de cada conexión del grafo
edges_renege_points <- edges_renege_test %>%
  st_centroid()

# Mapa de puntos RENeGe
p_renege_connection_points <- ggplot() +
  geom_sf(
    data = roads_test,
    linewidth = 0.7,
    color = "gray85"
  ) +
  geom_sf(
    data = edges_renege_points,
    aes(color = rho_mean, size = rho_abs),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "#2c7bb6",
    mid = "gray90",
    high = "#d7191c",
    midpoint = 0
  ) +
  scale_size_continuous(
    range = c(0.5, 3.5)
  ) +
  theme_minimal() +
  labs(
    title = "RENeGe Latent Effect at Network Connections",
    subtitle = "Point representation of posterior mean rho",
    color = "rho",
    size = "|rho|"
  )

# Guardar mapa en PNG
ggsave(
  filename = "figures/map_renege_connection_points.png",
  plot = p_renege_connection_points,
  width = 9,
  height = 7,
  dpi = 300
)



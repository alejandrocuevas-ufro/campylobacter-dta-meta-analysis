# ============================================================
# 01_prep_datos.R
# Preparación de datasets DTA Campylobacter (inmunoensayos)
# Adaptado de: Datasets_R_Campy_v2_estructura_nueva.R (autor del proyecto)
# Cambios: rutas de salida al workspace (la carpeta del proyecto es sólo lectura)
# ============================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr)
})

# Ruta del libro de extracción. Orden de preferencia:
#   1. variable de entorno CAMPY_DATA
#   2. data/DATA_METAANALYSIS.xlsx dentro del repositorio
#   3. la carpeta de trabajo original del proyecto
ROOT_DIR  <- "/Users/jano/Proyectos R/METAANALYSIS_CAMPY"
DATA_FILE <- Sys.getenv("CAMPY_DATA", unset = "")
if (!nzchar(DATA_FILE) || !file.exists(DATA_FILE)) {
  cand <- c(file.path("data", "DATA_METAANALYSIS.xlsx"),
            file.path(ROOT_DIR, "DATA_METAANALYSIS.xlsx"))
  DATA_FILE <- cand[file.exists(cand)][1]
  if (is.na(DATA_FILE)) stop("No se encontró DATA_METAANALYSIS.xlsx; define CAMPY_DATA.")
}
OUT_DIR   <- "salidas"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

INCLUDE_ELIGIBILITY_CAUTION <- TRUE
INCLUDE_ELIGIBILITY_NA      <- TRUE
COMPARISON_LEVELS <- c("LFA_Culture", "LFA_PCR", "EIA_Culture", "EIA_PCR")

quadas_cols <- c("rob_PS","rob_IT","rob_RS","rob_FT","ac_PS","ac_IT","ac_RS")

standardise_schema <- function(x) {
  if (!"Index_1" %in% names(x) && "Index" %in% names(x)) x$Index_1 <- x$Index
  if (!"Reference_1" %in% names(x) && "Reference" %in% names(x)) x$Reference_1 <- x$Reference
  req <- c("author","year","TP","FN","FP","TN","Index_1","Index_2","Reference_2")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("Faltan columnas en MetaDTA: ", paste(miss, collapse=", "))
  num <- intersect(c("year","TP","FN","FP","TN","n","sample_prev_%", quadas_cols), names(x))
  for (nm in num) x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
  x$Record_ID <- seq_len(nrow(x))
  x
}

make_study_id <- function(author, year, record_id) {
  base <- str_squish(str_remove(as.character(author), "\\s*\\([^)]*\\)\\s*$"))
  yr <- as.character(year)
  has_year <- !is.na(yr) & str_detect(base, paste0("\\b", yr, "\\b"))
  out <- ifelse(has_year | is.na(year), base, paste(base, year))
  bad <- is.na(base) | base == "" | is.na(out) | out == "NA"
  out[bad] <- paste0("row_", record_id[bad]); out
}

make_comparison <- function(index_2, reference_2) {
  case_when(
    index_2 == "lateral flow immunoassay" & reference_2 == "Culture" ~ "LFA_Culture",
    index_2 == "lateral flow immunoassay" & reference_2 == "PCR"     ~ "LFA_PCR",
    index_2 == "microplate assay"         & reference_2 == "Culture" ~ "EIA_Culture",
    index_2 == "microplate assay"         & reference_2 == "PCR"     ~ "EIA_PCR",
    TRUE ~ NA_character_)
}

market_rank_test <- function(index_1, index_2) {
  idx <- str_to_lower(as.character(index_1)); platform <- as.character(index_2)
  case_when(
    platform == "microplate assay" & str_detect(idx, "prospect") ~ 1L,
    platform == "microplate assay" & str_detect(idx, "ridascreen|rida screen") ~ 2L,
    platform == "microplate assay" & str_detect(idx, "premier") ~ 3L,
    platform == "microplate assay" ~ 99L,
    platform == "lateral flow immunoassay" & str_detect(idx, "immunocard|stat") ~ 1L,
    platform == "lateral flow immunoassay" & str_detect(idx, "ridaquick|rida quick") ~ 2L,
    platform == "lateral flow immunoassay" ~ 99L,
    TRUE ~ 99L)
}

eligible_filter <- function(x) {
  if (!"eligibility" %in% names(x)) return(rep(TRUE, nrow(x)))
  ok <- x$eligibility == 1
  if (INCLUDE_ELIGIBILITY_CAUTION) ok <- ok | x$eligibility == 2
  if (INCLUDE_ELIGIBILITY_NA)      ok <- ok | is.na(x$eligibility)
  ok[is.na(ok)] <- FALSE; ok
}

add_number_id <- function(x) {
  x %>% ungroup() %>% arrange(Study_ID, market_rank, Record_ID) %>%
    mutate(Number_ID = row_number(), study_label = paste0(Study_ID, " | ", Index_1))
}

select_one_per_study <- function(x) {
  x %>% group_by(Comparison, Study_ID) %>%
    arrange(market_rank, eligibility_rank, desc(n), Record_ID, .by_group = TRUE) %>%
    mutate(selected_for_global = row_number() == 1L,
           selection_order = row_number(),
           n_entries_same_study = n(),
           selection_rule = if_else(selected_for_global,
             paste0("Seleccionado por ranking de mercado: ", Index_1, " (rank=", market_rank, ")"),
             paste0("No seleccionado: existe un test con mayor prioridad de mercado en ", Study_ID))) %>%
    ungroup()
}

get_dataset <- function(x, comparison, selected = FALSE) {
  out <- filter(x, Comparison == comparison)
  if (selected) out <- filter(out, selected_for_global)
  add_number_id(out)
}

datos_raw <- data.frame(read_excel(DATA_FILE, sheet = "MetaDTA"), check.names = FALSE)

# ------------------------------------------------------------------
# Correcciones documentadas tras la conciliación de referencias
# (ver conciliacion_referencias.xlsx y notas_para_revision.md)
#  - Zaghloul 2012: se retira. El estándar de referencia del artículo es
#    compuesto (cultivo positivo, o cultivo negativo con PCR e inmunoensayo
#    positivos), de modo que incorpora la prueba índice; además la revista
#    no está indexada (sin DOI ni PMID) y el estudio nunca pasó la selección.
#  - Özcan 2022 y Tachikawa 2017: se incorporan desde patch_datos.csv.
#    Tachikawa queda marcado solo_sensibilidad = 1 (2x2 reconstruido del
#    resumen, QUADAS-2 no evaluable) y no entra en la síntesis principal.
# ------------------------------------------------------------------
datos_raw <- datos_raw[!grepl("Zaghloul", datos_raw$author), , drop = FALSE]
PATCH_FILE <- "patch_datos.csv"
if (file.exists(PATCH_FILE)) {
  # locale "C" en este entorno: leer sin fileEncoding y marcar UTF-8 a mano
  pt <- read.csv(PATCH_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in names(pt)) if (is.character(pt[[nm]])) Encoding(pt[[nm]]) <- "UTF-8"
  flag_pt <- ifelse(is.na(pt$solo_sensibilidad), 0, pt$solo_sensibilidad)
  for (nm in setdiff(names(datos_raw), names(pt))) pt[[nm]] <- NA
  pt <- pt[, names(datos_raw), drop = FALSE]
  n_base <- nrow(datos_raw)
  datos_raw <- rbind(datos_raw, pt)
  datos_raw$solo_sensibilidad <- c(rep(0, n_base), flag_pt)
} else {
  datos_raw$solo_sensibilidad <- 0
}

# ------------------------------------------------------------------
# Reglas de armonización acordadas con el equipo revisor (octubre 2026)
#
# (1) Muestras preseleccionadas -> riesgo de sesgo ALTO en selección de
#     pacientes. La planilla original calificaba de forma desigual el mismo
#     diseño: Granato 2010, Platts-Mills 2014 y Tolcin 2000 como incierto, y
#     Franco 2021, Kawatsu 2008 y Schnee 2018 como bajo, pese a que en todos
#     ellos el muestreo fue selectivo (preselected = 1, sampling_type = 2).
#     QUADAS-2 califica como riesgo alto el muestreo que no es consecutivo ni
#     aleatorio, de modo que se unifica el criterio para todas las entradas.
#     Afecta también a Lin 2026, que llega así desde patch_datos.csv.
#
# (2) Campo country: recoge en algunos registros el país del autor de
#     correspondencia y no el del estudio, y mezcla etiquetas en español y en
#     inglés. Se normaliza a inglés y se corrige Platts-Mills 2014, cuyas
#     muestras proceden de Tanzania, Bangladés y Perú.
# ------------------------------------------------------------------
pres <- !is.na(datos_raw$preselected) & datos_raw$preselected == 1
datos_raw$rob_PS[pres] <- 2
message("preselected -> rob_PS = alto: ", sum(pres), " entradas")

norm_pais <- function(x) {
  x <- as.character(x)
  x[grepl("^Alemania", x)]  <- "Germany"
  x[grepl("^Espa",     x)]  <- "Spain"
  x[grepl("^Finlandia", x)] <- "Finland"
  x[grepl("^Jap",      x)]  <- "Japan"
  x[grepl("^Turqu|^Turkey", x)] <- "Türkiye"
  x
}
datos_raw$country <- norm_pais(datos_raw$country)
datos_raw$country[grepl("Platts", datos_raw$author)] <- "Tanzania/Bangladesh/Peru"

datos_raw <- standardise_schema(datos_raw)

datos_campy <- datos_raw %>%
  mutate(Study_ID = make_study_id(author, year, Record_ID),
         Comparison = make_comparison(Index_2, Reference_2),
         market_rank = market_rank_test(Index_1, Index_2),
         eligibility_rank = case_when(eligibility == 1 ~ 1L, eligibility == 2 ~ 2L,
                                      is.na(eligibility) ~ 3L, TRUE ~ 9L),
         sensitivity = if_else((TP+FN) > 0, TP/(TP+FN), NA_real_),
         specificity = if_else((TN+FP) > 0, TN/(TN+FP), NA_real_),
         sample_prev_calc = if_else((TP+FN+FP+TN) > 0, (TP+FN)/(TP+FN+FP+TN), NA_real_))

datos_campy$eligibility_ok <- eligible_filter(datos_campy)

datos_sensibilidad <- datos_campy %>%
  filter(solo_sensibilidad == 1, !is.na(Comparison)) %>%
  add_number_id()

datos_comparaciones <- datos_campy %>%
  filter(solo_sensibilidad == 0,
         eligibility_ok, !is.na(Comparison),
         Reference_2 %in% c("Culture","PCR"),
         Index_2 %in% c("lateral flow immunoassay","microplate assay"),
         !is.na(TP), !is.na(FN), !is.na(FP), !is.na(TN)) %>%
  add_number_id()

global_selection_audit <- datos_comparaciones %>% select_one_per_study() %>%
  arrange(Comparison, Study_ID, selection_order)

LFA_Culture_all <- get_dataset(datos_comparaciones, "LFA_Culture")
LFA_PCR_all     <- get_dataset(datos_comparaciones, "LFA_PCR")
EIA_Culture_all <- get_dataset(datos_comparaciones, "EIA_Culture")
EIA_PCR_all     <- get_dataset(datos_comparaciones, "EIA_PCR")

LFA_Culture <- get_dataset(global_selection_audit, "LFA_Culture", TRUE)
LFA_PCR     <- get_dataset(global_selection_audit, "LFA_PCR", TRUE)
EIA_Culture <- get_dataset(global_selection_audit, "EIA_Culture", TRUE)
EIA_PCR     <- get_dataset(global_selection_audit, "EIA_PCR", TRUE)

datasets_global <- list(LFA_Culture = LFA_Culture, LFA_PCR = LFA_PCR,
                        EIA_Culture = EIA_Culture, EIA_PCR = EIA_PCR)
datasets_all <- list(LFA_Culture_all = LFA_Culture_all, LFA_PCR_all = LFA_PCR_all,
                     EIA_Culture_all = EIA_Culture_all, EIA_PCR_all = EIA_PCR_all)

comparison_summary <- data.frame(
  dataset = names(datasets_global),
  n_rows_global    = vapply(datasets_global, nrow, integer(1)),
  n_studies_global = vapply(datasets_global, function(z) n_distinct(z$Study_ID), integer(1)),
  n_rows_all       = vapply(datasets_all, nrow, integer(1)),
  n_studies_all    = vapply(datasets_all, function(z) n_distinct(z$Study_ID), integer(1)),
  row.names = NULL, stringsAsFactors = FALSE)

saveRDS(list(datos_campy = datos_campy, datos_comparaciones = datos_comparaciones,
             global = datasets_global, all = datasets_all,
             audit = global_selection_audit, sensibilidad = datos_sensibilidad),
        file.path(OUT_DIR, "datasets_dta_campy.rds"))

print(comparison_summary)
cat("\nExcluidos por eligibility==0:\n")
print(datos_campy %>% filter(eligibility == 0) %>%
        select(Study_ID, Index_1, Index_2, Reference_2, TP, FN, FP, TN))
cat("\nPlataformas inmunológicas presentes en la planilla pero fuera de las 4 comparaciones centrales:\n")
print(datos_campy %>% filter(!Index_2 %in% c("lateral flow immunoassay","microplate assay")) %>%
        count(Index_1, Index_2, Reference_2))


.script_dir <- "D:/Gene_analysis/test"
source(file.path(.script_dir, "step_1 define and dictionary.R"))
start_script("step_5 analysis_sets.R")

ax <- load_axis_scores()
meta <- read_csv(RE_GEO_META)
tm <- read_csv(RE_TCVGH_META)
source(TCVGH_OUTCOME_R)  # the 37 TCVGH in-hospital outcomes, listed in that file
tv <- TCVGH_IN_HOSPITAL_MORTALITY
cts <- load_cts_labels()
assert_unique_key(meta, c("dataset", "geo_accession"), "GEO metadata")

manifest <- ax[, c("dataset", "sample_id")]
manifest$patient_id <- paste(manifest$dataset, manifest$sample_id, sep = ":")
manifest$case_status <- "unknown"
manifest$time_from_baseline <- NA_real_
manifest$death <- NA_integer_
manifest$outcome_type <- "unknown"

geo <- function(ds) {
  i <- which(manifest$dataset == ds)
  g <- meta[match(manifest$sample_id[i], meta$geo_accession), , drop = FALSE]
  assert(!anyNA(g$geo_accession), ds, ": every profile needs deposited metadata")
  list(i = i, g = g)
}

x <- geo("GSE134347")
ds <- tolower(trimws(x$g$disease.state))
manifest$case_status[x$i] <- ifelse(ds == "sepsis", "sepsis",
                              ifelse(grepl("control|healthy", ds), "control",
                              ifelse(ds == "noninfectious", "nonseptic", "unknown")))

x <- geo("GSE26440")
manifest$case_status[x$i] <- ifelse(trimws(x$g$disease.state) == "septic shock patient", "septic_shock",
                              ifelse(trimws(x$g$disease.state) == "normal control", "control", "unknown"))
manifest$death[x$i] <- ifelse(trimws(x$g$outcome) == "Nonsurvivor", 1L,
                       ifelse(trimws(x$g$outcome) == "Survivor", 0L, NA_integer_))
manifest$outcome_type[x$i] <- "in_hospital"

x <- geo("GSE57065")
is_patient <- grepl("^Blood_P[0-9]+_H[0-9]+$", trimws(x$g$title))
manifest$case_status[x$i] <- ifelse(is_patient, "sepsis", "control")
pid <- sub("^Blood_(P[0-9]+)_H[0-9]+$", "\\1", trimws(x$g$title))
pid[!is_patient] <- trimws(x$g$title[!is_patient])
manifest$patient_id[x$i] <- paste("GSE57065", pid, sep = ":")
manifest$time_from_baseline[x$i] <-
  suppressWarnings(as.numeric(sub("[^0-9.].*$", "", trimws(x$g$collection.time)))) / 24

x <- geo("GSE65682")
is_healthy <- !is.na(x$g$icu_acquired_infection) & tolower(trimws(x$g$icu_acquired_infection)) == "healthy"
manifest$case_status[x$i] <- ifelse(is_healthy, "control", "sepsis")
manifest$death[x$i] <- suppressWarnings(as.integer(x$g$mortality_event_28days))
manifest$outcome_type[x$i] <- "28_day"

x <- geo("GSE95233")
is_sepsis <- trimws(x$g$time.point) %in% c("D01", "D02", "D03")
manifest$case_status[x$i] <- ifelse(is_sepsis, "sepsis", "control")
manifest$patient_id[x$i] <- paste("GSE95233", sub("_D[0-9]+$", "", trimws(x$g$title)), sep = ":")
manifest$time_from_baseline[x$i] <- match(trimws(x$g$time.point), c("D01", "D02", "D03")) - 1
manifest$death[x$i] <- ifelse(trimws(x$g$survival) == "Non Survivor", 1L,
                       ifelse(trimws(x$g$survival) == "Survivor", 0L, NA_integer_))
manifest$outcome_type[x$i] <- "in_hospital"

i <- which(manifest$dataset == "TCVGH_GSE216902")
g <- tm[match(manifest$sample_id[i], tm$sample_id), , drop = FALSE]
assert(!anyNA(g$sample_id), "every TCVGH profile needs a patient and visit")
manifest$patient_id[i] <- paste("TCVGH", g$patient_number, sep = ":")
manifest$case_status[i] <- "sepsis"
manifest$time_from_baseline[i] <- ifelse(g$visit_code == "Day1", 0, ifelse(g$visit_code == "Day8", 7, NA_real_))
manifest$death[i] <- as.integer(tv$in_hospital_death[match(as.character(g$patient_number),
                                                            as.character(tv$patient_number))])
manifest$outcome_type[i] <- "in_hospital"
assert(setequal(as.character(g$patient_number), as.character(tv$patient_number)),
       "the listed TCVGH outcomes must cover exactly the deposited GSE216902 patients")
record_value("tcvgh.outcome.patients", nrow(tv), "listed in lib/tcvgh_in_hospital_mortality.R")
record_value("tcvgh.outcome.in_hospital_deaths", sum(tv$in_hospital_death))
record_value("tcvgh.geo_condition_disagrees_with_death",
             sum((tv$geo_condition == "Non-Responder") != (tv$in_hospital_death == 1L)),
             "the deposited GEO Responder/Non-Responder label is not mortality")

# controls carry no prognostic endpoint
manifest$death[manifest$case_status %in% c("control", "nonseptic")] <- NA_integer_
manifest$outcome_available <- !is.na(manifest$death)
manifest$has_cts_label <- paste(manifest$dataset, manifest$sample_id) %in% paste(cts$cohort, cts$Sample)

# guards on the two case-status rules that were once wrong
assert(sum(manifest$dataset == "GSE26440" & manifest$case_status == "septic_shock") == 98L, "GSE26440 septic shock = 98")
assert(sum(manifest$dataset == "GSE26440" & manifest$case_status == "control") == 32L, "GSE26440 controls = 32")
assert(sum(manifest$dataset == "GSE95233" & manifest$case_status == "sepsis") == 102L, "GSE95233 sepsis samples = 102")
assert(sum(manifest$dataset == "GSE65682" & manifest$case_status == "control") == 42L, "GSE65682 controls = 42")
assert_unique_key(manifest, c("dataset", "sample_id"), "manifest")
write_csv(manifest, file.path(TABLE_DIR, "sample_manifest.csv"))

# ---- Table 1 --------------------------------------------------------------------------
sampling <- function(d) {
  t <- sort(unique(d$time_from_baseline[is.finite(d$time_from_baseline)]))
  if (length(t) < 2) return("Single")
  day <- t + 1
  if (length(day) == 2) sprintf("Repeated (baseline and %g days later)", t[2] - t[1])
  else sprintf("Repeated (Days %g-%g)", min(day), max(day))
}
mortality_def <- c(GSE65682 = "28-day mortality", GSE26440 = "In-hospital mortality",
                   GSE95233 = "In-hospital mortality", TCVGH_GSE216902 = "In-hospital mortality")
t1 <- do.call(rbind, lapply(COHORTS$dataset, function(ds) {
  d <- manifest[manifest$dataset == ds, ]
  data.frame(dataset = ds, accession = COHORTS$accession[COHORTS$dataset == ds],
             route = COHORTS$route[COHORTS$dataset == ds],
             platform = COHORTS$platform[COHORTS$dataset == ds],
             deposited_samples = nrow(d), sampling = sampling(d),
             mortality_definition = if (ds %in% names(mortality_def)) mortality_def[[ds]] else "N/A",
             represented_patients = length(unique(d$patient_id)),
             repeat_samples = nrow(d) - length(unique(d$patient_id)),
             controls_nonseptic = sum(d$case_status %in% c("control", "nonseptic")),
             stringsAsFactors = FALSE)
}))
write_csv(t1, file.path(TABLE_DIR, "Table1_cohort_inventory.csv"))
print(t1[, c("accession", "route", "platform", "deposited_samples", "sampling", "represented_patients")], row.names = FALSE)

record_value("cohorts.n", nrow(t1))
record_value("profiles.total", sum(t1$deposited_samples))
record_value("profiles.repeat_samples", sum(t1$repeat_samples))
record_value("patients.total", sum(t1$represented_patients))
for (r in seq_len(nrow(t1))) record_value(paste0("table1.deposited.", t1$dataset[r]), t1$deposited_samples[r])
for (rt in c("GEO", "ARCHS4", "recount3")) record_value(paste0("routes.", rt), sum(t1$route == rt))
record_value("table1.microarray_cohorts", sum(grepl("^Microarray", t1$platform)),
             "GSE134347 is GPL17586 (Affymetrix HTA 2.0), a microarray platform")
record_value("table1.is_microarray.GSE134347", as.numeric(grepl("^Microarray", t1$platform[t1$dataset == "GSE134347"])),
             "platform read from the GEO series matrix (!Series_platform_id = GPL17586)")
last_day <- function(ds) max(manifest$time_from_baseline[manifest$dataset == ds], na.rm = TRUE) + 1
record_value("table1.last_sampling_day.GSE57065", last_day("GSE57065"), "day 1 = first deposited sample")
record_value("table1.last_sampling_day.GSE95233", last_day("GSE95233"), "day 1 = first deposited sample")
record_value("table1.days_after_baseline.TCVGH", last_day("TCVGH_GSE216902") - 1,
             "GEO SOFT labels the two visits Day1 and Day8: baseline and 7 days later")

# ---- mortality analysis set: one sample per patient -----------------------------------
d <- merge(ax[, c("dataset", "sample_id", AXIS_CODES)], manifest, by = c("dataset", "sample_id"))
risk <- rbind(
  transform(d[d$dataset == "GSE65682" & d$case_status == "sepsis" & d$outcome_available, ], cohort = "GSE65682"),
  transform(d[d$dataset == "GSE26440" & d$case_status == "septic_shock" & d$outcome_available, ], cohort = "GSE26440"),
  transform(d[d$dataset == "GSE95233" & d$case_status == "sepsis" & d$outcome_available &
                d$time_from_baseline == 0, ], cohort = "GSE95233"),
  transform(d[d$dataset == "TCVGH_GSE216902" & d$outcome_available & d$time_from_baseline == 0, ], cohort = "TCVGH")
)
for (co in MORTALITY_COHORTS) {
  z <- risk[risk$cohort == co, ]
  assert(nrow(z) == length(unique(z$patient_id)), co, ": mortality set must be one sample per patient")
  record_value(paste0("mortality.n.", co), nrow(z))
  record_value(paste0("mortality.deaths.", co), sum(z$death))
}
record_value("mortality.n.total", nrow(risk))
record_value("mortality.deaths.total", sum(risk$death))
assert(nrow(risk) == 665L && sum(risk$death) == 156L, "mortality set must be 665 patients / 156 deaths")
write_csv(risk[, c("cohort", "dataset", "sample_id", "patient_id", "death", AXIS_CODES)],
          file.path(TABLE_DIR, "mortality_analysis_set.csv"))
cat(sprintf("   mortality analysis set: %d patients, %d deaths\n", nrow(risk), sum(risk$death)))


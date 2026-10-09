##補上TCVGH死亡狀態
TCVGH_IN_HOSPITAL_MORTALITY <- utils::read.table(header = TRUE, stringsAsFactors = FALSE, text = "
patient_number  in_hospital_death  geo_condition
13              0                  Responder
36              0                  Responder
41              1                  Non-Responder
43              0                  Non-Responder
44              0                  Non-Responder
93              0                  Responder
94              1                  Responder
97              1                  Non-Responder
100             0                  Non-Responder
106             0                  Responder
111             1                  Non-Responder
112             0                  Responder
115             0                  Responder
116             0                  Responder
117             0                  Non-Responder
122             0                  Responder
128             0                  Responder
130             0                  Responder
142             0                  Responder
148             0                  Responder
164             1                  Non-Responder
165             1                  Non-Responder
166             0                  Responder
168             1                  Non-Responder
169             0                  Responder
171             1                  Non-Responder
180             0                  Responder
183             0                  Responder
205             0                  Responder
208             0                  Responder
215             0                  Responder
220             0                  Responder
224             0                  Non-Responder
225             0                  Responder
230             0                  Non-Responder
237             0                  Non-Responder
241             0                  Responder
")

local({
  x <- TCVGH_IN_HOSPITAL_MORTALITY
  stopifnot(nrow(x) == 37L, !anyDuplicated(x$patient_number),
            all(x$in_hospital_death %in% c(0L, 1L)), sum(x$in_hospital_death) == 8L,
            all(x$geo_condition %in% c("Responder", "Non-Responder")))
})

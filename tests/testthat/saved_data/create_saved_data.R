# Set up a dataset for regular use

library(tidyverse)
load_all()

kd_signature_id <- "LINCSKD_28"
cp_signature_id <- "LINCSCP_5821"
oe_signature_id <- "LINCSOE_104"
inv_signature_id <- "LINCS_INV"

kd_signature <- get_signature(kd_signature_id)
cp_signature <- get_signature(cp_signature_id)
oe_signature <- get_signature(oe_signature_id)
inv_signature <- get_signature(inv_signature_id)

empty_signature <- inv_signature

filtered_at_0_any <-
    filter_signature(kd_signature, direction = "any", threshold = 0.0)

filtered_at_0_5_any <-
    filter_signature(kd_signature, direction = "any", threshold = 0.5)
filtered_at_0_5_up <-
    filter_signature(kd_signature, direction = "up", threshold = 0.5)
filtered_at_0_5_dn <-
    filter_signature(kd_signature, direction = "down", threshold = 0.5)

filtered_at_1_2_any <-
    filter_signature(kd_signature,
        direction = "any",
        threshold = c(-1L, 2L)
    )
filtered_at_1_2_up <-
    filter_signature(kd_signature,
        direction = "up",
        threshold = c(-1L, 2L)
    )
filtered_at_1_2_dn <-
    filter_signature(kd_signature,
        direction = "down",
        threshold = c(-1L, 2L)
    )

filtered_at_0_95_any <-
    filter_signature(kd_signature, direction = "any", prop = 0.05)
filtered_at_0_95_up <-
    filter_signature(kd_signature, direction = "up", prop = 0.05)
filtered_at_0_95_dn <-
    filter_signature(kd_signature, direction = "down", prop = 0.05)

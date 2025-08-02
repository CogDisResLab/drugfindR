library(vcr) # *Required* as vcr is set up on loading

invisible(vcr::vcr_configure(
    dir = vcr::vcr_test_path("cassettes"),
    record = "once",
    preserve_exact_body_bytes = TRUE,
    match_requests_on = c("method", "uri")
))

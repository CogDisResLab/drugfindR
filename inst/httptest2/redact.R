# Custom redacting functions for truncating iLINCS responses

truncateSignature <- function(response) {
    url <- url_parse(response[["url"]])

    if (url[["path"]] == "/api/ilincsR/downloadSignature") {
        print("Signature response")
        response
    } else {
        response
    }
}

truncateConcordants <- function(response) {
    url <- url_parse(response[["url"]])

    if (url[["path"]] == "/api/SignatureMeta/uploadAndAnalyze") {
        within_body_text(response, function(body) {
            body <- jsonlite::fromJSON(body)

            shortTable <- body[["status"]][["concordanceTable"]][1L:10L, ]
            body[["status"]][["concordanceTable"]] <- shortTable

            jsonlite::toJSON(body, auto_unbox = TRUE)
        })
    } else {
        response
    }
}

function(response) {
    require(jsonlite, quietly = TRUE)

    response |>
        redact_headers("Set-Cookie") |>
        truncateSignature() |>
        truncateConcordants() |>
        gsub_response("https\\://www.ilincs.org/api/", "ilincs/")
}

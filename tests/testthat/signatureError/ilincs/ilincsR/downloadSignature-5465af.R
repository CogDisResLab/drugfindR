structure(list(
    method = "GET", url = "ilincs/ilincsR/downloadSignature?sigID=&noOfTopGenes=Inf",
    status_code = 400L, headers = structure(list(
        Date = "Wed, 12 Nov 2025 20:57:55 GMT",
        Server = "Apache/2.4.65 (Unix) OpenSSL/3.5.1", `Strict-Transport-Security` = "max-age=63072000; includeSubDomains",
        Vary = "Origin,Accept-Encoding", `Access-Control-Allow-Credentials` = "true",
        `X-XSS-Protection` = "1; mode=block", `X-Download-Options` = "noopen",
        `X-Content-Type-Options` = "nosniff", `Content-Type` = "application/json; charset=utf-8",
        `X-Frame-Options` = "SAMEORIGIN", `Set-Cookie` = "REDACTED",
        Connection = "close", `Transfer-Encoding` = "chunked"
    ), class = "httr2_headers"),
    body = charToRaw("{\"error\":{\"statusCode\":400,\"name\":\"Error in iLincsR downloadSignature function. No genes or wrong data set picked.\",\"message\":\"argument \\\"sigID\\\" is missing, with no default\\n\\nIn call:\\ndownloadSignature(noOfGenes = c(\\\"Inf\\\"))\\n\",\"code\":\"ilincsR_ERROR\"}}"),
    timing = c(
        redirect = 0, namelookup = 0.000178, connect = 0.015839,
        pretransfer = 0.027281, starttransfer = 0.178171, total = 0.178199
    ), cache = new.env(parent = emptyenv())
), class = "httr2_response")

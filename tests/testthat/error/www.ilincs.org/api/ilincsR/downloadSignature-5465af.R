structure(list(
    method = "GET", url = "https://www.ilincs.org/api/ilincsR/downloadSignature?sigID=&noOfTopGenes=Inf",
    status_code = 400L, headers = structure(list(
        Date = "Tue, 26 Aug 2025 01:54:24 GMT",
        Server = "Apache/2.4.63 (Unix) OpenSSL/3.0.15", `Strict-Transport-Security` = "max-age=63072000; includeSubDomains",
        Vary = "Origin,Accept-Encoding", `Access-Control-Allow-Credentials` = "true",
        `X-XSS-Protection` = "1; mode=block", `X-Download-Options` = "noopen",
        `X-Content-Type-Options` = "nosniff", `Content-Type` = "application/json; charset=utf-8",
        `X-Frame-Options` = "SAMEORIGIN", `Set-Cookie` = "REDACTED",
        Connection = "close", `Transfer-Encoding` = "chunked"
    ), class = "httr2_headers"),
    body = charToRaw("{\"error\":{\"statusCode\":400,\"name\":\"Error in iLincsR downloadSignature function. No genes or wrong data set picked.\",\"message\":\"argument \\\"sigID\\\" is missing, with no default\\n\\nIn call:\\ndownloadSignature(noOfGenes = c(\\\"Inf\\\"))\\n\",\"code\":\"ilincsR_ERROR\"}}"),
    timing = c(
        redirect = 0, namelookup = 0.001464, connect = 0.037033,
        pretransfer = 0.128914, starttransfer = 0.314506, total = 0.314527
    ), cache = new.env(parent = emptyenv())
), class = "httr2_response")

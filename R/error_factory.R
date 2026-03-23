## Internal utility: wrap a function call and capture message/warning/error
## conditions as attributes on the returned value.
factory <- function (fun, debug=FALSE,
                     errval="An error occurred in the factory function",
                     types=c("message","warning","error")) {
    function(...) {
    errorOccurred <- FALSE
    warn <- err <- msg <- NULL
    res <- withCallingHandlers(tryCatch(fun(...),
                                        error = function(e) {
        if (debug) cat("error: ",conditionMessage(e),"\n")
        err <<- conditionMessage(e)
        errorOccurred <<- TRUE
        NULL
    }), warning = function(w) {
        if (!"warning" %in% types) {
            warning(conditionMessage(w))
        } else {
            warn <<- append(warn, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    },
    message = function(m) {
        if (debug) cat("message: ",conditionMessage(m),"\n")
        if (!"message" %in% types) {
            message(conditionMessage(m))
        } else {
            msg <<- append(msg, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    })
    if (errorOccurred) {
        if (!"error" %in% types) stop(err)
        res <- errval
    }

    setattr <- function(x, attrib, value) {
        attr(x,attrib) <- value
        x
    }

    attr_fun <- function(x,str,msg) {
        setattr(x,paste0("factory-",str), if(is.character(msg)) msg else NULL)
    }

    res <- attr_fun(res, "message", msg)
    res <- attr_fun(res, "warning", warn)
    res <- attr_fun(res, "error", err)

    return(res)
  }
}


    

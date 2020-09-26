#' @title check if package is available
#' @param package Package name
#' @noRd
check_package <- function(package){
    if (!requireNamespace(package, quietly = TRUE)) {
        stop(package, " package is needed for this function to work. Please install it.",
             call. = FALSE)
    }
}

#' @title register cores
#' @param cores A integer which defines the number of cores to be used in parallel
#' @noRd
register_cores <- function(cores){
    
    check_package("parallel")
    check_package("doParallel")
    
    parallel <- FALSE
    if (cores > 1){
        if (cores > parallel::detectCores()) cores <- parallel::detectCores()
        doParallel::registerDoParallel(cores)
        parallel = TRUE
    }
    return(parallel)
}


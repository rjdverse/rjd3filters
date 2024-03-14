#' @rdname moving_average
#' @export
setMethod(f = "show",
          signature = "moving_average",
          definition = function(object){
            print(.jcall(.ma2jd(object), "S", "toString"))
            invisible(object)
          })
#' @rdname finite_filters
#' @export
setMethod(f = "show",
          signature = c("finite_filters"),
          definition = function(object){
            x <- as.matrix(object)
            print(x)
          })

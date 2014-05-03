# Exemplo S3:
Person <- function(name, weight, height){
  x <- list(name = name, weight = weight, height = height)
  class(x) <- "Person"
  x
}

bmi <- function(x) UseMethod('bmi', x)
bmi.Person <- function(x) {
  x$weight/x$height^2
}

print.Person <- function(x) {
  cat("   Name: ", x$name, "\n",
      "Weight: ", x$weight, "Kg\n",
      "Height: ", x$height, "m\n",
      "   BMI: ", round(bmi(x), 2), "\n")
}


bob <- Person("Bob", 70, 1.76)
bob


# Exemplo S4:
setClass("PersonS4", 
         representation(name = 'character', 
                        weight = 'numeric', 
                        height = 'numeric'),
         prototype(name = NA_character_, 
                   weight = NA_real_, 
                   height = NA_real_),
         validity = function(object) {
           if(object@weight > object@height)
             TRUE
           else
             "Weight cannot be smaller than height!"
         })

setGeneric("bmiS4", function(object) standardGeneric("bmiS4"))
setMethod("bmiS4", signature("PersonS4"), 
          function(object) {
            object@weight/object@height^2
          })

setMethod("show", signature("PersonS4"),
          function(object) {
            cat("   Name: ", object@name, "\n",
                "Weight: ", object@weight, "Kg\n",
                "Height: ", object@height, "m\n",
                "   BMI: ", round(bmiS4(object), 2), "\n")
          })

bobS4 <- new("PersonS4", name = "Bob", weight = 70, height = 1.76)
bobS4

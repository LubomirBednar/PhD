# learn user defined functions
# The basic structure of a function
# Create the function my.mean()
my.mean <- function(x) {   # Single input called x
  
  output <- sum(x) / length(x) # Calculate output
  
  return(output)  # Return output to the user after running the function
  
}

my.mean(iris$Sepal.Length)
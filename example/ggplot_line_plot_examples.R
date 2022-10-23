# Import the tip library
library(tip)

# Create the variable that appears on the horizontal axis
x <- 1:10

# Create the variable that appears on the vertical axis
y <- rnorm(n = length(x), mean = 3, sd = 1)

# Create a label that appears on the horizontal axis
xlab <- "x"

# Create a label that appears on the vertical axis
ylab <- "y"

# Create the plot of y versus x with
ggplot_line_point(.x = x, .y = y, .xlab = xlab, .ylab = ylab)

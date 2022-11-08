# Import the tip library
library(tip)

# Choose an arbitrary random seed to generate the data
set.seed(4*8*15*16*23*42)

# Generate a symmetric posterior probability matrix
# Each element is the probability that the two subjects belong
# to the same cluster
n1 <- 10
posterior_prob_matrix <- matrix(NA, nrow = n1, ncol = n1)
for(i in 1:n1){
  for(j in i:n1){
    if(i != j){
      posterior_prob_matrix[i,j] <- runif(n=1,min=0,max=1)
      posterior_prob_matrix[j,i] <- posterior_prob_matrix[i,j]
    }else{
      posterior_prob_matrix[i,j] <- 1.0
    }
  }
}

# --- BEGIN GRAPH PLOT 1: NO COLORS OR NODE LABELS ---

# Set an arbitrary random seed
random_seed <- 815

# Set add_node_labels to FALSE
add_node_labels <- FALSE

# Set the node size
node_size <- 6

# Construct the graph plot
ggnet2_network_plot(.matrix_graph = posterior_prob_matrix,
                    .subject_names = vector(),
                    .subject_class_names = vector(),
                    .class_colors = vector(),
                    .class_shapes = vector(),
                    .random_seed = random_seed,
                    .node_size = node_size,
                    .add_node_labels = add_node_labels)

# --- END GRAPH PLOT 1: NO COLORS OR NODE LABELS ---

# --- BEGIN GRAPH PLOT 2: NO COLORS, BUT ADD NODE LABELS ---

# Set an arbitrary random seed
random_seed <- 815

# Add node labels to the plot
add_node_labels <- TRUE

# Set graph nodes (i.e. subject identifier) a larger size
node_size <- 6

# Add subject names to the plot
subject_names <- paste("S", 1:n1, sep = "")

# Construct the graph plot
ggnet2_network_plot(.matrix_graph = posterior_prob_matrix,
                    .subject_names = subject_names,
                    .subject_class_names = vector(),
                    .class_colors = vector(),
                    .class_shapes = vector(),
                    .random_seed = random_seed,
                    .node_size = 6,
                    .add_node_labels = TRUE)

# --- END GRAPH PLOT 2: NO COLORS, BUT ADD NODE LABELS ---


# --- BEGIN GRAPH PLOT 3: ADD COLORS AND NODE LABELS ---

# Set an arbitrary random seed
random_seed <- 815

# Add node labels to the plot
add_node_labels <- TRUE

# Set graph nodes (i.e. subject identifier) a larger size
node_size <- 10

# Add subject names to the plot
subject_names <- paste("S", 1:n1, sep = "")

# Create a vector of class labels
subject_class_names <- c("Class2","Class2","Class1","Class2","Class1",
                              "Class1","Class2","Class1","Class1","Class2")

# Assign a color to each class; this can be a character value or a hex value
class_colors <- c("Class1" = "skyblue", "Class2" = "#FF9999")

# Assign a pch integer value to each class
class_shapes <- c("Class1" = 16, "Class2" = 17)

# Generate the plot
ggnet2_network_plot(.matrix_graph = posterior_prob_matrix,
                    .subject_names = subject_names,
                    .subject_class_names = subject_class_names,
                    .class_colors = class_colors,
                    .class_shapes = class_shapes,
                    .random_seed = random_seed,
                    .node_size = node_size,
                    .add_node_labels = add_node_labels)

# --- END GRAPH PLOT 3: ADD COLORS AND NODE LABELS ---

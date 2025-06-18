# Code
rm(list = ls())
cat("\14")
while (dev.cur()>1) dev.off()

options(digits = 7)  # Default is usually 7

# Open required package libraries
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if(!require(data.table)) install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org")
if(!require(stringr)) install.packages("stringr", repos = "http://cran.us.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if(!require(knitr)) install.packages("knitr", repos = "http://cran.us.r-project.org")
if(!require(kableExtra)) install.packages("kableExtra", repos = "http://cran.us.r-project.org")
if(!require(scales)) install.packages("scales", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(readr)) install.packages("readr", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(rstudioapi)) install.packages("rstudioapi", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(funModeling)) install.packages("funModeling", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(gridExtra)) install.packages("gridExtra", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(doParallel)) install.packages("doParallel", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(nortest)) install.packages("nortest", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(broom)) install.packages("broom", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(rstatix)) install.packages("rstatix", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(corrplot)) install.packages("corrplot", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")
if(!require(pROC)) install.packages("corrplot", dependencies = c("Depends", "Suggests"), repos = "http://cran.us.r-project.org")

cat("downloading data from repository...\n ")

# Mark the script directory as R working directory
script_path <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(script_path))


#-------------------------load/download data set--------------------------
# Load data set if already saved in HD, if not, download and save
file_path <- file.path(getwd(), "my_data_BreastCancer_final.RData")

if (file.exists(file_path)) {
  load(file_path)
  message("Data loaded from: ", file_path)
} else{
  # downloads and process data from repository at University of California Irvine
name_data <- "breast_cancer_data"
download.file(
  url = "https://archive.ics.uci.edu/static/public/17/breast+cancer+wisconsin+diagnostic.zip",
  destfile = str_c(name_data,".zip"),  # Save with simpler name
  mode = "wb"
)

# Extract data to a sub-folder
dir.create(name_data, showWarnings = FALSE)
unzip(str_c(name_data,".zip"), exdir = name_data)
cat("data extracted to folder breast_cancer_data ...\n ")

# 6. Number of attributes: 32 (ID, diagnosis, 30 real-valued input features)
# 7. Attribute information
# 
# 1) ID number
# 2) Diagnosis (M = malignant, B = benign)
# 3-32)
# 
# Ten real-valued features are computed for each cell nucleus:
#   
# a) radius (mean of distances from center to points on the perimeter)
# b) texture (standard deviation of gray-scale values)
# c) perimeter
# d) area
# e) smoothness (local variation in radius lengths)
# f) compactness (perimeter^2 / area - 1.0)
# g) concavity (severity of concave portions of the contour)
# h) concave points (number of concave portions of the contour)
# i) symmetry 
# j) fractal dimension ("coastline approximation" - 1)

basic_features <- c(
  "radius",
  "texture",
  "perimeter",
  "area",
  "smoothness",
  "compactness",
  "concavity",
  "concave_points",
  "symmetry",
  "fractal_dim"
)

# Read data from unzipped cvs file "wdbc.data"
cancer_data <- read_csv(
  file = file.path(name_data, "wdbc.data"), 
  col_names = FALSE,  # No headers in original file
  show_col_types = FALSE # echo off column types
)

# Add column names (from wdbc.names). Specific feature names not used, just generic mean, 
# standard error and worst. The data is arranged in 10-column data blocks
colnames(cancer_data) <- c("iD", "diagnosis", 
                           paste0("avg_", basic_features),
                           paste0("se_", basic_features),
                           paste0("w_", basic_features))

# Convert Diagnosis to factor. More convenient for binary categorical data
cancer_data$diagnosis <- factor(cancer_data$diagnosis, 
                                levels = c("B", "M"),
                                labels = c("Benign", "Malignant"))


# Split into train and test sets (80:20)
set.seed(200, sample.kind = "Rounding")

indexes <- createDataPartition(cancer_data$diagnosis, times = 1, p = 0.2, list = FALSE)
test_set <- cancer_data[indexes,]
train_set <- cancer_data[-indexes,]

# save dataset to HD for other tests
save(cancer_data, train_set, test_set, file = "my_data_BreastCancer_final.RData")
cat("data successfully saved...\n ")
}
#---------------------------------------------------------------------
#-------------------------data exploration--------------------------


F_EXPLORATION <- FALSE
if(F_EXPLORATION){
  
  na_inf_summary <-funModeling::df_status(cancer_data, print_results = FALSE)
  
  kable(na_inf_summary, caption = "Summary of missing data") %>%
    kable_styling(latex_options = c("scale_down","HOLD_position"))
  
  prop_diag <- cancer_data %>%
    count(diagnosis) %>%  # counts occurrences of each diagnosis
    mutate(proportion = n / sum(n)) # calculates the proportion relative to the total
  
  # Plot proportions of "Malignant" and "Beningn" in the data set
  ggplot(prop_diag, aes(x = diagnosis, y = proportion, fill = diagnosis)) +
    geom_bar(stat = "identity") + # Use stat="identity" because 'y' now directly represents height
    geom_text(
      aes(label = scales::percent(proportion)), # Format as percentage
      vjust = -0.5, # Adjust vertical position slightly above the bar
      color = "black"
    ) +
    scale_fill_manual(values = c("Benign" = "steelblue", "Malignant" = "firebrick")) + 
    labs(
      x = "Diagnosis",
      y = "Proportion"
    ) +
    theme_light()
  
  # {r table-proportion-diag, echo=FALSE}
  prop_diag <- cancer_data %>%
    count(diagnosis) %>%  # counts occurrences of each diagnosis
    mutate(proportion = n / sum(n)) # calculates the proportion relative to the total
  
  kable(prop_diag, caption = "Summary of diagnosis proportions in data set") %>%
    kable_styling(latex_options = c("scale_down","HOLD_position"))
  
  
  # {r normality-test, include=FALSE}
  normality_results_AD <- cancer_data %>%
    select(-iD,-diagnosis) %>%
    map_df(~ broom::tidy(ad.test(.)), .id = "feature") %>% # ad.test() works on vectors, so a map must be used
    mutate(normal = ifelse(p.value > 0.05, "Normal", "Non-Normal")) %>%
    select(-method)
  
  # Anderson-Darling test for normality by features grouped by diagnosis
  normality_by_diagnosis <- cancer_data %>% 
    select(-iD) %>%
    group_by(diagnosis) %>% 
    dplyr::summarize(
      across(
        .cols = where(is.numeric),  # Auto-selects all numeric columns (cols 3-32)
        .fns = ~ broom::tidy(ad.test(.)),  # Anderson-Darling test
        .names = "{.col}"  # Column name format
      )
    ) %>%
    tidyr::pivot_longer(
      cols = -diagnosis,
      names_to = "feature",
      values_to = "test_results"
    ) %>%
    unnest(test_results)  # Unpack the tidy() output
  # View results
  normality_by_diagnosis <- normality_by_diagnosis %>%
    select(diagnosis, feature, statistic, p.value) %>%
    mutate(normal = ifelse(p.value > 0.05, "Normal", "Non-Normal")) %>%
    arrange(feature, diagnosis)
  
  # ```{r table-normality, echo=FALSE}
  # formatting for better table printing
  normality_results_AD <- normality_results_AD %>%
    mutate(
      p.value = sprintf("%.3e", p.value), # For 3 decimal places in scientific notation
      statistic = sprintf("%.4e", statistic) # For 4 decimal places
    )
  
  kable(normality_results_AD, caption = "Arlington-Darling normality test") %>%
    kable_styling(latex_options = c("scale_down","HOLD_position"))
  
  # {r plot-normality, fig.cap="Normality Status by Diagnosis Group" }
  # Normality Status by Diagnosis Group
  normality_by_diagnosis %>%
    ggplot(aes(x = diagnosis, y = feature, fill = normal)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(values = c("Normal" = "blue", "Non-Normal" = "firebrick")) +
    labs(
      x = "Diagnosis",
      y = "Feature",
      fill = "Normality"
    ) +
    theme_light()
  
  # ```{r wilcoxon-test, include=FALSE}
  # Wilcoxon test for all features (non-normal data)
  wilcoxon_results <- cancer_data %>%
    select(-iD) %>%
    pivot_longer(-diagnosis, names_to = "feature") %>%
    group_by(feature) %>%
    rstatix::wilcox_test(value ~ diagnosis) %>%  # From rstatix
    adjust_pvalue(method = "fdr") %>%  # Benjamini-Hochberg correction
    add_significance("p.adj")# %>%   arrange(p)
  
  # {r effect-size-rank, include=FALSE}
  effect_sizes_rank <- cancer_data %>%
    select(-iD) %>%
    pivot_longer(-diagnosis, names_to = "feature") %>%
    group_by(feature) %>%
    rstatix::wilcox_effsize(value ~ diagnosis) #%>%  

  effect_sizes_rank %>% arrange(desc(effsize))%>% print(n=30)
    # Merge results
  feature_importance <- wilcoxon_results %>%
    left_join(effect_sizes_rank, by = "feature") %>% arrange(desc(effsize))
  
  
  # {r feature-importance-tab, echo=FALSE}
  col_a<-which(feature_importance$effsize < 0.3 & feature_importance$p.adj<0.05) 
  col_b<-which(feature_importance$effsize < 0.3 & feature_importance$p.adj>0.05)
  
  feat_imp <- feature_importance %>% 
    select(feature,group1.x,group2.x,p,p.adj,effsize,magnitude) %>%
    rename("group 1" = group1.x, "group 2"= group2.x) %>%
    # formatting for better table printing
    mutate(
      p.adj = sprintf("%.4e", p.adj), 
      p = sprintf("%.4e", p) # For 4 decimal places
    )%>%
    arrange(desc(effsize))
  
  # Include data in table
  kable(feat_imp, caption = "Statistical Significance and size effect for the breast cancer dataset") %>%
    row_spec(
      row = col_a, # Add color to small effect rows
      color = "red") %>%
    row_spec(
      row = col_b, # Add color to small effect rows
      color = "blue",
    ) %>%
    kable_styling(latex_options = c("scale_down","HOLD_position"))
  
  # {r feature-importance-plot, fig.cap="Feature Significance Analysis", fig.pos="H"}
  feature_importance %>%
    ggplot(aes(x = effsize, y = -log10(p.adj), color = magnitude)) +
    geom_point() +
    geom_vline(xintercept = c(0.3, 0.5), linetype = "dashed") + # vertical lines for effects 0.10 - < 0.3 (small effect), 0.30 - < 0.5 (moderate effect) and >= 0.5 (large effect)
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") + # horiz line for adj.p-value
    ggrepel::geom_text_repel(
      aes(label = feature), size = 3,
      max.overlaps = Inf,
      box.padding = unit(0.5, "lines"), # Increase padding around labels
      point.padding = unit(0.5, "lines") # Increase padding around points
    ) +
    scale_color_manual(
      values = c(
        "small" = "#505050",
        "moderate" = "blue",
        "large" = "firebrick"  # Red for high importance
      ),
      breaks = c("small","moderate",  "large")  # Force legend order
    ) +
    labs(x = "Effect Size (rank-biserial correlation)", 
         y = "-log10(Adjusted p-value)") +
    theme_light()
  
  
  # {r correlation-matrix, include=FALSE}
  # Compute correlations for numeric features
  cor_matrix <- cor(cancer_data[, 3:32])
  
  
  # {r correlation-plot, fig.cap="Correlation between features"}
  ## Data dimensionality and Principal component analysis
  corrplot::corrplot(cor_matrix,
                     method = "square",      # Shape of the correlation indicators: "circle", "square", "ellipse", "number", "color", "pie"
                     type = "upper",         # Display only the "upper" triangle of the matrix to avoid redundancy
                     order = "hclust",       # Order variables by hierarchical clustering to group similar correlations
                     diag = FALSE,           # Do not display the diagonal (correlation of variable with itself is always 1)
                     tl.col = "black",       # Color of the text labels (variable names)
                     # tl.srt = 45,            # Rotate text labels for better readability if names are long
                     tl.cex = 0.6,           # Size of text labels
                     cl.pos = "r",           # Position of the color legend ("r" for right, "b" for bottom)
                     cl.cex = 0.7,           # Size of the color legend labels
  )
  
  # {r PCA-analisys, include=FALSE}
  # calculate the principal components
  pca_results <- prcomp(cancer_data[, 3:32], scale. = TRUE)
  
  pca_summary <- summary(pca_results)
  var_explained <- pca_summary$importance["Cumulative Proportion", ]
  
  data_plot_pca <- tibble(
    PC = 1:length(var_explained),
    Cumulative_Variance = as.numeric(var_explained)
  )
  
  
  # {r PCA-var-explained-plot, fig.cap="Cumulative Variance Explained by Principal Components"}
  # plot variance explained 
  ggplot(data_plot_pca, aes(x = PC, y = Cumulative_Variance)) +
    geom_line(color = "steelblue", size = 1.2) + # Line connecting the points
    geom_point(color = "steelblue", size = 3, shape = 21, fill = "lightblue") + # Points for each PC
    # horizontal line to indicate a target variance explained (90%)
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", alpha = 0.7) +
    annotate("text", x = max(data_plot_pca$PC) * 0.8, y = 0.95,
             label = "95% Variance Explained", vjust = -0.5, color = "red", size = 4) +
    labs(
      x = "Number of Principal Components",
      y = "Cumulative Proportion of Variance Explained"
    ) +
    # Format y-axis as percentages
    scale_y_continuous(labels = scales::percent, limits = c(data_plot_pca$Cumulative_Variance[1], 1)) +
    # x-axis breaks are at integer PC numbers
    scale_x_continuous(breaks = seq(1, length(var_explained), by = 1)) +
    theme_light()
  
  
  # {r loadings-pca-plot, fig.cap="Loadings of Original Variables on PC1-PC6"}
  X_pcs_to_consider <- 6
  # Access the loadings matrix
  loadings_matrix <- pca_results$rotation
  # Select the loadings for the first X principal components
  top_x_pcs_loadings <- loadings_matrix[, 1:X_pcs_to_consider]
  # 
  # Convert the loadings matrix to a "long" format for ggplot2
  # This creates columns for 'Variable', 'Principal_Component', and 'Loading'
  loadings_long_df <- top_x_pcs_loadings %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Variable") %>% # Convert row names (variable names) to a column
    pivot_longer(
      cols = starts_with("PC"), # Select all columns that start with "PC"
      names_to = "Principal_Component",
      values_to = "Loading"
    )
  # Create the Heatmap Plot
  ggplot(loadings_long_df, aes(x = Principal_Component, y = Variable, fill = Loading)) +
    geom_tile(color = "white", linewidth = 0.5) + # Create tiles with white borders
    scale_fill_gradient2(
      low = "blue",      # Color for negative loadings
      mid = "white",     # Color for loadings near zero
      high = "red",      # Color for positive loadings
      midpoint = 0,      # Center the color scale at zero
      #limit = max(abs(loadings_long_df$Loading)), # Ensure symmetric limits
      name = "Loading Value" # Legend title
    ) +
    labs(
      x = "Principal Component",
      y = "Original Variable"
    ) +
    theme_minimal() 
  
  # {r PC1-PC2-plot, fig.cap="PC1 vs PC2"}
  pca_df <- as.data.frame(pca_results$x)
  pca_df$diagnosis <- cancer_data$diagnosis
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = diagnosis)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Benign" = "lightblue", "Malignant" = "firebrick")) +
    theme_minimal() +
    stat_ellipse(type="norm", lwd = 1.)+
    theme_light()
  
}

################################################################################
################################################################################
################################################################################
# TRAINING AND CROSS VALIDATION
################################################################################
################################################################################
################################################################################

model_metrics_results <- function(m_method, tune_grid, params_tuned){
  # Inputs:
  # m_method: the tained method to be used ("knn", "rf", "svm", "xgb")
  # tune_grid: grid for optimizing each method-depending parameter (k, mtry, eta, sigma, etc...)
  # params_tuned: names of the optimized parameters (only for informatible, summarizing purposes)
  # Outputs:
  # 1: tibble summarizing metrics, 2: tuned models (ROC, Sens, Spec), 3: tibble ROC data, 4: predicted probabilities on test data (pred_probs)
  # some of them will be useful during ensemble modelling process  
  # List to store trained model objects
  trained_models_list <- list()
  
  # List to store ROC curves for plotting
  roc_curves_list <- list()
  
  # Create an empty tibble to store results
  metrics_model_results <- tibble(
    method = character(), # trained method
    metric = character(), # the metric used for selecting the "best" in CV
    Accuracy = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    AUC = numeric(), 
    F1 = numeric(),
    best_p_cv = list(), # store the optimal cross validation value(s)
    parameter = character() # 
  )
  
  # Metrics to optimize for during cross-validation
  # "ROC" is AUC, "Sens" is Sensitivity, "Spec" is Specificity
  optimization_metrics <- c("ROC", "Sens", "Spec")
  # iterate a for loop for training/testing a model by modifying the metrics for the cross validation/folding  
  for (i in 1:length(optimization_metrics)){
    
    opt_metric <- optimization_metrics[i]
    message(paste0("Training ", m_method, " model optimizing for: ", opt_metric))
    
    # Train the  model
    # The `metric` argument in `train` determines which metric `selectionFunction` uses
    # to pick the best tuning parameter (here, 'k').
    # Set an specific seed, for repeatibility results, can be neglected later
    set.seed(7, sample.kind = "Rounding")
    
    model_tuned <- train(
      xS, y,
      method = m_method,
      tuneGrid = tune_grid,
      trControl = fitControl,
      metric = opt_metric
    )
    
    # Store the trained model
    trained_models_list[[opt_metric]] <- model_tuned
    
    # Make predictions on the *test* set (xS_T), will be used in ensemble modelling
    pred_classes <- predict(model_tuned, newdata = xS_T, type = "raw")
    pred_probs <- predict(model_tuned, newdata = xS_T, type = "prob")
    
    # Calculate Confusion Matrix and metrics on the TEST SET
    results_cm <- confusionMatrix(
      pred_classes,
      y_T,
      positive = positive_class_label
    )
    
    # Calculate AUC for the test set
    roc_obj <- roc(
      response = y_T,
      predictor = pred_probs[[positive_class_label]], # probabilities for the positive class
      levels = levels(y_T),
    )
    # Store the best model AUC
    auc_value <- as.numeric(auc(roc_obj))
    roc_curves_list[[opt_metric]] <- roc_obj
    # Store the best model parameter (k_best, sigma_best, etc)
    best_param_value <- model_tuned$bestTune[params_tuned]
    
    # Store the results
    current_results <- tibble(
      method = m_method,
      metric = opt_metric,
      Accuracy = results_cm$overall["Accuracy"],
      Sensitivity = results_cm$byClass["Sensitivity"],
      Specificity = results_cm$byClass["Specificity"],
      AUC = auc_value,
      F1 = results_cm$byClass["F1"],
      best_p_cv = list(best_param_value), 
      parameter = paste(params_tuned, collapse = ", ")
    )
    
    metrics_model_results <- bind_rows(metrics_model_results, current_results)
  }
  
  # Return a 4 elements list
  metrics_model_results <- list(metrics_model_results,
                                trained_models_list,
                                roc_curves_list,
                                pred_probs)
} # end of function model_metrics_results()

############################################################################

cl <- makeCluster(detectCores() -1) # Create a cluster
registerDoParallel(cl)         # Register it with foreach/caret

# Train dataset
x <- train_set %>%
  select(-diagnosis,-iD)

y <- train_set$diagnosis

# Test dataset
x_T <- test_set %>%
  select(-diagnosis,-iD)

y_T <- test_set$diagnosis

# Preprocessing datasets (removing mean and scale by sd of training set)
preProc_params <- preProcess(x, method = c("center", "scale"))
xS <- predict(preProc_params, x)
xS_T <- predict(preProc_params, x_T)

fitControl <- trainControl(method = "repeatedcv",
                           number = 10, # 10-fold cross-validation
                           repeats = 5, # repeat each cross-validation 10 times
                           classProbs = TRUE, # class probabilities computed
                           returnResamp = "final", # only save the final resampled summary metrics (saves out-of-fold predictions for ensemble)
                           savePredictions = "final",
                           allowParallel = TRUE,    # allow parallel processing
                           summaryFunction = twoClassSummary, # For metrics like AUC, Sensitivity, Specificity
                           selectionFunction = "best"          # Select model with max AUC
)
# Assign the positive class label for analysis
positive_class_label <- "Malignant"


# set the grids for training each set of models
all_models_grid <-list(
  data.frame(k = seq(1, 60, 2)), # for KNN
  expand.grid(sigma = c(0.01, 0.1,0.25), C = seq(1,10,1)),  # for SVM
  data.frame(mtry = seq(1, 30, 2)), # for RF
  xgb_grid <- expand.grid(
    nrounds = c(100), #100
    max_depth = seq(2,10,2), # c(3, 6,9),
    eta = seq(0.1,1,0.25), #c(0.01, 0.1,0.5),
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8 )
)

# Training algorithms to be explored
all_models_method <- c("knn","svmRadial","rf","xgbTree")
# Map model methods to their primary tuning parameter names for display
param_names_map <- list(
  knn = "k",
  svmRadial = c("sigma","C"),
  rf = "mtry",
  xgbTree = c("nrounds", "max_depth", "eta")
)

try_method <- 1 # knn
# returns a 4 elements list: 
model_knn <- model_metrics_results(
  all_models_method[try_method], 
  all_models_grid[[try_method]], 
  param_names_map[[try_method]]
)

# update results in the overall tibble
all_model_results <- model_knn[[1]]


try_method <- 2 # svmRadial
model_svm <- model_metrics_results(
  all_models_method[try_method], 
  all_models_grid[[try_method]], 
  param_names_map[[try_method]]
)
all_model_results <- bind_rows(all_model_results, model_svm[[1]])

try_method <- 3 #"rf"
model_rf <- model_metrics_results(
  all_models_method[try_method],
  all_models_grid[[try_method]],
  param_names_map[[try_method]]
)

all_model_results <- bind_rows(all_model_results, model_rf[[1]])

try_method <- 4 # "xgbTree"
model_xgbTree <- model_metrics_results(
  all_models_method[try_method], 
  all_models_grid[[try_method]], 
  param_names_map[[try_method]]
)


all_model_results <- bind_rows(all_model_results, model_xgbTree[[1]])%>% arrange(desc(AUC))
kable(all_model_results)
stopCluster(cl)

#############################################################################
message("Training/testing PCA based models....")
#############################################################################

cl <- makeCluster(detectCores() -1) # Create a cluster
registerDoParallel(cl)         # Register it with foreach/caret
# Preprocessing datasets (removing mean and scale by sd of training set)

preProc_params <- preProcess(x, method = c("center", "scale","pca"), thresh = 0.95)
xS <- predict(preProc_params, x)
xS_T <- predict(preProc_params, x_T)

try_method <- 1 # knn
# function(m_method, tune_grid, params_tuned)

model_knn_PCA <- model_metrics_results(
  all_models_method[try_method], 
  all_models_grid[[try_method]], 
  param_names_map[[try_method]]
)

all_model_results_PCA <- model_knn_PCA[[1]]

try_method <- 2 # svmRadial
model_svm_PCA <- model_metrics_results(
  all_models_method[try_method], 
  all_models_grid[[try_method]], 
  param_names_map[[try_method]]
)
all_model_results_PCA <- bind_rows(all_model_results_PCA, model_svm_PCA[[1]])

try_method <- 3 #"rf"
model_rf_PCA <- model_metrics_results(
  all_models_method[try_method], 
  all_models_grid[[try_method]], 
  param_names_map[[try_method]]
)

all_model_results_PCA <- bind_rows(all_model_results_PCA, model_rf_PCA[[1]])

try_method <- 4 # "xgbTree"
model_xgbTree_PCA <- model_metrics_results(
  all_models_method[try_method], 
  all_models_grid[[try_method]], 
  param_names_map[[try_method]]
)

all_model_results_PCA <- bind_rows(all_model_results_PCA, model_xgbTree_PCA[[1]])%>% arrange(desc(AUC)) %>% mutate(method = paste0(method, "_PCA"))

kable(all_model_results_PCA)


###########################################################################
save(model_knn_PCA,model_svm_PCA,model_rf_PCA,model_xgbTree_PCA,
     model_knn,model_svm,model_xgbTree,model_rf,all_model_results_PCA,
     all_model_results,
     train_set,test_set, file="trained_models.RData")
############################################################################

################################################################################
################################################################################
################################################################################
# WEIGHTED PROBABILISTIC ENSEMBLE
################################################################################
################################################################################
################################################################################
compute_f12 <- function(w,p,yData) {
  w <- as.numeric(w)
  # ensemble_p <- w[1]*p1 + w[2]*p2 + w[n]*pn +...
  ensemble_p <- as.numeric(p%*%as.matrix(w))
  pred_class <- ifelse(ensemble_p >= 0.5, "Malignant", "Benign")
  pred_class <- factor(pred_class, levels = levels(yData))
  cm <- confusionMatrix(pred_class, yData, positive = "Malignant")
  roc_obj <- roc(yData, ensemble_p, levels = rev(levels(yData)))
  
  tibble(F1 = cm$byClass["F1"], AUC = as.numeric(auc(roc_obj)))
}

# Assign the positive class label for analysis
positive_class_label <- "Malignant"
##################################################################################
# OUT-OF-FOLD PREDICTIONS ENSEMBLE MODEL
message("Creating an ensemble model, from 3 all-features trained models....")
##################################################################################
# Selects the best metric from level-0 training ("knn","svmRadial","rf","xgbTree")
best_metric <- model_knn[[1]] %>% select(metric, AUC, F1) %>%   slice_max(AUC) %>% pull(metric)
all_best_metric <- tibble(knn = best_metric)
# get oof for knn model
knn_oof <- model_knn[[2]][[best_metric]]$pred %>%
  filter(k == model_knn[[2]][[best_metric]]$bestTune$k) %>% # filter by best optimizing parameter (k, sigma, mtry, etc)
  group_by(rowIndex) %>%
  dplyr::summarize(Malignant = mean(Malignant)) %>%  # Average predictions from multiple resamples (5 repetitions)
  arrange(rowIndex)

best_metric <- model_svm[[1]] %>% select(metric, AUC, F1) %>%   slice_max(AUC) %>% pull(metric)
all_best_metric <- all_best_metric %>% mutate(all_best_metric,svm = best_metric)
# get oof for svm model
svm_oof <- model_svm [[2]][[best_metric]]$pred %>% 
  filter(sigma == model_svm[[2]][[best_metric]]$bestTune$sigma) %>%
  group_by(rowIndex) %>%
  dplyr::summarize(Malignant = mean(Malignant)) %>%  
  arrange(rowIndex)

best_metric <- model_xgbTree[[1]] %>% select(metric, AUC, F1) %>% slice_max(AUC, with_ties = FALSE) %>% pull(metric)
all_best_metric <- all_best_metric %>% mutate(all_best_metric,xgb = best_metric)
# get oof for xgb model
xgb_oof <- model_xgbTree[[2]][[best_metric]]$pred %>% 
  filter(eta == model_xgbTree[[2]][[best_metric]]$bestTune$eta) %>%
  group_by(rowIndex) %>%
  dplyr::summarize(Malignant = mean(Malignant)) %>%  
  arrange(rowIndex)

meta_features_train <- knn_oof %>%
  left_join(svm_oof, by = "rowIndex", suffix = c(".knn", ".svm")) %>%
  left_join(xgb_oof, by = "rowIndex") %>%
  arrange(rowIndex)

yENS <- y[meta_features_train$rowIndex]  # Should match the length of meta_features_train

# Create grid
step <- 0.005  # adjust resolution as needed
grid <- expand.grid(
  w1 = seq(step, 1, step),
  w2 = seq(step, 1, step)
)
grid$w3 <- 1 - (grid$w1 + grid$w2)

# Keep valid rows
grid <- subset(grid, w3 >= 0)
grid$F1 <- numeric(nrow(grid))

cl <- makeCluster(detectCores() -1) # Create a cluster
registerDoParallel(cl)         # Register it with foreach/caret

# Prepare data set
p <- as.matrix(meta_features_train%>% select(-rowIndex))

# iterate models from predefined weights grid 
f1_AUC <- foreach(i = 1:nrow(grid), .combine = bind_rows,
                  .export = c("grid", "compute_f12", "p","yENS"),
                  .packages = c("caret","pROC","dplyr")
) %dopar% {
  # current set of weights
  current_weights <- grid[i, 1:3]
  # Obtain F1 and AUC for current weight
  compute_f12(current_weights, p, yENS)
}

stopCluster(cl)


aaa <- f1_AUC %>% summarise(b_F1 = which.max(F1),b_AUC = which.max(AUC))
grid$F1 <- f1_AUC$F1

# Best weights
best_row <- grid[aaa$b_F1, ]
best_weights <- as.numeric(best_row[1:3])

cat("Best weights:\n")
print(best_weights)
cat("Best F1:", best_row$F1, "\n")


# VERIFY MODEL ON TEST DATA PREDICTIONS
svm_p <- model_svm[[4]][, "Malignant"]
xgb_p <- model_xgbTree[[4]][, "Malignant"]
knn_p <- model_knn[[4]][, "Malignant"]
rf_p  <- model_rf[[4]][, "Malignant"]


# Use best weights to compute final ensemble prob
ensemble_p_best <- best_weights[1]*svm_p + best_weights[2]*xgb_p + best_weights[3]*knn_p
pred_class <- ifelse(ensemble_p_best >= 0.5, "Malignant", "Benign")
pred_class <- factor(pred_class, levels = levels(y_T))
cm <- confusionMatrix(pred_class, y_T, positive = "Malignant")

# ROC
roc_obj <- roc(y_T, ensemble_p_best, levels = rev(levels(y_T)))
# Store in a overal list for plotting
roc_ensemble <- list(roc_obj)

cat("AUC:", auc(roc_obj), "\n")
cat("Best F1:", best_row$F1, "\n")
cat("Best Acc:", cm$overall[["Accuracy"]])
cat("Best Sens:", cm$byClass[["Sensitivity"]])
cat("Best Spec:", cm$byClass[["Specificity"]])

current_results_w3 <- tibble(
  method = "ensemble wp3 P_OOF",
  metric = "F1",
  Accuracy = cm$overall[["Accuracy"]],
  Sensitivity = cm$byClass[["Sensitivity"]],
  Specificity = cm$byClass[["Specificity"]],
  AUC = as.numeric(auc(roc_obj)),
  F1 = best_row$F1,
  best_p_cv = list(round(best_weights,3)),# list(as.data.frame(t(best_weights))), 
  parameter = "w1,w2,w3"
)

all_model_results <- bind_rows(all_model_results, current_results_w3)
kable(all_model_results)
# kable(current_results_w3)

##################################################
# ENSEMBLE 4 MODELS
message("Creating an ensemble model, from 4 all-features trained models....")
##################################################
best_metric <- model_rf[[1]] %>% select(metric, AUC, F1) %>%   slice_max(AUC) %>% pull(metric)
all_best_metric <- all_best_metric %>% mutate(all_best_metric, rf = best_metric)
# get oof for xgb model
rf_oof <- model_rf[[2]][[best_metric]]$pred %>% 
  filter(mtry == model_rf[[2]][[best_metric]]$bestTune$mtry) %>%
  group_by(rowIndex) %>%
  dplyr::summarize(Malignant = mean(Malignant)) %>%  # Average predictions if multiple resamples
  arrange(rowIndex)

meta_features_train <- knn_oof %>%
  left_join(svm_oof, by = "rowIndex", suffix = c(".knn", ".svm")) %>%
  left_join(xgb_oof, by = "rowIndex") %>%
  left_join(rf_oof, by = "rowIndex") %>%
  arrange(rowIndex)

yENS <- y[meta_features_train$rowIndex]  # Should match the length of meta_features_train

# Create grid
step <- 0.05  # adjust resolution as needed
grid <- expand.grid(
  w1 = seq(step, 1, step),
  w2 = seq(step, 1, step),
  w3 = seq(step, 1, step)
)
grid$w4 <- 1 - (grid$w1 + grid$w2 + grid$w3)

# Keep valid rows
grid <- subset(grid, w4 >= 0)

grid$F1 <- numeric(nrow(grid))

cl <- makeCluster(detectCores() -1) # Create a cluster
registerDoParallel(cl)         # Register it with foreach/caret

p <- as.matrix(meta_features_train%>% select(-rowIndex))

f1_AUC <- foreach(i = 1:nrow(grid), .combine = bind_rows,
                  .export = c("grid", "compute_f12", "p","yENS"),
                  .packages = c("caret","pROC","dplyr")
) %dopar% {
  current_weights <- grid[i, 1:4]
  compute_f12(current_weights, p, yENS)
}

stopCluster(cl)

aaa <- f1_AUC %>% summarise(b_F1 = which.max(F1),b_AUC = which.max(AUC))
grid$F1 <- f1_AUC$F1

# Best weights
best_row <- grid[aaa$b_F1, ]
best_weights <- as.numeric(best_row[1:4])

cat("Best weights:\n")
print(best_weights)
cat("Best F1:", best_row$F1, "\n")

# Use best weights to compute final ensemble prob
ensemble_p_best <- best_weights[1]*svm_p + best_weights[2]*xgb_p + best_weights[3]*knn_p + best_weights[4]*rf_p
pred_class <- ifelse(ensemble_p_best >= 0.5, "Malignant", "Benign")
pred_class <- factor(pred_class, levels = levels(y_T))
cm <- confusionMatrix(pred_class, y_T, positive = "Malignant")

# ROC
roc_obj <- roc(y_T, ensemble_p_best, levels = rev(levels(y_T)))
roc_ensemble[[2]] <- roc_obj
cat("AUC:", auc(roc_obj), "\n")
cat("Best F1:", best_row$F1, "\n")
cat("Best Acc:", cm$overall[["Accuracy"]])
cat("Best Sens:", cm$byClass[["Sensitivity"]])
cat("Best Spec:", cm$byClass[["Specificity"]])

current_results_w4 <- tibble(
  method = "ensemble wp4 P_OOF",
  metric = "F1",
  Accuracy = cm$overall[["Accuracy"]],
  Sensitivity = cm$byClass[["Sensitivity"]],
  Specificity = cm$byClass[["Specificity"]],
  AUC = as.numeric(auc(roc_obj)),
  F1 = best_row$F1,
  best_p_cv = list(round(best_weights,3)),# list(as.data.frame(t(best_weights))), 
  parameter = "w1,w2,w3,w4"
)

all_model_results <- bind_rows(all_model_results, current_results_w4)
kable(all_model_results)
kable(current_results_w4)


##################################################################################
# ENSEMBLE OUT-OF-FOLD PREDICTIONS PCA
message("Creating an ensemble model, from 3 PCA-prepocessed trained models....")
##################################################################################

# best_metric <- model_knn_PCA[[1]]%>%select(metric,AUC,F1) %>% summarize(metric= metric[which.max(AUC)]) %>% pull()
best_metric <- model_knn_PCA[[1]] %>% select(metric, AUC, F1) %>%   slice_max(AUC) %>% pull(metric)
all_best_metric <- all_best_metric %>% mutate(all_best_metric,knn_pca = best_metric)
# get oof for knn model
knn_oof <- model_knn_PCA[[2]][[best_metric]]$pred %>%
  filter(k == model_knn_PCA[[2]][[best_metric]]$bestTune$k) %>%
  group_by(rowIndex) %>%
  dplyr::summarize(Malignant = mean(Malignant)) %>%  # Average predictions if multiple resamples
  arrange(rowIndex)

# best_metric <- model_svm_PCA[[1]]%>%select(metric,AUC,F1) %>% summarize(metric= metric[which.max(AUC)]) %>% pull()
best_metric <- model_svm_PCA[[1]] %>% select(metric, AUC, F1) %>%   slice_max(AUC) %>% pull(metric)
all_best_metric <- all_best_metric %>% mutate(all_best_metric,svm_pca = best_metric)

# get oof for svm model
svm_oof <- model_svm_PCA[[2]][[best_metric]]$pred %>% 
  filter(sigma == model_svm_PCA[[2]][[best_metric]]$bestTune$sigma) %>%
  group_by(rowIndex) %>%
  dplyr::summarize(Malignant = mean(Malignant)) %>%  # Average predictions if multiple resamples
  arrange(rowIndex)

# best_metric <- model_xgbTree_PCA[[1]]%>%select(metric,AUC,F1) %>% summarize(metric= metric[which.max(AUC)]) %>% pull()
best_metric <- model_xgbTree_PCA[[1]] %>% select(metric, AUC, F1) %>%   slice_max(AUC) %>% pull(metric)
all_best_metric <- all_best_metric %>% mutate(all_best_metric,xgb_pca = best_metric)

# get oof for xgb model
xgb_oof <- model_xgbTree_PCA[[2]][[best_metric]]$pred %>% 
  filter(eta == model_xgbTree_PCA[[2]][[best_metric]]$bestTune$eta) %>%
  group_by(rowIndex) %>%
  dplyr::summarize(Malignant = mean(Malignant)) %>%  # Average predictions if multiple resamples
  arrange(rowIndex)

meta_features_train <- knn_oof %>%
  left_join(svm_oof, by = "rowIndex", suffix = c(".knn", ".svm")) %>%
  left_join(xgb_oof, by = "rowIndex") %>%
  arrange(rowIndex)

yENS <- y[meta_features_train$rowIndex]  # Should match the length of meta_features_train

# Create grid
step <- 0.005  # adjust resolution as needed
grid <- expand.grid(
  w1 = seq(step, 1, step),
  w2 = seq(step, 1, step)
)
grid$w3 <- 1 - (grid$w1 + grid$w2)

# Keep valid rows
grid <- subset(grid, w3 >= 0)

grid$F1 <- numeric(nrow(grid))

cl <- makeCluster(detectCores() -1) # Create a cluster
registerDoParallel(cl)         # Register it with foreach/caret

p <- as.matrix(meta_features_train%>% select(-rowIndex))
f1_AUC <- foreach(i = 1:nrow(grid), .combine = bind_rows,
                  .export = c("grid", "compute_f12", "p","yENS"),
                  .packages = c("caret","pROC","dplyr")
) %dopar% {
  current_weights <- grid[i, 1:3]
  compute_f12(current_weights, p, yENS)
  
}
stopCluster(cl)

aaa <- f1_AUC %>% summarise(b_F1 = which.max(F1),b_AUC = which.max(AUC))
grid$F1 <- f1_AUC$F1

# Best weights
best_row <- grid[aaa$b_F1, ]
best_weights <- as.numeric(best_row[1:3])

cat("Best weights:\n")
print(best_weights)
cat("Best F1:", best_row$F1, "\n")

# VERIFY MODEL ON TEST DATA PREDICTIONS
svm_p <- model_svm_PCA[[4]][, "Malignant"]
xgb_p <- model_xgbTree_PCA[[4]][, "Malignant"]
knn_p <- model_knn_PCA[[4]][, "Malignant"]
rf_p  <- model_rf_PCA[[4]][, "Malignant"]


# Use best weights to compute final ensemble prob
ensemble_p_best <- best_weights[1]*svm_p + best_weights[2]*xgb_p + best_weights[3]*knn_p
pred_class <- ifelse(ensemble_p_best >= 0.5, "Malignant", "Benign")
pred_class <- factor(pred_class, levels = levels(y_T))
cm <- confusionMatrix(pred_class, y_T, positive = "Malignant")

# ROC
roc_obj <- roc(y_T, ensemble_p_best, levels = rev(levels(y_T)))
roc_ensemble[[3]] <- roc_obj
cat("AUC:", auc(roc_obj), "\n")
cat("Best F1:", best_row$F1, "\n")
cat("Best Acc:", cm$overall[["Accuracy"]])
cat("Best Sens:", cm$byClass[["Sensitivity"]])
cat("Best Spec:", cm$byClass[["Specificity"]])

current_results_w3_PCA <- tibble(
  method = "ensemble wp3 P_OOF PCA",
  metric = "F1",
  Accuracy = cm$overall[["Accuracy"]],
  Sensitivity = cm$byClass[["Sensitivity"]],
  Specificity = cm$byClass[["Specificity"]],
  AUC = as.numeric(auc(roc_obj)),
  F1 = best_row$F1,
  best_p_cv = list(round(best_weights,3)),# list(as.data.frame(t(best_weights))), 
  parameter = "w1,w2,w3"
)

all_model_results_PCA <- bind_rows(all_model_results_PCA, current_results_w3_PCA)
kable(all_model_results_PCA)
kable(current_results_w3_PCA)

##################################################
# ENSEMBLE 4 MODELS PCA
message("Creating an ensemble model, from 4 PCA-prepocessed trained models....")
##################################################
# best_metric <- model_rf_PCA[[1]]%>%select(metric,AUC,F1) %>% summarize(metric= metric[which.max(AUC)]) %>% pull()
best_metric <- model_rf_PCA[[1]] %>% select(metric, AUC, F1) %>%   slice_max(AUC) %>% pull(metric)
all_best_metric <- all_best_metric %>% mutate(all_best_metric,rf_pca = best_metric)

# get oof for xgb model
rf_oof <- model_rf_PCA[[2]][[best_metric]]$pred %>% 
  filter(mtry == model_rf_PCA[[2]][[best_metric]]$bestTune$mtry) %>%
  group_by(rowIndex) %>%
  dplyr::summarize(Malignant = mean(Malignant)) %>%  # Average predictions if multiple resamples
  arrange(rowIndex)


meta_features_train <- knn_oof %>%
  left_join(svm_oof, by = "rowIndex", suffix = c(".knn", ".svm")) %>%
  left_join(xgb_oof, by = "rowIndex") %>%
  left_join(rf_oof, by = "rowIndex") %>%
  arrange(rowIndex)

yENS <- y[meta_features_train$rowIndex]  # Should match the length of meta_features_train

# Create grid
step <- 0.05  # adjust resolution as needed
grid <- expand.grid(
  w1 = seq(step, 1, step),
  w2 = seq(step, 1, step),
  w3 = seq(step, 1, step)
)
grid$w4 <- 1 - (grid$w1 + grid$w2 + grid$w3)

# Keep valid rows
grid <- subset(grid, w4 >= 0)
grid$F1 <- numeric(nrow(grid))

cl <- makeCluster(detectCores() -1) # Create a cluster
registerDoParallel(cl)         # Register it with foreach/caret

p <- as.matrix(meta_features_train%>% select(-rowIndex))
f1_AUC <- foreach(i = 1:nrow(grid), .combine = bind_rows,
                  .export = c("grid", "compute_f12", "p","yENS"),
                  .packages = c("caret","pROC","dplyr")
) %dopar% {
  current_weights <- grid[i, 1:4]
  compute_f12(current_weights, p, yENS)
  
}
stopCluster(cl)

aaa <- f1_AUC %>% summarise(b_F1 = which.max(F1),b_AUC = which.max(AUC))
grid$F1 <- f1_AUC$F1

# Best weights
best_row <- grid[aaa$b_F1, ]
best_weights <- as.numeric(best_row[1:4])

cat("Best weights:\n")
print(best_weights)
cat("Best F1:", best_row$F1, "\n")

# Use best weights to compute final ensemble prob
ensemble_p_best <- best_weights[1]*svm_p + best_weights[2]*xgb_p + best_weights[3]*knn_p + best_weights[4]*rf_p
pred_class <- ifelse(ensemble_p_best >= 0.5, "Malignant", "Benign")
pred_class <- factor(pred_class, levels = levels(y_T))
cm <- confusionMatrix(pred_class, y_T, positive = "Malignant")

# ROC
roc_obj <- roc(y_T, ensemble_p_best, levels = rev(levels(y_T)))
roc_ensemble[[4]] <- roc_obj
cat("AUC:", auc(roc_obj), "\n")
cat("Best F1:", best_row$F1, "\n")
cat("Best Acc:", cm$overall[["Accuracy"]])
cat("Best Sens:", cm$byClass[["Sensitivity"]])
cat("Best Spec:", cm$byClass[["Specificity"]])

current_results_w4_PCA <- tibble(
  method = "ensemble wp4 P_OOF PCA",
  metric = "F1",
  Accuracy = cm$overall[["Accuracy"]],
  Sensitivity = cm$byClass[["Sensitivity"]],
  Specificity = cm$byClass[["Specificity"]],
  AUC = as.numeric(auc(roc_obj)),
  F1 = best_row$F1,
  best_p_cv = list(round(best_weights,3)),# list(as.data.frame(t(best_weights))), 
  parameter = "w1,w2,w3,w4"
)

all_model_results_PCA <- bind_rows(all_model_results_PCA, current_results_w4_PCA)
kable(all_model_results_PCA)
kable(current_results_w4_PCA)


#################
# summary of ensembles
all_ensembles <- bind_rows(
  current_results_w3,
  current_results_w4,
  current_results_w3_PCA,
  current_results_w4_PCA
)
kable(all_ensembles)

##############################################################################
# Storage ROC data for plotting
##############################################################################
roc_svm <- model_svm[[3]][[all_best_metric$svm]]
roc_xgb <- model_xgbTree[[3]][[all_best_metric$xgb]]
roc_knn <- model_knn[[3]][[all_best_metric$knn]]
roc_rf <- model_rf[[3]][[all_best_metric$rf]]

roc_svm_pca <- model_svm_PCA[[3]][[all_best_metric$svm_pca]]
roc_xgb_pca <- model_xgbTree_PCA[[3]][[all_best_metric$xgb_pca]]
roc_knn_pca <- model_knn_PCA[[3]][[all_best_metric$knn_pca]]
roc_rf_pca <- model_rf_PCA[[3]][[all_best_metric$rf_pca]]

roc_list <- list(SVM = roc_svm, xgbTree = roc_xgb, KNN=roc_knn, RF=roc_rf, ENS_W3 = roc_ensemble[[1]] , ENS_W4 = roc_ensemble[[2]])
roc_list_pca <- list(SVM = roc_svm_pca, xgbTree = roc_xgb_pca, KNN=roc_knn_pca, RF=roc_rf_pca, ENS_W3 = roc_ensemble[[3]] , ENS_W4 = roc_ensemble[[4]])

##############################################################################
# save data for report or temporary development
##############################################################################

save(model_knn_PCA,model_svm_PCA,model_rf_PCA,model_xgbTree_PCA,
     model_knn,model_svm,model_xgbTree,model_rf,all_model_results_PCA,
     all_model_results, all_ensembles, roc_list,roc_list_pca,
     train_set,test_set, file="trained_models_ALL.RData")

##############################################################################
# PLOT ROC CURVES: WITH OR WITHOUT PCA
##############################################################################
ggroc(roc_list, aes = c("color", "linetype")) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed") + # Add diagonal line
  labs(
    title = "ROC Curves for Breast Cancer Classification Models",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  scale_color_manual(values = c("blue", "red", "darkgreen","magenta","black","cyan")) + # Customize colors
  theme_minimal() +
  annotate("text", x = 0.75, y = 0.35, label = paste("AUC SVM:", round(auc(roc_svm), 4)), col = "blue") +
  annotate("text", x = 0.75, y = 0.30, label = paste("AUC XGB:", round(auc(roc_xgb), 4)), col = "red") +
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC KNN:", round(auc(roc_knn), 4)), col = "darkgreen")+
  annotate("text", x = 0.75, y = 0.20, label = paste("AUC RF:", round(auc(roc_rf), 4)), col = "magenta")+
  annotate("text", x = 0.75, y = 0.15, label = paste("AUC ENS W3:", round(as.numeric(roc_ensemble[[1]]$auc), 4)), col = "black")+
  annotate("text", x = 0.75, y = 0.10, label = paste("AUC ENS W4:", round(as.numeric(roc_ensemble[[2]]$auc), 4)), col = "cyan")+
  geom_line(linewidth = 1.2) +# Apply linewidth to the ROC curves themselves
  theme_light()+
  coord_cartesian(xlim = c(1.05, 0.65), ylim = c(0.65,1.05))

# PCA DATA

ggroc(roc_list_pca, aes = c("color", "linetype")) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed") + # Add diagonal line
  labs(
    title = "ROC Curves for Breast Cancer Classification Models (PCA data)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  scale_color_manual(values = c("blue", "red", "darkgreen","magenta","black","cyan")) + # Customize colors
  theme_minimal() +
  annotate("text", x = 0.75, y = 0.35, label = paste("AUC SVM PCA:", round(auc(roc_svm), 4)), col = "blue") +
  annotate("text", x = 0.75, y = 0.30, label = paste("AUC XGB PCA:", round(auc(roc_xgb), 4)), col = "red") +
  annotate("text", x = 0.75, y = 0.25, label = paste("AUC KNN PCA:", round(auc(roc_knn), 4)), col = "darkgreen")+
  annotate("text", x = 0.75, y = 0.20, label = paste("AUC RF PCA:", round(auc(roc_rf), 4)), col = "magenta")+
  annotate("text", x = 0.75, y = 0.15, label = paste("AUC ENS W3 PCA:", round(as.numeric(roc_ensemble[[1]]$auc), 4)), col = "black")+
  annotate("text", x = 0.75, y = 0.10, label = paste("AUC ENS W4 PCA:", round(as.numeric(roc_ensemble[[2]]$auc), 4)), col = "cyan")+
  geom_line(linewidth = 1.2) +# Apply linewidth to the ROC curves themselves
  theme_light() +
  coord_cartesian(xlim = c(1.05, 0.65), ylim = c(0.65,1.05))

# aaa<-bind_rows(all_model_results,all_model_results_PCA) %>% arrange(desc(AUC))

##############################################################################
# PLOT HEAT MAPS: WITH OR WITHOUT PCA
##############################################################################

all_model_results_H <- all_model_results %>% select(-best_p_cv,-parameter)
# Select best model per class (adjust criteria as needed)
best_models <- all_model_results_H %>%
  group_by(method) %>%
  slice_max(AUC, n = 1) %>%  # or another criterion like F1
  ungroup()

# Prepare data for heatmap
heatmap_data <- best_models %>%
  select(method, Accuracy, Sensitivity, Specificity, AUC, F1) %>%
  pivot_longer(cols = c(Accuracy, Sensitivity, Specificity, AUC, F1),
               names_to = "Metric",
               values_to = "Value")

# Plot heatmap
ggplot(heatmap_data, aes(x = Metric, y = method, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", Value)), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0.85, 1)) +
  theme_minimal() +
  labs(title = "Best Model per Class Performance Heatmap",
       x = "Metric",
       y = "Model",
       fill = "Score")

all_model_results_PCA_H <- all_model_results_PCA %>% select(-best_p_cv,-parameter)
# Select best model per class (adjust criteria as needed)
best_models <- all_model_results_PCA_H %>%
  group_by(method) %>%
  slice_max(AUC, n = 1) %>%  # or another criterion like F1
  ungroup()

# Prepare data for heatmap
heatmap_data <- best_models %>%
  select(method, Accuracy, Sensitivity, Specificity, AUC, F1) %>%
  pivot_longer(cols = c(Accuracy, Sensitivity, Specificity, AUC, F1),
               names_to = "Metric",
               values_to = "Value")

# Plot heatmap
ggplot(heatmap_data, aes(x = Metric, y = method, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", Value)), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0.85, 1)) +
  theme_minimal() +
  labs(title = "Best Model per Class Performance Heatmap (PCA)",
       x = "Metric",
       y = "Model",
       fill = "Score")



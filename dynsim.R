rm(list = ls())

# Isingのデータを使ってダイナミカルにシミュレーションする
library(stats)
library(qgraph)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(glmnet)
library(foreign)
library(bootnet)
library(IsingFit)

# データの読み込み
data <- read.spss("pone.0182162.s004.sav", 
                  to.data.frame=TRUE)

# データの整理(PHQ-9に絞って実行)
data_phq <- data %>% 
  mutate(phq9a = S_PHQ9_a, phq9b = S_PHQ9_b, 
         phq9c = S_PHQ9_c, phq9d = S_PHQ9_d, 
         phq9e = S_PHQ9_e, phq9f = S_PHQ9_f, 
         phq9g = S_PHQ9_g, phq9h = S_PHQ9_h, 
         phq9i = S_PHQ9_i) %>% 
  select(phq9a, phq9b, phq9c, phq9d, phq9e, phq9f, phq9g, phq9h, phq9i) %>%
  binarize(split = 1)

res <- estimateNetwork(data_phq, "IsingSampler")



simulate_treatment_network_dym <- function(W_init,
                                           b_init,
                                           target,
                                           connectivity = 1,
                                           edge_between_TC = 1,
                                           weight_bias = 1,
                                           TB = 1,
                                           trial = 10,
                                           baseline_iteration = 10,
                                           num_TC = 5,
                                           TC_iteration_per_component = 10,
                                           follow_up_iteration = 10,
                                           symptom_name = NULL) {
  
  
  
  # --- 内部関数: ロジスティック関数 ---
  sigmoid <- function(x) {
    x <- pmax(pmin(x, 700), -700)
    1 / (1 + exp(x)) 
  }
  
  # --- 内部関数: 症状ネットワークの推定 (正則化ロジスティック回帰) ---
  estimate_symptom_network <- function(X_series) {
    P <- ncol(X_series)
    T_time <- nrow(X_series)
    
    Y_data <- as.matrix(X_series[2:T_time, ])
    Z_data <- as.matrix(X_series[1:(T_time - 1), ])
    
    B_temporal <- matrix(0, nrow = P, ncol = P)
    Alpha_threshold <- numeric(P)
    
    for (i in 1:P) {
      Y_i <- Y_data[, i]
      
      if (sd(Y_i) == 0) {
        Alpha_threshold[i] <- 0 
        next
      }
      
      tryCatch({
        n_folds <- min(10, length(Y_i))
        if(n_folds < 3) n_folds <- length(Y_i)
        
        fit <- suppressWarnings(cv.glmnet(x = Z_data, y = Y_i, family = "binomial", alpha = 1, nfolds = n_folds))
        
        coefs <- coef(fit, s = "lambda.min")
        coefs_vec <- as.vector(coefs)
        
        Alpha_threshold[i] <- coefs_vec[1]
        B_temporal[i, ] <- coefs_vec[2:(P + 1)]
        
      }, error = function(e) {
      })
    }
    
    return(list(B = B_temporal, Alpha = Alpha_threshold))
  }
  
  # --- 内部関数: 1試行分のシミュレーション ---
  simulate_network <- function(W_init, b_init, target, connectivity, edge_between_TC,
                               weight_bias, TB, baseline_iteration, num_TC,
                               TC_iteration_per_component, follow_up_iteration, symptom_name) {
    
    W <- W_init
    b <- b_init
    
    # TCの重み生成のための平均・分散計算用
    if (isSymmetric(as.matrix(W))) {
      # 無向グラフ: 上三角成分のみ取得（重複回避）
      edges <- W[row(W) < col(W)]
    } else {
      # 有向グラフ: 対角成分以外の非ゼロ要素を全て取得
      edges <- W[row(W) != col(W) & W != 0]
    }
    
    if(length(edges) == 0) {
      mean_weight <- 0
      sd_weight <- 0.1 
    } else {
      mean_weight <- mean(edges, na.rm = TRUE)
      sd_weight <- sd(edges, na.rm = TRUE)
      if(is.na(sd_weight)) sd_weight <- 0.1
    }
    
    mean_threthold <- mean(b, na.rm = TRUE)
    sd_threthold <- sd(b, na.rm = TRUE)
    if(is.na(sd_threthold)) sd_threthold <- 0.1
    
    num_symptom <- nrow(W)
    X <- matrix(rep(0, length(b)), 1, length(b))
    
    D <- numeric(0)
    TC <- numeric(0)
    
    if (is.null(symptom_name)) {
      symptom_name <- letters[1:num_symptom]
    }
    
    safe_sample <- function(prob) {
      if (is.na(prob)) prob <- 0.5 
      if (prob < 0) prob <- 0
      if (prob > 1) prob <- 1
      sample(x = c(1, 0), size = 1, replace = T, prob = c(prob, 1 - prob))
    }
    
    # --- Baseline ---
    X_history_baseline <- matrix(0, nrow = baseline_iteration, ncol = length(b))
    
    for (s in 1:baseline_iteration) {
      A <- numeric(length(b))
      P_prob <- numeric(length(b))
      
      for (j in 1:length(b)) {
        
        input_val <- sum(connectivity * W[, j] * X)
        A[j] <- input_val
        exponent_val <- abs(b[j]) - A[j]
        exponent_val <- max(min(exponent_val, 100), -100)
        P_prob[j] <- 1 / (1 + exp(exponent_val))
        X[j] <- safe_sample(P_prob[j])
      }
      X_history_baseline[s, ] <- X
      D[s] <- sum(X[1:num_symptom]) / num_symptom
      TC[s] <- 0
    }
    
    last_phase_data <- X_history_baseline
    
    # --- TC Loop ---
    for (i in 1:num_TC) {
      
      #  次のTCが入る直前までの各症状の時系列データからネットワーク更新
      if (i > 1 && any(target == 1)) {
        symptom_data <- last_phase_data[, 1:num_symptom]
        
        if (sd(as.vector(symptom_data)) > 0) {
          estimated <- estimate_symptom_network(symptom_data)
          B_new <- estimated$B
          Alpha_new <- estimated$Alpha
          
          B_new[is.na(B_new)] <- 0
          Alpha_new[is.na(Alpha_new)] <- 0
          
          W[1:num_symptom, 1:num_symptom] <- t(B_new)
          b[1:num_symptom] <- Alpha_new
        }
      }
      
      # Wの拡張
      current_nodes_count <- nrow(W)
      new_dim <- current_nodes_count + 1
      W_expanded <- matrix(0, nrow = new_dim, ncol = new_dim)
      if (current_nodes_count > 0) {
        W_expanded[1:current_nodes_count, 1:current_nodes_count] <- W
      }
      W <- W_expanded
      
      # TCを入れるとtargetと負の相関を持つようになる
      for (j in 1:(num_symptom + i)) {
        if (j <= num_symptom) {
          if (target[j] == 1) {
            mw <- if(is.na(mean_weight)) 0 else mean_weight
            sw <- if(is.na(sd_weight)) 0.1 else sd_weight
            
            temp_weight <- -abs(rnorm(1, mean_weight, sw)) * weight_bias
            
            W[num_symptom + i, j] <- temp_weight
            W[j, num_symptom + i] <- temp_weight
          }
        } else if (j > num_symptom && j < (num_symptom + i)) {
          mw <- if(is.na(mean_weight)) 0 else mean_weight
          sw <- if(is.na(sd_weight)) 0.1 else sd_weight
          temp_weight <- rnorm(1, edge_between_TC * mw, sw)
          W[num_symptom + i, j] <- temp_weight
          W[j, num_symptom + i] <- temp_weight
        }
      }
      
      mt <- if(is.na(mean_threthold)) 0 else mean_threthold
      st <- if(is.na(sd_threthold)) 0.1 else sd_threthold
      b <- c(b, rnorm(1, mt * TB, st))
      
      TC_name <- paste0("TC", 1:i)
      colnames(W) <- c(symptom_name, TC_name)
      rownames(W) <- c(symptom_name, TC_name)
      
      X <- cbind(X, c(0))
      
      X_history_current <- matrix(0, nrow = TC_iteration_per_component, ncol = length(b))
      
      for (k in 1:TC_iteration_per_component) {
        A <- numeric(length(b))
        P_prob <- numeric(length(b))
        for (j in 1:length(b)) {
          input_val <- sum(connectivity * W[, j] * X)
          A[j] <- input_val
          exponent_val <- abs(b[j]) - A[j]
          exponent_val <- max(min(exponent_val, 100), -100)
          P_prob[j] <- 1 / (1 + exp(exponent_val))
          X[j] <- safe_sample(P_prob[j])
        }
        X_history_current[k, ] <- X
        tmp_D <- sum(X[1:num_symptom]) / num_symptom
        tmp_TC <- sum(X[(num_symptom + 1):(num_symptom + i)]) / i
        time_no <- baseline_iteration + (i - 1) * TC_iteration_per_component + k
        D[time_no] <- tmp_D
        TC[time_no] <- tmp_TC
      }
      last_phase_data <- X_history_current
    }
    
    # --- Follow-up ---
    for (k in 1:follow_up_iteration) {
      A <- numeric(length(b))
      P_prob <- numeric(length(b))
      for (j in 1:length(b)) {
        input_val <- sum(connectivity * W[, j] * X)
        A[j] <- input_val
        exponent_val <- abs(b[j]) - A[j]
        exponent_val <- max(min(exponent_val, 100), -100)
        P_prob[j] <- 1 / (1 + exp(exponent_val))
        X[j] <- safe_sample(P_prob[j])
      }
      tmp_D <- sum(X[1:num_symptom]) / num_symptom
      tmp_TC <- sum(X[(num_symptom + 1):(num_symptom + i)]) / i
      time_no <- length(D) + 1
      D[time_no] <- tmp_D
      TC[time_no] <- tmp_TC
    }
    
    return(list(D = D, TC = TC, W = W))
  }
  
  # --- メイン処理 ---
  target <- unlist(target)
  total_time <- baseline_iteration + num_TC * TC_iteration_per_component + follow_up_iteration
  D_iteration <- matrix(0, trial, total_time)
  TC_iteration <- matrix(0, trial, total_time)
  
  initial_num_symptoms <- nrow(W_init)
  final_nodes <- initial_num_symptoms + num_TC
  W_iteration <- matrix(0, trial, final_nodes^2)
  
  # 症状名の取得
  if (is.null(symptom_name)) {
    if (!is.null(colnames(W_init))) {
      symptom_name <- colnames(W_init)
    } else {
      symptom_name <- letters[1:initial_num_symptoms]
    }
  }
  
  for (l in 1:trial) {
    result <- simulate_network(W_init, b_init, target, connectivity, edge_between_TC,
                               weight_bias, TB, baseline_iteration, num_TC,
                               TC_iteration_per_component, follow_up_iteration,
                               symptom_name)
    D_iteration[l, ] <- result$D
    TC_iteration[l, ] <- result$TC
    W_iteration[l, ] <- as.vector(result$W)
  }
  
  W_names_full <- colnames(result$W) 
  D_plot <- colMeans(D_iteration, na.rm=TRUE)
  TC_plot <- colMeans(TC_iteration, na.rm=TRUE)
  D_sd <- apply(D_iteration, 2, sd, na.rm=TRUE)
  TC_sd <- apply(TC_iteration, 2, sd, na.rm=TRUE)
  
  # --- プロット用データの整形（症状のみ抽出） ---
  W_full_mean <- matrix(colMeans(W_iteration, na.rm=TRUE), final_nodes, final_nodes)
  colnames(W_full_mean) <- W_names_full
  rownames(W_full_mean) <- W_names_full
  
  W_symptoms_only <- W_full_mean[1:initial_num_symptoms, 1:initial_num_symptoms]
  
  color_node <- rep("lightblue", initial_num_symptoms)
  
  # 固定サイズ設定
  v_fixed <- 8
  
  # タイトル設定
  target_indices <- which(target == 1)
  if (length(target_indices) > 0) {
    targeted_names <- symptom_name[target_indices]
    plot_title <- paste("Intervention on:", paste(targeted_names, collapse = ", "))
  } else {
    plot_title <- "No Intervention"
  }
  
  temp_file_path <- file.path(tempdir(), "tmp.tiff")
  tiff(temp_file_path, width = 1200, height = 800, res = 400)
  tryCatch({
    qgraph(W_symptoms_only, layout = "spring", theme = "colorblind", 
           vsize = v_fixed,
           color = color_node, posCol = "blue", 
           negCol = "red", negDashed = F,
           title = plot_title) 
  }, error = function(e) {
    plot(1, type="n", axes=F, xlab="", ylab="", main="Network Plot Error")
  })
  dev.off()
  
  p2 <- tibble(time = 1:total_time,
               D_plot = D_plot,
               TC_plot = TC_plot,
               D_sd = D_sd,
               TC_sd = TC_sd) |>
    pivot_longer(cols = c(D_plot, TC_plot), names_to = "variable", values_to = "mean") |>
    mutate(sd = ifelse(variable == "D_plot", D_sd, TC_sd)) |>
    ggplot(aes(x = time, y = mean, color = variable, fill = variable)) +
    geom_line() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
    scale_fill_manual(name = "Variable", values = c("D_plot" = "blue", "TC_plot" = "orange")) +
    scale_color_manual(name = "Variable", values = c("D_plot" = "blue", "TC_plot" = "orange")) +
    scale_y_continuous(limits = c(0.0, 1.0), oob = scales::squish) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(x = "", y = "", title = "")
  
  qgraph_image <- ggdraw() + draw_image(temp_file_path)
  qgraph_image_grob <- ggplotGrob(qgraph_image)
  result_plot <- grid.arrange(qgraph_image_grob, p2, ncol = 2)
  
  result_text <- paste0(
    sprintf("The mean value of symptom at the final step(t=%d). = %.3f", total_time, D_plot[total_time]), "\n",
    sprintf("The mean value of treatment component at the final step(t=%d). = %.3f", total_time, TC_plot[total_time]), "\n",
    sprintf("The SD value of symptom at the final step(t=%d). = %.3f", total_time, D_sd[total_time]), "\n",
    sprintf("The SD value of treatment component at the final step(t=%d). = %.3f", total_time, TC_sd[total_time])
  )
  return(list(result_plot = result_plot, result_text = result_text))
}


# 推定したデータを使ってシミュレーション
simulate_treatment_network_dym(res$graph,
                               res$intercepts,
                               target = c(1,0,0,0,0,0,0,0,0),
                               connectivity = 1.3,
                               edge_between_TC = 1,
                               weight_bias = 1,
                               TB = 1,
                               trial = 10,
                               baseline_iteration = 10,
                               num_TC = 9,
                               TC_iteration_per_component = 10,
                               follow_up_iteration = 10)

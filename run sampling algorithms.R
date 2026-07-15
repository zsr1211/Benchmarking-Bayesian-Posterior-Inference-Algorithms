# This script is about running all posteirors in posteriordb with 4 methods 
# (1 MCMC and 3 PDMP) and 10 seeds, and saving the running reuslts as .rds files.



# load in packages
library(posteriordb)
library(cmdstanr)
library(PDMPSamplersR)
library(posterior)
library(bridgestan)
library(dplyr)



# update PDMPSamplersR and the Julia backbench if necessary
# PDMPSamplersR::pdmpsamplers_update()
# JuliaCall::julia_command("Pkg.status()")



# set posteriordb and select some posteriors
pdb <- pdb_local("/Users/zhangsr/Desktop/UvA/Thesis/posteriorDB/posteriordb")



# select posterior without reference draws from posteriorDB
summary_df_noref <- summary_df %>%
  filter(has_ref == 0) %>% 
  select(posterior, cond_laplace, cond_reference, has_unc_ref, has_ref, complexity_withNA, everything())

# Posterior groups (without posteriorDB reference)
posterior_groups_noref <- list(
  simple = summary_df_noref$posterior[summary_df_noref$complexity_withNA == "simple"],
  moderate = summary_df_noref$posterior[summary_df_noref$complexity_withNA == "moderate"],
  complex = summary_df_noref$posterior[summary_df_noref$complexity_withNA == "complex"]
)
posterior_groups_noref



# select posterior with reference draws from posteriorDB
summary_df_ref <- summary_df %>%
  filter(has_ref == 1) %>%
  select(posterior, cond_laplace, cond_reference, has_unc_ref, has_ref, complexity_withNA, everything())

# Posterior groups (with posteriorDB reference)
posterior_groups_ref <- list(
  simple = summary_df_ref$posterior[summary_df_ref$complexity_withNA == "simple"],
  moderate = summary_df_ref$posterior[summary_df_ref$complexity_withNA == "moderate"],
  complex = summary_df_ref$posterior[summary_df_ref$complexity_withNA == "complex"]
)
posterior_groups_ref



# select posterior with my own reference draws (from try_to_find_reference.R)
summary_df_myref <- summary_df %>%
  filter(has_unc_ref == 2) %>% 
  select(posterior, cond_laplace, cond_reference, has_unc_ref, has_ref, complexity_withNA, everything())

# Posterior groups (without posteriorDB reference)
posterior_groups_myref <- list(
  simple = summary_df_myref$posterior[summary_df_myref$complexity_withNA == "simple"],
  moderate = summary_df_myref$posterior[summary_df_myref$complexity_withNA == "moderate"],
  complex = summary_df_myref$posterior[summary_df_myref$complexity_withNA == "complex"]
)
posterior_groups_myref



# this posterior subset is used to first batch to run
posterior_groups_first <- list(
  simple = c(
    "arK-arK",
    "arma-arma11",
    "bball_drive_event_0-hmm_drive_0",
    "eight_schools-eight_schools_noncentered",
    "garch-garch11",
    "gp_pois_regr-gp_regr",
    "hudson_lynx_hare-lotka_volterra",
    "low_dim_gauss_mix-low_dim_gauss_mix",
    "one_comp_mm_elim_abs-one_comp_mm_elim_abs"
  ),
  moderate = c(
    "diamonds-diamonds",
    "gp_pois_regr-gp_pois_regr",
    "hmm_example-hmm_example",
    "nes2000-nes",
    "sblrc-blr"
  ),
  complex = c(
    "bball_drive_event_1-hmm_drive_1",
    "earnings-logearn_interaction",
    "eight_schools-eight_schools_centered",
    "mcycle_gp-accel_gp"
  )
)

# # the second batch
# posterior_groups_rest <- lapply(
#   names(posterior_groups_ref),
#   function(g) setdiff(posterior_groups[[g]], posterior_groups_first[[g]])
# )
# names(posterior_groups_rest) <- names(posterior_groups)
# posterior_groups_rest



# this posterior subset is used to check if the pipeline works
# posterior_groups <- list(
#   simple = c(
#     "arK-arK"),
#   moderate = c(
#     "diamonds-diamonds"
#   ),
#   complex = c(
#     "eight_schools-eight_schools_centered"
#   )
# )


# select some posteriors with smaller data set and no reference draws 
# neither posteriorDB reference, nor the reference we generated
posterior_groups_noref <- list(
  simple = c(
    "Mth_data-Mth_model",
    "Mh_data-Mh_model"
  ),
  moderate = c(
    "surgical_data-surgical_model",
    "sir-sir",
    "soil_carbon-soil_incubation",
    "seeds_data-seeds_model",
    "loss_curves-losscurve_sislob",
    "low_dim_gauss_mix_collapse-low_dim_gauss_mix_collapse"
  ),
  complex = c(
    "bball_drive_event_1-hmm_drive_1",
    "Survey_data-Survey_model",
    "GLMM_Poisson_data-GLMM_Poisson_model",
    "seeds_data-seeds_centered_model",
    "lsat_data-lsat_model",
    "pilots-pilots",
    "seeds_data-seeds_stanified_model",
    "rats_data-rats_model",
    "butterfly-multi_occupancy",
    "uk_drivers-state_space_stochastic_level_stochastic_seasonal"
  )
)




# # check if the worst BPS posterior
# BPS_check <- list(
#   complex = c(
#     "mesquite-mesquite",
#     "earnings-logearn_interaction"
#   )
# )




# fixed seed sets
set.seed(1234)
seeds <- sample(1:10000, 10, replace = FALSE) 
# seeds <- seeds[6]
seeds




# set methods
methods <- c("MCMC", "ZigZag", "BPS", "Boomerang")
# methods <- c("MCMC")
# methods <- c("BPS")




# set output folder and error file
results_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_final"
results_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_check_BPS"
dir.create(results_dir, showWarnings = FALSE)

error_dir <- file.path(results_dir, "error_logs")
dir.create(error_dir, showWarnings = FALSE, recursive = TRUE)




# main loops
# check if Julia setup works
PDMPSamplersR:::check_for_julia_setup()

# get pdmp dimension
get_pdmp_dimension <- function(model_path, data_path) {
  
  JuliaCall::julia_assign("_path_to_stan_model", normalizePath(model_path, mustWork = TRUE))
  JuliaCall::julia_assign("_path_to_stan_data", normalizePath(data_path, mustWork = TRUE))
  JuliaCall::julia_command("_tmp_pdmp_model = PDMPModel(_path_to_stan_model, _path_to_stan_data);")
  
  d <- JuliaCall::julia_eval("_tmp_pdmp_model.d")
  JuliaCall::julia_command("_tmp_pdmp_model = nothing; GC.gc()")
  d
}



# loops
# run each posterior groups
# change the posterior_groups_myref to other posterior subset
for (group_name in names(posterior_groups_myref)) {
  
  # run each posterior
  for (post_name in posterior_groups_myref[[group_name]]) {
    
    cat("\n==========\n")
    cat("Group:", group_name, "Posterior:", post_name, "\n")
    
    # get posterior
    po <- posterior(post_name, pdb)
    
    # get stan model and data path
    model_path <- stan_code_file_path(po)
    data_path <- data_file_path(po)
    
    cat("Model:", model_path, "\n")
    cat("Data:", data_path, "\n")
    
    # make posterior name safe for file paths
    safe_post_name <- gsub("[^A-Za-z0-9_\\-]", "_", post_name)
    
    # compile model
    mod <- cmdstan_model(model_path, force_recompile = TRUE)
    
    # get dimension
    d <- get_pdmp_dimension(model_path, data_path)
    cat("Dimension:", d, "\n")
    
    
    # run each seed
    for (seed in seeds) {
      
      file_x0 <- file.path(
        results_dir,
        paste0(group_name, "_", safe_post_name, "_seed", seed, "_x0.rds")
      )
      
      if (file.exists(file_x0)) {
        x0_fixed <- readRDS(file_x0)
        
        if (length(x0_fixed) != d) {
          stop(
            "Existing x0 has wrong dimension for posterior ", post_name,
            ". Expected d = ", d, ", but got length = ", length(x0_fixed)
          )
        }
        
        cat("Loaded existing x0:", file_x0, "\n")
      } else {
        set.seed(seed)
        x0_fixed <- rnorm(d)
        saveRDS(x0_fixed, file_x0)
        cat("Saved new x0:", file_x0, "\n")
      }
      
      
      # run each method (algorithm)
      for (method in methods) {
        
        cat("\nMethod:", method, "Seed:", seed, "\n")
        
        # output file for main result
        file_out <- file.path(
          results_dir,
          paste0(
            group_name, "_", safe_post_name, "_", method, "_seed", seed, ".rds"
          )
        )
        
        # output file for MCMC unconstrained draws
        file_draws <- file.path(
          results_dir,
          paste0(
            group_name, "_", safe_post_name, "_", method, "_seed", seed,
            "_unconstrained_matrix.rds"
          )
        )
        
        # error log file
        error_out <- file.path(
          error_dir,
          paste0(
            group_name, "_", safe_post_name, "_", method, "_seed", seed,
            "_ERROR.txt"
          )
        )
        
        # skip logic
        if (method == "MCMC") {
          if (file.exists(file_out) && file.exists(file_draws)) {
            cat("MCMC result and unconstrained draws already exist → skip\n")
            next
          } else if (file.exists(file_out) && !file.exists(file_draws)) {
            cat("MCMC result exists but unconstrained draws missing → rerun MCMC\n")
          }
        } else {
          if (file.exists(file_out)) {
            cat("Result already exists → skip\n")
            next
          }
        }
        
        tryCatch({
          
          # settings for each method
          if (method == "MCMC") {
            
            fit <- mod$sample(
              data = data_path,
              chains = 4,
              parallel_chains = 4,
              iter_warmup = 2000, 
              iter_sampling = 2000,
              seed = seed,
              adapt_delta = 0.9
            )
            
            result <- fit
            
            # save MCMC result first
            saveRDS(result, file_out)
            cat("Saved MCMC result:", file_out, "\n")
            
            # extract and save unconstrained draws
            unconstrained_matrix <- fit$unconstrain_draws(format = "matrix")
            saveRDS(unconstrained_matrix, file_draws)
            cat("Saved unconstrained draws:", file_draws, "\n")
          }
          
          
          if (method == "ZigZag") {
            result <- pdmp_sample_from_stanmodel(
              model_path,
              data_path,
              flow = "PreconditionedZigZag",
              algorithm = "GridThinningStrategy",
              T = 2000,
              grid_n = 30,
              t_warmup = 1000,
              sticky = FALSE,
              adaptive_scheme = "full",
              seed = seed,
              x0 = x0_fixed,
              support_boundary = support_boundary_control(
                mode = "line_search_truncated_refresh",
                max_refresh_attempts = 100L,
                refresh_probe_time = 1e-6,
                min_safe_time = 1e-8
              )
            )

            saveRDS(result, file_out)
            cat("Saved:", file_out, "\n")
          }

          
          if (method == "BPS") {
            result <- pdmp_sample_from_stanmodel(
              model_path,
              data_path,
              flow = "PreconditionedBPS",
              algorithm = "GridThinningStrategy",
              T = 2000,
              grid_n = 30,
              t_warmup = 1000,
              sticky = FALSE,
              adaptive_scheme = "full",
              seed = seed,
              x0 = x0_fixed,
              support_boundary = support_boundary_control(
                mode = "line_search_truncated_refresh",
                max_refresh_attempts = 100L,
                refresh_probe_time = 1e-6,
                min_safe_time = 1e-8
                )
            )

            saveRDS(result, file_out)
            cat("Saved:", file_out, "\n")
          }

          
          if (method == "Boomerang") {
            result <- pdmp_sample_from_stanmodel(
              model_path,
              data_path,
              flow = "AdaptiveBoomerang",
              algorithm = "GridThinningStrategy",
              T = 2000,
              grid_n = 30,
              t_warmup = 1000,
              sticky = FALSE,
              adaptive_scheme = "full",
              seed = seed,
              x0 = x0_fixed,
              support_boundary = support_boundary_control(
                mode = "line_search_truncated_refresh",
                max_refresh_attempts = 100L,
                refresh_probe_time = 1e-6,
                min_safe_time = 1e-8
              )
            )

            saveRDS(result, file_out)
            cat("Saved:", file_out, "\n")
          }
          
          
          # if successful, remove old error log if it exists
          if (file.exists(error_out)) {
            file.remove(error_out)
          }
          
        }, error = function(e) {
          
          cat("ERROR:", e$message, "\n")
          
          writeLines(
            c(
              paste("time:", Sys.time()),
              paste("group:", group_name),
              paste("posterior:", post_name),
              paste("method:", method),
              paste("seed:", seed),
              paste("model_path:", model_path),
              paste("data_path:", data_path),
              paste("file_out:", file_out),
              paste("file_draws:", file_draws),
              paste("error:", e$message)
            ),
            error_out
          )
          
          cat("Saved error log:", error_out, "\n")
        })
      }
    }
  }
}




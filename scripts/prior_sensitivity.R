### 
# Systematically determine sensitivity to prior 
library(tidyverse)
library(here)
library(furrr)

source("scripts/likelihood_functions.R")
source("scripts/run_pipeline.R")

plan(multisession(workers = 8))

theme_set(theme_classic())

#####
# Compare priors
geom_params <- expand_grid(ns = 1:100, param = seq(0.1,1,by= 0.1)) %>% 
  group_by(param) %>% 
  mutate(p = dgeom(ns, param)/sum(dgeom(1:100, param)))

geom_params %>% 
  # filter(param == 0.4 | param == 0.5) %>% 
  ggplot(aes(ns, p, color = paste0("geom(", param, ")"), group = param)) + 
  geom_line() + 
  geom_line(aes(ns, dpois(ns, 1), color = "poisson(1)"), color = "black") + 
  coord_cartesian(c(0,15)) + 
  labs(color = "prior for n", caption = "poisson(1) in black for reference", 
       x = "n_k", y= "prior probability") 

#####
# Simulate data similar to fixed 8TR conditions in the case when there are likely 2 episomes per cluster
set.seed(400)
real_mu <- 525 # based on observed data
real_sigma2 <- 140^2 # based on observed data
q <- 230
real_ns <- sample(1:5, q, replace = T, prob = dpois(1:5, 2)/sum(dpois(1:5, 2))) 
hist(real_ns)
I <- matrix(rnorm(q, real_mu*real_ns, sqrt(real_sigma2*real_ns)), ncol = q)
names(I) <- paste0("n", 1:q)

# save(real_ns, I, file =  here('results', 'trouble_shooting', 'simulated_data_clustered.RData'))
# load(here('results', 'trouble_shooting', 'simulated_data_clustered.RData'))


# modified_gibbs <- function(tau0, mu0, n_prior, n_param, I, n_iterations){
#   # print(n_prior)
#   n_prior <- list(n_prior, n_param)
#   df <- run_gibbs(tau0, mu0, I, n_iterations, n_prior = n_prior)
#   print(str_interp("done- tau0: ${tau0} mu0: ${mu0}, ${n_prior}, ${n_param}"))
#   save(df, file = here("results", "trouble_shooting", str_interp("MCMC_${n_prior}_${n_param}_tau0${tau_0}_mu0${mu0}.RData")))
#   return(df)
# }

n_iterations <- 10000
hyper_parameters <- expand_grid(tau0 = c(2/real_sigma2, 0.5/real_sigma2), 
                                mu0 = c(300, 1000)) %>% 
  mutate(chain = 1:nrow(.))

# hyper_parameters <- rbind(
#   expand_grid(hyper_parameters, n_prior = "geom", n_param = c(0.1, 0.2, 0.3, 0.4, 0.5)),
#   hyper_parameters %>% mutate(n_prior = "pois", n_param = 1)
# ) %>% 
#   arrange(desc(n_prior), n_param)
# 
# all_chains <- hyper_parameters %>% head %>% 
#   select(-chain) %>% 
#   future_pmap_dfr(modified_gibbs, I = I, n_iterations = n_iterations,
#                  .options = furrr_options(seed = T), .id = "run")

cat('starting')
baseline <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("pois", 1), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done with baseline\n')

geom0.1 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.1), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.2 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.2), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.3 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.3), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.4 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.4), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.5 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.5), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.6 <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I,  n_prior = list("geom", 0.6), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

# save(baseline, geom0.1, geom0.2, geom0.3, geom0.4, geom0.5, geom0.6, file = here('results', 'trouble_shooting', 'chains_for_clustered_sims.RData'))
# load(here('results', 'trouble_shooting', 'chains_for_clustered_sims.RData'))

all_runs <- rbind(
  baseline %>% mutate(prior = "poisson(1)"),
  geom0.1 %>% mutate(prior = "geom(0.1)"),
  geom0.2 %>% mutate(prior = "geom(0.2)"),
  geom0.3 %>% mutate(prior = "geom(0.3)"),
  geom0.4 %>% mutate(prior = "geom(0.4)"),
  geom0.5 %>% mutate(prior = "geom(0.5)"),
  geom0.6 %>% mutate(prior = "geom(0.6)")
) %>% 
  filter(iteration > 1000) 

inferred_param <- all_runs %>% group_by(prior) %>% 
  summarise(mu = DescTools::Mode(round(mu)), sigma = DescTools::Mode(round(sqrt(1/tau))))

all_runs %>% 
  ggplot(aes(iteration, mu, color = chain)) + 
  geom_line() + 
  facet_wrap(~prior) + 
  geom_hline(yintercept = real_mu) + 
  geom_text(data = data.frame(NA), aes(1100, real_mu, label = real_mu), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  geom_hline(data = inferred_param, aes(yintercept = mu), lty = "dashed") +
  geom_text(data = inferred_param, aes(1100, mu, label = mu), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  ylim(c(0, max(all_runs$mu))) + 
  labs(y = "mu")

all_runs %>% 
  ggplot(aes(iteration, 1/tau, color = chain)) + 
  geom_line() + 
  facet_wrap(~prior) + 
  geom_hline(yintercept = real_sigma2) + 
  ylim(c(0, max(1/all_runs$tau))) + 
  labs(y = "sigma^2")

all_runs %>% 
  ggplot(aes(iteration, sqrt(1/tau), color = chain)) + 
  geom_line() + 
  facet_wrap(~prior) + 
  geom_hline(yintercept = sqrt(real_sigma2)) + 
  geom_text(data = data.frame(NA), aes(1100, sqrt(real_sigma2), label = sqrt(real_sigma2)), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  geom_hline(data = inferred_param, aes(yintercept = sigma), lty = "dashed") +
  geom_text(data = inferred_param, aes(1100, sigma, label = sigma), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  ylim(c(0, max(sqrt(1/all_runs$tau)))) + 
  labs(y = "sigma")


convergence <- all_runs %>% 
  group_split(prior) %>% 
  future_map_dfr(function(x) convergence_results(select(x, -prior)) %>% mutate(prior = unique(x$prior)))

convergence %>% 
  filter(prior == "geom(0.1)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  labs(y = "parameter") + 
  
  convergence %>% 
  filter(prior == "geom(0.2)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence %>% 
  filter(prior == "geom(0.3)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence %>% 
  filter(prior == "geom(0.4)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  labs(y = "parameter") + 
  
  convergence %>% 
  filter(prior == "geom(0.5)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence %>% 
  filter(prior == "poisson(1)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  plot_layout(ncol = 3, byrow = T, guides = "collect") &
  geom_point(color="gray") &
  geom_point(data = . %>% filter(name %in% c("mu", "tau")) ,
             aes(color = as.character(name))) &
  facet_wrap(~prior, scales = "free_y") & 
  geom_vline(xintercept = 1, lty = "dashed") &
  xlim(c(0.99, 1.05)) &
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) &
  labs(color = "")

## Compare to simulated datas
cluster_modes <- all_runs %>%
  filter(chain == "1") %>% 
  select(-mu, -tau, -iteration, -chain) %>% 
  pivot_longer(starts_with("n"), values_to = "n_epi", names_to = "cluster_id") %>% 
  count(prior, cluster_id, n_epi) %>% 
  group_by(prior, cluster_id) %>% 
  filter(n == max(n)) %>% 
  select(-n) %>% 
  ungroup() %>% 
  left_join(data.frame(real_ns) %>% mutate(cluster_id = paste0("n", 1:nrow(.)))) 

cluster_modes %>% 
  count(prior, n_epi, real_ns) %>% 
  ggplot(aes(real_ns, n_epi, fill = n)) + 
  geom_tile() + 
  facet_wrap(~prior) + 
  labs(x = "real nk", y= "estimated nk") + 
  geom_abline(lty = "dashed")

cluster_modes %>% 
  mutate(correct = 
           factor(ifelse(n_epi == real_ns, "correct inference", "incorrect inference"),
                  levels = c("incorrect inference", "correct inference"))) %>% 
  ggplot(aes(n_epi, fill = correct)) + 
  geom_bar(position = "stack") +
  facet_wrap(~prior) + 
  geom_bar(aes(real_ns, color = "simulated data"), fill = NA, inherit.aes = F) + 
  scale_color_manual(values = c("black") )+ 
  scale_fill_manual(values = c("red", "lightgreen")) + 
  labs(color = "", fill = "") + 
  labs(x = "number of episomes per cluster", y = "count") + 
  geom_text(data = .  %>% count(prior, correct) %>% group_by(prior) %>% mutate(p = n/sum(n)) %>% 
              filter(correct == "correct inference"),
            aes(Inf, Inf, label = paste0("accuracy: ", round(p*100), "%")),
            hjust = 1, vjust = 1.1, inherit.aes = F, show.legend = F) 

sds <- all_runs %>% 
  filter(chain == "1") %>% 
  pivot_longer(starts_with("n"), names_to = "cluster_id") %>% 
  group_by(prior, cluster_id) %>% 
  summarise(sd = sd(value))

tibble(intensity = I[1,], cluster_id = names(I)) %>% 
  left_join(cluster_modes) %>% 
  left_join(sds) %>% 
  ggplot(aes(intensity, as.factor(n_epi))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = sd)) +
  facet_wrap(~prior, scales = "free_y") +
  scale_color_viridis_c() + 
  labs(x = "Simulated Intensity Data", 
       y = "Estimated number of episomes in cluster",
       color = "Standard deviation\nof the\nposterior distribution")



#####
# Simulate another dataset with larger variance
set.seed(400)
real_sigma2 <- 350^2 # based on observed data
I_lv <- matrix(rnorm(q, real_mu*real_ns, sqrt(real_sigma2*real_ns)), ncol = q)
names(I_lv) <- paste0("n", 1:q)

# save(real_ns, I_lv, file =  here('results', 'trouble_shooting', 'simulated_data_clustered_larger_variance.RData'))

cat('starting\n')
baseline_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("pois", 1), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done with baseline\n')

geom0.1_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.1), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.2_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.2), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.3_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.3), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.4_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.4), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.5_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.5), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

geom0.6_lv <- hyper_parameters %>%
  select(-chain) %>%
  future_pmap_dfr(run_gibbs, I_lv,  n_prior = list("geom", 0.6), n_iterations,
                  .options = furrr_options(seed = T), .id = "chain")

cat('done\n')

# save(baseline_lv, geom0.1_lv, geom0.2_lv, geom0.3_lv, geom0.4_lv, geom0.5_lv, geom0.6_lv, file = here('results', 'trouble_shooting', 'chains_for_clustered_sims_larger_variance.RData'))

all_runs_lv <- rbind(
  baseline_lv %>% mutate(prior = "poisson(1)"),
  geom0.1_lv %>% mutate(prior = "geom(0.1)"),
  geom0.2_lv %>% mutate(prior = "geom(0.2)"),
  geom0.3_lv %>% mutate(prior = "geom(0.3)"),
  geom0.4_lv %>% mutate(prior = "geom(0.4)"),
  geom0.5_lv %>% mutate(prior = "geom(0.5)"),
  geom0.6_lv %>% mutate(prior = "geom(0.6)")
) %>% 
  filter(iteration > 1000) 

inferred_param_lv <- all_runs_lv %>% group_by(prior) %>% 
  summarise(mu = DescTools::Mode(round(mu)), sigma = DescTools::Mode(round(sqrt(1/tau))))

all_runs_lv %>% 
  ggplot(aes(iteration, mu, color = chain)) + 
  geom_line() + 
  facet_wrap(~prior) + 
  geom_hline(yintercept = real_mu) + 
  geom_text(data = data.frame(NA), aes(1100, real_mu, label = real_mu), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  geom_hline(data = inferred_param_lv, aes(yintercept = mu), lty = "dashed") +
  geom_text(data = inferred_param_lv, aes(1100, mu, label = mu), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  ylim(c(0, max(all_runs_lv$mu))) + 
  labs(y = "mu")

all_runs_lv %>% 
  ggplot(aes(iteration, 1/tau, color = chain)) + 
  geom_line() + 
  facet_wrap(~prior) + 
  geom_hline(yintercept = real_sigma2) + 
  ylim(c(0, max(1/all_runs$tau))) + 
  labs(y = "sigma^2")

all_runs_lv %>% 
  ggplot(aes(iteration, sqrt(1/tau), color = chain)) + 
  geom_line() + 
  facet_wrap(~prior) + 
  geom_hline(yintercept = sqrt(real_sigma2)) + 
  geom_text(data = data.frame(NA), aes(1100, sqrt(real_sigma2), label = sqrt(real_sigma2)), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  geom_hline(data = inferred_param_lv, aes(yintercept = sigma), lty = "dashed") +
  geom_text(data = inferred_param_lv, aes(1100, sigma, label = sigma), 
            hjust = 0, vjust = -0.1, inherit.aes = F) +
  ylim(c(0, max(sqrt(1/all_runs_lv$tau))+10)) + 
  labs(y = "sigma")


cluster_modes_lv <- all_runs_lv %>%
  filter(chain == "1") %>% 
  select(-mu, -tau, -iteration, -chain) %>% 
  pivot_longer(starts_with("n"), values_to = "n_epi", names_to = "cluster_id") %>% 
  count(prior, cluster_id, n_epi) %>% 
  group_by(prior, cluster_id) %>% 
  filter(n == max(n)) %>% 
  select(-n) %>% 
  ungroup() %>% 
  left_join(data.frame(real_ns) %>% mutate(cluster_id = paste0("n", 1:nrow(.)))) 

cluster_modes_lv %>% 
  count(prior, n_epi, real_ns) %>% 
  ggplot(aes(real_ns, n_epi, fill = n)) + 
  geom_tile() + 
  facet_wrap(~prior) + 
  labs(x = "real nk", y= "estimated nk") + 
  geom_abline(lty = "dashed")

cluster_modes_lv %>% 
  mutate(correct = 
           factor(ifelse(n_epi == real_ns, "correct inference", "incorrect inference"),
                  levels = c("incorrect inference", "correct inference"))) %>% 
  ggplot(aes(n_epi, fill = correct)) + 
  geom_bar(position = "stack") +
  facet_wrap(~prior) + 
  geom_bar(aes(real_ns, color = "simulated data"), fill = NA, inherit.aes = F) + 
  scale_color_manual(values = c("black") )+ 
  scale_fill_manual(values = c("red", "lightgreen")) + 
  labs(color = "", fill = "") + 
  labs(x = "number of episomes per cluster", y = "count") + 
  geom_text(data = .  %>% count(prior, correct) %>% group_by(prior) %>% mutate(p = n/sum(n)) %>% 
              filter(correct == "correct inference"),
            aes(Inf, Inf, label = paste0("accuracy: ", round(p*100), "%")),
            hjust = 1, vjust = 1.1, inherit.aes = F, show.legend = F) 

convergence_lv <- all_runs_lv %>% 
  group_split(prior) %>% 
  future_map_dfr(function(x) convergence_results(select(x, -prior)) %>% mutate(prior = unique(x$prior)))

convergence_lv %>% 
  filter(prior == "geom(0.1)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  labs(y = "parameter") + 
  
  convergence_lv %>% 
  filter(prior == "geom(0.2)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence_lv %>% 
  filter(prior == "geom(0.3)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence_lv %>% 
  filter(prior == "geom(0.4)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  labs(y = "parameter") + 
  
  convergence_lv %>% 
  filter(prior == "geom(0.5)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  convergence_lv %>% 
  filter(prior == "poisson(1)") %>% 
  arrange(Rhat) %>% 
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(Rhat, name)) + 
  theme(axis.title.y = element_blank()) +
  
  plot_layout(ncol = 3, byrow = T, guides = "collect") &
  geom_point(color="gray") &
  geom_point(data = . %>% filter(name %in% c("mu", "tau")) ,
             aes(color = as.character(name))) &
  facet_wrap(~prior, scales = "free_y") & 
  geom_vline(xintercept = 1, lty = "dashed") &
  xlim(c(0.99, 1.05)) &
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) &
  labs(color = "")


sds_lv <- all_runs_lv %>% 
  filter(chain == "1") %>% 
  pivot_longer(starts_with("n"), names_to = "cluster_id") %>% 
  group_by(prior, cluster_id) %>% 
  summarise(sd = sd(value))

tibble(intensity = I_lv[1,], cluster_id = names(I_lv)) %>% 
  left_join(cluster_modes_lv) %>% 
  left_join(sds_lv) %>% 
  ggplot(aes(intensity, as.factor(n_epi))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = sd)) +
  facet_wrap(~prior, scales = "free_y") +
  scale_color_viridis_c() + 
  labs(x = "Simulated Intensity Data", 
       y = "Estimated number of episomes in cluster",
       color = "Standard deviation\nof the\nposterior distribution")


rbind(
tibble(intensity = I[1,], cluster_id = names(I)) %>%
  left_join(cluster_modes) %>% 
  left_join(sds) %>% mutate(data = "baseline simulation"),

tibble(intensity = I_lv[1,], cluster_id = names(I_lv)) %>% 
  left_join(cluster_modes_lv) %>% 
  left_join(sds_lv) %>% mutate(data = "simulation with larger variance")
) %>% 
  filter(prior == "geom(0.5)" | prior == "geom0.5") %>% 
  ggplot(aes(intensity, as.factor(n_epi))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = sd)) +
  facet_wrap(~data, ncol = 1) +
  scale_color_viridis_c() + 
  labs(x = "Simulated Intensity Data", 
       y = "Estimated number of episomes in cluster",
       color = "Standard deviation\nof the\nposterior distribution") +
  theme_bw()

rbind(
  tibble(intensity = I[1,], cluster_id = names(I)) %>%
    left_join(cluster_modes) %>% 
    left_join(sds) %>% mutate(data = "baseline simulation") %>% 
    filter(prior == "geom0.5") %>% 
    left_join(inferred_param),
  
  tibble(intensity = I_lv[1,], cluster_id = names(I_lv)) %>% 
    left_join(cluster_modes_lv) %>% 
    left_join(sds_lv) %>% mutate(data = "simulation with larger variance") %>% 
    filter(prior == "geom(0.5)") %>% 
    left_join(inferred_param_lv)
) %>% 
  mutate(intensity = intensity/mu) %>% 
  ggplot(aes(intensity, as.factor(n_epi))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = sd)) +
  facet_wrap(~data, ncol = 1) +
  scale_color_viridis_c() + 
  labs(x = "Simulated Intensity Data Normalized by Inferred Mean", 
       y = "Estimated number of episomes in cluster",
       color = "Standard deviation\nof the\nposterior distribution") +
  scale_x_continuous(breaks = 0:10) + 
  theme_bw()
  

#####
# Make silhouette plots
make_sil <- function(I, n_prior, modes){
  dist_matrix <- dist(t(I))  
  class <- modes %>% 
    filter(prior == n_prior) %>% 
    group_by(cluster_id) %>% 
    filter(n_epi == min(n_epi)) %>% 
    ungroup %>% 
    mutate(cluster = as.numeric(gsub("n", "", cluster_id))) %>% 
    arrange(cluster) %>% 
    pull(n_epi)
  sil <- silhouette(class, dist_matrix)
  fviz_silhouette(sil)
}

make_sil(I, "geom(0.1)", cluster_modes)
make_sil(I, "geom(0.2)", cluster_modes)
make_sil(I, "geom(0.3)", cluster_modes)
make_sil(I, "geom(0.4)", cluster_modes)
make_sil(I, "geom(0.5)", cluster_modes) + 
  coord_flip() + theme(axis.text.y = element_blank(), axis.text.x = element_text())
make_sil(I, "geom(0.6)", cluster_modes)
make_sil(I, "poisson(1)", cluster_modes)


make_sil(I_lv, "geom(0.1)", cluster_modes_lv)
make_sil(I_lv, "geom(0.2)", cluster_modes_lv)
make_sil(I_lv, "geom(0.3)", cluster_modes_lv)
make_sil(I_lv, "geom(0.4)", cluster_modes_lv)
make_sil(I_lv, "geom(0.5)", cluster_modes_lv)
make_sil(I_lv, "geom(0.6)", cluster_modes_lv)
make_sil(I_lv, "poisson(1)", cluster_modes_lv)

hist(real_ns, freq = F)
points(1:5, dgeom(1:5, 0.5)/sum(dgeom(1:100, 0.5)), col = "red")
points(1:5, dgeom(1:5, 0.6)/sum(dgeom(1:100, 0.6)), col = "blue")
points(1:5, dgeom(1:5, 0.4)/sum(dgeom(1:100, 0.4)), col = "green")

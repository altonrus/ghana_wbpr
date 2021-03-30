library(data.table)
library(ggplot2)
library(scales)

##
#
# Calibrating annual probability of AIDS
#
##



states_hiv <- c("noART_1", "noART_2", "noART_3",
                "ART_1", "ART_2", "ART_3", "ART_4", "AIDS", "HIV_death")

t_hiv <- matrix(0L, nrow = length(states_hiv), ncol = length(states_hiv))
rownames(t_hiv) <- states_hiv
colnames(t_hiv) <- states_hiv

t_hiv["noART_1", "noART_2"] <- 1
t_hiv["noART_2", "noART_3"] <- 1
t_hiv["noART_3", "noART_3"] <- 1
t_hiv["ART_1", "ART_2"] <- 1
t_hiv["ART_2", "ART_3"] <- 1
t_hiv["ART_3", "ART_4"] <- 1
t_hiv["ART_4", "ART_4"] <- 1
t_hiv["HIV_death", "HIV_death"] <- 1
# t_hiv["noART_4", "AIDS"] <- 0.12
# t_hiv["noART_4", "noART_4"] <- 1- t_hiv["noART_4", "AIDS"]
# t_hiv["AIDS", "HIV_death"] <- 0
t_hiv["AIDS", "AIDS"] <- 1-t_hiv["AIDS", "HIV_death"]

rowSums(t_hiv)



## HIV time to AIDS without death or all-cause mortality


sim_calib_pAIDS <- function(init, t_matrix, iters = 10){
  #t matrix is missing death, we add in at each step
  #i_age <- age_start
  
  markov_trace <- matrix( nrow = iters, ncol = nrow(t_matrix))
  colnames(markov_trace) <- colnames(t_matrix)
  
  markov_trace[1, ] <- init #%*% t_matrix
  
  for (i in 2:iters){
    
    markov_trace[i, ] <- markov_trace[i-1 , ] %*% t_matrix
  }
  return(markov_trace)
}

param_HIV_no_death <- function(t_hiv, p_AIDS){
  t_hiv["noART_3", "AIDS"] <- p_AIDS
  t_hiv["noART_3", "noART_3"] <- 1- p_AIDS
  return(t_hiv)
}


dt_targets_AIDS <- data.table(years = c(5, 7, 9),
                         peds_BC = c(0.123297855440826, 0.201760127084988, 0.252200158856235),
                         peds_LB = c(0.089671167593328, 0.151320095313741, 0.190551231135822),
                         peds_UB = c(0.173737887212073, 0.257804606830818, 0.319453534551231),
                         adult_BC = c(0.344875637425138, 0.564341952150226, 0.705427440187783),
                         adult_LB = c(0.250818645400101, 0.42325646411267, 0.532989621475214),
                         adult_UB = c(0.485961125462695, 0.721103605525289, 0.893541424237858))

dt_params <- data.table(
  name = c("p_AIDS", "p_AIDS_death", "p_HIV_death"),
  peds_BC = 0,
  peds_LB = 0,
  peds_UB = 0,
  adult_BC = 0,
  adult_LB = 0,
  adult_UB = 0
)
p_AIDS_to_try <- seq(0, .5, by = 0.0001)


for (param_version in colnames(dt_targets_AIDS)[2:7]){
  targets <- dt_targets_AIDS[, get(param_version) ]
  
  SSE_min <- 1

  init <- c(1, rep(0, nrow(t_hiv)-1))
  
  for (p_AIDS in p_AIDS_to_try){
        t_hiv["noART_3", "AIDS"] <- p_AIDS
        t_hiv["noART_3", "noART_3"] <- 1- p_AIDS
        
        
        #sum of squared error
        SSE <- sum((sim_calib_pAIDS(init, t_hiv)[c(6,8,10), "AIDS"] - targets)^2)
        
        # print(pAIDS)
        # print(SSE)
        
        if (SSE < SSE_min){
          SSE_min <- SSE
          p_AIDS_best <- p_AIDS
        }
  }
  dt_params[name == "p_AIDS" , c(eval(param_version)) := p_AIDS_best]
}


dt_params


init <- c(1, rep(0, nrow(t_hiv)-1))

AIDS_calib_peds <- data.table(
  age = 5:95,
  BC = sim_calib_pAIDS(init, 
                             param_HIV_no_death(t_hiv, dt_params[name == "p_AIDS", peds_BC]),iters= 95-5+1)[ , "AIDS"],
  LB = sim_calib_pAIDS(init, 
                            param_HIV_no_death(t_hiv, dt_params[name == "p_AIDS", peds_LB]),iters= 95-5+1)[ , "AIDS"],
  UB = sim_calib_pAIDS(init, 
                            param_HIV_no_death(t_hiv, dt_params[name == "p_AIDS", peds_UB]),iters= 95-5+1)[ , "AIDS"]
  
  
)



AIDS_calib_adult <- data.table(
  age = 40:95,
  BC = sim_calib_pAIDS(init, 
                             param_HIV_no_death(t_hiv, dt_params[name == "p_AIDS", adult_BC]),iters= 95-40+1)[ , "AIDS"],
  LB = sim_calib_pAIDS(init, 
                       param_HIV_no_death(t_hiv, dt_params[name == "p_AIDS", adult_LB]),iters= 95-40+1)[ , "AIDS"],
  UB = sim_calib_pAIDS(init, 
                       param_HIV_no_death(t_hiv, dt_params[name == "p_AIDS", adult_UB]),iters= 95-40+1)[ , "AIDS"]
  
)

AIDS_calib_all <- rbind(
  cbind(cohort = "Pediatric", melt(AIDS_calib_peds, id.vars = "age", variable.name = "Limit", value.name = "Proportion_with_AIDS")),
  cbind(cohort = "Adult", melt(AIDS_calib_adult, id.vars = "age", variable.name = "Limit", value.name = "Proportion_with_AIDS"))
)


ggplot()+
  geom_line(data = AIDS_calib_all, aes(x = age, y = Proportion_with_AIDS, color = cohort, linetype = Limit))+
  geom_pointrange(data = dt_targets_AIDS, aes(x = 5+years, y = peds_BC, ymin = peds_LB, ymax = peds_UB), color = "royalblue3")+
  geom_pointrange(data = dt_targets_AIDS, aes(x = 40+years, y = adult_BC, ymin = adult_LB, ymax = adult_UB), color = "red")+
  ylab("Proportion with AIDS")+
  xlab("Cohort age")+
  scale_y_continuous(labels = label_percent())





##
#
# Calibrating annual probability of death for each cohort
#
##

sim_markov <- function(init, t_matrix, age_start, all_cause_mort){
  #t matrix is missing death, we add in at each step
  
  #Simulate by year from age_start to age 95
  iters <- 95 - age_start + 1
  
  markov_trace <- matrix( nrow = iters, ncol = 1+nrow(t_matrix))
  colnames(markov_trace) <- c(colnames(t_matrix), "death")
  rownames(markov_trace) <- age_start:95
  
  # First row has initial condition
  markov_trace[as.character(age_start), ] <- init
  
  
  for (i in (age_start+1):95){
    
    #get probability of all-cause death
    p_death <- all_cause_mort[findInterval(i, all_cause_mort$age)]$p_death
    #Create temporary transition matrix that encorporates probability of all cause death 
    #  and down-weights all other transitions so it still sums to 1
    t_matrix_temp <- rbind(cbind(t_matrix*(1-p_death), death = p_death), death = c(rep(0, ncol(t_matrix)), 1))
    
    #Apply markovian transition
    markov_trace[as.character(i), ] <- markov_trace[as.character(i-1) , ] %*% t_matrix_temp
  }
  
  return(markov_trace)
}


all_cause_mort <- fread(".'results/ghana_allcause_mort.csv")



init <- c(1, rep(0, nrow(t_hiv)))
# p_AIDS_death_to_try <- seq(0, 1, by = 0.01)
# p_HIV_death_to_try <- seq(0, 0.1, by = 0.01)
p_AIDS_death_to_try <- seq(.2, .8, by = 0.001)
p_HIV_death_to_try <- seq(0, 0.1, by = 0.001)

dt_targets_survival <- data.table(
  years = c(4, 7),
  peds_BC = c(0.945454545454545, 0.79),
  peds_LB = c(0.8488, 0.63),
  peds_UB = c(0.9886, 0.88),
  adult_BC = c(0.902654867256637, 0.545132743362832),
  adult_LB = c(0.797607079646018, 0.391858407079646),
  adult_UB = c(0.964148672566372, 0.685398230088496)
)


for (param_version in colnames(dt_targets_AIDS)[2:7]){
  targets <- dt_targets_survival[, get(param_version) ]
  SSE_min <- 1
  
  init <- c(1, rep(0, nrow(t_hiv)))
  
  if (substr(param_version, 1, 4) == "peds"){
    p_AIDS <- dt_params[name == "p_AIDS", peds_BC]
    age_start <- 5
  } else {
    p_AIDS <- dt_params[name == "p_AIDS", adult_BC]
    age_start <- 40
  }
  t_hiv["noART_3", "AIDS"] <- p_AIDS
  t_hiv["noART_3", "noART_3"] <- 1- p_AIDS
  
  for (p_AIDS_death in p_AIDS_death_to_try){
    for (p_HIV_death in p_HIV_death_to_try){
      t_hiv["AIDS", "HIV_death"] <- p_AIDS_death
      t_hiv["AIDS", "AIDS"] <- 1 - p_AIDS_death
      t_hiv["noART_3", "HIV_death"] <- p_HIV_death
      t_hiv["noART_3", "noART_3"] <- 1- p_HIV_death - p_AIDS
      #sum of squared error
      SSE <- sum((rowSums(sim_markov(init, t_hiv, age_start, all_cause_mort)[c(5, 8), 1:8]) - targets)^2)
      # print(pAIDS)
      # print(SSE)
      if (SSE < SSE_min){
        p_HIV_death_best <- p_HIV_death
        p_AIDS_death_best <- p_AIDS_death
        SSE_min <- SSE
      }
    }
  }
  
  dt_params[name == "p_HIV_death" , c(eval(param_version)) := p_HIV_death_best]
  dt_params[name == "p_AIDS_death" , c(eval(param_version)) := p_AIDS_death_best]
}



SSE_min

dt_params
fwrite(dt_params, ".'results/HIV_params_calibrated.csv")

dt_params <- fread(".'results/HIV_params_calibrated.csv")

update_trans_matrix_death <- function(p_AIDS, p_AIDS_death, p_HIV_death){
  t_hiv["noART_3", "AIDS"] <- p_AIDS
  t_hiv["noART_3", "noART_3"] <- 1- p_AIDS
  t_hiv["AIDS", "HIV_death"] <- p_AIDS_death
  t_hiv["AIDS", "AIDS"] <- 1 - p_AIDS_death
  t_hiv["noART_3", "HIV_death"] <- p_HIV_death
  t_hiv["noART_3", "noART_3"] <- 1- p_HIV_death - p_AIDS
  return(t_hiv)
}














init <- c(1, rep(0, nrow(t_hiv)))

Death_calib_peds <- data.table(
  age = 5:95,
  BC = rowSums(sim_markov(init=init,
                          t_matrix = update_trans_matrix_death(dt_params[name == "p_AIDS", peds_BC],
                                                               dt_params[name == "p_AIDS_death", peds_BC],
                                                               dt_params[name == "p_HIV_death", peds_BC]), 
                          5, all_cause_mort)[, 1:8]),
  LB = rowSums(sim_markov(init,
                          update_trans_matrix_death(dt_params[name == "p_AIDS", peds_BC],
                                                    dt_params[name == "p_AIDS_death", peds_LB],
                                                    dt_params[name == "p_HIV_death", peds_LB]), 
                          5, all_cause_mort)[, 1:8]),
  UB = rowSums(sim_markov(init,
                          update_trans_matrix_death(dt_params[name == "p_AIDS", peds_BC],
                                                    dt_params[name == "p_AIDS_death", peds_UB],
                                                    dt_params[name == "p_HIV_death", peds_UB]), 
                          5, all_cause_mort)[, 1:8])
)



Death_calib_adult <- data.table(
  age = 40:95,
  BC = rowSums(sim_markov(init,
                          update_trans_matrix_death(dt_params[name == "p_AIDS", adult_BC],
                                                    dt_params[name == "p_AIDS_death", adult_BC],
                                                    dt_params[name == "p_HIV_death", adult_BC]), 
                          age_start = 40, all_cause_mort)[, 1:8]),
  LB = rowSums(sim_markov(init,
                          update_trans_matrix_death(dt_params[name == "p_AIDS", adult_BC],
                                                    dt_params[name == "p_AIDS_death", adult_LB],
                                                    dt_params[name == "p_HIV_death", adult_LB]), 
                          age_start = 40, all_cause_mort)[, 1:8]),
  UB = rowSums(sim_markov(init,
                          update_trans_matrix_death(dt_params[name == "p_AIDS", adult_BC],
                                                    dt_params[name == "p_AIDS_death", adult_UB],
                                                    dt_params[name == "p_HIV_death", adult_UB]), 
                          age_start = 40, all_cause_mort)[, 1:8])
)

Death_calib_all <- rbind(
  cbind(cohort = "Pediatric", melt(Death_calib_peds, id.vars = "age", variable.name = "Limit", value.name = "Proportion_alive")),
  cbind(cohort = "Adult", melt(Death_calib_adult, id.vars = "age", variable.name = "Limit", value.name = "Proportion_alive"))
)


ggplot()+
  geom_line(data = Death_calib_all, aes(x = age, y = Proportion_alive, color = cohort, linetype = Limit))+
  geom_pointrange(data = dt_targets_survival, aes(x = years+5, y = peds_BC, ymin = peds_LB, ymax = peds_UB), color = "royalblue3")+
  geom_pointrange(data = dt_targets_survival, aes(x = years+40, y = adult_BC, ymin = adult_LB, ymax = adult_UB), color = "red")+
  ylab("Proportion alive")+
  xlab("Cohort age")+
  scale_y_continuous(labels = label_percent())


Death_calib_all
AIDS_calib_all

calib_trace_combined <- cbind(Death_calib_all, "Proportion_with_AIDS" = AIDS_calib_all$Proportion_with_AIDS)
calib_trace_combined_melt <- melt(calib_trace_combined, measure.vars = c("Proportion_with_AIDS", "Proportion_alive"))

dt_targets_combined <- rbind(
  cbind(dt_targets_survival, "variable" = "Proportion_alive"),
  cbind(dt_targets_AIDS, "variable" = "Proportion_with_AIDS")
)

dt_targets_bc <- melt(dt_targets_combined, id.vars = c("years", "variable"), measure.vars =  c("peds_BC", "adult_BC"), value.name =  "Base case", variable.name = "cohort")
dt_targets_lb <- melt(dt_targets_combined, id.vars = c("years", "variable"), measure.vars =  c("peds_LB", "adult_LB"), value.name =  "Lower bound", variable.name = "cohort")
dt_targets_ub <- melt(dt_targets_combined, id.vars = c("years", "variable"), measure.vars =  c("peds_UB", "adult_UB"), value.name =  "Upper bound", variable.name = "cohort")
dt_targets_plot <- cbind(dt_targets_bc, "Lower bound" = dt_targets_lb$`Lower bound`, "Upper bound" = dt_targets_ub$`Upper bound`)
dt_targets_plot[ , cohort := ifelse(cohort == "peds_BC", "Pediatric", "Adult")]
dt_targets_plot[ , age := ifelse(cohort == "Pediatric", years+5, years+40)]

labs <- c("Proportion with AIDS (no ART, no deaths)",
          "Proportion alive (NO ART)")
names(labs)<- c("Proportion_with_AIDS", "Proportion_alive")
calib_trace_combined_melt[ , variable := factor(variable, levels = c("Proportion_with_AIDS", "Proportion_alive"))]
dt_targets_plot[ , variable := factor(variable, levels = c("Proportion_with_AIDS", "Proportion_alive"))]

ggplot()+
  geom_line(data = calib_trace_combined_melt, aes(x = age, y = value, color = cohort, linetype = Limit))+
  facet_grid(rows = vars(variable), 
             labeller = labeller(variable = labs))+
  geom_pointrange(data = dt_targets_plot, aes(x = age, y = `Base case`, ymin = `Lower bound`, ymax = `Upper bound`, color = cohort))+
  ylab("")+
  xlab("Cohort age")+
  scale_y_continuous(labels = label_percent())+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave(".'manuscript/figs/HIV_calib_plot.png",
       height = 7, width = 7, units = "in")


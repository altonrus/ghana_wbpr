library(data.table)
library(readxl)
library(ggplot2)
library(mc2d)
theme_set(theme_bw())

source("./Rscripts/model_functions.R")

## Data from file
all_cause_mort <- fread("./data/ghana_allcause_mort.csv")
dt_transitions <- data.table(read_excel("./data/Ghana_hea_parameters.xlsx", sheet = "Markov_probs"))
# transitions <- cbind(dt_transitions)
# transitions[ , value := `Base case`]
dt_microcost_params <- data.table(read_excel("./data/Ghana_hea_parameters.xlsx", sheet = "Microcost_params"))
dt_riskmod_params <- data.table(read_excel("./data/Ghana_hea_parameters.xlsx", sheet = "riskmod_params"))
dt_daly_params <- data.table(read_excel("./data/Ghana_hea_parameters.xlsx", sheet = "daly_params"))

# # # # #
# RUN BASECASE #################

#FORMAT DATA FOR MODEL
dt_transitions[ , value := `Base case`]
bc_microcost_params <- as.list(dt_microcost_params$`Base case`)
names(bc_microcost_params) <- dt_microcost_params$rname
#Get risk model params
params <- as.list(unlist(dt_riskmod_params[ , `Base case`]))
names(params) <- dt_riskmod_params$rname 
#get DALY params in list
bc_daly_params <- as.list(dt_daly_params[Disease=="acute"]$`Base case`)
names(bc_daly_params) <- dt_daly_params[Disease=="acute"]$rname
dw_hiv <- as.list(dt_daly_params[Disease=="hiv"]$`Base case`)
names(dw_hiv) <- dt_daly_params[Disease=="hiv"]$rname
dw_hbv <- as.list(dt_daly_params[Disease=="hbv"]$`Base case`)
names(dw_hbv) <- dt_daly_params[Disease=="hbv"]$rname
dw_hcv <- as.list(dt_daly_params[Disease=="hcv"]$`Base case`)
names(dw_hcv) <- dt_daly_params[Disease=="hcv"]$rname
bc_daly_params$hiv<- dw_hiv
bc_daly_params$hbv<- dw_hbv
bc_daly_params$hcv<- dw_hcv


#RUN MODEL
dalys_np <- calc_acute_dalys(bc_microcost_params, bc_daly_params, params)
costs_states <- calc_disease_state_costs(bc_microcost_params)
#Get outcomes for chronic viruses from running Markov models
hiv_outcomes <- hiv_calc_np_outcomes(costs_states, dt_transitions, params, 
                                            bc_daly_params, bc_microcost_params$p_IPmort)
hcv_outcomes <- hcv_calc_np_outcomes(costs_states, dt_transitions, params, 
                                            bc_daly_params, bc_microcost_params$p_IPmort)
hbv_outcomes <- hbv_calc_np_outcomes(costs_states, dt_transitions, params, 
                                            bc_daly_params, bc_microcost_params$p_IPmort)
#Format costs_np and dalys_np for outcome calculation
ae_acute <- c("sep", "syp", "mal", "ftr")
costs_np <- costs_states[ae_acute]
costs_np$hiv <- hiv_outcomes$np_cost
costs_np$hbv <- hcv_outcomes$np_cost
costs_np$hcv <- hbv_outcomes$np_cost
dalys_np$YLL$hiv <- hiv_outcomes$np_YLL
dalys_np$YLD$hiv <- hiv_outcomes$np_YLD
dalys_np$YLL$hbv <- hbv_outcomes$np_YLL
dalys_np$YLD$hbv <- hbv_outcomes$np_YLD
dalys_np$YLL$hcv <- hcv_outcomes$np_YLL
dalys_np$YLD$hcv <- hcv_outcomes$np_YLD
#Calc outcomes
dt_bc_outcomes <- calc_outcomes(costs_np, dalys_np, params)

dt_bc_costs <- rbind(unlist(costs_np))
dt_bc_dalys <- rbind(unlist(dalys_np))

## Scenario analysis
#Secondary infections
dt_bc_outcomes[ , sec_infections_burden_red := (
  sep.burden_reduced+mal.burden_reduced+ftr.burden_reduced+syp.burden_reduced+
    2*hiv.burden_reduced+2*hbv.burden_reduced+2*hcv.burden_reduced
)]
dt_bc_outcomes[ , sec_infections_net_cost := (
  prt_cost - sec_infections_burden_red
)]
dt_bc_outcomes[ , sec_infections_DALY_reduced := (
  sep.DALY_reduced+mal.DALY_reduced+ftr.burden_reduced+syp.burden_reduced+
    2*hiv.burden_reduced+2*hbv.burden_reduced+2*hcv.burden_reduced
)]

#No TTBI-sepsis
dt_bc_outcomes[ , no_sepsis_burden_red := (
  mal.burden_reduced+ftr.burden_reduced+syp.burden_reduced+
    hiv.burden_reduced+hbv.burden_reduced+hcv.burden_reduced
)]
dt_bc_outcomes[ , no_sepsis_net_cost := (
  prt_cost - no_sepsis_burden_red
)]
dt_bc_outcomes[ , no_sepsis_DALY_red := (
  mal.DALY_reduced+ftr.DALY_reduced+syp.DALY_reduced+
    hiv.DALY_reduced+hbv.DALY_reduced+hcv.DALY_reduced
)]
dt_bc_outcomes[ , no_sepsis_icer := (no_sepsis_net_cost)/no_sepsis_DALY_red]

fwrite(dt_bc_costs, "./results/BC_costs.csv")
fwrite(dt_bc_dalys, "./results/BC_dalys.csv")
fwrite(dt_bc_outcomes, "./results/BC_outcomes.csv")

# # # # #
# RUN PSA #################

n_iter = 10000
ae_acute <- c("sep", "syp", "mal", "ftr")

for(iter in 1:n_iter){
  
  #Generate microcost inputs from PERT distribution
  microcost_params <- as.list(unlist(dt_microcost_params[ , sample_params(mean = `Base case`, lb = Low, ub = High, dist = Distn, Param1=Param1, Param2=Param2), by=rname]$V1))
  names(microcost_params) <- dt_microcost_params$rname
  #Generate risk model params from indicated distribution
  params <- as.list(unlist(dt_riskmod_params[ , sample_params(mean = `Base case`, lb = Low, ub = High, dist = Distn, Param1=Param1, Param2=Param2), by=rname]$V1))
  names(params) <- dt_riskmod_params$rname 
  #Generate Markov prob inputs drawn from distributions
  markov_probs_not_varied <- dt_transitions[is.na(Low), ]
  markov_probs_not_varied[ , value := `Base case`]
  markov_probs_varied <- dt_transitions[!is.na(Low), ]
  markov_probs_varied[, `Base case` := as.numeric(`Base case`)]
  markov_probs_varied[, value := sample_params(mean = `Base case`, lb = Low, ub = High, dist = Distn, Param1=Param1, Param2=Param2), by=1:nrow(markov_probs_varied)]
  markov_probs <- rbind(markov_probs_not_varied, markov_probs_varied)
  #Generate DALY parameters from distributions
  daly_params <- as.list(unlist(dt_daly_params[Disease=="acute" , sample_params(mean = `Base case`, lb = Low, ub = High, dist = "PERT", Param1=Low, Param2=High), by=rname]$V1))
  names(daly_params) <- dt_daly_params[Disease=="acute"]$rname 
  dw_hiv <- as.list(unlist(dt_daly_params[Disease=="hiv" , sample_params(mean = `Base case`, lb = Low, ub = High, dist = "PERT", Param1=Low, Param2=High), by=rname]$V1))#
  names(dw_hiv) <- dt_daly_params[Disease=="hiv"]$rname
  dw_hbv <- as.list(unlist(dt_daly_params[Disease=="hbv" , sample_params(mean = `Base case`, lb = Low, ub = High, dist = "PERT", Param1=Low, Param2=High), by=rname]$V1))#
  names(dw_hbv) <- dt_daly_params[Disease=="hbv"]$rname
  dw_hcv <- as.list(unlist(dt_daly_params[Disease=="hcv" , sample_params(mean = `Base case`, lb = Low, ub = High, dist = "PERT", Param1=Low, Param2=High), by=rname]$V1))#
  names(dw_hcv) <- dt_daly_params[Disease=="hcv"]$rname
  daly_params$hiv<- dw_hiv
  daly_params$hbv<- dw_hbv
  daly_params$hcv<- dw_hcv
  
  
  #Calculate disease state costs and DALY outcomes for acute
  costs_states <- calc_disease_state_costs(microcost_params)
  dalys_np <- calc_acute_dalys(microcost_params, daly_params, params)
  
  #Get outcomes for chronic viruses from running Markov models
  hiv_outcomes <- hiv_calc_np_outcomes(costs_states, markov_probs, params, 
                                       daly_params, microcost_params$p_IPmort)
  hcv_outcomes <- hcv_calc_np_outcomes(costs_states, markov_probs, params, 
                                       daly_params, microcost_params$p_IPmort)
  hbv_outcomes <- hbv_calc_np_outcomes(costs_states, markov_probs, params, 
                                       daly_params, microcost_params$p_IPmort)
  
  
  #Format costs_np and dalys_np for outcome calculation
  costs_np <- costs_states[ae_acute]
  costs_np$hiv <- hiv_outcomes$np_cost
  costs_np$hbv <- hcv_outcomes$np_cost
  costs_np$hcv <- hbv_outcomes$np_cost
  dalys_np$YLL$hiv <- hiv_outcomes$np_YLL
  dalys_np$YLD$hiv <- hiv_outcomes$np_YLD
  dalys_np$YLL$hbv <- hbv_outcomes$np_YLL
  dalys_np$YLD$hbv <- hbv_outcomes$np_YLD
  dalys_np$YLL$hcv <- hcv_outcomes$np_YLL
  dalys_np$YLD$hcv <- hcv_outcomes$np_YLD
  
  dt_outcomes <- calc_outcomes(costs_np, dalys_np, params)
  
  
  #Add inputs and outputs to data.frames
  if(iter==1){
    dt_PSA_microcost <- as.data.table(microcost_params)
    dt_PSA_params <- as.data.table(params)
    dt_PSA_daly_params <- rbind(unlist(daly_params))
    dt_psa_costs <- rbind(unlist(costs_np))
    dt_psa_dalys <- rbind(unlist(dalys_np))
    dt_psa_outcomes <- dt_outcomes
    
    
  } else{
    dt_PSA_microcost <- rbind(dt_PSA_microcost,
                              as.data.table(microcost_params))
    dt_PSA_params <- rbind(dt_PSA_params,
                           as.data.table(params))
    dt_PSA_daly_params<-rbind(dt_PSA_daly_params,
                              rbind(unlist(daly_params)))
    dt_psa_costs <- rbind(dt_psa_costs,
                          rbind(unlist(costs_np)))
    dt_psa_dalys <- rbind(dt_psa_dalys,
                          rbind(unlist(dalys_np)))
    dt_psa_outcomes <- rbind(dt_psa_outcomes,
                             dt_outcomes)
  }
  if (iter %% 100 == 0) {print(paste0("Finished iter ", iter))}
}

##ADD SCENARIO ANALYSES TO OUTCOMES DATATABLES
## Sensitivity analysis
#Secondary infections
dt_psa_outcomes[ , sec_infections_burden_red := sep.burden_reduced+
                   mal.burden_reduced+ftr.burden_reduced+
                   syp.burden_reduced+2*hiv.burden_reduced+
                   2*hbv.burden_reduced+2*hcv.burden_reduced]
dt_psa_outcomes[ , sec_infections_net_cost := (
  prt_cost - sec_infections_burden_red
)]
dt_psa_outcomes[ , sec_infections_DALY_red := sep.burden_reduced+
                   mal.burden_reduced+ftr.burden_reduced+
                   syp.burden_reduced+2*hiv.burden_reduced+
                   2*hbv.burden_reduced+2*hcv.burden_reduced]

#No TTBI-sepsis
dt_psa_outcomes[ , no_sepsis_burden_red := (
  mal.burden_reduced+ftr.burden_reduced+syp.burden_reduced+
    hiv.burden_reduced+hbv.burden_reduced+hcv.burden_reduced
)]
dt_psa_outcomes[ , no_sepsis_net_cost := (
  prt_cost - no_sepsis_burden_red
)]
dt_psa_outcomes[ , no_sepsis_DALY_red := (
  mal.DALY_reduced+ftr.DALY_reduced+syp.DALY_reduced+
    hiv.DALY_reduced+hbv.DALY_reduced+hcv.DALY_reduced
)]
dt_psa_outcomes[ , no_sepsis_icer :=(no_sepsis_net_cost)/no_sepsis_DALY_red]




fwrite(dt_PSA_microcost, "./results/PSA_microcost.csv")
fwrite(dt_PSA_params, "./results/PSA_params.csv")
fwrite(dt_PSA_daly_params, "./results/PSA_daly_params.csv")
fwrite(dt_psa_costs, "./results/PSA_costs.csv")
fwrite(dt_psa_dalys, "./results/PSA_dalys.csv")
fwrite(dt_psa_outcomes, "./results/PSA_outcomes.csv")
#costs <- costs_states$HBV



# # # #
# Univariate sensitivity analysis ###########

#Form datatable with all parameters
dt_riskmod_params[ , Param_disp := ifelse(is.na(Disease), Parameter, paste0(Disease, " ", tolower(Parameter)))]
dt_microcost_params[ , Param_disp := ifelse(Disease=="General", paste0(Category,", ",Parameter), paste0(Disease, " ", tolower(Category),", ",Parameter))]
dt_transitions[ , rname := paste0("p_",Disease,"_",From,"_",To,"_",Cohort)]
dt_transitions[ , Param_disp := paste0(Disease, " ", From_disp," to ", To_disp)]
dt_transitions_params <- dt_transitions[!is.na(Low)]

dt_univ <- rbind(
  dt_riskmod_params[ , c("rname", "Param_disp", "Base case", "Low", "High")],
  dt_microcost_params[ , c("rname", "Param_disp", "Base case", "Low", "High")],
  dt_transitions_params[ , c("rname", "Param_disp", "Base case", "Low", "High")],
  dt_daly_params[, c("rname_long", "Description", "Base case", "Low", "High")],
  use.names=FALSE
)


dt_univ[, c("net_savings_low", "DALYs_averted_low", "nmb_low") :=  
          calc_univariate_change(rname,  Low, dt_transitions, 
                                 dt_microcost_params, 
                                 dt_riskmod_params,
                                 dt_daly_params),
        by = rname]
dt_univ[, c("net_savings_high", "DALYs_averted_high", "nmb_high") :=  
          calc_univariate_change(rname,  High, dt_transitions, 
                                 dt_microcost_params, 
                                 dt_riskmod_params,
                                 dt_daly_params),
        by = rname]
dt_univ[ , abs_net_savings_range := abs(net_savings_high - net_savings_low)]
dt_univ[ , abs_DALY_range := abs(DALYs_averted_high - DALYs_averted_low)]
dt_univ[ , abs_nmb_range := abs(nmb_high - nmb_low)]

dt_univ <- dt_univ[order(-abs_net_savings_range)]

fwrite(dt_univ, "./results/univ_sens_analysis.csv")






# # # # #
# Analysis of Markov models #############
# # # # #


#FORMAT DATA FOR MODEL
dt_transitions[ , value := `Base case`]
bc_microcost_params <- as.list(dt_microcost_params$`Base case`)
names(bc_microcost_params) <- dt_microcost_params$rname
#Get risk model params
params <- as.list(unlist(dt_riskmod_params[ , `Base case`]))
names(params) <- dt_riskmod_params$rname 
#get DALY params in list
bc_daly_params <- as.list(dt_daly_params[Disease=="acute"]$`Base case`)
names(bc_daly_params) <- dt_daly_params[Disease=="acute"]$rname
dw_hiv <- as.list(dt_daly_params[Disease=="hiv"]$`Base case`)
names(dw_hiv) <- dt_daly_params[Disease=="hiv"]$rname
dw_hbv <- as.list(dt_daly_params[Disease=="hbv"]$`Base case`)
names(dw_hbv) <- dt_daly_params[Disease=="hbv"]$rname
dw_hcv <- as.list(dt_daly_params[Disease=="hcv"]$`Base case`)
names(dw_hcv) <- dt_daly_params[Disease=="hcv"]$rname
bc_daly_params$hiv<- dw_hiv
bc_daly_params$hbv<- dw_hbv
bc_daly_params$hcv<- dw_hcv


#RUN MODEL
dalys_np <- calc_acute_dalys(bc_microcost_params, bc_daly_params, params)
costs_states <- calc_disease_state_costs(bc_microcost_params)

states_hiv <- c('HIV_SC_1', 'HIV_SC_2', 'HIV_SC_3', 'ART_1', 'ART_2', 'ART_3', 'ART_4', 
                'AIDS_SC', 'AIDS_ART_1', 'AIDS_ART_2', 'AIDS_RD', 'HIV_death')

#Create init vector with initial distribution
p_ART_init <- as.numeric(dt_transitions[Disease=="HIV" & From == "init" & To == "ART_1", value])
init_hiv <- rep(0, length(states_hiv)+1)
names(init_hiv) <- c(states_hiv,"death")
init_hiv["ART_1"] <- p_ART_init
init_hiv["HIV_SC_1"] <- 1-p_ART_init


t_hiv_peds <- gen_t_matrix(dt_transitions, "hiv", "Pediatric")
t_hiv_adult <- gen_t_matrix(dt_transitions, "hiv", "Adult")


HIV_peds <- sim_markov(init_hiv, t_hiv_peds,
                       age_start=5, 
                       all_cause_mort = all_cause_mort,
                       costs = costs_states$HIV,
                       DALY_weights = bc_daly_params$hiv
)
HIV_adult <- sim_markov(init_hiv, t_hiv_adult, 
                        age_start=40, 
                        all_cause_mort = all_cause_mort,
                        costs =  costs_states$HIV,
                        DALY_weights = bc_daly_params$hiv)

HIV_peds_trace <- data.table(HIV_peds$markov_trace, keep.rownames = TRUE)
HIV_peds_trace[ , Age := as.numeric(rn)]
HIV_peds_trace[ , HIV_SC := HIV_SC_1+HIV_SC_2+HIV_SC_3]
HIV_peds_trace[ , ART := ART_1+ART_2+ART_3+ART_4+AIDS_ART_1+AIDS_ART_2+AIDS_RD]
HIV_peds_trace <- HIV_peds_trace[ , c("Age", "HIV_SC", "AIDS_SC", "ART", "HIV_death", "death")]


HIV_peds_trace_melt <- melt(HIV_peds_trace, 
                            id.vars = c("Age"), variable.name = "Disease state",
                            value.name = "Percent of cohort")

HIV_adult_trace <- data.table(HIV_adult$markov_trace, keep.rownames = TRUE)
HIV_adult_trace[ , Age := as.numeric(rn)]
HIV_adult_trace[ , HIV_SC := HIV_SC_1+HIV_SC_2+HIV_SC_3]
HIV_adult_trace[ , ART := ART_1+ART_2+ART_3+ART_4+AIDS_ART_1+AIDS_ART_2+AIDS_RD]
HIV_adult_trace <- HIV_adult_trace[ , c("Age", "HIV_SC", "AIDS_SC", "ART", "HIV_death", "death")]


HIV_adult_trace_melt <- melt(HIV_adult_trace, 
                             id.vars = c("Age"), variable.name = "Disease state",
                             value.name = "Percent of cohort")

HIV_trace_both <- rbind(
  cbind(cohort = "Pediatric", HIV_peds_trace_melt),
  cbind(cohort = "Adult", HIV_adult_trace_melt)
)


dt_chronic_costs_by_state <- calc_costs_by_state(HIV_peds$cost_by_state_trace,
                                                 cohort="Pediatric",
                                                 disease_state_costs = costs_states[["HIV"]], "HIV")


dt_chronic_costs_by_state <- rbind(
  dt_chronic_costs_by_state,
  calc_costs_by_state(HIV_adult$cost_by_state_trace,
                      cohort="Adult",
                      disease_state_costs = costs_states[["HIV"]], "HIV")
)


ggplot(HIV_trace_both)+
  geom_line(aes(x = Age, y = `Percent of cohort`, color = `Disease state`))+
  geom_point(aes(x = Age, y = `Percent of cohort`, color = `Disease state`), size = 1.2)+
  facet_wrap(vars(cohort), ncol=1)+
  scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "bottom")


ggsave("manuscript/figs/hiv_trace.png",
       width=6.5,
       height = 6,
       units="in")

fwrite(HIV_peds_trace, "./results/markov_trace_HIV_peds.csv")
fwrite(HIV_adult_trace, "./results/markov_trace_HIV_peds.csv")





#  #  #  #
# HCV ###########
#  #  #  #

states_hcv <- c('acute_SC', 'acute_T', 'no_infection', 'chronic_SC', 'chronic_T', 
                'chronic_TF', 'CC_SC', 'CC_T', 'CC_TF', 'DCC_SC', 'DCC_T', 'DCC_TF', 
                'HCC', 'HCV_death')

#Create init vector with initial distribution
p_acute_T <- as.numeric(dt_transitions[Disease=="HCV" & From == "init" & To == "acute_T", value])
init_hcv <- c(1-p_acute_T, p_acute_T, rep(0, length(states_hcv)-1))
names(init_hcv) <- c(states_hcv,"death")

#BASE CASE
# transition matrics
t_hcv <- gen_t_matrix(dt_transitions, "hcv")

HCV_peds_bc <- sim_markov(init_hcv, t_hcv, 
                          age_start=5, 
                          all_cause_mort = all_cause_mort,
                          costs = costs_states$HCV,
                          DALY_weights = bc_daly_params$HCV
)


HCV_adult_bc <- sim_markov(init_hcv, t_hcv, 
                           age_start=40, 
                           all_cause_mort = all_cause_mort,
                           costs =  costs_states$HCV,
                           DALY_weights = bc_daly_params$HCV)


HCV_peds_trace <- data.table(HCV_peds_bc$markov_trace, keep.rownames = TRUE)
HCV_peds_trace[ , Age := as.numeric(rn)]
HCV_peds_trace[ , rn := NULL]

HCV_peds_trace_melt <- melt(HCV_peds_trace, 
                            id.vars = c("Age"), variable.name = "Disease state",
                            value.name = "Percent of cohort")

HCV_adult_trace <- data.table(HCV_adult_bc$markov_trace, keep.rownames = TRUE)
HCV_adult_trace[ , Age := as.numeric(rn)]
HCV_adult_trace[ , rn := NULL]


HCV_adult_trace_melt <- melt(HCV_adult_trace, 
                             id.vars = c("Age"), variable.name = "Disease state",
                             value.name = "Percent of cohort")


HCV_trace_both <- rbind(
  cbind(cohort = "Pediatric", HCV_peds_trace_melt),
  cbind(cohort = "Adult", HCV_adult_trace_melt)
)


dt_chronic_costs_by_state <- rbind(
  dt_chronic_costs_by_state,
  calc_costs_by_state(HCV_peds_bc$cost_by_state_trace,
                      cohort="Pediatric",
                      disease_state_costs = costs_states[["HCV"]], "HCV")
)


dt_chronic_costs_by_state <- rbind(
  dt_chronic_costs_by_state,
  calc_costs_by_state(HCV_adult_bc$cost_by_state_trace,
                      cohort="Adult",
                      disease_state_costs = costs_states[["HCV"]], "HCV")
)


ggplot(HCV_trace_both)+
  geom_line(aes(x = Age, y = `Percent of cohort`, color = `Disease state`))+
  geom_point(aes(x = Age, y = `Percent of cohort`, color = `Disease state`), size = 1.2)+
  facet_wrap(vars(cohort), ncol=1)+
  scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "bottom")

ggsave("./manuscript/figs/hcv_trace.png",
       width=6.5,
       height = 6,
       units="in")

fwrite(HCV_peds_trace, "./results/markov_trace_HCV_peds.csv")
fwrite(HCV_adult_trace, "./results/markov_trace_HCV_peds.csv")
#  #  #  #
# HBV ###########
#  #  #  #

states_hbv <- c('acute_SC', 'acute_CM', 'imm_tol_SC', 'imm_tol_CM', 'carrier_SC', 
                'carrier_CM', 'imm_react_SC', 'imm_react_T', 'chronic_SC', 
                'chronic_CM', 'chronic_T', 'CC_SC', 'CC_T', 'DCC_SC', 'DCC_T', 'HCC', 
                'no_infection', 'HBV_death')

#Create init vector with initial distribution
p_acute_CM <- as.numeric(dt_transitions[Disease=="HBV" & From == "init" & To == "acute_CM", value])
init_hbv <- c(1-p_acute_CM, p_acute_CM, rep(0, length(states_hbv)-1))





# transition matrics
t_hbv <- gen_t_matrix(dt_transitions, "hbv")



HBV_peds_bc <- sim_markov(init_hbv, t_hbv, 
                          age_start=5, 
                          all_cause_mort = all_cause_mort,
                          costs = costs_states$HBV,
                          DALY_weights = bc_daly_params$HBV
)


HBV_adult_bc <- sim_markov(init_hbv, t_hbv, 
                           age_start=40, 
                           all_cause_mort = all_cause_mort,
                           costs =  costs_states$HBV,
                           DALY_weights = bc_daly_params$HBV)

HBV_peds_trace <- data.table(HBV_peds_bc$markov_trace, keep.rownames = TRUE)
HBV_peds_trace[ , Age := as.numeric(rn)]
HBV_peds_trace[ , rn := NULL]

HBV_peds_trace_melt <- melt(HBV_peds_trace,
                            id.vars = c("Age"), variable.name = "Disease state",
                            value.name = "Percent of cohort")

HBV_adult_trace <- data.table(HBV_adult_bc$markov_trace, keep.rownames = TRUE)
HBV_adult_trace[ , Age := as.numeric(rn)]
HBV_adult_trace[ , rn := NULL]


HBV_adult_trace_melt <- melt(HBV_adult_trace,
                             id.vars = c("Age"), variable.name = "Disease state",
                             value.name = "Percent of cohort")


HBV_trace_both <- rbind(
  cbind(cohort = "Pediatric", HBV_peds_trace_melt),
  cbind(cohort = "Adult", HBV_adult_trace_melt)
)


dt_chronic_costs_by_state <- rbind(
  dt_chronic_costs_by_state,
  calc_costs_by_state(HBV_peds_bc$cost_by_state_trace,
                      cohort="Pediatric",
                      disease_state_costs = costs_states[["HBV"]], "HBV")
)


dt_chronic_costs_by_state <- rbind(
  dt_chronic_costs_by_state,
  calc_costs_by_state(HBV_adult_bc$cost_by_state_trace,
                      cohort="Adult",
                      disease_state_costs = costs_states[["HBV"]], "HBV")
)



ggplot(HBV_trace_both)+
  geom_line(aes(x = Age, y = `Percent of cohort`, color = `Disease state`))+
  geom_point(aes(x = Age, y = `Percent of cohort`, color = `Disease state`), size = 1.2)+
  facet_wrap(vars(cohort), ncol=1)+
  scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "bottom")

ggsave("manuscript/figs/hbv_trace.png",
       width=6.5,
       height = 6,
       units="in")
fwrite(HBV_peds_trace, "./results/markov_trace_HBV_peds.csv")
fwrite(HBV_adult_trace, "./results/markov_trace_HBV_peds.csv")

fwrite(dt_chronic_costs_by_state, "./results/chronic_costs_by_disease_state.csv")

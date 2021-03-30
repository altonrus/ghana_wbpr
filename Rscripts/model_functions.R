library(data.table)
library(mc2d)

sim_markov <- function(init, t_matrix, age_start, all_cause_mort,
                       costs, DALY_weights){
  # init <- init_hiv
  # age_start <- 5
  # t_matrix <- t_hiv_peds
  # costs<-state_costs$HIV
  #t matrix is missing death, we add in at each step
  
  #Simulate by year from age_start to age 95
  iters <- 95 - age_start + 1
  nstates <- nrow(t_matrix)+1
  state_names <- c(colnames(t_matrix), "death")
  
  markov_trace <- matrix( nrow = iters, ncol = length(state_names),
                          dimnames = list("Age" = age_start:95, 
                                          "State" = state_names))
  
  cost_by_state_trace <-  matrix( nrow = iters-1, ncol = length(state_names),
                                  dimnames = list("Age" = age_start:94, 
                                                  "State" = state_names))
  
  YLD_by_state_trace <-  matrix( nrow = iters-1, ncol = length(state_names),
                                  dimnames = list("Age" = age_start:94, 
                                                  "State" = state_names))
  # colnames(markov_trace) <- c(colnames(t_matrix), "death")
  # rownames(markov_trace) <- age_start:95
  
  # First row has initial condition
  markov_trace[as.character(age_start), ] <- init
  
  
  cost_mults <- structure(numeric(length(state_names)), names = state_names)
  YLD_mults <- structure(numeric(length(state_names)), names = state_names)
  
  for(state in state_names){
    cost_mults[state] <- unlist(costs[state])
    if(is.null(DALY_weights[[state]])){
      YLD_mults[state] <-0
    } else {
      YLD_mults[state] <- unlist(DALY_weights[state])
    }
    
  }
  

  
  
  for (i in 1:(iters-1)){
    curr_age = age_start+i
    #get probability of all-cause death
    p_death <- all_cause_mort[findInterval(curr_age, all_cause_mort$age)]$p_death
    #Create temporary transition matrix that encorporates probability of all cause death 
    #  and down-weights all other transitions so it still sums to 1
    t_matrix_temp <- rbind(cbind(t_matrix*(1-p_death), death = p_death), death = c(rep(0, ncol(t_matrix)), 1))
    #Correct the death for cause column (HIV_death, HCV_death, HBV death).
    # Must always be last row/col before death added
    t_matrix_temp[nrow(t_matrix_temp) - 1, ncol(t_matrix_temp) - 1 ] <- 1
    t_matrix_temp[nrow(t_matrix_temp) - 1, ncol(t_matrix_temp)] <- 0
    
    #Apply markovian transition
    current_dist <-  markov_trace[i , ]
    markov_trace[i+1, ] <- current_dist %*% t_matrix_temp
    
    prop_from_to <- (matrix(rep(current_dist, each = nstates), nrow = nstates, byrow = TRUE) *
                       t_matrix_temp)
    cost_by_state <- (0.5*rowSums(prop_from_to) + 0.5*colSums(prop_from_to))*cost_mults
    YLD_by_state <- (0.5*rowSums(prop_from_to) + 0.5*colSums(prop_from_to))*YLD_mults
    #Calculate costs at timestep using cycle tree correction
    cost_by_state_trace[i, ] <- cost_by_state
    YLD_by_state_trace[i, ] <- YLD_by_state
    
  }
  
  
  tot_cost_by_year <- rowSums(cost_by_state_trace)
  tot_YLD_by_year <- rowSums(YLD_by_state_trace)
  np_cost <- sum(tot_cost_by_year * 1.03^(0:(-1*(iters-2))))
  np_YLD <- sum(tot_YLD_by_year * 1.03^(0:(-1*(iters-2))))
  
  #YLL calculations
  normal_trace <- sim_normal_life(age_start)
  YLL_by_year <- (rowSums(markov_trace[, (ncol(markov_trace)-1):ncol(markov_trace)])-
                     normal_trace[,"dead"])
  #tot_YLL <- sum(YLL_by_year)
  np_YLL <- sum(YLL_by_year[-1] * 1.03^(0:(-1*(iters-2))))
  return(list("markov_trace" = markov_trace, 
              "tot_cost_by_year" = tot_cost_by_year,
              "np_cost" = np_cost, 
              "cost_by_state_trace" = cost_by_state_trace,
              "tot_YLD_by_year" = tot_YLD_by_year,
              "np_YLD" = np_YLD,
              "YLL_by_year"=YLL_by_year,
              "np_YLL"=np_YLL))
}

sim_normal_life <- function(age_start){
  #Simulate by year from age_start to age 95
  iters <- 95 - age_start + 1
  markov_trace <- matrix(nrow = iters, ncol=2,
                    dimnames = list("Age" = age_start:95, 
                                    "State" = c("alive","dead")))
  markov_trace[1, ] <- c(1,0)
  
  for (i in 1:(iters-1)){
    curr_age = age_start+i
    #get probability of all-cause death
    p_death <- all_cause_mort[findInterval(curr_age, all_cause_mort$age)]$p_death
    t_matrix_temp <- matrix(data=c(1-p_death, p_death, 0, 1),
                            nrow=2, byrow=TRUE)
    #Apply markovian transition
    markov_trace[i+1, ] <- markov_trace[i , ] %*% t_matrix_temp
  }
  return(markov_trace)
}





gen_t_matrix <- function(transitions, disease, cohort = "Adult"){
  transitions <- transitions[From != "init"]
  
  #State names (excludes all cause death; that's added in Sim_Markov)
  states <- list(hcv = c('acute_SC', 'acute_T', 'no_infection', 'chronic_SC', 'chronic_T', 
                         'chronic_TF', 'CC_SC', 'CC_T', 'CC_TF', 'DCC_SC', 'DCC_T', 'DCC_TF', 
                         'HCC', 'HCV_death'),
                 hbv = c('acute_SC', 'acute_CM', 'imm_tol_SC', 'imm_tol_CM', 'carrier_SC', 
                         'carrier_CM', 'imm_react_SC', 'imm_react_T', 'chronic_SC', 
                         'chronic_CM', 'chronic_T', 'CC_SC', 'CC_T', 'DCC_SC', 'DCC_T', 'HCC', 
                         'no_infection', 'HBV_death'),
                 hiv = c('HIV_SC_1', 'HIV_SC_2', 'HIV_SC_3', 'ART_1', 'ART_2', 'ART_3', 'ART_4', 
                         'AIDS_SC', 'AIDS_ART_1', 'AIDS_ART_2', 'AIDS_RD', 'HIV_death'))
  
  #Create empty transition vector
  t_matrix <- matrix(0L, nrow = length(states[[disease]]), ncol = length(states[[disease]]),
                  dimnames = list("From" = states[[disease]],
                              "To" = states[[disease]]))
  
  
  
  #Load in natural history parameters and treatment effectiveness
  t_nat_hist_params <- transitions[Disease == toupper(disease) & Category != "Treatment uptake" &
                                     `Base case` != "#" & Cohort %in% c("Both", cohort)]
  t_nat_hist_params[ , value := as.numeric(value)]
  for (row in 1:nrow(t_nat_hist_params)){
    t_matrix[t_nat_hist_params[row, From],
             t_nat_hist_params[row, To]] <- t_nat_hist_params[row, value]
  }
  
  #If HCV, handle special case of acute AVT treatment failure
  if(disease == "hcv"){
    t_matrix["acute_T", "no_infection"] <-  1 - (1-t_matrix["acute_SC", "no_infection"])*(1-t_matrix["acute_T", "no_infection"])
  }
  
  #Handle remainders for natural history model
  t_nat_hist_remainders <- transitions[Disease == toupper(disease) & Category != "Treatment uptake" & 
                                         `Base case` == "#" & Cohort %in% c("Both", cohort)]
  for (row in 1:nrow(t_nat_hist_remainders)){
    t_matrix[t_nat_hist_remainders[row, From],
             t_nat_hist_remainders[row, To]] <- 1 - sum(t_matrix[t_nat_hist_remainders[row, From], ])
  }
  
  #Factor in treatment uptake
  t_treatment <- transitions[Disease == toupper(disease) & Category == "Treatment uptake" &
                               Cohort %in% c("Both", cohort)]
  t_treatment[ , value := as.numeric(value)]
  for (row in 1:nrow(t_treatment)){
    t_matrix[t_treatment[row, From], ] <- t_matrix[t_treatment[row, From], ]*(1-t_treatment[row, value])
    t_matrix[t_treatment[row, From],
             t_treatment[row, To]] <- t_treatment[row, value]
  }
  

  

  
  t_matrix
  rowSums(t_matrix)
  
  
  return(t_matrix)
}
calc_disease_state_costs <- function(prm){
  cost_mult_acute_from_IPmort <- (1 - prm$p_IPmort + prm$p_IPmort*prm$p_acute_c_mort)
  
  
  c_disease_states <- list()
  c_disease_states$sep <- (prm$c_IP_day*prm$sep_d_IP_days + 
                               prm$sep_p_meds*prm$sep_c_meds)*cost_mult_acute_from_IPmort
  
  c_disease_states$mal <- (prm$c_IP_day*prm$mal_p_IP_days*prm$mal_d_IP_days +
                             prm$mal_p_OP*prm$mal_n_OP*prm$mal_c_OP+
                             prm$mal_c_diag+prm$mal_c_meds)*cost_mult_acute_from_IPmort
  c_disease_states$ftr <- (prm$ftr_c_meds+
                             prm$ftr_n_IP_days*prm$ftr_p_IP_days*prm$c_IP_day)*cost_mult_acute_from_IPmort
  c_disease_states$syp <- (prm$syp_p_diag*(prm$syp_c_diag+prm$syp_c_meds))*cost_mult_acute_from_IPmort
  
  #HIV
  c_HIV_states <- list()
  c_HIV_states$HIV_SC_1 <- 0
  c_HIV_states$HIV_SC_2 <- 0
  c_HIV_states$HIV_SC_3 <- 0
  c_HIV_states$ART_1 <- prm$hiv_cost_ART1
  c_HIV_states$ART_2 <- prm$hiv_cost_ART2
  c_HIV_states$ART_3 <- prm$hiv_cost_ART3
  c_HIV_states$ART_4 <- prm$hiv_cost_ART4
  c_HIV_states$AIDS_SC <- prm$hiv_cost_AIDS_SC
  c_HIV_states$AIDS_ART_1 <- prm$hiv_cost_AIDS_ART_1+prm$hiv_cost_ART1
  c_HIV_states$AIDS_ART_2 <- prm$hiv_cost_AIDS_ART_2+prm$hiv_cost_ART2
  c_HIV_states$AIDS_RD <- prm$hiv_cost_RD+prm$hiv_cost_ART4
  c_HIV_states$HIV_death <- 0
  c_HIV_states$death <- 0
  
  #HBV
  c_HBV_states <- list()
  c_HBV_states$acute_SC <- 0
  c_HBV_states$acute_CM <- (2*prm$hbv_c_HBsAg+prm$hbv_c_profile+
                           prm$hbv_c_DNAtest+prm$hep_c_OP_extensive+
                           prm$hep_c_OP_brief)
  c_HBV_states$no_infection <- 0
  c_HBV_states$imm_tol_CM <- (prm$hep_c_OP_extensive+prm$hbv_c_HBsAg+
                             prm$hbv_c_profile+prm$hbv_c_DNAtest+prm$hep_c_lft)
  c_HBV_states$imm_tol_SC <- 0
  c_HBV_states$carrier_CM <- 0
  c_HBV_states$carrier_SC <- 0
  c_HBV_states$imm_react_T <- (prm$hep_c_OP_extensive+prm$hbv_c_HBsAg+
                                 prm$hbv_c_profile+prm$hbv_c_antivirals+
                                 prm$hbv_c_DNAtest+prm$hep_c_lft)
  c_HBV_states$imm_react_SC <- 0
  c_HBV_states$chronic_SC <- 0
  c_HBV_states$chronic_CM <- (prm$hep_c_lft+prm$hep_c_BUNCE+
                             prm$hep_c_FBC+prm$hep_c_alphafeto+
                             prm$hep_c_ab_ultrasono+
                             prm$hep_c_OP_brief*prm$hbv_n_OP_chronic+
                             prm$hbv_c_profile + prm$hbv_c_DNAtest)
  c_HBV_states$chronic_T <- (prm$hep_c_lft+prm$hep_c_BUNCE+
                               prm$hep_c_FBC+prm$hep_c_alphafeto+
                               prm$hep_c_ab_ultrasono+
                               prm$hep_c_OP_brief*prm$hbv_n_OP_chronicT+
                               prm$hbv_c_antivirals+prm$hbv_c_profile + prm$hbv_c_DNAtest)
  c_HBV_states$CC_SC <- 0
  c_HBV_states$CC_T <- (prm$hep_c_lft + prm$hep_c_INR + prm$hep_c_BUNCE + 
                          prm$hep_c_FBC + 2*prm$hep_c_alphafeto + 
                          2*prm$hep_c_ab_ultrasono + 2*prm$hbv_c_DNAtest + 
                          prm$hep_c_OP_brief*prm$hbv_n_OP_CCT+
                          prm$hbv_c_antiviralsCirrhosis)
  c_HBV_states$DCC_SC <- 0
  c_HBV_states$DCC_T <- (prm$hep_c_lft + prm$hep_c_INR + prm$hep_c_BUNCE + 
                           prm$hep_c_FBC + 2*prm$hep_c_alphafeto + 
                           2*prm$hep_c_ab_ultrasono + 2*prm$hbv_c_DNAtest + 
                           2*prm$hep_c_band_lig + prm$hep_c_spiro + 
                           prm$hep_c_fluro + 
                           prm$hep_c_OP_brief*prm$hbv_n_OP_DCCT + 
                           prm$hbv_c_antiviralsCirrhosis)
  
  c_HBV_states$HCC <- (prm$hep_c_lft + prm$hep_c_INR + prm$hep_c_FBC + 
                         prm$hep_c_alphafeto + prm$hep_c_triphas_CT + 
                         prm$hep_c_band_lig + prm$hep_c_sorafenib +
                         prm$hep_c_chemoembo)
  c_HBV_states$HBV_death <- 0
  c_HBV_states$death <- 0
  
  c_HCV_states <- list()
  c_HCV_states$acute_SC <- 0
  c_HCV_states$acute_T <- (prm$hcv_c_screen + 2*prm$hcv_c_rna_test + 
                           prm$hcv_c_genotype + 
                             prm$hep_c_lft + prm$hep_c_BUNCE +
                             prm$hcv_c_antivirals+
                             prm$hep_c_OP_brief*prm$hcv_n_OP_acuteT)
                             
                            
  c_HCV_states$no_infection <- 0
  c_HCV_states$chronic_SC <- 0
  c_HCV_states$chronic_T <- (prm$hep_c_lft + prm$hep_c_BUNCE + 
                               prm$hep_c_FBC + prm$hep_c_alphafeto +
                               prm$hep_c_ab_ultrasono + 
                               prm$hcv_c_genotype + 2*prm$hcv_c_rna_test+
                               prm$hcv_n_OP_chronic*prm$hep_c_OP_brief+
                               prm$hcv_c_antivirals)
  c_HCV_states$chronic_TF <- (prm$hep_c_lft + prm$hep_c_BUNCE + 
                               prm$hep_c_FBC + prm$hep_c_alphafeto +
                               prm$hep_c_ab_ultrasono + 
                               prm$hcv_c_genotype + 2*prm$hcv_c_rna_test+
                               prm$hcv_n_OP_chronic*prm$hep_c_OP_brief)
  
  c_HCV_states$CC_SC <- 0
  c_HCV_states$CC_T <- (prm$hep_c_lft + prm$hep_c_INR + prm$hep_c_BUNCE + 
                          prm$hep_c_FBC + 2*prm$hep_c_alphafeto + 
                          2*prm$hep_c_ab_ultrasono + 2*prm$hcv_c_rna_test + 
                          prm$hcv_c_genotype+
                          prm$hcv_n_OP_cc*prm$hep_c_OP_brief+
                          prm$hcv_c_antivirals)
  c_HCV_states$CC_TF <- (prm$hep_c_lft + prm$hep_c_INR + prm$hep_c_BUNCE + 
                          prm$hep_c_FBC + 2*prm$hep_c_alphafeto + 
                          2*prm$hep_c_ab_ultrasono + 2*prm$hcv_c_rna_test + 
                          prm$hcv_c_genotype+
                          prm$hcv_n_OP_cc*prm$hep_c_OP_brief)
  c_HCV_states$DCC_SC <- 0
  
  c_HCV_states$DCC_T <- (prm$hep_c_lft + prm$hep_c_INR + prm$hep_c_BUNCE + 
                           prm$hep_c_FBC + 2*prm$hep_c_alphafeto + 
                           2*prm$hep_c_ab_ultrasono + 3*prm$hcv_c_rna_test + 
                           2*prm$hep_c_band_lig + prm$hep_c_spiro + 
                           prm$hcv_c_genotype+
                           prm$hep_c_fluro  + prm$hcv_n_OP_dcc*prm$hep_c_OP_brief+
                           prm$hcv_c_antiviral_dcc)
  c_HCV_states$DCC_TF <- (prm$hep_c_lft + prm$hep_c_INR + prm$hep_c_BUNCE + 
                           prm$hep_c_FBC + 2*prm$hep_c_alphafeto + 
                           2*prm$hep_c_ab_ultrasono + 3*prm$hcv_c_rna_test + 
                           2*prm$hep_c_band_lig + prm$hep_c_spiro + 
                           prm$hcv_c_genotype+
                           prm$hep_c_fluro  + prm$hcv_n_OP_dcc*prm$hep_c_OP_brief)
  
  c_HCV_states$HCC <- (prm$hep_c_lft + prm$hep_c_INR + prm$hep_c_FBC + 
                         prm$hep_c_alphafeto + prm$hep_c_triphas_CT + 
                         prm$hep_c_band_lig + prm$hep_c_sorafenib +
                         prm$hep_c_chemoembo)
  c_HCV_states$HCV_death <- 0
  c_HCV_states$death <- 0
  
  c_disease_states$HIV <- c_HIV_states
  c_disease_states$HCV <- c_HCV_states
  c_disease_states$HBV <- c_HBV_states
  return(c_disease_states)
}

calc_acute_dalys <- function(prm_c, prm_d, params){
  #Do not count DALYs from ip mortality unrelated to disease
  # prm_c <- bc_microcost_params
  # prm_d <- bc_daly_params
  
  
  
  YLL_np_peds <- sum(sim_normal_life(5)[-1,1] * 1.03^(0:(-1*(89))))
  YLL_np_adult <-sum(sim_normal_life(40)[-1,1] * 1.03^(0:(-1*(54))))
  YLL_per_death <- (params$p_peds*YLL_np_peds+
                      (1-params$p_peds)*YLL_np_adult)
  

  YLD <- list()
  YLD$sep <- (prm_d$sep_dw_ip*prm_c$sep_d_IP_days/365+#IP DALYs
                prm_d$sep_p_sequelae*prm_d$sep_dur_sequelae*prm_d$sep_dw_sequelae/365 #sequelae
              )*(1-prm_c$p_IPmort)*(1-prm_d$sep_p_death)
  YLD$mal <- (prm_d$mal_dur_dw*prm_d$mal_dw/365)*(1-prm_c$p_IPmort)*
    (1 - params$p_peds*prm_d$mal_peds_p_death)
  YLD$ftr <- (prm_d$ftr_dw*prm_d$ftr_dur_dw/365)*(1-prm_c$p_IPmort)
  YLD$syp <- (prm_d$syp_dw*prm_d$syp_dur_dw/365)*(1-prm_c$p_IPmort)
  
  YLL <- list()
  YLL$sep <- (prm_d$sep_p_death*YLL_per_death #mortality DALYs
              )*(1-prm_c$p_IPmort)
  YLL$mal <- params$p_peds*prm_d$mal_peds_p_death*YLL_np_peds*(1-prm_c$p_IPmort)
  YLL$ftr <- 0
  YLL$syp <- 0
  
  return(list("YLD" = YLD, "YLL"=YLL))
}

rtriang<- function(mean, lb, ub, n=1){
  mode = 3*mean - lb - ub
  if (mode < lb){
    print(paste0(
      "Cannot achieve desired expected value with lower bound of ", lb, "."
    ))
    lb = (3*mean - ub)/2
    mode=lb
    print(paste0(
      "Reset lb to ", lb, "."
    ))
  } else if (mode > ub) {
    print(paste0(
      "Cannot achieve desired expected value with upper bound of ", ub, "."
    ))
    ub = (3*mean - lb)/2
    mode=ub
    print(paste0(
      "Reset upper bound to ", ub, "."
    ))
  } 
  
  Fmode <- (mode - lb)/(ub - lb)
  U <- runif(n)
  rtris <- ifelse(U < Fmode, 
                  lb + (U*(ub - lb)*(mode - lb))^.5,
                  ub - ((1-U)*(ub - lb)*(ub - mode))^.5
                  )
  
  return(rtris)
}
rpert_mean<- function(mean, lb, ub, n=1, adjust_bounds = FALSE){
  mode = (6*mean - lb - ub)/4
  if (adjust_bounds == TRUE){
    if (mode < lb){
      print(paste0(
        #"Cannot achieve desired expected value with upper bound of ", ub, "."
      ))
      ub = 6*mean - 5*lb
      mode=lb
      print(paste0(
        #"Reset upper bound to ", ub, "."
      ))
    } else if (mode > ub) {
      print(paste0(
        #"Cannot achieve desired expected value with lower bound of ", lb, "."
      ))
      lb = 6*mean - 5*ub
      mode=ub
      print(paste0(
        #"Reset lower bound to ", lb, "."
      ))
    } 
  } else {
    if (mode < lb){
      mode = lb
      #print(paste0("Set mode to lower bound. Mean will be ", (5*lb+ub)/6,
      #             "instead of ", mean))
    } else if (mode > ub) {
      mode = ub
      #print(paste0("Set mode to lower bound. Mean will be ", (5*ub+lb)/6,
      #             "instead of ", mean))
    } 
  }
  return(rpert(n, lb, mode, ub))
}
sample_params<- function(mean, lb,ub, dist,Param1, Param2=Param2){
  if(dist == "PERT" & !is.na(dist)){
    return(rpert_mean(mean, lb, ub))
  } else if (dist == "Beta" & !is.na(dist) ) {
    return(rbeta(1, Param1, Param2))
  } else{
    return(mean)
  }
}
calc_outcomes <- function(costs_np, dalys_np, params,
                          WTP=2202){
  aes <- c("hiv", "sep", "hcv", "hbv",
          "syp", "mal", "ftr")
  outcomes <- list()
  
  outcomes["prt_cost"] <- params$c_PI*(1+params$p_waste)*params$n_comp
  
  for (ae in aes){
    cases_no_prt <-  (params$n_comp*
                       params[[paste0("bl_",ae)]]*
                       params[[paste0("p_clin_",ae)]])
    
    cases_prt <- cases_no_prt/params[[paste0("rred_",ae)]]
    
    cases_reduced <- cases_no_prt - cases_prt 
    
    burden_no_prt <- cases_no_prt*costs_np[[ae]]
    burden_prt <- cases_prt*costs_np[[ae]]
    burden_reduced <- cases_reduced*costs_np[[ae]]
    
    YLL_no_prt <- cases_no_prt*dalys_np$YLL[[ae]]
    YLL_prt <- cases_prt*dalys_np$YLL[[ae]]
    YLL_reduced <- cases_reduced*dalys_np$YLL[[ae]]
    
    YLD_no_prt <- cases_no_prt*dalys_np$YLD[[ae]]
    YLD_prt <- cases_prt*dalys_np$YLD[[ae]]
    YLD_reduced <- cases_reduced*dalys_np$YLD[[ae]]
    
    
    
    outcomes[[ae]]<- list(
      "cases_no_prt" = cases_no_prt,
      "cases_prt" = cases_prt,
      "cases_reduced" = cases_reduced,
      "burden_no_prt" = burden_no_prt,
      "burden_prt" = burden_prt,
      "burden_reduced" = burden_reduced,
      "YLL_no_prt" = YLL_no_prt,
      "YLL_prt" = YLL_prt,
      "YLL_reduced" = YLL_reduced,
      "YLD_no_prt" = YLD_no_prt,
      "YLD_prt" = YLD_prt,
      "YLD_reduced" = YLD_reduced,
      "DALY_no_prt" = YLD_no_prt+YLL_no_prt,
      "DALY_prt" = YLD_prt+YLL_prt,
      "DALY_reduced" = YLD_reduced+YLL_reduced      
    )
  }
  
  dt_outcomes <- data.table(rbind(unlist(outcomes)))
  #Cases
  dt_outcomes[ , all_ae.cases_no_prt := sum(.SD), .SDcols = paste0(aes,".cases_no_prt")]
  dt_outcomes[ , all_ae.cases_prt := sum(.SD), .SDcols = paste0(aes,".cases_prt")]
  dt_outcomes[ , all_ae.cases_reduced := sum(.SD), .SDcols = paste0(aes,".cases_reduced")]
  #Costs
  dt_outcomes[ , all_ae.burden_no_prt := sum(.SD), .SDcols = paste0(aes,".burden_no_prt")]
  dt_outcomes[ , all_ae.burden_prt := sum(.SD), .SDcols = paste0(aes,".burden_prt")]
  dt_outcomes[ , all_ae.burden_reduced := sum(.SD), .SDcols = paste0(aes,".burden_reduced")]
  dt_outcomes[ , net_savings := all_ae.burden_reduced - prt_cost]
  #YLL
  dt_outcomes[ , all_ae.YLL_no_prt := sum(.SD), .SDcols = paste0(aes,".YLL_no_prt")]
  dt_outcomes[ , all_ae.YLL_prt := sum(.SD), .SDcols = paste0(aes,".YLL_prt")]
  dt_outcomes[ , all_ae.YLL_reduced := sum(.SD), .SDcols = paste0(aes,".YLL_reduced")]
  #YLD
  dt_outcomes[ , all_ae.YLD_no_prt := sum(.SD), .SDcols = paste0(aes,".YLD_no_prt")]
  dt_outcomes[ , all_ae.YLD_prt := sum(.SD), .SDcols = paste0(aes,".YLD_prt")]
  dt_outcomes[ , all_ae.YLD_reduced := sum(.SD), .SDcols = paste0(aes,".YLD_reduced")]
  #DALYs
  dt_outcomes[ , all_ae.DALY_no_prt := sum(.SD), .SDcols = paste0(aes,".DALY_no_prt")]
  dt_outcomes[ , all_ae.DALY_prt := sum(.SD), .SDcols = paste0(aes,".DALY_prt")]
  dt_outcomes[ , all_ae.DALY_reduced := sum(.SD), .SDcols = paste0(aes,".DALY_reduced")]
  dt_outcomes[ , nmb := all_ae.DALY_reduced*WTP+net_savings]

  
  return(dt_outcomes)
}
calc_univariate_change <- function(param, value, dt_transitions, 
                                   dt_microcost_params, dt_riskmod_params,
                                   dt_daly_params){
  #FORMAT DATA FOR MODEL
  dt_transitions[ , value := `Base case`]
  bc_microcost_params <- as.list(dt_microcost_params$`Base case`)
  names(bc_microcost_params) <- dt_microcost_params$rname
  #Get risk model params
  params <- as.list(unlist(dt_riskmod_params[ , `Base case`]))
  names(params) <- dt_riskmod_params$rname 
  #get YLL and YLD for each 
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
  
  #Make parameter substitution
  if(param %in% dt_transitions$rname){
    dt_transitions[rname == param]$value <-  value
  } else if (param %in% names(bc_microcost_params)){
    bc_microcost_params[[param]] <- value
  } else if (param %in% names(params)){
    params[[param]] <- value
  } else if (param %in% dt_daly_params$rname_long){
    if(substr(param,1,3)=="hiv"){
      bc_daly_params[["hiv"]][[param]] <- value
    } else if(substr(param,1,3)=="hbv"){
      bc_daly_params[["hbv"]][[param]] <- value
    } else if(substr(param,1,3)=="hcv"){
      bc_daly_params[["hcv"]][[param]] <- value
    } else{
      bc_daly_params[[param]] <- value
    }
  } else {
    print(paste0("error: ", param, " not found."))
    return(NA)
  }
  
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
  outcomes <- calc_outcomes(costs_np, dalys_np, params)
  
  return(list("net_savings" = outcomes$net_savings,
         "DALY_averted" = outcomes$all_ae.DALY_reduced,
         "nmb" = outcomes$nmb))
}
hiv_calc_np_outcomes <- function(state_costs, markov_probs, params, daly_params, IPmort){
  # state_costs <- costs_states
  # markov_probs <- dt_transitions
  # daly_params <- bc_daly_params
  # IPmort <- bc_microcost_params$p_IPmort
  
  states_hiv <- c('HIV_SC_1', 'HIV_SC_2', 'HIV_SC_3', 'ART_1', 'ART_2', 'ART_3', 'ART_4', 
                  'AIDS_SC', 'AIDS_ART_1', 'AIDS_ART_2', 'AIDS_RD', 'HIV_death')
  
  #Create init vector with initial distribution
  #Create init vector with initial distribution
  p_ART_init <- as.numeric(markov_probs[Disease=="HIV" & From == "init" & To == "ART_1", value])
  init_hiv <- rep(0, length(states_hiv)+1)
  names(init_hiv) <- c(states_hiv,"death")
  init_hiv["ART_1"] <- p_ART_init
  init_hiv["HIV_SC_1"] <- 1-p_ART_init
  
  
  #Transition matrices
  t_hiv_peds <- gen_t_matrix(markov_probs, "hiv", "Pediatric")
  t_hiv_adult <- gen_t_matrix(markov_probs, "hiv", "Adult")
  #Run markov model
  HIV_peds <- sim_markov(init_hiv, t_hiv_peds,
                         age_start=5, 
                         all_cause_mort = all_cause_mort,
                         costs = state_costs$HIV,
                         DALY_weights = daly_params$hiv
  )
  HIV_adult <- sim_markov(init_hiv, t_hiv_adult, 
                             age_start=40, 
                             all_cause_mort = all_cause_mort,
                             costs =  state_costs$HIV,
                          DALY_weights = daly_params$hiv)
  #Take weighted average of peds and adult
  return(list(np_cost = (params$p_peds*HIV_peds$np_cost + 
                           (1-params$p_peds)*HIV_adult$np_cost)*(1-IPmort),
              np_YLL = (params$p_peds*HIV_peds$np_YLL + 
                            (1-params$p_peds)*HIV_adult$np_YLL)*(1-IPmort),
              np_YLD = (params$p_peds*HIV_peds$np_YLD + 
                          (1-params$p_peds)*HIV_adult$np_YLD)*(1-IPmort)))
  
}
hcv_calc_np_outcomes<- function(state_costs, markov_probs, params, daly_params, IPmort){
  states_hcv <- c('acute_SC', 'acute_T', 'no_infection', 'chronic_SC', 'chronic_T', 
                  'chronic_TF', 'CC_SC', 'CC_T', 'CC_TF', 'DCC_SC', 'DCC_T', 'DCC_TF', 
                  'HCC', 'HCV_death')
  
  #Create init vector with initial distribution
  p_acute_T <- as.numeric(markov_probs[Disease=="HCV" & From == "init" & To == "acute_T", value])
  init_hcv <- c(1-p_acute_T, p_acute_T, rep(0, length(states_hcv)-1))
  names(init_hcv) <- c(states_hcv,"death")
  
  #Transition matrices
  t_hcv <- gen_t_matrix(markov_probs, "hcv")
  #Run markov model
  HCV_peds <- sim_markov(init_hcv, t_hcv, 
                            age_start=5, 
                            all_cause_mort = all_cause_mort,
                            costs = costs_states$HCV,
                         DALY_weights = daly_params$hcv
  )
  HCV_adult <- sim_markov(init_hcv, t_hcv, 
                             age_start=40, 
                             all_cause_mort = all_cause_mort,
                             costs =  costs_states$HCV,
                          DALY_weights = daly_params$hcv)

  #Take weighted average of peds and adult
  return(list(np_cost = (params$p_peds*HCV_peds$np_cost + 
                           (1-params$p_peds)*HCV_adult$np_cost)*(1-IPmort),
              np_YLL = (params$p_peds*HCV_peds$np_YLL + 
                          (1-params$p_peds)*HCV_adult$np_YLL)*(1-IPmort),
              np_YLD = (params$p_peds*HCV_peds$np_YLD + 
                          (1-params$p_peds)*HCV_adult$np_YLD)*(1-IPmort)))
  
}
hbv_calc_np_outcomes<- function(state_costs, markov_probs, params, daly_params, IPmort){
  states_hbv <- c('acute_SC', 'acute_CM', 'imm_tol_SC', 'imm_tol_CM', 'carrier_SC', 
                  'carrier_CM', 'imm_react_SC', 'imm_react_T', 'chronic_SC', 
                  'chronic_CM', 'chronic_T', 'CC_SC', 'CC_T', 'DCC_SC', 'DCC_T', 'HCC', 
                  'no_infection', 'HBV_death')
  
  #Create init vector with initial distribution
  p_acute_CM <- as.numeric(markov_probs[Disease=="HBV" & From == "init" & To == "acute_CM", value])
  init_hbv <- c(1-p_acute_CM, p_acute_CM, rep(0, length(states_hbv)-1))
  #Transition matrices
  t_hbv <- gen_t_matrix(markov_probs, "hbv")
  #Run markov model
  HBV_peds <- sim_markov(init_hbv, t_hbv, 
                            age_start=5, 
                            all_cause_mort = all_cause_mort,
                            costs = costs_states$HBV,
                         DALY_weights = daly_params$hbv
  )
  
  
  HBV_adult <- sim_markov(init_hbv, t_hbv, 
                             age_start=40, 
                             all_cause_mort = all_cause_mort,
                             costs =  costs_states$HBV,
                          DALY_weights = daly_params$hbv)
  #Take weighted average of peds and adult
  return(list(np_cost = (params$p_peds*HBV_peds$np_cost + 
                           (1-params$p_peds)*HBV_adult$np_cost)*(1-IPmort),
              np_YLL = (params$p_peds*HBV_peds$np_YLL + 
                          (1-params$p_peds)*HBV_adult$np_YLL)*(1-IPmort),
              np_YLD = (params$p_peds*HBV_peds$np_YLD + 
                          (1-params$p_peds)*HBV_adult$np_YLD)*(1-IPmort)))
}
calc_costs_by_state <- function(cost_by_state_trace, cohort, disease_state_costs, disease){
  return(
    dt_costs_by_state <- data.table(
      "Disease" = disease,
      "Disease state" = colnames(cost_by_state_trace),
      "Annual cost" = disease_state_costs[colnames(cost_by_state_trace)],
      "Cohort" = cohort,
      "Undiscounted Lifetime cost" = colSums(cost_by_state_trace),
      "Net present lifetime cost" = colSums(cost_by_state_trace *  1.03^(0:(-1*(nrow(cost_by_state_trace)-1))))
    )
  )
}

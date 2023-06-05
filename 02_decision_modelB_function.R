### DECISION MODEL ###

#' \code{run_HPMBC_MX_ODE} implements the decision model
#' @param l_params_all List with all parameters of decision model

run_HPMBC_MX_ODE <- function(t, v_x, l_params_all) { 
  with(
    as.list(l_params_all),
    { 
      n_time <- floor(t)
      m_states <- matrix(data     = v_x, 
                         nrow     = n_ages,
                         ncol     = n_states, 
                         dimnames = list(v_names_ages, v_names_states))
      
      v_H          <- m_states[,"H"]
      v_preS1      <- m_states[,"preS1"]
      v_S1         <- m_states[,"S1"]
      v_S1T        <- m_states[,"S1T"]
      v_preS2      <- m_states[,"preS2"]
      v_S2         <- m_states[,"S2"]
      v_S2T        <- m_states[,"S2T"]
      v_preS3      <- m_states[,"preS3"]
      v_S3         <- m_states[,"S3"]
      v_S3T        <- m_states[,"S3T"]
      v_preS4      <- m_states[,"preS4"]
      v_S4         <- m_states[,"S4"]
      v_S4T        <- m_states[,"S4T"]
      v_BCD        <- m_states[,"BCD"]
      v_BCD_detect <- m_states[,"BCD_detect"]
      v_BCD_trt    <- m_states[,"BCD_trt"]
      v_OCD        <- m_states[,"OCD"]
      v_AllD       <- m_states[,"AllD"]
       v_i         <- m_states[,"I"]
       v_i_scr     <- m_states[,"I_Scr"]
       v_i_sym     <- m_states[,"I_Sym"]
       v_i_S1      <- m_states[,"I_S1"]
       v_i_S1_scr  <- m_states[,"S1_I_scr"]
       v_i_S1_sym  <- m_states[,"S1_I_sym"]
       v_i_S2      <- m_states[,"I_S2"]
       v_i_S2_scr  <- m_states[,"S2_I_scr"]
       v_i_S2_sym  <- m_states[,"S2_I_sym"]
       v_i_S3      <- m_states[,"I_S3"]
       v_i_S3_scr  <- m_states[,"S3_I_scr"]
       v_i_S3_sym  <- m_states[,"S3_I_sym"]
       v_i_S4      <- m_states[,"I_S4"]
       v_i_S4_scr  <- m_states[,"S4_I_scr"]
       v_i_S4_sym  <- m_states[,"S4_I_sym"]
      v_dtc        <- m_states[,"Dtc"]
      v_trt        <- m_states[,"Trt"]
      v_N          <- m_states[,"N"]
      
      
      #Create internal parameters
      v_mu <- l_params_all$m_r_demogmx[, "1"] # this is 1990, mortality 
      v_N <- (v_H + v_preS1 + v_S1 + v_S1T + v_preS2 + v_S2 + v_S2T + v_preS3 + v_S3 + v_S3T + v_preS4 + v_S4 + v_S4T)
      N = sum(v_N)
      v_beta_1   <- v_beta_1  #Breast cancer onset (lambda)
      v_beta_2   <- beta_2    #Progression form pre-clinical stage 1 to pre-clinical stage 2
      v_beta_3   <- beta_3    #Progression form pre-clinical stage 2 to pre-clinical stage 3
      v_beta_4   <- beta_4    #Progression form pre-clinical stage 3 to pre-clinical stage 4
      v_gamma_1  <- gamma_1   #Mortality of being detected for breast cancer but not treated
      v_gamma_2  <- gamma_2   #Mortality of being detected for breast cancer but not treated
      v_gamma_3  <- gamma_3   #Mortality of being detected for breast cancer but not treated
      v_gamma_4  <- gamma_4   #Mortality of being detected for breast cancer but not treated
      v_gammaT_1 <- gammaT_1  #Mortality of being detected for breast cancer and treated
      v_gammaT_2 <- gammaT_2  #Mortality of being detected for breast cancer and treated
      v_gammaT_3 <- gammaT_3  #Mortality of being detected for breast cancer and treated
      v_gammaT_4 <- gammaT_4  #Mortality of being detected for breast cancer and treated
      v_delta_1  <- delta_1   #Detection rate for stage 1
      v_delta_2  <- delta_2   #Detected rate for stage 2
      v_delta_3  <- delta_3   #Detected rate for stage 3
      v_delta_4  <- delta_4   #Detected rate for stage 4 
      v_phi      <- phi       #Time to event - time rate interval for waiting for healthcare services
      v_eta      <- v_eta     #Screening rate 
      v_rho      <- v_rho     #Sensitivity
      
      #Aging parameters in the model 
      # Vector for H pop growth 
      v_H_g <- c(N*r_birth, (v_r_aging*v_H)[-n_ages])
      #Vector of aging for preS1
      v_preS1_g <- c(0, (v_r_aging*v_preS1)[-n_ages])
      #Vector of aging for S1
      v_S1_g <- c(0, (v_r_aging*v_S1)[-n_ages])
      #Vector of aging for preS2
      v_preS2_g <- c(0, (v_r_aging*v_preS2)[-n_ages])
      #Vector of aging for S2
      v_S2_g <- c(0, (v_r_aging*v_S2)[-n_ages])
      #Vector of aging for preS3
      v_preS3_g <- c(0, (v_r_aging*v_preS3)[-n_ages])
      #Vector of aging for S3
      v_S3_g <- c(0, (v_r_aging*v_S3)[-n_ages])
      #Vector of aging for preS4
      v_preS4_g <- c(0, (v_r_aging*v_preS4)[-n_ages])
      #Vector of aging for S4
      v_S4_g <- c(0, (v_r_aging*v_S4)[-n_ages])
      #Vector of aging for S1 - trt
      v_S1T_g <- c(0, (v_r_aging*v_S1T)[-n_ages])
      #Vector of aging for S2 - trt
      v_S2T_g <- c(0, (v_r_aging*v_S2T)[-n_ages])
      #Vector of aging for S3 - trt
      v_S3T_g <- c(0, (v_r_aging*v_S3T)[-n_ages])
      #Vector of aging for S4 - trt
      v_S4T_g <- c(0, (v_r_aging*v_S4T)[-n_ages])
      
      ## Population Model
      ### Define differential equations
      ## Healthy 
      v_dH <- v_H_g -                                              # Incoming growth
        (v_r_aging + v_beta_1 + v_mu)*v_H                          # Leaving H

      ## Pre-clinical S1 
      v_dpreS1 <- v_preS1_g -                                                  # Incoming growth
        (v_r_aging + v_beta_2 + v_mu + v_delta_1 + (v_eta*v_rho))*v_preS1 +    # Leaving PreS1
        v_beta_1*v_H                                                           # From H

      ## Clinical S1
      v_dS1 <- v_S1_g -                                              # Incoming growth
        (v_r_aging + v_mu + v_gamma_1 + v_beta_2 + v_phi)*v_S1 +     # Leaving S1
        (v_delta_1 + v_eta)*v_preS1                                  # From PreS1

      ## Clinical S1 with trt
      v_dS1T <- v_S1T_g -                                # Incoming growth 
      (v_r_aging + v_mu + v_gammaT_1)*v_S1T +            # Leaving S1T
      (v_phi)*v_S1                                       # Incoming from S1
     
       ## Pre-clinical S2
      v_dpreS2 <- v_preS2_g -                                                 # Incoming growth
        (v_r_aging + v_beta_3 + v_mu + v_delta_2 + (v_eta*v_rho))*v_preS2 +   # Leaving PreS2
        v_beta_2*v_preS1                                                      # From PreS1

      ## Clinical S2
      v_dS2 <- v_S2_g -                                                                  # Incoming growth
        (v_r_aging + v_mu + v_gamma_2 + v_beta_3 + v_phi)*v_S2 +                         # Leaving S2
        (v_delta_2 + (v_eta*v_rho))*v_preS2 +                                            # From PreS2
        v_beta_2*v_S1                                                                    # From S1

      ## Clinical S2 with trt
      v_dS2T <- v_S2T_g -                                   # Incoming growth 
        (v_r_aging + v_mu + v_gammaT_2)*v_S2T +             # Leaving S2T
        (v_phi)*v_S2                                        # Incoming from S2
      
      ## Pre-clinical S3 
      v_dpreS3 <- v_preS3_g -                                                 # Incoming
        (v_r_aging + v_beta_4 + v_mu + v_delta_3 + (v_eta*v_rho))*v_preS3 +   # Leaving PreS3
        v_beta_3*v_preS2                                                      # From PreS2

      ## Clinical S3
      v_dS3 <- v_S3_g -                                                   # Incoming
        (v_r_aging + v_mu + v_gamma_3 + v_beta_4 + v_phi)*v_S3 +          # Leaving S3
        (v_delta_3 + (v_eta*v_rho))*v_preS3 +                             # From PreS3
        beta_3*v_S2                                                       # From S2
        
      ## Clinical S3 with trt
      v_dS3T <- v_S3T_g -                                   # Incoming growth 
        (v_r_aging + v_mu + v_gammaT_3)*v_S3T +             # Leaving S3T
        (v_phi)*v_S3                                        # Incoming from S3
      
      ## Pre-clinical S4  
      v_dpreS4 <- v_preS4_g -                                       # Incoming
        (v_r_aging + v_mu + v_delta_4 + (v_eta*v_rho))*v_preS4 +    # Leaving PreS4
        v_beta_4*v_preS3                                            # From PreS3

      ## Clinical S4
      v_dS4 <- v_S4_g -                                      # Incoming
        (v_r_aging + v_mu + v_gamma_4 + v_phi)*v_S4 +        # Leaving S4
        (v_delta_4 + (v_eta*v_rho))*v_preS4 +                # From PreS4
        v_S3*beta_4                                          # From S3

      ## Clinical S4 with trt
      v_dS4T <- v_S4T_g -                                   # Incoming growth 
        (v_r_aging + v_mu + v_gammaT_4)*v_S4T +             # Leaving S4T
        (v_phi)*v_S4                                        # Incoming from S4
      
      
      ## Breast Cancer Death
      v_BCD <- v_gamma_1*v_S1 +                              # Death from S1
        v_gamma_2*v_S2 +                                     # Death from S2
        v_gamma_3*v_S3 +                                     # Death from S3
        v_gamma_4*v_S4 +                                     # Death from S4
        v_gammaT_1*v_S1T +                                   # Death from S1T
        v_gammaT_2*v_S2T +                                   # Death from S2T
        v_gammaT_3*v_S3T +                                   # Death from S3T
        v_gammaT_4*v_S4T                                     # Death from S4T
      
      v_BCD_detect <- 
        v_gamma_1*v_S1 +                                    # Death from S1
        v_gamma_2*v_S2 +                                    # Death from S2
        v_gamma_3*v_S3 +                                    # Death from S3
        v_gamma_4*v_S4                                      # Death from S4
        
        v_BCD_trt <- 
        v_gammaT_1*v_S1T +                                   # Death from S1T
        v_gammaT_2*v_S2T +                                   # Death from S2T
        v_gammaT_3*v_S3T +                                   # Death from S3T
        v_gammaT_4*v_S4T                                     # Death from S4T

      ## Death from All Causes  
      v_OCD <- v_mu*v_H +                                 # Death from Health
        v_mu*v_preS1 +                                    # Death from PreS1
        v_mu*v_preS2 +                                    # Death from PreS2
        v_mu*v_preS3 +                                    # Death from PreS3
        v_mu*v_preS4 +                                    # Death from PreS4
        v_mu*v_S1 +                                       # Death from S1
        v_mu*v_S2 +                                       # Death from S2
        v_mu*v_S3 +                                       # Death from S3
        v_mu*v_S4 +                                       # Death from S4
        v_mu*v_S1T +                                      # Death from S1T
        v_mu*v_S2T +                                      # Death from S2T
        v_mu*v_S3T +                                      # Death from S3T
        v_mu*v_S4T                                        # Death from S4T

      ## Total Death
      v_AllD <- v_BCD + v_OCD                             # All Death

      ## Total Incidence
      v_i <- v_preS1*((v_eta*v_rho) + v_delta_1) +
        v_preS2*((v_eta*v_rho) + v_delta_2) +
        v_preS3*((v_eta*v_rho) + v_delta_3) +
        v_preS4*((v_eta*v_rho) + v_delta_4)
    
      
      ## Screening Incidence - incidence detected by screening
      v_i_scr <- (v_preS1 + v_preS2 + v_preS3 + v_preS4)*(v_eta*v_rho)
      
      ##Symptom Incidence - incidence detected by symptoms
      v_i_sym <-v_preS1*(v_delta_1) + v_preS2*(v_delta_2) + v_preS3*(v_delta_3) + v_preS4*(v_delta_4)
      
      ## Incidence by stage  
      #Stage 1 
      v_i_S1     <- v_preS1*((v_eta*v_rho) + v_delta_1)
      v_i_S1_scr <- v_preS1*(v_eta*v_rho)
      v_i_S1_sym <- v_preS1*(v_delta_1)
      
      #Stage 2  
      v_i_S2        <- v_preS2*((v_eta*v_rho) + v_delta_2) 
      v_i_S2_scr    <- v_preS2*(v_eta*v_rho) 
      v_i_S2_sym    <- v_preS2*(v_delta_2) 
      
      #Stage 3
      v_i_S3     <- v_preS3*((v_eta*v_rho) + v_delta_3) 
      v_i_S3_scr <- v_preS3*(v_eta*v_rho)
      v_i_S3_sym <- v_preS3*(v_delta_3) 
      
      #Stage 4  
      v_i_S4 <- v_preS4*((v_eta*v_rho) + v_delta_4)
      v_i_S4_scr <- v_preS4*(v_eta*v_rho)
      v_i_S4_sym <- v_preS4*(v_delta_4)
      
      # The number of women treated 
      v_trt <-
        v_S1*v_phi +
        v_S2*v_phi +
        v_S3*v_phi +
        v_S4*v_phi
      
      #The number of women transitioning in detected and no treated 
      v_dtc <-    
        v_S1*beta_2 +
        v_S2*beta_3 +
        v_S3*beta_4
      
      # Total Population 
      v_N <- (v_H + v_preS1 + v_S1 + v_preS2 + v_S2 + v_preS3 + v_S3 + v_preS4 + v_S4 + v_S1T + v_S2T + v_S3T + v_S4T)
      
      v_dx <- c(v_dH, v_dpreS1, v_dS1, v_dS1T, v_dpreS2, v_dS2, v_dS2T,
                v_dpreS3, v_dS3, v_dS3T, v_dpreS4, v_dS4, v_dS4T,
                v_BCD, v_BCD_detect, v_BCD_trt, v_OCD, v_AllD, 
                v_i, v_i_scr, v_i_sym, v_i_S1, v_i_S1_scr, v_i_S1_sym, 
                v_i_S2, v_i_S2_scr, v_i_S2_sym, v_i_S3, v_i_S3_scr, 
                v_i_S3_sym, v_i_S4, v_i_S4_scr, v_i_S4_sym, 
                v_trt, v_dtc, v_N)      #Combine results into a single vector dx
      list(v_dx)                        #Return result as a list
      

    }           
    
  )
  
}

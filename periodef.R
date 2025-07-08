library(dplyr); library(haven); library(rlang)

# function to get periodontal data from the NHANES website from the 2009-10, 2011-12 or
# 2013-14 paired years
# year : year in which the NHANES survey was conducted, takes the values 2009, 2011 or 2013
# site : (optional) vector denoting the names of the six sites of each tooth
nhanes_PD = function(year = 2011, site = NULL){
  
  website = paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", year, "/DataFiles/")
  
  if(year == 2009){
    ext = "_F.XPT"
  } else if(year == 2011){
    ext = "_G.XPT"
  } else if(year == 2013){
    ext = "_H.XPT"
  }
  
  # read the data from NHANES website
  perio = read_xpt(paste0(website, "OHXPER", ext))
  n_dat = nrow(perio)
  
  # The respective tooth codes for 'DF', 'MdF','MF','DL', 'MdL', 'ML' 
  # as given in the NHANES data
  code_sites = c('D', 'M', 'S', 'P', 'L', 'A')
  n_site = length(code_sites)
  
  # construct the tooth vector
  # third molars (1, 16, 17, 31) are excluded in NHANES data
  tooth = c(paste0("0", 2:9), paste0(c(10:15, 18:31)))
  n_tooth = length(tooth)
  
  # PD - pocket depth
  names_PD = rep(paste0("OHX", tooth, "PC"), each = n_site)
  names_PD = paste0(names_PD, rep(code_sites, times = n_tooth))
  mat_PD = perio[ , names_PD]
  
  # CAL - clinical attachment loss
  names_CAL = rep(paste0("OHX", tooth, "LA"), each = n_site)
  names_CAL = paste0(names_CAL, rep(code_sites, times = n_tooth))
  mat_CAL = perio[ , names_CAL]
  
  # if the site names are not mentioned, use the default names
  if(is.null(site)){
    site = c('DF', 'MdF', 'MF', 'DL', 'MdL', 'ML')
  }
  
  # construct the tooth-site level data
  df = tibble(SEQN = rep(perio$SEQN, each = n_site*n_tooth),
              tooth = rep(tooth, each = n_site, times = n_dat), 
              site = rep(site, times = n_dat*n_tooth), 
              PD = as.vector(t(mat_PD)), CAL = as.vector(t(mat_CAL))) 
  
  # group the data by SEQN, tooth 
  df = df %>%
    group_by(SEQN, tooth)
  df$tooth = as.numeric(df$tooth)
  
  return(df)
}


# function to extract variables either from a dataframe or from the R environment
# x : the variable that should be extracted
# data : (optional) data frame, if the variable is a column from the dataframe
extract_variable = function(x, data = NULL) {
  
  # tidy evaluation for variable extraction in data or parent frame
  x_quo = enquo(x)
  if (!is.null(data)) {
    return(eval_tidy(x_quo, data))
  } else {
    return(eval_tidy(x_quo, parent.frame()))
  }
}


# function to construct the tooth-site level dataframe from the user-input
# SEQN : unique sequence number for each subject
# tooth : variable denoting the tooth for each subject
# site : variable denoting the site for each tooth
# PD : value of measured pocket depth (in mm) for each tooth-site
# CAL : value of measured clinical attachment loss (in mm) for each tooth-site
# data : (optional) dataframe containing the above variables
get_data = function(SEQN, tooth, site, PD, CAL, data = NULL) {
  
  # extract the variables using the extract_variable() function
  SEQN.val = extract_variable({{SEQN}}, data)
  tooth.val = extract_variable({{tooth}}, data)
  site.val = extract_variable({{site}}, data)
  PD.val = extract_variable({{PD}}, data)
  CAL.val = extract_variable({{CAL}}, data)
  
  # store the data as a tibble 
  df = tibble(SEQN = SEQN.val, 
              tooth = as.numeric(tooth.val), 
              site = site.val, 
              PD = PD.val, 
              CAL = CAL.val)
  
  # group the data by SEQN, tooth 
  df = df %>%
    group_by(SEQN, tooth)
  
  return(df)
}


# function to check for pairs of non-adjacent teeth
# x : vector containing the teeth for a subject
check_non_ajc = function(x){
  
  n = length(x)
  
  if(n < 2)
    return(FALSE)
  else {
    # Check all unique pairs (i < j)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        if (abs(x[i] - x[j]) > 1) return(TRUE)
      }
    }
    return(FALSE)
  }
}


# function to count the number of opposing pairs of teeth
# teeth_vector : vector containing the teeth for a subject
# exclude_third_molars : FALSE/TRUE denoting whether third molars are present or not
count_opposing_pairs = function(teeth_vector, exclude_third_molars = TRUE) {
  
  pairs = 0
  
  if(exclude_third_molars == TRUE)
    upper_teeth = 2:15
  else upper_teeth = 1:16
  
  for (u in upper_teeth) {
    # for each tooth in the upper jaw, get the corresponding tooth in the lower jaw
    # 1->31, 2->31, 3->30, 4->29, ..., 15->18, 16->17 if including third molars
    # 2->31, 3->30, 4->29, ..., 15->18 if excluding third molars
    l = 33 - u  
    
    # check for the presence of pairs
    if (u %in% teeth_vector && l %in% teeth_vector) {
      pairs = pairs + 1
    }
  }
  return(pairs)
}


# function to get the periodontal status classification according to the 2018 definition
# SEQN : unique sequence number for each subject
# tooth : variable denoting the tooth for each subject
# site : variable denoting the site for each tooth
# PD : value of measured pocket depth (in mm) for each tooth-site
# CAL : value of measured clinical attachment loss (in mm) for each tooth-site
# data : (optional) dataframe containing the above variables
# interdental_sites : vector specifying the site names which are interdental
# buccal_oral_sites : vector specifying the site names which are buccal or oral
# exclude_third_molars : FALSE/TRUE denoting whether third molars are present or not
# coded_values : (optional) vector specifying PD and/or CAL measurements that could not be assessed 
#               and may be coded (such as 99 in NHANES)
periostage = function(SEQN, tooth, site, PD, CAL, data = NULL, 
                      interdental_sites, buccal_oral_sites,
                      exclude_third_molars = TRUE, coded_values = NULL){
  
  # construct the dataframe using user-input
  df = get_data(SEQN = {{SEQN}}, tooth = {{tooth}}, site = {{site}}, 
                PD = {{PD}}, CAL = {{CAL}}, data = data)
  
  # extract the SEQN for each subject in the data
  out = tibble(SEQN = unique(df$SEQN))
  
  #################################################################################
  # Check for Periodontitis Case Definition Assessable
  # >= 2 teeth with buccal/oral CAL and PD measurements 
  #                     OR
  # >= 2 non-adjacent teeth with interdental CAL measurements 
  #################################################################################
  
  # Check for >= 2 teeth with buccal/oral CAL and PD measurements
  
  buccal_oral_dat = df %>%
    filter(site %in% buccal_oral_sites, 
           !is.na(CAL), !(CAL %in% coded_values), 
           !is.na(PD), !(PD %in% coded_values))
  
  buccal_oral_seqn = buccal_oral_dat %>%
    group_by(SEQN) %>%
    filter(length(unique(tooth)) >= 2) %>%
    pull(SEQN) %>%
    unique()
  
  # Check for >= 2 non adjacent interdental CAL measurements
  interdental_dat = df %>%
    filter(site %in% interdental_sites,
           !is.na(CAL), !(CAL %in% coded_values)) 
  
  interdental_seqn = interdental_dat %>%
    group_by(SEQN) %>%
    filter(check_non_ajc(sort(unique(tooth)))) %>%
    pull(SEQN)%>%
    unique()
  
  buccal_interdental_seqn = sort(unique(c(interdental_seqn, buccal_oral_seqn)))
  
  out$assess = ifelse(out$SEQN %in% buccal_interdental_seqn, 'Yes', 'No')
  
  #################################################################################
  # Check for Periodontitis Case 
  # >= 2 teeth with buccal/oral CAL >= 3 mm and PD > 3 mm
  #                     OR
  # >= 2 non-adjacent teeth with interdental CAL >= 1 mm
  #################################################################################
  
  # Check for >= 2 buccal or oral CAL >= 3mm and PD > 3 mm
  
  buccal_oral_PD = buccal_oral_dat %>%
    filter(SEQN %in% buccal_oral_seqn, 
           !is.na(CAL), !(CAL %in% coded_values), 
           !is.na(PD), !(PD %in% coded_values),
           CAL >= 3, PD > 3) 
  
  buccal_oral_PD_seqn = buccal_oral_PD %>%
    group_by(SEQN) %>%
    filter(length(unique(tooth)) >= 2) %>%
    pull(SEQN) %>%
    unique()
  
  # Check for >= 2 non adjacent interdental CAL >= 1mm 
  
  interdental_PD = interdental_dat %>%
    filter(SEQN %in% interdental_seqn, 
           !is.na(CAL), !(CAL %in% coded_values), CAL >= 1) 
  
  interdental_PD_seqn = interdental_PD %>%
    group_by(SEQN) %>%
    filter(check_non_ajc(sort(unique(tooth)))) %>%
    pull(SEQN)%>%
    unique()
  
  # collecting SEQN that satisfy either the interdental or buccal/oral condition
  buccal_interdental_PD_seqn = sort(unique(c(interdental_PD_seqn, buccal_oral_PD_seqn)))
  
  out$case = NA
  out$case[out$assess == 'Yes'] = ifelse(buccal_interdental_seqn %in% buccal_interdental_PD_seqn,
                                         'Yes', 'No')
  
  #################################################################################
  # Determine the stage of Periodontitis Cases 
  #################################################################################
  
  out$stage = NA
  
  # Stage 0 refers to subjects for which the periodontitis case criterion are not fulfilled
  # meaning they are not classified as having periodontitis
  out$stage[out$case == 'No'] = 0 
  
  # create a tooth-site level dataframe for the PD cases
  df_PD_stage = df %>%
    filter(SEQN %in% buccal_interdental_PD_seqn)
  
  # Check for Stage 1
  # Maximum interdental CAL 1 - 2 mm
  PD_stage1_seqn = df_PD_stage %>%
    filter(site %in% interdental_sites,
           !is.na(CAL), !(CAL %in% coded_values)) %>%
    group_by(SEQN) %>%
    filter(max(CAL) >= 1, max(CAL) <= 2) %>%
    pull(SEQN) %>%
    unique()
  
  out$stage[out$SEQN %in% PD_stage1_seqn] = 1
  
  # Check for Stages 2 and 3
  # Maximum interdental CAL 3 - 4 mm
  PD_stage23_seqn = df_PD_stage %>%
    filter(site %in% interdental_sites,
           !is.na(CAL), !(CAL %in% coded_values)) %>%
    group_by(SEQN) %>%
    filter(max(CAL) >= 3, max(CAL) <= 4) %>%
    pull(SEQN) %>%
    unique()
  
  # Stage = 3, if PD >= 6 mm at >=2 non adjacent teeth
  # Stage = 2, otherwise
  nonadjc_PD_stage23 = df_PD_stage %>%
    filter(SEQN %in% PD_stage23_seqn, 
           !is.na(PD), !(PD %in% coded_values),
           PD >=6) 
  
  nonadjc_PD_stage23_seqn = nonadjc_PD_stage23 %>%
    group_by(SEQN) %>%
    filter(check_non_ajc(sort(unique(tooth)))) %>%
    pull(SEQN) %>%
    unique()
  
  if(length(nonadjc_PD_stage23_seqn)==0){
    nonadjc_PD_stage23_seqn = NULL
  } 
  
  PD_stage2_seqn = PD_stage23_seqn[!(PD_stage23_seqn %in% nonadjc_PD_stage23_seqn)]
  
  out$stage[out$SEQN %in% PD_stage2_seqn] = 2
  out$stage[out$SEQN %in% nonadjc_PD_stage23_seqn] = 3
  
  # Check for Stages 3 and 4
  # Maximum interdental CAL >= 5 mm 
  PD_stage34_seqn = df_PD_stage %>%
    filter(site %in% interdental_sites,
           !is.na(CAL), !(CAL %in% coded_values)) %>%
    group_by(SEQN) %>%
    filter(max(CAL) >= 5) %>%
    pull(SEQN)%>%
    unique()
  
  PD_stage34 = df_PD_stage %>%
    filter(SEQN %in% PD_stage34_seqn) 
  
  # Stage = 4, if # opposing pairs of natural teeth < 10
  PD_stage4_seqn = PD_stage34 %>%
    group_by(SEQN) %>%
    filter(!is.na(CAL), !is.na(PD)) %>%
    filter(count_opposing_pairs(teeth_vector = tooth, 
                                exclude_third_molars = exclude_third_molars) < 10) %>%
    pull(SEQN)%>%
    unique()
  
  PD_stage3_seqn = PD_stage34_seqn[!(PD_stage34_seqn %in% PD_stage4_seqn)]
  out$stage[out$SEQN %in% PD_stage3_seqn] = 3
  out$stage[out$SEQN %in% PD_stage4_seqn] = 4  
  
  return(out)
}
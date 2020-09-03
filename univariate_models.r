require(gamm4)
require(plyr)
require(parallel)
require(stringr)
require(MuMIn)
require(pracma)
require(tableone)
library(dplyr)
library(broom)
library(tidyr)
library(feather)

gamma_variables = c('bisbas_ss_bas_drive',
'cbcl_aggressive',
'cbclulebreak',
'cbcl_attention',
'cbcl_social',
'cbcl_thought',
'cbcl_somatic',
'cbcl_withdep',
'cbcl_anxdep',
'mania',
'prosociality_p',
'prosociality_y',
'psychosis_severity_score', 
'upps_ss_lack_of_perseverance',
'upps_ss_positive_urgency',
"KSADS_depression_S_p",
"KSADS_bipolar_S_p",           
"KSADS_psychosis_S_p",          
"KSADS_anxiety_S_p",           
"KSADS_ocd_S_p",                
"KSADS_eatdis_S_p",            
"KSADS_adhd_S_p",               
"KSADS_oppcond_S_p",           
"KSADS_devdis_S_p",             
"KSADS_ptsd_S_p",              
"KSADS_insomnia_S_p",           
"KSADS_suicidality_S_p",       
"KSADS_homicidality_S_p",       
"KSADS_ANY_S_p",               
"KSADS_depression_S_y",         
"KSADS_bipolar_S_y",           
"KSADS_anxiety_S_y",            
"KSADS_insomnia_S_y",          
"KSADS_suicidality_S_y",        
"KSADS_ANY_S_y")    
    # KSADS_symptoms)

binomial_variables = c('KSADS_ANY_D_y', 'KSADS_ANY_D_p')

# Helper function to reformat ids from PRS file
reformat_id <- function(id){
    ans = unlist(strsplit(id, 'NDAR_'))[2]
    ans = paste0('NDAR_', unlist(strsplit(ans, '_'))[1])
    return(ans)
}

read_in_ABCD = function(nda_file=character(), prs_file=character(), pcs_file=character(), event_filt='baseline_year_1_arm_1'){
    if (!endsWith(nda_file, '.feather')){
      # Read in nda file
      df = readRDS(nda_file)
    }else{
      df = read_feather(nda_file)
    }
    df = df[df$eventname==event_filt, ]

    # Read in PRS if they are given 
    if (length(prs_file) != 0){
        # Read Polygenic Scores
        prs = read.table(prs_file, header=TRUE)
        colnames(prs)[1] = 'src_subject_id'
        reformated_rows = lapply(as.character(prs$src_subject_id), reformat_id)
        prs = prs[!duplicated(reformated_rows),]
        prs$src_subject_id = reformated_rows[!duplicated(reformated_rows)]
        # Concatonate _PRS to cols
        # colnames(prs)[2:dim(prs)[2]] = paste0(colnames(prs)[2:dim(prs)[2]], '_PRS')
        # Join with data frame
        df = join(df, prs)
    }

    # Read in Genetic PCs if they are given 
    if (length(pcs_file) != 0){
        # Read PCs
        pcs = read.table(pcs_file, header=TRUE)
        colnames(pcs)[colnames(pcs)=='IID'] = 'src_subject_id'
        df = join(df, pcs, by='src_subject_id')
    }
    return(df)
}

save_table_1 = function(df, outfile, na_rm_vars){
    varList = c('age', 'sex', 'high.educ', 'household.income', 'race_ethnicity')
    ind = rowSums(is.na(df[, c(na_rm_vars)]))==0
    tableone = CreateTableOne(vars = varList, data = df[ind, ])
    tab1Mat <- print(tableone, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
    ## Save to a CSV file
    write.table(tab1Mat, outfile, quote=FALSE, sep='\t')
}

create_model_table = function(DVs, 
IVs, 
covariates, 
random_effects = '~(1|abcd_site/rel_family_id)', 
reduced_models=TRUE,
combined_models=FALSE){

    # Create Model Table
    covariate_str = paste(covariates, collapse=' + ')
    if (!combined_models){
        model_table = data.frame(DV=rep(DVs, length(IVs)), 
            IV=rep(IVs, each=length(DVs)))
        model_table$full = paste0(model_table$DV, ' ~ ', model_table$IV, ' + ', covariate_str)
    }else{
        model_table = data.frame(DV=DVs)
        model_table[, paste0('IV_', seq(length(IVs)))] = t(replicate(length(DVs), IVs))
        model_table[, 'full'] = paste0(model_table$DV, ' ~ ', paste0(IVs, collapse=' + '), ' + ', covariate_str)
    }

    if (reduced_models){
            model_table$reduced = paste0(model_table$DV, ' ~ ', covariate_str)
    }
    if (random_effects != FALSE){
        model_table$random_effects = random_effects
    }

    gamma_ind = model_table$DV %in% gamma_variables
    binomial_ind = model_table$DV %in% binomial_variables
    model_table$distribution = 'normal'
    model_table$distribution[gamma_ind] = 'gamma'
    model_table$distribution[binomial_ind] = 'binomial'
    return(model_table)
}

gamm4_link = function(form, distribution, random, dat){
    print(paste0('Fitting ', form))
    if (distribution=='normal'){
        return <- gamm4(formula=formula(form), random=formula(random), data=dat)
    }else if (distribution=='gamma'){
        return <- gamm4(formula=formula(form), random=formula(random), data=dat, family=Gamma(link = "log"))
    }else if (distribution=='binomial'){
        return <- gamm4(formula=formula(form), random=formula(random), data=dat, family=binomial())
    }
}

glm_link = function(form, distribution, dat){
  print(paste0('Fitting ', form))
  tryCatch({
  if (distribution=='normal'){
    model <- glm(formula(form), family=gaussian(link="identity"), data=dat, control=list(maxit=10000))
  }else if (distribution=='gamma'){
    model <- glm(formula(form), family=Gamma(link="log"), data=dat, control=list(maxit=10000))
  }else if (distribution=='binomial'){
    model <- glm(formula(form), family=binomial(link = "logit"), data=dat, control=list(maxit=10000))
  }
  return(model)
  }, warning = function(war) {
      print(paste("warning:  ",war))
      return(model)
  }, error = function(err) {
      print(paste("error:  ",err))
      return(NULL)
  }, finally = {

  })
}

extract_reg_results <- function(full_models, reduced_models, model_table, model_type, suffix=''){
  if (model_type=='GAMM4'){
    for (i in 1:dim(model_table)[1]){
        sumry = summary(full_models[,i]$gam)
        model_table[i, paste0('coeff', suffix)] = as.numeric(sumry$p.coeff[2])
        model_table[i, paste0('coeff_SE', suffix)] = as.numeric(sumry$se[2])
        model_table[i, paste0('upper', suffix)] = as.numeric(sumry$p.coeff[2]) + 1.96*as.numeric(sumry$se[2])
        model_table[i, paste0('lower', suffix)] = as.numeric(sumry$p.coeff[2]) - 1.96*as.numeric(sumry$se[2])
        DF = as.numeric(sumry$residual.df)
        t =  as.numeric(sumry$p.t[2])
        model_table[i, paste0('t', suffix)] = t
        model_table[i, paste0('r2_tstat', suffix)] = (t^2/(t^2 + DF))
        model_table[i, paste0('pval', suffix)] = sumry$p.pv[[2]]
        model_table[i, paste0('r2_full', suffix)] = sumry$r.sq
        # Likelihood ratio test
        model_table[i, paste0('n_full', suffix)] = length(full_models[,i]$gam$residuals)
        model_table[i, paste0('full_converge_flag', suffix)] = as.numeric(isempty(full_models[,i]$mer@optinfo$warnings))
        if (!is.null(reduced_models)){
            model_table[i, paste0('r2_delta', suffix)] = r.squaredLR(full_models[,i]$mer, reduced_models[,i]$mer)
            lrttest_table = anova(full_models[,i]$mer, reduced_models[,i]$mer)
            model_table[i, paste0('lrt_chisq', suffix)] = lrttest_table['full_models[, i]$mer', 'Chisq']
            model_table[i, paste0('lrt_pval', suffix)] = lrttest_table['full_models[, i]$mer', 'Pr(>Chisq)']
            model_table[i, paste0('n_reduced', suffix)] = length(reduced_models[,i]$gam$residuals)
            model_table[i, paste0('reduced_converge_flag', suffix)] = as.numeric(isempty(reduced_models[,i]$mer@optinfo$warnings))
        }
    }
  }else if (model_type=='glm'){
    emptylist<-sapply(X=full_models, FUN=isempty)
    if (!isempty(which(emptylist==FALSE))){ #at least one model converged - if no models converged will return model_table same as input
      if (is.null(reduced_models)){
        NAmat<-matrix(NA,nrow = dim(model_table)[1],ncol = 9)
        colnames(NAmat)<-c(paste0('coeff', suffix),paste0('coeff_SE', suffix),paste0('upper', suffix),paste0('lower', suffix),paste0('t', suffix),paste0('r2_tstat', suffix),paste0('pval', suffix),paste0('n_full', suffix),paste0('full_converge_flag', suffix))
      } else if (!is.null(reduced_models)){
        NAmat<-matrix(NA,nrow = dim(model_table)[1],ncol = 14)
        colnames(NAmat)<-c(paste0('coeff', suffix),paste0('coeff_SE', suffix),paste0('upper', suffix),paste0('lower', suffix),paste0('t', suffix),paste0('r2_tstat', suffix),paste0('pval', suffix),paste0('n_full', suffix),paste0('full_converge_flag', suffix),
                           paste0('r2_delta', suffix),paste0('lrt_deviance', suffix),paste0('lrt_pval', suffix),paste0('n_reduced', suffix),paste0('reduced_converge_flag', suffix))
      }
      model_table<-cbind(model_table,NAmat)
      for (i in 1:dim(model_table)[1]){
        if (!is.null(full_models[[i]])){ #if model converged then read out stats
          IV = as.character(model_table[i, 'IV'])
          sumry = summary(full_models[[i]])
          model_table[i, paste0('coeff', suffix)] = as.numeric(sumry$coef[IV,1])
          model_table[i, paste0('coeff_SE', suffix)] = as.numeric(sumry$coef[IV,2])
          model_table[i, paste0('upper', suffix)] = as.numeric(sumry$coef[IV,1]) + 1.96*as.numeric(sumry$coef[IV,2])
          model_table[i, paste0('lower', suffix)] = as.numeric(sumry$coef[IV,1]) - 1.96*as.numeric(sumry$coef[IV,2])
          t = as.numeric(sumry$coef[IV,3])
          DF = as.numeric(sumry$df.residual)
          model_table[i, paste0('t', suffix)] = t
          model_table[i, paste0('r2_tstat', suffix)] = (t^2/(t^2 + DF))
          model_table[i, paste0('pval', suffix)] = sumry$coef[IV,4]
          model_table[i, paste0('n_full', suffix)] = length(full_models[[i]]$residuals)
          model_table[i, paste0('full_converge_flag', suffix)] = as.numeric(full_models[[i]]$converged)
          if (!is.null(reduced_models)){
            model_table[i, paste0('r2_delta', suffix)] = r.squaredLR(full_models[[i]], reduced_models[[i]])
            # Likelihood ratio test
            lrttest_table = anova(full_models[[i]], reduced_models[[i]], test="LRT")
            model_table[i, paste0('lrt_deviance', suffix)] = lrttest_table[2,4]
            model_table[i, paste0('lrt_pval', suffix)] = lrttest_table[2,5]
            model_table[i, paste0('n_reduced', suffix)] = length(reduced_models[[i]]$residuals)
            model_table[i, paste0('reduced_converge_flag', suffix)] = as.numeric(reduced_models[[i]]$converged)
          }
        }
      }
    }
  }
    return(model_table)
  }

purrr_extract <- function(full_models, model_table, model_type, suffix=''){
  if (model_type=='GAMM4'){
  var_names = names(summary(full_models[,1]$gam)$p.coeff) #assumes all models have same covariates
  NAmat<-matrix(NA,nrow = dim(model_table)[1],ncol = length(var_names)*5)
  colnames(NAmat)<-c(paste0(var_names, '_coeff', suffix),paste0(var_names, '_coeffSE', suffix),paste0(var_names, '_t', suffix), paste0(var_names, '_r2_tstat', suffix) ,paste0(var_names, '_pval', suffix))
  model_table<-cbind(model_table,NAmat)
  for (i in 1:dim(model_table)[1]){
      sumry = summary(full_models[,i]$gam)
      model_table[i, paste0(var_names, '_coeff', suffix)] = as.numeric(sumry$p.coeff)
      model_table[i, paste0(var_names, '_coeffSE', suffix)] = as.numeric(sumry$se)
      #model_table[i, paste0(var_names, '_upper', suffix)] = as.numeric(sumry$p.coeff) + 1.96*as.numeric(sumry$se)
      #model_table[i, paste0(var_names, '_lower', suffix)] = as.numeric(sumry$p.coeff) - 1.96*as.numeric(sumry$se)
      ts = as.numeric(sumry$p.t)
      DF = as.numeric(sumry$residual.df)
      model_table[i, paste0(var_names, '_t', suffix)] = ts
      model_table[i, paste0(var_names, '_r2_tstat', suffix)] = (ts^2/(ts^2 + DF))
      model_table[i, paste0(var_names, '_pval', suffix)] = sumry$p.pv
      model_table[i, paste0('n_full', suffix)] = length(full_models[,i]$gam$residuals)
      model_table[i, paste0('r2_full', suffix)] = sumry$r.sq
      model_table[i, paste0('full_converge_flag', suffix)] = as.numeric(isempty(full_models[,i]$mer@optinfo$warnings))
    }
  }else if (model_type=='glm'){
    emptylist<-sapply(X=full_models, FUN=isempty)
    if (!isempty(which(emptylist==FALSE))){
      firstind<-min(which(emptylist==FALSE))
    var_names = rownames(summary(full_models[[firstind]])$coef) #assumes all models have same covariates
    NAmat<-matrix(NA,nrow = dim(model_table)[1],ncol = length(var_names)*5)
    colnames(NAmat)<-c(paste0(var_names, '_coeff', suffix),paste0(var_names, '_coeffSE', suffix),paste0(var_names, '_t', suffix), paste0(var_names, '_r2_tstat', suffix), paste0(var_names, '_pval', suffix))
    model_table<-cbind(model_table,NAmat)
    for (i in 1:dim(model_table)[1]){
      if (!is.null(full_models[[i]])){
        sumry = summary(full_models[[i]])
        model_table[i, paste0(var_names, '_coeff', suffix)] = as.numeric(sumry$coef[,1])
        model_table[i, paste0(var_names, '_coeffSE', suffix)] = as.numeric(sumry$coef[,2])
        #model_table[i, paste0(var_names, '_upper', suffix)] = as.numeric(sumry$coef[,1]) + 1.96*as.numeric(sumry$coef[,2])
        #model_table[i, paste0(var_names, '_lower', suffix)] = as.numeric(sumry$coef[,1]) - 1.96*as.numeric(sumry$coef[,2])
        t = as.numeric(sumry$coef[,3])
        DF = as.numeric(sumry$df.residual)
        model_table[i, paste0(var_names, '_t', suffix)] = t
        model_table[i, paste0(var_names, '_r2_tstat', suffix)] = (t^2)/(t^2 + DF)
        model_table[i, paste0(var_names, '_pval', suffix)] = sumry$coef[,4]
        model_table[i, paste0('n_full', suffix)] = length(full_models[[i]]$residuals)
        model_table[i, paste0('full_converge_flag', suffix)] = as.numeric(full_models[[i]]$converged)
      }
    }
    }
  }
    return(model_table)
    # model_table$coef_table =  full[2, ] %>% map(summary) %>% map(c('coefficients'))
    # model_table$tidied %>% map(~ select(., term, estimate)) %>% map(~ spread(., term, estimate)) %>% map(data.frame) %>% ldply(rbind)
    # model_table$tidied %>% map(~ select(., term, std.error)) %>% map(~ spread(., term, std.error)) %>% map(data.frame) %>% ldply(rbind)
    # model_table$tidied %>% map(~ select(., term, statistic)) %>% map(~ spread(., term, statistic)) %>% map(data.frame) %>% ldply(rbind)
}

quant_norm = function(df_col){
    n = length(df_col)
    Fn = ecdf(df_col)
    return(qnorm(Fn(df_col)-0.5/n))
}

pre_proc_vars = function(df, model_table, QN_IVs=TRUE, QN_gauss_DVs=TRUE){
    gamma_vars = unique(as.character(model_table[model_table$distribution=='gamma', 'DV']))
    binomial_vars = unique(as.character(model_table[model_table$distribution=='binomial', 'DV']))
    gaussian_vars = unique(as.character(model_table[model_table$distribution=='normal', 'DV']))

    ## Print out which distributions are being assumed for which variables
    print('Gamma variables:')
    print(gamma_vars)
    print('Binomial variables:')
    print(binomial_vars)
    if (isempty(gamma_vars)==FALSE){
        # Ensure that all entries are greater than zero for gamma variables
        if (length(gamma_vars)>1){
            min_vals = apply(df[, gamma_vars], 2, FUN=min)
            df[, gamma_vars] = sweep(df[, gamma_vars],  2, min_vals - 0.1)
        }else{
            min_val = min(df[, gamma_vars])
            df[, gamma_vars] = df[, gamma_vars] - (min_val - 0.1)
        }
    }

    if (QN_IVs){
        IV_cols = colnames(model_table)[startsWith(colnames(model_table), 'IV')]    
        IVs = as.character(unique(unlist(model_table[, IV_cols])))
        if (length(IVs)>1){
            IVs_numeric = IVs[sapply(df[, IVs], is.numeric)]
        }else{
            IVs_numeric = IVs[is.numeric(df[, IVs])]
        }
        print('Qunatile normalizing independant variables:')
        print(IVs_numeric)
        if (length(IVs_numeric)>1){
            df[, IVs_numeric] = as.numeric(apply(df[, IVs_numeric], 2, quant_norm))
        }else{
            df[, IVs_numeric] = as.numeric(quant_norm(df[, IVs_numeric]))
        }
    }
    
    if (QN_gauss_DVs & (isempty(gaussian_vars)==FALSE)){
      print('Qunatile normalizing gaussian dependant variables:')
      print(gaussian_vars)
      if (length(gaussian_vars)>1){
        df[, gaussian_vars] = as.numeric(apply(df[, gaussian_vars], 2, quant_norm))
      }else{
        df[, gaussian_vars] = as.numeric(quant_norm(df[, gaussian_vars]))
      }
    }
    # for (binomial_var in binomial_vars){
        
    # }
    return(df)
}

subselect_df = function(df, subset_codition){
    subset_col = strsplit(subset_condition, split="[-+*/)( ]|[^x][0-9]+|^[0-9]+")[[1]][1]
    assign(subset_col, df[, subset_col])
    selection_ind = eval(parse(text=subset_condition))
    return(df[selection_ind, ])
}

run_models = function(model_table, df, model_type='GAMM4', cores=10, QN_IVs=TRUE, QN_gauss_DVs=TRUE, save_all_coeffs=FALSE, return_models=FALSE, suffix=''){
  df = pre_proc_vars(df, model_table, QN_IVs=QN_IVs, QN_gauss_DVs=QN_gauss_DVs)  
  if (model_type=='GAMM4'){
      if (cores>1){
        cores=1
        print('Parallel processing not currently available with gamm4. Number of cores set to 1.')
      }
        
        if ('reduced' %in% colnames(model_table)){unique_reduced = unique(model_table$reduced)}
        if (cores==1){
            if ('reduced' %in% colnames(model_table)){
                reduced = mapply(gamm4_link, form=unique_reduced, distribution=model_table$distribution[1:length(unique_reduced)], random=model_table$random[1:length(unique_reduced)], MoreArgs=list(dat=df))
            }
            full = mapply(gamm4_link, form=model_table$full, distribution=model_table$distribution, model_table$random, MoreArgs=list(dat=df))
        }else{
            if ('reduced' %in% colnames(model_table)){
                reduced = mcmapply(gamm4_link, form=model_table$full, distribution=model_table$distribution[1:length(unique_reduced)], random=model_table$random[1:length(unique_reduced)], mc.cores=cores, MoreArgs=list(dat=df))
            }
            full = mcmapply(gamm4_link, form=model_table$full, distribution=model_table$distribution, random=model_table$random, mc.cores=cores, MoreArgs=list(dat=df))
        }
        if ('reduced' %in% colnames(model_table)){
            repetitions = dim(model_table)[1]/length(unique_reduced)
            reduced = do.call(cbind, replicate(repetitions, reduced, simplify=FALSE))
        }else{
            reduced = NULL
        }
        
    }else if (model_type=='permutation'){
        # Permutation code goes here
    }else if (model_type=='glm'){
      if ('reduced' %in% colnames(model_table)){unique_reduced = unique(model_table$reduced)}
      if (cores==1){
        if ('reduced' %in% colnames(model_table)){
          reduced = mapply(glm_link, form=unique_reduced, distribution=model_table$distribution[1:length(unique_reduced)], MoreArgs=list(dat=df), SIMPLIFY=FALSE)
        }
        full = mapply(glm_link, form=model_table$full, distribution=model_table$distribution, MoreArgs=list(dat=df), SIMPLIFY=FALSE)
      }else{
        if ('reduced' %in% colnames(model_table)){
          reduced = mcmapply(glm_link, form=unique_reduced, distribution=model_table$distribution[1:length(unique_reduced)], mc.cores=cores, MoreArgs=list(dat=df), SIMPLIFY=FALSE)
        }
        full = mcmapply(glm_link, form=model_table$full, distribution=model_table$distribution, mc.cores=cores, MoreArgs=list(dat=df), SIMPLIFY=FALSE)
      }
      if ('reduced' %in% colnames(model_table)){
        repetitions = dim(model_table)[1]/length(unique_reduced)
        reduced = do.call(cbind, replicate(repetitions, reduced, simplify=FALSE))
      }else{
        reduced = NULL
      }
      
    }

    if (save_all_coeffs){
      model_table = purrr_extract(full, model_table, model_type, suffix=suffix)
    }else{
      model_table = extract_reg_results(full, reduced, model_table, model_type, suffix=suffix)
    }

    if (return_models==TRUE){
          if ('reduced' %in% colnames(model_table)){
            output<-list(model_table=model_table,full=full,reduced=reduced)
            return(output)
          } else {
            output<-list(model_table=model_table,full=full)
            return(output)
          }
          
    }else{
        return(model_table)
    }
}

select_singletons<-function(df){
  df<-df[sample(nrow(df)),]
  df<-distinct(df,rel_family_id, .keep_all=TRUE)
  return(df)
}





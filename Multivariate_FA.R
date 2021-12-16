# -----  Multivariate twin models  ----- #

library(umx)
library(tidyverse)
library(rlist)

## -- Load data
data_path = "~/data"
load(paste(data_path,"twin_data.R",sep="/"))

## -- Functions
# - Table with model results
show_models <- function(model) {
  iden=mxCheckIdentification(model)
  stats=summary(model)
  name=model@name
  aic=stats[["AIC.Mx"]]
  m2ll=stats[["Minus2LogLikelihood"]]
  bic=stats[["BIC.Mx"]]
  identified=iden[["status"]]
  non_id_parm=paste(iden[["non_identified_parameters"]],collapse = ",")
  parm_to_drop=iden[["non_identified_parameters"]]
  out = list(data.frame(Model=name,AIC=aic,LogLike=m2ll,BIC=bic,identification=identified,non_id=non_id_parm),
             parm_to_drop)
  return(out)
}
# - Modifies models by dropping non-identified parameters and returns table with model results
modify_models <- function(model_list,df,parm2drop,net) {
  round=2
  df= df %>% add_column(Iterations = 1)
  while (all(df$non_id=="None") == FALSE){
    vec=which(df$non_id != "None")
    sublist=model_list[vec]
    paramstmp=params2drop[vec]
    namestmp=df$Model[vec]
    for (m in 1:length(sublist)){
      paths=as.vector(paramstmp[[m]])
      sublist[[m]]=assign(paste("cor",m,round, sep="_"), umxModify(sublist[[m]], update = paths, name = paste(namestmp[m],paste(paths,collapse = ","))))
    }
    #--- Data frame showing created models
    new_res = map(sublist,show_models)
    
    #--- Visualize res
    modified_info=map(new_res, 1) %>% reduce(rbind) %>% add_column(iterations=round)
    modified_params2drop=map(new_res, 2) 
    
    df[vec,] = modified_info[,]
    params2drop[vec] = modified_params2drop
    model_list[vec] = sublist
    round=round+1
  }
  df = df %>% add_column(Fiber_gr = net, .before = "Model") %>% add_column(ind = c(1:nrow(.)))
  out=list(df,model_list)
  return(out)
}
# - Returns Table with confidence intervals
cis_table <- function(model) {
  ccii=umxConfint(model,parm="smart",run = TRUE)
  table=as.data.frame(ccii@output[["confidenceIntervals"]]) %>% 
    add_column(parameter=row.names(.),.before = "lbound")
  if (NA %in% table$lbound | NA %in% table$ubound) {
    parameters = table$parameter
    estimates = table$estimate
    model=mxRun(model, intervals = TRUE)
    CIdetail=summary(model, verbose= TRUE)$CIdetail[, c("parameter", "value", "side")]
    tr=CIdetail %>% 
      filter(str_detect(parameter, "1,1|2,2|3,3|4,4|5,5|6,6|7,7|8,8|9,9|10,10|cp_"))
    
    table_com = cbind(tr %>% filter(side == "lower") %>% select(all_of(c("parameter","value"))),
                      tr %>% filter(side == "upper") %>% select(value)) %>%
      setNames(c("parameter","lbound","ubound")) %>% 
      mutate_at(.vars = c("parameter"), ~ ifelse(.=="a_cp_r1c1",
                                                 "top.a_cp[1,1]",ifelse(.=="c_cp_r1c1",
                                                                        "top.a_cp[1,1]",ifelse(.=="e_cp_r1c1",
                                                                                               "top.e_cp[1,1]",.)))) %>%
      filter(parameter %in% parameters) %>% 
      arrange(match(parameter, table$parameter)) %>% 
      add_column(estimate=estimates, .before = "ubound")
    
    table = table_com
    
  }
  table = as_tibble(table) %>% mutate_if(is.numeric,round,2) %>% mutate_at(.vars=c("lbound","ubound"), as.character) %>%
    mutate(cis = paste("(",lbound,"-",ubound,")",sep=""))
  
  return(table)
}
# - Set optimizer
umx_set_optimizer("SLSQP")

### --- Common pathway models --- ###

## -- Define models
# - Model names
names=c("Drop cg cs", "Drop ag as", "Drop cs", "Drop as", "Drop cg", "Drop ag", "Drop cg cs as")
# - Model syntax
models=c("(^cs_)|(^c_cp_)","(^as_)|(^a_cp_)","(cs_)","(as_)","(^c_cp_)","(^a_cp_)","(^cs_)|(^c_cp_)|(as_)")


### --- Make models for all fiber tract groups

tracing=c(rep("func",4),rep("atlas",3))
groups=list(gl_fibers,fa_fibers,core_fibers,ext_fibers,fa_fibers,core_fibers,ext_fibers)
var_names=c("Global","Face","CFN","EFC","Face_atlas","CFN_atlas","EFC_atlas")

for (c in 1:length(tracing)) {
  Models_list = list()
  group=groups[[c]]
  tr=tracing[c]
  var_name=var_names[c]
  print(paste("---------------------------", var_name, "---------------------------",sep=" "))
  if (tr=="func"){
    cp_1 = umxCP("All_paths", selDVs = group, sep = "_T", nFac = 1, dzData = DZ, mzData = MZ)
  }
  else{
    cp_1 = umxCP("All_paths", selDVs = group, sep = "_T", nFac = 1, dzData = DZ_atlas, mzData = MZ_atlas)
  }
  # - Add to list
  Models_list=Models_list %>% list.append(cp_1)
  # - Fit all models
  for (m in 1:length(models)){
    model = umxModify(cp_1, regex = models[m], name = names[m])
    Models_list=Models_list %>% list.append(model)
  }
  # - Inspect results
  models_output = map(Models_list,show_models)
  models_info=map(models_output, 1) %>% reduce(rbind)
  params2drop=map(models_output, 2)
  assign(paste("model_info",var_name,sep="_"),models_info)
  assign(paste("model_list",var_name,sep="_"),Models_list)
  
  # - Modify models
  print(paste("---------------------------", "Modify models", "---------------------------",sep=" "))
  out=modify_models(Models_list,models_info,parm2drop,var_name)
  models_info_final = out[[1]]
  Models_list_final = out[[2]]
  assign(paste("models_info_final",var_name,sep="_"),models_info_final)
  assign(paste("model_list_final",var_name,sep="_"),Models_list_final)
}

# - Joint table with fitted models
model_summary=rbind(models_info_final_Global,models_info_final_Face,
                   models_info_final_CFN,models_info_final_EFC,
                   models_info_final_Face_atlas,models_info_final_CFN_atlas,models_info_final_EFC_atlas) 
  
# - Joint table with final models 
final_models = model_summary %>%
  group_by(Fiber_gr) %>% 
  filter(AIC == min(AIC)) %>% distinct(AIC, .keep_all = T)

## -- Estimate confidence intervals for final models
all_mods=list(model_list_final_Global[[4]],model_list_final_Face[[5]],
              model_list_final_CFN[[2]],model_list_final_EFC[[8]],
              model_list_final_Face_atlas[[4]],model_list_final_CFN_atlas[[6]],
              model_list_final_EFC_atlas[[2]])

CIs=map(all_mods,cis_table)


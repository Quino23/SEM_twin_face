
# -----  Univariate twin models  ----- #

library(umx)
library(tidyverse)
library(rlist)
library(cowplot)

## -- Load data
data_path = "~/data"
load(paste(data_path,"twin_data.R",sep="/")) 


## -- Functions and setup
# - Output model parameter estimates, model fit estimates and confidence intervals
output <- function(model,dig) {
  #model=mxRun(model, intervals = TRUE)
  model=mxTryHard(model,intervals = TRUE)
  # - Estimates
  name=model@name
  stats=summary(model)
  bic=stats[["BIC.Mx"]]
  aic=stats[["AIC.Mx"]]
  fitLogLike=model@fitfunction$result[1,1]
  a=model@output$algebras$top.a_std[1,1]
  c=model@output$algebras$top.c_std[1,1]
  e=model@output$algebras$top.e_std[1,1]
  # - CIS
  cis=model@output[["confidenceIntervals"]] %>% as.vector()
  cis=cis[-c(4,5,6)]
  if (NA %in% cis){
    CIs_to_set = names(omxGetParameters(model, free = TRUE))
    model_summary = summary(mxModel(model, mxCI(CIs_to_set, interval = 0.95)), verbose = TRUE)
    CIdetail = model_summary$CIdetail[, c("parameter", "value", "side")]
    cis_v = CIdetail[,"value"] %>% as.vector()
    cis_match_order = c(cis_v[1],cis_v[3],cis_v[5],
                        cis_v[2],cis_v[4],cis_v[6])
    
    pos=which(is.na(cis))
    cis[pos] = cis_match_order[pos]
  }
  cis.ch=cis %>% round(dig) %>% as.character()
  aci = paste("(",cis.ch[1],"-",cis.ch[4],")",sep="")
  cci = paste("(",cis.ch[2],"-",cis.ch[5],")",sep="")
  eci = paste("(",cis.ch[3],"-",cis.ch[6],")",sep="")   
  
  # - se CIS
  if (length(model@output[["standardErrors"]]) == 4){
    aSE=model@output[["standardErrors"]][2]
    cSE=model@output[["standardErrors"]][3]
    eSE=model@output[["standardErrors"]][4]
  } 
  else if (length(model@output[["standardErrors"]]) == 3) {
    if (model@output$algebras$top.a_std[1,1]==0) {
      aSE=0
      cSE=model@output[["standardErrors"]][2]
      eSE=model@output[["standardErrors"]][3]
    } else {
      aSE=model@output[["standardErrors"]][2]
      cSE=0
      eSE=model@output[["standardErrors"]][3]}
  }
  else {
    aSE=0
    cSE=0
    eSE=model@output[["standardErrors"]][2]
  }
  
  aupper=a + (1.96 * aSE)
  alower=a - (1.96 * aSE)
  cupper=c + (1.96 * cSE)
  clower=c - (1.96 * cSE)
  eupper=e + (1.96 * eSE)
  elower=e - (1.96 * eSE)
  
  se.cis=c(alower,clower,elower,aupper,cupper,eupper)
  se.cis[se.cis>1]=1
  se.cis[se.cis<0]=0
  
  se.cis.ch=se.cis %>% round(dig) %>% as.character()
  se.aci = paste("(",se.cis.ch[1],"-",se.cis.ch[4],")",sep="")
  se.cci = paste("(",se.cis.ch[2],"-",se.cis.ch[5],")",sep="")
  se.eci = paste("(",se.cis.ch[3],"-",se.cis.ch[6],")",sep="") 
  
  
  res=data.frame(Model=name,
                 a=a,aci_l=cis[1],aci_u=cis[4],aci=aci,
                 se.aci_l=se.cis[1],se.aci_u=se.cis[4],se.aci=se.aci,
                 
                 c=c,cci_l=cis[2],cci_u=cis[5],cci=cci,
                 se.cci_l=se.cis[2],se.cci_u=se.cis[5],se.cci=se.cci,
                 
                 e=e,eci_l=cis[3],eci_u=cis[6],eci=eci,
                 se.eci_l=se.cis[3],se.eci_u=se.cis[6],se.eci=se.eci,
                 
                 AIC=aic,LogLike=fitLogLike,BIC=bic)
  
  print(res)
}

# - Set proper fiber tract names for figures
name_changer=function(x,tracts,tracts_new){
  x$Fiber=as.character(x$Fiber)
  for (i in 1:length(tracts)){
    name=tracts[i]
    new_name=tracts_new[i]
    pos=x$Fiber==name
    x$Fiber[pos]=new_name
  }
  x$Fiber=as.factor(x$Fiber)
  return(x)
}

# - set optimizer
umx_set_optimizer("SLSQP")

### --- Univariate models for Global fiber tracts and face fibers traced following functional localization

## - Lists to store models
Models_ace=list()
Models_ae=list()
Models_ce=list()
Models_e=list()

## - Make models
for(i in 1:length(fibers)){
  # -  Current fiber tract
  var = fibers[i]
  # - Make competing models
  mod_ace=umxACE(selDVs = var, dzData = DZ, mzData = MZ, sep = "_T", autoRun = F)
  mod_ae=umxModify(mod_ace, update = "c_r1c1", name = "AE", autoRun = F)
  mod_ce=umxModify(mod_ace, update = "a_r1c1", name = "CE", autoRun = F)
  mod_e=umxModify(mod_ae, update = "a_r1c1", name = "E",    autoRun = F)
  # - Pack models into list
  Models_ace=Models_ace %>% list.append(mod_ace)
  Models_ae=Models_ae %>% list.append(mod_ae)
  Models_ce=Models_ce %>% list.append(mod_ce)
  Models_e=Models_e %>% list.append(mod_e)
}

results_ace = map(Models_ace,output,2) %>% reduce(rbind)
results_ae = map(Models_ae,output,2) %>% reduce(rbind)
results_ce = map(Models_ce,output,2) %>% reduce(rbind)
results_e = map(Models_e,output,2) %>% reduce(rbind)

# - Model parameter table
net=c(rep("Global",10),rep("Core",3),rep("Ext",4))
univariate_func=rbind(results_ace,results_ae,results_ce,results_e) %>%
  add_column(Fiber=rep(fibers,4), .before="Model") %>%
  add_column(Net=rep(net,4)) %>% add_column(Net_gen=rep(c(rep("Global",10),rep("Face",7)),4))  %>%
  mutate_at(.vars = "Fiber", ~ (factor(., levels = fibers)))

### --- Univariate models for face fibers traced following atlas ROIs

## - Lists to store models
Models_ace=list()
Models_ae=list()
Models_ce=list()
Models_e=list()

## - Make models
for(i in 1:length(fa_fibers)){
  # -  Current fiber tract
  var = fa_fibers[i]
  # - Make competing models
  mod_ace=umxACE(selDVs = var, dzData = DZ_atlas, mzData = MZ_atlas, sep = "_T", autoRun = F)
  mod_ae=umxModify(mod_ace, update = "c_r1c1", name = "AE", autoRun = F)
  mod_ce=umxModify(mod_ace, update = "a_r1c1", name = "CE", autoRun = F)
  mod_e=umxModify(mod_ae, update = "a_r1c1", name = "E",    autoRun = F)
  # - Pack models into list
  Models_ace=Models_ace %>% list.append(mod_ace)
  Models_ae=Models_ae %>% list.append(mod_ae)
  Models_ce=Models_ce %>% list.append(mod_ce)
  Models_e=Models_e %>% list.append(mod_e)
}

results_ace = map(Models_ace,output,2) %>% reduce(rbind)
results_ae = map(Models_ae,output,2) %>% reduce(rbind)
results_ce = map(Models_ce,output,2) %>% reduce(rbind)
results_e = map(Models_e,output,2) %>% reduce(rbind)

# - Model parameter table
univariate_atlas=rbind(results_ace,results_ae,results_ce,results_e) %>%
  add_column(Fiber=rep(fa_fibers,4), .before="Model") %>%
  add_column(Net=rep(c(rep("Core",3),rep("Ext",4)),4)) %>% 
  add_column(Net_gen="Face_atlas") %>%
  mutate_at(.vars = "Fiber", ~ (factor(., levels = fa_fibers)))

### --- Complete result tables
fibers_plot=c("ATR","CGC","UNC","Fminor","Fmajor","CGH","ILF","SLF","CST","IFO",
              "FFA-OFA","FFA-V1V2","OFA-V1V2",
              "FFA-ATL","OFA-ATL","pSTS-ATL","ATL-V1V2")

# - Add proportion of explained variance 
univariate_results = rbind(univariate_func,univariate_atlas) %>% 
  mutate(ast=a^2/(a^2+c^2+e^2),
         cst=c^2/(a^2+c^2+e^2),
         est=e^2/(a^2+c^2+e^2)) %>%
  name_changer(fibers,fibers_plot)


## - Select best-fitting models based on AIC
best_gl=univariate_results %>% filter(Net_gen=="Global") %>%
  group_by(Fiber) %>%
  slice_min(n = 1, AIC) %>% 
  arrange(factor(Fiber, levels = c("ATR","CGC","UNC","Fminor","Fmajor","CGH","ILF","SLF","CST","IFO")))

best_fa=univariate_results %>% filter(Net_gen=="Face") %>%
  group_by(Fiber) %>%
  slice_min(n = 1, AIC) %>% 
  arrange(factor(Fiber, levels = c("FFA-OFA","FFA-V1V2","OFA-V1V2",
                                   "FFA-ATL","OFA-ATL","pSTS-ATL","ATL-V1V2")))

best_fa_atlas=univariate_results %>% filter(Net_gen=="Face_atlas") %>%
  group_by(Fiber) %>% 
  filter(AIC == min(AIC)) %>% 
  arrange(factor(Fiber, levels = c("FFA-OFA","FFA-V1V2","OFA-V1V2",
                                   "FFA-ATL","OFA-ATL","pSTS-ATL","ATL-V1V2")))


final_univariate_models = rbind(best_gl,best_fa,best_fa_atlas) %>%
  mutate_if(is.numeric,round,2) %>% select(Fiber, Model, a,aci,c,cci,e,eci,AIC,LogLike,BIC)


## - Figures
gl_plot = best_gl %>% 
  pivot_longer(cols=c("ast","cst","est"), names_to = "Source", values_to = "Variance") %>%
  mutate(Effect=ifelse(Source == "ast",
                       "Additive genetic (a)", ifelse(Source == "cst",
                                                      "Shared environment (c)","Unique influence (e)"))) %>%
  ggplot(aes(fill=Effect, y=Variance, x=Fiber)) +
  coord_flip() +
  geom_bar(stat="identity") + 
  theme(axis.title.x = element_text(color="black", size=11),axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_text(color="black", size=11),
        legend.position="bottom",
        panel.background = element_blank())

fa_plot = best_fa %>% 
  pivot_longer(cols=c("ast","cst","est"), names_to = "Source", values_to = "Variance") %>%
  mutate(Effect=ifelse(Source == "ast",
                       "Additive genetic (a)", ifelse(Source == "cst",
                                                      "Shared environment (c)","Unique influence (e)"))) %>%
  ggplot(aes(fill=Effect, y=Variance, x=Fiber)) +
  coord_flip() +
  geom_bar(stat="identity") + 
  theme(axis.title.x = element_text(color="black", size=11),axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_text(color="black", size=11),
        legend.position="bottom",
        panel.background = element_blank())

fa_plot_atlas = best_fa_atlas %>% 
  pivot_longer(cols=c("ast","cst","est"), names_to = "Source", values_to = "Variance") %>%
  mutate(Effect=ifelse(Source == "ast",
                       "Additive genetic (a)", ifelse(Source == "cst",
                                                      "Shared environment (c)","Unique influence (e)"))) %>%
  ggplot(aes(fill=Effect, y=Variance, x=Fiber)) +
  coord_flip() +
  geom_bar(stat="identity") + 
  theme(axis.title.x = element_text(color="black", size=11),axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_text(color="black", size=11),
        legend.position="bottom",
        panel.background = element_blank()) 

legend = get_legend(fa_plot_atlas)

# - Multiplota
prow <- plot_grid( gl_plot + theme(legend.position="none"),
                   fa_plot + theme(legend.position="none"),
                   fa_plot_atlas + theme(legend.position="none"),
                   align = 'vh',
                   labels = c("A", "B","C"),
                   hjust = -1,
                   nrow = 1)

prow2 <- plot_grid(prow + theme(legend.position="none"),
                   legend,
                   align = 'v',
                   hjust = -1,
                   nrow = 2,
                   rel_heights = c(1,0.1))


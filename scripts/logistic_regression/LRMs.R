library(PRROC)
library(caret)
library(plotROC)
library(pROC)

working_dir=<Path here>


ml_models=data.frame()
k=1
for (j in c("XL","EC","XLEC","iXLEC_docking","XLEC_Alphafold","Alphafold")){
  tryCatch({
    
    # specify the type of training method used $ 5 of folds 
    ctrlspecs= trainControl(method="cv", number=5,  classProbs = TRUE, savePredictions = TRUE)
    
    
    if(j=="XL"){
      all_models_stat=data.frame(read.delim(paste(working_dir,"models_xlcom_retrospect_clean_85.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
      
      
      for (i in 1:nrow(all_models_stat)){
        tryCatch({
          if(all_models_stat[i,"na_contacts"]>=10){all_models_stat[i,"interface"]=1}else{all_models_stat[i,"interface"]=0}
        },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
      }
      
      all_models_stat$interface[all_models_stat$interface==1]="one"
      all_models_stat$interface[all_models_stat$interface==0]="zero"
      all_models_stat[,"interface"]=as.factor(all_models_stat[,"interface"])
      class(all_models_stat[,"interface"])
      nrow( all_models_stat$interface[all_models_stat$interface==0])
      model=NULL
      model=train(interface ~n_xl_identified, data=all_models_stat, method="glm", family=binomial, trControl=ctrlspecs)
      saveRDS(model, paste(working_dir,"model_xl.rds",sep="/"))
    }
    
   
    
    if(j=="EC"){all_models_stat=data.frame(read.delim(paste(working_dir, "about_dataset","models_ecs_retrospect_clean_85.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
   
     for (i in 1:nrow(all_models_stat)){
      tryCatch({
        if(all_models_stat[i,"na_contacts"]>=10){all_models_stat[i,"interface"]=1}else{all_models_stat[i,"interface"]=0}
      },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
    }
    
    all_models_stat$interface[all_models_stat$interface==1]="one"
    all_models_stat$interface[all_models_stat$interface==0]="zero"
    all_models_stat[,"interface"]=as.factor(all_models_stat[,"interface"])
    class(all_models_stat[,"interface"])
    
    model=NULL
    model=train(interface ~N_eff +precent_ecs+ max_ecs+ kurtosis+ skewness, data=all_models_stat, method="glm", family=binomial, trControl=ctrlspecs)
    saveRDS(model, paste(working_dir,"model_ecs.rds",sep="/"))
    
    }
    
  
    
    if(j=="XLEC"){all_models_stat=data.frame(read.delim(paste(working_dir, "about_dataset","models_xlec2_retrospect_clean_85.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    for (i in 1:nrow(all_models_stat)){
      tryCatch({
        if(all_models_stat[i,"na_contacts"]>=10){all_models_stat[i,"interface"]=1}else{all_models_stat[i,"interface"]=0}
      },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
    }
    
    all_models_stat$interface[all_models_stat$interface==1]="one"
    all_models_stat$interface[all_models_stat$interface==0]="zero"
    all_models_stat[,"interface"]=as.factor(all_models_stat[,"interface"])
    class(all_models_stat[,"interface"])
    
    
    model=NULL
    model=train(interface ~ n_xl_identified+   N_eff +precent_ecs+ max_ecs+ kurtosis+ skewness,data=all_models_stat, method="glm", family=binomial, trControl=ctrlspecs)
    saveRDS(model, paste(working_dir,"model_xlec.rds",sep="/"))
    
    }
    

    if(j=="iXLEC_docking"){all_models_stat=data.frame(read.delim(paste(working_dir, "about_dataset","models_xlec2_retrospect_clean_85.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    for (i in 1:nrow(all_models_stat)){
      tryCatch({
        if(all_models_stat[i,"na_contacts"]>=10){all_models_stat[i,"interface"]=1}else{all_models_stat[i,"interface"]=0}
      },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
    }
    
    all_models_stat$interface[all_models_stat$interface==1]="one"
    all_models_stat$interface[all_models_stat$interface==0]="zero"
    all_models_stat[,"interface"]=as.factor(all_models_stat[,"interface"])
    class(all_models_stat[,"interface"])
    
    
    model=NULL 
    model=train(interface ~ Haddock_score+ + n_Haddock_cluster+ n_xl_identified+ n_xl_30Ang+ n_ecs_15Ang_0.999percentile+  N_eff +precent_ecs+ max_ecs+ kurtosis+ skewness,data=all_models_stat, method="glm", family=binomial, trControl=ctrlspecs)
    saveRDS(model, paste(working_dir,"model_xlec_dock.rds",sep="/"))

    
      }  
    
    
  
    if(j=="XLEC_Alphafold"){all_models_stat=data.frame(read.delim(paste(working_dir, "about_dataset","models_alphafold_retrospect_clean_85.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    for (i in 1:nrow(all_models_stat)){
      tryCatch({
        if(all_models_stat[i,"na_contacts"]>=10){all_models_stat[i,"interface"]=1}else{all_models_stat[i,"interface"]=0}
      },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
    }
    
    all_models_stat$interface[all_models_stat$interface==1]="one"
    all_models_stat$interface[all_models_stat$interface==0]="zero"
    all_models_stat[,"interface"]=as.factor(all_models_stat[,"interface"])
    class(all_models_stat[,"interface"])
    
    
    model=NULL 
    model =train(interface ~ AF_ptm +  n_xl_identified + n_xl_30Ang+ n_ecs_10Ang_0.999percentile + N_eff + precent_ecs + max_ecs + kurtosis +skewness ,data=all_models_stat, method="glm", family=binomial, trControl=ctrlspecs)
    saveRDS(model, paste(working_dir,"model_xlec_AF.rds",sep="/"))
    
    }
    
    
    if(j=="Alphafold"){all_models_stat=data.frame(read.delim(paste(working_dir, "about_dataset","models_alphafold_retrospect_clean_85.csv",sep="/"),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    for (i in 1:nrow(all_models_stat)){
      tryCatch({
        if(all_models_stat[i,"na_contacts"]>=10){all_models_stat[i,"interface"]=1}else{all_models_stat[i,"interface"]=0}
      },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
    }
    
    all_models_stat$interface[all_models_stat$interface==1]="one"
    all_models_stat$interface[all_models_stat$interface==0]="zero"
    all_models_stat[,"interface"]=as.factor(all_models_stat[,"interface"])
    class(all_models_stat[,"interface"])
    
    
    model=NULL 
    model =train(interface ~ AF_ptm ,data=all_models_stat, method="glm", family=binomial, trControl=ctrlspecs)
    saveRDS(model, paste(working_dir,"model_AF.rds",sep="/"))
    

    
    }
    
    
    
    names(all_models_stat)[40]="Displacement(Ang)"
    names(all_models_stat)[41]="Angle(deg)"
    
  
    
    
    #predict outcome using model from layer of cross-validation 
    # Out-of-sample prediction:
    predict=NULL
    predict=data.frame(model[["pred"]])
    
    
    
    # GLM variable inportance (predictor variables )
    varImp(model)
    
    ml_models[k,"ml_model"]=j
    # summarize results
    confmat=caret::confusionMatrix( predict$pred, predict$obs )
    confmat
    ml_models[k,"accuracy"]=round(confmat[["overall"]][["Accuracy"]],digits=3)
    ml_models[k,"balanced_accuracy"]=round(confmat[["byClass"]][["Balanced Accuracy"]],digits = 3)
    ml_models[k,"precision"]= round(confmat[["byClass"]][["Precision"]],digits=3)
    ml_models[k,"recall"]=round(confmat[["byClass"]][["Recall"]],digits=3)
    #F1=round((2 * precision * recall) / (precision + recall),digits=4)
    F1 = round( (2 * confmat[["byClass"]][["Precision"]] * confmat[["byClass"]][["Recall"]] ) / (confmat[["byClass"]][["Precision"]] + confmat[["byClass"]][["Recall"]] ) , digits=3)
    ml_models[k,"F1"]=round(F1, digits = 3)
    ml_models[k,"sensitivity"]=round(as.numeric(confmat[["byClass"]][["Sensitivity"]]), digits = 3)
    ml_models[k,"specificity"]=round(as.numeric(confmat[["byClass"]][["Specificity"]]), digits = 3)
    ml_models[k,"kappa"]=round(confmat[["overall"]][["Kappa"]], digits = 3)
    
    
    
    
    score1=NULL
    score1=predict[predict[,"obs"]=="one",]$one
    score0=NULL
    score0=predict[predict[,"obs"]=="zero",]$one
    ##################################
    #ROC plot
    ##################################
    roc=NULL
    #roc= roc.curve(predict$obs, predict$obs, curve = T)
    
    roc= roc.curve(score1, score0, curve = T)
    plot(roc, main = 'ROC curve',cex.axis=1, cex.main=1.5, cex.lab=1.4)
    ml_models[k,"roc_auc"]=round(roc$auc, digits = 3)
    roc_df=data.frame()
    roc_df=data.frame(roc[["curve"]])
    write.table(roc_df, file=paste(working_dir,"/roc_df_",j,".csv",sep=""),sep="\t",row.names=FALSE,quote=F)
    
    roc_plot=NULL
    roc_plot = ggplot(roc_df, aes(x=X1,y=X2,color=X3))
    roc_plot = roc_plot + geom_line(size=1.5) + theme(axis.text = element_text(colour = "blue")) 
    roc_plot = roc_plot + geom_abline(slope = 1, intercept = 0, size = 0.4)
    roc_plot = roc_plot + ggtitle(paste("ROC, AUC=",round(roc$auc,digits = 3),sep=" ")) 
    roc_plot = roc_plot +  xlab("Specificity")+ ylab("Senstivity")+ labs(color="Thresolds")
    roc_plot = roc_plot + theme_bw() +  theme(plot.background=element_blank(),panel.grid.minor= element_blank(),panel.grid.major = element_blank())
    roc_plot = roc_plot + theme(strip.text = element_text(size=16, color="black"),axis.text=element_text(size=10),axis.title=element_text(size=20),axis.text.x = element_text(size=14,color="black",  hjust = 1),axis.text.y = element_text(size=14,color="black"), legend.text = element_text(size = 14, color = "black"),  legend.title  = element_text(size = 16, color = "black"))
    roc_plot
    
    ggsave(paste( paste("roc",j, sep = "_"),".pdf",sep = ""), plot = last_plot(), path = working_dir , scale = 1, width = NA, height = NA, dpi=300)
    ggsave(paste( paste("roc",j, sep = "_"),".png",sep = ""), plot = last_plot(), path = working_dir , scale = 1, width = NA, height = NA, dpi=300)
    
    
    ##################################
    #ploting Precision-Recall Curve
    ##################################
    pr=NULL
    pr= pr.curve(score1, score0, curve = T)
    ml_models[k,"pr_auc"]=round(pr$auc.integral, digits = 3)
    pr_df=data.frame()
    pr_df=data.frame(pr[["curve"]])
    write.table(pr_df, file=paste(working_dir,"/pr_df_",j,".csv",sep=""),sep="\t",row.names=FALSE,quote=F)
    
    plot(pr, main="Out-of-sample PR curve")
    pr_plot =NULL
    pr_plot = ggplot(pr_df, aes(x = X1, y= X2, color=X3)) +  geom_path(size=1.5) 
    #pr_plot = pr_plot +  geom_abline(slope = 1, intercept = 0, size = 0.4)
    pr_plot = pr_plot +  xlab("Recall")+ ylab("Precision")+ labs(color="Thresolds")
    pr_plot = pr_plot + ggtitle(paste("PR curve, AUC=",round(pr$auc.integral,digits = 3),sep=" ")) 
    pr_plot = pr_plot + theme_bw() +  theme(plot.background=element_blank(),panel.grid.minor= element_blank(),panel.grid.major = element_blank())
    pr_plot = pr_plot + theme(strip.text = element_text(size=16, color="black"),axis.text=element_text(size=10),axis.title=element_text(size=20),axis.text.x = element_text(size=14,color="black",  hjust = 1),axis.text.y = element_text(size=14,color="black"), legend.text = element_text(size = 14, color = "black"),  legend.title  = element_text(size = 16, color = "black"))
    pr_plot
    
    ggsave(paste(paste("pr",j, sep = "_"),".pdf",sep = ""), plot = pr_plot, path = working_dir , scale = 1, width = NA, height = NA, dpi=300)
    ggsave(paste(paste("pr",j, sep = "_"),".png",sep = ""), plot = last_plot(), path = working_dir , scale = 1, width = NA, height = NA, dpi=300)
    
    assign( paste("roc",j, sep = "_"), roc)
    assign( paste("roc_df",j, sep = "_"), roc_df )
    assign( paste("pr_df",j, sep = "_"), pr_df )
    assign( paste("pr",j, sep = "_"), pr)
    assign( paste( "model",j, sep = "_"), model )
    k=k+1

   
    
    },error=function(e){cat("step=",j ,"ERROR :",conditionMessage(e), "\n")})
}

write.table(ml_models, file=paste(working_dir,"LR_models_statistics.txt",sep="/"),sep="\t",row.names=FALSE,quote=F)




library(rgl)
library(Rpdb)
library(bio3d)
library(reshape2)
library(tidyr)
library(reshape2)
library(pryr)
library(pracma)
library(seqinr)
library(fitdistrplus)
library(ggstance)
library(ggpubr)
library(tseries)

working_dir="<define working directory>"
protein_pair_list="<define path to list of protein paris>"
EC= "<define path to EVcomplex inter coupling file *CouplingScores_inter.csv>" 
EC_haddock="<define path to EC docking results directory>" 
EC_haddock_stat="<define path to statistics file of EC docking results>" 
EC_haddock_dist_rest= "<define path to EC restraints file with distances from docking>" 
python3="<define path to python3>/bin/python3"
haddock_ec_restraints="<define path to haddock_ec_restraints.py"

align_stat= "<define path to EVcomplex job summery file *job_statistics_summary.csv>"
Xl= "<define path to inter XL file>"
XL_haddock="<define path to XL docking results directory>" 
XL_haddock_stat="<define path to statistics file of XL docking results>" 
XL_haddock_dist_rest= "<define path to XL restraints file with distances from docking>" 
haddock_xl_restraints="<define path to haddock_xl_restraints.py"
stretcher= "<define path to stretcher excuatble>"
new_html= "<define path to haddock intialization file template>" 
cluster_dir="<define cluster directory where docking is taking place>"
submit_slurm= "<define path to cluster submission file template>" 

dir.create(paste(working_dir,"dataset", sep="/"))


list=data.frame(read.delim(paste(protein_pair_list),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
#In list file every ptotein pair should have 
#prefix: index number of every protein pairs 
#uid1: Uniprot ID for protein 1 as <ID_species>
#uid2: Uniprot ID for protein 2 as <ID_species>



#This will create a subdirectory for every protein pair under dataset directory named as prefix_protein1_protein2
#and downlaod the protein sequnce for uniprot webserver 
# station should be connected to the internet 

for (i in 1:nrow(list)){
  tryCatch({
    
    print(i)
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")


 },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}



########################################
#                                      #
#    calculate pairwise seq identity   #
#                                      #
########################################

#Run pair wise seqeunce similarity between pairs
for(i in 1:nrow(list)){
  tryCatch({
    print(i)

    #create a directory for the every protein
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
 
    # run structure for global sequence Identity 
    system(paste("cd", complex_dir, "&& ", stretcher,  "-asequence", paste(protein1_dir, paste(list[i,"uid1"],".fasta",sep = ""), sep='/') , "-bsequence", paste(protein2_dir, paste(list[i,"uid2"],".fasta",sep = ""), sep='/'), "-out" ,paste(list[i,"uid1"],"vs",list[i,"uid2"], "stretcher.out",sep="_"), sep=" "))
    
    temp_file=data.frame()
    temp_file=data.frame(read.delim(paste(complex_dir,paste(list[i,"uid1"],"vs",list[i,"uid2"], "stretcher.out",sep="_"), sep='/'), header=FALSE, stringsAsFactors = FALSE))
    temp_file=data.frame(grep("# Identity: ", temp_file[,"V1"],value = TRUE))
    list[i,"pairwise_identity"]=gsub(".*# Identity:    ", "", temp_file[,1])
 
  }, error=function(e){cat("ERROR: i=",i,conditionMessage(e), "\n")})
  
} 


list[,"pairwise_identity"]=gsub(".*\\(", "", list[,"pairwise_identity"])
list[,"pairwise_identity"]=gsub("\\).*", "", list[,"pairwise_identity"])
list[,"pairwise_identity"]=gsub("\\%.*", "", list[,"pairwise_identity"])



#################################
#                               #
#           Check EC            #
#                               #
#################################
list[,"ecs_complete"]= " "

for (i in 1:nrow(list)){
  tryCatch({
    
    print(i)
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
    
    ec_list=data.frame()
    ec_list=data.frame(read.delim(paste(EC),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    
    if(nrow(ec_list)!=0){list[i,"ecs_complete"]="Complete"}
    
    ec_list= ec_list[order( ec_list[,"cn"], decreasing = TRUE),]
    
    list[i,"kurtosis"]=kurtosis(ec_list[,"cn"])
    list[i,"skewness"]=skewness(ec_list[,"cn"])
    list[i,"squared_skewness"]=( list[i,"skewness"])^2
    list[i,"max_ecs"]=as.numeric(max(ec_list[,"cn"]))
    list[i,"n_ecs"]=as.numeric(nrow(ec_list))
    list[i,"precent_ecs"]=round(as.numeric(nrow(ec_list))/as.numeric(list[i,"seq1_len"]*list[i,"seq2_len"]),digits = 4)
    
    list[i,"Jarque-Bera"]=(list[i,"n_ecs"]/6)*(((list[i,"skewness"])^2) + (((list[i,"kurtosis"])^2)/4))
    
    align_list=data.frame()
    align_list=data.frame(read.delim(paste(align_stat),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    list[i,"N_eff"]=as.numeric(align_list[1,"N_eff"])
    
  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}


list=list[list[,"aln1"]!="" & list[,"aln2"]!="",]
list=list[list[,"ecs_complete"]=="Complete",]




#################################
#                               #
#           Check XL            #
#                               #
#################################


for (i in 1:nrow(list)){
  tryCatch({
    
    print(i)
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
    
    
    xl_list=data.frame()
    xl_list=data.frame(read.delim(paste(XL),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    xl_list=xl_list[!is.na(xl_list[,"uid1"]),]
    list[i,"n_xl_identified"]=as.numeric(nrow(xl_list))
    xl_list=xl_list[!is.na(xl_list[,"A_i"]) & !is.na(xl_list[,"A_j"]),]
    list[i,"n_xl_found"]=as.numeric(nrow(xl_list))
    
    
  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}

list=list[list[,"n_xl_found"]>=1,]




##########################################
#                                        #
#   Check XL & EC Dockings are complete  #
#                                        #
##########################################

list[,"dock_xl"]=" "
list[,"dock_ec"]=" "


for (i in 1:nrow(list)){
  tryCatch({
    
    print(i)
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
    
    haddock_check=data.frame()
    haddock_check=data.frame(read.delim(paste(XL_haddock,"run_xl/structures/it1/water/file.nam",sep="/"),header=FALSE, row.name=NULL, sep=" ",comment.char='#',stringsAsFactors = FALSE))
    if(nrow(haddock_check)!=0){list[i,"dock_xl"]="Complete"}

  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}



for (i in 1:nrow(list)){
  tryCatch({
    
    print(i)
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
    
    haddock_check=data.frame()
    haddock_check=data.frame(read.delim(paste(EC_haddock,"run_ec/structures/it1/water/file.nam",sep="/"),header=FALSE, row.name=NULL, sep=" ",comment.char='#',stringsAsFactors = FALSE))
    if(nrow(haddock_check)!=0){list[i,"dock_ec"]="Complete"}
    
  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}


#Remove cases when either of the docking protocols is not complete

list=list[ list[,"dock_xl"]=="Complete" &  list[,"dock_ec"]=="Complete",]
list=list[ !is.na(list[,"prefix"]),]




##########################################
#                                        #
#    Create XLEC hybrid restaints file   #
#                                        #
##########################################


#1. EC restraints 

for (i in 1:nrow(list)){
  tryCatch({
    print(i)
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")   
    
    dock_dir=" "
    dock_dir=paste(complex_dir,"haddock_xlec",sep="/")
    dir.create(dock_dir)
 
    #Upload EC docking statistics file 
    models_data=data.frame()
    models_data=data.frame(read.delim(paste(EC_haddock_stat),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    
     #filter to the model that achieved the highest number of EC restarints within less than 15 Angestrom 
     models_data=  models_data[ models_data[,"dock_protocol"] %in% c("ecs"),]
     models_data= models_data[order(as.numeric(as.character( models_data[,"n_ecs_15Ang"])), decreasing = TRUE),]
     models_data[,"test"]=paste( models_data[,"prefix"], models_data[,"dock_protocol"],sep=" ")
     models_data= models_data[!duplicated( models_data[,"test"]),]
     models_data= models_data[ models_data[,"prefix"]==list[i,"prefix"],]
    
    # upload the restraints distance file and select the restraints
    ec_list=data.frame()
    ec_list=data.frame(read.delim(paste(EC_haddock_dist_rest),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    
 
    if(models_data[1,"model_rank"] %in% c("1","2","3","4","5","6","7","8","9","10")){
      colname=paste(models_data[1,"model"],".pdb",sep="")
    }else{
      colname=paste("haddock_ecs_clust",sub(".*clust_","",models_data[1,"model_rank"]),"_top1_",sub(".*haddock_ecs_","",models_data[1,"model_no."]),".pdb",sep="")
    }
    
    ec_list= ec_list[ec_list[, colname]<=15,]
    ec_list= ec_list[!is.na(ec_list[, colname]),]
    write.csv(ec_list, file =paste(dock_dir, paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),"ecslist_distmodels_15A.csv",sep="_"),sep="/"),row.names=FALSE,  quote = TRUE)
   
    # generate restraints file
    system(paste("cd ",dock_dir," && python3 ", paste(haddock_ec_restraints) ,paste(dock_dir, paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),"ecslist_distmodels_15A.csv",sep="_"),sep="/"),paste(dock_dir,"unambig_ecs.tbl",sep="/"),nrow(ec_list),sep=" ")) 

  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}




#2. XL restraints 

for (i in 1:nrow(list)){
  tryCatch({
    print(i)
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
   
 
    dock_dir=" "
    dock_dir=paste(complex_dir,"haddock_xlec",sep="/")
    dir.create(dock_dir)
  
    #Upload XL docking statistics file 
    models_data=data.frame()
    models_data=data.frame(read.delim(paste(XL_haddock_stat),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
  
    #filter to the model that achieved the highest number of XL restarints within less than 35 Angestrom 
     models_data=  models_data[ models_data[,"dock_protocol"] %in% c("xl"),]
     models_data=  models_data[order(as.numeric(as.character(   models_data[,"n_xl_35Ang"])), decreasing = TRUE),]
     models_data[,"test"]=paste( models_data[,"prefix"], models_data[,"dock_protocol"],sep=" ")
     models_data=  models_data[!duplicated( models_data[,"test"]),]
     models_data=  models_data[ models_data[,"prefix"]==list[i,"prefix"],]
 

    # upload the restraints distance file and select the restraints
     xl_list=data.frame()
     xl_list=data.frame(read.delim(paste(XL_haddock_dist_rest),header=TRUE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    
    if(models_data[1,"model_rank"] %in% c("1","2","3","4","5","6","7","8","9","10")){
      colname=paste(models_data[1,"model"],".pdb",sep="")
    }else{
      colname=paste("haddock_xlcom_clust",sub(".*clust_","",models_data[1,"model_rank"]),"_top1_",sub(".*haddock_xlcom_","",models_data[1,"model_no."]),".pdb",sep="")
    }
    
    xl_list= xl_list[xl_list[, colname]<=35,]
    xl_list= xl_list[!is.na(xl_list[, colname]),]
    
    xl_list=xl_list[!xl_list[,"A_i"]==" ",]
    xl_list=xl_list[!xl_list[,"A_j"]==" ",]
    xl_list=xl_list[!is.na(xl_list[,"uid1"]),]
    write.csv(xl_list, file =paste(dock_dir , paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),"xlinks_distmodels_35A.csv",sep="_"),sep="/"),row.names=FALSE,  quote = TRUE)
  
   # generate restraints file
    system(paste("cd ",dock_dir," && python3 ", paste(haddock_xl_restraints),paste(dock_dir , paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),"xlinks_distmodels_35A.csv",sep="_"),sep="/"),paste(dock_dir,"unambig_xl.tbl",sep="/"),as.numeric(nrow(xl_list)),sep=" ")) 


  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}



#3. combining XL and EC restraints files in one file

for (i in 1:nrow(list)){
  tryCatch({
    
    print(i)
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
     
    dock_dir=" "
    dock_dir=paste(complex_dir,"haddock_xlec",sep="/")
    dir.create(dock_dir)
    
    xl_unambig=data.frame()
    xl_unambig=data.frame(read.delim(paste(dock_dir,"unambig_xl.tbl",sep="/"),header=FALSE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    
    ecs_unambig=data.frame()
    ecs_unambig=data.frame(read.delim(paste(dock_dir,"unambig_ecs.tbl",sep="/"),header=FALSE, row.name=NULL, sep=",",comment.char='#',stringsAsFactors = FALSE))
    
    unambig=data.frame()   
    unambig=rbind(ecs_unambig,xl_unambig)
    names(unambig)[1]=" "
    write.csv(unambig, file=paste(dock_dir,"unambig.tbl",sep="/"), row.names=FALSE,  quote = FALSE)

    
  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}




###############################################
#                                             #
#    Create XLEC docking intialization file   #
#                                             #
###############################################

for (i in 1:nrow(list)){
  
  tryCatch({
    
    print(i)
    
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
    
    
    dock_dir=" "
    dock_dir=paste(complex_dir,"haddock_xlec",sep="/")

    new.html=data.frame()
    new.html=data.frame(read.delim(paste(new_html), stringsAsFactors = FALSE))
    
    if(list[i,"protein1_struc.model"]=="structure"){
      new.html[10,]=paste("PDB_FILE1=",  paste(paste(sub("\\_.*", "",list[i,"uid1"]),list[i,"protein1_pdb"],list[i,"protein1_chain"],sep='_'),"pdb", sep = '.'),"<BR>",sep="")
    }else{
      new.html[10,]=paste("PDB_FILE1=",  paste(paste(sub("\\_.*", "",list[i,"uid1"]),"hm",list[i,"protein1_model_pdb"],list[i,"protein1_model_chain"],sep='_'),"pdb", sep = '.'),"<BR>",sep="")
    }
    
    if(list[i,"protein2_struc.model"]=="structure"){
      new.html[11,]=paste("PDB_FILE2=",   paste(paste(sub("\\_.*", "",list[i,"uid2"]),list[i,"protein2_pdb"],list[i,"protein2_chain"],sep='_'),"pdb", sep = '.'), "<BR>",sep="")
    }else{
      new.html[11,]=paste("PDB_FILE2=",  paste(paste(sub("\\_.*", "",list[i,"uid2"]),"hm",list[i,"protein2_model_pdb"],list[i,"protein2_model_chain"],sep='_'),"pdb", sep = '.'),"<BR>",sep="")
    }
    
    new.html[12,]=paste("PROJECT_DIR=",cluster_dir,"<BR>",sep="")
    new.html[13,]=paste("RUN_NUMBER=_xlec","<BR>",sep="")
    new.html[14,]=paste("UNAMBIG_TBL=",paste(cluster_dir,"unambig.tbl",sep="/"),"<BR>",sep="")
    write.csv(new.html, file=paste(dock_dir,"new.html",sep='/'),row.names=FALSE,  quote = FALSE)
    
  
  },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}




###############################################
#                                             #
#       Create XLEC Slurm submission file     #
#                                             #
###############################################


for (i in 1:nrow(list)){
  #for (i in 1:10){
  tryCatch({
    #if(list[,"protein1_pdbtaxonomy"]=="Homo sapiens" &&list[,"protein2_pdbtaxonomy"]=="Homo sapiens"){
    
    complex_dir=" "
    complex_dir=paste(working_dir,"dataset/complexes", paste(list[i,"prefix"],sub("\\_.*", "",list[i,"uid1"]),sub("\\_.*", "",list[i,"uid2"]),sep='_'),sep="/")
   
    dock_dir=" "
    dock_dir=paste(complex_dir,"haddock_xlec",sep="/")
    submit_file=data.frame()
    submit_file=data.frame(read.delim(submit_slurm, stringsAsFactors = FALSE))
    submit_file[12,]=" "
    submit_file[13,]=" "
    submit_file[13,]=paste("cd ", cluster_dir ,"&& python /n/groups/marks/software/haddock2.2/Haddock/RunHaddock.py", sep=" ")
    submit_file[14,]=" "
    submit_file[15,]=paste("cd ", paste(cluster_dir,"run_xlec",sep="/") ,"&& python /n/groups/marks/software/haddock2.2/Haddock/RunHaddock.py > run_xlec.log", sep=" ")                   
    submit_file[16,]=paste("cd ", paste(cluster_dir,"run_xlec",sep="/") ,"&& rm -r structures/it0", sep=" ")                   
    submit_file[17,]=paste("cd ", paste(cluster_dir,"run_xlec/structures/it1/water",sep="/") ,"&& $HADDOCKTOOLS/ana_structures.csh", sep=" ")
    submit_file[18,]=paste("cd ", paste(cluster_dir,"run_xlec/structures/it1/water/analysis/",sep="/") ,"&& $HADDOCKTOOLS/cluster_struc -f  haddock_xlec_rmsd.disp 5 4 >cluster.out", sep=" ")
    submit_file[19,]=paste("cd ", paste(cluster_dir,"run_xlec/structures/it1/water/",sep="/") ,"&& $HADDOCKTOOLS/ana_clusters.csh -best 4 analysis/cluster.out", sep=" ")
    names(submit_file)[1]="#!/bin/bash" 
    write.csv(submit_file, file=paste(dock_dir,"submit_haddock.sh",sep='/'),row.names=FALSE,  quote = FALSE)
    
    },error=function(e){cat("step=", i ,"ERROR :",conditionMessage(e), "\n")})
}







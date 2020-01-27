 ```
>  library(hyperoverlap)
>  library(tictoc)
>  library(foreach)
>  library(doParallel)
>  library(doParallel)
>  library(hypervolume)
>  
>  df = read.csv("conifer_masterdata_fixed.csv")
>  
>  ##############################################
>  # HYPEROVERLAP DEFAULT 3D
>  ##############################################
>  no_cores <- 10
>  cl<-makeCluster(no_cores)
>  registerDoParallel(cl)
>  
>  foreach(i=1:10, .packages = c("hyperoverlap", "tictoc")) %dopar% {
>    
>    tic()
>    x = hyperoverlap_set(df[,6:8], df$Genus, cost=1000, kernel.degree = 5, stoppage.threshold = 0.4, write.to.file=TRUE,path=paste0("hyperoverlap_", Sys.time(),"_",i, "/"))
>    write.csv(x, paste0("results_",i,".csv"))
>    toc()
>  }
>  
>  stopCluster(cl)
>  
>  ##############################################
>  # LINEAR ONLY HYPEROVERLAP 3D
>  ##############################################
>  no_cores <- 10
>  cl<-makeCluster(no_cores)
>  registerDoParallel(cl)
>  
>  foreach(i=1:10, .packages = c("hyperoverlap", "tictoc")) %dopar% {
>  
>  tic()
>  x = hyperoverlap_set(df[,6:8], df$Genus, cost=1000, kernel="linear", stoppage.threshold = 0.4, write.to.file=FALSE)
>  write.csv(x, paste0("results_",i,".csv"))
>  toc()
>  }
>  
>  stopCluster(cl)
>  #######################################################
>  # HYPERVOLUME 3D:
>  #######################################################
>  
>  taxa <- as.character(unique(df$Genus))
>  which(taxa=="Wollemia")
>  
>  taxa=taxa[-60] #remove Wollemia
>  taxa=taxa[-28] #remove
>  dir.create("hv")
>  
>  ######################
>  # DEFAULT PARAMETERS 3D
>  ######################
>  
>  no_cores <- 10
>  cl<-makeCluster(no_cores)
>  registerDoParallel(cl)
>  
>  foreach(k=1:10, .packages = c("hypervolume", "tictoc")) %dopar% {
>    
>    tic()
>    #build hypervolumes
>    for (i in 1:length(taxa)){
>      name = taxa[i]
>      dat = df[which(df$Genus==taxa[i]),6:10]
>      hv = hypervolume_svm(dat,name=name, verbose=F)
>      saveRDS(hv, file=paste0("./hv/",name,k,".rds"))
>    }
>    
>    results = list(NULL,NULL,NULL,NULL)
>    names(results) <-  c("entity1", "entity2", "result", "jaccard")
>    
>    
>    for (i in 1:(length(taxa))){
>      hv1 = readRDS(paste0("./hv/",taxa[i],k,".rds"))
>      
>      
>      for (j in 1:length(taxa)){
>        
>        if(j==i){
>          results[[1]] <- c(results[[1]], taxa[i])
>          results[[2]] <- c(results[[2]], taxa[i])
>          results[[3]] <- c(results[[3]], " ")
>          results[[4]] <- c(results[[4]], NA)
>        }
>          
>          
>          if (j>  i){
>            
>            hv2 = readRDS(paste0("./hv/",taxa[j],k,".rds"))
>            hvset=hypervolume_set(hv1,hv2, check.memory=F, verbose=F)
>            
>            results[[1]] <- c(results[[1]], hv1@Name)
>            results[[2]] <- c(results[[2]], hv2@Name)
>            
>            if(hypervolume_overlap_statistics(hvset)[1]==0){
>              results[[3]] <- c(results[[3]], "non-overlap")
>              results[[4]] <- c(results[[4]], hypervolume_overlap_statistics(hvset)[1])
>            } else {
>              results[[3]] <- c(results[[3]], "overlap")
>              results[[4]] <- c(results[[4]], hypervolume_overlap_statistics(hvset)[1])
>            }
>            
>          }
>        }
>      }
>      
>      results <-  data.frame(results)
>      results$entity1 <- factor(results$entity1,levels=taxa)
>      results$entity2 <-  factor(results$entity2,levels=taxa)
>      results <- results[order(results$entity1,results$entity2),]
>      write.csv(results, paste0("results_hv_",k,".csv"))
>      
>    toc()
>  } 
>  stopCluster(cl)
>  
>  ######################
>  # SAMPLES PER POINT INCREASED BY FACTOR OF 100 3D
>  ######################
>  
>  
>  no_cores <- 10
>  cl<-makeCluster(no_cores)
>  registerDoParallel(cl)
>  
>  foreach(k=1:10, .packages = c("hypervolume", "tictoc")) %dopar% {
>    
>    tic()
>    #build hypervolumes
>    for (i in 1:length(taxa)){
>      name = taxa[i]
>      dat = df[which(df$Genus==taxa[i]),6:10]
>      hv = hypervolume_svm(dat,name=name, verbose=F, samples.per.point=ceiling((10^(5 + sqrt(ncol(dat))))/nrow(dat)))
>      saveRDS(hv, file=paste0("./hv/",name,k,".rds"))
>    }
>    
>    results = list(NULL,NULL,NULL,NULL)
>    names(results) <-  c("entity1", "entity2", "result", "jaccard")
>    
>    #calculate overlap
>    for (i in 1:(length(taxa))){
>      hv1 = readRDS(paste0("./hv/",taxa[i],k,".rds"))
>      
>      
>      for (j in 1:length(taxa)){
>        
>        if(j==i){
>          results[[1]] <- c(results[[1]], taxa[i])
>          results[[2]] <- c(results[[2]], taxa[i])
>          results[[3]] <- c(results[[3]], " ")
>          results[[4]] <- c(results[[4]], NA)
>        }
>        
>        
>        if (j>  i){
>          
>          hv2 = readRDS(paste0("./hv/",taxa[j],k,".rds"))
>          hvset=hypervolume_set(hv1,hv2, check.memory=F, verbose=F)
>          
>          results[[1]] <- c(results[[1]], hv1@Name)
>          results[[2]] <- c(results[[2]], hv2@Name)
>          
>          if(hypervolume_overlap_statistics(hvset)[1]==0){
>            results[[3]] <- c(results[[3]], "non-overlap")
>            results[[4]] <- c(results[[4]], hypervolume_overlap_statistics(hvset)[1])
>          } else {
>            results[[3]] <- c(results[[3]], "overlap")
>            results[[4]] <- c(results[[4]], hypervolume_overlap_statistics(hvset)[1])
>          }
>          
>        }
>      }
>    }
>    
>    results <-  data.frame(results)
>    results$entity1 <- factor(results$entity1,levels=taxa)
>    results$entity2 <-  factor(results$entity2,levels=taxa)
>    results <- results[order(results$entity1,results$entity2),]
>    write.csv(results, paste0("results_hv5_",k,".csv"))
>    
>    toc()
>  } 
>  stopCluster(cl)
>  
>  
>  # FIVE DIMENSIONS
>  
>  
>  ##############################################
>  # HYPEROVERLAP DEFAULT 5D
>  ##############################################
>  no_cores <- 10
>  cl<-makeCluster(no_cores)
>  registerDoParallel(cl)
>  
>  foreach(i=1:10, .packages = c("hyperoverlap", "tictoc")) %dopar% {
>    
>    tic()
>    x = hyperoverlap_set(df[,6:10], df$Genus, cost=1000, kernel.degree = 5, stoppage.threshold = 0.4, write.to.file=TRUE,path=paste0("hyperoverlap_", Sys.time(),"_",i, "/"))
>    write.csv(x, paste0("results_",i,".csv"))
>    toc()
>  }
>  
>  stopCluster(cl)
>  
>  ##############################################
>  # LINEAR ONLY HYPEROVERLAP 5D
>  ##############################################
>  no_cores <- 10
>  cl<-makeCluster(no_cores)
>  registerDoParallel(cl)
>  
>  foreach(i=1:10, .packages = c("hyperoverlap", "tictoc")) %dopar% {
>    
>    tic()
>    x = hyperoverlap_set(df[,6:10], df$Genus, cost=1000, kernel="linear", stoppage.threshold = 0.4, write.to.file=FALSE)
>    write.csv(x, paste0("results_",i,".csv"))
>    toc()
>  }
>  
>  stopCluster(cl)
>  #######################################################
>  # HYPERVOLUME 5D
>  #######################################################
>  
>  taxa=taxa[-28] #remove Metasequoia (too few points)
>  
>  ######################
>  # DEFAULT PARAMETERS 5D
>  ######################
>  
>  no_cores <- 10
>  cl<-makeCluster(no_cores)
>  registerDoParallel(cl)
>  
>  foreach(k=1:10, .packages = c("hypervolume", "tictoc")) %dopar% {
>    
>    tic()
>    #build hypervolumes
>    for (i in 1:length(taxa)){
>      name = taxa[i]
>      dat = df[which(df$Genus==taxa[i]),6:10]
>      hv = hypervolume_svm(dat,name=name, verbose=F)
>      saveRDS(hv, file=paste0("./hv/",name,k,".rds"))
>    }
>    
>    results = list(NULL,NULL,NULL,NULL)
>    names(results) <-  c("entity1", "entity2", "result", "jaccard")
>    
>    
>    for (i in 1:(length(taxa))){
>      hv1 = readRDS(paste0("./hv/",taxa[i],k,".rds"))
>      
>      
>      for (j in 1:length(taxa)){
>        
>        if(j==i){
>          results[[1]] <- c(results[[1]], taxa[i])
>          results[[2]] <- c(results[[2]], taxa[i])
>          results[[3]] <- c(results[[3]], " ")
>          results[[4]] <- c(results[[4]], NA)
>        }
>        
>        
>        if (j>  i){
>          
>          hv2 = readRDS(paste0("./hv/",taxa[j],k,".rds"))
>          hvset=hypervolume_set(hv1,hv2, check.memory=F, verbose=F)
>          
>          results[[1]] <- c(results[[1]], hv1@Name)
>          results[[2]] <- c(results[[2]], hv2@Name)
>          
>          if(hypervolume_overlap_statistics(hvset)[1]==0){
>            results[[3]] <- c(results[[3]], "non-overlap")
>            results[[4]] <- c(results[[4]], hypervolume_overlap_statistics(hvset)[1])
>          } else {
>            results[[3]] <- c(results[[3]], "overlap")
>            results[[4]] <- c(results[[4]], hypervolume_overlap_statistics(hvset)[1])
>          }
>          
>        }
>      }
>    }
>    
>    results <-  data.frame(results)
>    results$entity1 <- factor(results$entity1,levels=taxa)
>    results$entity2 <-  factor(results$entity2,levels=taxa)
>    results <- results[order(results$entity1,results$entity2),]
>    write.csv(results, paste0("results_hv_",k,".csv"))
>    
>    toc()
>  } 
>  stopCluster(cl)
>  
>  
>  ###########################################

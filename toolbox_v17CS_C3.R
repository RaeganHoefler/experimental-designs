superRunner<-function(task_scheduler_path=NULL,task_index=NULL,task_geno_dictionary_path=NULL,phenodata_path=NULL,yd_monitor_path=NULL,n_geno_iter=NULL,h2=NULL,step_ar1=0.1,step_ar1ar1=0.1,fieldbook.return=FALSE,model.return=FALSE, ge="ds",cores=0, Amat_path = 0){
  tmp<-numeric()
  fbb<-list()
  modell<-list()
  for (ki in 1:n_geno_iter){
    
    task_scheduler<-read.csv(task_scheduler_path,stringsAsFactors = FALSE)
    task_geno_dictionary<-readRDS(task_geno_dictionary_path)
    pheno0<-readRDS(phenodata_path)[[task_scheduler$pheno_set[task_index]]]
    pheno<-pheno0[task_geno_dictionary[[task_scheduler$map[task_index]]][[task_scheduler$size[task_index]]][task_scheduler$stpoint[task_index],],]
    
    ##geno<-readRDS('core_base/AA.rds')
    yd_monitor<-readRDS(yd_monitor_path)
    
    fieldbook<-list(fb=layoutGenerator(design = task_scheduler$design[task_index],
                                       n_row = task_scheduler$n_row[task_index],
                                       n_col = task_scheduler$n_col[task_index],
                                       n_env = task_scheduler$n_env[task_index],  ## 5 environments here
                                       trt = task_scheduler$trt[task_index],
                                       rep = task_scheduler$rep[task_index],  ## for the alpha design, this is equal to the number of complete blocks
                                       inc_blocks = task_scheduler$inc_blocks[task_index],
                                       checks = task_scheduler$checks[task_index],
                                       plot = FALSE),
                    id=task_scheduler$id[task_index],   
                    trial_size=task_scheduler$size[task_index],
                    map=task_scheduler$map[task_index],
                    pheno_set=task_scheduler$pheno_set[task_index],
                    x=task_scheduler$col[task_index],
                    y=task_scheduler$row[task_index])
    
 
    
    
    ## here, we update the fieldbook$fb with real phenotypic data and yield monitor information
    fieldbook$fb<-layoutFiller(fieldbook = fieldbook$fb,
                               pheno = pheno, ## this dataframe must contain the same number of columns as specificed in the fieldbook[[i]]$fb$env
                               yd_monitor = yd_monitor[[fieldbook$map]],  ## the yield monitor to use depends on the size of the design 
                               x = fieldbook$x, ## this are the starting points. in the future, we are going to randomize this number
                               y = fieldbook$y, ## this are the starting points. in the future, we are going to randomize this number
                               h2 = h2) ## some extra noise
    
    
    cat(paste0('\n-----------------------------------------------------------\n',
               'TASK SCHEDULER: model with id=',fieldbook$id,'\n',
               'yield map ',fieldbook$map,', size ',fieldbook$trial_size,
               ', starting points row=',fieldbook$y,', col=',fieldbook$x,'\n*ITERATION: ',ki,'/',n_geno_iter,
               '\nNOTE 1: printing fieldbook is ',ifelse(fieldbook.return,'ACTIVE','FALSE'),', printing models is ',ifelse(model.return,'ACTIVE','FALSE'),
               '\nNOTE 2: parallelization is ',ifelse(cores>0,paste0('ACTIVE and output printing is limited [used cores=',cores,']'),'INACTIVE and output printing is normald'),
               '\n-----------------------------------------------------------\n'))
    
    
    
    ## USING THIS K MATRIX FOR NOW. IN THE FUTURE, WE WILL PROVIDE MARKER DATA AS AN ARGUMENT IN THE SUPERRUNNER FUNCTION
    if (Amat_path==0){
      K <- diag(length(unique(fieldbook$fb$trt)))
      colnames(K) <- rownames(K) <- as.character(unique(fieldbook$fb$trt))
    }else{
      K <- as.matrix(read.table(Amat_path))
      colnames(K)<-rownames(K)
      K <- K[as.character(unique(fieldbook$fb$trt_names)),as.character(unique(fieldbook$fb$trt_names))]
      colnames(K) <- rownames(K) <- as.character(unique(fieldbook$fb$trt))
    }
    
    E <- diag(length(unique(fieldbook$fb$env)))
    colnames(E) <- rownames(E) <- as.character(unique(fieldbook$fb$env))
    
    ## now, invoking the coreRunner function
    if (model.return){
      
      tmp0<-coreRunner(design = task_scheduler$design[task_index],fieldbookInput = fieldbook$fb,E2=E,K2 = K,debug=model.return,step_ar1a=step_ar1,step_ar1ar1a = step_ar1ar1, ge=ge,ncore=cores)
      tmp<-rbind(tmp,cbind(fieldbook$id,fieldbook$x,fieldbook$y,fieldbook$map,fieldbook$trial_size,task_scheduler$stpoint[task_index],
                           tmp0[[1]]))
      modell[[ki]]<-tmp0[[2]]
    }else{
      tmp<-rbind(tmp,cbind(fieldbook$id,fieldbook$x,fieldbook$y,fieldbook$map,fieldbook$trial_size,task_scheduler$stpoint[task_index],
                           coreRunner(design = task_scheduler$design[task_index],fieldbookInput = fieldbook$fb,E2=E,K2 = K,debug=model.return,step_ar1a=step_ar1,step_ar1ar1a = step_ar1ar1, ge=ge,ncore=cores))) 
    }
    
    cat(paste0('\nTASK SCHEDULER: ITERATION ',ki,'/',n_geno_iter,' done!'))
    fbb[[ki]]<-fieldbook
  }
  aindex<-ifelse(Amat_path==0,'no','yes')
  tmp2<-data.frame(id=tmp[,1],row=tmp[,3],col=tmp[,2],map=tmp[,4],trial_size=tmp[,5],st_point=tmp[,6],pheno_set=task_scheduler$pheno_set[task_index],markerdata=aindex,her=h2,tmp[,-c(1:6)])
  if (model.return & fieldbook.return){
    return(list(tmp2,fbb,modell))
  }
  if (fieldbook.return){
    return(list(tmp2,fbb))
  }
  if (model.return){
    return(list(tmp2,modell))
  }
  if (model.return==FALSE & fieldbook.return==FALSE){
    return(tmp2)
  }
  
}



superRunnerSpatial<-function(spatial_fieldbook_dictionary_path=NULL,
                             task_scheduler_spatial_path=NULL,
                             task_index=NULL,
                             phenodata_path=NULL,
                             task_geno_dictionary_path=NULL,
                             yd_monitor_path=NULL,
                             n_geno_iter=NULL,
                             h2=NULL,
                             step_ar1=0.1,
                             step_ar1ar1=0.1,
                             fieldbook.return=FALSE,
                             model.return=FALSE, 
                             ge="ds",
                             cores=0, 
                             Amat_path = 0){
  tmp<-numeric()
  fbb<-list()
  modell<-list()
  for (ki in 1:n_geno_iter){
    spatial_fb_dir<-readRDS(spatial_fieldbook_dictionary_path)
    task_scheduler_spatial<-read.csv(task_scheduler_spatial_path)
    
    ## grabbing variables
    yds<-task_scheduler_spatial$yds[task_index]
    size<-task_scheduler_spatial$size[task_index]
    starting_point<-task_scheduler_spatial$st_point[task_index]
    design<-task_scheduler_spatial$spatial_model[task_index]
    pheno_set<-task_scheduler_spatial$pheno[task_index]
    
    task_geno_dictionary<-readRDS(task_geno_dictionary_path)
    pheno0<-readRDS(phenodata_path)[[pheno_set]]
    pheno<-pheno0[task_geno_dictionary[[yds]][[size]][starting_point,],]
    
    ##geno<-readRDS('core_base/AA.rds')
    yd_monitor<-readRDS(yd_monitor_path)
    
    #spatial_fb_dir<-tmp
    
    fieldbook<-list(fb=spatial_fb_dir[[yds]][[size]][[starting_point]][[design]]$fieldbook,
                    id=spatial_fb_dir[[yds]][[size]][[starting_point]][[design]]$id,   
                    trial_size=size,
                    map=yds,
                    pheno_set=starting_point,
                    x=spatial_fb_dir[[yds]][[size]][[starting_point]][[design]]$x,
                    y=spatial_fb_dir[[yds]][[size]][[starting_point]][[design]]$y)
    
    
    
    
    ## here, we update the fieldbook$fb with real phenotypic data and yield monitor information
    
    fieldbook$fb<-layoutFiller(fieldbook = fieldbook$fb,
                               pheno = pheno, ## this dataframe must contain the same number of columns as specificed in the fieldbook[[i]]$fb$env
                               yd_monitor = yd_monitor[[fieldbook$map]],  ## the yield monitor to use depends on the size of the design 
                               x = 1, ## this are the starting points. in the future, we are going to randomize this number
                               y = 1, ## this are the starting points. in the future, we are going to randomize this number
                               h2 = h2) ## some extra noise
    
    
    
    cat(paste0('\n-----------------------------------------------------------\n',
               'TASK SCHEDULER: model with id=',fieldbook$id,'\n',
               'yield map ',fieldbook$map,', size ',fieldbook$trial_size,
               ', starting points row=',fieldbook$y,', col=',fieldbook$x,'\n*ITERATION: ',ki,'/',n_geno_iter,
               '\nNOTE 1: printing fieldbook is ',ifelse(fieldbook.return,'ACTIVE','FALSE'),', printing models is ',ifelse(model.return,'ACTIVE','FALSE'),
               '\nNOTE 2: parallelization is ',ifelse(cores>0,paste0('ACTIVE and output printing is limited [used cores=',cores,']'),'INACTIVE and output printing is normald'),
               '\n-----------------------------------------------------------\n'))
    
    
    
    ## USING THIS K MATRIX FOR NOW. IN THE FUTURE, WE WILL PROVIDE MARKER DATA AS AN ARGUMENT IN THE SUPERRUNNER FUNCTION
    if (Amat_path==0){
      K <- diag(length(unique(fieldbook$fb$trt)))
      colnames(K) <- rownames(K) <- as.character(unique(fieldbook$fb$trt))
    }else{
      K <- as.matrix(read.table(Amat_path))
      colnames(K)<-rownames(K)
      K <- K[as.character(unique(fieldbook$fb$trt_names)),as.character(unique(fieldbook$fb$trt_names))]
      colnames(K) <- rownames(K) <- as.character(unique(fieldbook$fb$trt))
    }
    
    E <- diag(length(unique(fieldbook$fb$env)))
    colnames(E) <- rownames(E) <- as.character(unique(fieldbook$fb$env))
    
    ## now, invoking the coreRunner function
    if (model.return){
      
      tmp0<-coreRunner_sp(fieldbookInput = fieldbook$fb,E2=E,K2 = K,debug=model.return,step_ar1a=step_ar1,step_ar1ar1a = step_ar1ar1, ge=ge,ncore=cores)
      tmp<-rbind(tmp,cbind(fieldbook$id,fieldbook$x,fieldbook$y,fieldbook$map,fieldbook$trial_size,
                           tmp0[[1]]))
      modell[[ki]]<-tmp0[[2]]
    }else{
      tmp<-rbind(tmp,cbind(fieldbook$id,fieldbook$x,fieldbook$y,fieldbook$map,fieldbook$trial_size,
                           coreRunner_sp(fieldbookInput = fieldbook$fb,E2=E,K2 = K,debug=model.return,step_ar1a=step_ar1,step_ar1ar1a = step_ar1ar1, ge=ge,ncore=cores))) 
    }
    
    cat(paste0('\nTASK SCHEDULER: ITERATION ',ki,'/',n_geno_iter,' done!'))
    fbb[[ki]]<-fieldbook
  }
  
  aindex<-ifelse(Amat_path==0,'no','yes')
  tmp2<-data.frame(id=tmp[,1],row=tmp[,3],col=tmp[,2],map=tmp[,4],trial_size=tmp[,5],stpoint=starting_point,pheno_set=pheno_set,markerdata=aindex,her=h2,tmp[,-c(1:5)])
  if (model.return & fieldbook.return){
    return(list(tmp2,fbb,modell))
  }
  if (fieldbook.return){
    return(list(tmp2,fbb))
  }
  if (model.return){
    return(list(tmp2,modell))
  }
  if (model.return==FALSE & fieldbook.return==FALSE){
    return(tmp2)
  }
  
}


## task generator

taskGenerator<-function(scheduler_path,yd_dim,n_starting_points,total_geno){
  ## master scheduler
  master_scheduler<-read.csv(scheduler_path,header=TRUE,stringsAsFactors = FALSE)
  
  stp<-stp2<-list()
  for (i in 1:length(unique(master_scheduler$size))){
    xx<-master_scheduler[which(master_scheduler$size==i),]
    maxtrt<-max(xx$trt)
    tm<-which(is.na(xx$checks)==FALSE)
    maxtrt<-maxtrt+max(xx$checks[tm])
    maxrow<-max(xx$n_row)
    maxcol<-max(xx$n_col)
    tmp<-matrix(NA,n_starting_points,4)
    tmp2<-matrix(NA,n_starting_points,maxtrt)
    
    colrange<-1:(yd_dim[2]-maxcol+1)
    rowrange<-1:(yd_dim[1]-maxrow+1)
    
    
    for (k in 1:n_starting_points){
      y<-sample(rowrange,1)
      x<-sample(colrange,1)
      id_geno<-sample(1:total_geno,maxtrt)
      tmp[k,]<-c(k,y,x,k)   #map, startig point, x and y 
      tmp2[k,]<-id_geno
    }
    stp[[i]]<-tmp
    stp2[[i]]<-tmp2
  }
  
  tasks<-numeric()
  for (i in 1:nrow(master_scheduler)){
    tt<-do.call('rbind',replicate(n_starting_points,
                                  master_scheduler[i,],simplify = FALSE))
    task<-data.frame(stpoint=stp[[master_scheduler$size[i]]][,1],
                     row=stp[[master_scheduler$size[i]]][,2],
                     col=stp[[master_scheduler$size[i]]][,3],
                     tt)
    tasks<-rbind(tasks,task)
  } 
  return(list(tasks=tasks,geno_dictionary=stp2,unique_starting_points=stp))
}


## spatial fieldbook creator
#starting_points_dic<-starting_points_dictionary
#spatial_master_scheduler<-read.csv('~/Dropbox/inifap/proyectos/lucia_gxe2_2019/core_base/spatial_master_scheduler.csv',sep = ',')
#yd_monitor<-readRDS('~/Dropbox/inifap/proyectos/lucia_gxe2_2019/core_base/inuse/yd_monitor.rds')

spatial_fb_creator<-function(starting_points_dic=NULL,spatial_master_scheduler=NULL,yd_monitor=NULL,yds='all',size='all',sp='all',rc=NULL){
  
  ## do all for everything
  if (yds=='all' & size=='all' & sp[1]=='all'){
    st_dic<-starting_points_dic
  }
  ## all yds, a specific size, all starting points
  if (yds == 'all' & size != 'all' & sp[1] == 'all'){
    st_dic<-list()
    for (i in 1:length(starting_points_dic)){
      st_dic[[i]]<-list()
      st_dic[[i]][[1]]<-starting_points_dic[[i]][[size]]
    }
  }
  ## all yds, a specific size, a specific starting point
  if (yds == 'all' & size != 'all' & sp[1] != 'all'){
    st_dic<-list()
    for (i in 1:length(starting_points_dic)){
      st_dic[[i]]<-list()
      st_dic[[i]][[1]]<-matrix(starting_points_dic[[i]][[size]][sp,],ncol = 4)
    }
  }
  
  ## all yds, all sizes, a specific starting point
  if (yds == 'all' & size == 'all' & sp[1] != 'all'){
    st_dic<-list()
    for (i in 1:length(starting_points_dic)){
      st_dic[[i]]<-list()
      for (j in 1:length(starting_points_dic[[i]])){
        st_dic[[i]][[j]]<-matrix(starting_points_dic[[i]][[j]][sp,],ncol = 4)
      }
    }
  }
  
  ## a specific yds, a specific size, a specific starting point
  if (yds != 'all' & size != 'all' & sp[1] != 'all'){
    st_dic<-list()
    st_dic[[1]]<-list()
    st_dic[[1]][[1]]<-matrix(starting_points_dic[[yds]][[size]][sp,],ncol = 4)
  }

  print(length(st_dic))
  print(length(st_dic[[1]]))
  print(dim(st_dic[[1]][[1]]))
  
  
  spatial_fieldbook_dictionary<-list()
  for (yds in 1:length(st_dic)){
    spatial_fieldbook_dictionary[[yds]]<-list()
    for (si in 1:length(st_dic[[yds]])){
      spatial_fieldbook_dictionary[[yds]][[si]]<-list()
      for (i in 1:nrow(st_dic[[yds]][[si]])){
        spatial_fieldbook_dictionary[[yds]][[si]][[i]]<-list()
        xx<-spatial_master_scheduler[which(spatial_master_scheduler$size==si),]
        print(dim(xx))
        maxrow<-max(xx$nrow)
        maxcol<-max(xx$ncol)
        for (j in 1:nrow(xx)){
          spatial_fieldbook_dictionary[[yds]][[si]][[i]][[j]]<-list()
          x0<-st_dic[[yds]][[si]][i,3]
          y0<-st_dic[[yds]][[si]][i,2]
          x1<-x0+maxcol-1
          y1<-y0+maxrow-1
          selected<-yd_monitor[[yds]][y0:y1,x0:x1]
          sel_long <- melt(selected,id.vars = rownames(selected),
                           variable.name = "c")
          colnames(sel_long)<-c("Row","Column","Yld")
          fit.ar1 <- gls(Yld ~ 1, data = sel_long,corr = corAR1(form= ~ 1|Row))
          fit.ar2 <- gls(Yld ~ 1, data = sel_long,corr = corAR1(form= ~1|Column))
          
          if (xx$rho1[j]==0){
            rho1<-as.numeric(coef(fit.ar1$modelStruct$corStruct,unconstrained = FALSE))
            rho2<-as.numeric(coef(fit.ar2$modelStruct$corStruct,unconstrained = FALSE))
          }else{
            rho1<-xx$rho1[j]
            rho2<-xx$rho2[j]
          }
          
          a1 <- rcDiGGer(numberOfTreatments=xx$trt[j],
                         rowsInDes=xx$nrow[j], columnsInDes=xx$ncol[j],
                         rowsInRep=xx$rowsInRep[j], columnsInRep=xx$colsInRep[j],
                         twoPhase = FALSE,
                         rngSeeds = c(2351, 8342), #random number generator seeds
                         maxInterchanges = 50000) 
          if (rc==FALSE){
          b1 <- Block(rowsInBlock=xx$rowsInBlock[j], columnsInBlock=xx$columnsInBlock[j],
                      blockVarianceRatio=2.0)
          
          c1 <- Correlation(rowModel="AR", rowParameter=rho1,   ##### CAMBIAR ESTO
                            columnModel="AR", columnPar=rho2,    ##### CAMBIAR ESTO
                            spatialVarianceRatio=1.0) 
          
          o1 <- Objective(weight=1.0, linearCovariate="BOTH", # options "ROWS", "COLUMNS" or "BOTH" to model linear trends associated with rows and columns
                          numberOfBlocks=1, blocks=list(b1),
                          corr=c1)
          
          p1 <- Phase(rowsInSwapBlock=xx$rowsInSwapBlock[j], columnsInSwapBlock=xx$columnsInSwapBlock[j],
                      rowsInCorrelationBlock=xx$rowsInCorrelationBlock[j],
                      columnsInCorrelationBlock=xx$columnsInCorrelationBlock[j], aType="A++",
                      targetAValue=0, maxInterchanges=500000,
                      searchIntensity=100, numberOfObjectives=1,
                      objectives=list(o1))
          
          addPhase(a1, p1)
          
          
          a1 <- run(a1)
          
          }
          tmp<-data.frame(row=a1$dlist$ROW,col=a1$dlist$RANGE,trt=a1$dlist$TRT,c_block=a1$dlist$REP,i_block=a1$dlist$B111,type=NA)
          
          fb<-do.call('rbind',replicate(xx$n_env[j],
                                        tmp,simplify = FALSE))
          
          fb$env<-rep(1:xx$n_env[j],each=nrow(tmp))
          
          spatial_fieldbook_dictionary[[yds]][[si]][[i]][[j]]<-list(yd_monitor=yds,
                                                                    size=si,
                                                                    starting_point=i,
                                                                    spatial_design=j,
                                                                    id=xx$id[j],
                                                                    x=x0,
                                                                    y=y0,
                                                                    fieldbook=fb,
                                                                    rhos=c(rho1,rho2))
        }
      }
    }
  }
  return(spatial_fieldbook_dictionary)
}


## provided maps are in coordinates, this function transform them to a matrix

ind2map<-function(y=NULL,row=NULL,col=NULL,center=TRUE){
  
  map<-matrix(NA,nrow=length(unique(col)),ncol=length(unique(row)))
  
  if (center){
    y<-y-mean(y)
  }
  
  for (i in 1:length(y)){
    map[i]<-y[i]
  }
  
  map<-t(map)
  rownames(map)<-paste0('r',1:nrow(map))
  colnames(map)<-paste0('c',1:ncol(map))
  return(map)
}


## layout creator

layoutGenerator<-function(design=NULL,n_row=NULL,n_col=NULL,n_env=NULL,trt=NULL,rep=NULL,inc_blocks=NULL,checks=NULL,plot=FALSE){
  design<-tolower(design)
  
  ## UNREP
  
  if (design=='unrep'){
    #cat('\nNOTE: unrep design is filled row-wise')
    if (n_col*n_row != trt){
      stop(paste0('\nERROR: map dimensions are rows=',n_row,', cols=',n_col,', for a total of ',n_row*n_col,' e.u. You are requesting trt=',trt,', please correct...\n'))
    }else{
      
      book1<-data.frame(row=numeric(length = trt*n_env)*NA,
                        col=numeric(length = trt*n_env)*NA,
                        trt=numeric(length = trt*n_env)*NA,
                        c_block=numeric(length = trt*n_env)*NA,
                        i_block=numeric(length = trt*n_env)*NA,
                        type=numeric(length = trt*n_env)*NA,
                        env=numeric(length=trt*n_env)*NA)
      cont<-1
      
      for (env in 1:n_env){
        
        map0<-matrix(sample(1:trt),n_row,n_col)
        for (i in 1:trt){
          f<-which(i==map0,arr.ind = TRUE)
          for (j in 1:nrow(f)){
            book1$trt[cont]<-i
            book1$row[cont]<-f[j,1]
            book1$col[cont]<-f[j,2]
            book1$env[cont]<-env
            cont<-cont+1
          }
        }
      }
      f<-which(book1$env==1)
      if (plot){
        xx<-book1[f,]
        print(ggplot(xx,aes(x=row,y=col,fill=as.factor(trt)))+geom_tile()+geom_text(aes(label = round(trt, 1))) + theme(legend.position = 'none')+coord_flip())
      }
    }
  }
  
  
  
  
  ## CRD
  
  if (design=='crd'){
    #cat('\nNOTE: crd is filled column-wise')
    if (n_col*n_row != trt*rep){
      stop(paste0('\nERROR: map dimensions are rows=',n_row,', cols=',n_col,', for a total of ',n_row*n_col,' e.u. You are requesting trt*rep=',trt*rep,', please correct...\n'))
    }else{
     
        
        book1<-data.frame(row=numeric(length = trt*rep*n_env)*NA,
                          col=numeric(length = trt*rep*n_env)*NA,
                          trt=numeric(length = trt*rep*n_env)*NA,
                          c_block=numeric(length = trt*rep*n_env)*NA,
                          i_block=numeric(length = trt*rep*n_env)*NA,
                          type=numeric(length = trt*rep*n_env)*NA,
                          env=numeric(length = trt*rep*n_env)*NA)
        cont<-1
        for (env in 1:n_env){
          map0<-matrix(sample(rep(1:trt,rep)),nrow = n_row,byrow = FALSE)
          for (i in 1:trt){
            f<-which(i==map0,arr.ind = TRUE)
            for (j in 1:nrow(f)){
              book1$trt[cont]<-i
              book1$row[cont]<-f[j,1]
              book1$col[cont]<-f[j,2]
              book1$env[cont]<-env
              cont<-cont+1
            }
          }
          
        }
        if (plot){
          f<-which(book1$env==1)
          xx<-book1[f,]
          print(ggplot(xx,aes(x=row,y=col,fill=as.factor(trt)))+geom_tile()+
                  geom_text(aes(label = round(trt, 1))) + 
                  theme(legend.position = 'none')+coord_flip())
          
        }
      }
    }

  ## ALPHA
  
  if (design=='alpha'){
    #cat('\nNOTE: alpha is filled column-wise')
    if (n_col*n_row != trt*rep){
      stop(paste0('\nERROR: map dimensions are rows=',n_row,', cols=',n_col,', for a total of ',n_row*n_col,' e.u. You are requesting trt*rep=',trt*rep,', please correct...\n'))
    }else{
      if (is.null(inc_blocks)){
        stop(paste0('\nERROR: argument <inc_blocks> is mandatory. It refers to the number of incomplete blocks within the whole experiment\n'))
      }else{
        
        
        book1<-data.frame(row=numeric(length = trt*rep*n_env)*NA,
                          col=numeric(length = trt*rep*n_env)*NA,
                          trt=numeric(length = trt*rep*n_env)*NA,
                          c_block=numeric(length = trt*rep*n_env)*NA,
                          i_block=numeric(length = trt*rep*n_env)*NA,
                          type=numeric(length = trt*rep*n_env)*NA,
                          env=numeric(length = trt*rep*n_env)*NA)
        cont<-1
        for (env in 1:n_env){
          
          
          ty<-design.alpha(trt = 1:trt, k = inc_blocks, r = rep, serie = 2, seed = 0, kinds = "Super-Duper",randomization=TRUE)
          
          map0<-cbind(ty$sketch$rep1,ty$sketch$rep2)
         
          ty2<-ty$book
          
          for (i in 1:nrow(ty2)){
            
              book1$trt[cont]<-ty2$`1:trt`[i]
              book1$c_block[cont]<-ty2$replication[i]
              book1$i_block[cont]<-ty2$block[i]
              f<-which(book1$trt[cont]==map0,arr.ind = TRUE)
              book1$row[cont]<-f[ty2$replication[i],1]
              book1$col[cont]<-f[ty2$replication[i],2]
              book1$env[cont]<-env
              cont<-cont+1
            
          }
          
        }
        if (plot){
          f<-which(book1$env==1)
          xx<-book1[f,]
          print(ggplot(xx,aes(x=row,y=col,fill=as.factor(trt)))+geom_tile()+
                  geom_text(aes(label = round(trt, 1))) + 
                  geom_hline(yintercept = seq(0.5,n_col,by = n_col/rep))+
                  geom_vline(xintercept = seq(0.5,n_row,by = n_row/inc_blocks),color='green',lty=3,lwd=2)+
                  theme(legend.position = 'none')+coord_flip())
          
        }
      }
    }
  }
  
  
  ## PREP_C
  
  if (design=='prep_c'){
    if (n_col*n_row != trt+(checks*rep)){
      stop(paste0('\nERROR: map dimensions are rows=',n_row,', cols=',n_col,', for a total of ',n_row*n_col,' e.u. You are requesting trt+checks*rep=',trt+(checks*rep),', please correct...\n'))
    }else{
      if (is.null(checks)){
        stop(paste0('\nERROR: argument <checks< is mandatory. It refers to the number of checks repeated <rep> number of times\n'))
      }else{
        
        #cat(paste0('\nNOTE: prep_c is filled column-wise. The rep/unrep proportion is ',round(checks/trt,2)))
        
        book1<-data.frame(row=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          col=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          trt=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          c_block=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          i_block=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          type=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          env=numeric(length = (trt+(rep*checks))*n_env)*NA)
        cont<-1
        
        for (env in 1:n_env){
          ct<-numeric()
          ct2<-numeric()
          for (i in 1:rep){
            ct<-c(ct,sample(c((trt+1):(trt+checks),rep(NA,trt/rep))))
            ct2<-c(ct2,rep(i,checks+(trt/rep)))
          }
          ct[which(is.na(ct))]<-sample(1:trt)
          map0<-matrix(ct,n_row,n_col,byrow = FALSE)
          bbi<-matrix(ct2,n_row,n_col,byrow = FALSE)
          
          for (i in 1:(trt+checks)){
            f<-which(i==map0,arr.ind = TRUE)
            for (j in 1:nrow(f)){
              book1$trt[cont]<-i
              book1$c_block[cont]<-bbi[f[j,1],f[j,2]]
              book1$type[cont]<-ifelse(i<=trt,'trt','check')
              book1$row[cont]<-f[j,1]
              book1$col[cont]<-f[j,2]
              book1$env[cont]<-env
              cont<-cont+1
            }
          }
        }
        if (plot){
          f<-which(book1$env==1)
          xx<-book1[f,]
          print(ggplot(xx,aes(x=row,y=col,fill=as.factor(trt)))+geom_tile()+
                  geom_text(aes(label = round(trt, 1))) + 
                  geom_hline(yintercept = seq(0.5,n_col,by = n_col/rep),lwd=3)+
                  theme(legend.position = 'none')+coord_flip())
          cat('\nLocation 1 is shown in the map')
        }
      }
    }
  }
  
  ## PREP_M
  
  if (design=='prep_m'){
    if (n_col*n_row != trt+(checks*rep)){
      stop(paste0('\nERROR: map dimensions are rows=',n_row,', cols=',n_col,', for a total of ',n_row*n_col,' e.u. You are requesting trt+checks*rep=',trt+(checks*rep),', please correct...\n'))
    }else{
      if (is.null(checks)){
        stop(paste0('\nERROR: argument <checks< is mandatory. It refers to the number of checks repeated <rep> number of times\n'))
      }else{
        
        #cat(paste0('\nNOTE: prep_m is filled column-wise. The rep/unrep proportion is ',round(checks/trt,2)))
        
        book1<-data.frame(row=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          col=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          trt=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          c_block=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          i_block=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          type=numeric(length = (trt+(rep*checks))*n_env)*NA,
                          env=numeric(length = (trt+(rep*checks))*n_env)*NA)
        cont<-1
        
        for (env in 1:n_env){
          ct<-numeric()
          ct2<-numeric()
          for (i in 1:rep){
            ct<-c(ct,sample(c((trt+1):(trt+checks),rep(NA,trt/rep))))
            ct2<-c(ct2,rep(i,checks+(trt/rep)))
          }
          ct[which(is.na(ct))]<-sample(1:trt)
          map0<-matrix(ct,n_row,n_col,byrow = FALSE)
          bbi<-matrix(ct2,n_row,n_col,byrow = FALSE)
          
          for (i in 1:(trt+checks)){
            f<-which(i==map0,arr.ind = TRUE)
            for (j in 1:nrow(f)){
              book1$trt[cont]<-i
              book1$c_block[cont]<-bbi[f[j,1],f[j,2]]
              book1$type[cont]<-ifelse(i<=trt,'trt','check')
              book1$row[cont]<-f[j,1]
              book1$col[cont]<-f[j,2]
              book1$env[cont]<-env
              cont<-cont+1
            }
          }
        }
        ## converting prep_c to prep_m
        
        
        idx<-seq(1,10000,checks)
        #idx<-seq(1,trt+checks,checks)
        tmp<-numeric()
        basek<-1:(trt+checks)
        chbase<-rep(basek,3)
        for (i in 1:env){
          ch<-idx[i]:(idx[i]+checks-1)
          ch<-chbase[ch]
          tmp<-c(tmp,c(basek[-ch],rep(ch,each=rep)))
        }
        
        book1$trt<-tmp
        
        if (plot){
          f<-which(book1$env==1)
          xx<-book1[f,]
          print(ggplot(xx,aes(x=row,y=col,fill=as.factor(trt)))+geom_tile()+
                  geom_text(aes(label = round(trt, 1))) + 
                  geom_hline(yintercept = seq(0.5,n_col,by = n_col/rep),lwd=3)+
                  theme(legend.position = 'none')+coord_flip())
          cat('\nLocation 1 is shown in the map')
        }
      }
    }
  }
  
  return(book1) 
}


## layout filler

layoutFiller<-function(fieldbook=NULL,pheno=NULL,yd_monitor=NULL,x=0,y=0,h2=0.5){
  
  
  nenv<-ncol(pheno)
  nt<-length(unique(fieldbook$trt))
  
  pheno_s<-sample(rownames(pheno),nt)
  
  fieldbook$trt_names<-pheno_s[fieldbook$trt]
  
  for (j in 1:nrow(fieldbook)){
    fieldbook$ygorro[j]<-as.numeric(pheno[fieldbook$trt_names[j],paste0('L',fieldbook$env[j])])
    fieldbook$field_noise[j]<-yd_monitor[y+fieldbook$row[j]-1,x+fieldbook$col[j]-1] 
  }
  vg <- var(fieldbook$ygorro)
  vextra <- var(fieldbook$field_noise)
  
  ##################change negative Vextra_noise ###
  vextranoise <- (2*((1-h2)/h2)*vg)-vextra
  
  if (vextranoise < 0) {vextranoise <- 0} else {vextranoise <- vextranoise}
  fieldbook$extra_noise <- rnorm(nrow(fieldbook),mean=0,sd=sqrt(vextranoise)) #replications r=2
  
   ## here we add extra noise based on h2
  fieldbook$y<-fieldbook$ygorro+
    fieldbook$field_noise+
    fieldbook$extra_noise
  
  
  return(fieldbook)
}

################################### FOR SUPERRUNNER  ############################################
myNULLmodel <- function(design,trait,data,E3, K3, soft="sommer",debug=FALSE, ge="ds",re="ds"){
  K3<<-K3
  E3<<-E3
  data[,"trtf"] <- as.character(data[,"trt"])
  data[,"rowf"] <- as.character(data[,"row"])
  data[,"colf"] <- as.character(data[,"col"])
  data[,"i_blockf"] <- as.character(data[,"i_block"])
  data[,"c_blockf"] <- as.character(data[,"c_block"])
  data[,"ci_blockf"] <- as.character(paste(data[,"c_block"],data[,"i_block"]))
  data[,"ci_block"] <- as.numeric(as.factor(data[,"ci_blockf"]))
  data[,"envf"] <- as.character(as.factor(data[,"env"]))
  ###Change - calculate the compound symetry structure for genetic effects
  data[,"EnvName"] <- paste(data$envf,data$trtf,sep = ":")
  
  # myfixedformula <- as.formula(paste(trait,"~ envf"))
  # 
  # coln <- c("c_block","i_block")
  # form <- NULL
  # for(icol in coln){
  #   levs <- table(data[,icol])
  #   if(length(levs) > 1){
  #     toadd <- paste0(" vs(",ge,"(envf,diag(5)),",icol,")")
  #     if(is.null(form)){form <- toadd}else{form <- paste(form,toadd, sep=" + ")}
  #   }
  # }
  
  basefixed <- paste0(trait,"~ envf")
  toadd_f <- NULL
  levs_c_blockf <- table(data[,"c_blockf"])
  if(length(levs_c_blockf) > 1){toadd_f <- paste0("envf:c_blockf")}
  if(!is.null(toadd_f)){basefixed <- paste(basefixed,toadd_f, sep=" + ")}
  
  myfixedformula <- as.formula(basefixed)
  
  toadd <- NULL
  levs_i_blockf <- table(data[,"i_blockf"])
  if(length(levs_i_blockf) > 1){toadd <- paste0(" vs(ds(envf), i_blockf)")}

  ###Individual trial level####
  basefixed_i <- paste0(trait,"~ 1")
  toadd_f_i <- NULL
  if(length(levs_c_blockf) > 1){toadd_f_i <- paste0("c_blockf")}
  if(!is.null(toadd_f_i)){basefixed_i <- paste(basefixed_i,toadd_f_i, sep="+")}
  myfixedformula_i <- as.formula(basefixed_i)
  
  toadd_i <- NULL
  if(length(levs_i_blockf) > 1){toadd_i <- paste0("i_blockf")}
  ###
  
  if(design == 'unrep'){
  if(soft=="asreml"){
      
    }else{
      
      myrandomformula <- as.formula(paste0("~ vs(trtf,Gu=K3)"))
      myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))

      mod.s <- try(
        mmer(fixed=myfixedformula,
             random=myrandomformula,
             rcov=myresidualformula,
             data=data, verbose = FALSE), silent = TRUE)
      if(class(mod.s) != "try-error"){
        if(length(mod.s) > 0){
          pp <- predict(mod.s, RtermsToForce = 1, FtermsToForce = 1:2)
          #pp <- predict(mod.s)
        }else{pp<-list()}
      }else{pp<-list()}
      if (class(mod.s) != "try-error"  & length(pp) > 0) {
        # print(summary(mod.s)$varcomp)
        VCP <- summary(mod.s)$varcomp
        # vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
        # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # if(length(vgmain) > 0){vg <- vg+vgmain}
        # ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
        # mve <- mean(ve)
        # h2 <-  mean(vg/(vg+ve))#average across locations
        # 
        # ######cambios
        # mvg<-mean(vg)
        # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
        # vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
        # vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
        # mvgxe<-mean(vgxe)
        # 
        vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
        vg <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
        
        mvg <- mean(vg)
        mve <- mean(ve)
        
        #h2 <-  mean(vg/(vg + ve)) #average across locations
        
        # Global H2 Cullis ####
        n.g      <- length(unique(mod.s$data$trtf))
        C22.g    <- mod.s$PevU[["u:trtf"]][[1]]       #
        trC22.g  <- sum(diag(as.matrix(C22.g)))         # trace
        vdBLUP.g <- 2/n.g*(trC22.g-(sum(C22.g)-trC22.g)/(n.g-1)) # Mean variance of a difference of two genotypic BLUPs
      
        h2 <- 1-(vdBLUP.g / 2 / mvg)
        
        h2L1 <- NA
        h2L2 <- NA
        h2L3 <- NA
        h2L4 <- NA
        h2L5 <- NA
        h2L <- NA
        mvgxe <- NA
        ratioGGE <- NA
        
        nenv <- unique(pp$fitted$env)
        withinCor <- numeric()
        for(ienv in nenv){
          useEnv <- which(pp$fitted$env == ienv)
          withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
        }
        pa <- mean(withinCor)
        bic <- mod.s$BIC
        
        ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
        orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
        orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
        orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
        orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
        best <- numeric()
        #new change realized selection gain
        rsg <- numeric()
        for(ienv in nenv){
          a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
          b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
          best[ienv] <- length(intersect(a,b))/length(a)
          rsg[ienv] <- (mean(orderPred[which(orderPred$env == ienv)[1:ngeno],"ygorro"]) - mean(orderPred[which(orderPred$env == ienv),"ygorro"]))/mean(orderPred[which(orderPred$env == ienv),"ygorro"])*100
        }
        
        mbest <- mean(best)
        mrsg <- mean(rsg)
        outp<-c(h2,h2L1,h2L2,h2L3,h2L4,h2L5,h2L,mve,pa,bic,mbest,mvg,mvgxe,ratioGGE,mrsg)
        
        ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
        
      }else{
        ifelse(debug==FALSE,return(rep(NA,15)),return(list(rep(NA,15),'model-error')))
      }
    }
  }else{
    if(soft=="asreml"){
      
    }else{
      
      # if(ge == "cs"){cs <- TRUE; ge="us"}else{cs<-FALSE}
      # baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3)")
      # if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3) ")}  #if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3)")} ##
      # if(!is.null(form)){baserandom <- paste(baserandom,form, sep=" + ")}
      #if(ge == "us"){baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3)")}
  
      if(ge == "ds"){
        baserandom <- paste0("~ vs(EnvName,Gu=kronecker(E3,K3, make.dimnames=TRUE))"," + vs(trtf,Gu=K3)")
      } else {if(ge == "us"){baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3)")}}
      
      if(!is.null(toadd)){baserandom <- paste(baserandom,toadd, sep=" + ")}
      
      myrandomformula <- as.formula(baserandom)
      myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))
      
      #Individual trial level####
      baserandom_i <- paste0("~ vs(trtf,Gu=K3)")
      if(!is.null(toadd_i)){baserandom_i <- paste(baserandom_i,toadd_i, sep=" + ")}
      myrandomformula_i <- as.formula(baserandom_i)
      myresidualformula_i <- as.formula(paste0("~ units"))
      
      h2L_i <- numeric()
      env_i <- unique(data$envf)
      for(ii in env_i){
        mod.s_i <- try(
          mmer(fixed=myfixedformula_i,
               random=myrandomformula_i,
               rcov=myresidualformula_i,
               data=data[data$envf==ii,], verbose = FALSE), silent = TRUE)
        if(class(mod.s_i) != "try-error"){
          if(length(mod.s_i) > 0){
        n.g_i      <- length(unique(data$trtf))  # number of genotypes
        vc.g_i     <- as.numeric(mod.s_i$sigma[["u:trtf"]]) # genetic variance component
        C22.g_i    <- mod.s_i$PevU[["u:trtf"]][[1]]       #
        trC22.g_i  <- sum(diag(as.matrix(C22.g_i)))         # trace
        vdBLUP.g_i <- 2/n.g_i*(trC22.g_i-(sum(C22.g_i)-trC22.g_i)/(n.g_i-1)) # Mean variance of a difference of two genotypic BLUPs
        # H2 Cullis at individual trial level ###
        h2L_i[ii] <- 1-(vdBLUP.g_i / 2 / vc.g_i)
          }
        }
      }
      h2L1 <- h2L_i[1]
      h2L2 <- h2L_i[2]
      h2L3 <- h2L_i[3]
      h2L4 <- h2L_i[4]
      h2L5 <- h2L_i[5]
      h2L <- mean(h2L_i,na.rm=T)
      ###
      
      mod.s <- try(
        mmer(fixed=myfixedformula,
             random=myrandomformula,
             rcov=myresidualformula,
             data=data, verbose = FALSE), silent = TRUE)
      if(class(mod.s) != "try-error"){
        if(length(mod.s) > 0){
          pp <- predict(mod.s, RtermsToForce = 1:2, FtermsToForce = 1:2)
          #pp <- predict(mod.s)
        }else{pp<-list()}
      }else{pp<-list()}
      if (class(mod.s) != "try-error"  & length(pp) > 0) {
        # print(summary(mod.s)$varcomp)
        VCP <- summary(mod.s)$varcomp
        # vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
        # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # if(length(vgmain) > 0){vg <- vg+vgmain}
        # ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
        # mve <- mean(ve)
        # h2 <-  mean(vg/(vg+ve))#average across locations
        # 
        # ######cambios##
        # mvg<-mean(vg)
        # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
        # vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
        # vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
        # mvgxe<-mean(vgxe)
        # 
        vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
        mvg <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        mvgxe<-mean(VCP[intersect(grep("u:EnvName",rownames(VCP)),vcpIndex),"VarComp"])
        
        #vg_gxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
        #mvg_ge <- mean(vg_gxe)
        mve <- mean(ve)
        
        #h2 <- 1 - (sum(diag(mod.s$PevU$`u:trtf`[[trait]]))/nrow(mod.s$PevU$`u:trtf`[[trait]])/(mvg))
        # h2 <-  mean(vg_gxe/(vg_gxe + ve/2)) #r=2 
        #h2 <-  mean((mvg+mvgxe)/(mvg + mvgxe + ve/2)) #r=2 #average across locations
        
        # Global H2 Cullis ####
        n.g      <- length(unique(mod.s$data$trtf))
        C22.g    <- mod.s$PevU[["u:trtf"]][[1]]
        trC22.g  <- sum(diag(as.matrix(C22.g)))         # trace
        vdBLUP.g <- 2/n.g*(trC22.g-(sum(C22.g)-trC22.g)/(n.g-1)) # Mean variance of a difference of two genotypic BLUPs
        h2 <- 1-(vdBLUP.g / 2 / mvg)        
        
        # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
        # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
        # mvg <- mean(vg)
        # mvgxe <- mvg_ge - mvg
        ratioGGE <- mvgxe/mvg
        
        nenv <- unique(pp$fitted$env)
        withinCor <- numeric()
        for(ienv in nenv){
          useEnv <- which(pp$fitted$env == ienv)
          withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
        }
        pa <- mean(withinCor)
        bic <- mod.s$BIC
        
        ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
        orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
        orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
        orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
        orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
        best <- numeric()
        #new change realized selection gain
        rsg <- numeric()
        for(ienv in nenv){
          a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
          b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
          best[ienv] <- length(intersect(a,b))/length(a)
          rsg[ienv] <- (mean(orderPred[which(orderPred$env == ienv)[1:ngeno],"ygorro"]) - mean(orderPred[which(orderPred$env == ienv),"ygorro"]))/mean(orderPred[which(orderPred$env == ienv),"ygorro"])*100
        }
        
        mbest <- mean(best)
        mrsg <- mean(rsg)
        outp<-c(h2,h2L1,h2L2,h2L3,h2L4,h2L5,h2L,mve,pa,bic,mbest,mvg,mvgxe,ratioGGE,mrsg)
        
        ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
        
      }else{
        ifelse(debug==FALSE,return(rep(NA,15)),return(list(rep(NA,15),'model-error')))
      }
    } 
  }
}

# my2Dmodel <- function(trait,data, K3, soft="sommer",debug=FALSE, ge="ds",re="ds"){
#   K3<<-K3
#   data[,"trtf"] <- as.character(data[,"trt"])
#   data[,"rowf"] <- as.character(data[,"row"])
#   data[,"colf"] <- as.character(data[,"col"])
#   data[,"i_blockf"] <- as.character(data[,"i_block"])
#   data[,"c_blockf"] <- as.character(data[,"c_block"])
#   data[,"ci_blockf"] <- as.character(paste(data[,"c_block"],data[,"i_block"]))
#   data[,"ci_block"] <- as.numeric(as.factor(data[,"ci_blockf"]))
#   data[,"envf"] <- as.character(as.factor(data[,"env"]))
#   
#   myfixedformula <- as.formula(paste(trait,"~ envf"))
#   
#   coln <- c("c_block","i_block")
#   form <- NULL
#   for(icol in coln){
#     levs <- table(data[,icol])
#     if(length(levs) > 1){
#       toadd <- paste0(" vs(",ge,"(envf,diag(5)),",icol,")")
#       if(is.null(form)){form <- toadd}else{form <- paste(form,toadd, sep=" + ")}
#     }
#   }
#   
#   if(soft=="asreml"){
#     
#   }else{
#     if(ge == "cs"){cs <- TRUE; ge="us"}else{cs<-FALSE}
#     baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3) + vs(ds(envf),rowf) + vs(ds(envf),colf) + vs(ds(envf),spl2D(row,col))")  # 
#     if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3) ")}
#     
#     if(!is.null(form)){baserandom <- paste(baserandom,form, sep=" + ")}
#     myrandomformula <- as.formula(baserandom)
#     myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))
#     print(myrandomformula)
#     mod.s <- try(
#       mmer(fixed=myfixedformula,
#            random=myrandomformula,
#            rcov=myresidualformula,
#            data=data, verbose = FALSE), silent = TRUE)
#     if(class(mod.s) != "try-error"){
#       if(length(mod.s) > 0){
#         pp <- predict(mod.s)
#       }else{pp<-list()}
#     }else{pp<-list()}
#     if (class(mod.s) != "try-error"  & length(pp) > 0) {
#       VCP <- summary(mod.s)$varcomp
#       vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
#       vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       if(length(vgmain) > 0){vg <- vg+vgmain}
#       ve <- abs(VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"])
#       mve <- mean(ve)
#       h2 <-  mean(vg/(vg+ve))#average across locations
#       ######cambios###
#       mvg<-mean(vg)
#       vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
#       vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
#       vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
#       mvgxe<-mean(vgxe)
#       #
#       nenv <- unique(pp$fitted$env)
#       withinCor <- numeric()
#       for(ienv in nenv){
#         useEnv <- which(pp$fitted$env == ienv)
#         withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
#       }
#       pa <- mean(withinCor)
#       bic <- mod.s$BIC
#       
#       ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
#       orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
#       orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
#       orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
#       orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
#       best <- numeric()
#       for(ienv in nenv){
#         a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
#         b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
#         best[ienv] <- length(intersect(a,b))/length(a)
#       }
#       mbest <- mean(best)
#       outp<-c(h2,mve,pa,bic,mbest,mvg,mvgxe)
#       # h2 <- 1 - (mean(mod.s$PevU$`u:trtf`[[trait]])/(vg))
#       ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
#     }else{
#       ifelse(debug==FALSE,return(rep(NA,7)),return(list(rep(NA,7),'model-error')))
#     }
#   }
# }
# 
# myAR1model <- function(trait,data, rho,direction,K3,soft="sommer",debug=FALSE, ge="ds",re="ds"){
#   K3<<-K3
#   data[,"trtf"] <- as.character(data[,"trt"])
#   data[,"rowf"] <- as.character(data[,"row"])
#   data[,"colf"] <- as.character(data[,"col"])
#   data[,"i_blockf"] <- as.character(data[,"i_block"])
#   data[,"c_blockf"] <- as.character(data[,"c_block"])
#   data[,"ci_blockf"] <- as.character(paste(data[,"c_block"],data[,"i_block"]))
#   data[,"ci_block"] <- as.numeric(as.factor(data[,"ci_blockf"]))
#   data[,"envf"] <- as.character(as.factor(data[,"env"]))
#   
#   oposite <- setdiff(c("rowf","colf"),c(direction,paste0(direction,"f")))
#   touse <- intersect(c("rowf","colf"),c(direction,paste0(direction,"f")))
#   
#   myfixedformula <- as.formula(paste(trait,"~ envf"))
#   myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))
#   
#   
#   coln <- c("c_block","i_block")
#   form <- NULL
#   for(icol in coln){
#     levs <- table(data[,icol])
#     if(length(levs) > 1){
#       toadd <- paste0("vs(",ge,"(envf,diag(5)),",icol,")")
#       if(is.null(form)){form <- toadd}else{form <- paste(form,toadd, sep=" + ")}
#     }
#   }
#   
#   if(soft=="asreml"){
#     
#   }else{
#     
#     if(ge == "cs"){cs <- TRUE; ge="us"}else{cs<-FALSE}
#     baserandom <- paste0("~  vs(",ge,"(envf), trtf,Gu=K3) + vs(ds(envf),",touse,", Gu=AR1(",touse,", rho=",rho,")) + vs(ds(envf),",oposite,")")
#     if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3)")}
#     
#     if(!is.null(form)){baserandom <- paste(baserandom,form, sep=" + ")}
#     myrandomformula <- as.formula(baserandom)
#     print(myrandomformula)
#     mod.s <- try(
#       mmer(fixed=myfixedformula,
#            random=myrandomformula,
#            rcov=myresidualformula,
#            data=data, verbose = FALSE), silent = TRUE)
#     if(class(mod.s) != "try-error"){
#       if(length(mod.s) > 0){
#         pp <- predict(mod.s)
#       }else{pp<-list()}
#     }else{pp<-list()}
#     if (class(mod.s) != "try-error"  & length(pp) > 0) {
#       VCP <- summary(mod.s)$varcomp
#       vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
#       vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       if(length(vgmain) > 0){vg <- vg+vgmain}
#       ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
#       mve <- mean(ve)
#       h2 <-  mean(vg/(vg+ve))#average across locations
#       ######cambios##
#       mvg<-mean(vg)
#       vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
#       vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
#       vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
#       mvgxe<-mean(vgxe)
#       
#       nenv <- unique(pp$fitted$env)
#       withinCor <- numeric()
#       for(ienv in nenv){
#         useEnv <- which(pp$fitted$env == ienv)
#         withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
#       }
#       pa <- mean(withinCor)
#       bic <- mod.s$BIC
#       
#       
#       ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
#       orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
#       orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
#       orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
#       orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
#       best <- numeric()
#       for(ienv in nenv){
#         a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
#         b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
#         best[ienv] <- length(intersect(a,b))/length(a)
#       }
#       mbest <- mean(best)
#       outp<-c(h2,mve,pa,bic,mbest,mvg,mvgxe)
#       # h2 <- 1 - (mean(mod.s$PevU$`u:trtf`[[trait]])/(vg))
#       ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
#     }else{
#       ifelse(debug==FALSE,return(rep(NA,7)),return(list(rep(NA,7),'model-error')))
#     }
#   }
#   
# }

myAR1AR1model <- function(design,trait,data, rho1, rho2,E3,K3,soft="sommer",debug=FALSE, ge="ds",re="ds"){
  K3<<-K3
  E3<<-E3
  data[,"trtf"] <- as.character(data[,"trt"])
  data[,"rowf"] <- as.character(data[,"row"])
  data[,"colf"] <- as.character(data[,"col"])
  data[,"i_blockf"] <- as.character(data[,"i_block"])
  data[,"c_blockf"] <- as.character(data[,"c_block"])
  data[,"ci_blockf"] <- as.character(paste(data[,"c_block"],data[,"i_block"]))
  data[,"ci_block"] <- as.numeric(as.factor(data[,"ci_blockf"]))
  data[,"envf"] <- as.character(as.factor(data[,"env"]))
  data[,"EnvName"] <- paste(data$envf,data$trtf,sep = ":")
  
  # myfixedformula <- as.formula(paste(trait,"~ envf"))
  # 
  # coln <- c("c_block","i_block")
  # form <- NULL
  # for(icol in coln){
  #   levs <- table(data[,icol])
  #   if(length(levs) > 1){
  #     toadd <- paste0(" vs(",ge,"(envf,diag(5)),",icol,")")
  #     if(is.null(form)){form <- toadd}else{form <- paste(form,toadd, sep=" + ")}
  #   }
  # }
  basefixed <- paste0(trait,"~ envf")
  toadd_f <- NULL
  levs_c_blockf <- table(data[,"c_blockf"])
  if(length(levs_c_blockf) > 1){toadd_f <- paste0("envf:c_blockf")}
  if(!is.null(toadd_f)){basefixed <- paste(basefixed,toadd_f, sep=" + ")}
  
  myfixedformula <- as.formula(basefixed)
  
  toadd <- NULL
  levs_i_blockf <- table(data[,"i_blockf"])
  if(length(levs_i_blockf) > 1){toadd <- paste0("vs(ds(envf), i_blockf)")}

  ##Individual trial level####
  basefixed_i <- paste0(trait,"~ 1")
  toadd_f_i <- NULL
  if(length(levs_c_blockf) > 1){toadd_f_i <- paste0("c_blockf")}
  if(!is.null(toadd_f_i)){basefixed_i <- paste(basefixed_i,toadd_f_i, sep="+")}
  myfixedformula_i <- as.formula(basefixed_i)

  toadd_i <- NULL
  if(length(levs_i_blockf) > 1){toadd_i <- paste0("i_blockf")}
  ##
  
  if(design == 'unrep'){
    if(soft=="asreml"){
      
    }else{

      myrandomformula <- as.formula(paste0("~ vs(trtf,Gu=K3) + vs(ds(envf),rowf,Gu=AR1(rowf, rho=",rho1,")) + vs(ds(envf),colf, Gu=AR1(colf, rho=",rho2,"))"))
      myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))

      print(myrandomformula)
      mod.s <- try(
        mmer(fixed=myfixedformula,
             random=myrandomformula,
             rcov=myresidualformula,
             data=data, verbose = FALSE), silent = TRUE)
      if(class(mod.s) != "try-error"){
        if(length(mod.s) > 0){
          pp <- predict(mod.s, RtermsToForce = 1, FtermsToForce = 1:2)
          #pp <- predict(mod.s)
        }else{pp<-list()}
      }else{pp<-list()}
      if (class(mod.s) != "try-error"  & length(pp) > 0) {
        # print(summary(mod.s)$varcomp)
        VCP <- summary(mod.s)$varcomp
        # vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
        # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # if(length(vgmain) > 0){vg <- vg+vgmain}
        # ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
        # mve <- mean(ve)
        # h2 <-  mean(vg/(vg+ve))#average across locations
        # 
        # ######cambios###
        # mvg<-mean(vg)
        # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
        # vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
        # vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
        # mvgxe<-mean(vgxe)
        # ##
        # vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
        # vg <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
        vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
        vg <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]

        mvg <- mean(vg)
        mve <- mean(ve)
        
        #h2 <-  mean(vg/(vg + ve)) #average across locations
        
        # Global H2 Cullis ###
        n.g      <- length(unique(mod.s$data$trtf))
        C22.g    <- mod.s$PevU[["u:trtf"]][[1]]       #
        trC22.g  <- sum(diag(as.matrix(C22.g)))         # trace
        vdBLUP.g <- 2/n.g*(trC22.g-(sum(C22.g)-trC22.g)/(n.g-1)) # Mean variance of a difference of two genotypic BLUPs
        
        h2 <- 1-(vdBLUP.g / 2 / mvg)
        
        h2L1 <- NA
        h2L2 <- NA
        h2L3 <- NA
        h2L4 <- NA
        h2L5 <- NA
        h2L <- NA
        mvgxe <- NA
        ratioGGE <- NA
        
        nenv <- unique(pp$fitted$env)
        withinCor <- numeric()
        for(ienv in nenv){
          useEnv <- which(pp$fitted$env == ienv)
          withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
        }
        pa <- mean(withinCor)
        bic <- mod.s$BIC
        
        ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
        orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
        orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
        orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
        orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
        best <- numeric()
        #new change realized selection gain
        rsg <- numeric()
        for(ienv in nenv){
          a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
          b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
          best[ienv] <- length(intersect(a,b))/length(a)
          rsg[ienv] <- (mean(orderPred[which(orderPred$env == ienv)[1:ngeno],"ygorro"]) - mean(orderPred[which(orderPred$env == ienv),"ygorro"]))/mean(orderPred[which(orderPred$env == ienv),"ygorro"])*100
        }
        
        mbest <- mean(best)
        mrsg <- mean(rsg)
        outp<-c(h2,h2L1,h2L2,h2L3,h2L4,h2L5,h2L,mve,pa,bic,mbest,mvg,mvgxe,ratioGGE,mrsg)
        # h2 <- 1 - (mean(mod.s$PevU$`u:trtf`[[trait]])/(vg))
        ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
        
      }else{
        ifelse(debug==FALSE,return(rep(NA,15)),return(list(rep(NA,15),'model-error')))
      }
    }
  }else{
    if(soft=="asreml"){
      
    }else{
      # if(ge == "cs"){cs <- TRUE; ge="us"}else{cs<-FALSE}
      #baserandom <- paste0("~  vs(",ge,"(envf), trtf,Gu=K3) + vs(ds(envf),rowf, Gu=AR1(rowf, rho=",rho1,")) + vs(ds(envf),colf, Gu=AR1(colf, rho=",rho2,"))")
      # if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3)")}
      # 
      # if(!is.null(form)){baserandom <- paste(baserandom,form, sep=" + ")}
      #if(ge == "us"){baserandom <- paste0("~  vs(",ge,"(envf), trtf,Gu=K3) + vs(ds(envf),rowf,Gu=AR1(rowf, rho=",rho1,")) + vs(ds(envf),colf, Gu=AR1(colf, rho=",rho2,"))")}
      if(ge == "ds"){
        baserandom <- paste0("~ vs(EnvName,Gu=kronecker(E3,K3, make.dimnames=TRUE)) + vs(trtf,Gu=K3) + vs(ds(envf),rowf,Gu=AR1(rowf, rho=",rho1,")) + vs(ds(envf),colf, Gu=AR1(colf, rho=",rho2,"))")
      } else {if(ge == "us"){baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3)")}}
      
      if(!is.null(toadd)){baserandom <- paste(baserandom,toadd, sep=" + ")}
      
      myrandomformula <- as.formula(baserandom)
      myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))
      
      ##Individual trial level####
      baserandom_i <- paste0("~ vs(trtf,Gu=K3) + vs(rowf,Gu=AR1(rowf, rho=",rho1,")) + vs(colf, Gu=AR1(colf, rho=",rho2,"))")
      if(!is.null(toadd_i)){baserandom_i <- paste(baserandom_i,toadd_i, sep=" + ")}
      myrandomformula_i <- as.formula(baserandom_i)
      myresidualformula_i <- as.formula(paste0("~ units"))

      h2L_i <- numeric()
      env_i <- unique(data$envf)
      for(ii in env_i){
        mod.s_i <- try(
          mmer(fixed=myfixedformula_i,
               random=myrandomformula_i,
               rcov=myresidualformula_i,
               data=data[data$envf==ii,], verbose = FALSE), silent = TRUE)
        if(class(mod.s_i) != "try-error"){
          if(length(mod.s_i) > 0){
            n.g_i      <- length(unique(data$trtf))  # number of genotypes
            vc.g_i     <- as.numeric(mod.s_i$sigma[["u:trtf"]]) # genetic variance component
            C22.g_i    <- mod.s_i$PevU[["u:trtf"]][[1]]       #
            trC22.g_i  <- sum(diag(as.matrix(C22.g_i)))         # trace
            vdBLUP.g_i <- 2/n.g_i*(trC22.g_i-(sum(C22.g_i)-trC22.g_i)/(n.g_i-1)) # Mean variance of a difference of two genotypic BLUPs
            # H2 Cullis at individual trial level ###
            h2L_i[ii] <- 1-(vdBLUP.g_i / 2 / vc.g_i)
          }
        }
      }
      h2L1 <- h2L_i[1]
      h2L2 <- h2L_i[2]
      h2L3 <- h2L_i[3]
      h2L4 <- h2L_i[4]
      h2L5 <- h2L_i[5]
      h2L <- mean(h2L_i,na.rm=T)
      ##
      
      mod.s <- try(
        mmer(fixed=myfixedformula,
             random=myrandomformula,
             rcov=myresidualformula,
             data=data, verbose = FALSE), silent = TRUE)
      if(class(mod.s) != "try-error"){
        if(length(mod.s) > 0){
          #       pp <- predict(mod.s)
          pp <- predict(mod.s, RtermsToForce = 1:2, FtermsToForce = 1:2)
        }else{pp<-list()}
      }else{pp<-list()}
      if (class(mod.s) != "try-error"  & length(pp) > 0) {
        VCP <- summary(mod.s)$varcomp
        # vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
        # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # if(length(vgmain) > 0){vg <- vg+vgmain}
        # ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
        # mve <- mean(ve)
        # h2 <-  mean(vg/(vg+ve))#average across locations
        # ######cambios###
        # mvg<-mean(vg)
        # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
        # vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
        # vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        # if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
        # mvgxe<-mean(vgxe)
        ####
        vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
        mvg <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        mvgxe<-VCP[intersect(grep("u:EnvName",rownames(VCP)),vcpIndex),"VarComp"]
        # vg_gxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
        ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
        # mvg_ge <- mean(vg_gxe)
        mve <- mean(ve)
        
        #h2 <- 1 - (sum(diag(mod.s$PevU$`u:trtf`[[trait]]))/nrow(mod.s$PevU$`u:trtf`[[trait]])/(mvg))
        #h2 <-  mean(vg_gxe/(vg_gxe + ve/2)) #r=2 
        #h2 <-  mean((mvg+mvgxe)/(mvg + mvgxe + ve/2)) #r=2 
        
        # Global H2 Cullis ###
        n.g      <- length(unique(mod.s$data$trtf))
        C22.g    <- mod.s$PevU[["u:trtf"]][[1]]
        trC22.g  <- sum(diag(as.matrix(C22.g)))         # trace
        vdBLUP.g <- 2/n.g*(trC22.g-(sum(C22.g)-trC22.g)/(n.g-1)) # Mean variance of a difference of two genotypic BLUPs
        h2 <- 1-(vdBLUP.g / 2 / mvg)        
        
        # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
        # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
        # mvg <- mean(vg)
        # mvgxe <- mvg_ge - mvg
        ratioGGE <- mvgxe/mvg
        
        nenv <- unique(pp$fitted$env)
        withinCor <- numeric()
        for(ienv in nenv){
          useEnv <- which(pp$fitted$env == ienv)
          withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
        }
        pa <- mean(withinCor)
        bic <- mod.s$BIC
        
        ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
        orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
        orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
        orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
        orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
        best <- numeric()
        #new change realized selection gain
        rsg <- numeric()
        for(ienv in nenv){
          a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
          b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
          best[ienv] <- length(intersect(a,b))/length(a)
          rsg[ienv] <- (mean(orderPred[which(orderPred$env == ienv)[1:ngeno],"ygorro"]) - mean(orderPred[which(orderPred$env == ienv),"ygorro"]))/mean(orderPred[which(orderPred$env == ienv),"ygorro"])*100
        }
        
        mbest <- mean(best)
        mrsg <- mean(rsg)
        outp<-c(h2,h2L1,h2L2,h2L3,h2L4,h2L5,h2L,mve,pa,bic,mbest,mvg,mvgxe,ratioGGE,mrsg)
        #outp<-c(h2,mve,pa,bic,mbest,mvg,mvgxe)
        # h2 <- 1 - (mean(mod.s$PevU$`u:trtf`[[trait]])/(vg))
        ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
        
      }else{
        ifelse(debug==FALSE,return(rep(NA,15)),return(list(rep(NA,15),'model-error')))
      }
    }  
  } 
}


## core runner
coreRunner<-function(design = NULL,fieldbookInput=NULL,E2=NULL,K2=NULL,debug=FALSE,ge="ds",re="ds",step_ar1a=0.1,step_ar1ar1a=0.1,ncore=0){   ## input must be a list with 4 experimental designs
  
  K2<<-K2
  E2<<-E2
  
  spatialModels<-c('none','spl2D','AR1','AR1xAR1')
  
  bucket<-data.frame(spatial=numeric(length(4)),
                     h2G=numeric(length(4)),
                     h2L1=numeric(length(4)),
                     h2L2=numeric(length(4)),
                     h2L3=numeric(length(4)),
                     h2L4=numeric(length(4)),
                     h2L5=numeric(length(4)),
                     h2L=numeric(length(4)),
                     ve=numeric(length(4)),
                     pa=numeric(length(4)),
                     bic=numeric(length(4)),
                     best=numeric(length(4)),
                     vg=numeric(length(4)),
                     vgxe=numeric(length(4)),
                     ratioGE=numeric(length(4)),
                     rsg=numeric(length(4)),
                     rho1=numeric(length(4))*NA,
                     rho2=numeric(length(4))*NA
                     
  )
  
  
  model_bucket<-list()
  
  
  ## no adjustment
  
  cat('\n**running with no spatial adjustment...')
  bucket[1,1]<-spatialModels[1]
  if (debug){
    invisible(capture.output(tmpm<-try(myNULLmodel(design = design, trait="y", data=fieldbookInput,E3=E2, K3=K2,debug=TRUE,ge=ge,re=re))))
    #tmpm<-try(myNULLmodel(trait="y", data=fieldbookInput, K=K,debug=TRUE,ge=ge,re=re))
    
    if(class(tmpm) != "try-error"){
      bucket[1,2:16]<-tmpm[[1]]
      model_bucket[[1]]<-tmpm[[2]]}else{
        bucket[1,2:16]<-rep(NA,15)
      }
  }else{
    invisible(capture.output(rtrt<-try(myNULLmodel(design = design,trait="y", data=fieldbookInput,E3=E2, K3=K2,debug=FALSE,ge=ge,re=re))))
    if(class(rtrt) != "try-error"){
      bucket[1,2:16]<-rtrt
    }else{
      bucket[1,2:16]<-rep(NA,15)
    }
  }
  cat(' done!')
  
  ## spl2D models
  # cat('\n**running spl2D correction...')
  # bucket[2,1]<-spatialModels[2]
  # if (debug){
  #   invisible(capture.output(tmpm<-try(my2Dmodel(trait="y", data=fieldbookInput, K3=K2,debug = TRUE,ge=ge,re=re))))
  #   #tmpm<-try(my2Dmodel(trait="y", data=fieldbookInput, K3=K2,debug = TRUE,ge=ge,re=re))
  #   if(class(tmpm) != "try-error"){
  #     bucket[2,2:8]<-tmpm[[1]]
  #     model_bucket[[2]]<-tmpm[[2]]}else{
  #       bucket[2,2:8]<-rep(NA,7)
  #     }
  #   
  # }else{
  #   invisible(capture.output(rtrt<-my2Dmodel(trait="y", data=fieldbookInput, K3=K2,debug = FALSE,ge=ge,re=re)))
  #   if(class(rtrt) != "try-error"){
  #     
  #     bucket[2,2:8]<-rtrt
  #   }else{
  #     bucket[2,2:8]<-rep(NA,7)
  #     
  #   }
  #   
  # }
  # cat(' done!')
  # 
  # 
  # 
  # ## AR1 models
  # ars <- seq(-1,1,step_ar1a)
  # cat('\n**running AR1 adjustment...\n')
  # h2ars <- matrix(NA,nrow=length(ars),ncol=7);colnames(h2ars) <- c("h2","ve","pa","bic","best","vg","vgxe"); rownames(h2ars) <- ars
  # 
  # if (debug){
  #   arsmodel<-list()
  #   
  #   if (ncore>0){
  #     n.cores <- ncore 
  #     cl <- snow::makeCluster(n.cores)
  #     registerDoSNOW(cl)
  #     
  #     arsmodel<-foreach(k= 1:length(ars),.packages = 'sommer',.export = 'myAR1model') %dopar% myAR1model(trait="y", data=fieldbookInput,rho=ars[k], direction = "col", K3=K2,debug = TRUE,ge=ge,re=re)
  #     snow::stopCluster(cl)
  #     
  #     for (k in 1:length(arsmodel)){
  #       if(class(arsmodel[[k]]) != "try-error"){
  #         h2ars[k,]<- arsmodel[[k]][[1]]
  #       }else{
  #         h2ars[k,]<- rep(NA,7)
  #         
  #       }
  #     }
  #     cat('done!\n')
  #     
  #   }else{
  #     
  #     for(k in 1:length(ars)){
  #       cat(paste0('      rho=',ars[k],'... '))
  #       invisible(capture.output(arsmodel[[k]]<-try(myAR1model(trait="y", data=fieldbookInput,rho=ars[k], direction = "col", K3=K2,debug = TRUE,ge=ge,re=re))))
  #       if(class(arsmodel[[k]]) != "try-error"){
  #         h2ars[k,]<- arsmodel[[k]][[1]]
  #       }else{
  #         h2ars[k,]<- rep(NA,7)
  #         
  #       }
  #       cat('done!\n')
  #     }
  #     
  #   }
  #   if(length(which(!is.na(h2ars))) > 0){
  #     bucket[3,1]<-spatialModels[3]
  #     maxi<-which(h2ars[,3] == max(h2ars[,3], na.rm = TRUE))
  #     bucket[3,2:8]<-h2ars[maxi[1],]
  #     bucket[3,9]<-ars[maxi]
  #     model_bucket[[3]]<-arsmodel[[maxi[1]]][[2]]
  #   }else{
  #     bucket[3,1]<-spatialModels[3]
  #     bucket[3,2:8]<-rep(NA,7)
  #     model_bucket[[3]]<-'error'
  #   }
  #   
  # }else{
  #   
  #   
  #   
  #   if (ncore>0){
  #     n.cores <- ncore 
  #     cl <- snow::makeCluster(n.cores)
  #     registerDoSNOW(cl)
  #     
  #     rtrt<-foreach(k= 1:length(ars),.packages = 'sommer',.export = 'myAR1model') %dopar% myAR1model(trait="y", data=fieldbookInput,rho=ars[k], direction = "col", K3=K2,debug = FALSE,ge=ge,re=re)
  #     snow::stopCluster(cl)
  #     
  #     for(k in 1:length(rtrt)){
  #       if(class(rtrt[[k]]) != "try-error"){
  #         h2ars[k,] <- rtrt[[k]]
  #       }else{
  #         h2ars[k,]<-rep(NA,7)
  #       }
  #       cat('done!\n')
  #     }
  #     
  #   }else{
  #     
  #     
  #     for(k in 1:length(ars)){
  #       cat(paste0('      rho=',ars[k],'... '))
  #       invisible(capture.output(rtrt <- try(myAR1model(trait="y", data=fieldbookInput,rho=ars[k], direction = "col", K3=K2,debug=FALSE, ge=ge,re=re))))
  #       if(class(rtrt) != "try-error"){
  #         h2ars[k,] <- rtrt
  #       }else{
  #         h2ars[k,]<-rep(NA,7)
  #       }
  #       cat('done!\n')
  #     }
  #     
  #   }
  #   
  #   
  #   # store the one that gave the maximum PA
  #   if(length(which(!is.na(h2ars))) > 0){
  #     bucket[3,1]<-spatialModels[3]
  #     maxi<-which(h2ars[,3] == max(h2ars[,3], na.rm = TRUE))
  #     bucket[3,2:8]<-h2ars[maxi,]
  #     bucket[3,9]<-ars[maxi]
  #   }else{
  #     bucket[3,1]<-spatialModels[3]
  #     bucket[3,2:8]<-rep(NA,7)
  #     model_bucket[[3]]<-'error'
  #   }
  #   
  # }
  
  ## AR1xAR1 models
  ars <- seq(-1,1,step_ar1ar1a)
  ars2 <- expand.grid(ars,ars)
  cat('**running AR1AR1 adjustment...\n')
  h2ars2 <- matrix(NA,nrow=nrow(ars2),ncol=15);colnames(h2ars2) <- c("h2G","h2L1","h2L2","h2L3","h2L4","h2L5","h2L","ve","pa","bic","best","vg","vgxe","ratioGE","mrsg"); rownames(h2ars2) <- apply(ars2,1,function(x){paste(x,collapse = "-")})
  
  if (debug){
    
    if (ncore>0){
      n.cores <- ncore 
      cl <- snow::makeCluster(n.cores)
      registerDoSNOW(cl)
      
      arsmodel<-foreach(k= 1:nrow(ars2),.packages = 'sommer',.export = 'myAR1AR1model') %dopar% myAR1AR1model(design=design,trait="y", data=fieldbookInput,rho1=ars2[k,1],rho2 = ars2[k,2],E3=E2,K3=K2,debug=TRUE, ge=ge,re=re)
      snow::stopCluster(cl)
      
      for (k in 1:length(arsmodel)){
        if(class(arsmodel[[k]]) != "try-error"){
          h2ars2[k,]<- arsmodel[[k]][[1]]
        }else{
          h2ars2[k,]<- rep(NA,15)
          
        }
      }
      cat('done!\n')
      
    }else{
      
      arsmodel<-list()
      for(k in 1:nrow(ars2)){
        cat(paste0('      rho1=',format(round(ars2[k,1],2),nsmall = 2),', rho2=',format(round(ars2[k,2],2),nsmall = 2),'... '))
        invisible(capture.output(arsmodel[[k]]<-try(myAR1AR1model(design=design,trait="y", data=fieldbookInput,rho1=ars2[k,1],rho2 = ars2[k,2],E3=E2,K3=K2,debug=TRUE, ge=ge,re=re))))
        if(class(arsmodel[[k]]) != "try-error"){
          h2ars2[k,]<- arsmodel[[k]][[1]]
        }else{
          h2ars2[k,]<- rep(NA,15)
          
        }
        cat('done!\n')
      }
      
    }
    if(length(which(!is.na(h2ars2))) > 0){
      bucket[4,1]<-spatialModels[4]
      maxi<-which(h2ars2[,9] == max(h2ars2[,9], na.rm = TRUE))[1]
      #maxi<-which(h2ars2[,3] == max(h2ars2[,3], na.rm = TRUE))[1]
      bucket[4,2:16]<-h2ars2[maxi,] #max(h2ars2)
      bucket[4,17:18]<-ars2[maxi,]
      model_bucket[[4]]<-arsmodel[[maxi]][[2]]
    }else{
      bucket[4,1]<-spatialModels[4]
      bucket[4,2:16]<-rep(NA,15)
      model_bucket[[4]]<-'error'
    }
  }else{
    
    
    if (ncore>0){
      n.cores <- ncore 
      cl <- snow::makeCluster(n.cores)
      registerDoSNOW(cl)
      
      rtrt<-foreach(k= 1:length(ars),.packages = 'sommer',.export = 'myAR1AR1model') %dopar% myAR1AR1model(design=design,trait="y", data=fieldbookInput,rho1=ars2[k,1],rho2 = ars2[k,2],E3=E2,K3=K2,debug = FALSE,ge=ge,re=re)
      snow::stopCluster(cl)
      
      for(k in 1:length(rtrt)){
        if(class(rtrt[[k]]) != "try-error"){
          h2ars2[k,] <- rtrt[[k]]
        }else{
          h2ars2[k,]<-rep(NA,15)
        }
        cat('done!\n')
      }
      
    }else{
      
      for(k in 1:nrow(ars2)){
        cat(paste0('      rho1=',format(round(ars2[k,1],2),nsmall = 2),', rho2=',format(round(ars2[k,2],2),nsmall = 2),'... '))
        invisible(capture.output(rtrt <- try(myAR1AR1model(design=design,trait="y", data=fieldbookInput,rho1=ars2[k,1],rho2 = ars2[k,2],E3=E2,K3=K2,debug = FALSE,ge=ge,re=re))))
        if(class(rtrt) != "try-error"){
          h2ars2[k,] <- rtrt
        }else{
          h2ars2[k,]<-rep(NA,15)
        }
        cat('done!\n')
      }
      
    }
    
    
    if(length(which(!is.na(h2ars2))) > 0){
      bucket[4,1]<-spatialModels[4]
      maxi<-which(h2ars2[,9] == max(h2ars2[,9], na.rm = TRUE))
      bucket[4,2:16]<-h2ars2[maxi,] #max(h2ars2)
      bucket[4,17:18]<-ars2[maxi,]
    }else{
      bucket[4,1]<-spatialModels[4]
      bucket[4,2:16]<-rep(NA,15)
      model_bucket[[4]]<-'error'
    }
    
  }
  
  if (debug){
    return(list(bucket,model_bucket))
  }else{
    return(bucket)
  }
  
  
}

################################### FOR SUPERRUNNER SPATIAL  ############################################
## spatial models
# 15% best genotypes say in % how many best were identified by the model

myNULLmodel_sp <- function(trait,data,E3, K3, soft="sommer",debug=FALSE, ge="ds",re="ds"){
  K3<<-K3
  E3<<-E3
  data[,"trtf"] <- as.character(data[,"trt"])
  data[,"rowf"] <- as.character(data[,"row"])
  data[,"colf"] <- as.character(data[,"col"])
  data[,"i_blockf"] <- as.character(data[,"i_block"])
  data[,"c_blockf"] <- as.character(data[,"c_block"])
  data[,"ci_blockf"] <- as.character(paste(data[,"c_block"],data[,"i_block"]))
  data[,"ci_block"] <- as.numeric(as.factor(data[,"ci_blockf"]))
  data[,"envf"] <- as.character(as.factor(data[,"env"]))
  ###Change - calculate the compound symetry structure for genetic effects
  data[,"EnvName"] <- paste(data$envf,data$trtf,sep = ":")
  
  levs_c_blockf <- table(data[,"c_blockf"])
  myfixedformula <- as.formula(paste(trait,"~ envf + envf:c_block"))

  myfixedformula_i <- as.formula(paste(trait,"~ c_block"))
  
#   coln <- c("c_block","i_block")
#   form <- NULL
#  for(icol in coln){
#    levs <- table(data[,icol])
#    if(length(levs) > 1){
#      toadd <- paste0(" vs(",ge,"(envf,diag(5)),",icol,")")
#      if(is.null(form)){form <- toadd}else{form <- paste(form,toadd, sep=" + ")}
#    }
#  }
  
  if(soft=="asreml"){
    
  }else{
    
    # if(ge == "cs"){cs <- TRUE; ge="us"}else{cs<-FALSE}
    # baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3)")
    # if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3)  +vs(ds(envf),rowf) + vs(ds(envf),colf) ")}  #if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3)")} ##
    # if(!is.null(form)){baserandom <- paste(baserandom,form, sep=" + ")}
  
    # if(ge == "us"){baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3) + vs(ds(envf),rowf) + vs(ds(envf),colf)")}
    if(ge == "ds"){
      baserandom <- paste0("~ vs(EnvName,Gu=kronecker(E3,K3, make.dimnames=TRUE))"," + vs(trtf,Gu=K3) + vs(ds(envf),rowf) + vs(ds(envf),colf)")
    } else {if(ge == "us"){baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3) + vs(ds(envf),rowf) + vs(ds(envf),colf)")}}
    
    myrandomformula <- as.formula(baserandom)
    myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))
    
    # Individual trial level####
    baserandom_i <- paste0("~ vs(trtf,Gu=K3) + rowf + colf")
    myrandomformula_i <- as.formula(baserandom_i)
    myresidualformula_i <- as.formula(paste0("~ units"))
    
    h2L_i <- numeric()
    env_i <- unique(data$envf)
    for(ii in env_i){
      mod.s_i <- try(
        mmer(fixed=myfixedformula_i,
             random=myrandomformula_i,
             rcov=myresidualformula_i,
             data=data[data$envf==ii,], verbose = FALSE), silent = TRUE)
      if(class(mod.s_i) != "try-error"){
        if(length(mod.s_i) > 0){
      n.g_i      <- length(unique(data$trtf))  # number of genotypes
      vc.g_i     <- as.numeric(mod.s_i$sigma[["u:trtf"]]) # genetic variance component
      C22.g_i    <- mod.s_i$PevU[["u:trtf"]][[1]]       #
      trC22.g_i  <- sum(diag(as.matrix(C22.g_i)))         # trace
      vdBLUP.g_i <- 2/n.g_i*(trC22.g_i-(sum(C22.g_i)-trC22.g_i)/(n.g_i-1)) # Mean variance of a difference of two genotypic BLUPs
      # H2 Cullis at individual trial level ###
      h2L_i[ii] <- 1-(vdBLUP.g_i / 2 / vc.g_i)
        }
      }
    }
    h2L1 <- h2L_i[1]
    h2L2 <- h2L_i[2]
    h2L3 <- h2L_i[3]
    h2L4 <- h2L_i[4]
    h2L5 <- h2L_i[5]
    h2L <- mean(h2L_i,na.rm=T)
    ###
    
    mod.s <- try(
      mmer(fixed=myfixedformula,
           random=myrandomformula,
           rcov=myresidualformula,
           data=data, verbose = FALSE), silent = TRUE)
    if(class(mod.s) != "try-error"){
      if(length(mod.s) > 0){
        # pp <- predict(mod.s)
        pp <- predict(mod.s, RtermsToForce = 1:2, FtermsToForce = 1:2)
      }else{pp<-list()}
    }else{pp<-list()}
    if (class(mod.s) != "try-error"  & length(pp) > 0) {
      # print(summary(mod.s)$varcomp)
      VCP <- summary(mod.s)$varcomp
      # vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
      # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      # vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      # if(length(vgmain) > 0){vg <- vg+vgmain}
      # ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
      # mve <- mean(ve)
      # h2 <-  mean(vg/(vg+ve))#average across locations
      # ######cambios###
      # mvg<-mean(vg)
      # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
      # vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
      # vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      # if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
      # mvgxe<-mean(vgxe)
      # 
      
      vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
      mvg <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      mvgxe<-mean(VCP[intersect(grep("u:EnvName",rownames(VCP)),vcpIndex),"VarComp"])
      
      # vg_gxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
      # mvg_ge <- mean(vg_gxe)
      mve <- mean(ve)
      
      #h2 <- 1 - (sum(diag(mod.s$PevU$`u:trtf`[[trait]]))/nrow(mod.s$PevU$`u:trtf`[[trait]])/(mvg))
      #h2 <-  mean(vg_gxe/(vg_gxe + ve/length(levs_c_blockf)))#average across locations
      #h2 <-  mean((mvg+mvgxe)/(mvg + mvgxe + ve/2)) #r=2 
      
      # Global H2 Cullis ###
      n.g      <- length(unique(mod.s$data$trtf))
      C22.g    <- mod.s$PevU[["u:trtf"]][[1]]
      trC22.g  <- sum(diag(as.matrix(C22.g)))         # trace
      vdBLUP.g <- 2/n.g*(trC22.g-(sum(C22.g)-trC22.g)/(n.g-1)) # Mean variance of a difference of two genotypic BLUPs
      h2 <- 1-(vdBLUP.g / 2 / mvg)        
      
      # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
      # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
      # mvg <- mean(vg)
      # mvgxe <- mvg_ge - mvg
      ratioGGE <- mvgxe/mvg
      
      nenv <- unique(pp$fitted$env)
      withinCor <- numeric()
      for(ienv in nenv){
        useEnv <- which(pp$fitted$env == ienv)
        withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
      }
      pa <- mean(withinCor)
      bic <- mod.s$BIC
      
      ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
      orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
      orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
      orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
      orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
      best <- numeric()
      #new change realized selection gain
      rsg <- numeric()
      for(ienv in nenv){
        a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
        b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
        best[ienv] <- length(intersect(a,b))/length(a)
        rsg[ienv] <- (mean(orderPred[which(orderPred$env == ienv)[1:ngeno],"ygorro"]) - mean(orderPred[which(orderPred$env == ienv),"ygorro"]))/mean(orderPred[which(orderPred$env == ienv),"ygorro"])*100
      }
      mbest <- mean(best)
      mrsg <- mean(rsg)
      outp<-c(h2,h2L1,h2L2,h2L3,h2L4,h2L5,h2L,mve,pa,bic,mbest,mvg,mvgxe,ratioGGE,mrsg)
      #outp<-c(h2,mve,pa,bic,mbest,mvg,mvgxe)
      
      ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
      
    }else{
      ifelse(debug==FALSE,return(rep(NA,15)),return(list(rep(NA,15),'model-error')))
    }
  }
}

# my2Dmodel_sp <- function(trait,data, K3, soft="sommer",debug=FALSE, ge="ds",re="ds"){
#   K3<<-K3
#   data[,"trtf"] <- as.character(data[,"trt"])
#   data[,"rowf"] <- as.character(data[,"row"])
#   data[,"colf"] <- as.character(data[,"col"])
#   data[,"i_blockf"] <- as.character(data[,"i_block"])
#   data[,"c_blockf"] <- as.character(data[,"c_block"])
#   data[,"ci_blockf"] <- as.character(paste(data[,"c_block"],data[,"i_block"]))
#   data[,"ci_block"] <- as.numeric(as.factor(data[,"ci_blockf"]))
#   data[,"envf"] <- as.character(as.factor(data[,"env"]))
#   
#   myfixedformula <- as.formula(paste(trait,"~ envf"))
#   
#   coln <- c("c_block","i_block")
#   form <- NULL
# #  for(icol in coln){
# #    levs <- table(data[,icol])
# #    if(length(levs) > 1){
# #      toadd <- paste0(" vs(",ge,"(envf,diag(5)),",icol,")")
# #      if(is.null(form)){form <- toadd}else{form <- paste(form,toadd, sep=" + ")}
# #    }
# #  }
#   
#   if(soft=="asreml"){
#     
#   }else{
#     if(ge == "cs"){cs <- TRUE; ge="us"}else{cs<-FALSE}
#     baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3) +  vs(ds(envf),rowf) + vs(ds(envf),colf) + vs(ds(envf),spl2D(row,col))")  #+ vs(ds(envf),rowf) + vs(ds(envf),colf) 
#     if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3)")}
#     
#     if(!is.null(form)){baserandom <- paste(baserandom,form, sep=" + ")}
#     myrandomformula <- as.formula(baserandom)
#     myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))
#     print(myrandomformula)
#     mod.s <- try(
#       mmer(fixed=myfixedformula,
#            random=myrandomformula,
#            rcov=myresidualformula,
#            data=data, verbose = FALSE), silent = TRUE)
#     if(class(mod.s) != "try-error"){
#       if(length(mod.s) > 0){
#         pp <- predict(mod.s)
#       }else{pp<-list()}
#     }else{pp<-list()}
#     if (class(mod.s) != "try-error"  & length(pp) > 0) {
#       VCP <- summary(mod.s)$varcomp
#       vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
#       vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       if(length(vgmain) > 0){vg <- vg+vgmain}
#       ve <- abs(VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"])
#       mve <- mean(ve)
#       h2 <-  mean(vg/(vg+ve))#average across locations
#       ######cambios###
#       mvg<-mean(vg)
#       vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
#       vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
#       vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
#       mvgxe<-mean(vgxe)
#       
#       nenv <- unique(pp$fitted$env)
#       withinCor <- numeric()
#       for(ienv in nenv){
#         useEnv <- which(pp$fitted$env == ienv)
#         withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
#       }
#       pa <- mean(withinCor)
#       bic <- mod.s$BIC
#       
#       ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
#       orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
#       orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
#       orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
#       orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
#       best <- numeric()
#       for(ienv in nenv){
#         a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
#         b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
#         best[ienv] <- length(intersect(a,b))/length(a)
#       }
#       mbest <- mean(best)
#       outp<-c(h2,mve,pa,bic,mbest,mvg,mvgxe)
#       # h2 <- 1 - (mean(mod.s$PevU$`u:trtf`[[trait]])/(vg))
#       ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
#     }else{
#       ifelse(debug==FALSE,return(rep(NA,7)),return(list(rep(NA,7),'model-error')))
#     }
#   }
# }
# 
# myAR1model_sp <- function(trait,data, rho,direction,K3,soft="sommer",debug=FALSE, ge="ds",re="ds"){
#   K3<<-K3
#   data[,"trtf"] <- as.character(data[,"trt"])
#   data[,"rowf"] <- as.character(data[,"row"])
#   data[,"colf"] <- as.character(data[,"col"])
#   data[,"i_blockf"] <- as.character(data[,"i_block"])
#   data[,"c_blockf"] <- as.character(data[,"c_block"])
#   data[,"ci_blockf"] <- as.character(paste(data[,"c_block"],data[,"i_block"]))
#   data[,"ci_block"] <- as.numeric(as.factor(data[,"ci_blockf"]))
#   data[,"envf"] <- as.character(as.factor(data[,"env"]))
#   
#   oposite <- setdiff(c("rowf","colf"),c(direction,paste0(direction,"f")))
#   touse <- intersect(c("rowf","colf"),c(direction,paste0(direction,"f")))
#   
#   myfixedformula <- as.formula(paste(trait,"~ envf"))
#   myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))
#   
#   
#   coln <- c("c_block","i_block")
#   form <- NULL
# #  for(icol in coln){
# #    levs <- table(data[,icol])
# #    if(length(levs) > 1){
# #      toadd <- paste0("vs(",ge,"(envf,diag(5)),",icol,")")
# #      if(is.null(form)){form <- toadd}else{form <- paste(form,toadd, sep=" + ")}
# #    }
# #  }
#   
#   if(soft=="asreml"){
#     
#   }else{
#     
#     if(ge == "cs"){cs <- TRUE; ge="us"}else{cs<-FALSE}
#     baserandom <- paste0("~  vs(",ge,"(envf), trtf,Gu=K3) + vs(ds(envf),",touse,", Gu=AR1(",touse,", rho=",rho,")) + vs(ds(envf),",oposite,")")
#     if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3)  ")}
#     
#     if(!is.null(form)){baserandom <- paste(baserandom,form, sep=" + ")}
#     myrandomformula <- as.formula(baserandom)
#     print(myrandomformula)
#     mod.s <- try(
#       mmer(fixed=myfixedformula,
#            random=myrandomformula,
#            rcov=myresidualformula,
#            data=data, verbose = FALSE), silent = TRUE)
#     if(class(mod.s) != "try-error"){
#       if(length(mod.s) > 0){
#         pp <- predict(mod.s)
#       }else{pp<-list()}
#     }else{pp<-list()}
#     if (class(mod.s) != "try-error"  & length(pp) > 0) {
#       VCP <- summary(mod.s)$varcomp
#       vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
#       vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       if(length(vgmain) > 0){vg <- vg+vgmain}
#       ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
#       mve <- mean(ve)
#       h2 <-  mean(vg/(vg+ve))#average across locations
#       ######cambios###
#       mvg<-mean(vg)
#       vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
#       vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
#       vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
#       if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
#       mvgxe<-mean(vgxe)
#       
#       nenv <- unique(pp$fitted$env)
#       withinCor <- numeric()
#       for(ienv in nenv){
#         useEnv <- which(pp$fitted$env == ienv)
#         withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
#       }
#       pa <- mean(withinCor)
#       bic <- mod.s$BIC
#       
#       
#       ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
#       orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
#       orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
#       orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
#       orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
#       best <- numeric()
#       for(ienv in nenv){
#         a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
#         b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
#         best[ienv] <- length(intersect(a,b))/length(a)
#       }
#       mbest <- mean(best)
#       outp<-c(h2,mve,pa,bic,mbest,mvg,mvgxe)
#       # h2 <- 1 - (mean(mod.s$PevU$`u:trtf`[[trait]])/(vg))
#       ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
#     }else{
#       ifelse(debug==FALSE,return(rep(NA,7)),return(list(rep(NA,7),'model-error')))
#     }
#   }
#   
# }

myAR1AR1model_sp <- function(trait,data, rho1, rho2,E3,K3,soft="sommer",debug=FALSE, ge="ds",re="ds"){
  K3<<-K3
  E3<<-E3
  data[,"trtf"] <- as.character(data[,"trt"])
  data[,"rowf"] <- as.character(data[,"row"])
  data[,"colf"] <- as.character(data[,"col"])
  data[,"i_blockf"] <- as.character(data[,"i_block"])
  data[,"c_blockf"] <- as.character(data[,"c_block"])
  data[,"ci_blockf"] <- as.character(paste(data[,"c_block"],data[,"i_block"]))
  data[,"ci_block"] <- as.numeric(as.factor(data[,"ci_blockf"]))
  data[,"envf"] <- as.character(as.factor(data[,"env"]))
  data[,"EnvName"] <- paste(data$envf,data$trtf,sep = ":")
  
  levs_c_blockf <- table(data[,"c_blockf"])
  myfixedformula <- as.formula(paste(trait,"~ envf + envf:c_block"))

  myfixedformula_i <- as.formula(paste(trait,"~ c_block"))
  
  coln <- c("c_block","i_block")
  form <- NULL
#  for(icol in coln){
#    levs <- table(data[,icol])
#    if(length(levs) > 1){
#      toadd <- paste0(" vs(",ge,"(envf,diag(5)),",icol,")")
#      if(is.null(form)){form <- toadd}else{form <- paste(form,toadd, sep=" + ")}
#    }
#  }
  
  if(soft=="asreml"){
    
  }else{
    
    # if(ge == "cs"){cs <- TRUE; ge="us"}else{cs<-FALSE}
    # baserandom <- paste0("~  vs(",ge,"(envf), trtf,Gu=K3) + vs(ds(envf),rowf, Gu=AR1(rowf, rho=",rho1,")) + vs(ds(envf),colf, Gu=AR1(colf, rho=",rho2,"))")
    # if(cs){baserandom <- paste0(baserandom," + vs(trtf,Gu=K3) ")}
    # if(!is.null(form)){baserandom <- paste(baserandom,form, sep=" + ")}
    
    # if(ge == "us"){baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3) + vs(ds(envf),rowf, Gu=AR1(rowf, rho=",rho1,")) + vs(ds(envf),colf, Gu=AR1(colf, rho=",rho2,"))")}
    if(ge == "ds"){
      baserandom <- paste0("~ vs(EnvName,Gu=kronecker(E3,K3, make.dimnames=TRUE)) + vs(trtf,Gu=K3) + vs(ds(envf),rowf,Gu=AR1(rowf, rho=",rho1,")) + vs(ds(envf),colf, Gu=AR1(colf, rho=",rho2,"))")
    } else {if(ge == "us"){baserandom <- paste0("~ vs(",ge,"(envf), trtf,Gu=K3)+ vs(ds(envf),rowf, Gu=AR1(rowf, rho=",rho1,")) + vs(ds(envf),colf, Gu=AR1(colf, rho=",rho2,"))")}}
    
    myrandomformula <- as.formula(baserandom)
    myresidualformula <- as.formula(paste0("~ vs(",re,"(envf), units)"))
    
    #Individual trial level####
    baserandom_i <- paste0("~ vs(trtf,Gu=K3) + vs(rowf,Gu=AR1(rowf, rho=",rho1,")) + vs(colf, Gu=AR1(colf, rho=",rho2,"))")
    myrandomformula_i <- as.formula(baserandom_i)
    myresidualformula_i <- as.formula(paste0("~ units"))

    h2L_i <- numeric()
    env_i <- unique(data$envf)
    for(ii in env_i){
      mod.s_i <- try(
        mmer(fixed=myfixedformula_i,
             random=myrandomformula_i,
             rcov=myresidualformula_i,
             data=data[data$envf==ii,], verbose = FALSE), silent = TRUE)
      if(class(mod.s_i) != "try-error"){
        if(length(mod.s_i) > 0){
      n.g_i      <- length(unique(data$trtf))  # number of genotypes
      vc.g_i     <- as.numeric(mod.s_i$sigma[["u:trtf"]]) # genetic variance component
      C22.g_i    <- mod.s_i$PevU[["u:trtf"]][[1]]       #
      trC22.g_i  <- sum(diag(as.matrix(C22.g_i)))         # trace
      vdBLUP.g_i <- 2/n.g_i*(trC22.g_i-(sum(C22.g_i)-trC22.g_i)/(n.g_i-1)) # Mean variance of a difference of two genotypic BLUPs
      # H2 Cullis at individual trial level ###
      h2L_i[ii] <- 1-(vdBLUP.g_i / 2 / vc.g_i)
        }
      }
    }
    h2L1 <- h2L_i[1]
    h2L2 <- h2L_i[2]
    h2L3 <- h2L_i[3]
    h2L4 <- h2L_i[4]
    h2L5 <- h2L_i[5]
    h2L <- mean(h2L_i, na.rm = T)
    ##

    mod.s <- try(
      mmer(fixed=myfixedformula,
           random=myrandomformula,
           rcov=myresidualformula,
           data=data, verbose = FALSE), silent = TRUE)
    if(class(mod.s) != "try-error"){
      if(length(mod.s) > 0){
        #pp <- predict(mod.s)
        pp <- predict(mod.s, RtermsToForce = 1:2, FtermsToForce = 1:2)
      }else{pp<-list()}
    }else{pp<-list()}
    if (class(mod.s) != "try-error"  & length(pp) > 0) {
      VCP <- summary(mod.s)$varcomp
      # vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
      # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      # vgmain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      # if(length(vgmain) > 0){vg <- vg+vgmain}
      # ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
      # mve <- mean(ve)
      # h2 <-  mean(vg/(vg+ve))#average across locations
      # ######cambios###
      # mvg<-mean(vg)
      # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
      # vgxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
      # vgxemain <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      # if(length(vgxemain) > 0){vgxe <- vgxe+vgxemain}
      # mvgxe<-mean(vgxe)
      ####

      vcpIndex <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==2)
      mvg <- VCP[intersect(grep("u:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      mvgxe<-VCP[intersect(grep("u:EnvName",rownames(VCP)),vcpIndex),"VarComp"]
      
      # vg_gxe <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex),"VarComp"]
      ve <- VCP[intersect(grep("units",rownames(VCP)),vcpIndex),"VarComp"]
      # mvg_ge <- mean(vg_gxe)
      mve <- mean(ve)
      
      #h2 <- 1 - (sum(diag(mod.s$PevU$`u:trtf`[[trait]]))/nrow(mod.s$PevU$`u:trtf`[[trait]])/(mvg))
      # h2 <-  mean(vg_gxe/(vg_gxe + ve/length(levs_c_blockf)))#average across locations
      #h2 <-  mean((mvg+mvgxe)/(mvg + mvgxe + ve/2)) #r=2 
      
      # Global H2 Cullis ###
      n.g      <- length(unique(mod.s$data$trtf))
      C22.g    <- mod.s$PevU[["u:trtf"]][[1]]
      trC22.g  <- sum(diag(as.matrix(C22.g)))         # trace
      vdBLUP.g <- 2/n.g*(trC22.g-(sum(C22.g)-trC22.g)/(n.g-1)) # Mean variance of a difference of two genotypic BLUPs
      h2 <- 1-(vdBLUP.g / 2 / mvg)        
      
      # vcpIndex2 <- which(unlist(lapply(rownames(VCP), function(x){length(strsplit(x,":")[[1]])}))==3)
      # vg <- VCP[intersect(grep("[0-9]:trtf",rownames(VCP)),vcpIndex2),"VarComp"]
      # mvg <- mean(vg)
      # mvgxe <- mvg_ge - mvg
      ratioGGE <- mvgxe/mvg
      
      nenv <- unique(pp$fitted$env)
      withinCor <- numeric()
      for(ienv in nenv){
        useEnv <- which(pp$fitted$env == ienv)
        withinCor[ienv]<- cor(pp$fitted$predicted.value.y[useEnv], pp$fitted$ygorro[useEnv])
      }
      pa <- mean(withinCor)
      bic <- mod.s$BIC
      
      ngeno <- round(length(unique(pp$fitted$trtf))*.15) # find best 10%
      orderGorro <- unique(pp$fitted[with(pp$fitted, order(envf,-ygorro)), ])
      orderGorro<-orderGorro[!duplicated(orderGorro[,c("trt_names","env")]),]
      orderPred <- pp$fitted[with(pp$fitted, order(envf,-predicted.value.y)), ]
      orderPred<-orderPred[!duplicated(orderPred[,c("trt_names","env")]),]
      best <- numeric()
      #new change realized selection gain
      rsg <- numeric()
      for(ienv in nenv){
        a <- orderGorro[which(orderGorro$env == ienv)[1:ngeno],"trtf"]
        b <- orderPred[which(orderPred$env == ienv)[1:ngeno],"trtf"]
        best[ienv] <- length(intersect(a,b))/length(a)
        rsg[ienv] <- (mean(orderPred[which(orderPred$env == ienv)[1:ngeno],"ygorro"]) - mean(orderPred[which(orderPred$env == ienv),"ygorro"]))/mean(orderPred[which(orderPred$env == ienv),"ygorro"])*100
      }
      
      mbest <- mean(best)
      mrsg <- mean(rsg)
      outp<-c(h2,h2L1,h2L2,h2L3,h2L4,h2L5,h2L,mve,pa,bic,mbest,mvg,mvgxe,ratioGGE,mrsg)
      #outp<-c(h2,mve,pa,bic,mbest,mvg,mvgxe)
      # h2 <- 1 - (mean(mod.s$PevU$`u:trtf`[[trait]])/(vg))
      ifelse(debug==FALSE,return(outp),return(list(outp,mod.s)))
      
    }else{
      ifelse(debug==FALSE,return(rep(NA,15)),return(list(rep(NA,15),'model-error')))
    }
  }
}


## core runner
coreRunner_sp<-function(fieldbookInput=NULL,E2=NULL,K2=NULL,debug=FALSE,ge="ds",re="ds",step_ar1a=0.1,step_ar1ar1a=0.1,ncore=0){   ## input must be a list with 4 experimental designs
  
  K2<<-K2
  spatialModels<-c('none','spl2D','AR1','AR1xAR1')
  
  bucket<-data.frame(spatial=numeric(length(4)),
                     h2G=numeric(length(4)),
                     h2L1=numeric(length(4)),
                     h2L2=numeric(length(4)),
                     h2L3=numeric(length(4)),
                     h2L4=numeric(length(4)),
                     h2L5=numeric(length(4)),
                     h2L=numeric(length(4)),
                     ve=numeric(length(4)),
                     pa=numeric(length(4)),
                     bic=numeric(length(4)),
                     best=numeric(length(4)),
                     vg=numeric(length(4)),
                     vgxe=numeric(length(4)),
                     ratioGE=numeric(length(4)),
                     rsg=numeric(length(4)),
                     rho1=numeric(length(4))*NA,
                     rho2=numeric(length(4))*NA
                     
  )
  
  
  model_bucket<-list()
  
  
  ## no adjustment
  
  cat('\n**running with no spatial adjustment...')
  bucket[1,1]<-spatialModels[1]
  if (debug){
    invisible(capture.output(tmpm<-try(myNULLmodel_sp(trait="y", data=fieldbookInput,E3=E2, K3=K2,debug=TRUE,ge=ge,re=re))))
    #tmpm<-try(myNULLmodel(trait="y", data=fieldbookInput, K=K,debug=TRUE,ge=ge,re=re))
    
    if(class(tmpm) != "try-error"){
      bucket[1,2:16]<-tmpm[[1]]
      model_bucket[[1]]<-tmpm[[2]]}else{
        bucket[1,2:16]<-rep(NA,15)
      }
  }else{
    invisible(capture.output(rtrt<-try(myNULLmodel_sp(trait="y", data=fieldbookInput,E3=E2, K3=K2,debug=FALSE,ge=ge,re=re))))
    if(class(rtrt) != "try-error"){
      bucket[1,2:16]<-rtrt
    }else{
      bucket[1,2:16]<-rep(NA,15)
    }
  }
  cat(' done!')
  
  # ## spl2D models
  # cat('\n**running spl2D correction...')
  # bucket[2,1]<-spatialModels[2]
  # if (debug){
  #   invisible(capture.output(tmpm<-try(my2Dmodel_sp(trait="y", data=fieldbookInput, K3=K2,debug = TRUE,ge=ge,re=re))))
  #   #tmpm<-try(my2Dmodel(trait="y", data=fieldbookInput, K3=K2,debug = TRUE,ge=ge,re=re))
  #   if(class(tmpm) != "try-error"){
  #     bucket[2,2:8]<-tmpm[[1]]
  #     model_bucket[[2]]<-tmpm[[2]]}else{
  #       bucket[2,2:8]<-rep(NA,7)
  #     }
  #   
  # }else{
  #   invisible(capture.output(rtrt<-my2Dmodel_sp(trait="y", data=fieldbookInput, K3=K2,debug = FALSE,ge=ge,re=re)))
  #   if(class(rtrt) != "try-error"){
  #     
  #     bucket[2,2:8]<-rtrt
  #   }else{
  #     bucket[2,2:8]<-rep(NA,7)
  #     
  #   }
  #   
  # }
  # cat(' done!')
  # 
  # 
  # 
  # ## AR1 models
  # ars <- seq(-1,1,step_ar1a)
  # cat('\n**running AR1 adjustment...\n')
  # h2ars <- matrix(NA,nrow=length(ars),ncol=7);colnames(h2ars) <- c("h2","ve","pa","bic","best","vg","vgxe"); rownames(h2ars) <- ars
  # 
  # if (debug){
  #   arsmodel<-list()
  #   
  #   if (ncore>0){
  #     n.cores <- ncore 
  #     cl <- snow::makeCluster(n.cores)
  #     registerDoSNOW(cl)
  #     
  #     arsmodel<-foreach(k= 1:length(ars),.packages = 'sommer',.export = 'myAR1model_sp') %dopar% myAR1model_sp(trait="y", data=fieldbookInput,rho=ars[k], direction = "col", K3=K2,debug = TRUE,ge=ge,re=re)
  #     snow::stopCluster(cl)
  #     
  #     for (k in 1:length(arsmodel)){
  #       if(class(arsmodel[[k]]) != "try-error"){
  #         h2ars[k,]<- arsmodel[[k]][[1]]
  #       }else{
  #         h2ars[k,]<- rep(NA,7)
  #         
  #       }
  #     }
  #     cat('done!\n')
  #     
  #   }else{
  #     
  #     for(k in 1:length(ars)){
  #       cat(paste0('      rho=',ars[k],'... '))
  #       invisible(capture.output(arsmodel[[k]]<-try(myAR1model_sp(trait="y", data=fieldbookInput,rho=ars[k], direction = "col", K3=K2,debug = TRUE,ge=ge,re=re))))
  #       if(class(arsmodel[[k]]) != "try-error"){
  #         h2ars[k,]<- arsmodel[[k]][[1]]
  #       }else{
  #         h2ars[k,]<- rep(NA,7)
  #         
  #       }
  #       cat('done!\n')
  #     }
  #     
  #   }
  #   if(length(which(!is.na(h2ars))) > 0){
  #     bucket[3,1]<-spatialModels[3]
  #     maxi<-which(h2ars[,3] == max(h2ars[,3], na.rm = TRUE))
  #     bucket[3,2:8]<-h2ars[maxi[1],]
  #     bucket[3,9]<-ars[maxi]
  #     model_bucket[[3]]<-arsmodel[[maxi[1]]][[2]]
  #   }else{
  #     bucket[3,1]<-spatialModels[3]
  #     bucket[3,2:8]<-rep(NA,7)
  #     model_bucket[[3]]<-'error'
  #   }
  #   
  # }else{
  #   
  #   
  #   
  #   if (ncore>0){
  #     n.cores <- ncore 
  #     cl <- snow::makeCluster(n.cores)
  #     registerDoSNOW(cl)
  #     
  #     rtrt<-foreach(k= 1:length(ars),.packages = 'sommer',.export = 'myAR1model_sp') %dopar% myAR1model_sp(trait="y", data=fieldbookInput,rho=ars[k], direction = "col", K3=K2,debug = FALSE,ge=ge,re=re)
  #     snow::stopCluster(cl)
  #     
  #     for(k in 1:length(rtrt)){
  #       if(class(rtrt[[k]]) != "try-error"){
  #         h2ars[k,] <- rtrt[[k]]
  #       }else{
  #         h2ars[k,]<-rep(NA,7)
  #       }
  #       cat('done!\n')
  #     }
  #     
  #   }else{
  #     
  #     
  #     for(k in 1:length(ars)){
  #       cat(paste0('      rho=',ars[k],'... '))
  #       invisible(capture.output(rtrt <- try(myAR1model_sp(trait="y", data=fieldbookInput,rho=ars[k], direction = "col", K3=K2,debug=FALSE, ge=ge,re=re))))
  #       if(class(rtrt) != "try-error"){
  #         h2ars[k,] <- rtrt
  #       }else{
  #         h2ars[k,]<-rep(NA,7)
  #       }
  #       cat('done!\n')
  #     }
  #     
  #   }
  #   
  #   
  #   # store the one that gave the maximum PA
  #   if(length(which(!is.na(h2ars))) > 0){
  #     bucket[3,1]<-spatialModels[3]
  #     maxi<-which(h2ars[,3] == max(h2ars[,3], na.rm = TRUE))
  #     bucket[3,2:8]<-h2ars[maxi,]
  #     bucket[3,9]<-ars[maxi]
  #   }else{
  #     bucket[3,1]<-spatialModels[3]
  #     bucket[3,2:8]<-rep(NA,7)
  #     model_bucket[[3]]<-'error'
  #   }
  #   
  # }
  
  ## AR1xAR1 models
  ars <- seq(-1,1,step_ar1ar1a)
  ars2 <- expand.grid(ars,ars)
  cat('**running AR1AR1 adjustment...\n')
  h2ars2 <- matrix(NA,nrow=nrow(ars2),ncol=15);colnames(h2ars2) <- c("h2G","h2L1","h2L2","h2L3","h2L4","h2L5","h2L","ve","pa","bic","best","vg","vgxe", "ratioGE", "rsg"); rownames(h2ars2) <- apply(ars2,1,function(x){paste(x,collapse = "-")})
  
  if (debug){
    
    
    if (ncore>0){
      n.cores <- ncore 
      cl <- snow::makeCluster(n.cores)
      registerDoSNOW(cl)
      
      arsmodel<-foreach(k= 1:nrow(ars2),.packages = 'sommer',.export = 'myAR1AR1model_sp') %dopar% myAR1AR1model_sp(trait="y", data=fieldbookInput,rho1=ars2[k,1],rho2 = ars2[k,2],E3=E2,K3=K2,debug=TRUE, ge=ge,re=re)
      snow::stopCluster(cl)
      
      for (k in 1:length(arsmodel)){
        if(class(arsmodel[[k]]) != "try-error"){
          h2ars2[k,]<- arsmodel[[k]][[1]]
        }else{
          h2ars2[k,]<- rep(NA,15)
          
        }
      }
      cat('done!\n')
      
    }else{
      
      arsmodel<-list()
      for(k in 1:nrow(ars2)){
        cat(paste0('      rho1=',format(round(ars2[k,1],2),nsmall = 2),', rho2=',format(round(ars2[k,2],2),nsmall = 2),'... '))
        invisible(capture.output(arsmodel[[k]]<-try(myAR1AR1model_sp(trait="y", data=fieldbookInput,rho1=ars2[k,1],rho2 = ars2[k,2],E3=E2,K3=K2,debug=TRUE, ge=ge,re=re))))
        if(class(arsmodel[[k]]) != "try-error"){
          h2ars2[k,]<- arsmodel[[k]][[1]]
        }else{
          h2ars2[k,]<- rep(NA,15)
          
        }
        cat('done!\n')
      }
      
    }
    if(length(which(!is.na(h2ars2))) > 0){
      bucket[4,1]<-spatialModels[4]
      maxi<-which(h2ars2[,9] == max(h2ars2[,9], na.rm = TRUE))[1]
      #maxi<-which(h2ars2[,3] == max(h2ars2[,3], na.rm = TRUE))[1]
      bucket[4,2:16]<-h2ars2[maxi,] #max(h2ars2)
      bucket[4,17:18]<-ars2[maxi,]
      model_bucket[[4]]<-arsmodel[[maxi]][[2]]
    }else{
      bucket[4,1]<-spatialModels[4]
      bucket[4,2:16]<-rep(NA,15)
      model_bucket[[4]]<-'error'
    }
  }else{
    
    
    if (ncore>0){
      n.cores <- ncore 
      cl <- snow::makeCluster(n.cores)
      registerDoSNOW(cl)
      
      rtrt<-foreach(k= 1:length(ars),.packages = 'sommer',.export = 'myAR1AR1model_sp') %dopar% myAR1AR1model_sp(trait="y", data=fieldbookInput,rho1=ars2[k,1],rho2 = ars2[k,2],E3=E2,K3=K2,debug = FALSE,ge=ge,re=re)
      snow::stopCluster(cl)
      
      for(k in 1:length(rtrt)){
        if(class(rtrt[[k]]) != "try-error"){
          h2ars2[k,] <- rtrt[[k]]
        }else{
          h2ars2[k,]<-rep(NA,15)
        }
        cat('done!\n')
      }
      
    }else{
      
      for(k in 1:nrow(ars2)){
        cat(paste0('      rho1=',format(round(ars2[k,1],2),nsmall = 2),', rho2=',format(round(ars2[k,2],2),nsmall = 2),'... '))
        invisible(capture.output(rtrt <- try(myAR1AR1model_sp(trait="y", data=fieldbookInput,rho1=ars2[k,1],rho2 = ars2[k,2],E3=E2,K3=K2,debug = FALSE,ge=ge,re=re))))
        if(class(rtrt) != "try-error"){
          h2ars2[k,] <- rtrt
        }else{
          h2ars2[k,]<-rep(NA,15)
        }
        cat('done!\n')
      }
      
    }
    
    
    if(length(which(!is.na(h2ars2))) > 0){
      bucket[4,1]<-spatialModels[4]
      maxi<-which(h2ars2[,9] == max(h2ars2[,9], na.rm = TRUE))
      bucket[4,2:16]<-h2ars2[maxi,] #max(h2ars2)
      bucket[4,17:18]<-ars2[maxi,]
    }else{
      bucket[4,1]<-spatialModels[4]
      bucket[4,2:16]<-rep(NA,15)
      model_bucket[[4]]<-'error'
    }
    
  }
  
  if (debug){
    return(list(bucket,model_bucket))
  }else{
    return(bucket)
  }
  
  
}



############################# EXTRA FUNCTIONS FOR AGRICOLAE #########################################

design.alpha <-function (trt, k, r, serie = 2, seed = 0, kinds = "Super-Duper",randomization = TRUE){
  number <- 10
  if (serie > 0) 
    number <- 10^serie
  name.trt <- c(paste(deparse(substitute(trt))))
  ntr <- length(trt)
  if (seed == 0) {
    genera <- runif(1)
    seed <- .Random.seed[3]
  }
  set.seed(seed, kinds)
  s <- ntr/k
  if (ntr%%k != 0) {
   # cat("\nThe size of the block is not appropriate", "\nthe number of treatments must be multiple of k (size block) \n")
  }else {
    serie <- ""
    if (r == 2 & k <= s) {
      alpha <- matrix(0, nrow = k, ncol = r)
      alpha[2, 2] <- 1
      for (i in 3:k) {
        alpha[i, 2] <- alpha[i - 1, 2] + 1
      }
      serie <- "I"
    }
    if (r == 3 & s%%2 != 0 & k <= s) {
      alpha <- matrix(0, nrow = k, ncol = r)
      alpha[2, 2] <- 1
      alpha[2, 3] <- s - 1
      for (i in 3:k) {
        alpha[i, 2] <- alpha[i - 1, 2] + 1
        alpha[i, 3] <- alpha[i - 1, 3] - 1
      }
      serie <- "II"
    }
    if (r == 3 & s%%2 == 0 & k < s) {
      s1 <- s/2
      alpha <- matrix(0, nrow = k, ncol = r)
      alpha[2, 2] <- 1
      alpha[2, 3] <- s1
      for (i in 3:k) {
        alpha[i, 2] <- alpha[i - 1, 2] + 1
        alpha[i, 3] <- alpha[i - 2, 3] + 1
      }
      serie <- "III"
    }
    if (r == 4 & s%%2 != 0 & s%%3 != 0 & k <= s) {
      s2 <- (s + 1)/2
      alpha <- matrix(0, nrow = k, ncol = r)
      alpha[2, 2] <- 1
      alpha[2, 3] <- s - 1
      alpha[2, 4] <- s2
      for (i in 3:k) {
        alpha[i, 2] <- alpha[i - 1, 2] + 1
        alpha[i, 3] <- alpha[i - 1, 3] - 1
        alpha[i, 4] <- alpha[i - 2, 4] + 1
      }
      serie <- "IV"
    }
    if (serie == "") {
 #     cat("\nhelp(design.alpha): to see the series of alpha generators\n")
      stop
    }
    else {
      nf <- nrow(alpha)
      nc <- ncol(alpha)
      cc <- rep(alpha[, 1], s)
      for (i in 2:r) {
        cc <- c(cc, rep(alpha[, i], s))
      }
      dim(cc) <- c(nf, s, r)
      for (m in 1:r) cc[, 1, m] <- alpha[, m]
      for (i in 2:s) {
        for (j in 1:nf) {
          for (m in 1:r) {
            cc[j, i, m] <- cc[j, i - 1, m] + 1
            if (cc[j, i, m] >= s) 
              cc[j, i, m] <- 0
          }
        }
      }
      for (j in 1:nf) {
        cc[j, , ] <- cc[j, , ] + (j - 1) * s
      }
      intermediate <- cc
    #  cat("\nAlpha Design (0,1) - Serie ", serie, "\n")
      E <- (ntr - 1) * (r - 1)/((ntr - 1) * (r - 1) + r * 
                                  (s - 1))
    #  cat("\nParameters Alpha Design\n=======================")
    #  cat("\nTreatmeans :", ntr)
    #  cat("\nBlock size :", k)
    #  cat("\nBlocks     :", s)
    #  cat("\nReplication:", r, "\n")
    #  cat("\nEfficiency factor\n(E )", E, "\n\n<<< Book >>>\n")
      parameters <- list(design = "alpha", trt = trt, k = k, 
                         r = r, serie = serie, seed = seed, kinds = kinds)
      statistics <- data.frame(treatments = ntr, blocks = s, 
                               Efficiency = E)
      rownames(statistics) <- "values"
      for (m in 1:r) {
        for (j in 1:s) {
          aleatorio <- 1:k
          if (randomization) 
            aleatorio <- sample(1:k, k)
          cc[, j, m] <- cc[aleatorio, j, m]
        }
      }
      for (m in 1:r) {
        aleatorio <- 1:s
        if (randomization) 
          aleatorio <- sample(1:s, s)
        cc[, , m] <- cc[, aleatorio, m]
      }
      cc <- cc + 1
      block <- gl(s, k)
      md <- as.numeric(cc[, , 1])
      bp <- 1:ntr
      if (randomization) 
        bp <- sample(1:ntr, ntr)
      trt <- trt[bp]
      mtr <- trt[md]
      book <- data.frame(block = as.factor(block), trt = as.factor(mtr), 
                         replication = 1)
      for (i in 2:r) {
        md <- as.numeric(cc[, , i])
        mtr <- trt[md]
        book1 <- data.frame(block = as.factor(block), 
                            trt = as.factor(mtr), replication = i)
        book <- rbind(book, book1)
      }
      Rep <- book$replication
      plots <- Rep * number + (1:ntr)
      cols <- as.numeric(rep(gl(k, 1), s * r))
      book <- data.frame(plots = plots, cols = cols, book)
      book <- data.frame(row.names = NULL, book)
      book$block <- gl(s * r, k)
      book[, 2] <- as.factor(book[, 2])
      book[, 5] <- as.factor(book[, 5])
      names(book)[4] <- name.trt
      tr <- as.character(book[, 4])
      dim(tr) <- c(k, s, r)
      if (r == 2) 
        design <- list(rep1 = t(tr[, , 1]), rep2 = t(tr[, 
                                                        , 2]))
      if (r == 3) 
        design <- list(rep1 = t(tr[, , 1]), rep2 = t(tr[, 
                                                        , 2]), rep3 = t(tr[, , 3]))
      if (r == 4) 
        design <- list(rep1 = t(tr[, , 1]), rep2 = t(tr[, 
                                                        , 2]), rep3 = t(tr[, , 3]), rep4 = t(tr[, , 
                                                                                                4]))
      outdesign <- list(parameters = parameters, statistics = statistics, 
                        sketch = design, book = book)
      return(outdesign)
    }
  }
} 



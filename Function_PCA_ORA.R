PCA_calc<-function(EXPRESSIONTABLE=NULL,
                   MAPPINGFILE=NULL,
                   NAME_TABELLE=NULL,
                   ALL_GENES=NULL,
                   EXPRESSIONTABLE_MEAN=NULL,
                   TREATMENT_TABLE=NULL,
                   GROUPED_TABLE=NULL,
                   PLOT_COLOURS=NULL,
                   PCA_SCALING=NULL,
                   PCA_SCORE_TABLE=NULL,
                   ORA_FILTERING=NULL,
                   TOP=NULL,
                   TPM=NULL){
  
  ##################################### PACKAGES ###########################################################################
  
  library(dplyr)
  library(stringr)
  library(RColorBrewer)
  library(ggplot2)
  
  
  ################################################ INFO FILE/ INFO PDF #####################################################
  ## Create an INFO file as output for PCA summary etc.
  
  Name_Info<-paste0("Info_",NAME_TABELLE,".txt") 
  
  # add Text/ header to info file
  cat("Info file: Output PCA/ORA calculations and plot", file = Name_Info)
  # add 2 newlines
  cat("\n\n\n", file = Name_Info, append = TRUE)
  
  num.plots<-5
  
  Plots_output<-vector(num.plots, mode='list')
  
###################################################################################################################################
################################################## CHECK INPUT ####################################################################
###################################################################################################################################  

################################################## EXPRESSIONTABLE ################################################################
  
  if (EXPRESSIONTABLE[1,1] != rownames(EXPRESSIONTABLE)[1]){
    print("EXPRESSIONTABLE is valid.")
    
  }else if (EXPRESSIONTABLE[1,1] != rownames(EXPRESSIONTABLE)[1] && is.element(EXPRESSIONTABLE[1,1],ALL_GENES)==T){
    rownames(EXPRESSIONTABLE)<-EXPRESSIONTABLE[,1]
    print("First column set as rownames (Gene IDs).")
  
  }
  #else if (EXPRESSIONTABLE[1,1] != rownames(EXPRESSIONTABLE)[1] && is.element(EXPRESSIONTABLE[1,1],ALL_GENES)==F){

  #  EXPRESSIONTABLE1<-cbind.data.frame(rownames(EXPRESSIONTABLE),EXPRESSIONTABLE,stringsAsFactors=FALSE)
   # rownames(EXPRESSIONTABLE1)<-EXPRESSIONTABLE1[,1]
    #EXPRESSIONTABLE1<-data.matrix(EXPRESSIONTABLE1)
  #  print("Rownames set as First column (Gene IDs).")

  #    }
  #else{
  #  stop("First column of EXPRESSIONTABLE and rownames of EXPRESSIONTABLE don't match ALL_GENES. Class character or GeneIDs are needed.")
 # }
  
################################################## MAPPINGFILE ####################################################################
  
  if (is.element(rownames(EXPRESSIONTABLE)[1],MAPPINGFILE[,3])==FALSE | is.element(rownames(EXPRESSIONTABLE)[1],MAPPINGFILE[,3])==FALSE){
    
    stop("EXPRESSIONTABLE GeneIDs are not equal to MAPPINGFILE Identifiern/GeneIDs. Please check! or 
         EXPRESSIONTABLE_MEAN GeneIDs are not equal to MAPPINGFILE Identifiern/GeneIDs. Please check!")
  
  }else if (is.element(rownames(EXPRESSIONTABLE)[1],MAPPINGFILE[,3])==FALSE && is.element(rownames(EXPRESSIONTABLE)[1],MAPPINGFILE[,3])==FALSE){
    
    stop("EXPRESSIONTABLE GeneIDs are not equal to MAPPINGFILE Identifiern/GeneIDs. Please check! and
         EXPRESSIONTABLE_MEAN GeneIDs are not equal to MAPPINGFILE Identifiern/GeneIDs. Please check!")
  
  }else{
    print("MAPPINGFILE is valid.")
  }
  
################################################## ALL_GENES ######################################################################
   
 if (is.null(ALL_GENES)==TRUE && is.character(EXPRESSIONTABLE[,1])==TRUE){
    ALL_GENES<-as.data.frame(EXPRESSIONTABLE[,1])
    print("All GeneIDs of EXPRESSIONTABLE used as ALL_GENES.")
    
  }else if (is.null(ALL_GENES)==TRUE && is.character(EXPRESSIONTABLE[,1])==FALSE){
    stop("First column of EXPRESSIONTABLE is numeric. Class character is needed (Gene IDs).")
    
  }else if (is.null(ALL_GENES)==FALSE && is.character(EXPRESSIONTABLE[,1])==TRUE){
    print("ALL_GENES was set manually.")
  }
  
################################################## EXPRESSIONTABLE_MEAN ############################################################ 
 
  if(is.character(EXPRESSIONTABLE_MEAN[,1])==TRUE && levels(as.factor(EXPRESSIONTABLE_MEAN[,1] != rownames(EXPRESSIONTABLE_MEAN)))==FALSE){ ## hier noch ein kleiner Fehler
    rownames(EXPRESSIONTABLE_MEAN)<-EXPRESSIONTABLE_MEAN[,1]
    print("First column set as rownames for EXPRESSIONTABLE_MEAN")
    
  }else if (is.character(EXPRESSIONTABLE_MEAN[,1])==FALSE){
    stop("First column of EXPRESSIONTABLE_MEAN is numeric. Class character is needed (Gene IDs)")
  }else{
    print("EXPRESSIONTABLE_MEAN is valid.")
  }
 
################################################## GROUPED_TABLE ###################################################################   
  
  
  if (is.null(GROUPED_TABLE)==TRUE && (ncol(TREATMENT_TABLE))==2){
    GROUPED_TABLE<-factor(GROUPED_TABLE[,1])
    print("GROUPED_TABLE was not specified. First column of TREATMENT_TABLE was used for Treatment.")
    
  }else if(is.null(GROUPED_TABLE)==TRUE && (ncol(TREATMENT_TABLE))==3){
    GROUPED_TABLE<- factor(paste(TREATMENT_TABLE[,2],TREATMENT_TABLE[,3],sep=".")) 
    print("GROUPED_TABLE was not specified. Second and third column of TREATMENT_TABLE were grouped.")
    
  }else if(is.null(GROUPED_TABLE)==TRUE && (ncol(TREATMENT_TABLE))==4){
    GROUPED_TABLE<- factor(paste(TREATMENT_TABLE[,2],TREATMENT_TABLE[,3],TREATMENT_TABLE[,4],sep=".")) 
    print("GROUPED_TABLE was not specified. Second, third and fourth column of TREATMENT_TABLE were grouped.")
    
  }else if(is.null(GROUPED_TABLE)==TRUE && (ncol(TREATMENT_TABLE))==5 && is.element(TREATMENT_TABLE[,1],rownames(TREATMENT_TABLE))==T){
    GROUPED_TABLE<- factor(paste(TREATMENT_TABLE[,2],TREATMENT_TABLE[,3],TREATMENT_TABLE[,4],TREATMENT_TABLE[,5],sep=".")) 
    print("GROUPED_TABLE was not specified. Second, third, fourth and fifth column of TREATMENT_TABLE were grouped.")
  }else if(is.null(GROUPED_TABLE)==TRUE && (ncol(TREATMENT_TABLE))>4  && is.element(TREATMENT_TABLE[,1],rownames(TREATMENT_TABLE))==F){
    stop("There are to many columns in GROUPED_TABLE. Please check the table. Should consist of GeneID and max. 3 variable/parameter columns") 
    
  }else if (is.null(GROUPED_TABLE)==FALSE){
    print("GROUPED_TABLE was set manually")
  }
################################################## TPM_TABLE ###################################################################   
  
## This table will be needed to create Boxplots for the Top genes along the PC axes
## Need to look like this: rownames= Gene IDs; Colnames: Sample names, columns should be numeric TPM values
  if (is.null(TPM)==F && is.numeric(TPM[,1])==T ){  
    print("TPM was set correctly. First column of TREATMENT_TABLE was used for Treatment.")
 } else {
    stop('No table for TPM was found!')
  }
  
  
############################################ CHECK FOR PLOT_COLOURS TABLE ##########################################################


if (is.null(PLOT_COLOURS)==TRUE ){
  print("PLOT_COLOURS will be specified automatically.")
  
    ## Create the PLOT_COLOURS Table
    symb<-c("21","22","23","24","25")
    
    #Farben max 12
    cols<-brewer.pal(12, "Paired")
    
    #Bordercolour max 3
    
    oco<-c("black","gray77","grey100")
    
    # vetcors for loop
    co<-vector()
    sy<-vector()
    oc<-vector()
   
    if (ncol(TREATMENT_TABLE)==2){ 
      a<-TREATMENT_TABLE[,2]       
      af<-as.factor(a)             
      al<-levels(af)               
      la<-length(al)               
      
      
      if (la<=12){                
        for (t in 1:a){
          for (c in 1:la){
            assign(get(al)[c],cols[c])
          }
          
          co<-rbind(co,get(get(sorted[1,3])[t])) 
          sy<-rbind(sy,NA)                       
          oc<-rbind(oc,NA)                       
        }
      }else{
        stop("Too many variables are present. Please check your TREATMENT_TABLE. Max. 10 and/or 5 variables (3 columns)!")
      } 
      PLOT_COLOURS<-cbind(co,sy,oc) # final colours for plotting
      colnames(PLOT_COLOURS)<-c("Symbolcolour","Symbols","Bordercolour")
      TREATMENT_TABLE<-cbind(TREATMENT_TABLE,PLOT_COLOURS)
    }
    else if (ncol(TREATMENT_TABLE)==3){ # if table columns = 3
      a<-TREATMENT_TABLE[,2]
      b<-TREATMENT_TABLE[,3]
      
      af<-as.factor(a)
      bf<-as.factor(b)
      
      al<-levels(af)
      bl<-levels(bf)
      
      la<-length(al)
      lb<-length(bl) 
      
      n<-as.data.frame(cbind(c(la,lb),c("la","lb"),c("a","b"),c("al","bl")))
      sorted<-n[order(n[,1],decreasing = T),]
      
      if (as.numeric(sorted[1,1])<=12 && as.numeric(sorted[2,1])<=5){ 
        for (t in 1:length(get(as.vector(sorted[1,3])))){ 
          for (c in 1:get(as.character(sorted[1,2]))){
            assign(get(as.character(sorted[1,4]))[c],cols[c]) 
          }
          for (s in 1:get(as.character(sorted[2,2]))){
            assign(get(as.character(sorted[2,4]))[s],symb[s]) 
          }
          
          co<-rbind(co,get(as.character(get(as.vector(sorted[1,3]))[t])))
          sy<-rbind(sy,get(as.character(get(as.vector(sorted[2,3]))[t])))
          oc<-rbind(oc,NA)
        }
      }else{
        stop("Too many variables are present. Please check your TREATMENT_TABLE. Max. 10 and/or 5 variables! (3 columns)")
      } 
      
      PLOT_COLOURS<-cbind(co,sy,oc) # final colours for plotting
      colnames(PLOT_COLOURS)<-c("Symbolcolour","Symbols","Bordercolour")
      TREATMENT_TABLE<-cbind(TREATMENT_TABLE,PLOT_COLOURS)
      
    }
    else if (ncol(TREATMENT_TABLE)==4){ # if table columns = 4
      a<-TREATMENT_TABLE[,2]
      b<-TREATMENT_TABLE[,3]
      c<-TREATMENT_TABLE[,4]
      
      af<-as.factor(a)
      bf<-as.factor(b)
      cf<-as.factor(c)
      
      al<-levels(af)
      bl<-levels(bf)
      cl<-levels(cf)
      
      la<-length(al)
      lb<-length(bl)
      lc<-length(cl)
      
      n<-as.data.frame(cbind(c(la,lb,lc),c("la","lb","lc"),c("a","b","c"),c("al","bl","cl")))
      sorted<-n[order(n[,1],decreasing = T),]
      
      
      if (as.numeric(sorted[1,1])<=12 && as.numeric(sorted[2,1])<=5 && as.numeric(sorted[3,1])<=3){
        for (t in 1:length(get(as.vector(sorted[1,3])))){ 
          for (c in 1:get(as.character(sorted[1,2]))){
            assign(get(as.character(sorted[1,4]))[c],cols[c]) # colours
          }
          for (s in 1:get(as.character(sorted[2,2]))){
            assign(get(as.character(sorted[2,4]))[s],symb[s]) # symbols
          }
          for (o in 1:get(as.character(sorted[3,2]))){
            assign(get(as.character(sorted[3,4]))[o],oco[o]) # border colours
          }
          
          co<-rbind(co,get(as.character(get(as.vector(sorted[1,3]))[t])))
          sy<-rbind(sy,get(as.character(get(as.vector(sorted[2,3]))[t])))
          oc<-rbind(oc,get(as.character(get(as.vector(sorted[3,3]))[t])))
          }
        }else{
          stop("Too many variables are present. Please check your TREATMENT_TABLE. Max. 10 and/or 5, 3 variables (3 columns)!")
        } 
      
      PLOT_COLOURS<-cbind(co,sy,oc) # final colours for plotting
      colnames(PLOT_COLOURS)<-c("Symbolcolour","Symbols","Bordercolour")
      TREATMENT_TABLE<-cbind(TREATMENT_TABLE,PLOT_COLOURS)
      
      
    
  }
    else if (is.null(PLOT_COLOURS)==FALSE && ncol(TREATMENT_TABLE)>4){
      stop("There are too many columns in TREATMENT_TABLE (Max=4) for automatic colour and symbol selection.
            Please provide a seperate PLOT_COULOURS table with less columns.")
    
  
  }
  
  }else if (is.null(PLOT_COLOURS)==FALSE){
  
      print("PLOT_COLOURS was specified manually.")
    }
    
####################################################################################################################################
############################################## START CALCULATIONS  #################################################################
####################################################################################################################################
  
  TREATMENT_TABLE<-as.data.frame(TREATMENT_TABLE)
  rownames(TREATMENT_TABLE)<-TREATMENT_TABLE[,1]
  keep<-intersect(colnames(EXPRESSIONTABLE),rownames(TREATMENT_TABLE))
  TABLE<-EXPRESSIONTABLE[,keep ,drop=F]        
  bincodes<-as.vector(MAPPINGFILE$BINCODE)             
  identifiers<-as.vector(MAPPINGFILE$IDENTIFIER)
   
  MAPPINGFILE->Map1
  Map1[,3]<-as.vector(Map1[,3])
  #print(head(Map1))
  Map_ORA<-Map1[,c(1,3)]
  #print(head(Map_ORA))

  print('Start PCA calculations')
  
  PCA_SCALING
  
  if (is.null("PCA_SCALING")==FALSE && is.character(PCA_SCALING)==TRUE && (PCA_SCALING=="ROWS")==TRUE){
    
    print("Scaling on ROWS")
    PCA<-prcomp(t(TABLE), scale=TRUE, center=TRUE)
    Summary_PCA<-summary(PCA) 
    Scores_PCA<-PCA$rotation             
    Loadings_PCA<-PCA$x   
    
    
  }else if (is.null("PCA_SCALING")==FALSE && is.character(PCA_SCALING)==TRUE && (PCA_SCALING=="COLUMNS")==TRUE){
    print("Scaling on COLUMNS")
    PCA <- prcomp(TABLE, scale=TRUE, center=TRUE) 
    Summary_PCA<-summary(PCA)
    Scores_PCA<-PCA$x                                                              
    Loadings_PCA<-PCA$rotation                                                     
 
  
  }else if(is.null("PCA_SCALING")==TRUE) {
    print("Default scaling on COLUMNS")
    PCA <- prcomp(TABLE, scale=TRUE, center=TRUE) 
    Summary_PCA<-summary(PCA)
    Scores_PCA<-PCA$x                                                              
    Loadings_PCA<-PCA$rotation 
    
  } else {
    print("Default scaling on COLUMNS")
    PCA <- prcomp(TABLE, scale=TRUE, center=TRUE) 
    Summary_PCA<-summary(PCA)
    # Ausgabe von Scores und Loadings 
    Scores_PCA<-PCA$x                                                              
    Loadings_PCA<-PCA$rotation 
  }
  
  print("PCA calculations completed. :)")
  
   
  if(is.null(PCA_SCORE_TABLE)==TRUE || PCA_SCORE_TABLE==1){
    print("Start annotation Scores PCA and build Outputfile.")
    
    Annotation_Genes<-function(GENELIST){
      MAP=MAPPINGFILE
      PCA_Scores=Scores_PCA[,1:3]
      annotated_genes<-c()
      
      for (g in GENELIST){
        map<-MAP[MAP$IDENTIFIER==g,]
        genes<-PCA_Scores[rownames(PCA_Scores)==g,]
  
        if (nrow(map)==0){
          map<-rbind(map,c(rep('No Mapman Annotation',2),g,'No Mapman Annotation'))
          colnames(map)<-c("BINCODE","NAME","IDENTIFIER","DESCRIPTION")
          both<-cbind(map, t(as.data.frame(genes)))
          annotated_genes<-rbind(annotated_genes,both)
        
        }else{
          both<-cbind(map,matrix(genes, nrow=nrow(map), ncol=length(genes), byrow=TRUE))
          colnames(both)<-c(colnames(map),'PC1','PC2','PC3')
          annotated_genes<-rbind(annotated_genes,both)
         
        }
      }
      return(annotated_genes)
    }
  
    annotated_Scores_list<-lapply(rownames(Scores_PCA),function(x) Annotation_Genes(x)) 
    
    annotated_Scores<-bind_rows(annotated_Scores_list)
    
    name_save<-paste0('Scores_PCA_PC1-PC3_Annotated_',NAME_TABELLE,'.txt')
    
    write.table(file=name_save, annotated_Scores,sep="\t",row.names = F)
    
    rm(name_save)
    
    print("END annotation Scores PCA.")
  
  }else if(PCA_SCORE_TABLE==2){
    print("No Annotation of PCA genes was build.")
  }
  
  # SCORE-PLOT: 
  
  cat("Principal Component Analysis\n\n", file = Name_Info, append = TRUE) 
  #capture.output(Name_Info, file = Name_Info, append = TRUE)
  cat("Summary PCA:\n", file = Name_Info, append = TRUE)
  capture.output(Summary_PCA,file = Name_Info, append = TRUE)
  cat("\n", file = Name_Info, append = TRUE)
  
  Importance<-Summary_PCA$importance
  Importance[2,1]-> Importance_PC1
  Importance[2,2]-> Importance_PC2
  Importance[2,3]-> Importance_PC3
 
  cat(paste0('PC1 trägt ', Importance[2,1], ' zur Erklärung der Gesamtvarianz bei.\n'),file = Name_Info, append = TRUE)
  cat(paste0('PC2 trägt ',Importance[2,2],' zur Erklärung der Gesamtvarianz bei.\n'),file = Name_Info, append = TRUE)
  cat(paste0('PC1 und PC2 zusammen tragen ',Importance[3,2],' zur Erklärung der Gesamtvarianz bei.\n'),file = Name_Info, append = TRUE)
  cat(paste0('PC3 trägt ',Importance[2,3],' zur Erklärung der Gesamtvarianz bei.\n'),file = Name_Info, append = TRUE)
  cat(paste0('PC1, PC2 und PC3 zusammen tragen ',Importance[3,3],' zur Erklärung der Gesamtvarianz bei.\n'),file = Name_Info, append = TRUE)

  
  Importance_PC1 <- round(Importance_PC1*100, digits=1)
  xlab1 <-paste('PC1 ','(',Importance_PC1,' %)')
  
  Importance_PC2 <- round(Importance_PC2*100, digits=1)
  ylab1 <-paste('PC2','(',Importance_PC2,' %)')
  
  Importance_PC3 <- round(Importance_PC3*100, digits=1)
  ylab2 <-paste('PC3','(',Importance_PC3,' %)')
  
  
  xlab<-c(xlab1,ylab1)
  ylab<-c(ylab1,ylab2)
  
  Name_Info_pdf<-paste0(NAME_TABELLE,'.pdf')
  pdf(Name_Info_pdf)
  
 
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste(NAME_TABELLE), cex = 1.6, col = "black")
  
  
  plot(PCA, type="l")
  plot(Scores_PCA,main="Scores PCA: Genes in experimental space")
  biplot(PCA, cex=0.9, las=1)
  
  # PC1 x PC2
  plot(Loadings_PCA, type='n',xlab=xlab1, ylab=ylab1,main='PC1 x PC2')
  text(Loadings_PCA, labels=rownames(Loadings_PCA), cex=0.8)
  
  # PC1 x PC3
  plot(Loadings_PCA[,c(1,3)], type='n',xlab=xlab1, ylab=ylab2,main='PC1 x PC3')
  text(Loadings_PCA[,c(1,3)], labels=rownames(Loadings_PCA[,c(1,3)]), cex=0.8)

  # PC2 x PC3
  plot(Loadings_PCA[,c(2,3)], type='n',xlab=ylab1, ylab=ylab2,main='PC2 x PC3')
  text(Loadings_PCA[,c(2,3)], labels=rownames(Loadings_PCA[,c(2,3)]), cex=0.8)
  
  dev.off()
  
  
  par(mfrow=c(1,1))
  cat("\n", file = Name_Info, append = TRUE)
  
  
  #Maximal and Minimal PC1:
  cat('Maximal und Minimalwerte der PC1:\n',file = Name_Info, append = TRUE)
  PC1<-summary(sort(PCA$rotation[,1]))
  capture.output(PC1,file = Name_Info, append = TRUE)
  cat("\n", file = Name_Info, append = TRUE)
  
  #Maximal and Minimal PC2:')
  cat('Maximal und Minimalwerte der PC2:\n',file = Name_Info, append = TRUE)
  PC2<-summary(sort(PCA$rotation[,2]))
  capture.output(PC2,file = Name_Info, append = TRUE)
  cat("\n", file = Name_Info, append = TRUE)
  
  #Maximal and Minimal PC3:')
  cat('Maximal und Minimalwerte der PC3:\n',file = Name_Info, append = TRUE)
  PC3<-summary(sort(PCA$rotation[,3]))
  capture.output(PC3,file = Name_Info, append = TRUE)
  
  
  
  rm(PC1,PC2,PC3)
  
  name_tab1a<-paste0(NAME_TABELLE,'_PC1_pos','.txt')                    
  name_tab1b<-paste0(NAME_TABELLE,'_PC1_neg','.txt')
  name_tab2a<-paste0(NAME_TABELLE,'_PC2_pos','.txt')
  name_tab2b<-paste0(NAME_TABELLE,'_PC2_neg','.txt')
  name_tab3a<-paste0(NAME_TABELLE,'_PC3_pos','.txt')
  name_tab3b<-paste0(NAME_TABELLE,'_PC3_neg','.txt')
  
  names_tabs<-c(name_tab1a,name_tab1b,name_tab2a,name_tab2b,name_tab3a,name_tab3b)
  
      
  name_tab1<-paste0(NAME_TABELLE,'_PC1_pos_neg','.txt')
  name_tab2<-paste0(NAME_TABELLE,'_PC2_pos_neg','.txt')
  name_tab3<-paste0(NAME_TABELLE,'_PC3_pos_neg','.txt')
  
  names_tabs2<-c(name_tab1,name_tab2,name_tab3)
  

  if (is.null(TOP)==TRUE){
    TOP=50
  }else if(is.null(TOP)==FALSE && TOP<=1000){
    TOP=TOP
    print(paste("You've choosen ",TOP, " highly involved genes in PC1-PC3 for ORA calculation.", sep="" ))
  }else if(is.null(TOP)==FALSE && TOP>1000 || is.null(TOP)==FALSE && TOP>nrow(Scores_PCA)){
    stop(paste('Number of Genes for TOP is too high. Maximum number of genes for PCA is ',nrow(Scores_PCA),'. Please choose a value less/equals 1000', sep=""))
  }

  tabs<-vector()
  tabs_plot<-vector()
  cnams<-vector()
  cnams2<-vector()
  for (x in 1:3){
      n<-as.data.frame(sort(Scores_PCA[,x],decreasing=T))
      n1<-as.data.frame(sort(Scores_PCA[,x],decreasing=F))
      tab<-unique(rownames(n)[1:TOP])
      tab2<-unique(rownames(n1)[1:TOP])
      
      tabs<-cbind(tabs,tab,tab2)
      cnams<-c(cnams,paste0('Identifier_PC',x,'_pos'),paste0('Identifier_PC',x,'_neg'))
      
      tab_p<-Map1[Map1[,3]%in%tab,3:4] 
      tab_p2<-Map1[Map1[,3]%in%tab2,3:4]
      
      tabs_plot<-c(tabs_plot,tab_p[!duplicated(tab_p$IDENTIFIER),],tab_p2[!duplicated(tab_p2$IDENTIFIER),])
      cnams2<-c(cnams2,paste0('Identifier_PC',x,'_pos'),paste0('Description_PC',x,'_pos'),paste0('Identifier_PC',x,'_neg'),paste0('Description_PC',x,'_neg'))
      names(tabs_plot)<-cnams2
     }
      
 colnames(tabs)<-cnams 

  #### Boxplots Gene and TPM 
  TPM<-cbind(rownames(TPM),TPM)
  TPM_tabs<-vector() 
  
  ##subset der TPM Tabelle: TPM
  
  if (is.null(TPM)==FALSE){
    for (j in 1:ncol(tabs)){
	  subset_TPM<-TPM[TPM[,1] %in% unique(tabs[,j]),]
	  
	 
	  Annotation_Genes2<-function(GENELIST){
	      MAP=MAPPINGFILE
	      TPM_sub=subset_TPM
	      annotated_genes2<-c()
	      
	   for (g in GENELIST){
    		map<-MAP[MAP$IDENTIFIER==g,]
    		genes<-TPM_sub[TPM_sub[,1]==g,-1]
    	  
    		if (nrow(map)==0){
    		  map<-rbind(map,c(rep('No Mapman Annotation',2),g,'No Mapman Annotation'))
    		  colnames(map)<-c("BINCODE","NAME","IDENTIFIER","DESCRIPTION")
    		  both<-cbind(map, as.data.frame(genes))
    		  colnames(both)<-c(colnames(map),colnames(subset_TPM[,-1]))
    		  annotated_genes2<-rbind(annotated_genes2,both)
    		}else{
    		  both<-cbind(map,as.data.frame(t(replicate(nrow(map), unlist(genes[1,])))))
    		  colnames(both)<-c(colnames(map),colnames(subset_TPM[,-1]))
    		  annotated_genes2<-rbind(annotated_genes2,both)
    		}
	      }
	      return(annotated_genes2)
	    }
	  
	  annotated_TPM_list<-lapply(tabs[,j],function(x) Annotation_Genes2(x))
	  annotated_TPM<-data.frame(Reduce(rbind, annotated_TPM_list))
	  annotated_TPM<-annotated_TPM[!duplicated(annotated_TPM[,3]),]
	  #annotated_TPM<-as.data.frame(bind_rows(annotated_TPM_list))
	  write.table(file=paste0('Top',TOP,'_Genes_',cnams[j],'_Annotation+TPM.txt'),sep="\t",annotated_TPM, row.names = F)
	  
	  ##### Boxplots Top Gene#####
	  Name_Boxplots_pdf<-paste0(NAME_TABELLE,'_Boxplots_Top',TOP,'_Genes_',cnams[j],'.pdf')
	  pdf(Name_Boxplots_pdf)
	  
	  for (n in 1:nrow(annotated_TPM)){
	    gene<-annotated_TPM[n,3]
	    description<-str_split(annotated_TPM[n,4],'&')[[1]][1]
	    
	    TPM_tab<-as.data.frame(cbind(colnames(TPM)[-1], GROUPED_TABLE,as.matrix(t(annotated_TPM[,c(-1:-4)]))))
	    colnames(TPM_tab)<-c('samples','group',annotated_TPM[,3])
	    TPM_tab<-TPM_tab[,-1]
	    #for (m in 2:ncol(TPM_tab)){
	  
	    boxplot(as.numeric(TPM_tab[,gene])~as.character(TPM_tab[,1]),
	            data = TPM_tab,
	            xaxt = "n",
	            ylab="log2CPM", 
	            xlab="Treatment", main=paste0(gene, '\n',description),cex.main=0.7,
	            col=c(rep(c("green","red"),length(unique(as.character(TPM_tab[,2]))))))

	    text(x = 1:length(levels(GROUPED_TABLE)),
	         ## Move labels to just below bottom of chart.
	         y = par("usr")[3] - 0.45,
	         ## Use names from the data list.
	         labels = levels(GROUPED_TABLE),
	         ## Change the clipping region.
	         xpd = NA,
	         ## Rotate the labels by 35 degrees.
	         srt = 35,
	         ## Adjust the labels to almost 100% right-justified.
	         adj = 0.965,
	         ## Increase label size.
	         cex = 0.7)
	   
	   # }
	    
	    
	  }
	  dev.off()
	
	  
	  TPM_tabs<-rbind(TPM_tabs,annotated_TPM)
	  
  	
	}
  }else{
  
  stop('There is no TPM table provided!')
  }

  direc<-getwd()
  
 print('END PCA CALCULATIONS')
 #################################################### OVERREPRESENTATION ANALYSIS #####################################################################  

print('START ORA CALCULATIONS')  

 Read_line<-function(bincode,identifier,LISTE2){        
  if (identifier == ''){
    return(LISTE2);}
  if (bincode %in% names(LISTE2)){                  
    LISTE2[[bincode]]<-c(identifier, LISTE2[[bincode]]) 
  } else {
    LISTE2[[bincode]]<-c(identifier);  
  }
  newbincode<-sub("\\.[^\\.]+$","", bincode)  
  if (bincode!=newbincode){
    LISTE2<-Read_line(newbincode,identifier,LISTE2)    
  } 
  return(LISTE2);
}

LISTE<-list()
for (i in 1:nrow(Map_ORA)){
  LISTE<-Read_line(bincodes[i],identifiers[i], LISTE)
}
#print(head(LISTE))

Counts_elements<- function(MAPPINGLISTE,GENES,RESULT){     
  rownames(RESULT)<-names(MAPPINGLISTE)                   
  for (h in 1:length(MAPPINGLISTE)){                      
    ol<-intersect(MAPPINGLISTE[[h]], GENES)               
    RESULT[h,1]<-length(GENES)                            
    RESULT[h,2]<-length(ol)                               
    RESULT[h,3]<-RESULT[h,1]-RESULT[h,2]                  
    
  }
  return(RESULT)
}

dataframe<-data.frame(matrix(0,nrow=length(LISTE),ncol=3))

contigenz<-function(x){                                   
  x<-as.vector(x)
  res<-fisher.test(matrix(c(x[2],x[3],x[4],x[5]),nrow=2)) 
  return(res[[1]])                                        
}

print('ORA Erstellen der Tabelle für ALLE Gene')

ALLE<-Counts_elements(LISTE, ALL_GENES[,1],dataframe)           

if (is.null(ALLE)==TRUE){
  stop("ALLE table is eqals zero. Please check your input.")
}


print('Build endtable ORA')

tabs_ORA<-tabs


pvalues<-data.frame(matrix(0,nrow=length(LISTE),ncol=ncol(EXPRESSIONTABLE))) 
ERF_TAB<-data.frame(matrix(0,nrow=length(LISTE),ncol=1)) 
for (k in 1:ncol(tabs_ORA)){                              
  rownames(pvalues)<-names(LISTE)                         
  colnames(pvalues)<-names_tabs                           
  rownames(ERF_TAB)<-names(LISTE)                         
  Tabs<-as.vector(tabs_ORA[,k])                        
  MEINS<-Counts_elements(LISTE, Tabs,dataframe)     
  NMEINS<-ALLE$X2-MEINS$X2                                
  NMEINSSUM<-ALLE$X3-MEINS$X3                             
  MEINS_ERF<-as.data.frame(MEINS$X2/MEINS$X1)                                   
  ALLE_ERF<-as.data.frame(ALLE$X2/ALLE$X1)                             
  ERF<-MEINS_ERF[,1]/ALLE_ERF[,1]                                
  FINAL<-cbind(MEINS,NMEINS,NMEINSSUM)                   
  pvalues[k]<-apply(FINAL,1,contigenz)                    
  ERF_TAB<-cbind(ERF_TAB,ERF)
  
}
ERF_TAB<-ERF_TAB[,-1]                                     
colnames(ERF_TAB)<-names_tabs  

print("Enrichmentfactor was calculated")


ORA_results<-data.frame(matrix(0,nrow=length(LISTE),ncol=1))
ORA_names<-vector()
for (g in 1:ncol(ERF_TAB)){ #)
  rownames(ORA_results)<-names(LISTE)
  ORA_results<-cbind(ORA_results,ERF_TAB[,g],pvalues[,g])
  ORA_names<-c(ORA_names,paste0(colnames(ERF_TAB[g]),'_',colnames(pvalues[g])),paste0('Pvalue_',colnames(pvalues[g])))
}
ORA_results<-ORA_results[,-1]              
colnames(ORA_results)<-ORA_names           
Bincodes<-c(rownames(ORA_results))         
ORA_results<-cbind(Bincodes,ORA_results)  

print('END: ORA CALCULATIONS')

print('START: MERGE TABLES')

Finaldata_TPM_over<-data.frame(matrix(0,nrow=6,ncol=1)) 
m=1 

for (bi in c(2,4,6,8,10,12)){
  name_PC<-c('PC1_pos','PC1_neg','PC2_pos','PC2_neg','PC3_pos','PC3_neg')
  name_table<-c("BINCODES","NAME", "IDENTIFIER", "DESCRIPTION", paste0("ERF_",name_PC[m]),paste0("pvalue_",name_PC[m]))
  
  #get only the PC column and their ERF/Pvalues
  Tabelle<-ORA_results[,c(1,bi,bi+1)]
  
  #Filter the dataframe by ERF and pvalue
  Sorting_ERF<-Tabelle[which(Tabelle[,2]>=1),]
  Sorting_pval<-Sorting_ERF[which(Sorting_ERF[,3]<=0.05),] 
  
  
  Bins_sort<-Sorting_pval[,1]
  Mapping<-subset(Map1, Map1[,1]%in% Bins_sort)
  
  tot<-data.frame()
  Map_odupli<-Mapping[!duplicated(Mapping[,2]),]
  for (n in  Map_odupli[,1]){
    binc<-Mapping[which(Mapping[,1]==n),]
    if (nrow(binc)==1){
      binc<-binc[,-5]
      colnames(binc)<-colnames(Map_odupli[,-5])
      tot<-rbind(tot,binc)
    }else{
      genes<-paste(binc[,3], collapse = ", ") 
      co<-cbind(binc[1,1:2],genes,binc[1,4])
      colnames(co)<-colnames(Map_odupli[,-5])
      tot<-rbind(tot,co)
    }
  }
  
  #adding ERF and pvalue column
  tot3<-data.frame()
  for (p in tot[,1]){
    t<-Sorting_pval[p,] # ERF and pvalue for p
    z<-tot[which(tot[,1]==p),]
    
    com<-cbind(z[,1:4],t[,2:3])
    tot3<-rbind(tot3,com)
  }
  
  frame<-cbind(name_table,t(tot3))
  Finaldata_TPM_over<-cbind(Finaldata_TPM_over,frame)
  
  m<-m+1
}


# Output for Plot
save_tab<-t(Finaldata_TPM_over)[-1,]
name_save<-paste0('ORA_Toplist_overrepresentation_',NAME_TABELLE,'.txt')

write.table(file=name_save, save_tab,sep="\t",row.names = F)
print('ENDE: ZUSAMMENF?GEN DER ORA SPALTEN, ERSTELLEN EINER GESAMTDATEI')


##########################################################################################################################################################
#########################################################   PLOT einzeln, ENDPLOT   ######################################################################
##########################################################################################################################################################

print('START PLOT')
GROUPED_TABLE ->Treatment
#Treatment<-as.matrix(Treatment1)
#Treatment<-as.factor(t(Treatment$Treatment))
Treatment_Levels <- levels(Treatment)

cat("\n", file = Name_Info, append = TRUE)
cat("Treatment levels:\n", file = Name_Info, append = TRUE)
capture.output(levels(Treatment),file = Name_Info, append = TRUE)

Symbols<-as.numeric(as.vector(PLOT_COLOURS$Symbols))
Symbolcolour<-as.vector(PLOT_COLOURS$Symbolcolour)
Bordercolour<-as.vector(PLOT_COLOURS$Bordercolour)

#Legend_parameters<-TREATMENT_TABLE[!duplicated(TREATMENT_TABLE[,2:8]),]
Legend<-cbind(as.character(Treatment1),PLOT_COLOURS)  
Legend_parameters<-as.data.frame(Legend[!duplicated(Legend[,1]),])

Symbols_leg <-as.numeric(as.vector(Legend_parameters$Symbols))
Colour_leg <-as.vector(Legend_parameters$Symbolcolour) 
Border_leg<-as.vector(Legend_parameters$Bordercolour)
labels_leg<-as.vector(Legend_parameters[,1])


####################################### BARPLOT Parameter ##################################################################################
############################################################################################################################################
ORA_results<-ORA_results[order(ORA_results[,1]),]
Bin_descripton<-Map1[Map1[,1]%in%ORA_results[,1],]
Bin_descripton<-Bin_descripton[order(Bin_descripton[,1]),]

Bin_descripton<-Bin_descripton[!duplicated(Bin_descripton[,1]),]
ORA_results<-ORA_results[match(Bin_descripton[,1],ORA_results[,1]),]
DESCRIPTION<-as.vector(Bin_descripton[,2])

END<-cbind(ORA_results[,1],Bin_descripton[,c(-1,-5)],ORA_results[,-1])
ERF_Liste<-list()


for (y in 5:ncol(END)){
  if ((y %% 2 ==0) ==FALSE){
    #print(y)
    #y=5
    Test_y<-END[,c(1:4,y,y+1)]
    Test<-filter(Test_y, Test_y[,5]>=1,) #[c(1:4,y,y+1)]
    Test<-list(Test[Test[,6]<=0.05,]) 
    ERF_Liste<-cbind(ERF_Liste,Test) 
  }
}

total_plot<-function(PCxpos,PCxneg,PCypos,PCyneg,PCx,PCy,ERF_COLUMN_ypos,ERF_COLUMN_yneg,ERF_COLUMN_xpos,ERF_COLUMN_xneg,TOP_NUMBER,XLAB,YLAB){  
  
 # PC1xPC2<-total_plot(1,2,3,4,1,2,5,xlab[1],ylab[1])
  #PCxpos=tabs_plot[1:2]
  #PCxneg=tabs_plot[3:4]
  #PCypos=tabs_plot[5:6]
  #PCyneg=tabs_plot[7:8]
  #PCx=1
  #PCy=2
  #ERF_COLUMN_ypos=3
  #ERF_COLUMN_yneg=4
  #ERF_COLUMN_xpos=1
  #ERF_COLUMN_xneg=2
  #TOP_NUMBER=5
  #XLAB=xlab[1]
  #YLAB=ylab[1]
  
  top<-TOP_NUMBER
  
  colfunc_pos_y <-colorRampPalette(c('indianred2','red4'))  
  colfunc_neg_y<-colorRampPalette(c('navyblue','lightskyblue1'))
  colfunc_pos_x<- colorRampPalette(c('red4','indianred2')) 
  colfunc_neg_x<-colorRampPalette(c('lightskyblue1','navyblue'))
  
  topgenes_table<-function(PCXPOS,PCXNEG,PCYPOS,PCYNEG){
    PCxp=as.data.frame(PCXPOS)
    PCxn=as.data.frame(PCXNEG)
    PCyp=as.data.frame(PCYPOS)
    PCyn=as.data.frame(PCYNEG)
    
    l1<-cbind(PCxp[1:5,],PCxn[1:5,],PCyp[1:5,],PCyn[1:5,]) 
    l<-l1[,c(2,4,6,8)]
    l[,1]<-as.character(l[,1])
    l[,2]<-as.character(l[,2])
    l[,3]<-as.character(l[,3])
    l[,4]<-as.character(l[,4])
    
    top_axes<-vector()
    for (s in 1:ncol(l)){
      top_plot<-vector()
      for (u in 1:nrow(l)){                                   
        nam<-strsplit(l[u,s], " ",fixed=T)
        if  (grepl(nam[[1]][1], "(oriGinal",fixed=T) | grepl(nam[[1]][1],"(original",fixed=T) | grepl(nam[[1]][2],"(oriGinal",fixed=T) | grepl(nam[[1]][2],"(original",fixed=T)){
          top_plot<-rbind(top_plot,nam[[1]][1])
        }else if (grepl(nam[[1]][1],"(oriGinal",fixed=T) | grepl(nam[[1]][1],"(original",fixed=T)){
          top_plot<-rbind(top_plot,"no hits")
        }else if(grepl(nam[[1]][1], "(oriGinal",fixed=T)==F | grepl(nam[[1]][1],"(original",fixed=T)==F | grepl(nam[[1]][2],"(oriGinal",fixed=T)==F | grepl(nam[[1]][2],"(original",fixed=T)==F){
          
          top_plot<-rbind(top_plot,paste(nam[[1]][1],nam[[1]][2]))
        }
      }
      top_axes<-cbind(top_axes,top_plot)
    }
    colnames(top_axes)<-colnames(l1[,c(1,3,5,7)])
    
    top_x_axis<-as.vector(cbind(top_axes[1:5,2],top_axes[c(5,4,3,2,1),1]))#c(5,4,3,2,1)
    top_y_axis<-as.vector(cbind(top_axes[1:5,4],top_axes[c(5,4,3,2,1),3]))#c(5,4,3,2,1)
    
    pos_neg_top_data<-cbind(top_x_axis,top_y_axis)
    pos_neg_top_data<-cbind(top_x_axis,top_y_axis)
    
    #print(pos_neg_top_vec)
    return(pos_neg_top_data)
    
  }
  barplot_table<-function(table2,TRUEorFALSE){
    Porder<-list()
    p_ordered<-table2[order(table2[,6],decreasing = TRUE),]    
    p_ordered[,2]<-as.character(p_ordered[,2])                        
    names_plot<-list()
    names1<-vector()
    bins<-data.frame()
    for (u in 1:nrow(p_ordered)){                                 
      nam<-strsplit(p_ordered[u,2], " " ,fixed=T)
      names_plot<-rbind(names_plot,nam)
    }
    
    for (l in 1:nrow(names_plot)){
      if (is.na(names_plot[[l]][4])==F){
        
        newnames<-paste0(names_plot[[l]][1]," ",names_plot[[l]][2]," ",names_plot[[l]][3],"*")
        
        if (l==1){
          
          names1<-rbind(names1,newnames)
          bins<-rbind(bins, p_ordered[l,])
          
        }else if(grepl(newnames,names1[l-1],fixed = T)==FALSE){
          
          names1<-rbind(names1,newnames)
          bins<-rbind(bins, p_ordered[l,])
          
        }else{
          print(paste0("Shortened Bin name was duplicated and removed but Sub-Bin might be different. Check output file for more information about"," Bin: ",paste(names_plot[[l]][1],names_plot[[l]][2],names_plot[[l]][3],names_plot[[l]][4],names_plot[[l]][5])))
        }
      }else if(is.na(names_plot[[l]][3])==F && is.na(names_plot[[l]][4])==T){
        
        names1<-rbind(names1,paste(names_plot[[l]][1],names_plot[[l]][2],names_plot[[l]][3]))
        
        bins<-rbind(bins, p_ordered[l,])
        
      }else if (is.na(names_plot[[l]][3])==T && is.na(names_plot[[l]][2])==F){
        
        names1<-rbind(names1,paste(names_plot[[l]][1],names_plot[[l]][2]))
        bins<-rbind(bins, p_ordered[l,])
        
      }else if(is.na(names_plot[[l]][3])==T && is.na(names_plot[[l]][2])==T){
        names1<-rbind(names1,names_plot[[l]][1])
        bins<-rbind(bins, p_ordered[l,])
      }
    }
    
    
    new_ordered<-cbind(bins[,1],names1,bins[,3:6])
    
    
    names_plot<-as.data.frame(names_plot)
    all_plot<-new_ordered[,c(2,5,6)]
    all_without_dupl<-all_plot[!duplicated(all_plot[,1]),]
    all_without_dupl[,1]<-as.character(all_without_dupl[,1])
    for (o in TRUEorFALSE){
      if (o == TRUE){
        all_without_dupl<-all_without_dupl[(nrow(all_without_dupl)-(top-1)):nrow(all_without_dupl),]
      }else{
        all_without_dupl<-all_without_dupl[1:top,]
      }
    }
    return(all_without_dupl)
  }
  
  n<-layout(matrix(c(0,1,1,2,
                     1,1,1,3,
                     4,5,5,6,
                     4,5,5,6,
                     0,7,7,6),4,5),
            widths = c(1.7,1.7,2.5,2.5,1.7), 
            heights = c(2,2.5,2.5,2.7)) 
  
   cat('\n\n',file = Name_Info, append = TRUE)
  cat('Plot ORA and PCA results',file = Name_Info, append = TRUE)
  
  ################################### ORA_BARPLOT y-axes

  
  #1 Barplot y-axes
  p<-as.data.frame(ERF_Liste[ERF_COLUMN_ypos])
  #1 Barplot y-Achse
  p2<-as.data.frame(ERF_Liste[ERF_COLUMN_yneg])
  
  
  ypos<-barplot_table(p,TRUE)    
  ypos[,1]<-sapply(ypos[,1], function(x) gsub( "\\."," ", x))
  ypos[,2]
  cat('Pathways y-axis positive direction\n', file = Name_Info, append = TRUE)
  capture.output(ypos[,c(1,3)], file = Name_Info, append = TRUE) 
  
  
  yneg<-barplot_table(p2,TRUE)
  yneg<-yneg[order(yneg[,3],decreasing = F),]
  yneg[,1]<-sapply(yneg[,1], function(x) gsub( "\\."," ", x))
  
  cat('Pathways y-axis negative direction\n', file = Name_Info, append = TRUE)
  capture.output(yneg[,c(1,3)], file = Name_Info, append = TRUE)
  
 
  par(mar=c(3,23,14,0)+.3)
  m<-barplot(c(-log10(yneg[,2]),-log10(ypos[,2])), ylim = rev(range(c(0,log10(ypos[,2])))),yaxt='n',horiz = T,
             col =c(colfunc_pos_y(top),colfunc_neg_y(top)),axes=F,space=0.1,width=.17)
  
  mtext(side = 3, text = "log10 Enrichmentfactor", line = 1.4,cex = 0.8,adj = 0.7)
  axis(3,at=-round(max(c(log10(ypos[(nrow(ypos)-4):nrow(ypos),2]),log10(yneg[1:top,2])))):0,labels = round(max(c(log10(ypos[(nrow(ypos)-4):nrow(ypos),2]),log10(yneg[1:top,2])))):0, line = -1.)
  
  text(y=m, (par("usr")[1]-0.5),labels = as.vector(c(ypos[,1],yneg[,1])) , offset=-7.5,srt = 15, pos = 2, xpd =T,  cex = 1,font = 3)
  mtext(side = 3, text = "ORA Analysis", line = 5, cex = 1, font=2,adj=0.7)
  
  print("Sorted boxplots y-axis: smallest pvalues should be at the beginning and end!")
  print(c(yneg[,3],ypos[,3]))
  
  ############################################# pvalues legend ###############################################################################
  p3<-as.data.frame(ERF_Liste[ERF_COLUMN_xpos]) 
  p4<-as.data.frame(ERF_Liste[ERF_COLUMN_xneg]) 
  

  xpos<-barplot_table(p3,TRUE)
  xpos[,1]<-sapply(xpos[,1], function(x) gsub( "\\."," ", x))
  cat('Pathways x-axis positive direction\n', file = Name_Info, append = TRUE)
  capture.output(xpos[,c(1,3)], file = Name_Info, append = TRUE) 
  
  xneg<-barplot_table(p4,FALSE)
  xneg[,1]<-sapply(xneg[,1], function(x) gsub( "\\."," ", x))
  #einf?gen der Werte von xneg etc uz der INFO file
  cat('Pathways x-axis negative direction\n', file = Name_Info, append = TRUE)
  capture.output(xneg[,c(1,3)], file = Name_Info, append = TRUE)
  
  
  pos_pval<-c(format(min(c(ypos[1,3],xpos[1,3])),scientific = TRUE, digits = 3),format(min(c(ypos[3,3]),xpos[3,3]),scientific = TRUE, digits = 3),
              format(min(c(ypos[5,3]),xpos[5,3]),scientific = TRUE, digits = 3))
  
  neg_pval<-c(format(min(c(yneg[5,3]),xneg[5,3]),scientific = TRUE, digits = 3),format(min(c(yneg[3,3]),xneg[3,3]),scientific = TRUE, digits = 3),
              format(min(c(yneg[1,3],xneg[1,3])),scientific = TRUE, digits = 3))
  
  
  #####################################  Legend der pvalues
  #2 Legenden
  #Legende pvalue
  par(mar=c(7,5,2,0))
  legend_image <- as.raster(matrix(colfunc_neg_x(top), ncol=1))
  legend_image2<- as.raster(matrix(colfunc_pos_y(top), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'pvalue',cex.main=1.5,adj=0.5,bty='l')
  rasterImage(legend_image, 1.0,0,1.4,1) 
  rasterImage(legend_image2, 0.6,0,1,1)
  text(x=1.5, y = seq(0,1,l=3),cex=1, labels = neg_pval,adj = -0.05)
  text(x=0.2, y = seq(0,1,l=3),cex=1, labels = pos_pval,adj =)
  box(lty = '21', col = 'black')
  
  ######################################  Legend PCA
  
  
  if ((length(labels_leg) >= 20)== TRUE){
    par(mar=c(15,0,1,0))
    plot(c(1,2,3,4),c(1,1,1,1), axes = F, ylab = '',xlab = '',col='white')
    legend(x=1.3,y=1.9,box.lty=1, pch=Symbols_leg[1:20],xpd=T,cex=0.8,pt.lwd=2,pt.cex = 2,col=Border_leg[1:20], legend=labels_leg[1:20],
           pt.bg = Colour_leg[1:20],text.width=1) #,bty = '0'
    legend(x=2.7,y=1.9,box.lty=1, pch=Symbols_leg[21:length(Border_leg)],xpd=T,cex=0.8,pt.lwd=2,pt.cex = 2,col=Border_leg[21:length(Border_leg)],
           legend=labels_leg[21:length(Border_leg)],
           pt.bg = Colour_leg[21:length(Border_leg)],text.width=1) #,bty = '0'
    
    mtext(side = 3, text = "PCA", cex = 1, font=2,adj=0.54, line = 5)# 3 steht für die PC1
  }else if(( length(labels_leg)>=12 & length(labels_leg) <20)==TRUE){
    par(mar=c(8,0,4,1.5))
    plot(c(1,2,3,4),c(1,1,1,1), axes = F, ylab = '',xlab = '',col='white')
    legend(x=1.3,y=1.55,box.lty=1, pch=Symbols_leg[1:10],xpd=T,cex=1,pt.lwd=2.5,pt.cex = 2,col=Border_leg[1:10], legend=labels_leg[1:10],
           pt.bg = Colour_leg[1:10],text.width=1) #,bty = '0'
    legend(x=2.7,y=1.55,box.lty=1, pch=Symbols_leg[11:length(Border_leg)],xpd=T,cex=1,pt.lwd=2.5,pt.cex = 3,col=Border_leg[11:length(Border_leg)],
           legend=labels_leg[11:length(Border_leg)],
           pt.bg = Colour_leg[11:length(Border_leg)],text.width=1) #,bty = '0'
    mtext(side = 3, text = "PCA", cex = 1, font=2,adj=0.54, line = 2.5)
  }else{
    par(mar=c(6,0,0,0))
    plot(c(1,2,3,4),c(1,1,1,1), axes = F, ylab = '',xlab = '',col='white')
    legend('center',box.lty=1, pch=Symbols_leg,xpd=T,cex=1,pt.lwd=2.5,pt.cex = 3,col=Border_leg, legend=labels_leg,pt.bg = Colour_leg,text.width=1) #,bty = '0'
    mtext(side = 3, text = "PCA", cex = 1, font=2,adj=0.5, line = 1)
  }
  
  
  ###################################### Top5 Gene x-axes PC
   #3 Top5 gene PC1/PC2
  # noch zu bearbeiten
  
  PC_genes_xaxis<-topgenes_table(PCxpos,PCxneg,PCypos,PCyneg)  
  
  first_col<-as.data.frame(lapply(PC_genes_xaxis[,1], function(x) gsub( "\\(original description:","no hits", x)))
  first_col<-as.data.frame(lapply(first_col, function(x) gsub( "\\(original","", x)))
  first_col<-as.data.frame(lapply(first_col, function(x) gsub( "component of","", x)))
  
  second_col<-as.data.frame(lapply(PC_genes_xaxis[,2], function(x) gsub( "\\(original description:","no hits", x)))
  second_col<-as.data.frame(lapply(second_col, function(x) gsub( "\\(original","", x)))
  second_col<-as.data.frame(lapply(second_col, function(x) gsub( "component of","", x)))
  
  PC_genes_xaxis<-cbind(t(first_col),t(second_col))
  
  
  
  #par(mar=c(1,6,8,5))
  par(mar=c(0.5,6,5,5))
  x<-c(0,1,2,3,4,5,6,7,8,9,10,11)
  y<-c(0,0,0,0,0,0,0,0,0,0,0,0)
  
  
  plot(x,y, pch=1, cex= c(4.2,4.2,4.0,3.7,3.4,3.1,3.1,3.4,3.7,4.0,4.2,4.2),yaxt='n',axes=F,ylab='',xlab='',col =c('white',rep('navyblue',5),rep('red',5),'white'))
  text(x,y,adj = 0, labels =c('',PC_genes_xaxis[1:10,1],''), offset=2.6,srt =-12, pos = 3, cex = 0.95,font =3)
   
  ####################################### PCA PLot
  
  par(mar=c(6,5,0,2))
  plot(Loadings_PCA[,c(PCx,PCy)], xlab=XLAB, ylab=YLAB,pch=Symbols, cex.lab=1.7,cex=4, lwd=2.5,bg=Symbolcolour,col=Bordercolour,cex.axis=1.5)
  mtext(side = 3, text = paste0('PCA ',NAME_TABELLE), line = 0.5, cex = 1, font=2)#axis(4,at=5:1,labels=names_plot[1:5,1],line = -1.3)
  
  ####################################### Barplot entlang der x_Achse
 

  par(mar=c(13,1,2.5,18))
  m<-barplot(c(log10(xneg[,2]),log10(xpos[,2])), ylim = rev(range(c(0,log10(xpos[,2])))),yaxt='n',
             col =c(colfunc_neg_x(top),colfunc_pos_x(top)),space=0.1,width=.1)  
  text(x= m, (par("usr")[3]+0.4) ,adj = 0.5, labels =as.vector(c(xneg[,1],xpos[,1])), offset=-1,srt = 15, pos = 2, xpd = TRUE, cex = 1,font =3)
  #axis(4,at=round(max(c(log10(xneg[,2]),log10(xpos[,2])))):0,labels = round(max(c(log10(xneg[,2]),log10(xpos[,2])))):0, line = -2.5)
  axis(4,at=round(max(c(log10(xneg[,2]),log10(xpos[,2])))):0,labels = round(max(c(log10(xneg[,2]),log10(xpos[,2])))):0, line = -1.5)
  
  mtext(side = 4, text = "log10 Enrichmentfactor", line = 1,cex = 0.8)
  mtext(side = 4, text = "ORA Analysis", line = 5.0, cex = 1, font=2)
  
  print("Sorted boxplots y-axis: smallest pvalues should be at the beginning and end!")
  print(c(xneg[,3],xpos[,3]))
  
  ######################################### Top5 Gene entlang der y _Achse
  
  par(mar=c(8,0,0,0))
  y<-c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5)
  x<-c(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0)
  
  n<-plot(x,y, pch=1, cex= c(4.5,4.2,3.9,3.6,3.3,3.3,3.6,3.9,4.2,4.5,4.5),yaxt='n',axes=F,ylab='',xlab='',col=c(rep('navyblue',5),rep('red',5),'white'))
  
  text(x,y,adj = 2, labels = c(PC_genes_xaxis[1:10,2],''), offset=1,srt =-5, pos = 4, cex = 1,font =3)
  
  return(PC_genes_xaxis)  
}


PC1xPC2<-total_plot(tabs_plot[1:2],tabs[3:4],tabs_plot[5:6],tabs_plot[7:8],1,2,3,4,1,2,5,xlab[1],ylab[1])
PC1xPC3<-total_plot(tabs_plot[1:2],tabs[3:4],tabs_plot[9:10],tabs_plot[11:12],1,3,5,6,1,2,5,xlab[1],ylab[2])
PC2XPC3<-total_plot(tabs_plot[5:6],tabs_plot[7:8],tabs_plot[9:10],tabs_plot[11:12],2,3,5,6,3,4,5,xlab[2],ylab[2])

print('ENDE PLOT')

}

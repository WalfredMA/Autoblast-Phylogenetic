library(XML)
library(RCurl)

####parameters#################################################
	 x<-"NP_034283.1"
	 program<-"blastp"
	 database<-"nr"
	 maxnumtarget<-1000
	 runpsi<-'True'
	 psithresh<-0.005
	 runtime<- 10
	 findsimilar<-'FALSE'
###############################################################

	 
	 dir.create("result")
	 similar <- 0
	 i <- 1
	 l<-0
     numberofalignmenttable <-c()
     numberofalignmenttable[0] <- 0
     numberofalignmentrecord <-c()
     numberofalignmentrecord[1] <- 0

     
	 if (findsimilar == 'True') {runtime <- -1} 
	 
	 if (program == "blastn" | program=="tblastn" | program =="tblastx") 
	 {page <-"Nucleotides"} else {page="Proteins"}
	 	 
	 if (runpsi == 'True') {runpsi <- "psiBlast"} else{runpsi <- "Blastp"}
	 
	 if(page=="Proteins"){
	     	type2<-"PROTEIN"} else{type2<-"DNA"}

############## main loop ######################################	 


     while (runtime != 0 && similar == 0 ) {
     	

##############call for blast###################################    


	 
	     seq <- x
	      if (i == 1 ){
	     url1 <- paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?QUERY=",seq,"&DATABASE=",database,"&BLAST_PROGRAMS=",runpsi,"&MAX_NUM_SEQ=",maxnumtarget,"&I_THRESH=",psithresh,"&FILTER=L&EXPECT=10&PROGRAM=",program,"&CLIENT=web&SERVICE=plain&NCBI_GI=on&PAGE=",page,"&CMD=Put",sep="")} else {
	     url1 <- paste("http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Put&VIEW_RESULTS=FromRes&NEXT_I=Go&CDD_RID=data_cache_seq:",seq,"&PREV_RID=",IDout,"&RUN_PSIBLAST=on&STEP_NUMBER=",i-1, "&UNIQ_OBJ_NAME=A_SearchResults_1bekt0_1dwT_draFWZwu3F7_23to07_1q8Pon&QUERY_INDEX=0", sep="")}
	     
	     
	     
###############scan page and get quest ID#######################	  	   
  
	     t<-0
	     repeat {
	     if (t == 10) {break}
	     t <- t + 1
	     blast1 <- try(scan(url1,what="raw"))
	     if ('try-error' %in% class(blast1)) 
	     {Sys.sleep(5)	
	     next}else {break}
	     }
   
	     IDout <- blast1[rev(grep("RID",blast1))[1] + 2]
	     s1<- as.character('\"')
	     if(grepl(s1,IDout)==TRUE){
	     	
	     	IDout<-strsplit(IDout, s1)[[1]][2]     	
	     }
	     
	     
##############loop to request download until blast finish and get reall result############# 

    
	     
         result <- 0
         name <- paste ("./result/result_",i,".xml",sep = "")
         l<-0
         t<-0
         repeat {
         t<-t+1
         if (l>300) {break}
         if (t>100) {break}
         Sys.sleep(10)
         try<-try(download.file (paste("http://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&RID=",IDout,"&FORMAT_TYPE=XML&FORMAT_OBJECT=Alignment&CMD=Get",sep = "") , name))  
         if('try-error' %in% class(try)){next}
         read<- read.csv(name)
         l<-length(read[[1]])
         }
         
############creat a child folder and extract data need################  

       
         folder <- paste ("./result/blast_",i,sep = "")
         dir.create(folder) 
         
         result <- xmlParse(name, useInternalNodes = TRUE)
         
         BlastOutput_program<- xmlValue(result[["/BlastOutput/BlastOutput_program"]])
         BlastOutput_version<- xmlValue(result[["/BlastOutput/BlastOutput_version"]])
         BlastOutput_reference<- xmlValue(result[["/BlastOutput/BlastOutput_reference"]])
         BlastOutput_db<- xmlValue(result[["/BlastOutput/BlastOutput_db"]])
         BlastOutput_queryID<- xmlValue(result[["/BlastOutput/BlastOutput_query-ID"]])
         BlastOutput_querydef<- xmlValue(result[["/BlastOutput/BlastOutput_query-def"]])
         BlastOutput_querylen<- xmlValue(result[["/BlastOutput/BlastOutput_query-len"]])
         
         data <- result[["/BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits"]]
         numberofalignment <- xmlSize(data)
         dataframe<-xmlToDataFrame(data)
         dataframe<-as.vector(dataframe)
         
         file1<- c("\n***********program**********", BlastOutput_program, "\n**********version**********",BlastOutput_version, "\n**********reference**********",BlastOutput_reference,"\n**********database**********",BlastOutput_db,"\n**********queryID**********",BlastOutput_queryID,"\n**********querydef**********",BlastOutput_querydef,"\n**********querylen**********",BlastOutput_querylen,"\n**********number of alignment**********", numberofalignment)
         name <-paste(folder,"/overview.txt",sep="") 
         
         cat(file1, file=name,sep= "\n")
                          
         name<- paste(folder,"/Hit_id.txt",sep="")
       write.table(paste("\t\t\t",as.vector(unlist(dataframe[2])), "\n", sep=""), sep = "", na = "NA",file = name)

        
         name<- paste(folder,"/Hit_accession.txt",sep="")
         write.table(paste("\t\t\t",as.vector(unlist(dataframe[4])), "\n", sep=""), sep = "", na = "NA",file = name)
         
         name<- paste(folder,"/Hit_len.txt",sep="")
         write.table(paste("\t\t\t",as.vector(unlist(dataframe[5])), "\n", sep=""), sep ="", na = "NA",file = name)
         
         name<- paste(folder,"/Hit_def.txt",sep="")
write.table(paste(as.vector(unlist(dataframe[3])), "\n\n", sep=""), sep ="\n", na = "NA",file = name)         
         Hit_hsps<-c()
         hit<-result[["/BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits"]]
         
         Hit_def<- as.vector(unlist(dataframe[3]))
         
         species<- c()
         Hsp_hseq<-c()
         Hsp_qseq<-c()
         Hsp_midline<-c()
         acce<-as.vector(unlist(dataframe[2]))
         
         for (k in 1:numberofalignment) {
         	
         	split <- c()
         	split <- strsplit(Hit_def[k], "\\[")
            split <- strsplit(split[[1]][length(split[[1]])], "\\]")
            species[k]<-split[[1]][1]
         	
         	Hsp_hseq0<-hit[[k]][["Hit_hsps"]][[1]][["Hsp_hseq"]]
         	Hsp_hseq1<-xmlValue(Hsp_hseq0)
         	Hsp_hseq[k]<-paste("\n>",acce[k]," ",Hit_def[k], "\n",Hsp_hseq1, sep="")
         	
         	Hsp_qseq0 <- hit[[k]][["Hit_hsps"]][[1]][["Hsp_qseq"]]
         	Hsp_qseq1 <- xmlValue(Hsp_qseq0)
         	Hsp_qseq[k]<-paste("\n>",acce[k]," ",Hit_def[k], "\n",Hsp_qseq1, sep="")

         	Hsp_midline0 <- hit[[k]][["Hit_hsps"]][[1]][["Hsp_midline"]]
         	Hsp_midline1 <- xmlValue(Hsp_midline0)
         	Hsp_midline[k]<-paste("\n>",acce[k]," ",Hit_def[k], "\n",Hsp_midline1, sep="")
         	
         }
         
         name<- paste(folder,"/Hit_species.txt",sep="")
         write.table(as.matrix(unlist(species),ncol=1,row.names=c(1:numberofalignment))
, file = name, sep = "\n", na = "NA")

         name<- paste(folder,"/Hsp_hseq.txt",sep="")
         cat(Hsp_hseq, file=name,sep="\n")
         
  	     name<- paste(folder,"/Hsp_qseq.txt",sep="")
         cat(Hsp_qseq, file=name, sep="\n")

         name<- paste(folder,"/Hsp_midline.txt",sep="")
         cat(Hsp_midline, file=name,sep="\n")
         
         name<-paste(folder,"/completeseq.fasta",sep="")
         accesionlist<-paste(as.vector(unlist(dataframe[4])),collapse=";")
         
         urlfasta<-paste("http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&dopt=fasta&db=",database, sep="")
         
         
         
         
##########use extracted accesion number to call for full sequence fasta file##############


         
         t<-0
         repeat{
         t<-t+1
         if(t>10){break}
         fasta<-postForm(urlfasta, val= accesionlist)
         if ('try-error' %in% class(fasta)){
         	Sys.sleep(10)
         	next} else {break}
         }
         cat(fasta, file=name)
                  
##################use fasta file to build phylogenetic tree################################




  	
	     phylogenetree<-paste("./clustalw/clustalw2 -INFILE=./result/blast_", i , "/completeseq.fasta -TYPE=",type2 ," -TREE -QUICKTREE -OUTPUT=PHYLIP -OUTFILE=phylo", sep="")
	     try(system(as.character(phylogenetree)))


#################substract accesion ID to species name, get another tree########################

	     
         treefile<-paste(folder,"/completeseq.ph",sep="")
         read2<-as.vector(unlist(readLines(treefile)))
         read1<- c()
         read3<- c()
         read1<- read2
         ID2<-as.vector(unlist(dataframe[2]))
         species1<-as.vector(unlist(species))
         species2<-c()
         
         for(p in 1:numberofalignment){
             n<-1
             if (length(strsplit(species1[p], " ")[[1]]) > 4) {next}
             while(species1[p]  %in%  species2 == TRUE){
             	species1[p] <- paste(n,as.vector(unlist(species))[p],sep="")
             	n<-n+1
             }
             species2[p]<-species1[p]
             repn<-grep(ID2[p], read2,fixed=TRUE)
             read1[repn]<-sub(ID2[p],species2[p],read2[repn],fixed=TRUE)
         }
         
         cat(read1, file=paste(folder,"/speciesphylotree.ph",sep=""))

#################record number of hit and into new turn#############
         
	     numberofalignmenttable[i] <-numberofalignment
	     numberofalignmentrecord[i+1] <- 	numberofalignment     
	    
	     
	     if (numberofalignment == numberofalignmentrecord[i]) {
	     	similar <- 1
	     }
	     
	     
	     i <- i+1
	     runtime <- runtime-1 
	     cat(paste("\nnumberofalignment=", numberofalignment,"\n\n\n next run#",i,sep=""))
	     
 }
 ##########end of main loop########################################
 
 ##########hit numbers of each term#################################
	     write.table(paste(numberofalignmenttable,"\n\n",sep=""), file="./result/number of alignments.txt")
	     

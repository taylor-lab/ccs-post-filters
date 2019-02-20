library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(tidyr)

specify_decimal = function(x, k) format(round(x, k), nsmall=k)
"%ni%" = Negate("%in%")
curr_date = format(Sys.Date(),"%d-%b-%y")

amaf_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_analysis_mafs'
amaf_files=list.files(amaf_dir,pattern="*muts.maf$",full.names=TRUE)
B = length(unique(amaf_files)); B 

#mapping_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_sample_mapping'
#mapping_files=list.files(mapping_dir,pattern="*sample_mapping.txt$",full.names=TRUE)
#M = length(unique(mapping_files)); M 

##Need to get flags from VarDict so load vcfs
#vcf_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_vcfs'
#vardict_vcf_files = list.files(vcf_dir,pattern="*vardict.vcf$",full.names=TRUE)

#Results directory
batch_ccs_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/cs_deliveries'

#For each batch
for(j in 1:B)
{
  print(j);

  ##Load in the maf file(s) to get the batch id since we are looping by the batch ids 
  maf_file_base_name = str_match(amaf_files[j],"(Proj.*.muts.maf)")[,1]; 
  batchid = str_replace(maf_file_base_name,".muts.maf","");
  batch_out_dir = paste0(batch_ccs_dir,'/res_',batchid)
  
  #Final result maf
  final_ccs_maf = paste0('/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/cs_deliveries/final/',batchid,'.muts.ccs_tagged_123118.maf')

  #if(batchid %in% c("Proj_06058_CFO","Proj_06058_CFB","Proj_set_06058_UT_1"))
  if(!file.exists(final_ccs_maf))
  {
        print(batchid)
        print(final_ccs_maf)
        #nn = paste(j,batchid)
        #write.table(nn,"/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/batch_number_mapping_138.txt",row.names=F,quote=F,sep="\t",append=T)
       
        ##Load in the mapping file(s)
        #batchid_mfile = paste0(batchid,'_sample_mapping.txt')
        #mfile_name = gsub(",", "",toString(mapping_files[as.numeric(grep(batchid_mfile,mapping_files))]))
        #mfile = fread(mfile_name,header=FALSE); head(mfile); dim(mfile)

        ##Load in the analysis maf(s)
        maf0 = fread(amaf_files[j])
        maf = maf0
        print(dim(distinct(maf)))
        #table(maf0$Tumor_Sample_Barcode)
        maf = maf %>% mutate(t_var_freq = t_alt_count/t_depth)
        print(dim(distinct(maf)))
        #########################################################  SNVS
        ##Read in the subsetted vcfs for all samples
        vcffiles <- list.files(path=batch_out_dir,pattern="^VardictVcf_s")
        loadandmerge <- function(text) {
          x <- strsplit(text,split="_")
          name <- paste(x[[1]][2],x[[1]][3],x[[1]][4],x[[1]][5],x[[1]][6],sep="_")
          name = str_replace(name,".txt","")
          file <- paste0(batch_out_dir,"/",text)
          print(file)
          if(!is.na(file.size(file))) {
            aux <- read.table(file,header=F,sep="\t",as.is=T)
            aux$samp <- name
            aux
          } else {
            cat(paste0("Read.table didn't work for ",text,"!\n"))
          }
        }
        print(batch_out_dir)
        pantog <- do.call(rbind,lapply(vcffiles,loadandmerge))
        print(dim(pantog))
        colnames(pantog) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE","NORMAL","samp")
        pantog$id <- paste(pantog$samp,pantog$CHROM,pantog$POS,sep=";")

        ##Separate to get the tags associated with the tumor sample and merge to the original vcf
        pantog_2 <- pantog %>% select(samp,CHROM,POS,id,FORMAT,SAMPLE) %>% separate_rows(FORMAT,SAMPLE, sep = "\\:") %>% spread(FORMAT,SAMPLE);   dim(pantog_2)
        pantog_3 <- pantog %>% select(samp,CHROM,POS,id,INFO) %>% separate_rows(INFO,sep=";") %>% separate(INFO,c("INFO","VALUE"),sep="=") %>% spread(INFO,VALUE); dim(pantog_3)
        pantog2 <- merge(pantog_2,pantog_3,by=c("samp","CHROM","POS","id"))
        dim(pantog2)
        print("SNV's done!!")
        #########################################################  DELETIONS
        
        ##Load these new vcfs for deletions back in
        vcffiles3 <- list.files(path=batch_out_dir,pattern="^VardictVcf_deletions")
        loadandmerge3 <- function(text) {
          x <- strsplit(text,split="_")
          name <- paste(x[[1]][3],x[[1]][4],x[[1]][5],x[[1]][6],x[[1]][7],sep="_")
          name = str_replace(name,".txt","")
          file <- paste0(batch_out_dir,"/",text)
          print(file); 
          if( file.size(file)>0 & !is.na(file.size(file)) ) {
            aux <- read.table(file,header=F,sep="\t",as.is=T)
            aux$samp <- name
            aux
          } else {
            cat(paste0("Read.table didn't work for ",text,"!\n"))
          }
        }
        
        deltog <- do.call(rbind,lapply(vcffiles3,loadandmerge3))
        dim(deltog)
        if(!is.null(deltog))
        {
        colnames(deltog) <- c("CHROM","POSwrong","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE","NORMAL","samp")
        deltog$POS <- deltog$POSwrong + 1
        deltog$id <- paste(deltog$samp,deltog$CHROM,deltog$POS,sep=";"); dim(deltog)
        deltog_2 <- deltog %>% select(samp,CHROM,POS,id,FORMAT, SAMPLE) %>% separate_rows(FORMAT,SAMPLE, sep = '\\:') %>% spread(FORMAT,SAMPLE); dim(deltog_2)
        deltog_3 <- deltog %>% select(samp,CHROM,POS,id,INFO) %>% separate_rows(INFO,sep=";") %>% separate(INFO,c("INFO","VALUE"),sep="=") %>% spread(INFO,VALUE); dim(deltog_3)
        deltog2 <- merge(deltog_2,deltog_3,by=c("samp","CHROM","POS","id")); dim(deltog2)
        
        print("DEL's done!!")
       
        #########################################################
        
        alltog <- rbind(pantog2,deltog2)
        
        print("Combining done!!")
        dim(alltog)
        } else {print ("No deletions!!"); alltog = pantog2}
        #########################################################
        #########################################################  RESCUE MUTECT CALLS

        ##Load the mutect vcfs in
        vcffiles2 <- list.files(path=batch_out_dir,pattern="^MutectVcf")
        loadandmerge2 <- function(text) {
          x <- strsplit(text,split="_")
          name <- paste(x[[1]][3],x[[1]][4],x[[1]][5],x[[1]][6],x[[1]][7],sep="_")
          name = str_replace(name,".txt","")
          file <- paste0(batch_out_dir,"/",text)
          if( file.size(file)>0 & !is.na(file.size(file)) ) {
            aux <- read.table(file,sep="\t",as.is=TRUE,header=FALSE)
            aux$samp <- name
            aux
          } else {
            cat(paste0("Read.table didn't work for ",text,"!\n"))
          }
        }
        mutvcf <- lapply(vcffiles2,loadandmerge2)
        if(!is.null(mutvcf))
        {
          mutvcf2_ <- do.call(rbind,mutvcf)
          #mutvcf2 <- lapply(mutvcf2_,setNames,nm=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE","NORMAL","samp"))
          mutvcf2 <- setnames(mutvcf2_,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE","NORMAL","samp"))
          mutvcf3 <- mutvcf2[mutvcf2$FILTER=="PASS",] ##Take only those that have pass in the filter column, so assuming this is accurate, these should be found in the mutect mafs
          mutvcf3$id <- paste(mutvcf3$samp,mutvcf3$CHROM,mutvcf3$POS,sep=";")
          ##Note whether these sites are found in the mutect maf and would have been "rescued" anyways. As in they would have been found in the maf anyways but the "set" would have just said Mutect instead of Vardict
          dim(mutvcf3)
          print("MUTECT calls done!!")
        } else{ print("NO Mutect calls!!")}
        ##Removed things that aren't in the vcf because these filters only really apply to sites found by Vardict
        maf$id <- paste(maf$Tumor_Sample_Barcode,maf$Chromosome,maf$Start_Position,sep=";")
        dim(maf) 
        names(maf)
        
        
        #########################################################    WHY?!! AND IF NO MUTECT THEN WE CAN SKIP THE NEXT SECTION?
        maf3 <- merge(maf,alltog,by="id") 
        dim(maf3) #671
        print("MERGE with Mutect done!!")
        #########################################################    
        
        ##Not sure, but sometimes there are duplicate rows for certain IDs. some of them the original maf has double lines and sometimes it's in the vcf...not too often
        maf3$LMQ <- ifelse(maf3$TYPE.y!="SNV",ifelse(maf3$NM > 2, "remove",ifelse(maf3$MQ < 55, "remove","ok")),"ok")
        
        ##Need to break down the columns further to get strand bais
        ##Here focus on the strand bias for reads supporting the alternate allele
        maf3$forA <- sapply(strsplit(maf3$ALD,","),function(x){x[1]})
        maf3$revA <- sapply(strsplit(maf3$ALD,","),function(x){x[2]})
        ##If a site (of any type) has 10 reads, then there MUST be at least 1 read in each direction or it will be flagged for removal
        maf3$strandfilt2 <- ifelse(maf3$TYPE.y=="SNV",ifelse(as.numeric(maf3$VD) >= 10,ifelse(maf3$forA==0,"remove",ifelse(maf3$revA==0,"remove","ok")),"ok"),ifelse(as.numeric(maf3$VD) >= 10, ifelse(maf3$forA==0,"remove",ifelse(maf3$revA==0,"remove","ok")),"ok"))
        
        ##Flag those that would likely be "rescued" by Mutect, meaning that it would be in the final maf anyways even in Vardict never found it
        maf3$mutres <- ifelse(maf3$id %in% mutvcf3$id,"rescue","not")
        
        ##How many are removed in that they don't have strand bias, pass the mapping quality filter or are rescued by mutect
        maf3$remove <- ifelse(maf3$mutres=="rescue","ok",ifelse(maf3$strandfilt2=="remove","remove",ifelse(maf3$LMQ=="remove","remove","ok")))
        
        ##Don't want to filter out sites where ALL reads in the tumor go in 1 direction. This generally occurs when the location is at the end of what can be caught by the bait
        maf3$forR <- sapply(strsplit(maf3$RD,","),function(x){x[1]})
        maf3$revR <- sapply(strsplit(maf3$RD,","),function(x){x[2]})
        maf3$fortot <- as.numeric(maf3$forA) + as.numeric(maf3$forR)
        maf3$revtot <- as.numeric(maf3$revA) + as.numeric(maf3$revR)
        maf3$refbias <- ifelse(maf3$fortot==0,"keep",ifelse(maf3$revtot==0,"keep","not"))
        
        ##Keep even if there is strand bias if it's all biased
        maf3$strandfiltnew <- ifelse(maf3$strandfilt2=="remove",ifelse(maf3$refbias=="keep","ok","remove"),maf3$strandfilt2)
        ##Update the remove column
        maf3$removenew <- ifelse(maf3$mutres=="rescue","ok",ifelse(maf3$strandfiltnew=="remove","remove",ifelse(maf3$LMQ=="remove","remove","ok")))
        
        ##Just renaming the variable
        maf4 <- maf3 
        
        print("maf4 is ready!!")
        
        #########################################################        
        ##Read in tumor fillouts
        #########################################################
        
        fillfile = list.files(path=paste0(batch_out_dir,"/fillout/"),pattern="tumor.per_sample_output_fo.maf")
        #Pull the files now
        loadandmerge4 <- function(text) {
          x <- strsplit(text,split="_")
          name <- paste(x[[1]][1],x[[1]][2],x[[1]][3],x[[1]][4],x[[1]][5],sep="_")
          name <- gsub("bam","",name)
          file <- paste0(batch_out_dir,"/fillout/",text)
          if(!is.na(file.size(file))) {
            aux <- fread(file,data.table=FALSE)
            aux$samp <- name
            aux
          } else {
            cat(paste0("Read.table didn't work for ",text,"!\n"))
          }
        }
        fill <- do.call(rbind,lapply(fillfile,loadandmerge4))
        fill$forVAF <- fill$t_alt_count_forward / fill$t_alt_count
        fill$forMAF <- ifelse(fill$forVAF > 0.5, 1-fill$forVAF, fill$forVAF)
        fill$FavStrand <- ifelse(fill$forVAF > 0.5, "For", "Rev")
        fill$Filt <- FALSE
        fill$Filt[is.na(fill$forMAF) | fill$forMAF==0] <- TRUE ##Now marking ones without any reads as TRUE because I'm using this to rescue vardict. And I don't want to rescue something without any supporting reads in Abra...then we just go with whatever VarDict says
        print("Tumor fillout is done!!")
        dim(fill)
        
        #########################################################        
        ##Read in normal fillouts
        #########################################################
        
        Nfillfile = list.files(path=paste0(batch_out_dir,"/fillout/"),pattern="normal.per_sample_output_fo.maf")
        Nfill <- do.call(rbind,lapply(Nfillfile,loadandmerge4))
        Nfill$forVAF <- Nfill$t_alt_count_forward / Nfill$t_alt_count
        Nfill$forMAF <- ifelse(Nfill$forVAF > 0.5, 1-Nfill$forVAF, Nfill$forVAF)
        Nfill$FavStrand <- ifelse(Nfill$forVAF > 0.5, "For", "Rev")
        Nfill$Filt <- FALSE
        Nfill$Filt[is.na(Nfill$forMAF) | Nfill$forMAF==0] <- TRUE ##Don't really use this
        dim(Nfill)
        
        print("Normal fillout is done!!")
        #########################################################        
        
        ##Ok, let's try to combine the files now
        maf4$Ttag <- paste(maf4$Tumor_Sample_Barcode,maf4$Chromosome,maf4$Start_Position,maf4$Reference_Allele,maf4$Tumor_Seq_Allele2,sep=";")
        maf4$Ntag <- paste(maf4$Matched_Norm_Sample_Barcode,maf4$Chromosome,maf4$Start_Position,maf4$Reference_Allele,maf4$Tumor_Seq_Allele2,sep=";")
        print("Tags are ready!!")
        

        fill$Ttag <- paste(sapply(strsplit(as.character(fill$Tumor_Sample_Barcode),"[.]"),function(x){x[1]}),fill$Chromosome,fill$Start_Position,fill$Reference_Allele,fill$Tumor_Seq_Allele1,sep=";")
        fill <- fill[!duplicated(fill$Ttag),]
        rownames(fill) <- fill$Ttag
        dim(fill)
        colnames(fill) <- c(colnames(fill)[1:32],"Tfill_refcount","Tfill_altcount",colnames(fill)[35:37],"Tfill_totcount","Tfill_VF","Tfill_totcountFor","Tfill_refcountFor","Tfill_altcountFor",colnames(fill)[43],"Tfill_forVAF","Tfill_forMAF","Tfill_FavStrand","Tfill_Filt",colnames(fill)[48])


        Nfill$Ntag <- paste(sapply(strsplit(as.character(Nfill$Tumor_Sample_Barcode),"[.]"),function(x){x[1]}),Nfill$Chromosome,Nfill$Start_Position,Nfill$Reference_Allele,Nfill$Tumor_Seq_Allele1,sep=";")
        Nfill <- Nfill[!duplicated(Nfill$Ntag),]
        rownames(fill) <- fill$Ntag
        colnames(Nfill) <- c(colnames(Nfill)[1:32],"Nfill_refcount","Nfill_altcount",colnames(fill)[35:37],"Nfill_totcount","Nfill_VF","Nfill_totcountFor","Nfill_refcountFor","Nfill_altcountFor",colnames(Nfill)[43],"Nfill_forVAF","Nfill_forMAF","Nfill_FavStrand","Nfill_Filt",colnames(Nfill)[48])
        
        print("About to merge fillout columns and maf columns!!")
        maf5.1 <- merge(maf4,fill[,c(33,34,38:42,44:48)],by="Ttag",all.x=TRUE)
        maf5 <- merge(maf5.1,Nfill[,c(33,34,38:42,44:48)],by="Ntag",all.x=TRUE)
        dim(maf5)
        
        ##OK, now toss any site that has >3 reads in the matched normal
        maf5$normreads <- ifelse(maf5$Nfill_altcount > 3,"remove","ok")
        ##Right now, looking for normal reads is only being done on those found in vcfs. This may want to be changed later to include mutect variants as well
        
        ##Tag strand bias properly
        ##strandtog: includes the 10 reads and alternate-supporting strand bias in strandfilt2, rescues strandbiased sites if there is a totalbias as well (as long as the mapping quality is high enough), rescues if it can be found in the mutect vcf or if Abra doesn't think there is strand bias
        maf5$strandtog <- ifelse(maf5$strandfilt2=="ok","ok",
                                 ifelse(maf5$refbias=="keep",ifelse(maf5$MQ > 40,"ok",
                                                                    ifelse(maf5$mutres=="rescue","ok",
                                                                           ifelse(maf5$Tfill_Filt=="TRUE","remove","ok"))),
                                        ifelse(maf5$mutres=="rescue","ok",
                                               ifelse(maf5$Tfill_Filt=="TRUE","remove","ok"))))
        ##refbiasfix: just tagging how many would still be rescued by considering total bias when continue to remove those with low mapping quality
        maf5$refbiasfix <- ifelse(maf5$refbias=="keep",ifelse(maf5$MQ > 40,"keep","not"),"not")
        ##removefinal: put everything together for all sites found in the vardict vcfs. This follows the flowchart
        ##NOTE: If you want to add in those only found in Mutect, change maf5 here to something else
        maf5$removefinal <- ifelse(maf5$set=="VarDict",
                                   ifelse(maf5$hotspot_whitelist=="TRUE","ok",
                                          ifelse(maf5$Nfill_altcount > 3, "remove",
                                                 ifelse(maf5$TYPE.x=="SNV",ifelse(maf5$strandfilt2=="ok","ok",
                                                                                  ifelse(maf5$refbiasfix=="keep","ok",
                                                                                         ifelse(maf5$mutres=="rescue","ok",
                                                                                                ifelse(is.na(maf5$Tfill_forMAF),"remove",
                                                                                                       ifelse(maf5$Tfill_forMAF==0,"remove","ok"))))),
                                                        ifelse(maf5$LMQ=="remove","remove",
                                                               ifelse(maf5$strandfilt2=="ok","ok",
                                                                      ifelse(maf5$refbiasfix=="keep","ok",
                                                                             ifelse(maf5$mutres=="rescue","ok",
                                                                                    ifelse(is.na(maf5$Tfill_forMAF),"remove",
                                                                                           ifelse(maf5$Tfill_forMAF==0,"remove","ok"))))))))),"ok")

        #write.table(maf5,paste0(batch_out_dir,"/",batchid,"_FullFinalMaf.txt"),quote=FALSE,sep="\t",row.names=FALSE)
        table(maf5$refbiasfix, maf5$strandtog)
        maf5 = maf5 %>% mutate(var_alli = str_c(Tumor_Sample_Barcode,';',Chromosome,';',Start_Position),
                                             var_tag = str_c(Tumor_Sample_Barcode,':',Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2))
        #########################################################
        
        mutect_only = subset(maf, set=="MuTect"); dim(mutect_only);
        mutect_only$Ttag <- paste(mutect_only$Tumor_Sample_Barcode,mutect_only$Chromosome,mutect_only$Start_Position,mutect_only$Reference_Allele,mutect_only$Tumor_Seq_Allele2,sep=";")
        mutect_only$Ntag = paste(mutect_only$Matched_Norm_Sample_Barcode,mutect_only$Chromosome,mutect_only$Start_Position,mutect_only$Reference_Allele,mutect_only$Tumor_Seq_Allele2,sep=";")
        mutect_only = merge(mutect_only,Nfill[,c(33,34,38:42,44:48)],by="Ntag",all.x=TRUE)
        mutect_only$normreads = ifelse(mutect_only$Nfill_altcount > 3,"remove","ok")
        mutect_only$removefinal = ifelse(mutect_only$Nfill_altcount > 3,"remove","ok")
        write.table(mutect_only,paste0(batch_out_dir,"/",batchid,"_FullFinalMaf_Mutect_Only.txt"),quote=FALSE,sep="\t",row.names=FALSE)
        print(table(maf5$mutres,maf5$strandfilt2))

        #########################################################          
        ###Backfill NOT NOW BUT RIGHT NOW
        
        head(maf5)
        table(maf5$removefinal)
        mutect_only = mutect_only %>% mutate(var_alli = str_c(Tumor_Sample_Barcode,';',Chromosome,';',Start_Position),
                                             var_tag = str_c(Tumor_Sample_Barcode,':',Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2))
        print("MUTECT Stuff")
        table(mutect_only$removefinal)
        df = rbind.fill(maf5,mutect_only)
        dim(distinct(df)); #696
        df = distinct(df)
        df = df %>% mutate(var_alli = str_c(Tumor_Sample_Barcode,';',Chromosome,';',Start_Position))
        head(df)
        table(df$removefinal)
        
        ##Check what all exists in the original maf's FILTER column
        names(maf0); dim(maf0) #677
        table(maf0$FILTER)
        
        ##Make a var_tag in the original maf as well as the new maf with the removefinal column
        maf_ = distinct(maf0) %>% mutate(var_tag = str_c(Tumor_Sample_Barcode,':',Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2),var_alli = str_c(Tumor_Sample_Barcode,';',Chromosome,';',Start_Position))
        dim(maf_)
        df = distinct(df) %>% mutate(var_tag = str_c(Tumor_Sample_Barcode,':',Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2),
                                     var_alli = str_c(Tumor_Sample_Barcode,';',Chromosome,';',Start_Position))
        dim(df)
        table(df$removefinal)
        table(df$removefinal,df$FILTER)
        
        maf_$removefinal <- ifelse(maf_$var_tag %in% maf5$var_tag[maf5$removefinal=="remove" & !is.na(maf5$removefinal)],"remove",
                                   ifelse(maf_$var_tag %in% mutect_only$var_tag[mutect_only$removefinal=="remove"],"remove","ok"))
        maf_$FILTER <- ifelse(maf_$FILTER=="PASS",ifelse(maf_$removefinal=="remove","ccs_filter",maf_$FILTER),
                                   ifelse(maf_$removefinal=="remove",paste0(maf_$FILTER,",ccs_filter"),maf_$FILTER))
        
        print("OVERALL Stuff")
        print(table(maf_$FILTER))
        dim(maf_)
        write.table(maf_,final_ccs_maf,quote=FALSE,row.names = FALSE,sep="\t",append=FALSE)
        write.table(maf_,paste0(batch_out_dir,"/",batchid,".muts.ccs_tagged_123118.maf"),quote=FALSE,row.names = FALSE,sep="\t",append=FALSE)
        
        #Candidates for IGV review
        #table(maf_$FILTER, maf_$Variant_Type)

  }
}

#DEL DNP INS ONP SNP
#ccs_manual_filter 120   1  43  16  29
#common_variant      2   0   0   0  21
#PASS               28   2  20   0 395
#head(filter(maf_,FILTER=="ccs_manual_filter", Variant_Type == "SNP") %>% select(var_tag, t_alt_count, t_depth))
#filter(maf_,FILTER=="ccs_manual_filter", Variant_Type == "DEL") %>% select(var_tag, t_alt_count, t_depth)
#filter(maf_,FILTER=="ccs_manual_filter", Variant_Type == "INS") %>% select(var_tag, t_alt_count, t_depth)
#filter(maf_,FILTER=="ccs_manual_filter", Variant_Type == "DNP") %>% select(var_tag, t_alt_count, t_depth)

#ls /ifs/res/taylorlab/chavans/roslin_2.4_deliveries/Proj_07058_G/bam/*.bam | grep "s_C_006195_M001_d"
#ls /ifs/res/taylorlab/chavans/roslin_2.4_deliveries/Proj_07058_G/bam/*.bam | grep "s_C_006195_P001_d"

#cmd = 'Rscript /ifs/res/taylorlab/chavans/roslin_2.4_deliveries/AllisonFilters_ApplyFilters_Step3.R'
#cmd_bsub = paste("bsub -e",'/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/cs_deliveries',"-n 1 -R rusage[mem=200] -We 12:00",cmd)

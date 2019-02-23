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
vcf_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_vcfs'
vardict_vcf_files = list.files(vcf_dir,pattern="*vardict.vcf$",full.names=TRUE)
mutect_vcf_files = list.files(vcf_dir,pattern="*mutect.vcf$",full.names=TRUE)

#Results directory
batch_ccs_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/cs_deliveries'

#For each batch
for(j in 1:B)
{
    print(j)
    
    ##Load in the maf file(s) to get the batch id since we are looping by the batch ids 
    maf_file_base_name = str_match(amaf_files[j],"(Proj.*.muts.maf)")[,1]; 
    batchid = str_replace(maf_file_base_name,".muts.maf",""); 
    print(batchid)
    batch_out_dir = paste0(batch_ccs_dir,'/res_',batchid)
    if(!dir.exists(batch_out_dir)) 
          { 
            print(paste0(batch_out_dir," was missing!")) 
          }
          x = paste0('mkdir -p ',batch_out_dir); system(x);
    
    #if( batchid == "Proj_07250_X")
    {
      ##Load in the mapping file(s)
      #batchid_mfile = paste0(batchid,'_sample_mapping.txt')
      #mfile_name = gsub(",", "",toString(mapping_files[as.numeric(grep(batchid_mfile,mapping_files))]))
      #mfile = fread(mfile_name,header=FALSE); head(mfile); dim(mfile)
    
      ##Load in the analysis maf(s)
      maf0 = fread(amaf_files[j])
      maf = maf0
      #table(maf0$Tumor_Sample_Barcode)
      maf = maf %>% mutate(t_var_freq = t_alt_count/t_depth)
    
      vcfs = list.files(batch_out_dir,pattern="^MutectVcf_VardictOverlap*",full.names=TRUE)
    
      if(length(vcfs) == 0)
      {
        #########################################################  SNVS
        
        bashvcf = do.call(rbind,lapply(unique(maf$Tumor_Sample_Barcode),function(x){
          file = gsub(",", "",toString(vardict_vcf_files[as.numeric(grep(x,vardict_vcf_files))]))
          if(is.na(file)){print(paste0("!!!VCF not found for ",x))}
          mafsearch = paste(unique(maf$Start_Position[maf$Tumor_Sample_Barcode==x]),collapse=" -e ")
        
          cmd = paste0("'grep -e ",mafsearch," ",file," > ",batch_out_dir,"/","VardictVcf_snvs_",x,".txt'")
          cmd_bsub = paste("bsub -e",batch_out_dir,"-n 1 -R rusage[mem=10] -We 1:57",cmd)
          write.table(cmd,paste0(batch_out_dir,"/",batchid,".sh"),quote=FALSE,row.names = FALSE,sep="\t",append=FALSE)
          write.table(cmd_bsub,paste0(batch_out_dir,"/",batchid,".bsub.sh"),quote=FALSE,row.names = FALSE,sep="\t",append=FALSE)
          system(cmd_bsub)
          }))
    
        #########################################################  DELETIONS
    
        ##Grab the deletions. This is because deletions are coded differently between mafs and vcfs
        bashvcf2 = do.call(rbind,lapply(unique(maf$Tumor_Sample_Barcode),function(x){
          file = gsub(",", "",toString(vardict_vcf_files[as.numeric(grep(x,vardict_vcf_files))]))
          if(is.na(file)){print(paste0("!!!VCF not found for ",x))}
          
          mafsearch = paste(unique(maf$Start_Position[maf$Tumor_Sample_Barcode==x & maf$Variant_Type=="DEL"])-1,collapse=" -e ") ##subtract 1 so it can find it
          flag = ifelse(mafsearch == "",0,1); 
            if(flag==1)
              {
              cmd2 = paste0("'grep -e ",mafsearch," ",file," > ",batch_out_dir,"/","VardictVcf_deletions_",x,".txt'")
              cmd2_bsub = paste("bsub -e",batch_out_dir,"-n 1 -R rusage[mem=10] -We 1:57",cmd2)
              write.table(cmd2,paste0(batch_out_dir,"/",batchid,".sh"),quote=FALSE,row.names = FALSE,sep="\t",append=FALSE)
              write.table(cmd2_bsub,paste0(batch_out_dir,"/",batchid,".bsub.sh"),quote=FALSE,row.names = FALSE,sep="\t",append=FALSE)
              system(cmd2_bsub)
              }
          }))
    
        #########################################################  RESCUE MUTECT CALLS
        
        ##Grab the vcfs for mutect to see how many of these filtered out ones may also be found in mutect ( how many of the SNVs would be rescued by mutect?)
        ##CAUTION: THIS IS USING THE VCF AND THESE MAY BE FILTERED OUT BY MUTECT BY THE TIME THEY GET TO THE MAF
        
        readmut <- do.call(rbind,lapply(unique(maf$Tumor_Sample_Barcode),function(x){
          cat(paste0(x,"\n"))
          file = gsub(",", "",toString(mutect_vcf_files[as.numeric(grep(x,mutect_vcf_files))]))
          if(is.na(file)){print(paste0("!!!VCF not found for ",x))}
          mafsearch <- paste(unique(maf$Start_Position[maf$Tumor_Sample_Barcode==x & maf$TYPE=="SNV"]),collapse=" -e ")
          cmd3 = paste0("'grep -e ",mafsearch," ",file," > ",batch_out_dir,"/","MutectVcf_VardictOverlap_",x,".txt'")
          cmd3_bsub = paste("bsub -e",batch_out_dir,"-n 1 -R rusage[mem=10] -We 1:57",cmd3)
          write.table(cmd3,paste0(batch_out_dir,"/",batchid,".sh"),quote=FALSE,row.names = FALSE,sep="\t",append=FALSE)
          write.table(cmd3_bsub,paste0(batch_out_dir,"/",batchid,".bsub.sh"),quote=FALSE,row.names = FALSE,sep="\t",append=FALSE)
          system(cmd3_bsub)     
          }))
      } 
    }   
}

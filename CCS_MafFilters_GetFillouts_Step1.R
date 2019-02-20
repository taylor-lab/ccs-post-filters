#libraries & custom functions
suppressWarnings(library(data.table))
suppressWarnings(library(plyr))
suppressWarnings(library(dplyr))
suppressWarnings(library(stringr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))

specify_decimal = function(x, k) format(round(x, k), nsmall=k)
"%ni%" = Negate("%in%")
curr_date = format(Sys.Date(),"%d-%b-%y")

bam_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_bams'
bam_files=list.files(bam_dir,pattern="*.bam$",full.names=TRUE)
print(length(unique(bam_files)))

sample_pairing_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_sample_pairing'
sample_pairing_files=list.files(sample_pairing_dir,pattern="*.txt$",full.names=TRUE)
print(length(unique(sample_pairing_files))) 

#N_bam_files = list.files(bam_dir,pattern=".*N00.*.bam$",full.names=TRUE) #temp exception for the 6058 batches

amaf_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_analysis_mafs'
amaf_files=list.files(amaf_dir,pattern="*muts.maf$",full.names=TRUE)
B = length(unique(amaf_files)); B 

#For each batch
for(j in 1:B)
{
  
  #Patient level maf_fillout across exomes of same patient
  big_exome_maf=fread(amaf_files[j],skip = 1 )
  N = length(unique(substring(big_exome_maf$Tumor_Sample_Barcode,1,10))); N
  pat_ids = unique(substring(big_exome_maf$Tumor_Sample_Barcode,1,10))
  S = length(unique(big_exome_maf$Tumor_Sample_Barcode)); S
  sample_ids = unique(big_exome_maf$Tumor_Sample_Barcode); sample_ids
  maf_file_base_name = str_match(amaf_files[j],"(Proj.*.muts.maf)")[,1]; 
  batchid = str_replace(maf_file_base_name,".muts.maf","");
  
  	#if(batchid %in% c("Proj_00940"))
  	{
    
    	fo_dir = paste0('/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/cs_deliveries/res_',batchid,'/fillout')
    	curr_amaf_TNpairs = fread(amaf_files[j],skip = 1) %>% select(Tumor_Sample_Barcode,Matched_Norm_Sample_Barcode) %>% distinct(.)

        if(!dir.exists(fo_dir)) 
        	{ 
        		print(paste0(fo_dir,"was missing!")) 
        	}
        	x = paste0('mkdir -p ',fo_dir); system(x);
        	#print(fo_dir)
      
        	fo_cmds_N = paste0(fo_dir,'/fo_cmds_N.txt')
        	fo_cmds_T = paste0(fo_dir,'/fo_cmds_T.txt')
  
            #For each tumor sample
            for(i in 1:S)
            {  
	              pat = substring(curr_amaf_TNpairs$Tumor_Sample_Barcode[i],1,10) #10 normally, 12 for exception batches
	              pat_maf = filter(big_exome_maf, grepl(pat,Tumor_Sample_Barcode))
	              print(pat)
	              print(dim(pat_maf))
	        
	              pat_maf_in = paste0(fo_dir,'/',pat,'.per_patient_input_fo.maf')
	              print(pat_maf_in)      
	              write.table(pat_maf, pat_maf_in, sep="\t",row.names=FALSE,append=FALSE,quote=FALSE)
	            
	              T_sample_bam = gsub(",", "",toString(bam_files[as.numeric(grep(curr_amaf_TNpairs$Tumor_Sample_Barcode[i],bam_files))]))
	              N_sample_bam = gsub(",", "",toString(bam_files[as.numeric(grep(curr_amaf_TNpairs$Matched_Norm_Sample_Barcode[i],bam_files))]))
	              
	              #Setting filenames
	              T_sample_maf_out = paste0(fo_dir,'/',curr_amaf_TNpairs$Tumor_Sample_Barcode[i],'.tumor.per_sample_output_fo.maf')
	              N_sample_maf_out = paste0(fo_dir,'/',curr_amaf_TNpairs$Tumor_Sample_Barcode[i],'.normal.per_sample_output_fo.maf')
	              
	              print(curr_amaf_TNpairs$Tumor_Sample_Barcode[i])
	              T_sample_fo = paste('/home/chavans/git/ngs-filters/maf_fillout.py','-m',pat_maf_in,'-b',T_sample_bam,'-g GRCh37 -o',T_sample_maf_out,'-mo')
	              N_sample_fo = paste('/home/chavans/git/ngs-filters/maf_fillout.py','-m',pat_maf_in,'-b',N_sample_bam,'-g GRCh37 -o',N_sample_maf_out,'-mo')
	              
	              T_fo_bsub = paste("bsub -e",fo_dir,"-n 1 -R rusage[mem=10] -We 1:57",T_sample_fo)
	              N_fo_bsub = paste("bsub -e",fo_dir,"-n 1 -R rusage[mem=10] -We 1:57",N_sample_fo)
	        
	              system(T_fo_bsub)
	              system(N_fo_bsub)
	        
	              write.table(T_fo_bsub,fo_cmds_T,sep="\t",row.names=FALSE,append=TRUE,quote=FALSE)
	              write.table(N_fo_bsub,fo_cmds_N,sep="\t",row.names=FALSE,append=TRUE,quote=FALSE)
            }
    }
}


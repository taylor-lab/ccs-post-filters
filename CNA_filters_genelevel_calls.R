library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(jsonlite)
specify_decimal = function(x, k) format(round(x, k), nsmall=k)
"%ni%" = Negate("%in%")

sample_pairing = fread('/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_sample_pairing/master_sample_pairing_123118.txt')
names(sample_pairing)=c("N","T")
sample_pairing<-filter(sample_pairing,N!="na" & T!="na")
total_tumors=dim(sample_pairing)[1]
print(paste0(total_tumors," samples"))

amaf_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_analysis_mafs'
amaf_files=list.files(amaf_dir,pattern="*muts.maf$",full.names=TRUE)
B = length(unique(amaf_files)); B 

facets_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_facets'
fc_lu_table0 = fread("/ifs/res/taylorlab/richara4/gitstuff/facets-suite/FACETS_CALL_table.tsv")
fc_lu_table = fc_lu_table0 %>% mutate(emtag = str_c(WGD,';',mcn,';',lcn))
WGDcut = 0.5

oncokb = fromJSON(readLines('http://oncokb.org/api/v1/genes', warn=F)); head(oncokb)
oncokb_tsg = filter(oncokb, tsg=="TRUE") %>% select(hugoSymbol) %>% distinct(.)

#For each batch
for(j in 1:1)
{
  ##Load in the maf file(s) to get the batch id since we are looping by the batch ids 
  maf_file_base_name = str_match(amaf_files[j],"(Proj.*.muts.maf)")[,1]; 
  batchid = str_replace(maf_file_base_name,".muts.maf","");
  batch_out_dir = paste0(batch_ccs_dir,'/res_',batchid)
  curr_amaf = fread(amaf_files[j]) %>% select(Tumor_Sample_Barcode,Matched_Norm_Sample_Barcode) %>% distinct(.)
  print(curr_amaf)
  S = length(unique(curr_amaf$Tumor_Sample_Barcode))
  print(batchid)
  
  ### Step1 Extract EM-based purity calls (from the selected fit purity.out file) -- per sample
  for(i in 1:S)
  {  
    sample_dir_name = paste0(curr_amaf$Tumor_Sample_Barcode[i],'__',curr_amaf$Matched_Norm_Sample_Barcode[i])
    fit_dir = 'facets_R0.5.6c100p500'
    facets_sample_dir = paste0(facets_dir,"/",sample_dir_name,"/",fit_dir)
    facets_out_file = list.files(facets_sample_dir,pattern="*_purity.out",full.name=TRUE)
    purity___ = grep("Purity", readLines(facets_out_file), value = TRUE); 
    purity__ = gsub(" ","",purity___); purity_ = gsub("#Purity=","",purity__); 
    purity = purity_[1]; curr_amaf$PurityEM[i] = purity
    curr_amaf$CFcut[i] <- 0.6 * as.numeric(purity)
  }
  print(curr_amaf)
  
  ## To handle purity = NA's in Step 3
  curr_amaf = curr_amaf %>% mutate(CFcut = ifelse(is.na(curr_amaf$CFcut),10,CFcut)) 
  
  ### Step 2.1 CONVERT GENELEVEL CALLS CNCF-based to EM-based -- per batch
  genelevelcalls_ = fread(paste0('/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/',batchid,'/analysis/',batchid,'.gene.cna.txt')) 
  genelevelcalls0 =  genelevelcalls_ %>% 
                     mutate(segid = str_c(Tumor_Sample_Barcode,';',chr,';',seg.start,';',seg.end)) %>%
                     mutate(mcn.em = tcn.em - lcn.em, seg.len = seg.end - seg.start) 

  u.genelevelcalls = genelevelcalls0 %>% 
                   select(Tumor_Sample_Barcode,tcn,lcn,tcn.em,lcn.em,chr,seg.start,seg.end,frac_elev_major_cn,WGD,mcn,FACETS_CNA,FACETS_CALL) %>% distinct(.) %>%
                   mutate(mcn.em = tcn.em - lcn.em, seg.len = seg.end - seg.start)
  
  frac_elev_major_cn.em = u.genelevelcalls %>%  
                          group_by(Tumor_Sample_Barcode) %>% 
                          summarize(frac_elev_major_cn.em = sum(as.numeric(mcn.em >= 2)*as.numeric(seg.len), na.rm = T) / sum(as.numeric(seg.len)))
  
  genelevelcalls0 = left_join(genelevelcalls0,frac_elev_major_cn.em, by="Tumor_Sample_Barcode")
  
          ### Step 2.2 Going back to original un-unique genelevel calls, just carry forward frac_elev_major_cn.em
          genelevelcalls0 = genelevelcalls0 %>% 
                            mutate(WGD.em = ifelse(frac_elev_major_cn.em > WGDcut, "WGD", "no WGD")) %>%
                            mutate(emtag = str_c(WGD.em,';',mcn.em,';',lcn.em))
          
          genelevelcalls0 = genelevelcalls0 %>% 
                            mutate(FACETS_CALL.em = plyr::mapvalues(emtag, fc_lu_table$emtag, fc_lu_table$FACETS_CALL)) %>% ##ASK
                            mutate(FACETS_CALL.em = ifelse(tcn.em >= 6,"AMP",FACETS_CALL.em)) %>%
                            mutate(FACETS_CALL.em = ifelse( (!is.na(tcn.em) & !is.na(lcn.em) & FACETS_CALL.em %ni% c(unique(fc_lu_table$FACETS_CALL)) ), "ILLOGICAL", FACETS_CALL.em))
          table(genelevelcalls0$FACETS_CALL.em)
          
          seg.count = plyr::count(genelevelcalls0, vars="segid") %>% 
                            mutate(count = freq) %>% select(-c(freq))
          
          genelevelcalls0 = inner_join(genelevelcalls0,seg.count,by='segid')
  
  ### Step 3 SET columns needed for filters & APPLY filters
          genelevelcalls0 = genelevelcalls0 %>% 
                    mutate(CFcut = plyr::mapvalues(Tumor_Sample_Barcode, curr_amaf$Tumor_Sample_Barcode, curr_amaf$CFcut))
  
        ### CHECK FOR *AMP/AMP (LOH)* 1. seg size <= 10 Mb AND 2. tcn.em >8 OR gene.count <= 10 OR cf.em >= 0.6*Purity
        ### CHECK FOR *HOMDEL* 1. seg size <= 10 Mb AND 2. gene.count <= 10 
        ### SUPRESS ALL OTHER CALLS FACETS_CNA = "SuppressedCall"
          genelevelcalls0 = genelevelcalls0 %>%
                    mutate(FACETS_CALL.ori = FACETS_CALL.em, 
                           FACETS_CALL.em = ifelse( FACETS_CALL.em %in% c("AMP","AMP (LOH)","HOMDEL"), 
                                              ifelse( (FACETS_CALL.em %in% c("AMP","AMP (LOH)") & seg.len < 10000000 & (tcn.em > 8 | count <=10 | cf.em > CFcut )), FACETS_CALL.em, 
                                                ifelse( (FACETS_CALL.em == "HOMDEL" & seg.len < 10000000 & count <= 10), FACETS_CALL.em, "ccs_filter")), FACETS_CALL.em )) 
          table(genelevelcalls0$FACETS_CALL.em)
          write.table(genelevelcalls3,paste0(batch_out_dir,"/",batchid,"_ccs_filtered_genelevel_cna.txt"),row.names=FALSE,quote=FALSE,sep="\t")
          
          homdeltsg_review = filter(genelevelcalls0, FACETS_CALL.em == "ccs_filter", FACETS_CALL.ori == "HOMDEL", Hugo_Symbol %in% unique(oncokb_tsg$hugoSymbol), seg.len < 25000000)
       
  table(homdeltsg_review$FACETS_CALL.ori)
  write.table(homdeltsg_review,paste0(batch_out_dir,"/",batchid,"_ccs_manual_review_genelevel_candidate_tsg_homdels.txt"),row.names=FALSE,quote=FALSE,sep="\t")

}

  
  
  
  
  
  

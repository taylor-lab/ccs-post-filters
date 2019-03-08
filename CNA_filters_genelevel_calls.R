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

pairing_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_sample_pairing'
pairing_files=list.files(pairing_dir,pattern="*sample_pairing.txt$",full.names=TRUE)
P = length(unique(pairing_files)); P 

facets_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_facets'
batch_ccs_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/cs_deliveries'

fc_lu_table0 = fread("/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/FACETS_CALL_table.tsv")
fc_lu_table = fc_lu_table0 %>% mutate(emtag = str_c(WGD,';',mcn,';',lcn))
WGDcut = 0.5

oncokb = fromJSON(readLines('http://oncokb.org/api/v1/genes', warn=F)); head(oncokb)
oncokb_tsg = filter(oncokb, tsg=="TRUE") %>% select(hugoSymbol) %>% distinct(.)

#For each batch
for(j in 1:P)
{
  ##Load in the maf file(s) to get the batch id since we are looping by the batch ids 
    pairing_file_base_name = str_match(pairing_files[j],"(Proj.*_sample_pairing.txt)")[,1]; 
    batchid = str_replace(pairing_file_base_name,"_sample_pairing.txt",""); 
    print(batchid)
    print(j)

    sample_pairing = fread(pairing_files[j],header=FALSE); head(sample_pairing); dim(sample_pairing)
    names(sample_pairing)=c("N","T")
    sample_pairing<-filter(sample_pairing,N!="na" & T!="na")
    total_tumors=dim(sample_pairing)[1]
    print(paste0(total_tumors," samples"))

    batch_out_dir = paste0(batch_ccs_dir,'/res_',batchid)
    batch_raw_cna = paste0(batch_out_dir,'/',batchid,'_genelevel_cna.txt')
  
  ### Step1 Extract EM-based purity calls (from the selected fit purity.out file) -- per sample
  for(i in 1:total_tumors)
  {  
    sample_dir_name = paste0(sample_pairing$T[i],'__',sample_pairing$N[i])
    fit_dir = 'facets_R0.5.6c100p500'
    facets_sample_dir = paste0(facets_dir,"/",sample_dir_name,"/",fit_dir)
    facets_out_file = list.files(facets_sample_dir,pattern="*_purity.out",full.name=TRUE)
    print(facets_out_file)
    purity___ = grep("Purity", readLines(facets_out_file[1]), value = TRUE); 
    purity__ = gsub(" ","",purity___); 
    purity_ = gsub("#Purity=","",purity__); 
    purity = purity_[1]; 
    sample_pairing$PurityEM[i] = purity
    print(purity)
    sample_pairing$CFcut[i] <- 0.6 * as.numeric(purity)
  }
  print(sample_pairing)
  
  ## To handle purity = NA's in Step 3
  sample_pairing = sample_pairing %>% mutate(CFcut = ifelse(is.na(sample_pairing$CFcut),10,CFcut)) 
  
  ### Step 2.1 CONVERT GENELEVEL CALLS CNCF-based to EM-based -- per batch
    #curr_genelevel = genelevel_files[as.numeric(grep(batchid,genelevel_files))] 
    #genelevelcalls_ = fread(paste0('/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/',batchid,'/analysis/',batchid,'.gene.cna.txt')) 
  curr_genelevel = fread(batch_raw_cna)
  genelevelcalls0 =  curr_genelevel %>% 
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
                    mutate(CFcut = plyr::mapvalues(Tumor_Sample_Barcode, sample_pairing$T, sample_pairing$CFcut))
  
          genelevelcalls0 = genelevelcalls0 %>%
                    mutate(FACETS_CALL.ori = FACETS_CALL.em, 
                           FACETS_CALL.em = ifelse( FACETS_CALL.em %in% c("AMP","AMP (LOH)","AMP (BALANCED)","HOMDEL"), 
                                              ifelse( (FACETS_CALL.em %in% c("AMP","AMP (LOH)","AMP (BALANCED)") & seg.len < 10000000 & (tcn.em > 8 | count <=10 | cf.em > CFcut )), FACETS_CALL.em, 
                                                ifelse( (FACETS_CALL.em == "HOMDEL" & seg.len < 10000000 & count <= 10), FACETS_CALL.em, "ccs_filter")), FACETS_CALL.em )) 
          table(genelevelcalls0$FACETS_CALL.em)
          homdeltsg_review = filter(genelevelcalls0, FACETS_CALL.em == "ccs_filter", FACETS_CALL.ori == "HOMDEL", Hugo_Symbol %in% unique(oncokb_tsg$hugoSymbol), seg.len < 25000000)
          
          genelevelcalls0 = genelevelcalls0 %>% 
                            mutate(FACETS_CNA.em = plyr::mapvalues(FACETS_CALL.em, fc_lu_table$FACETS_CALL, fc_lu_table$FACETS_CNA)) %>%
                            mutate(FACETS_CNA.em = ifelse(FACETS_CALL.em=="ccs_filter",0,FACETS_CNA.em)) %>%
                            mutate(FACETS_CNA.em = ifelse(FACETS_CALL.em=="ILLOGICAL",NA,FACETS_CNA.em)) %>%
                            select(-c(FACETS_CNA, FACETS_CALL, CFcut, FACETS_CALL.ori, count, segid, seg.len, emtag, WGD.em, mcn, mcn.em, frac_elev_major_cn.em ))
          
          genelevelcalls_final = genelevelcalls0 %>% 
                            mutate(FACETS_CNA = FACETS_CNA.em, FACETS_CALL = FACETS_CALL.em) %>%
                            select(-c(FACETS_CALL.em,FACETS_CNA.em))
          
          table(genelevelcalls_final$FACETS_CNA, genelevelcalls_final$FACETS_CALL)       
          head(genelevelcalls_final)
          write.table(genelevelcalls_final,paste0(batch_out_dir,"/",batchid,"_ccs_filtered_genelevel_cna.txt"),row.names=FALSE,quote=FALSE,sep="\t")

  table(homdeltsg_review$FACETS_CALL.ori)
  write.table(homdeltsg_review,paste0(batch_out_dir,"/",batchid,"_ccs_homdeltsg_review_candidates.txt"),row.names=FALSE,quote=FALSE,sep="\t")

}

## Compare with IMPACT
  
  
  
  
  
  

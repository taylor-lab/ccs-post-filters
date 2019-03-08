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

facets_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_facets'
#cmo_facets = '/opt/common/CentOS_6-dev/python/python-2.7.10/bin/cmo_facets'
batch_ccs_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/cs_deliveries'

target_file = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/Homo_sapiens.GRCh37.75.canonical_exons.bed.ilist.facets_format_specific.CORRECTED'

pairing_dir = '/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_sample_pairing'
pairing_files=list.files(pairing_dir,pattern="*sample_pairing.txt$",full.names=TRUE)
P = length(unique(pairing_files)); P 

for(j in 1:P)
{
		#j = 1
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

		list_cncf_batch = NA
		for(i in 1:total_tumors)
		{
			sample_dir_name = paste0(facets_dir,'/',sample_pairing$T[i],'__',sample_pairing$N[i])
			sample_hisens_cncf = list.files(sample_dir_name, pattern="*_hisens.cncf.txt$",full.names=TRUE, recursive = TRUE)
			list_cncf_batch[i] = sample_hisens_cncf
		}		

			list_cncf_batch = str_replace_all(toString(list_cncf_batch),",","")
			print(list_cncf_batch)

			cmd = paste('cmo_facets --suite-version 1.5.5 geneLevel -f',list_cncf_batch,'-o',batch_raw_cna,'-t', target_file)
			system(cmd)   
			perl_cmd = paste("perl -i -pe  's/^.*all_facets\\/(.*)__.*facets.*txt/$1/g' "," ",batch_raw_cna)
			system(perl_cmd)
			#cmd_bsub = paste("bsub -e",batch_out_dir,"-n 1 -R rusage[mem=10] -We 1:57",cmd)	
}
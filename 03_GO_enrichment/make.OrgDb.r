library(AnnotationForge) #加载
library(GO.db)
args <- commandArgs(T)
egg<-read.csv(args[1],header=T,sep="\t") 
egg[egg==""] <- NA

###提取注释的GO信息###
library(dplyr)
library(stringr)
gene_info <- egg %>%dplyr::select(GID = query, GENENAME = seed_ortholog) %>% na.omit() #根据egg第一和第二列的标题提取前两列。
goterms <- egg %>%dplyr::select(query, GOs) %>% na.omit() %>% filter(str_detect(GOs,"GO"))
all_go_list=str_split(goterms$GOs,",") #分隔goterms的GOs值，逗号是分隔符
gene2go <- data.frame(GID = rep(goterms$query, times = sapply(all_go_list, length)), GO = unlist(all_go_list), EVIDENCE = "IEA") %>% filter(str_detect(GO,"GO"))

###提取注释的KEGG信息###
koterms <- egg %>%dplyr::select(GID = query, KO=KEGG_ko)%>%na.omit()%>% filter(str_detect(KO,"ko")) #根据egg第一列和KEGG_ko列的标题提取基因的KEGG注释。

#load("/dellfsqd2/ST_OCEAN/USER/liyanan2/03_job/single_cell_seq_analysis/make.OrgDb/kegg_info.RData")
load(args[2])
colnames(ko2pathway)=c("KO",'Pathway')
koterms$KO=str_replace_all(koterms$KO,"ko:","") #把koterms的KO值的前缀ko:去掉，与ko2pathway格式一致。
gene2pathway <- koterms %>% left_join(ko2pathway, by = "KO") %>%dplyr::select(GID, Pathway) %>%na.omit() #合并koterms和ko2pathway到gene2pathway，将基因与pathway的对应关系整理出来
gene2pathway_name<-left_join(gene2pathway,pathway2name,by="Pathway") #合并gene2pathway和pathway2name
write.table(gene2pathway_name,file="gene2pathway_name.txt",sep="\t",row.names=F,quote=F) 

makeOrgPackage(gene_info=gene_info, go=gene2go, ko=koterms,  pathway=gene2pathway, version="0.0.1", maintainer='liyanan2', author='liyanan2',outputDir=".", tax_id=args[3], genus=args[4], species=args[5],goTable="go")



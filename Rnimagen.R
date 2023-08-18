library(tidyverse)
library(data.table)
library(plotly)
library(tidytext)
library(dplyr)
# # everything with original read
# mydata_in <- fread("/home/snm0205/nimagen/HIS/29962_22_results.txt", sep = "\t", header = F)
# mydata_in <- fread("/home/snm0205/nimagen/8898.results_03242022.txt", sep = "\t", header = F)
# mydata_in <- fread("/home/snm0205/nimagen/112_CAU.results_03242022.txt", sep = "\t", header = F)
# mydata_in <- fread("/home/snm0205/nimagen/112_CAU_results_04082022.txt", sep = "\t", header = F)
# mydata_in <- fread("/home/snm0205/nimagen/2486_AFA_results_04082022.txt", sep = "\t", header = F)
# setnames(mydata_in, c("leftp", "leftp_seq",
#                       "leftp_ham", "tseq",
#                       "rightp", "rightp_seq", 
#                       "rightp_ham", "readcount"))
# 
# #write_delim(mydata_in, file = "8898.results.R1_1.txt", delim = "\t", col_names = F)
# 
# # mydata_in_1 is everything except original read but derived original read by paste0
# mydata_in %>% mutate(lentseq = nchar(tseq),
#                      origRead = paste0(leftp_seq, tseq, rightp_seq ), 
#                      lenOrigRead = nchar(origRead))-> mydata_in_1
# 
# mydata_in_1 %>% distinct(rightp) %>% count()
# 
# primer_in <- fread("Updated_IDseek_SNP85_Primer_Seq_hg38_03232021_SM020922.txt", sep = "\t")
# 
# primer_in %>% anti_join(mydata_in_1, by=c("Left primer"="leftp")) -> foo
# 
# # just complete read
# mydataR <- fread("8898.results.R1_2.txt", sep = "\t", header = F, skip = 2,
#                  fill = T)
# 
# #everything except read
# mydata_exR <- fread("8898.results.R1.txt", sep = "\t", header = F)
# 
# 
# mydata_in_1 %>% filter(leftp_seq=="ATCCTTTCCTGCTTTTCCTCCTCCC", 
#                        rightp_seq=="GGGGGCCTGTCCTTTACTCACTATA"
#                        )
# 
# mydata_in %>% filter(leftp_seq=="ATCCTTTCCTGCTTTTCCTCCTCCC", 
#                      rightp_seq=="GGGGGCCTGTCCTTTACTCACTATA"
#                      )
# 
# 
# # just leftpsseq, tseq, and rightpseq
# ltr <- fread("8898.results.R1_3.txt", sep = "\t", header = F)
# 
# #corrected code
# mydata_correct <- fread("8898.results.R1_4.txt", sep = "\t", header = F, skip = 2,
#                         fill = T)


readFile <- function(x){
  tib <- fread(x, sep = "\t", header = F, skip = 2)
  colnames(tib) <- c("leftp", "leftp_seq",
                   "leftp_ham", "tseq",
                   "rightp", "rightp_seq", 
                   "rightp_ham", "readcount")
  
  tib$Filename <- x
  return(tib)
}


inargs<-commandArgs(trailingOnly = TRUE)
#print(inargs)
cau11<-"CAU/11/*results.txt"
cau22<-"CAU/22/*results.txt"
cau33<-"CAU/33/*results.txt"
afa11<-"AFA/11/*results.txt"
afa22<-"AFA/22/*results.txt"
afa33<-"AFA/33/*results.txt"
his11<-"HIS/11/*results.txt"
his22<-"HIS/22/*results.txt"
his33<-"HIS/33/*results.txt"

filescau11<- Sys.glob(cau11)
filescau22<- Sys.glob(cau22)
filescau33<- Sys.glob(cau33)
filesafa11<- Sys.glob(afa11)
filesafa22<- Sys.glob(afa22)
filesafa33<- Sys.glob(afa33)
fileshis11<- Sys.glob(his11)
fileshis22<- Sys.glob(his22)
fileshis33<- Sys.glob(his33)

all2gethercau11 <- rbindlist(
  map(filescau11, readFile)
)

all2gethercau22 <- rbindlist(
  map(filescau22, readFile)
)

all2gethercau33 <- rbindlist(
  map(filescau33, readFile)
)

all2getherafa11 <- rbindlist(
  map(filesafa11, readFile)
)

all2getherafa22 <- rbindlist(
  map(filesafa22, readFile)
)


all2getherafa33 <- rbindlist(
  map(filesafa33, readFile)
)

all2getherhis11 <- rbindlist(
  map(fileshis11, readFile)
)

all2getherhis22 <- rbindlist(
  map(fileshis22, readFile)
)

all2getherhis33 <- rbindlist(
  map(fileshis33, readFile)
)


all2gether <- bind_rows(all2gethercau11, all2gethercau22, all2gethercau33,
                        all2getherafa11, all2getherafa22, all2getherafa33,
                        all2getherhis11, all2getherhis22, all2getherhis33
                        )



all2gether %>% separate(Filename, into = c("Pop","Fuzz","Sample"),sep = "/") %>%
  separate(Sample, into = c("Sample", NA), sep = "_") %>% 
  #filter(Sample==112) -> foo
  group_by(Fuzz, Pop, Sample, leftp, tseq) %>% 
  summarise(ReadCTotal=sum(readcount)) %>%
  arrange(desc(ReadCTotal), .by_group = TRUE) -> foo1

# there is warning message about additional pieces discarded because for fuzzes 22 and 33
# the sample names are <samplename>_<fuzz>_reuslts.txt

# read in the primer data to get the name of the locus and left join with dataframe
# above to properly label x-axis

primer_in <- fread("Updated_IDseek_SNP85_Primer_Seq_hg38_03232021_SM04082022_primerMod.txt")
primer_in %>% mutate(leftp = ifelse(New_switch==1,`Left primer`,`Right primer`))-> foo

foo1 %>% left_join(foo, by = "leftp") -> foo2

foo2 %>%
  mutate(rank = rank(-ReadCTotal, ties.method = "first")) %>%
  top_n(-2, rank) %>%
  mutate(logReadCTotal = log10(ReadCTotal)) -> foo2.1

foo1 %>% filter(Fuzz==11) %>% group_by(leftp) %>% count() -> bar2

#faceting across fuzzes
foo2 %>% 
  mutate(Fuzz=case_when(Fuzz=="11"~1,
                               Fuzz=="22"~2,
                               Fuzz=="33"~3)) %>%
  ggplot(aes(x=reorder(RefSNP,ReadCTotal),y=ReadCTotal)) +
  geom_boxplot(color="deepskyblue4")+
  facet_grid(Fuzz~.)+
  labs(y="Read count",
       x="SNPs") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



#Looking at all 85 markers across all fuzzs and pops
foo2.1 %>% filter(Fuzz==11) %>%
  ggplot(aes(x=reorder(RefSNP,ReadCTotal),y=ReadCTotal)) + 
  geom_violin(scale="width",trim = FALSE)+
  geom_boxplot(color="deepskyblue4", width=0.09) + 
  #geom_hline(yintercept = 50, color="red") + coord_cartesian(ylim=c(0,2500)) + 
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +
  #geom_text(aes(0,50,label = 50, vjust = 1.5, hjust=-0.25, color="red")) +
  labs(y="Read count",
       x="SNPs") + facet_grid(Pop~.)# -> q

ggplotly(q)

# barplot CAU
foo2.1 %>% filter( Fuzz == "11", Pop=="CAU") %>% 
  ggplot(aes(x=reorder(RefSNP,ReadCTotal),y=ReadCTotal)) +
  geom_bar(aes(fill=factor(rank)),stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_wrap(~Sample )

# barplot AFA
foo2.1 %>% filter( Fuzz == "11", Pop=="AFA") %>% 
  ggplot(aes(x=reorder(RefSNP,ReadCTotal),y=ReadCTotal)) +
  geom_bar(aes(fill=factor(rank)),stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_wrap(~Sample )

# barplot HIS
foo2.1 %>% filter( Fuzz == "11", Pop=="HIS") %>% 
  ggplot(aes(x=reorder(RefSNP,ReadCTotal),y=ReadCTotal)) +
  geom_bar(aes(fill=factor(rank)),stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_wrap(~Sample)

foo2.1 %>% filter( Fuzz == "11", Pop=="HIS") %>%
  filter(Sample=="29711" | Sample== "24507" | Sample == "27300" | Sample == "19971" | Sample == "13007") %>%
  ggplot(aes(x=reorder(RefSNP,ReadCTotal),y=ReadCTotal)) +
  geom_bar(aes(fill=factor(rank)),stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) + #ylim(0, 2000)
  facet_wrap(~Sample) #+ ylim(0, 2000)




# threshold set at 50 reads
foo2 %>% filter(Fuzz==11) %>%
  ggplot(aes(x=reorder(RefSNP,ReadCTotal),y=ReadCTotal)) + 
  geom_boxplot(color="deepskyblue4") + 
  geom_hline(yintercept = 50, color="red") + coord_cartesian(ylim=c(0,2500)) + 
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none") +
  geom_text(aes(0,50,label = 50, vjust = 1.5, hjust=-0.25, color="red")) +
  labs(y="Analytical threshold (Read count)",
       x="SNPs") + facet_grid(Pop~.)
#scale_y_continuous(breaks = sort(c(seq(0,2500, length.out=500), 50)))


foo2 %>% filter(Fuzz==11) %>% group_by(RefSNP) %>%
  mutate(y10=quantile(ReadCTotal,0.1),
         y15=quantile(ReadCTotal,0.15)
         ) -> foo4

foo1 %>%
  ggplot(aes(x=))
#-> foo3

foo4 %>% summarise(y0 = min(ReadCTotal), y25 = quantile(ReadCTotal, 0.25), 
                   y50 = median(ReadCTotal), y75 = quantile(ReadCTotal, 0.75),
                   y100 = max(ReadCTotal), y10 = unique(y10), y15 = unique(y15), 
                   # Pop=unique(Pop)
                   ) -> foo5

foo5 %>%
  ggplot(aes(x=RefSNP, ymin=y0, lower=y25, middle = y50, upper = y75, ymax = y100)) +
  geom_boxplot(stat = "identity") +
  geom_segment(aes(x=1:n-0.5,xend=1:n+0.5,y=y10,yend=y10),
               inherit.aes = FALSE, color="orange") +
  geom_segment(aes(x=1:n-0.5,xend=1:n+0.5,y=y15,yend=y15),
               color="red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


foo5 %>% mutate(row_id=row_number()) %>%
  ggplot(aes(x=RefSNP, ymin=y0, lower=y25, middle = y50, upper = y75, ymax = y100)) +
  geom_boxplot(stat = "identity") +
  geom_segment(aes(x=row_id-0.5,xend=row_id+0.5,y=y10,yend=y10),
                inherit.aes = FALSE, color="orange") +
  geom_segment(aes(x=row_id-0.5,xend=row_id+0.5,y=y15,yend=y15),
                color="red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#%>%
  foo2 %>% ggplot(aes(x=RefSNP,y=ReadCTotal)) + 
  geom_boxplot(aes(color=Pop),outlier.colour = "black", outlier.shape = 1) + 
  #geom_boxplot(outlier.colour = "black", outlier.shape = 1) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_segment(x=x0s, xend=x1s,y=y0s, yend=y0s)
  #facet_grid(Fuzz~., scales = "free")
  # geom_hline(yintercept = 4000) +
  # geom_hline(yintercept = 3000) +
  # geom_hline(yintercept = 2000) +
  # geom_hline(yintercept = 1000)  #-> q


































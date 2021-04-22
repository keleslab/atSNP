## ----eval=TRUE, echo=TRUE, results = "markup"-------------------------------------------------------------
library(atSNP)


## ----eval=TRUE, echo=TRUE, results = "markup"-------------------------------------------------------------
data(encode_library)
length(encode_motif)
encode_motif[1]


## ----eval=TRUE, echo=TRUE, results="markup",tidy=TRUE-----------------------------------------------------
encode_motif[[1]]
GetIUPACSequence(encode_motif[[1]])


## ----eval=TRUE, echo=TRUE, results = "markup",tidy=TRUE---------------------------------------------------
length(encode_motifinfo)
head(encode_motifinfo)


## ----eval=TRUE, echo=TRUE, results="markup",tidy=TRUE-----------------------------------------------------
encode_motifinfo[names(encode_motif[1])]


## ----eval=TRUE, echo = TRUE, results = "markup", tidy = TRUE----------------------------------------------
data(jaspar_library)
jaspar_motif[[1]]
jaspar_motifinfo[names(jaspar_motif[1])]


## ----eval=FALSE, echo=TRUE, results="hide"----------------------------------------------------------------
## ## Source: http://meme.nbcr.net/meme/doc/examples/sample-dna-motif.meme-io
## pwms <- LoadMotifLibrary(
##  urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/sample-dna-motif.meme-io.txt")
## 
## ## Source: http://compbio.mit.edu/encode-motifs/motifs.txt
## pwms <- LoadMotifLibrary(
##  urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/motifs.txt",
##  tag = ">", transpose = FALSE, field = 1,
##  sep = c("\t", " ", ">"), skipcols = 1,
##  skiprows = 1, pseudocount = 0)
## 
## ## Source: http://johnsonlab.ucsf.edu/mochi_files/JASPAR_motifs_H_sapiens.txt
## pwms <- LoadMotifLibrary(
##  urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/JASPAR_motifs_H_sapiens.txt",
##  tag = "/NAME",skiprows = 1, skipcols = 0, transpose = FALSE,
##  field = 2)
## 
## ## Source: http://jaspar.genereg.net/html/DOWNLOAD/ARCHIVE/JASPAR2010/all_data/matrix_only/matrix.txt
## pwms <- LoadMotifLibrary(
##  urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/matrix.txt",
##  tag = ">", skiprows = 1, skipcols = 1, transpose = TRUE,
##  field = 1, sep = c("\t", " ", "\\[", "\\]", ">"),
##  pseudocount = 1)
## 
## ## Source: http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt
## pwms <- LoadMotifLibrary(
##  urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/pfm_vertebrates.txt",
##  tag = ">", skiprows = 1, skipcols = 0, transpose = TRUE, field = 1,
##  sep = c(">", "\t", " "), pseudocount = 1)


## ----eval=FALSE, echo=TRUE, results="markup",tidy=FALSE---------------------------------------------------
## data(example)
## write.table(snp_tbl, file = "test_snp_file.txt",
##             row.names = FALSE, quote = FALSE)
## snp_info <- LoadSNPData("test_snp_file.txt", genome.lib = "BSgenome.Hsapiens.UCSC.hg38", half.window.size = 30, default.par = TRUE, mutation = FALSE)
## ncol(snp_info$sequence) == nrow(snp_tbl)
## snp_info$rsid.rm
## 


## ----eval=FALSE, echo=TRUE, results="markup",tidy=FALSE---------------------------------------------------
## mutation_info <- LoadSNPData("test_snp_file.txt", genome.lib = "BSgenome.Hsapiens.UCSC.hg38", half.window.size = 30, default.par = TRUE, mutation = TRUE)
## ncol(mutation_info$sequence) == nrow(snp_tbl)
## file.remove("test_snp_file.txt")


## ----eval=FALSE, echo=TRUE, results="markup", tidy=FALSE--------------------------------------------------
## snp_info1 <- LoadSNPData(snpids = c("rs5050", "rs616488", "rs11249433", "rs182799", "rs12565013", "rs11208590"), genome.lib ="BSgenome.Hsapiens.UCSC.hg38", snp.lib = "SNPlocs.Hsapiens.dbSNP144.GRCh38", half.window.size = 30, default.par = TRUE, mutation = FALSE)


## ----eval=TRUE, echo = TRUE, results="markup",tidy=FALSE--------------------------------------------------
snp_info2 <- LoadFastaData(ref.urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/sample_1.fasta", snp.urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/sample_2.fasta", default.par = TRUE)


## ----eval=TRUE, echo=TRUE, results="markup",tidy=TRUE-----------------------------------------------------
data(example)
names(motif_library)
str(snpInfo)
## to look at the motif information
data(encode_library)
encode_motifinfo[names(motif_library)]


## ----eval=TRUE, echo=TRUE, results="markup"---------------------------------------------------------------
atsnp.scores <- ComputeMotifScore(motif_library, snpInfo, ncores = 1)
atsnp.scores$snp.tbl
atsnp.scores$motif.scores


## ----eval=TRUE,echo=TRUE,results="markup"-----------------------------------------------------------------
atsnp.result <- ComputePValues(motif.lib = motif_library, snp.info = snpInfo,
                               motif.scores = atsnp.scores$motif.scores, ncores = 1, testing.mc=TRUE)
atsnp.result


## ----eval=TRUE, echo = TRUE, results="markup"-------------------------------------------------------------
head(atsnp.result[order(atsnp.result$pval_rank), c("snpid", "motif", "pval_ref", "pval_snp", "pval_rank")])


## ----eval=TRUE, echo = FALSE, results="hide"--------------------------------------------------------------
pval_rank_bh = p.adjust(atsnp.result$pval_rank, method = "BH")
atsnp.result = cbind(atsnp.result, pval_rank_bh)


## ----eval=TRUE, echo = FALSE, results="markup"------------------------------------------------------------
atsnp.result[c("snpid", "motif", "pval_rank", "pval_rank_bh")]


## ----eval=FALSE, echo =TRUE,results="markup"--------------------------------------------------------------
## library(qvalue)
## qval_rank = qvalue(atsnp.result$pval_rank, pi0=0.1)$qvalues
## atsnp.result = cbind(atsnp.result, qval_rank)


## ----eval=TRUE,echo=TRUE,results="markup"-----------------------------------------------------------------
match.subseq_result <- MatchSubsequence(snp.tbl = atsnp.scores$snp.tbl, motif.scores = atsnp.result, motif.lib = motif_library, snpids = c("rs53576", "rs7412"), motifs = names(motif_library)[1], ncores = 1)
match.subseq_result[c("snpid", "motif", "IUPAC", "ref_match_seq", "snp_match_seq")]


## ----eval=TRUE, echo=TRUE, fig.align="center", fig.height=6, fig.width=6, warning=TRUE, dpi=600, include=TRUE, results="markup"----
match.seq<-dtMotifMatch(atsnp.scores$snp.tbl, atsnp.scores$motif.scores, snpids="rs53576", motifs="SIX5_disc1", motif.lib = motif_library)
plotMotifMatch(match.seq,  motif.lib = motif_library)


## ----eval=TRUE,echo=FALSE,results="markup",cache=FALSE----------------------------------------------------
print(sessionInfo())


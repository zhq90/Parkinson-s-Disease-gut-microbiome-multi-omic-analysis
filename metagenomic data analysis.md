1. quality control

    We used trimmomatic for quality control

    ```bash
    java -jar trimmomatic-0.39.jar \
    PE -threads 2 -phred 33 reads1.fq.gz reads2.fq.gz \
    clean_reads1.fq.gz clean_reads2.fq.gz unpaired.fq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:150 ILLUMINACLIP:TruSeq3-PE.fa:0:0:0
    ```

2. reads contamination filter

    We used bowtie2 to align reads to hg38 and phiX genome to filter out human and virus contamination reads.

    ```bash
    bowtie2 -p 10 -x hg38_phiX_index -1 reads1.fq -2 reads2.fq\
    -S aligned.sam --un-conc-gz output_%.fastq.gz
    ```

3. reads classification

    We used kraken2 to classify the reads by using UHGG2 database ([ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gu...](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0/kraken2_db_uhgg_v2/))

    step1:

    ```bash
    kraken2 --db uhgg2_v2 --output sampleID.kraken2output \
    --report sampleID.kraken2report --report-zero-counts --paired \
    --gzip-compressed reads1.fq.gz reads2.fq.gz
    ```

    step2:

    ```bash
    bracken -d uhgg2_v2 -i sampleID.kraken2report -r 150 -t 0 -l S
    bracken -d uhgg2_v2 -i sampleID.kraken2report -r 150 -t 0 -l G
    ```

4. metacyc function level analysis

    step1: convert fq to fa and merge paired reads

    ```bash
    seqkit fq2fa reads1.fq.gz -o reads1.fa.gz
    seqkit fq2fa reads2.fq.gz -o reads2.fa.gz

    zcat reads1.fa.gz reads2.fa.gz > reads.fa.gz
    ```

    step2: 

    ```bash
    humann -i read.fa.gz -o ./ --input-format fasta.gz \
    --search-mode uniref90 --diamond /data/home/zhq/program/diamond-2.0.14/bin/
    ```

5. Taxonomic analysis  
    differential analysis

    ```R
    library(ALDEx2)
    library(tidyverse)

    df <- read.csv("S_level_counts_PDHP_matchedID_filtered.csv", header = TRUE, sep = ",", check.names = FALSE)
    dim(df)
    df <- dplyr::filter(df, name != "Unknown species")
    dim(df)
    metadata <- read.csv("metadata_matchedID.csv", header = TRUE, sep = ",", check.names = FALSE)
    pairinfo <- read.csv("paired_info.csv", header = TRUE, sep = ",", check.names = FALSE)

    df <- column_to_rownames(df, var = "name")
    colnames(metadata)
    metadata <- column_to_rownames(metadata, var = "sampleID")
    length(pairinfo$PD[1:30])

    PDHP_list <- append(pairinfo$PD[1:30], pairinfo$HP[1:30])
    dfPDHP <- df[,PDHP_list]
    mdPDHP <- metadata[PDHP_list,]

    PD_HP_da <- ALDEx2::aldex(reads = dfPDHP,
                              conditions = mdPDHP$Group,
                              test = "t", effect = TRUE, denom = "all",
                              mc.samples = 1000,
                              paired.test = TRUE,
                              useMC = TRUE)
    write.csv(file = "PDvsHP_result.csv", x = PD_HP_da)

    aldex.plot(PD_HP_da, type = "MW")

    save.image("DF_aldex2_result.RData")

    ```

    α-diversity analysis

    ```R
    library(vegan)
    library(ggplot2)
    library(tidyverse)
    library(ggordiplots)
    library(ggpubr)
    library(reshape2)

    df <- read.csv("S_level_counts_PDHP_matchedID_filtered.csv", header = TRUE, sep = ",", check.names = FALSE)
    metadata <- read.csv("metadata_matchedID.csv", header = TRUE, sep = ",", check.names = FALSE)
    pairinfo <- read.csv("paired_info.csv", header = TRUE, sep = ",", check.names = FALSE)

    df <- column_to_rownames(df, var = "name")
    colnames(metadata)
    metadata <- column_to_rownames(metadata, var = "sampleID")
    length(pairinfo$PD[1:30])

    PDHP_list <- append(pairinfo$PD[1:30], pairinfo$HP[1:30])
    dfPDHP <- df[,sort(PDHP_list, decreasing = TRUE)]
    dfPDHP <- t(dfPDHP)
    mdPDHP <- metadata[sort(PDHP_list, decreasing = TRUE),]
    mdPDHP$Group <- factor(mdPDHP$Group, levels = c("PD","HP"))
    rownames(dfPDHP) == rownames(mdPDHP)

    raremax <- min(rowSums(dfPDHP))
    dfPDHP_rarefy <- rrarefy(dfPDHP, sample = raremax)

    mdPDHP$alpha_shanno <- vegan::diversity(t(dfPDHP_rarefy),
                                            index = "shannon",
                                            MARGIN = 2)
    mdPDHP$alpha_invsimpson <- vegan::diversity(t(dfPDHP_rarefy),
                                                index = "invsimpson",
                                                MARGIN = 2)
    mdPDHP$alpha_simpson <- vegan::diversity(t(dfPDHP_rarefy),
                                             index = "simpson",
                                             MARGIN = 2)

    simpson_pvalue <- wilcox.test(alpha_simpson ~ Group, mdPDHP)
    shanno_pvalue <- wilcox.test(alpha_shanno ~ Group, mdPDHP)
    invsimpson_pvalue <- wilcox.test(alpha_invsimpson ~ Group, mdPDHP)

    simpson_pvalue$p.value

    alpha_plot <- function(alpha_type, y_axis, outputname){
      temp <- mdPDHP[c(alpha_type,"Group")]
      temp <- reshape2::melt(temp)
      y_lab <- substring(alpha_type, 7)
      p <- ggboxplot(data = temp, x = "Group", y = "value", #orientation = "horizontal",
                     color = "Group", linetype = 1, size = 0.8, point.size = 0.8, add = "jitter", add.params = list(alpha = 0.6))+
        theme_bw()+
        scale_color_manual(values = c("#FFB889","#9BADB2"))+
        stat_compare_means(method = "wilcox.test", label.y = y_axis, label.x = 1.4, label = "p.signif")+
        scale_alpha_discrete(0.6)+
        labs(x = "", y = y_lab)+
        theme(axis.title.y = element_text(size = 14),
              axis.text.x = element_text(size = 14))
      
      ggsave(filename = outputname, plot = p, limitsize = FALSE,
             width = 2.5, height = 5, dpi = 300)
    }
    mdPDHP$alpha_shanno

    alpha_plot(alpha_type = "alpha_invsimpson", y_axis = 60, "invsimpson_boxplot.pdf")
    alpha_plot(alpha_type = "alpha_shanno", y_axis = 5.4, "shanno_boxplot.pdf")
    alpha_plot(alpha_type = "alpha_simpson", y_axis = 1.0, "simpson_boxplot.pdf")


    ```

    β-diversity

    ```R
    library(vegan)
    library(tidyverse)
    library(ggplot2)
    library(ggpubr)

    df <- read.csv("S_level_counts_PDHP_matchedID_filtered.csv", header = TRUE, sep = ",", check.names = FALSE)
    metadata <- read.csv("metadata_matchedID.csv", header = TRUE, sep = ",", check.names = FALSE)
    pairinfo <- read.csv("paired_info.csv", header = TRUE, sep = ",", check.names = FALSE)

    df <- column_to_rownames(df, var = "name")
    colnames(metadata)
    metadata <- column_to_rownames(metadata, var = "sampleID")
    length(pairinfo$PD[1:30])

    PDHP_list <- append(pairinfo$PD[1:30], pairinfo$HP[1:30])
    dfPDHP <- df[,sort(PDHP_list, decreasing = TRUE)]
    dfPDHP <- t(dfPDHP)
    mdPDHP <- metadata[sort(PDHP_list, decreasing = TRUE),]

    rownames(dfPDHP) == rownames(mdPDHP)

    raremax <- min(rowSums(dfPDHP))
    dfPDHP_rarefy <- rrarefy(dfPDHP, sample = raremax)

    data.bray <- vegdist(dfPDHP_rarefy, method = "bray")
    mean(data.bray)

    beta <- adonis2(data.bray ~ Group, data = mdPDHP, permutations = 1000)

    df <- data.frame(distance = append(mod$distances[pairinfo$锘縋D[1:30]], mod$distances[pairinfo$HP[1:30]]), Group = mdPDHP$Group)

    p <- ggboxplot(data = df, x = "Group", y = "distance",
                   linetype = 1, size = 1, point.size = 1.5, add = "jitter",
                   add.params = list(alpha = 0.6), color = "Group", orientation = "horizontal")+
      theme_bw()+
      scale_color_manual(values = c("#FFB889","#9BADB2"))+
      stat_compare_means(method = "wilcox.test", comparisons = list(c("PD","HP")), paired = TRUE)

    ggsave(filename = "horizontal_beta_boxplot.pdf", plot = p, limitsize = FALSE,
            dpi = 300, width = 10, height = 4)

    ```

6. differential abundance species association construction

    We used R package NetCoMi to construct species association network.

    ```R
    library(NetCoMi)
    library(tidyverse)
    library(phyloseq)

    df <- read.csv("sig_species_counts.csv", header = T, sep = ",", check.names = T)
    md <- read.csv("metadata_PDHP.csv", header = T, sep = ",")
    df <- column_to_rownames(df, var = "name")
    md <- column_to_rownames(md, var = "sampleID")
    rownames(df)
    rownames(md) <- make.names(rownames(md))

    taxmat <- matrix("-", nrow = nrow(df), ncol = 7)
    rownames(taxmat) <- rownames(df)
    colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    taxmat[,7] <- rownames(taxmat)
    taxmat[,7]

    otu <- otu_table(df, taxa_are_rows = T)
    tax <- tax_table(taxmat)
    pseq <- phyloseq(otu, tax, sample_data(md))                           

    net_spacc <- netConstruct(t(df),
                              filtTax = "none",
                              filtTaxPar = "none",
                              filtSamp = "none",
                              measure = "sparcc",
                              measurePar = NULL,
                              normMethod = "none",
                              zeroMethod = "none",
                              sparsMethod = "t-test",
                              adjust = "adaptBH",
                              alpha = 0.05,
                              dissFunc = "signed",
                              lfdrThresh = 0.05,
                              verbose = 2,
                              seed = 99)

    prop_spacc <- netAnalyze(net_spacc,
                             centrLCC = TRUE,
                             clustMethod = "cluster_optimal",
                             hubPar = "eigenvector",
                             weightDeg = F,
                             normDeg = F)

    p <- plot(prop_spacc,
              nodeColor = "cluster",
              labelFont = 1,labelScale = F
              )
    legend(0.4,1.1, title = "estimated association:",
           legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
           bty = "n", horiz = TRUE)

    plot(prop_spacc,
         labels = FALSE,
         nodeColor = "cluster")
    legend(0.4,1.1,
           legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
           bty = "n", horiz = TRUE)

    write.csv(file = "cluster_info.csv", x = prop_spacc$clustering$clust1)

    save.image("file.RData")


    write.csv(file = "association_ma")
    ```

7. Gut micrbiome gene feature selection

    step1: we aligned cleaned reads to GMGC ([https://gmgc.embl.de/](https://gmgc.embl.de/)) gut microboome gene catagolue.

    ```bash
    salmon quant -l A -i GMGC10_humangut95nr_norare_salmon_index \
    -1 reads1.fq.gz -2 reads2.fq.gz -o sampleID.salmonoutput --meta \
    -p 20

    ```
    step2: We removed gene which was not expressed in more than half of the cohort, and used Wilcoxon signed-rank test to performed differential analysis.

    ```python
    import pandas as pd
    from scipy import stats
    from statsmodels.stats.multitest import fdrcorrection

    df = pd.read_csv("tpm_PDHP_filtered.csv", header = 0, sep = ',')
    pairinfo = pd.read_csv("PDHP_paired_info.csv", header = 0, sep = ",")

    df = df.set_index("ID")
    df = df.fillna(0)

    PD_samples = list(pairinfo["PD"])
    HP_samples = list(pairinfo["HP"])

    def wilcox(dm, x, y, i):
        w, p = stats.wilcoxon(dm.loc[i,x], dm.loc[i,y])
        #print(p)
        return p

    df["pvalue"] = -1
    for i in df.index:
        df.loc[i,"pvalue"] = wilcox(df, PD_samples, HP_samples, i)

    df["FDR"] = fdrcorrection(df["pvalue"], method = "indep")[1]

    df[["pvalue","FDR"]].to_csv("wilcoxon_results.csv", index = True, index_label = "ID")

    dm = df[df["FDR"] < 0.05]
    dm.shape
    dm[["pvalue","FDR"]].to_csv("sig_wilcoxon_results.csv", index = True, index_label = "ID")

    dm.index

    matrix = pd.read_csv("tpm_PDHP_filtered.csv", header = 0, sep = ",")
    matrix = matrix.set_index("ID")
    matrix = matrix.loc[dm.index, :]
    matrix.to_csv("sig_wilcox_matrix.csv", index = True, index_label = "ID")

    ```
    step3: We used lassonet to perform feature selection and select top 25 gene as features to construc disease prediction model in the next.

    ```python
    import pandas as pd
    from lassonet import LassoNetClassifierCV
    from sklearn.preprocessing import LabelEncoder
    import numpy as np

    df = pd.read_csv("sig_wilcox_matrix.csv", header = 0, sep = ",")
    md = pd.read_csv("metadataPDHP.csv", header = 0, sep = ",")

    df = df.set_index("ID")
    df = df.T
    md = md.set_index("sampleID")
    md = md[["Group"]]
    df = df.loc[md.index,:]
    md = LabelEncoder().fit_transform(md)

    model = LassoNetClassifierCV(cv = 5,
                                 lambda_start = 0.001,
                                 random_state=1)
    path = model.fit(np.array(df), np.array(md))

    len(model.feature_importances_)
    df.columns
    len(model.best_selected_)

    results = pd.DataFrame(np.full((len(df.columns),2),""), index = df.columns, columns=["feature_importance","best_selected"])
    results["feature_importance"] = model.feature_importances_
    results["feature_selected"] = model.best_selected_
    results.to_csv("lassonet_feature_results.csv", index = True, index_label = "ID")

    ```

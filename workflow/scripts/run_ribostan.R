input_gtf <- "resources/gencode.v37.primary_assembly.annotation.CCDS.gtf"
input_fasta <- "resources/GRCh38.primary_assembly.genome.fa"
input_bam <- "results/split_bam/transcriptome/human/ribo_4E_minus_rep1_01_R1.bam"

{input.gtf}
{input.fasta}
{output.anno}
{input.bam}
{output.rpfs}
{output.offsets_df}
{output.psites}
{output.filtered_anno}
{output.uorf_anno} -> uorf_grl
{output.quant}
{output.tsv}

{output.morf_grl}

library(Ribostan)

# Generate unfiltered ORF annotation (including potential uORFs)
anno <- load_annotation(gtf = input_gtf, fafile = input_fasta, add_uorfs = TRUE)
# saveRDS(anno, file = "{output.anno}")

# Ribosome footprints (ribosome protected fragments)
rpfs <- get_readgr(ribobam = input_bam, anno)
# saveRDS(rpfs, file = "{output.rpfs}")

# Offsets
offsets_df <- get_offsets(rpfs, anno)
# saveRDS(offsets_df, file = "{output.offsets_df}")

# P-sites
psites <- get_psite_gr(rpfs, offsets_df, anno)
# saveRDS(psites, file = "{output.psites}")

# ORF annotation filtered for periodicity (multitaper test)
filtered_anno <- periodicity_filter_uORFs(psites, anno, remove=TRUE)
# saveRDS(filtered_anno, file = "{output.filtered_anno}")

# Annotation with only mORFs
morf_grl <- filtered_anno$cdsgrl[!filtered_anno$uORF]
# saveRDS(morf_grl, file = "{output.morf_grl}")

# Annotation with only uORFs
uorf_grl <- filtered_anno$cdsgrl[filtered_anno$uORF]
# saveRDS(uorf_grl, file = "{output.uorf_anno}")

# Calculate ribosomal TPMs from filtered annotation
ritpms <- get_ritpms(psites, filtered_anno)
# saveRDS()

# Total read count from all RPFs
n_reads <- dplyr::n_distinct(names(rpfs))

# ORF lengths
orf_lengths <- GenomicRanges::width(filtered_anno$trspacecds)
names(orf_lengths) <- names(filtered_anno$trspacecds)

# 
nucfracs <- (ritpms[names(orf_lengths)] * orf_lengths[names(orf_lengths)])
nucfracs <- nucfracs / sum(nucfracs)



# Split ribosomal TPMs by mORF and uORF
morf_ritpms <- ritpms[names(morf_grl)]
uorf_ritpms <- ritpms[names(uorf_grl)]




# uorf_lengths <- GenomicRanges::width(filtered_anno$trspacecds[names(uorf_ritpms)])



nucfracs <- (uorf_ritpms * lengths)
nucfracs <- nucfracs / sum(nucfracs)
counts <- n_reads * nucfracs
counts <- tibble::enframe(counts, "Name", "NumReads")
ritpmdf <- tibble::enframe(ritpms, "Name", "ritpm")

uorflens <- filtered_anno$trspacecds[names(uorf_ritpms)] %>%
  GenomicRanges::width() %>%
  setNames(names(filtered_anno$trspacecds[names(uorf_ritpms)])) %>%
  tibble::enframe("Name", "Length") %>%
  dplyr::mutate(EffectiveLength = .data$Length)

output <- uorflens %>%
  dplyr::left_join(ritpmdf) %>%
  dplyr::left_join(counts)

saveRDS(uorf_ritpms, file = "{output.quant}")

output %>% readr::write_tsv("{output.tsv}")







"results/ribostan/{sample}/{sample}.ribostan.human_morf_quant.rds"
"results/ribostan/{sample}/{sample}.ribostan.human_morf_quant.tsv"

"results/ribostan/{sample}/{sample}.ribostan.human_offsets.rds"

"results/ribostan/{sample}/{sample}.ribostan.yeast_morf_quant.rds"
"results/ribostan/{sample}/{sample}.ribostan.yeast_morf_quant.tsv"

        anno = "results/ribostan/{sample}/{sample}.ribostan.human_anno.rds",
        rpfs = "results/ribostan/{sample}/{sample}.ribostan.human_rpfs.rds",
        offsets_df = "results/ribostan/{sample}/{sample}.ribostan.human_offsets_df.rds",
        psites = "results/ribostan/{sample}/{sample}.ribostan.human_psites.rds",
        filtered_anno = "results/ribostan/{sample}/{sample}.ribostan.human_filtered_anno.rds",
        uorf_anno = "results/ribostan/{sample}/{sample}.ribostan.human_uorf_anno.rds",
        quant = "results/ribostan/{sample}/{sample}.ribostan.human_uorf_quant.rds",
        tsv = "results/ribostan/{sample}/{sample}.ribostan.human_uorf_quant.tsv"






rule ribostan_human_morfs:
    input:
        bam = "results/split_bam/transcriptome/human/{sample}.bam",
        fasta = config["human_genome"],
        gtf = config["human_gtf"]
    output:
        quant = "results/ribostan/{sample}/{sample}.ribostan.human_morf_quant.rds",
        tsv = "results/ribostan/{sample}/{sample}.ribostan.human_morf_quant.tsv"
    threads: 2

        library(Ribostan);
        anno <- load_annotation(gtf = "{input.gtf}", fafile = "{input.fasta}", add_uorfs = FALSE);
        rpfs <- get_readgr(ribobam = "{input.bam}", anno);
        offsets_df <- get_offsets(rpfs, anno);
        psites <- get_psite_gr(rpfs, offsets_df, anno);
        ritpms <- get_ritpms(psites, anno);
        n_reads <- dplyr::n_distinct(names(rpfs));
        lengths <- GenomicRanges::width(anno$trspacecds[names(ritpms)]);
        nucfracs <- (ritpms * lengths);
        nucfracs <- nucfracs / sum(nucfracs);
        counts <- n_reads * nucfracs;
        counts <- tibble::enframe(counts, "Name", "NumReads");
        ritpmdf <- tibble::enframe(ritpms, "Name", "ritpm");
        cdslens <- anno$trspacecds %>%
            GenomicRanges::width() %>%
            setNames(names(anno$trspacecds)) %>%
            tibble::enframe("Name", "Length") %>%
            dplyr::mutate(EffectiveLength = .data$Length);
        output <- cdslens %>%
            dplyr::left_join(ritpmdf) %>%
            dplyr::left_join(counts);
        saveRDS(ritpms, file = "{output.quant}");
        output %>% readr::write_tsv("{output.tsv}")



rule get_offsets:
    input:
        bam = "results/split_bam/transcriptome/human/{sample}.bam",
        fasta = config["human_genome"],
        gtf = config["human_gtf"]
    output:
        "results/ribostan/{sample}/{sample}.ribostan.human_offsets.rds"
    threads: 2
    shell: r"""
        R -e 'library(Ribostan);
        anno <- load_annotation(gtf = "{input.gtf}", fafile = "{input.fasta}", add_uorfs = FALSE);
        rpfs <- get_readgr(ribobam = "{input.bam}", anno);
        offsets_df <- get_offsets(rpfs, anno);
        saveRDS(offsets_df, file = "{output}")'
        """


rule ribostan_yeast_morfs:
    input:
        bam = "results/split_bam/transcriptome/yeast/{sample}.bam",
        fasta = config["yeast_genome"],
        gtf = config["yeast_gtf"],
        offsets = "results/ribostan/{sample}/{sample}.ribostan.human_offsets.rds"
    output:
        quant = "results/ribostan/{sample}/{sample}.ribostan.yeast_morf_quant.rds",
        tsv = "results/ribostan/{sample}/{sample}.ribostan.yeast_morf_quant.tsv"
    threads: 2

        library(Ribostan);
        anno <- load_annotation(gtf = "{input.gtf}", fafile = "{input.fasta}", add_uorfs = FALSE);
        rpfs <- get_readgr(ribobam = "{input.bam}", anno);
        offsets_df <- readRDS(file = "{input.offsets}");
        psites <- get_psite_gr(rpfs, offsets_df, anno);
        ritpms <- get_ritpms(psites, anno);
        n_reads <- dplyr::n_distinct(names(rpfs));
        lengths <- GenomicRanges::width(anno$trspacecds[names(ritpms)]);
        nucfracs <- (ritpms * lengths);
        nucfracs <- nucfracs / sum(nucfracs);
        counts <- n_reads * nucfracs;
        counts <- tibble::enframe(counts, "Name", "NumReads");
        ritpmdf <- tibble::enframe(ritpms, "Name", "ritpm");
        cdslens <- anno$trspacecds %>%
            GenomicRanges::width() %>%
            setNames(names(anno$trspacecds)) %>%
            tibble::enframe("Name", "Length") %>%
            dplyr::mutate(EffectiveLength = .data$Length);
        output <- cdslens %>%
            dplyr::left_join(ritpmdf) %>%
            dplyr::left_join(counts);
        saveRDS(ritpms, file = "{output.quant}");
        output %>% readr::write_tsv("{output.tsv}")



rule ribostan_human_uorfs:
    input:
        bam = "results/split_bam/transcriptome/human/{sample}.bam",
        fasta = config["human_genome"],
        gtf = config["human_gtf"]
    output:
        anno = "results/ribostan/{sample}/{sample}.ribostan.human_anno.rds",
        rpfs = "results/ribostan/{sample}/{sample}.ribostan.human_rpfs.rds",
        offsets_df = "results/ribostan/{sample}/{sample}.ribostan.human_offsets_df.rds",
        psites = "results/ribostan/{sample}/{sample}.ribostan.human_psites.rds",
        filtered_anno = "results/ribostan/{sample}/{sample}.ribostan.human_filtered_anno.rds",
        uorf_anno = "results/ribostan/{sample}/{sample}.ribostan.human_uorf_anno.rds",
        quant = "results/ribostan/{sample}/{sample}.ribostan.human_uorf_quant.rds",
        tsv = "results/ribostan/{sample}/{sample}.ribostan.human_uorf_quant.tsv"
    threads: 2

        library(Ribostan);
        anno <- load_annotation(gtf = "{input.gtf}", fafile = "{input.fasta}", add_uorfs = TRUE);
        saveRDS(anno, file = "{output.anno}");
        rpfs <- get_readgr(ribobam = "{input.bam}", anno);
        saveRDS(rpfs, file = "{output.rpfs}");
        offsets_df <- get_offsets(rpfs, anno);
        saveRDS(offsets_df, file = "{output.offsets_df}");
        psites <- get_psite_gr(rpfs, offsets_df, anno);
        saveRDS(psites, file = "{output.psites}");
        filtered_anno <- periodicity_filter_uORFs(psites, anno, remove=TRUE);
        saveRDS(filtered_anno, file = "{output.filtered_anno}");
        uorfgrl <- filtered_anno$cdsgrl[filtered_anno$uORF];
        saveRDS(uorfgrl, file = "{output.uorf_anno}");
        ritpms <- get_ritpms(psites, filtered_anno);
        uorf_ritpms <- ritpms[names(uorfgrl)];
        n_reads <- dplyr::n_distinct(names(rpfs));
        lengths <- GenomicRanges::width(filtered_anno$trspacecds[names(uorf_ritpms)]);
        nucfracs <- (uorf_ritpms * lengths);
        nucfracs <- nucfracs / sum(nucfracs);
        counts <- n_reads * nucfracs;
        counts <- tibble::enframe(counts, "Name", "NumReads");
        ritpmdf <- tibble::enframe(ritpms, "Name", "ritpm");
        uorflens <- filtered_anno$trspacecds[names(uorf_ritpms)] %>%
            GenomicRanges::width() %>%
            setNames(names(filtered_anno$trspacecds[names(uorf_ritpms)])) %>%
            tibble::enframe("Name", "Length") %>%
            dplyr::mutate(EffectiveLength = .data$Length);
        output <- uorflens %>%
            dplyr::left_join(ritpmdf) %>%
            dplyr::left_join(counts);
        saveRDS(uorf_ritpms, file = "{output.quant}");
        output %>% readr::write_tsv("{output.tsv}")

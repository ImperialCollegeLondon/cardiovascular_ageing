make_coloc_tests <- function(gene_traits,
                             gene_ebi_eqtl_locus,
                             myxlim,
                             sdY,
                             sdY_gtex=8.744022,
                             debug=FALSE) {
  if(debug) {
    browser()
  }
  # prepare for coloc
  C1 <- gene_traits %>%
    dplyr::filter(chr == my_snps$chr,
                  pos >= myxlim[1] & pos <= myxlim[2]) %>%
    dedupgw(snpcol = "rs_id", pcol = "p") %>% # implicit rename of snpcol and pcol to snp and p
    dplyr::select(snp, BETA, SE, p) %>%
    dplyr::rename(pvalues = p, beta = BETA, varbeta = SE) %>%
    dplyr::mutate(varbeta = varbeta * varbeta)
  # test coloc for all eqtl
  untrait <- unique(gene_ebi_eqtl_locus$trait)
  # k <- 1
  coloc_all <- list()
  coloc_2 <- list()
  for(k in 1:length(untrait)) {
    unk <- untrait[k]
    C2 <- gene_ebi_eqtl_locus %>%
      dplyr::filter(chr == my_snps$chr,
                    trait == unk,
                    pos >= myxlim[1] & pos <= myxlim[2]) %>%
      dedupgw(snpcol = "rs_id", pcol = "p") %>%
      dplyr::ungroup() %>%
      dplyr::select(snp, beta, se, p, maf, an) %>%
      dplyr::rename(pvalues = p, varbeta = se) %>%
      dplyr::mutate(varbeta = varbeta * varbeta) %>%
      dplyr::filter(!is.na(beta)) %>%
      unique()
    # get colocsnps, i.e. snps in both tables
    colocsnps <- intersect(C1$snp, C2$snp)
    #message(sprintf("Running coloc analysis - N snps C1 / C2 / intersect: %i / %i / %i",
    #                length(unique(C1$snp)),
    #                length(unique(C2$snp)),
    #                length(unique(colocsnps))))
    coloc_input <- C1 %>%
       dplyr::inner_join(C2, by="snp") %>%
       dplyr::rename(beta1 = beta.x, se1 = varbeta.x,
                     beta2= beta.y, se2 = varbeta.y) %>%
       dplyr::filter(!is.na(beta1) & !is.na(beta2) & !is.na(se1) & !is.na(se2)) %>%
       dplyr::mutate(type = "quant")
    coloc_2[[unk]] <- UKBRlib::coloc_wrapper(coloc_input)

    D1 <- make_coloc_dataset(C1, sdY = sdY)
    #D1x <- D1[c("pvalues", "beta", "sdY","type","snp")]
    #coloc::check_dataset(D1) #, warn.minp=1e-10) #, type = "quant")
    D2 <- make_coloc_dataset(C2, sdY = sdY_gtex, debug=F) # use maf and an
    #coloc::check_dataset(D2) #, warn.minp=1e-10) #, type = "quant")
    suppressWarnings(coloc_all[[unk]] <- coloc::coloc.abf(dataset1 = D1, dataset2 = D2))
  }
  coloctests <- dplyr::bind_rows(lapply(coloc_all, function(x) {
    x$summary
  })) %>%
    dplyr::mutate(trait = untrait) %>%
    dplyr::relocate(trait)
  #return(coloctests)
  return(list(coloctests=coloctests, coloctests_ukbrlib = coloc_2))
}


dedupgw <- function(gwtab, snpcol = "rs_id", pcol = "p",
                    graceful_duplicates = TRUE,
                    debug=FALSE) {
  if(debug) {
    browser()
  }
  gwtab <- gwtab %>%
    dplyr::rename(snp = `snpcol`,
                  p = `pcol`)
  tmp1 <- gwtab %>%
    dplyr::group_by(snp) %>%
    dplyr::summarise(p = min(p),
                     .groups = "drop")

  out <- tmp1 %>%
    dplyr::left_join(gwtab, by=c("snp", "p"))
  print(sprintf("nrow before dedup: %i", nrow(gwtab)))
  print(sprintf("nrow after  dedup: %i", nrow(out)))
  dupsleft <- which(table(out$snp)>1)
  if(length(dupsleft)>0) {
    msg <- sprintf("Ohjemine, still duplicate snp ids left. Something is wrong. Ids: %s",
            paste(names(dupsleft), collapse=","))
    if(graceful_duplicates) {
      warning(msg)
      out <- out %>%
        dplyr::filter(!snp %in% names(dupsleft))
    } else {
      stop(msg)
    }
  }
  return(out)
}
make_coloc_dataset <- function(C1x, sdY=NULL, debug=FALSE) {
  if(debug) {
    browser()
  }
  stopifnot(all(c("pvalues","beta","varbeta","snp") %in% colnames(C1x)))
  D1x <- C1x %>%
    #dplyr::rename(pvalues = p, varbeta = se) %>%
    #dplyr::rename(pvalues = p, beta = beta1, varbeta = se1, snp = snp) %>%
    dplyr::filter(!is.na(pvalues) & !is.na(beta) & !is.na(varbeta) & !is.na(snp)) %>%
    #dplyr::filter(!is.na(maf) & !is.na(an)) %>%
    dplyr::filter(beta != Inf)
  # deduplicate
  tmp1 <- D1x %>%
    dplyr::group_by(snp) %>%
    dplyr::summarise(pvalues = min(pvalues))
  D1x <- tmp1 %>%
    dplyr::left_join(D1x, by=c("snp", "pvalues"))

  if(is.null(sdY)) {
    # requires maf and n to estimate sdY
    # sdY used to interprete scale of betas
    stopifnot(all(c("maf", "an") %in% colnames(C1x)))
    out <- list(snp = D1x$snp,
                pvalues = D1x$pvalues,
                beta = D1x$beta,
                varbeta = D1x$varbeta,
                maf = D1x$maf,
                an = D1x$an,
                type = "quant")
  } else {
    out <- list(snp = D1x$snp,
       pvalues = D1x$pvalues,
       beta = D1x$beta,
       varbeta = D1x$varbeta,
       sdY = sdY,
       type = "quant")
  }
  return(out)
}

# functions for finemapping and coloc
myliftover_hg_37_to_38 <- function(df,
                                   chrom_col="Chromosome",
                                   pos_col="Pos_start",
                                   debug=F) {
  if(debug) {
    browser()
  }
  stopifnot(all(c(chrom_col, pos_col) %in% colnames(df)))
  #df<-df %>% dplyr::mutate(Chromosome=!!sym(chrom_col),Pos=!!sym(pos_col))
  # make sure the input Chromosome for makeGRangesFromDataFrame is in format chrX, chrY
  df <- df %>%
    dplyr::rename(chr = !!sym(chrom_col),
                  start = !!sym(pos_col)) %>%
    dplyr::mutate(end = as.numeric(start),
                  start = as.numeric(start),
                  chr = as.character(chr)) %>%
    dplyr::mutate(chr_temp=if_else(chr=='23','X',chr)) %>%
    dplyr::mutate(chr_temp=if_else(chr_temp=='24','Y',chr_temp))
  # load chaining info
  #chain <- rtracklayer::import.chain(config::get("19_38_chain_file"))
  chain <- rtracklayer::import.chain("/Users/gfkqx/gitlab/bullseye/r_data/hg19ToHg38.over.chain")

  # maybe construct from this string the GR object?
  #df <- df %>%
  #  dplyr::mutate(lostr = sprintf("%s:%s-%s", Chromosome, Pos, Pos))
  gr <- GenomicRanges::makeGRangesFromDataFrame(df %>%
                                                  dplyr::mutate(seqname=paste0('chr',as.character(chr_temp))),
                                                seqnames.field='seqname',
                                                start.field='start',
                                                end.field='end',
                                                strand='*', keep.extra.columns=TRUE)
  gr_up <- rtracklayer::liftOver(gr,chain) # gives info back which sequences are not found!!
  dfout <- as.data.frame(unlist(gr_up)) %>%
    dplyr::mutate(pos=as.numeric(start)) %>%
    dplyr::relocate(chr,pos) %>%
    dplyr::select(-seqnames,-start,-end,-width,-strand,-chr_temp) %>% as_tibble() %>%
    dplyr::mutate(chr=as.character(chr))
  # reset column names?
  #colnames(dfout)[colnames(dfout)=='Chromosome']<-chrom_col
  #colnames(dfout)[colnames(dfout)=='Pos']<-pos_col
  return(dfout)
}

liftover_ss <- function(filename,
                        mytrait,
                        debug=FALSE) {
  if(debug) {
    browser()
  }
  stopifnot(file.exists(filename))
  # read the hg37 data
  #new_names = c(rs_id="ID", chr="seqnames", pos="start", ref="REF", alt="ALT", p="P") # for consistency with open targets
  gene_traits <- readr::read_delim(filename, show_col_types = FALSE) %>%
    dplyr::rename(rs_id=ID, chr=Chromosome,
                  pos = Position, ref = REF, alt = ALT, p = P) %>%
    dplyr::relocate(rs_id, chr, pos, ref, alt, p, BETA, SE, T_STAT) %>%
    dplyr::arrange(p) %>%
    #dplyr::filter(p <= 1e-5) %>%
    dplyr::mutate(trait = mytrait) %>%
    dplyr::relocate(trait)

  Sys.setenv(R_CONFIG_FILE="~/gitlab/ukb-data-pipelines/config.yml")
  Sys.setenv(R_CONFIG_ACTIVE = "standard")
  ddd <- myliftover_hg_37_to_38(gene_traits, chrom_col = "chr", pos="pos", debug=F)
  return(ddd)
}

compute_ldmats <- function(gwas, pcut = 5e-8,
                           debug=FALSE) {
  if(debug) {
    browser()
  }

  # select significant
  gwassig <- gwas %>%
    dplyr::filter(p<=pcut) %>%
    dplyr::filter(stringr::str_detect(string = rs_id, pattern = "^rs.*$"))

  ldmats <- list()
  allclust <- list()
  uchr <- sort(unique(gwassig$chr))
  for(chrx in uchr) {
    xn <- paste0("chr",chrx)
    print(sprintf("Compute ldmat for %s", xn))
    # chrx <- uchr[1]
    gwchr <- gwassig %>%
      dplyr::filter(chr == chrx) %>%
      dplyr::arrange(p)
    mysnps <- gwchr$rs_id
    ttx <- try(ldx <- TwoSampleMR::ld_matrix(mysnps, with_alleles = TRUE, pop = "EUR"))
    if("try-error" %in% class(ttx)) {
      res <- ieugwasr::api_query("ld/matrix", query = list(rsid = mysnps, pop = "EUR"),
                                 access_token = NULL) %>% ieugwasr::get_query_content()
      length(res$snplist)
      print(sprintf("Error in fetching LDinfo, using lowest pvalue snp."))
      ldx <- matrix(1, nrow=length(mysnps), ncol=length(mysnps), dimnames=list(mysnps, mysnps))
      # ldx[diag(ldx)] <- 1 # set to 'no LD', to put them in the same ld bucket
    }
    ldmats[[xn]] <- ldx
  }
  return(ldmats)
}


identify_lead_snps_byproximity <- function(gwas,
                                           ldmats, # ldmatrices
                                           pcut = 5e-8,
                                           ldcut = 0.9,
                                           plot_clusters = TRUE,
                                           debug=FALSE) {
  if(debug) {
    browser()
  }

  # select significant
  gwassig <- gwas %>%
    dplyr::filter(p<=pcut) %>%
    dplyr::filter(stringr::str_detect(string = rs_id, pattern = "^rs.*$"))

  allclust <- list()
  uchr <- sort(unique(gwassig$chr))
  for(chrx in uchr) {
    xn <- paste0("chr",chrx)
    message(sprintf("Finding leads for %s", xn))
    # chrx <- uchr[1]
    gwchr <- gwassig %>%
      dplyr::filter(chr == chrx) %>%
      dplyr::arrange(p)
    mysnps <- gwchr$rs_id
    if(length(mysnps) == 1) {
      ldclust <- list()
      ldclust[[mysnps]] <- mysnps
    } else {
      if(is.null(ldmats)) {
        ldx <- TwoSampleMR::ld_matrix(mysnps, with_alleles = TRUE, pop = "EUR")
      } else {
        ldx <- ldmats[[xn]]
      }
      # make clusters
      if(all(ldx==1)) {
        # no ld info was present, take snps as they come
        ldsnps <- colnames(ldx)
      } else {
        pc <- pheatmap::pheatmap(ldx, show_colnames = FALSE, show_rownames = FALSE)
        # get ordering of snps according to hclust object
        #browser()
        hc <- pc$tree_row
        ldsnps <- labels(hc)
      }
      # identify clusters via ldmat
      ldclust <- list()
      # init while loop
      ldsnp <- ldsnps[1] # current cluster id
      ldsnps_remain <- setdiff(ldsnps, ldsnp)
      ldclust[[ldsnp]] <- ldsnp
      #for(ldsnp in ldsnps) {
      cnt <- 1
      while(length(ldsnps_remain)>0) {
        # ldsnp <- ldsnps_remain[1]
        # go trough remaining ldsnps, check LD
        # make a new cluster, named by initial lead snp
        ldsx <- ldsnps_remain[1]
        ldsnps_remain <- setdiff(ldsnps_remain, ldsx)
        # check ld, if <= 0.8 => next cluster, else in cluster
        myld <- ldx[ldsnp, ldsx]
        if(myld >= 0.8) {
          ldclust[[ldsnp]] <- c(ldclust[[ldsnp]], ldsx)
        } else {
          # new cluster
          ldclust[[ldsx]] <- ldsx
          ldsnp <- ldsx
        }
        cnt <- cnt + 1
        if(cnt == 10000) {
          break
        }
      }
    }

    # set representative as snp with lowest p!
    # i <- 1
    ldclustrep <- list()
    for(i in 1:length(ldclust)) {
      xnx <- names(ldclust)[i]
      ldobj <- ldclust[[xnx]]
      ldobjsnp <- gsub("_.*$", "", ldobj) # formating from ldmat computation needs to be reverted
      gwsel <- gwassig %>%
        dplyr::filter(rs_id %in% ldobjsnp) %>%
        dplyr::arrange(p) %>%
        dplyr::slice(1) # take first as rep
      ldclustrep[[xnx]] <- gwsel
    }
    # save representative ld cluster snps -> can be more than one!
    allclust[[xn]] <- ldclustrep
  }
  # make dataframes for each chromosome, holding all potential independent loci
  allclust <- lapply(allclust, dplyr::bind_rows)
  return(allclust)
}





identify_lead_snps <- function(gwas,
                               pcut = 5e-8,
                               k = 9,
                               ldx = NULL,
                               plot_clusters = TRUE,
                               debug=FALSE) {
  if(debug) {
    browser()
  }

  # select significant
  gwassig <- gwas %>%
    dplyr::filter(p<=pcut) %>%
    dplyr::filter(stringr::str_detect(string = rs_id, pattern = "^rs.*$"))
  # get ldmatrix, if not provided
  if(is.null(ldx)) {
    stopifnot(require(TwoSampleMR))
    if(nrow(gwassig)>500) {
      # only take 500 top snps
      warning(sprintf("Found %s signficant snps, trimming to 500 for LD computation ...", nrow(gwassig)))
      gwassig <- gwassig %>%
        dplyr::arrange(p) %>%
        dplyr::slice(500)
    }
    mysnps <- gwassig$rs_id
    ldx <- TwoSampleMR::ld_matrix(mysnps, with_alleles = TRUE, pop = "EUR")
  }
  # align snps in gwas table with ldx snps
  #else {
    # must have matching snps in gwassig and ldx
  ldsnps <- rownames(ldx)
  ldsnps <- gsub("_.*$", "", ldsnps)
  gwassig2 <- tibble::tibble(rs_id = ldsnps) %>%
    dplyr::inner_join(gwassig, by="rs_id")
  stopifnot(all(gwassig2$rs_id == ldsnps))
  if(nrow(gwassig2) < nrow(gwassig)) {
    warning(sprintf("Found ld info for %s snps, removing all others.", nrow(gwassig2)))
  }
  gwassig <- gwassig2
  rm(gwassig2)

  stopifnot(require(pheatmap))
  stopifnot(require(dendextend))
  pc <- pheatmap::pheatmap(ldx, show_colnames = FALSE, show_rownames = FALSE)
  # get hclust object
  hc <- pc$tree_row
  dd <- as.dendrogram(hc)

  # plot clusters with colored branches if you wish
  if(plot_clusters) {
    d1 <- color_branches(dd, k = k)
    plot(d1)
  }

  # choose 9 groups, best clustering, get labels for each group
  # then select best snp with lowest p for each cluster
  ct <- dendextend::cutree(dd, k = k)

  myleads <- list()
  for(ulx in unique(ct)) {
    # ulx <- unique(ct)[1]
    clsnpsall <- names(ct[ct == ulx])
    clsnps <- gsub("_.*$", "", clsnpsall)
    minp <- gwassig %>%
      dplyr::filter(rs_id %in% clsnps) %>%
      dplyr::group_by(chr) %>%
      dplyr::summarise(p = min(p))

    gwascl <- gwassig %>%
      dplyr::filter(rs_id %in% clsnps) %>%
      dplyr::filter(p %in% minp$p,
                    chr %in% minp$chr) %>%
      dplyr::arrange(p)
    xn <- paste0("chr",paste(unique(gwascl$chr), collapse="_"))
    myleads[[xn]] <- gwascl
  }

  return(list(ldmatrix = ldx, gwas_significant = gwassig, myleads = myleads))
}


find_closest_genes <- function(SNP_list,
                               plusi=300000,
                               plusf=300000,
                               genes = NULL,
                               debug=FALSE) { #ecg_31k_tab) {
  if(debug) {
    browser()
  }
  stopifnot(require(biomaRt))
  stopifnot(require(stringr))
  if(is.null(genes)) {
    Hs_Dataset = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host = "www.ensembl.org")
    attributes = listAttributes(Hs_Dataset)
    genes = getBM(attributes = c("chromosome_name","external_gene_name", "ensembl_gene_id", "start_position", "end_position", "description") ,mart = Hs_Dataset)
    genes$chromosome_name <- as.numeric(genes$chromosome_name) # please Ignore warning message
    genes <- genes[!is.na(genes$chromosome_name),]
    genes <- genes[order(genes$chromosome_name, genes$start_position),]
  }
  # Your SNP list
  #SNP_list <- ecg_31k_tab

  # Loop
  closest_genes <- list()
  # k <- 2
  for(k in 1:nrow(SNP_list)) {
    # kth row
    SNP <- SNP_list %>%
      dplyr::slice(k)
    cat("Chromosome_Number:", SNP$chr,"\n")
    gene = tibble::as_tibble(subset(genes, genes$chromosome_name == SNP$chr))

    # exactly one row
    posleft <- SNP$pos - plusi
    posright <- SNP$pos + plusf
    gx <- gene %>%
      dplyr::mutate(snp_pos = SNP$pos) %>%
      dplyr::filter(start_position >= posleft, # genes in range
                    end_position <= posright,
                    external_gene_name != "") %>%
      dplyr::filter(!stringr::str_detect(description, pattern="seudogene")) %>% # no pseudogenes
      dplyr::filter(!stringr::str_detect(external_gene_name, pattern="-AS1")) %>% # no antisense genes
      dplyr::mutate(distance_start = abs(snp_pos - start_position),
                    distance_end = abs(snp_pos - end_position)) %>%
      dplyr::mutate(distance_min = pmin(distance_start, distance_end)) %>%
      dplyr::arrange(distance_min) %>%
      dplyr::mutate(afterstart = snp_pos >= start_position,
                    beforeend = snp_pos <= end_position) %>%
      dplyr::mutate(ingene = afterstart & beforeend) %>%
      dplyr::mutate(rs_id = SNP$rs_id) %>%
      dplyr::relocate(rs_id, external_gene_name)
    if(any(gx$ingene)) {
      gx <- gx %>%
        dplyr::filter(ingene==TRUE) %>%
        dplyr::mutate(distance_min = as.character(distance_min)) %>% # need character to combine with non-ingene
        dplyr::group_by(rs_id, chromosome_name, ingene) %>%
        dplyr::summarise(external_gene_name = paste(unique(external_gene_name), collapse=","),
                         distance_min =  paste(unique(distance_min), collapse=","),
                         .groups = "drop") %>%
        dplyr::select(rs_id, external_gene_name, chromosome_name, distance_min, ingene)
    } else {
      # find next hit by distance, take top 2
      gx <- gx  %>%
        dplyr::slice_head(n=2) %>%
        dplyr::group_by(chromosome_name) %>%
        dplyr::summarise(external_gene_name = paste(external_gene_name, collapse=","),
                         distance_min = paste(distance_min, collapse=",")) %>%
        dplyr::mutate(rs_id = SNP$rs_id,
                      ingene = FALSE) %>% # all not in gene, but somewhere aside
        dplyr::select(rs_id, external_gene_name, chromosome_name, distance_min, ingene)
        #dplyr::relocate(rs_id, external_gene_name)
    }
    closest_genes[[SNP$rs_id]] <- gx
  }# m loop
  SNP_list <- SNP_list %>%
    dplyr::inner_join(dplyr::bind_rows(closest_genes),
                      by="rs_id") %>%
    dplyr::arrange(chr, external_gene_name, p, -abs(BETA)) %>%
    dplyr::select(-chromosome_name, -TEST) %>%
    dplyr::relocate(rs_id, external_gene_name, chr, pos, p, BETA, SE, ref, alt, T_STAT, A1, OBS_CT)

  return(SNP_list)
}





#######################################
#####    chromo Documentation    ######
#######################################

chromoDocumentation <- function() {
  file <- system.file("doc", "documentation.html", package = "chromo")
  if (file == "") stop("Documentation HTML not found. Try reinstalling the package.")
  utils::browseURL(file)
}


#################################
#####    chromo Vignette   ######
#################################

chromoVignette <- function() {
  file <- system.file("doc", "vignette.html", package = "chromo")
  if (file == "") stop("Vignette HTML not found. Try reinstalling the package.")
  utils::browseURL(file)
}


###############################
#####    chromo class    ######
###############################

setClass(
  Class = "Chromo",
  slots = list(
    data = "data.frame",
    columns = "ANY",
    classification = "ANY",
    genome = "ANY",
    composition = "ANY",
    density = "ANY",
    ora = "ANY",
    interactions = "ANY"
  )
)


################################
#####    chromoInitiate    #####
################################

chromoInitiate <- function(
    DEdf,
    pval_cutoff = 0.05,
    log2fc_cutoff = 1,
    gene_col,
    fc_col,
    p_col,
    keep_rep_feat = F,
    transcript_type = c("Mt_tRNA","Mt_rRNA","protein_coding","lncRNA","snRNA","rRNA","pseudogene","misc_RNA","processed_pseudogene",
                        "transcribed_unprocessed_pseudogene","rRNA_pseudogene","unprocessed_pseudogene","TEC","miRNA","transcribed_processed_pseudogene",
                        "snoRNA","unitary_pseudogene","transcribed_unitary_pseudogene","sRNA","IG_V_gene","IG_C_pseudogene","TR_J_gene",
                        "TR_V_gene","TR_J_pseudogene","TR_D_gene","TR_C_gene","IG_V_pseudogene","scaRNA","ribozyme","artifact","IG_C_gene",
                        "IG_J_gene","IG_J_pseudogene","IG_D_gene","translated_processed_pseudogene","TR_V_pseudogene","IG_pseudogene","vault_RNA","scRNA")
){

  all_features <- read.delim(system.file("extdata", "hsapiens_gene_ensembl.tsv", package = "chromo")) # from Ensembl
  cytobands <- read.delim(system.file("extdata", "hsapiens_cytogenicbands.tsv", package = "chromo")) # from UCSC's table browser

  chromoObject <- new("Chromo", genome = list(inp_annot = "hsapiens", cytobands = cytobands))

  # Annotate and filter
  DEdf <- DEdf %>%
    dplyr::left_join(all_features, by = stats::setNames("external_gene_name", gene_col)) %>%
    dplyr::filter(
      complete.cases(dplyr::select(., gene_col, "gene_biotype", "chromosome_name", "start_position", "end_position")),
      gene_biotype %in% transcript_type,
      chromosome_name %in% c(as.character(seq(1, 22)), "X", "Y", "MT")
    ) %>%
    {
      if (!keep_rep_feat) {
        dplyr::filter(., !duplicated(!!sym(gene_col)) | !duplicated(!!sym(gene_col), fromLast = TRUE))
      } else {.}
    } %>%
    dplyr::mutate(
      chromosome_name = factor(chromosome_name, levels = c(as.character(seq(1, 22)), "X", "Y", "MT")),
      gene_length = end_position - start_position,
      avg_position = (end_position + start_position)/2
    ) %>%
    dplyr::mutate(
      DEG = dplyr::case_when(
        !!rlang::sym(fc_col) > log2fc_cutoff & !!rlang::sym(p_col) < pval_cutoff ~ "UP",
        !!rlang::sym(fc_col) < -log2fc_cutoff & !!rlang::sym(p_col) < pval_cutoff ~ "DOWN",
        TRUE ~ "NO"
      )
    ) %>%
    dplyr::arrange(chromosome_name, start_position)

  chromoObject@columns <- list(
    gene_col = gene_col,
    fc_col = fc_col,
    p_col = p_col,
    chromosome = "chromosome_name",
    start_position = "start_position",
    end_position = "end_position",
    avg_position = "avg_position",
    gene_length = "gene_length",
    DEG = "DEG"
  )
  chromoObject@classification <- list(
    pval_cutoff = pval_cutoff,
    log2fc_cutoff = log2fc_cutoff
  )
  chromoObject@data <- DEdf

  return(chromoObject)
}


###################################
#####    chromoReclassify    ######
###################################

chromoReclassify <- function(
    chromoObject,
    pval_cutoff = 0.05,
    log2fc_cutoff = 1
){
  chromoObject@data <- chromoObject@data %>%
    mutate(
      DEG = case_when(
        !!sym(chromoObject@columns$fc_col) > log2fc_cutoff & !!sym(chromoObject@columns$p_col) < pval_cutoff ~ "UP",
        !!sym(chromoObject@columns$fc_col) < -log2fc_cutoff & !!sym(chromoObject@columns$p_col) < pval_cutoff ~ "DOWN",
        TRUE ~ "NO"
      )
    )

  chromoObject@classification <- list(pval_cutoff = pval_cutoff, log2fc_cutoff = log2fc_cutoff)

  return(chromoObject)
}


###################################
###      chromo Composition     ###
###################################

chromoComposition <- function(
    chromoObject,
    separate_by = "chromosome_name", #separate group
    only_expr_features = F,
    pct_expr_cols = c("pct.1", "pct.2"), # seurat
    score_method = "hyp_padj", # "pct", "hyp", "hyp_padj", "n_DEGs"
    padj_method = "BH" # existing params of p.adjust() function # only if using "hyp_padj"
){

  aux <- chromoObject@data %>%
    {if (only_expr_features) filter(., !!sym(pct_expr_cols[[1]]) + !!sym(pct_expr_cols[[2]]) != 0) else .}

  if(score_method %in% c("pct", "n_DEGs")){
    compo_df <- aux %>%
      group_by(!!sym(separate_by), DEG) %>%
      summarise(total = n()) %>%
      group_by(!!sym(separate_by)) %>%
      mutate(
        compo = (total / sum(total)) * 100
      ) %>%
      ungroup() %>%
      filter(
        DEG != "NO"
      )

    if(score_method == "n_DEGs"){
      compo_df$compo <- compo_df$total
    }

  }else if(score_method %in% c("hyp","hyp_padj")){
    compo_df <- aux %>%
      group_by(!!sym(separate_by), DEG) %>%
      summarise(total = n())

    totals <- list()
    totals$UP <- aux %>% filter(DEG == "UP") %>% nrow()
    totals$DOWN <- aux %>% filter(DEG == "DOWN") %>% nrow()
    totals$NO <- aux %>% filter(DEG == "NO") %>% nrow()

    # hypergeometric test
    compo_df$compo <- apply(compo_df, 1, function(x){
      phyper(
        as.numeric(x["total"]) - 1,
        totals[[x["DEG"]]],
        nrow(aux) - totals[[x["DEG"]]],
        sum(compo_df$total[compo_df[[separate_by]] == x[separate_by]]),
        lower.tail = FALSE
      )
    })

    compo_df <- compo_df %>%
      filter(DEG != "NO")

    if(score_method == "hyp_padj"){
      compo_df$compo <- p.adjust(compo_df$compo, method = padj_method)
    }
  }

  chromoObject@composition <- list(
    compo_df = compo_df,
    separate_by = separate_by,
    only_expr_features = only_expr_features,
    score_method = score_method
  )

  if (score_method == "hyp_padj") {
    chromoObject@composition$padj_method <- padj_method
  }
  if(only_expr_features){
    chromoObject@columns$pct_expr_cols <- pct_expr_cols
  }

  return(chromoObject)
}


#######################################
###      chromo Composition Plot    ###
#######################################

chromoCompositionPlot <- function(
    chromoObject,
    highlight_features = "top_by_separator", # "top_by_separator", "top_overall" or a list of names
    show_if_not_deg = T, # only if providing a list of names
    n_top_features = 1, # only if using "top_by_separator" or "top_overall"
    fc_line = T,
    title_xaxis = "Chromosome",
    title_yaxis = "Log2 fold change",

    color_dot_down = "#3771c8aa",
    color_dot_up = "#ff2200aa",
    color_dot_no = "#dddddd33",
    color_bar_down = "#3771c866",
    color_bar_up = "#ff220066",
    color_gene_name_down = "#2D5EAA",
    color_gene_name_up = "#aa0000",
    color_gene_name_no = "#555555",
    color_score_down = "#2D5EAA",
    color_score_up = "#aa0000",
    color_line_up = "#aa0000",
    color_line_down = "#2D5EAA",
    color_xaxis_text = "black",
    color_xaxis_label = "black",
    color_yaxis_text = "black",
    color_yaxis_label = "black",

    size_dot_alt = 1.2,
    size_dot_no = 0.8,
    size_bar = 0.8,
    size_gene_name = 3.2,
    size_score = ifelse(chromoObject@composition$score_method %in% c("hyp", "hyp_padj"), 7, 4),
    size_line = 0.4,
    size_xaxis_text = 14,
    size_xaxis_label = 12,
    size_yaxis_text = 12,
    size_yaxis_label = 12,

    style_gene_name = "plain",
    style_score = "bold",
    style_line = 2, # 1 - continuous, 2- dashed, 3 - dotted, etc.
    style_xaxis_text = "bold",
    style_xaxis_label = "plain",
    style_yaxis_text = "bold",
    style_yaxis_label = "plain"
){

  gene_col <- chromoObject@columns$gene_col
  fc_col <- chromoObject@columns$fc_col
  separate_by <- chromoObject@composition$separate_by
  compo_df <- chromoObject@composition$compo_df
  pct_expr_cols <- chromoObject@columns$pct_expr_cols

  aux <- chromoObject@data %>%
    {if (chromoObject@composition$only_expr_features) filter(., !!sym(pct_expr_cols[[1]]) + !!sym(pct_expr_cols[[2]]) != 0) else .}

  compo_df <- compo_df %>%
    mutate(
      compo = case_when(
        compo == 0 ~ min(compo[compo > 0], na.rm = TRUE),
        TRUE ~ compo
      ),
      compo = case_when(
        chromoObject@composition$score_method %in% c("hyp","hyp_padj") ~ -log10(compo),
        TRUE ~ compo
      )
    )

  max_compo <- max(compo_df$compo, na.rm = T)

  compo_df <- compo_df %>%
    mutate(
      proportion = (compo/max_compo * (max(abs(aux[[fc_col]])))),
      proportion = case_when(
        DEG == "UP" ~ proportion,
        DEG == "DOWN" ~ -proportion,
        TRUE ~ proportion
      ),
      compo = case_when(
        chromoObject@composition$score_method %in% c("hyp","hyp_padj") ~ ifelse(compo > -log10(0.001),"***",ifelse(compo > -log10(0.01),"**",ifelse(compo > -log10(0.05),"*",""))),
        TRUE ~ paste0(round(compo, 1), "%")
      ),
      y_axis = case_when(
        DEG == "UP" ~ 1.2 * max(aux[[fc_col]]),
        DEG == "DOWN" ~ 1.2 * min(aux[[fc_col]]),
        TRUE ~ NA
      ),
      color = case_when(
        DEG == "UP" ~ color_score_up,
        DEG == "DOWN" ~ color_score_down,
        TRUE ~ NA
      )
    )

  # highlight features
  if(highlight_features %in% c("top_by_separator", "top_overall")){
    max_genes <- aux %>%
      filter(DEG == "UP") %>%
      { if (highlight_features == "top_by_separator") group_by(., !!sym(separate_by)) else . } %>%
      slice_max(order_by = !!sym(fc_col), n = n_top_features, with_ties = FALSE) %>%
      { if (highlight_features == "top_by_separator") ungroup(.) else . } %>%
      dplyr::select(!!sym(gene_col), !!sym(separate_by))

    min_genes <- aux %>%
      filter(DEG == "DOWN") %>%
      { if (highlight_features == "top_by_separator") group_by(., !!sym(separate_by)) else . } %>%
      slice_min(order_by = !!sym(fc_col), n = n_top_features, with_ties = FALSE) %>%
      { if (highlight_features == "top_by_separator") ungroup(.) else . } %>%
      dplyr::select(!!sym(gene_col), !!sym(separate_by))

    hlf <- rbind(max_genes, min_genes) %>%
      mutate(gene_sep = paste0(!!sym(gene_col), "_", !!sym(separate_by))) %>%
      pull(gene_sep)
  }else{
    hlf <- highlight_features
  }

  local_aux <- aux %>%
    mutate(
      highlight = case_when(
        highlight_features %in% c("top_by_separator", "top_overall") & paste0(!!sym(gene_col), "_", !!sym(separate_by)) %in% hlf ~ !!sym(gene_col),
        !highlight_features %in% c("top_by_separator", "top_overall") & show_if_not_deg & !!sym(gene_col) %in% hlf ~ !!sym(gene_col),
        !highlight_features %in% c("top_by_separator", "top_overall") & !show_if_not_deg & !!sym(gene_col) %in% hlf & DEG %in% c("UP", "DOWN") ~ !!sym(gene_col),
        TRUE ~ NA
      ),
      color = case_when(
        !is.na(highlight) & DEG == "UP" ~ color_gene_name_up,
        !is.na(highlight) & DEG == "DOWN" ~ color_gene_name_down,
        !is.na(highlight) & DEG == "NO" ~ color_gene_name_no,
        TRUE ~ NA
      )
    )

  set.seed(42) # gene name position

  deg_plot <- ggplot()+
    scale_x_discrete(limits = levels(chromoObject@data[[separate_by]])) + # fixes order of separate_by levels. NOTE: requires user to change separate_by into factor!
    geom_bar(
      data = compo_df,
      aes(x = !!sym(separate_by), y = proportion, fill = DEG),
      stat = "identity",
      width = size_bar
    ) +
    scale_fill_manual(values = c("UP" = color_bar_up, "DOWN" = color_bar_down))

  if(fc_line){
    deg_plot <- deg_plot +
      geom_hline(yintercept = chromoObject@classification$log2fc_cutoff, color = color_line_up, linewidth = size_line, linetype = style_line) +
      geom_hline(yintercept = -chromoObject@classification$log2fc_cutoff, color = color_line_down, linewidth = size_line, linetype = style_line)
  }

  deg_plot <- deg_plot +
    geom_jitter( # not DEGs
      data = local_aux %>% filter(DEG == "NO"),
      aes(x = !!sym(separate_by), y = !!sym(fc_col), color = DEG),
      position = position_jitterdodge(
        jitter.width = 0.4,
        jitter.height = 0,
        dodge.width = 0,
        seed = 42
      ),
      size = size_dot_no
    )+
    geom_jitter( # DEGs
      data = local_aux %>% filter(DEG %in% c("UP", "DOWN")),
      aes(x = !!sym(separate_by), y = !!sym(fc_col), color = DEG),
      position = position_jitterdodge(
        jitter.width = 0.4,
        jitter.height = 0,
        dodge.width = 0,
        seed = 42
      ),
      size = size_dot_alt
    )+
    labs(
      x = title_xaxis,
      y = title_yaxis
    )+
    scale_color_manual(values = c("DOWN" = color_dot_down, "NO" = color_dot_no, "UP" = color_dot_up)) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.x = element_text(color = color_xaxis_text, size = size_xaxis_text, face = style_xaxis_text),
      axis.title.x = element_text(color = color_xaxis_label, size = size_xaxis_label, face = style_xaxis_label),
      axis.text.y = element_text(color = color_yaxis_text, size = size_yaxis_text, face = style_yaxis_text),
      axis.title.y = element_text(color = color_yaxis_label, size = size_yaxis_label, face = style_yaxis_label),
      legend.position = "none"
    )+
    geom_text_repel(
      data = local_aux,
      aes(x = !!sym(separate_by), y = !!sym(fc_col), color = DEG, label = highlight),
      max.overlaps = Inf,
      color = local_aux$color,
      size = size_gene_name,
      fontface = style_gene_name
    ) +
    annotate( #score
      geom = "text",
      x = compo_df[[separate_by]],
      y = compo_df$y_axis,
      label = compo_df$compo,
      color = compo_df$color,
      size = size_score,
      fontface = style_score
    )

  return(deg_plot)
}


#################################
#####    chromo Density    ######
#################################

chromoDensity <- function(
    chromoObject,
    bandwidth = "nrd0",
    cluster_threshold = 20, # 20%
    DEG_type = list("UP", "DOWN"), # c("UP", "DOWN"), "UP, "DOWN"
    padj_method = "none" # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", fdr", "none")
){

  calculate_density <- function(subset, chro) {
    dens <- density(
      subset[[avg_position]],
      kernel = "gaussian",
      from = 1,
      to = max(cytobands[cytobands$chr == chro,"baseEnd"]),
      bw = bandwidth
    )
    return(data.frame(x = dens$x, y = dens$y))
  }

  gene_col <- chromoObject@columns$gene_col
  fc_col <- chromoObject@columns$fc_col
  chromosome <- chromoObject@columns$chromosome
  start_position <- chromoObject@columns$start_position
  end_position <- chromoObject@columns$end_position
  avg_position <- chromoObject@columns$avg_position
  gene_length <- chromoObject@columns$gene_length
  DEG <- chromoObject@columns$DEG
  cytobands <- chromoObject@genome$cytobands
  DEdf <- chromoObject@data %>%
    filter(!!sym(chromosome) != "MT") %>%
    mutate(!!sym(chromosome) := factor(!!sym(chromosome), levels = c(as.character(seq(1, 22)), "X", "Y")))

  if(!is.list(DEG_type)){
    DEG_type <- list(DEG_type)
  }

  # DEGs density
  DEG_density <- DEdf %>%
    filter(!!sym(DEG) %in% DEG_type) %>%
    group_by(!!sym(chromosome)) %>%
    do({
      dens_df <- calculate_density(., chro = unique(.data[[chromosome]]))
      dens_df
    }) %>%
    mutate(
      y = if_else(y < max(y)*(cluster_threshold/100), 0, y)
    ) %>%
    ungroup() %>%
    mutate(!!sym(chromosome) := factor(!!sym(chromosome)))

  # adding 0s on borders
  for (i in levels(DEG_density[[chromosome]])) {
    DEG_density <- DEG_density %>%
      add_row(!!sym(chromosome) := i, x = 0, y = 0) %>%
      add_row(!!sym(chromosome) := i, x = max(cytobands[cytobands$chr == i,"baseEnd"]), y = 0)
  }

  # ordering
  DEG_density <- DEG_density %>%
    mutate(!!sym(chromosome) := factor(!!sym(chromosome), levels = levels(DEdf[[chromosome]]))) %>%
    arrange(!!sym(chromosome), x)

  # Clustering
  DEG_clusters <- DEG_density[DEG_density$y != 0 | (DEG_density$y == 0 & (c(TRUE, DEG_density$y[-length(DEG_density$y)] != 0) | c(DEG_density$y[-1] != 0, TRUE))), , drop = FALSE]
  DEG_clusters <- as.data.frame(DEG_clusters)
  if(DEG_clusters$y[1] == 0 & DEG_clusters$y[2] == 0){DEG_clusters <- DEG_clusters[-1,]}
  if(DEG_clusters$y[nrow(DEG_clusters)] == 0 & DEG_clusters$y[nrow(DEG_clusters)-1] == 0){DEG_clusters <- DEG_clusters[-nrow(DEG_clusters),]}

  # list
  DEG_density_list <- list()
  boundaries <- which(DEG_clusters$y == 0)
  for (i in seq(1, length(boundaries), by=2)) {
    DEG_density_list[[(i+1)/2]] <- DEG_clusters[boundaries[i]:boundaries[i + 1], ]
  }

  # remaking clusters df
  DEG_clusters <- map2_dfr(DEG_density_list, seq_along(DEG_density_list), ~ {
    df <- .x
    cluster_num <- .y
    data.frame(
      chromosome = unique(df[[chromosome]]),
      cluster_num = cluster_num,
      end_position = max(df$x),
      start_position = min(df$x)
    )
  }) %>%
    mutate(size = end_position - start_position)

  if(nrow(DEG_clusters) > 0){

    DEG_clusters$all_features <- NA
    DEG_clusters$DEGs <- NA
    for (i in 1:nrow(DEG_clusters)) {
      DEG_clusters$all_features[i] <- DEdf %>%
        filter(!!sym(chromosome) == DEG_clusters$chromosome[i], !!sym(avg_position) < DEG_clusters$end_position[i], !!sym(avg_position) > DEG_clusters$start_position[i]) %>%
        pull(!!sym(gene_col)) %>%
        paste(collapse = ";")
      DEG_clusters$DEGs[i] <- DEdf %>%
        filter(!!sym(DEG) %in% DEG_type) %>%
        filter(!!sym(chromosome) == DEG_clusters$chromosome[i], !!sym(avg_position) < DEG_clusters$end_position[i], !!sym(avg_position) > DEG_clusters$start_position[i]) %>%
        pull(!!sym(gene_col)) %>%
        paste(collapse = ";")
    }

    DEG_clusters <- DEG_clusters %>%
      mutate(
        n_features = str_count(all_features, ";") + 1,
        n_DEG = str_count(DEGs, ";") + 1
      ) %>%
      filter(n_DEG > 1) # removing clusters with only 1 DEG

    if(nrow(DEG_clusters) > 0){
      # Hypergeometric test
      total_DEG <- DEdf %>% filter(!!sym(DEG) %in% DEG_type) %>% nrow()
      total_no_DEG <- DEdf %>% filter(!!sym(DEG) == "NO") %>% nrow()

      DEG_clusters$pval <- apply(DEG_clusters, 1, function(x){
        phyper(
          as.numeric(x[["n_DEG"]]) - 1,
          total_DEG,
          total_no_DEG,
          as.numeric(x[["n_features"]]),
          lower.tail = FALSE
        )
      })

      # padj and score
      DEG_clusters <- DEG_clusters %>%
        mutate(
          pval = p.adjust(DEG_clusters$pval, method = padj_method),
          pval = case_when(
            pval == 0 ~ min(pval[pval > 0], na.rm = TRUE),
            TRUE ~ pval
          ),
          score = -log10(pval)
        ) %>%
        arrange(-score)

      # Numbering clusters
      DEG_clusters <- DEG_clusters %>%
        mutate(cluster_num = seq(1, nrow(DEG_clusters)))

      # Bands affected by each cluster
      DEG_clusters <- DEG_clusters %>%
        rowwise() %>%
        mutate(
          bands = {
            chr_val   <- as.character(.data[["chromosome"]])
            start_val <- as.numeric(.data[[start_position]])
            end_val   <- as.numeric(.data[[end_position]])

            aux <- cytobands %>%
              filter(
                chr == chr_val,
                (baseStart <= end_val & baseStart >= start_val) |
                  (baseEnd <= end_val & baseEnd >= start_val) |
                  (baseStart <= start_val & baseEnd >= end_val)
              )

            paste(aux$band, collapse = ";")
          }
        ) %>%
        ungroup() %>%
        dplyr::select(cluster_num, chromosome, start_position, end_position, size, all_features, DEGs, n_features, n_DEG, pval, score, bands)
    }
  }

  chromoObject@density[[paste(DEG_type, collapse = ifelse(length(DEG_type) > 1, "_", ""))]] = list(
    DEG_clusters = DEG_clusters,
    bandwidth = bandwidth,
    threshold = cluster_threshold,
    padj_method = padj_method
  )

  return(chromoObject)
}


#####################################
#####    chromo Density Plot   ######
#####################################

chromoDensityPlot <- function(
    chromoObject,
    DEG_type = list("UP", "DOWN"), # list("UP", "DOWN"), "UP", "DOWN"
    n_top_clusters = nrow(DEG_clusters),
    include_genes = F,
    include_density = T,
    all_chr = T,

    color_enrich = "#990099",
    color_enrich_up = "#dd2200",
    color_enrich_down = "#0022dd"
){

  gene_col <- chromoObject@columns$gene_col
  fc_col <- chromoObject@columns$fc_col
  chromosome <- chromoObject@columns$chromosome
  start_position <- chromoObject@columns$start_position
  end_position <- chromoObject@columns$end_position
  avg_position <- chromoObject@columns$avg_position
  gene_length <- chromoObject@columns$gene_length
  DEG <- chromoObject@columns$DEG
  cytobands <- chromoObject@genome$cytobands
  DEdf <- chromoObject@data %>%
    filter(!!sym(chromosome) != "MT") %>%
    mutate(!!sym(chromosome) := factor(!!sym(chromosome), levels = c(as.character(seq(1, 22)), "X", "Y")))

  if(!is.list(DEG_type)){
    DEG_type <- list(DEG_type)
  }

  density_type <- paste(DEG_type, collapse = ifelse(length(DEG_type) > 1, "_", ""))

  plot_color <- ifelse(length(DEG_type) == 2, color_enrich, ifelse("UP" %in% DEG_type, color_enrich_up, color_enrich_down))

  DEG_clusters <- chromoObject@density[[density_type]][["DEG_clusters"]]

  if(nrow(DEG_clusters) == 0){
    return(ggplot() + theme_classic()) # return empty ggplot
  }

  # top clusters and bands to include
  top_clusters <- DEG_clusters %>% arrange(-score) %>% head(n_top_clusters)

  if(all_chr){
    chr_with_clu <- c(as.character(seq(1,22)), "X", "Y")
  }else{
    chr_with_clu <- unique(top_clusters$chromosome)
  }

  bands_to_keep <- cytobands %>%
    group_by(chr) %>%
    summarize(
      baseStart = min(baseStart),
      baseEnd   = max(baseEnd),
      .groups   = "drop"
    ) %>%
    mutate(
      band = row_number()
    ) %>%
    rbind(cytobands %>% filter(sub("^[^_]*_", "", band) %in% c("p11.1", "p11", "q11.1", "q11")))

  # plot
  plots_list <- c()
  for(i in 1:length(chr_with_clu)){

    # plotting
    plots_list[[i]] <- ggplot() +
      labs(
        title = NULL,
        x = NULL,
        y = paste0("chr", chr_with_clu[i])) +
      theme_minimal() +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, max(as.numeric(cytobands %>% filter(chr %in% chr_with_clu) %>% pull(baseEnd))))) +
      coord_cartesian(clip = "off")+
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.title.y = element_text(face = "bold", size = 20, angle = 0, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank()
      )

    if(include_density & chr_with_clu[i] %in% top_clusters$chromosome){

      dens <- density(
        DEdf %>% filter(!!sym(chromosome) == chr_with_clu[i], DEG %in% DEG_type) %>% pull(avg_position),
        kernel = "gaussian",
        from = 1,
        to = max(cytobands[cytobands$chr == chr_with_clu[i], "baseEnd"]),
        bw = chromoObject@density[[density_type]]$bandwidth
      )
      dens_deg <- data.frame(x = dens$x, y_deg = dens$y) %>%
        mutate(y_deg = y_deg/max(y_deg))

      dens <- density(
        DEdf %>% filter(!!sym(chromosome) == chr_with_clu[i]) %>% pull(avg_position),
        kernel = "gaussian",
        from = 1,
        to = max(cytobands[cytobands$chr == chr_with_clu[i], "baseEnd"]),
        bw = chromoObject@density[[density_type]]$bandwidth
      )
      dens_all <- data.frame(x = dens$x, y_all = dens$y) %>%
        mutate(y_all = y_all/max(y_all))

      dens <- dens_deg %>%
        left_join(dens_all, by = "x") %>%
        mutate(
          y = y_deg - y_all,
          y = case_when(
            y < 0 ~ 0,
            TRUE ~ y
          ),
          y = 1.5*y/max(y)
        ) %>%
        add_row(x = 0, y = 0) %>%
        add_row(x = max(cytobands[cytobands$chr == chr_with_clu[i], "baseEnd"]), y = 0) %>%
        arrange(x)

      plots_list[[i]] <- plots_list[[i]] +
        geom_polygon(
          data = dens,
          aes(x = x, y = y),
          color = plot_color,
          linewidth = 0.5,
          fill = paste0(plot_color, "77")
        )
    }

    # fixes chr sizes
    if(include_density & !(chr_with_clu[i] %in% top_clusters$chromosome)){
      plots_list[[i]] <- plots_list[[i]] +
        geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1.5), alpha = 0)
    }

    if(include_genes){
      plots_list[[i]] <- plots_list[[i]] +
        geom_point( # not DEGs
          data = DEdf %>% filter(!!sym(chromosome) == chr_with_clu[i], !!sym(DEG) == "NO"),
          aes(x = !!sym(avg_position)), y = 0.1, color = "#00000022"
        ) +
        geom_point( # DEGs
          data = DEdf %>% filter(!!sym(chromosome) == chr_with_clu[i], !!sym(DEG) %in% DEG_type),
          aes(x = !!sym(avg_position)), y = 0.1, color = paste0(plot_color, "aa")
        )
    }

    # cytogenetic bands
    for (j in bands_to_keep[bands_to_keep$chr == chr_with_clu[i], "band", drop = TRUE]) {

      plots_list[[i]] <- plots_list[[i]] +
        annotate(
          "rect",
          xmin = bands_to_keep[bands_to_keep$band == j, "baseStart", drop = T],
          xmax = bands_to_keep[bands_to_keep$band == j, "baseEnd", drop = T],
          ymin = -1,
          ymax = 0,
          fill = ifelse(sub("^[^_]*_", "", j) %in% c("p11.1", "p11", "q11.1", "q11"), "black", "#ffffff"),
          color = "black"
        )
    }

    # clusters
    for (j in top_clusters[top_clusters$chromosome == chr_with_clu[i], "cluster_num", drop = T]) {

      plots_list[[i]] <- plots_list[[i]] +
        annotate(
          "rect",
          xmin = top_clusters[top_clusters$cluster_num == j, "start_position", drop = T],
          xmax = top_clusters[top_clusters$cluster_num == j, "end_position", drop = T],
          ymin = -1,
          ymax = 0,
          fill = colorRampPalette(c("white", plot_color))(256)[
            as.integer(rescale(
              top_clusters[top_clusters$cluster_num == j, "score", drop = TRUE],
              to = c(0, 1),
              from = c(min(top_clusters$score), max(top_clusters$score))
            ) * 255) + 1
          ],
          color = plot_color,
          size = 0.8,
          alpha = 0.8
        )
    }

    # cluster name
    for (j in top_clusters[top_clusters$chromosome == chr_with_clu[i], "cluster_num", drop = T]) {
      plots_list[[i]] <- plots_list[[i]] +
        annotate(
          "text",
          x = (top_clusters[top_clusters$cluster_num == j,"end_position", drop = T] + top_clusters[top_clusters$cluster_num == j,"start_position", drop = T])/2,
          y = -0.5,
          label = j,
          color = "black",
          size = 6,
          fontface = "bold"
        )
    }
  }

  combined_plot <- purrr::reduce(plots_list, `/`)

  return(combined_plot)
}


###############################
#####     chromoZoom     ######
###############################

chromoZoom <- function(
    chromoObject,
    DEG_type = list("UP", "DOWN"), # list("UP", "DOWN"), "UP", "DOWN"
    mt_density = F,
    cluster = 1,

    color_line_up = "#aa0000",
    color_line_down = "#2D5EAA",
    color_xaxis_text = "black",
    color_xaxis_label = "black",
    color_yaxis_text = "black",
    color_yaxis_label = "black",

    size_gene_name = 3,
    size_line = 0.4,
    size_xaxis_text = 12,
    size_xaxis_label = 12,
    size_yaxis_text = 12,
    size_yaxis_label = 12,

    style_gene_name = "bold",
    style_line = 2, # 1 - continuous, 2- dashed, 3 - dotted, etc.
    style_xaxis_text = "bold",
    style_xaxis_label = "plain",
    style_yaxis_text = "bold",
    style_yaxis_label = "plain"
){

  gene_col <- chromoObject@columns$gene_col
  fc_col <- chromoObject@columns$fc_col
  chromosome <- chromoObject@columns$chromosome
  start_position <- chromoObject@columns$start_position
  end_position <- chromoObject@columns$end_position
  avg_position <- chromoObject@columns$avg_position
  gene_length <- chromoObject@columns$gene_length
  DEG <- chromoObject@columns$DEG
  cytobands <- chromoObject@genome$cytobands

  custom_labels <- function(x) {
    ifelse(x >= 1e9, paste0(round(x / 1e9, 2), " Gb"),
           ifelse(x >= 1e6, paste0(round(x / 1e6, 2), " Mb"),
                  ifelse(x >= 1e3, paste0(round(x / 1e3, 2), " kb"), as.character(x))))
  }

  DEdf <- chromoObject@data

  if(!is.list(DEG_type)){
    DEG_type <- list(DEG_type)
  }

  density_type <- paste(DEG_type, collapse = ifelse(length(DEG_type) > 1, "_", ""))

  if(mt_density){
    cluster_df <- chromoObject@density$MT[[density_type]]$DEG_clusters %>% filter(cluster_num == cluster)
  }else{
    cluster_df <- chromoObject@density[[density_type]]$DEG_clusters %>% filter(cluster_num == cluster)
  }

  features_in_cluster <- unlist(strsplit(cluster_df[["all_features"]], ";"))
  DEGs_in_cluster <- unlist(strsplit(cluster_df[["DEGs"]], ";"))
  not_DEGs <- setdiff(features_in_cluster, DEGs_in_cluster)

  # min and max log2fc
  fc_vector <- DEdf %>% filter(!!sym(gene_col) %in% features_in_cluster) %>% pull(!!sym(fc_col))
  max_fc <- max(fc_vector)
  min_fc <- min(fc_vector)
  size_adj <- (max_fc - min_fc)/30

  zoom_plot <- ggplot() +
    labs(
      title = paste0(
        "Density: ", density_type,
        ", Cluster: ", cluster,
        ", Chr", cluster_df[["chromosome"]],
        ", ", custom_labels(cluster_df[["start_position"]]), " - ", custom_labels(cluster_df[["end_position"]])
      ),
      x = NULL,
      y = "Log2 fold change"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(min_fc - 4*size_adj, max_fc), breaks = seq(ceiling(min_fc), floor(max_fc), 1)) +
    scale_x_continuous(expand = c(0, 0), labels = custom_labels) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.text.x = element_text(color = color_xaxis_text, size = size_xaxis_text, face = style_xaxis_text),
      axis.title.x = element_text(color = color_xaxis_label, size = size_xaxis_label, face = style_xaxis_label),
      axis.text.y = element_text(color = color_yaxis_text, size = size_yaxis_text, face = style_yaxis_text),
      axis.title.y = element_text(color = color_yaxis_label, size = size_yaxis_label, face = style_yaxis_label),
      legend.position = "none" # Remove legends
    ) +
    geom_vline(xintercept = cluster_df[["start_position"]], color = "#dddd00", size = size_line, linetype = style_line) +
    geom_vline(xintercept = cluster_df[["end_position"]], color = "#dddd00", size = size_line, linetype = style_line) +
    geom_hline(yintercept = chromoObject@classification$log2fc_cutoff, color = color_line_up, size = size_line, linetype = style_line) +
    geom_hline(yintercept = -chromoObject@classification$log2fc_cutoff, color = color_line_down, size = size_line, linetype = style_line) +
    geom_hline(yintercept = 0, color = "#444444", size = size_line, linetype = style_line)

  # Bands
  if(!mt_density){
    bands <- unlist(strsplit(cluster_df[["bands"]], ";"))
    for(i in bands){
      zoom_plot <- zoom_plot +
        annotate(
          "rect",
          xmin = cytobands[cytobands$band == i, "baseStart"],
          xmax = cytobands[cytobands$band == i, "baseEnd"],
          ymin = min_fc - 4*size_adj,
          ymax = min_fc - 1.5*size_adj,
          fill = "#ffffff00",
          color = "black"
        ) +
        annotate(
          "text",
          x = (cytobands[cytobands$band == i, "baseStart"] + cytobands[cytobands$band == i, "baseEnd"])/2,
          y = min_fc - 2.75*size_adj,
          label = sub("^[^_]*_", "", i),
          color = "black",
          size = 3.5,
          fontface = "bold"
        )
    }
  }

  # not DEGs
  for (i in not_DEGs){
    zoom_plot <- zoom_plot +
      annotate("rect",
               xmin = DEdf %>% filter(!!sym(gene_col) == i) %>% pull(!!sym(start_position)),
               xmax = DEdf %>% filter(!!sym(gene_col) == i) %>% pull(!!sym(end_position)),
               ymin = DEdf %>% filter(!!sym(gene_col) == i) %>% pull(!!sym(fc_col)) - size_adj,
               ymax = DEdf %>% filter(!!sym(gene_col) == i) %>% pull(!!sym(fc_col)),
               fill = '#77777777'
      )
  }

  # DEGs
  for (i in DEGs_in_cluster){
    zoom_plot <- zoom_plot +
      annotate("rect",
               xmin = DEdf %>% filter(!!sym(gene_col) == i) %>% pull(!!sym(start_position)),
               xmax = DEdf %>% filter(!!sym(gene_col) == i) %>% pull(!!sym(end_position)),
               ymin = DEdf %>% filter(!!sym(gene_col) == i) %>% pull(!!sym(fc_col)) - size_adj,
               ymax = DEdf %>% filter(!!sym(gene_col) == i) %>% pull(!!sym(fc_col)),
               fill = ifelse(DEdf %>% filter(!!sym(gene_col) == i) %>% pull(!!sym(DEG)) == "DOWN", "#0033ffaa", "#ff3300aa")
      )
  }

  set.seed(42)

  # gene name
  zoom_plot <- zoom_plot +
    geom_text_repel(
      data = DEdf %>% filter(!!sym(gene_col) %in% DEGs_in_cluster),
      aes(x = avg_position, y = !!sym(fc_col), label = !!sym(gene_col), color = DEG), ######
      hjust = 0.5,
      vjust = 0.5,
      size = size_gene_name,
      fontface = style_gene_name,
      inherit.aes = F
    )+
    scale_color_manual(values = c("DOWN" = "#002277", "UP" = "#772200"))

  return(zoom_plot)
}


##############################
######    chromo ORA   #######
##############################

chromoORA <- function(
    chromoObject,
    DEG_type = list("UP", "DOWN"), # list("UP", "DOWN"), "UP", "DOWN"
    cluster = 1,
    ont_type = "BP", # "BP", "CC", "MF", "ALL"
    mt_density = F
){

  if(!is.list(DEG_type)){
    DEG_type <- list(DEG_type)
  }

  density_type <- paste(DEG_type, collapse = ifelse(length(DEG_type) > 1, "_", ""))

  if(mt_density){
    DEGs_in_cluster <- chromoObject@density$MT[[density_type]][["DEG_clusters"]] %>%
      filter(cluster_num == cluster) %>%
      pull(DEGs)
  }else{
    DEGs_in_cluster <- chromoObject@density[[density_type]][["DEG_clusters"]] %>%
      filter(cluster_num == cluster) %>%
      pull(DEGs)
  }

  DEGs_in_cluster <- strsplit(DEGs_in_cluster, ";")[[1]]

  aux <- enrichGO(
    keyType = "SYMBOL",
    gene = DEGs_in_cluster,
    universe = chromoObject@data[, chromoObject@columns$gene_col, drop = T],
    OrgDb = org.Hs.eg.db,
    ont = ont_type,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  aux <- aux@result %>% mutate(score = -log10(p.adjust))

  if(ont_type != "ALL"){
    aux$ONTOLOGY <- ont_type
  }

  if(mt_density){
    chromoObject@ora$MT[[density_type]][[cluster]] <- aux
  }else{
    chromoObject@ora[[density_type]][[cluster]] <- aux
  }

  return(chromoObject)
}


###################################
######    chromo ORA Plot   #######
###################################

chromoORAPlot <- function(
    chromoObject,
    DEG_type = list("UP", "DOWN"), # list("UP", "DOWN"), "UP", "DOWN"
    cluster = 1,
    mt_density = F,
    number_of_onto = 5,
    highlight = "ONTOLOGY",
    color_list = c("CC" = "#FBB4AEaa", "MF"  = "#B3CDE3aa", "BP" = "#CCEBC5aa"),
    text_size = 4,
    title = paste0(DEG_type, "_cluster: ", cluster)
){

  if(!is.list(DEG_type)){
    DEG_type <- list(DEG_type)
  }

  density_type <- paste(DEG_type, collapse = ifelse(length(DEG_type) > 1, "_", ""))

  if(mt_density){
    df <- chromoObject@ora$MT[[density_type]][[cluster]] %>%
      arrange(-score) %>%
      head(number_of_onto)
  }else{
    df <- chromoObject@ora[[density_type]][[cluster]] %>%
      arrange(-score) %>%
      head(number_of_onto)
  }

  df <- df %>%
    mutate(
      Description = factor(Description, levels = unique(df %>% arrange(score) %>% pull(Description)))
    )

  bar_plot <- ggplot(df) +
    geom_col(
      aes(x = score, y = Description, fill = !!sym(highlight)),
      width = 0.6
    ) +
    geom_vline(
      xintercept = -log10(0.05),
      linetype = "dashed",
      color = "#cc2222"
    ) +
    geom_vline(
      xintercept = -log10(0.01),
      linetype = "dashed",
      color = "#cc2222"
    ) +
    geom_vline(
      xintercept = -log10(0.001),
      linetype = "dashed",
      color = "#cc2222"
    ) +
    annotate(
      "text",
      x = -log10(0.05),
      y = nrow(df) + 0.6,
      label = "*",
      color = "#cc2222",
      size = 6
    )+
    annotate(
      "text",
      x = -log10(0.01),
      y = nrow(df) + 0.6,
      label = "**",
      color = "#cc2222",
      size = 6
    )+
    annotate(
      "text",
      x = -log10(0.001),
      y = nrow(df) + 0.6,
      label = "***",
      color = "#cc2222",
      size = 6
    )+
    labs(
      title = title,
      y = NULL,
      x = expression("-log"[10]*"(padjust)")
    ) +
    scale_x_continuous(
      limits = c(0, max(df$score)),
      breaks = seq(0, max(df$score), by = 1),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(0, "mm"),
      axis.line.y.left = element_line(color = "black"),
      axis.text = element_blank(),
      axis.title.x = element_text(size = 12)
    ) +
    scale_fill_manual(
      values = color_list
    ) +
    geom_text(
      aes(0, y = Description, label = Description),
      hjust = 0,
      nudge_x = 0.05,
      color = "black",
      size = text_size
    )

  return(bar_plot)
}


# ###################################
# #####    chromo Density MT   ######
# ###################################
#
# chromoDensityMT <- function(
#     chromoObject,
#     bandwidth = 1.06 * sd(DEG_density) * length(DEG_density)^(-1/5)/4, # Silverman's rule: nrd0
#     cluster_threshold = 20, # 20%
#     DEG_type = list("UP", "DOWN"), # c("UP", "DOWN"), "UP, "DOWN"
#     padj_method = "none" # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", fdr", "none")
# ){
#
#   gene_col <- chromoObject@columns$gene_col
#   fc_col <- chromoObject@columns$fc_col
#   chromosome <- chromoObject@columns$chromosome
#   start_position <- chromoObject@columns$start_position
#   end_position <- chromoObject@columns$end_position
#   avg_position <- chromoObject@columns$avg_position
#   gene_length <- chromoObject@columns$gene_length
#   DEG <- chromoObject@columns$DEG
#   DEdf <- chromoObject@data %>%
#     filter(!!sym(chromosome) == "MT")
#
#   if(nrow(DEdf) == 0){
#     stop(paste0("There are no ", paste(DEG_type, collapse = ifelse(length(DEG_type) > 1, "_", "")), " in MT! Try using reclassify."))
#   }
#
#   if(!is.list(DEG_type)){
#     DEG_type <- list(DEG_type)
#   }
#
#   DEG_density <- DEdf %>%
#     filter(!!sym(DEG) %in% DEG_type) %>%
#     pull(!!sym(avg_position))
#
#   DEG_density <- (DEG_density / 16023) * 2 * pi
#   DEG_density <- circular(DEG_density, units = "radians", modulo = "asis")
#
#   # Kernel Density
#   DEG_density <- density.circular(
#     DEG_density,
#     bw = bandwidth,
#     kernel = "wrappednormal"
#   )
#
#   DEG_density <- data.frame(
#     x = DEG_density$x,
#     y = DEG_density$y
#   )
#
#   # clustering
#   DEG_density <- DEG_density %>%
#     mutate(
#       y = if_else(y < max(y)*(cluster_threshold/100), 0, y)
#     ) %>%
#     arrange(x) %>%
#     mutate(
#       border = case_when(
#         y > 0 & lag(y, default = dplyr::last(y)) == 0 ~ "start",      # potential bug here. with 1 base clusters!
#         y > 0 & lead(y, default = dplyr::first(y)) == 0 ~ "end",
#         TRUE ~ NA
#       )
#     ) %>%
#     filter(complete.cases(.))
#
#   DEG_clusters <- data.frame(start_position = numeric(), end_position = numeric())
#
#   nrow_degdensity <- nrow(DEG_density)
#   for (i in seq_len(nrow_degdensity / 2)) {
#
#     if(i == nrow_degdensity / 2){ #last iteration
#       aux <- data.frame(
#         start_position = min(DEG_density$x),
#         end_position = max(DEG_density$x)
#       )
#       DEG_clusters <- rbind(DEG_clusters, aux)
#
#     }else{
#       start_idx <- which(DEG_density$border == "start")[1]
#       end_idx <- which(DEG_density$border == "end" & seq_along(DEG_density$border) > start_idx)[1]
#
#       aux <- data.frame(
#         start_position = DEG_density$x[start_idx],
#         end_position = DEG_density$x[end_idx]
#       )
#
#       DEG_density <- DEG_density[-c(start_idx, end_idx), ]
#       DEG_clusters <- rbind(DEG_clusters, aux)
#     }
#   }
#
#   # Size
#   DEG_clusters <- DEG_clusters %>%
#     mutate(
#       start_position = round(16023*start_position/(2*pi)),
#       end_position = round(16023*end_position/(2*pi)),
#       size = case_when(
#         start_position < end_position ~ end_position - start_position,
#         end_position < start_position ~ 16023 - start_position + end_position
#       )
#     )
#
#   # Features
#   DEG_clusters$all_features <- NA
#   DEG_clusters$DEGs <- NA
#
#   for (i in 1:nrow(DEG_clusters)) {
#     DEG_clusters$all_features[i] <- DEdf %>%
#       filter(!!sym(avg_position) < DEG_clusters$end_position[i], !!sym(avg_position) > DEG_clusters$start_position[i]) %>%
#       pull(!!sym(gene_col)) %>%
#       paste(collapse = ";")
#
#     DEG_clusters$DEGs[i] <- DEdf %>%
#       filter(!!sym(DEG) %in% DEG_type) %>%
#       filter(!!sym(avg_position) < DEG_clusters$end_position[i], !!sym(avg_position) > DEG_clusters$start_position[i]) %>%
#       pull(!!sym(gene_col)) %>%
#       paste(collapse = ";")
#   }
#
#   # counting
#   DEG_clusters <- DEG_clusters %>%
#     mutate(
#       n_features = str_count(all_features, ";") + 1,
#       n_DEG = str_count(DEGs, ";") + 1
#     )
#
#   # Hypergeometric test
#   total_DEG <- DEdf %>% filter(!!sym(DEG) %in% DEG_type) %>% nrow()
#   total_no_DEG <- DEdf %>% filter(!!sym(DEG) == "NO") %>% nrow()
#
#   DEG_clusters$pval <- apply(DEG_clusters, 1, function(x){
#     phyper(
#       as.numeric(x[["n_DEG"]]) - 1,
#       total_DEG,
#       total_no_DEG,
#       as.numeric(x[["n_features"]]),
#       lower.tail = FALSE
#     )
#   })
#
#   # padj and score
#   DEG_clusters <- DEG_clusters %>%
#     mutate(
#       pval = p.adjust(DEG_clusters$pval, method = padj_method),
#       score = -log10(pval)
#     ) %>%
#     arrange(-score) %>%
#     filter(n_DEG > 1) # removing clusters with only 1 DEG
#
#   DEG_clusters <- DEG_clusters %>%
#     ungroup() %>%
#     mutate(
#       cluster_num = seq(1,nrow(DEG_clusters)),
#       start_position = as.numeric(start_position),
#       end_position = as.numeric(end_position),
#       size = as.numeric(size)
#     ) %>%
#     dplyr::select(cluster_num, start_position, end_position, size, all_features, DEGs, n_features, n_DEG, pval, score)
#
#   chromoObject@density$MT[[paste(DEG_type, collapse = ifelse(length(DEG_type) > 1, "_", ""))]] = list(
#     DEG_clusters = DEG_clusters,
#     bandwidth = bandwidth,
#     threshold = cluster_threshold,
#     padj_method = padj_method
#   )
#
#   return(chromoObject)
# }
#
#
# #######################################
# #####    chromo Density Plot MT  ######
# #######################################
#
# chromoDensityPlotMT <- function(
#     chromoObject,
#     DEG_type = list("UP", "DOWN"), # list("UP", "DOWN"), "UP", "DOWN"
#     n_top_clusters = 3,
#     include_genes = F,
#     include_density = T,
#
#     color_enrich = "#990099",
#     color_enrich_up = "#dd2200",
#     color_enrich_down = "#0022dd"
# ){
#
#   gene_col <- chromoObject@columns$gene_col
#   fc_col <- chromoObject@columns$fc_col
#   chromosome <- chromoObject@columns$chromosome
#   start_position <- chromoObject@columns$start_position
#   end_position <- chromoObject@columns$end_position
#   avg_position <- chromoObject@columns$avg_position
#   gene_length <- chromoObject@columns$gene_length
#   DEG <- chromoObject@columns$DEG
#   DEdf <- chromoObject@data %>%
#     filter(!!sym(chromosome) == "MT")
#
#   if(!is.list(DEG_type)){
#     DEG_type <- list(DEG_type)
#   }
#
#   density_type <- paste(DEG_type, collapse = ifelse(length(DEG_type) > 1, "_", ""))
#
#   plot_color <- ifelse(length(DEG_type) == 2, color_enrich, ifelse("UP" %in% DEG_type, color_enrich_up, color_enrich_down))
#
#   top_clusters <- chromoObject@density$MT[[density_type]][["DEG_clusters"]] %>%
#     arrange(-score) %>%
#     head(n_top_clusters)
#
#   r_min <- 1.6
#   r_max <- 2
#
#   # plot
#   mt_plot <- ggplot() +
#     labs(
#       title = NULL,
#       x = NULL,
#       y = NULL
#     ) +
#     theme_minimal() +
#     scale_y_continuous(expand = c(0, 0)) +
#     scale_x_continuous(expand = c(0, 0)) +
#     coord_cartesian(clip = "off") +
#     theme(
#       plot.background = element_rect(fill = "white", color = NA),
#       panel.background = element_rect(fill = "white", color = NA),
#       panel.border = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       axis.line = element_line(color = "black", linewidth = 0.5),
#       axis.title.y = element_text(face = "bold", size = 20, angle = 0, vjust = 0.5),
#       axis.text.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.title.x = element_blank(),
#       axis.line.x = element_blank(),
#       axis.line.y = element_blank(),
#       legend.position = "none"
#     ) +
#     annotate(
#       "text",
#       x = 0,
#       y = 0,
#       label = "MT",
#       color = "black",
#       size = 40
#     )
#
#   # Density
#   if(include_density){
#
#     # DEGs density
#     DEG_density <- DEdf %>%
#       filter(!!sym(DEG) %in% DEG_type) %>%
#       pull(!!sym(avg_position))
#
#     DEG_density <- (DEG_density / 16023) * 2 * pi
#     DEG_density <- circular(DEG_density, units = "radians", modulo = "asis")
#
#     dens <- density.circular(
#       DEG_density,
#       bw = chromoObject@density$MT[[density_type]]$bandwidth,
#       kernel = "wrappednormal"
#     )
#
#     dens_deg <- data.frame(x = dens$x, y_deg = dens$y) %>%
#       mutate(y_deg = y_deg/max(y_deg))
#
#     # # All features density
#     # all_density <- DEdf %>%
#     #   pull(!!sym(avg_position))
#     #
#     # all_density <- (all_density / 16023) * 2 * pi
#     # all_density <- circular(all_density, units = "radians", modulo = "asis")
#     #
#     # dens <- density.circular(
#     #   all_density,
#     #   bw = chromoObject@density$MT[[density_type]]$bandwidth,
#     #   kernel = "wrappednormal"
#     # )
#     #
#     # dens_all <- data.frame(x = dens$x, y_all = dens$y) %>%
#     #   mutate(y_all = y_all/max(y_all))
#
#     # # Difference
#     # dens <- dens_deg %>%
#     #   left_join(dens_all, by = "x") %>%
#     #   mutate(
#     #     y = y_deg - y_all,
#     #     y = case_when(
#     #       y < 0 ~ 0,
#     #       TRUE ~ y
#     #     ),
#     #     y = 1.5*y/max(y),
#     #     x_cart = r_max * cos(x) * (1 + y),
#     #     y_cart = r_max * sin(x) * (1 + y)
#     #   ) %>%
#     #   arrange(x)
#
#     dens <- dens_deg %>%
#       mutate(
#         y = y_deg,
#         y = y/max(y),
#         x_cart = r_max * sin(x) * (1 + y*0.4),
#         y_cart = r_max * cos(x) * (1 + y*0.4)
#       ) %>%
#       arrange(x)
#
#     mt_plot <- mt_plot +
#       geom_polygon(
#         data = dens,
#         aes(x = x_cart, y = y_cart),
#         color = plot_color,
#         linewidth = 0.5,
#         fill = paste0(plot_color, "77")
#       )
#   }
#
#   # Base arc
#   mt_plot <- mt_plot +
#     geom_arc_bar(
#       aes(x0 = 0,
#           y0 = 0,
#           r0 = r_min,
#           r  = r_max,
#           start = 0,
#           end = pi * 2),
#       fill  = "white",
#       color = "black",
#       size  = 1.2
#     )
#
#   # Clusters
#   mt_plot <- mt_plot +
#     geom_arc_bar(
#       data = top_clusters,
#       aes(
#         x0    = 0,
#         y0    = 0,
#         r0    = r_min,
#         r     = r_max,
#         start = 2 * pi * start_position / 16023,
#         end   = 2 * pi * end_position   / 16023,
#         group = cluster_num,
#         fill  = colorRampPalette(c("white", plot_color))(256)[
#           as.integer(rescale(score, to = c(0, 1), from = c(min(score), max(score)))) * 255 + 1
#         ]
#       ),
#       color = plot_color,
#       alpha = 0.8,
#       size  = 1.2
#     ) +
#     geom_text(
#       data = top_clusters,
#       aes(
#         x = (r_max + r_min) / 2 * sin(((start_position + end_position) / 2) / 16023 * 2 * pi),
#         y = (r_max + r_min) / 2 * cos(((start_position + end_position) / 2) / 16023 * 2 * pi),
#         label = cluster_num
#       ),
#       color = "black",
#       size  = 10,
#       fontface = "bold"
#     )
#
#   # Genes
#   if(include_genes){ ######################### ta cropando as bolinhas dos genes! FIX!!! coord_cartesian(clip = "off")?
#     aux <- DEdf %>%
#       mutate(
#         x = r_max * sin((avg_position/16023) * 2 * pi),
#         y = r_max * cos((avg_position/16023) * 2 * pi)
#       )
#
#     mt_plot <- mt_plot +
#       geom_point( # not DEGs
#         data = aux %>% filter(!!sym(DEG) == "NO"),
#         aes(x = x, y = y),
#         color = "#00000044",
#         size = 8
#       ) +
#       geom_point( # DEGs
#         data = aux %>% filter(!!sym(DEG) %in% DEG_type),
#         aes(x = x, y = y),
#         color = paste0(plot_color, "aa"),
#         size = 8
#       )
#   }
#
#   return(mt_plot)
# }


#######################################
######    chromo Interactions   #######
#######################################             for integrated HI-C data across several cell types, tissues and conditions



#####################################
######    chromo Regulation   #######
#####################################             for integrated chip-seq data across several cell types, tissues and conditions




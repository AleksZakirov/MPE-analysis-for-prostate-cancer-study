subset_cols <- function(df, pattern1, pattern2) {
  sample_names = colnames(df)[2:ncol(df)] # get sample names
  columns_pattern1 = sample_names[grep(pattern1, sample_names)]
  columns_pattern2 = sample_names[grep(pattern2, sample_names)]
  data_subset = df[, c('mz', columns_pattern1, columns_pattern2)] # subset the data
  data_subset = data_subset %>% column_to_rownames("mz") # convert mz to rownames
  # convert to numeric
  data_subset = data_subset %>% 
    mutate_if(is.character, as.numeric) %>% 
    mutate_if(is.factor, as.numeric)
  
  return(data_subset)
}

# function to calculate the fc, log2fc, pvalue, adjpvalue
get_stats = function(subset_df, pattern1, pattern2, pval_type){
  group1_df = subset_df %>% 
    select(dplyr::starts_with(pattern1))
  group2_df = subset_df %>%
    select(dplyr::starts_with(pattern2))
  
  
  # define empty data frame
  stats_df = data.frame(group1 = rowMeans(group1_df), group2 = rowMeans(group2_df))
  stats_df$fc = stats_df$group1/stats_df$group2
  stats_df$log2fc = log2(stats_df$fc)
  stats_df$pvalue = sapply(1:nrow(subset_df), function(x) {
    wilcox.test(as.numeric(group1_df[x,]), as.numeric(group2_df[x,]), paired = FALSE, exact = FALSE)$p.value
  })

  # rename columns
  colnames(stats_df)[1] = paste0('mean_',pattern1)
  colnames(stats_df)[2] = paste0('mean_',pattern2)
  
  # add adjusted pvalue
  
  stats_df$pvalue_adj = p.adjust(stats_df$pvalue, method = "BH")
  
  stats_df = stats_df %>% tibble::rownames_to_column('mz')
  
  # convert to numeric
  stats_df = stats_df %>% 
    mutate_if(is.character, as.numeric) %>% 
    mutate_if(is.factor, as.numeric)
  
  
  
  return(stats_df)
}

# function to get fc and pval results
get_group_fc_results = function(data, pattern1, pattern2, comparison_name, pval_type, save, cust_column_name, cust_column_value) {
  
  df_subset = subset_cols(data, pattern1 = pattern1, pattern2 = pattern2)
  print(head(df_subset))
  
  # calculate the fc, log2fc, pvalue/adjpvalue
  fc_results = get_stats(df_subset, pattern1 = pattern1, pattern2 = pattern2, pval_type = pval_type)

  # check if there are any NA values
  if (any(is.na(fc_results$pvalue))){
    print(paste("NA values found in pvalue column for", comparison_name))
  }
  if (any(is.na(fc_results$fc))){
    print(paste("NA values found in fc column for", comparison_name))
  }

  # drop NA values
  print(paste("Number of rows before dropping NA values:", nrow(fc_results)))

  fc_results = fc_results %>% drop_na()
  print(paste("Number of rows after dropping NA values:", nrow(fc_results)))
  
  # add custom column if needed
  if (!is.null(cust_column_name)){
    fc_results[cust_column_name] = cust_column_value
  }
  
  # write to csv
  if (save){
    save_name = paste0('Data/', comparison_name, '_fc_', pval_type, '_', 'results.csv')
    write.csv(fc_results, save_name, row.names = FALSE)
    print(paste("Saved", save_name))
  }
  return(fc_results)
}

# function to calculate molecule+adduct weights for supplied adducts and a list of mw
compute_adduct_weights = function (adduct_table, mw) {
  # adduct table should contain only adducts for which the mw should be calculated
  ion.name <- adduct_table$Ion_Name
  ion.mass <- adduct_table$Ion_Mass
  mass.list <- as.list(ion.mass)
  mass.user <- lapply(mass.list, function(x) eval(parse(text = paste(gsub("PROTON", 1.00727646677, x)))))
  mass.df = as.data.frame(mass.user)
  colnames(mass.df) = ion.name
  mass.df = cbind(mw, mass.df)
  
  return(mass.df)
}

# function to extract compounds and corresponding associated adducts

get_matched_comps = function(adducts_table, mz_list, ms.type) {
  
  # get unique entries from KEGG peaks
  adducts_dedup = distinct(adducts_table)
  
  #get peaks
  peaks_vector = adducts_dedup[,1]
  
  # drop mz column and get adducts
  adducts_dedup = adducts_dedup[,-1]
  adducts_vector = colnames(adducts_dedup)
  
  # calculate the tolerance for all peaks
  tolerance = mz_list*ms.type*1e-06
  
  # check which experimental peaks are within the tolerance; get the coordinates of the compounds that are within the tolerance
  matched_coords = lapply(1:length(mz_list),
                          function(x) {
                            suppressWarnings(cbind(which((abs(adducts_dedup-mz_list[x])<= tolerance[x]) == TRUE, arr.ind = TRUE), exp_peak = mz_list[x]))
                          })
  
  matched_coords = data.frame(do.call(rbind, matched_coords))
  
  # convert to a dataframe
  matched_comps_df = data.frame(mz = peaks_vector[matched_coords$row], adduct = adducts_vector[matched_coords$col], exp_peak = matched_coords$exp_peak)
  
  # drop duplicated entries
  matched_comps_df = distinct(matched_comps_df)
  
  return(matched_comps_df)
  
}

# function to merge the fc table with annotation
merge_to_annot = function(fc_results_table, annotation_table, keep_sig = T, keep_annot = F, save = F){
  # left join tables
  fc_res_table = read.csv(fc_results_table)
  fc_annot_table = left_join(fc_res_table, annotation_table, by = c('mz' = 'mzs'))
  
  save_name_pattern = 'results_annot'
  if (keep_sig != FALSE){
    
    if (keep_sig == TRUE){
      keep_sig = 0.05
    }
    
    fc_annot_table = fc_annot_table %>% filter(pvalue <= keep_sig) %>% arrange(pvalue)
    save_name_pattern = paste0(save_name_pattern, '_sig')
  }
  
  if (keep_annot) {
    fc_annot_table = fc_annot_table %>% filter(!is.na(compound))
    save_name_pattern = paste0(save_name_pattern, '_mapped_only')
  }
  
  if (save) {
    save_name = gsub("results", save_name_pattern, fc_results_table)
    write.csv(fc_annot_table, save_name, row.names = F)
  }
  return(fc_annot_table)
}

compute_ORA = function(cpd_vector = NA, path_set = NA, hit_filt = NA) 
{
  q.size <- length(cpd_vector)
  if (is.na(cpd_vector) || q.size == 0) {
    print("No valid compound names found!")
    return(0)
  }
  
  uniq.count <- length(unique(unlist(path_set, use.names = FALSE)))
  set.size <- length(path_set)
  
  if (set.size == 1) {
    print("Cannot perform enrichment analysis on a single metabolite set!")
    return(0)
  }
  hits <- lapply(path_set, function(x) {
    x[x %in% cpd_vector]
  })
  hit.num <- unlist(lapply(hits, function(x) length(x)), use.names = FALSE)
  if (sum(hit.num > 0) == 0) {
    print("No match was found to the selected metabolite set library!")
    return(0)
  }
  set.num <- unlist(lapply(path_set, length), use.names = FALSE)
  res.mat <- matrix(NA, nrow = set.size, ncol = 6)
  rownames(res.mat) <- names(path_set)
  colnames(res.mat) <- c("total", "expected", "hits", "Raw p", 
                         "Holm p", "FDR")
  for (i in 1:set.size) {
    res.mat[i, 1] <- set.num[i]
    res.mat[i, 2] <- q.size * (set.num[i]/uniq.count)
    res.mat[i, 3] <- hit.num[i]
    res.mat[i, 4] <- phyper(hit.num[i] - 1, set.num[i], 
                            uniq.count - set.num[i], q.size, lower.tail = F)
  }
  res.mat[, 5] <- p.adjust(res.mat[, 4], "holm")
  res.mat[, 6] <- p.adjust(res.mat[, 4], "fdr")
  res.mat <- res.mat[hit.num > 0, ]
  ord.inx <- order(res.mat[, 4])
  res.mat <- signif(res.mat[ord.inx, ], 3)
  res.mat = data.frame(res.mat)
  
  
  if (!is.na(hit_filt)){
    res.mat = res.mat[res.mat['hits'] >= hit_filt,]
    print(paste0('Keeping pathways with >= ', hit_filt, ' hits'))
  }
  
  return(res.mat)
}




calc_perm = function(input_cpdlist, total_matched_cpds, pathway_set, matches.res, input_mzlist){
  ora.vec <- input_cpdlist
  query_set_size <- length(ora.vec)
  current.mset <- pathway_set
  total_cpds <- unique(total_matched_cpds)
  total_feature_num <- length(total_cpds)
  size <- negneg <- vector(mode = "list", length = length(current.mset))
  cpds <- lapply(current.mset, function(x) intersect(x, total_cpds))
  feats <- lapply(current.mset, function(x) intersect(x, ora.vec))
  feat_len <- unlist(lapply(feats, length))
  set.num <- unlist(lapply(cpds, length))
  negneg <- sizes <- vector(mode = "list", length = length(current.mset))
  for (i in seq_along(current.mset)) {
    cpd_mzs = unique(matches.res[matches.res['compound_id'] == feats[i], 'mzs'])
    sign_mz_count = length(intersect(cpd_mzs, sign_mz))
    sizes[[i]] <- min(feat_len[i], sign_mz_count)
    negneg[[i]] <- total_feature_num + sizes[[i]] - set.num[i] - 
      query_set_size
  }
  unsize <- as.integer(unlist(sizes))
  res.mat <- matrix(0, nrow = length(current.mset), ncol = 1)
  fishermatrix <- cbind(unsize - 1, set.num, (query_set_size + 
                                                unlist(negneg)), query_set_size)
  res.mat[, 1] <- apply(fishermatrix, 1, function(x) phyper(x[1], 
                                                            x[2], x[3], x[4], lower.tail = FALSE))
  perm_records <- list(res.mat, as.matrix(unsize))
  return(perm_records)
}

customORA = function(kegg_path = NULL, path_size = NULL, matched_peaks = NULL, pval_threshold = 0.05, hit_threshold = 2){
  
  if (is.null(kegg_path)) {
    print('Provide kegg_path')
    stop()
  }
  
  if (is.null(path_size)) {
    print('Provide path_size')
    stop()
  }
  print(
    paste0(
      'Using path_size = ',
      path_size,
      ', i.e keeping only pathways with at least ',
      path_size,
      ' compounds.'
    )
  )
  
  if (is.null(hit_threshold)) {
    print('Provide hit_threshold')
    stop()
  }
  
  # read in updated KEGG database
  kegg_db = read.csv(kegg_path)
  
  # keep only data for pathways with at least path_size compounds under it
  filt_kegg_db = kegg_db %>%
    filter(pathway_id %in% filter(plyr::count(kegg_db$pathway_id), freq >=
                                    path_size)$x)
  
  print(paste(
    'Dropped',
    n_distinct(kegg_db$pathway_id) - n_distinct(filt_kegg_db$pathway_id),
    'pathways with only one compound'
  ))
  
  # drop rows without given mass
  nbefore = n_distinct(filt_kegg_db$compound_id)
  filt_kegg_db = filt_kegg_db %>% drop_na(complete_compound_mass)
  nafter = n_distinct(filt_kegg_db$compound_id)
  
  print(paste('Dropped', nbefore - nafter, 'compounds without mass.'))
  
  
  # prepare a pathway:compounds dictionary from dataframe for ORA
  path_set_dict = vector(mode = 'list', length = 0)
  
  path_list = filt_kegg_db %>% pull(pathway_name) %>% unique()
  
  path_set_dict = lapply(path_list, function (x) {
    filt_kegg_db[filt_kegg_db$pathway_name == x, 'compound_id']
  })
  names(path_set_dict) = path_list
  
  # remove pathways with 1000+ entries - these are not useful
  path_set_dict_filt = path_set_dict[which(lapply(path_set_dict, length) < 1000)]
  
  # filter out entries without pathways mapped
  matched_peaks = matched_peaks %>% filter(pathway_name != '')
  
  # get a list of significant compounds for ORA and filter
  p005 = matched_peaks %>% filter(pvalue <= pval_threshold)
  
  p005up = p005 %>% filter(mean_fc > 1)
  p005down = p005 %>% filter(mean_fc <= 1)
  
  sign_cpd_up = p005up %>% pull(compound_id) %>% unique()
  sign_cpd_down = p005down %>% pull(compound_id) %>% unique()
  
  print(paste0('Significantly UPregulated compounds: ', length(sign_cpd_up)))
  print(paste0('Significantly DOWNregulated compounds: ', length(sign_cpd_down)))
  
  # run ORA
  print('Running ORA...')
  ora_results_up = compute_ORA(cpd_vector = sign_cpd_up, path_set = path_set_dict_filt, hit_filt = hit_threshold)
  ora_results_down = compute_ORA(cpd_vector = sign_cpd_down, path_set = path_set_dict_filt, hit_filt = hit_threshold)
  
  ora_results_up = ora_results_up %>% mutate(direction = 'up') %>% rownames_to_column('pathway')
  ora_results_down = ora_results_down %>% mutate(direction = 'down') %>% rownames_to_column('pathway')
  ora_results = rbind(ora_results_up ,ora_results_down) %>% 
    dplyr::mutate(FC = hits/expected, .after = FDR) %>% 
    dplyr::mutate(log2FC = log2(hits/expected), .after = FC) %>% 
    mutate(FC = ifelse(direction == 'down', FC*-1, FC)) %>% 
    mutate(log2FC = ifelse(direction == 'down', log2FC*-1, log2FC))
  print('Running ORA...done')
  
  return(ora_results)
}


customPeakMap = function(kegg_path = NULL, results_path = NULL, path_size = 2, ppm = 5, mode = NULL){
  
  if (is.null(kegg_path)) {
    print('Provide kegg_path')
    stop()
  }
  
  if (is.null(results_path)) {
    print('Provide results_path')
    stop()
  }
  
  if (is.null(path_size)) {
    print('Provide path_size')
    stop()
  }
  print(
    paste0(
      'Using path_size = ',
      path_size,
      ', i.e keeping only pathways with at least ',
      path_size,
      ' compounds.'
    )
  )
  
  if (is.null(ppm)) {
    print('Provide ppm')
    stop()
  }
  print(paste0('Using ppm = ', ppm))
  
  if (is.null(mode)) {
    print('Provide mode')
    stop()
  }
  if (!mode %in% c('negative', 'positive')) {
    print('Mode has to be either "positive" or "negative".')
    stop()
  }
  print(paste0('Using mode = ', mode))
  
  
  # read in updated KEGG database
  kegg_db = read.csv(kegg_path)
  
  # read the fold-change pval results table
  fc_table = read.csv(results_path)
  
  exp_peak_list = fc_table$mz
  
  # select adducts of interest
  if (mode == 'negative') {
    adduct_formulas = read.csv('Resource_tables/neg_adduct_table.csv')
    my_adduct_formulas = adduct_formulas[c(1, 4, 8), ]
    
  } else {
    adduct_formulas = read.csv('Resource_tables/pos_adduct_table.csv')
    my_adduct_formulas = adduct_formulas[c(2, 5, 13), ]
  }
  
  
  # keep only data for pathways with at least path_size compounds under it
  filt_kegg_db = kegg_db %>% filter(pathway_id %in% filter(plyr::count(kegg_db$pathway_id), freq >= path_size)$x)
  
  print(paste(
    'Dropped',
    n_distinct(kegg_db$pathway_id) - n_distinct(filt_kegg_db$pathway_id),
    'pathways with only one compound'
  ))
  
  # drop rows without given mass
  nbefore = n_distinct(filt_kegg_db$compound_id)
  filt_kegg_db = filt_kegg_db %>% drop_na(complete_compound_mass)
  nafter = n_distinct(filt_kegg_db$compound_id)
  
  print(paste('Dropped', nbefore - nafter, 'compounds without mass.'))
  
  print('Computing adduct masses...')
  # compute theoretical masses for all compounds+adducts in the kegg dataset
  adducts_computed = compute_adduct_weights(my_adduct_formulas,filt_kegg_db$complete_compound_mass)
  
  print('Computing adduct masses...done')
  
  print('Matching peaks to computed adducts...')
  
  # match experimental peaks to calculated adducts
  matched_comps = get_matched_comps(adducts_table = adducts_computed,
                                    mz_list = exp_peak_list,
                                    ms.type = ppm)
  print('Matching peaks to computed adducts...done')
  
  # merge with kegg data for complete annotation (join on theoretical mass)
  matched_comps_annot = left_join(matched_comps, filt_kegg_db, by = c('mz' = 'complete_compound_mass'))
  
  # convert mz to numeric
  matched_comps_annot$mz = as.numeric(matched_comps_annot$mz)
  
  # add pval and fc data (join on experimental peak)
  matched_comps_annot = matched_comps_annot %>% left_join(fc_table, by = c('exp_peak' =
                                                                             'mz'))
  
  # add mean fc - averaging FC across the same compounds with different adducts
  mean_fc_df = matched_comps_annot %>% select(compound_id, compound_name, fc) %>% group_by(compound_id, compound_name) %>% dplyr::summarise(mean_fc = mean(fc)) %>% ungroup() %>% filter(compound_id != '')
  
  matched_comps_annot = matched_comps_annot %>% left_join(mean_fc_df)
  
  # print some basic stats
  print(paste("Unique peaks mapped:", n_distinct(matched_comps_annot$exp_peak)))
  print(paste(
    "Unique compounds mapped:",
    n_distinct(matched_comps_annot$compound_id)
  ))
  print(paste(
    "Unique compounds in LIPIDMAPS:",
    n_distinct(
      matched_comps_annot %>% filter(is_lm == 'True') %>% pull(compound_id)
    )
  ))
  
  # NOTE: some compounds have exactly the same molecular formula, hence they match the same peak.
  # This can inflate some pathway representation.
  
  # reorder columns
  matched_comps_annot = matched_comps_annot[, c(
    'exp_peak',
    'adduct',
    'pvalue',
    'pvalue_adj',
    colnames(fc_table)[2:3],
    'fc',
    'mean_fc',
    'mode',
    'mz',
    'compound_id',
    'compound_name',
    'compound_synonyms',
    'compound_formula',
    'is_hmdb',
    'is_lm',
    'lm_id',
    'lm_name',
    'lm_abbrev',
    'lm_formula',
    'pathway_id',
    'pathway_name',
    'human_pathway',
    'rat_pathway',
    'mouse_pathway'
  )]
  
  return(matched_comps_annot)
}


plotPathVolcano = function(df, min_hits = 2, pval = 0.05, fold_change = 2, grp1 = 'Group 1', grp2 = 'Group 2', x_max = NULL, y_max = NULL, title = 'Plot title', subtitle = 'Plot subtitle', caption = 'Plot caption', colour_sec = c('firebrick3','dodgerblue3','gray50')){
  require(tidyverse)
  require(ggplot2)
  # process dataframe
  label_grp1 = paste0("Up-regulated in ", grp1)
  label_grp2 = paste0("Up-regulated in ", grp2)
  
  df = df %>% mutate(colour = case_when(FC > fold_change & Raw.p <= pval ~ label_grp1,
                                        FC < -fold_change & Raw.p <= pval ~ label_grp2,
                                        TRUE ~ "Unchanged")) %>% 
    mutate(colour = factor(colour, levels = c(label_grp1, label_grp2, 'Unchanged'))) %>% 
    mutate(label = ifelse(Raw.p <= pval & abs(FC) >= fold_change, Pathway, '')) %>%
    mutate(fdr_sig = ifelse(FDR <= pval & abs(FC) >= fold_change, 'sig', 'non-sig')) %>% 
    filter(hits >= min_hits)
  
  p = ggplot(df, aes(FC, -log10(Raw.p))) + # -log10 conversion
    
    geom_point(aes(color = colour, size = hits/total)) +
    geom_point(data = df[df$fdr_sig == 'sig',], aes(size=hits/total), shape = 21, stroke=1.5) +
    
    xlab(expression("Pathway enrichment (FC)")) +
    ylab(expression("-log"[10]*"p.val"))+
    
    scale_color_manual(values = colour_sec) +
    scale_shape_manual(values = c(5,0))+
    
    guides(colour = guide_legend(override.aes = list(size=3)))+
    
    geom_vline(xintercept = c(-fold_change, fold_change), linetype = 'dashed', alpha = 0.4)+
    geom_hline(yintercept = -log10(pval), linetype = 'dashed', alpha = 0.4)+
    
    geom_text_repel(data = df[df$direction == 'up',],
                    aes(label = label), 
                    xlim = c(fold_change,NA),
                    size = 3, 
                    box.padding = 0.75, 
                    point.padding =0.5, 
                    min.segment.length = 0,
                    max.overlaps = Inf, 
                    max.time = 2)+
    geom_text_repel(data = df[df$direction == 'down',],
                    aes(label = label),
                    xlim = c(NA, -fold_change),
                    size = 3, 
                    box.padding = 0.75, 
                    point.padding =0.5, 
                    min.segment.length = 0,
                    max.overlaps = Inf, 
                    max.time = 2)+
    
    theme_minimal()+
    labs(title = title,
         subtitle = subtitle,
         caption = caption,
         color = NULL, size = 'Pathway coverage')+
    coord_cartesian(clip = "off")
  
  
  if (is.null(x_max)) {
    p = p + xlim(-max(abs(df$FC)), max(abs(df$FC)))
  } else {
    p = p + xlim(-x_max, x_max)
  }
  
  if (!is.null(y_max)) {
    p = p + ylim(0, y_max)
  }
  
  
  return(p)
  
}

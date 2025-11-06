server <- function(input, output, session){
  results <- reactiveVal(NULL); displayed_data <- reactiveVal(data.frame())
  seqs <- reactiveVal(NULL); using_sample <- reactiveVal(FALSE)
  
  observeEvent(input$fasta, { seqs(readAAStringSet(input$fasta$datapath)); using_sample(FALSE) })
  
  observeEvent(input$use_sample, {
    validate(need(length(sample_fasta) > 0, "Sample data file 'mhc_db/sample_data.fasta' could not be loaded."))
    seqs(sample_fasta)
    using_sample(TRUE)
    updateTextInput(session, "uniprot_id", value = "P18272")
  })
  
  observeEvent(input$submit, {
    disable("submit"); on.exit(enable("submit"), add = TRUE)
    validate(need(seqs(), "Upload or use sample FASTA"))
    validate(need(length(input$alleles) > 0, "Select ≥1 HLA allele"))
    
    cur <- isolate(seqs()); uid <- toupper(isolate(input$uniprot_id %||% ""))
    
    message("\n[Step] Computing multiple sequence alignment...")
    withProgress(message = "Computing MSA & consensus...", value = 0.2, {
      if(length(cur) > 1){
        msa_res <- suppressMessages(msa(cur, method = "ClustalW"))
        cons_vec <- get_consensus_vec(msa_res)
        cons_aln_len <- length(cons_vec)
        cons_ungap <- paste(cons_vec[cons_vec != "-"], collapse = "")
        map_ug2aln <- which(cons_vec != "-")
        cons_aln <- get_conserved_regions_aln(msa_res, input$min_conserved_len)
        cons_aln$aln_end <- cons_aln$aln_start + cons_aln$aln_length - 1L
        message("[MSA] Conserved blocks (alignment coords): ",
                paste0("[", cons_aln$aln_start, "–", cons_aln$aln_end, "]", collapse = " ; "))
      } else {
        s <- clean_aa_string(as.character(cur[[1]]))
        cons_vec <- strsplit(s, "")[[1]]
        cons_aln_len <- length(cons_vec)
        cons_ungap <- gsub("-", "", s)
        map_ug2aln <- seq_along(cons_vec[cons_vec != "-"])
        cons_aln <- data.frame(aln_start = 1L, aln_length = cons_aln_len)
      }
      incProgress(0.6, detail = "Consensus & conserved regions ready")
    })
    
    message("[Step] Searching local database for published epitopes...")
    # Note: 'cons_ungap' is passed to the function's 'consensus_ungapped' argument
    pub_df <- tryCatch({
      if(nzchar(uid)) find_published_epitopes_local(cons_ungap, uid) else data.frame()
    }, error = function(e){ message("  Error: ", e$message); data.frame() })
    
    pep_lengths <- seq(input$peptide_length[1], input$peptide_length[2])
    all_hits <- list()
    if(!isTRUE(input$debug_skip_api)){
      message("[Step] Querying IEDB for binding predictions...")
      withProgress(message = "Querying IEDB...", value = 0, {
        total <- length(input$alleles) * length(pep_lengths); step <- 0
        for(a in input$alleles){
          for(l in pep_lengths){
            step <- step + 1; incProgress(step/total, detail = paste(a, l))
            message("  → ", a, " (", l, "mer)")
            peps <- generate_peptides(cons_ungap, l); if(!length(peps)) next
            df <- query_iedb(peps, a, l)
            if(is_nonempty_df(df)){
              df$Start <- vapply(df$Peptide, function(p){
                
                # --- THIS WAS THE ERROR ---
                # It should be 'cons_ungap', not 'consensus_ungapped'
                pos <- as.integer(regexpr(p, cons_ungap, fixed = TRUE))
                # --- END FIX ---
                
                if(!is.na(pos) && pos > 0L) pos else NA_integer_
              }, integer(1))
              df <- df[!is.na(df$Start),]
              if(!nrow(df)) next
              df$End <- df$Start + nchar(df$Peptide) - 1L
              df$Length <- l
              df$Allele <- a
              all_hits[[paste(a, l)]] <- df
            }
          }
        }
      })
    }
    pred_df <- if(length(all_hits)) bind_rows(all_hits) else data.frame()
    if(nrow(pred_df)) pred_df$RowID <- seq_len(nrow(pred_df))
    binders <- if(nrow(pred_df)) filter(pred_df, Rank < 2) else pred_df
    
    results(list(predicted = binders, published = pub_df,
                 conserved_aln = cons_aln, consensus_vec = cons_vec,
                 consensus_ungapped = cons_ungap, map_ungap_to_aln = map_ug2aln,
                 consensus_aln_len = cons_aln_len))
    
    message("[✓] Pipeline complete. Results ready.")
  })
  
  output$top_table <- renderDT({
    r <- results(); df <- r$predicted
    if(!is_nonempty_df(df)) return(datatable(data.frame()))
    df <- arrange(df, Rank); displayed_data(df)
    datatable(df, selection = list(mode = "multiple", target = "row"),
              options = list(pageLength = 10), rownames = FALSE) %>%
      formatStyle("AffinityClass",
                  backgroundColor = styleEqual(c("Strong","Intermediate","Weak"),
                                               c("#e6ffe6","#ffffe6","#ffe6e6")))
  })
  
  output$published_table <- renderDT({
    r <- results(); df <- r$published
    if(!is_nonempty_df(df)) return(datatable(data.frame()))
    
    if(!inherits(df, "data.table")) setDT(df)
    
    df <- df[, .(
      `Main Epitope ID` = Main_Epitope_ID,
      Peptide,
      Allele,
      `# Assays` = `#_Assays`,
      Source = case_when(
        grepl("tcell", SourceFiles, ignore.case = TRUE) & 
          grepl("mhc", SourceFiles, ignore.case = TRUE) ~ "MHC / T Cell",
        grepl("mhc", SourceFiles, ignore.case = TRUE) ~ "MHC",
        grepl("tcell", SourceFiles, ignore.case = TRUE) ~ "T Cell",
        TRUE ~ SourceFiles
      ),
      Start, End, Assays
    )]
    
    datatable(
      df,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE,
      escape = FALSE
    ) %>%
      formatStyle(
        columns = names(df),
        whiteSpace = "nowrap"
      )
  })
  
  output$sequence_plot <- renderPlotly({
    r <- results()
    if (is.null(r))
      return(plot_ly() %>% layout(xaxis = list(title = "Position")))
    
    pred <- r$predicted
    pub  <- r$published
    cons_aln <- r$conserved_aln
    map  <- r$map_ungap_to_aln
    aln_len <- length(r$consensus_vec)
    sel <- input$top_table_rows_selected
    df_shown <- displayed_data()
    
    if (is_nonempty_df(pub)) {
      pub$Source <- case_when(
        grepl("tcell", pub$SourceFiles, ignore.case = TRUE) & 
          grepl("mhc", SourceFiles, ignore.case = TRUE) ~ "MHC / T Cell",
        grepl("mhc", pub$SourceFiles, ignore.case = TRUE) ~ "MHC",
        grepl("tcell", pub$SourceFiles, ignore.case = TRUE) ~ "T Cell",
        TRUE ~ pub$SourceFiles
      )
    }
    
    add_aln <- function(df) {
      if (!is_nonempty_df(df) || !"Start" %in% names(df)) return(df)
      df <- df[df$Start > 0 & df$Start <= length(map), , drop = FALSE]
      if (!nrow(df)) return(df)
      df$AlnStart <- map[df$Start]
      df
    }
    
    predA <- add_aln(pred)
    pubA  <- add_aln(pub)
    
    plt <- plot_ly()
    
    if (is_nonempty_df(predA)) {
      plt <- plt %>% add_trace(
        data = predA, x = ~AlnStart, y = ~Allele,
        type = "scatter", mode = "markers", key = ~RowID,
        marker = list(size = 10, opacity = 0.85, color = ~Rank,
                      colorscale = "Reds", reversescale = TRUE,
                      colorbar = list(title = list(text = "<b>Binding Rank</b>"))),
        text = ~paste0("<b>Predicted:</b> ", Peptide,
                      "<br>Consensus start: ", Start, "–", End,
                      "<br>Allele: ", Allele,
                      "<br>Rank: ", round(Rank, 2), "%",
                      "<br>Affinity: ", round(Affinity, 1), " nM"),
        hoverinfo = "text", name = "Predicted"
      )
    }
    
    if (!is.null(sel) && length(sel) > 0 && is_nonempty_df(df_shown)) {
      selected <- df_shown[sel, , drop = FALSE]
      selected <- add_aln(selected)
      if (nrow(selected) > 0) {
        plt <- plt %>% add_trace(
          data = selected, x = ~AlnStart, y = ~Allele,
          type = "scatter", mode = "markers", key = ~RowID,
          marker = list(size = 14, color = "#1f77b4", symbol = "circle",
                        line = list(color = "#1f77b4", width = 1.5)),
          hoverinfo = "text",
          text = ~paste0("<b>Selected:</b> ", Peptide,
                        "<br>Allele: ", Allele,
                        "<br>Start: ", Start, "–", End),
          name = "Selected"
        )
      }
    }
    
    if (is_nonempty_df(pubA)) {
      plt <- plt %>% add_trace(
        data = pubA, x = ~AlnStart, y = ~Allele,
        type = "scatter", mode = "markers",
        marker = list(color = "#00B894", size = 9, symbol = "circle"),
        text = ~paste0("<b>Published:</b> ", Peptide,
                      "<br>Consensus start: ", Start, "–", End,
                      "<br>Main ID: ", Main_Epitope_ID %||% "",
                      "<br>Allele: ", Allele,
                      "<br># Assays: ", `#_Assays` %||% "",
                      "<br>Source: ", Source %||% ""),
        hoverinfo = "text", name = "Published"
      )
    }
    
    shape_list <- list()
    if (isTRUE(input$highlight_regions) && is_nonempty_df(cons_aln)) {
      cons_aln$aln_start <- as.numeric(cons_aln$aln_start)
      cons_aln$aln_length <- as.numeric(cons_aln$aln_length)
      cons_aln$aln_end <- cons_aln$aln_start + cons_aln$aln_length - 1
      shape_list <- lapply(seq_len(nrow(cons_aln)), function(i) {
        reg <- cons_aln[i, ]
        list(type = "rect", xref = "x", yref = "paper",
             x0 = reg$aln_start - 0.5, x1 = reg$aln_end + 0.5,
             y0 = 0, y1 = 1, fillcolor = "rgba(30,144,255,0.25)",
             line = list(color = "transparent"), layer = "below")
      })
    }
    
    plt %>%
      layout(
        xaxis = list(title = "Position",
                     range = c(0.5, aln_len + 0.5),
                     tickmode = "linear", dtick = 20, showgrid = FALSE),
        yaxis = list(title = "HLA Allele", type = "category",
                     gridcolor = "rgba(0,0,0,0.1)", gridwidth = 0.5),
        legend = list(orientation = "h", y = -0.2),
        plot_bgcolor = "rgba(240,240,240,0.5)",
        shapes = shape_list
      )
  }) %>%
    bindCache(results(), input$top_table_rows_selected, input$highlight_regions)
}

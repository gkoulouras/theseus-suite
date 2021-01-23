function(input, output, session) { 
  
  # disable tab2 on page load
  shinyjs::js$disableTab("digestion_results_tab")
  shinyjs::disable("no_of_misceavages")
  shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
  
  observeEvent(input$switchButtonDigest, {
    if (input$switchButtonDigest == TRUE) {
      shinyjs::show("fastaFile", anim=TRUE, animType="slide")
      shinyjs::show("fastaformat", anim=TRUE, animType="slide")
      shinyjs::hide("proteinListTextInput", anim=TRUE, animType="fade")
      shinyjs::hide("uploadExampleDT", anim=FALSE, animType="slide")
      shinyjs::hide("submitProtList", anim=TRUE, animType="fade")
    }
    else {
      shinyjs::hide("fastaFile", anim=TRUE, animType="slide")
      shinyjs::hide("fastaformat", anim=TRUE, animType="slide")
      shinyjs::show("proteinListTextInput", anim=TRUE, animType="fade")
      shinyjs::show("uploadExampleDT", anim=FALSE, animType="slide")
      shinyjs::show("submitProtList", anim=TRUE, animType="fade")
    }
  })
  
  observeEvent(input$uploadExampleDT, {
    updateTextInput(session, "proteinListTextInput", value="")
    updateTextInput(session, "proteinListTextInput", value=example_dataset)
  })

  
  observeEvent(input$tabs, {
    if (input$tabs == "resetsmenuitem") {
      js$reset()
      shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
    }
  })
  
  rv <- reactiveValues() 
  observeEvent(input$Id049, {
    rv$choice=input$Id049
    #print(rv$choice)
  })
  
  protease <- reactiveValues()
  observeEvent(input$protease, {
    protease$choice=input$protease
  })
  
  miscleavage <- reactiveValues()
  observeEvent(input$no_of_misceavages, {
    miscleavage$choice=input$no_of_misceavages
  })
  
  max_peptide_mass <- reactiveValues()
  observeEvent(input$mol_weight, {
    max_peptide_mass$choice=input$mol_weight
  })
  
  min_peptide_length <- reactiveValues()
  observeEvent(input$peptide_length, {
    min_peptide_length$choice <- input$peptide_length
  })
  
  observeEvent(input$fastaFile, {
    shinyjs::js$enableTab("digestion_results_tab")
    updateTabsetPanel(session, "digestion_tabSet", selected = "digestion_results_tab")
    updateVerticalTabsetPanel(session, "vertical_results_tabSet", selected = "vertical_entries_tab")
    shinyjs::addClass(selector = "body", class = "sidebar-collapse")
  })
  
  output$selected_enzyme <- renderText({ 
    HTML(paste0("Proteolutic enzyme: ", tools::toTitleCase(protease$choice) ) )
  })
  
  output$minimum_pep_len <- renderText({ 
    HTML(paste0("Min. peptide length: ", min_peptide_length$choice) )
  })
  
  output$maximum_pep_mass <- renderText({ 
    HTML(paste0("Max. peptide mass: ", max_peptide_mass$choice) )
  })
  
  output$missed_cleavages <- renderText({ 
    paste0("Allowed miscleavages: ", miscleavage$choice)
  })
  
  protGroups_digest <- reactive({
    req(input$protGroupsDigestion)
    ext <- tools::file_ext(input$protGroupsDigestion$datapath)
    validate(need(ext == "txt", "Please upload a txt file"))
    dt.seq <- read.table(input$protGroupsDigestion$datapath, header = T, sep = "\t", fill = TRUE, quote = "", check.names = F, stringsAsFactors = F)
    required_columns <- c("id","Protein IDs", "Sequence coverage [%]")
    if(all(required_columns %in% names(dt.seq)) == F) {
      tbl <- NULL
      shinyjs::runjs("$('#proteinGroups_progress').css('visibility','hidden')")
      js$swalErrorAlert()
    }
    else{
      columns2keep <- names(dt.seq)
      protGroups_keep <- dt.seq
      x <- protGroups_keep$`Protein IDs`
      protGroups_keep["count_of_prot_ID"] <- data.frame(matrix(unlist(lengths(regmatches(x, gregexpr(";", x))) + 1), nrow=length(lengths(regmatches(x, gregexpr(";", x)))), byrow=T))
      prot_IDs <- unlist(strsplit(as.character(x), split=";"))
      protGroups_keep <- protGroups_keep[rep(row.names(protGroups_keep), protGroups_keep$count_of_prot_ID), 1:ncol(protGroups_keep) - 1]
      protGroups_keep["UniProtID"] <- data.frame(prot_IDs)
      
      identifiable_proteins_subset <- identifiable_proteins()[,c("UniProtID", "Theoritical max. coverage [%]")]
      names(identifiable_proteins_subset) <- c("UniProtID", "Theoritical max coverage per protein ID [%]")
      # nonidentifiable_proteins_selected_cols <- non_identifiable_proteins()[,c("UniProtID", "Theoritical max. coverage [%]")]
      
      merged_dt <- merge(protGroups_keep, identifiable_proteins_subset, by = "UniProtID", all.x = TRUE)
      
      tbl <- as.data.frame(merged_dt[,c(columns2keep, "UniProtID", "Theoritical max coverage per protein ID [%]")] %>% group_by(id) %>%
        mutate(`Protein IDs` = paste(`UniProtID`, collapse=";"), `Theoritical max coverage per protein ID [%]` = paste(`Theoritical max coverage per protein ID [%]`, collapse=";")) %>%
        select(all_of(columns2keep), `Theoritical max coverage per protein ID [%]`) %>% distinct())
      
    
      # dt_out <- protGroups_digest()
      
      temp_dt <- tbl[,c("id", "Protein IDs", "Sequence coverage [%]", "Theoritical max coverage per protein ID [%]")]
      # temp_dt$`Theoritical max. coverage [%]` <- sapply(strsplit(as.character(temp_dt$`Theoritical max coverage per protein ID [%]`),";"), function(x, na.strings = NA) if(all(is.na(as.numeric(x)))) NA_real_ else max(as.numeric(x), na.rm = TRUE) )
      temp_dt$`Theoritical max. coverage [%]` <- sapply(strsplit(as.character(temp_dt$`Theoritical max coverage per protein ID [%]`),";"), function(x) {
        x1 <- suppressWarnings(as.numeric(x))
        if(all(is.na(x1))) NA_real_ else max(x1, na.rm = TRUE)
      })
      # write.table(temp_dt[order(as.numeric(temp_dt$id)),c("id", "Protein IDs", "Sequence coverage [%]", "Theoritical max coverage per protein ID [%]", "Theoritical max. coverage [%]")], "C:\\Users\\grigo\\Desktop\\SICyLIA_v.2.0\\myShinyApp\\tbl.txt", quote=F, sep = "\t", row.names=F, col.names=T)
      temp_dt[,c("id", "Protein IDs", "Sequence coverage [%]", "Theoritical max coverage per protein ID [%]", "Theoritical max. coverage [%]")]
      tbl <- temp_dt
      
      }
    return(tbl)
  })
  

  non_identifiable_proteins <- reactive({
    update_busy_bar(20)
    df <- digestion_init()
    update_busy_bar(30)
    df2 <- proteinStringSet_splitter()
    df_pep <- identifiable_peptides()
    update_busy_bar(40)
    showNotification(paste0("Applying digestion parameters to ", nrow(df2), " proteins..."), type="message", duration = 3)

    df$identifiable <- ifelse(df$length >= min_peptide_length$choice & df$`mass [Da]` < max_peptide_mass$choice, 1, 0)

    all_uniProtIDs <- unique(df[,c("UniProtID")])
    identifiable_uniProtIDs <- unique(df_pep$UniProtID)
    non_identifiable_uniProtIDs <- setdiff(all_uniProtIDs, identifiable_uniProtIDs)
    non_identifiable_uniProtIDs <- as.data.frame(non_identifiable_uniProtIDs)
    names(non_identifiable_uniProtIDs) <- c("UniProtID")
    tbl <- non_identifiable_uniProtIDs
    update_busy_bar(50)
    coord <- df[,c("UniProtID","start","end","identifiable")]
    dt_coordinates <- group_by(coord, UniProtID) %>%
      summarise(start = paste(start, collapse = ","),
                end = paste(end, collapse = ","),
                identifiable = paste(identifiable, collapse = ","))
    update_busy_bar(60)
    showNotification("Calculating sequence coverage ...", type="message", duration = 3)

    grouped_cleaved <- as.data.frame(table(df[,c("UniProtID", "identifiable")]))
    update_busy_bar(70)
    grouped_non_identifiable <- grouped_cleaved[grouped_cleaved$identifiable ==0,]
    grouped_non_identifiable <- grouped_non_identifiable[,c("UniProtID","Freq")]
    names(grouped_non_identifiable) <- c("UniProtID", "Non-identifiable")
  
    grouped_identifiable <- grouped_cleaved[grouped_cleaved$identifiable ==1,]
    grouped_identifiable <- grouped_identifiable[,c("UniProtID","Freq")]
    names(grouped_identifiable) <- c("UniProtID", "Identifiable")
    update_busy_bar(80)
    showNotification("Generating final output ...", type="message", duration = 3)

    merged_dt <- merge(tbl, df2[,c("sequence","UniProtID","GeneName","ProteinName","ProteinLength")], by = "UniProtID")
    # merged_dt <- merged_dt[,c("sequence","UniProtID","GeneName","ProteinName","ProteinLength")]

    merged_dt <- merge(merged_dt, grouped_non_identifiable, by = "UniProtID")
    # merged_dt <- merge(merged_dt, grouped_identifiable, by = "UniProtID")
    merged_dt$`Total peptides` <- merged_dt$`Non-identifiable`

    merged_dt <- merge(merged_dt, dt_coordinates, by = "UniProtID")
    update_busy_bar(90)
    
    if (nrow(merged_dt)>0) {
      merged_dt$links4seq <- createLink4seq(merged_dt[, c("UniProtID")], merged_dt[, c("sequence")], merged_dt[, c("GeneName")], merged_dt[, c("ProteinName")], merged_dt[, c("start")], merged_dt[, c("end")], merged_dt[, c("identifiable")])
      merged_dt$coverage <- 0
      merged_dt$Identifiable <- 0
      ret <- merged_dt[, c("UniProtID","GeneName","ProteinName", "Identifiable", "Total peptides", "ProteinLength", "coverage", "links4seq")]
      names(ret) <- c("UniProtID","Gene name","Protein name","Identifiable peptides", "Total peptides", "Protein length", "Theoritical max. coverage [%]","")
    }
    else {
      ret <- data.table(UniProtID=character(), `Gene name`=character(), `Protein name`=character(), `Identifiable peptides`=character(),  `Total peptides`=character(),  `Protein length`=character(),  `Theoritical max. coverage [%]`=character(), `Sequence`=character() )
    }
    update_busy_bar(100)
    return(ret)
  })
  

    identifiable_proteins <- reactive({
    update_busy_bar(20)
    dt <- digestion_init()
    update_busy_bar(30)
    
    update_busy_bar(40)
    dt2 <- proteinStringSet_splitter()
    
    showNotification(paste0("Applying digestion parameters to ", nrow(dt2), " proteins..."), type="message", duration = 3)
    dt$identifiable <- ifelse(dt$length >= min_peptide_length$choice & dt$`mass [Da]` < max_peptide_mass$choice, 1, 0)
    
    update_busy_bar(50)
    coord <- dt[,c("UniProtID","start","end","identifiable")]
    dt_coordinates <- group_by(coord, UniProtID) %>%
        summarise(start = paste(start, collapse = ","), 
                end = paste(end, collapse = ","),
                identifiable = paste(identifiable, collapse = ","))
    
    update_busy_bar(60)
    showNotification("Calculating sequence coverage ...", type="message", duration = 3)
    grouped_cleaved <- as.data.frame(table(dt[,c("UniProtID","identifiable")]))
    
    update_busy_bar(70)
    grouped_non_identifiable <- grouped_cleaved[grouped_cleaved$identifiable ==0,]
    grouped_non_identifiable <- grouped_non_identifiable[,c("UniProtID","Freq")]
    names(grouped_non_identifiable) <- c("UniProtID", "Non-identifiable")
    
    grouped_identifiable <- grouped_cleaved[grouped_cleaved$identifiable ==1,]
    grouped_identifiable <- grouped_identifiable[,c("UniProtID","Freq")]
    names(grouped_identifiable) <- c("UniProtID", "Identifiable")
    
    update_busy_bar(80)
    num <- aggregate(dt[dt$identifiable ==1,]$length, by=list(Category=dt[dt$identifiable ==1,]$UniProtID), FUN=sum)
    names(num) <- c("UniProtID","length")
    unq_identifiable_prot_ids <- unique(dt[dt$identifiable ==1,]$UniProtID)
    # # denom <- aggregate(cleaved_dt[cleaved_dt$UniProtID %in% unq_identifiable_prot_ids,]$length, by=list(Category=cleaved_dt[cleaved_dt$UniProtID %in% unq_identifiable_prot_ids,]$UniProtID), FUN=sum)
    denom <- dt2[dt2$UniProtID %in% unq_identifiable_prot_ids, c("UniProtID","ProteinLength")]
    names(denom) <- c("UniProtID","length")
    seq_coverage <- merge(num, denom, by= c("UniProtID"))
    seq_coverage$coverage <- round(seq_coverage$length.x/seq_coverage$length.y*100 ,2)
    seq_coverage <- seq_coverage[,c("UniProtID","coverage")]
    
    update_busy_bar(90)
    showNotification("Generating final output ...", type="message", duration = 3)
    tbl <- subset(dt, identifiable==1)
    tbl <- as.data.frame(unique(tbl[, c("UniProtID")]))
    names(tbl) <- c("UniProtID")
    merged_dt <- merge(tbl, dt2, by = "UniProtID")
    
    merged_dt <- merge(merged_dt, grouped_non_identifiable, by = "UniProtID")
    merged_dt <- merge(merged_dt, grouped_identifiable, by = "UniProtID")
    merged_dt$`Total peptides` <- merged_dt$Identifiable + merged_dt$`Non-identifiable`
    
    merged_dt <- merge(merged_dt, dt_coordinates, by = "UniProtID")
    
    merged_dt <- merge(merged_dt, seq_coverage, by = "UniProtID")
    
    if (nrow(merged_dt)>0) {

      merged_dt$links4seq <- createLink4seq(merged_dt[, c("UniProtID")], merged_dt[, c("sequence")], merged_dt[, c("GeneName")], merged_dt[, c("ProteinName")], merged_dt[, c("start")], merged_dt[, c("end")], merged_dt[, c("identifiable")])
      
      ret <- merged_dt[, c("UniProtID","GeneName","ProteinName", "Identifiable", "Total peptides", "ProteinLength", "coverage", "links4seq")]
      names(ret) <- c("UniProtID","Gene name","Protein name","Identifiable peptides", "Total peptides", "Protein length", "Theoritical max. coverage [%]","")
    }
    else {
      ret <- data.table(UniProtID=character(), `Gene name`=character(), `Protein name`=character(), `Identifiable peptides`=character(),  `Total peptides`=character(),  `Protein length`=character(),  `Theoritical max. coverage [%]`=character(), `Sequence`=character() )
    }
    update_busy_bar(100)
    return(ret)
  })
  
  
  identifiable_peptides <- reactive({
    dt <- digestion_init()
    dt$identifiable <- ifelse(dt$length >= min_peptide_length$choice & dt$`mass [Da]` < max_peptide_mass$choice, 1, 0)
    tbl <- subset(dt, identifiable==1)
    tbl <- tbl[, c("UniProtID", "peptide", "start", "end", "length", "mass [Da]")]
    names(tbl) <- c("UniProtID", "Peptide", "Start", "End", "Length", "Peptide mass [Da]")
    return(tbl)
  })
  
  
  frequent_peptides <- reactive({
    update_busy_bar(0)
    dt <- identifiable_peptides()
    update_busy_bar(50)
    tbl <- as.data.frame(dt[,c("UniProtID", "Peptide")] %>%
      group_by(Peptide) %>%
      summarise(Frequency=n(), `Unique Proteins`=n_distinct(UniProtID)) 
      # %>% arrange(desc(Frequency)) 
      )
    update_busy_bar(100)
    return(tbl)
  })
  
  
  dg <- reactive({
    update_busy_bar(0)
    showNotification("Digesting proteins ...", type="message", duration = 3)
    if (input$switchButtonDigest == TRUE) {
      dt.seq <- fastafile_reader()
    }
    else if (input$switchButtonDigest == FALSE) {
      dt.seq <- proteinList_reader()
    }
    digested <- digestedPep(dt.seq)
    digested$UniProtID <- stri_match_first_regex(substr(digested$group_name, 4, nchar(digested$group_name)), "(.*?)\\|")[,2] 
    return(digested)
  })
  
  
  cl <- reactive({
    update_busy_bar(0)
    if (input$switchButtonDigest == TRUE) {
      dt.seq <- fastafile_reader()
    }
    else if (input$switchButtonDigest == FALSE) {
      dt.seq <- proteinList_reader()
    }
    cl <- clRanges(dt.seq) 
    return(cl)
  })
  
  
  digestion_init <- reactive({
    cl_dt <- merge(dg()[,c("ID","UniProtID","value")], cl()[,c("ID","start","end","width")], by = "ID")
    names(cl_dt) <- c("ID", "UniProtID", "peptide", "start", "end", "length")
    showNotification("Calculating molecular mass of peptides ...", type="message", duration = 3)
    cl_dt$`mass [Da]` <- sapply(strsplit(cl_dt$peptide, ""), function(x) round(sum(mol_weightAA[x]) + mol_weightAA["H2O"], 4))
    return(cl_dt)
  })
  
  
  fastafile_reader <- reactive({
    req(input$fastaFile)
    ext <- tools::file_ext(input$fastaFile$datapath)
    validate(need(ext == "fasta", "Please upload a fasta file"))
    dt.seq <- readAAStringSet(input$fastaFile$datapath)
    return(dt.seq)
  })
  
  #################################
  #################################
  observeEvent(input$submitProtList, {
    if (input$proteinListTextInput != "") {
      inputIDs <- strsplit(input$proteinListTextInput, "\\s+")[[1]]
      if (length(inputIDs) > 0 & length(inputIDs) <= 200) {
        proteinList_reader()
        digestion_init()
        proteinStringSet_splitter()
        shinyjs::js$enableTab("digestion_results_tab")
        updateTabsetPanel(session, "digestion_tabSet", selected = "digestion_results_tab")
        updateVerticalTabsetPanel(session, "vertical_results_tabSet", selected = "vertical_entries_tab")
        shinyjs::addClass(selector = "body", class = "sidebar-collapse")
      }
      else if (length(inputIDs) > 200) {
        showNotification("Up to 200 proteins can be submitted simultaneously. Please upload a fasta file instead", type="warning", duration = 3)
      }
    }
    else {
      showNotification("Please submit at least one protein ID", type="warning", duration = 3)
    }
  })
  
  
  proteinList_reader <- reactive({
    AAset<-AAStringSet("")
    header_c <- c("")
    j <- 1
    inputIDs <- unique(unlist(strsplit(input$proteinListTextInput, "\\s+")[[1]]))
    show_modal_progress_line(value=0, easing = "linear", text = "Retrieving data from UniProt. This step might take a while, please wait...")
    for (i in 1:length(inputIDs)) {
      update_modal_progress(value = i / length(inputIDs),
                            text = paste0("Fetching protein information from UniProt: ", inputIDs[i], " (", i, "/",length(inputIDs), ")" ))
      r <- GET(paste0("https://www.uniprot.org/uniprot/", inputIDs[i], ".fasta"))
      if (r$status_code==200) {
        content <- content(r, as="text", encoding="UTF-8")
        header <- substr(content, 1, regexpr("\n", content)-1)
        sequen <- substr(content, regexpr("\n", content)[1], nchar(content))
        sequen <- gsub("\n", "", sequen)
        sequen <- gsub("\\|", "", sequen)
        header <- gsub(">", "", header)
        AAset[j] <- AAStringSet(sequen)
        header_c[j] <- header
        j = j + 1
      }
      else{
        showNotification(paste0("Accession ID ", inputIDs[i], " returned a bad request"), type="error", duration = 4)
      }
    }
    names(AAset) <- header_c
    remove_modal_progress()
    return(AAset)
  })
  
  
  
  proteinStringSet_splitter <- reactive({
    if (input$switchButtonDigest == TRUE) {
      dt.seq <- fastafile_reader()
    }
    else if (input$switchButtonDigest == FALSE) {
      dt.seq <- proteinList_reader()
    }
    df.seq <- as.data.frame(dt.seq)
    df.seq$header <- names(dt.seq)
    df.seq$UniProtID <- stri_match_first_regex(substr(df.seq$header, 4, nchar(df.seq$header)), "(.*?)\\|")[,2]
    df.seq$GeneName <- sub(" .*", "", gsub(".*GN=","",df.seq$header))
    df.seq$ProteinName <- substr(stri_match_first_regex(sub(".*? ", "", df.seq$header), "(.*?)\\=")[,2],1,nchar(stri_match_first_regex(sub(".*? ", "", df.seq$header), "(.*?)\\=")[,2])-3)
    df.seq$ProteinLength <- nchar(as.character(df.seq$x))
    names(df.seq) <- c("sequence","header","UniProtID","GeneName","ProteinName","ProteinLength")
    df.seq <- df.seq[ , !(names(df.seq) %in% "header")]
    rownames(df.seq) <- NULL
    return(df.seq)
  })
  
  
  bulk_digestion <- reactive({
    update_busy_bar(0)
    show_modal_progress_line(value=0, easing = "linear", text = "Bulk digestion in progress. This step might take a while, please wait...")
    protNum <- rep(NA, length(enzyme_list))
    peptNum <- rep(NA, length(enzyme_list))
    seqCover <- rep(NA, length(enzyme_list))
    if (input$switchButtonDigest == TRUE) {
      dt.seq <- fastafile_reader()
    }
    else if (input$switchButtonDigest == FALSE) {
      dt.seq <- proteinList_reader()
    }
    if (length(dt.seq) == 0) {
      remove_modal_progress()
      update_busy_bar(100)
    }
    for (i in 1:length(enzyme_list)) {
      update_modal_progress(value = i / length(enzyme_list),
                            text = paste0("Enzyme evaluation: ", names(enzyme_list)[i], " (", i, "/", length(enzyme_list), ")" ))
      dg2 <- cleave(dt.seq, enzym=enzyme_list[[i]], missedCleavages = 0, custom = NULL, unique = FALSE)
      dg2 <- as.data.frame(dg2)
      dg2$mol.weight <- sapply(strsplit(dg2$value, ""), function(x) sum(mol_weightAA[x]))
      dg2$identifiable <- ifelse(nchar(dg2$value) >= min_peptide_length$choice & dg2$mol.weight < max_peptide_mass$choice, 1, 0)
      v <- subset(dg2, identifiable==1)
      
      PepLen <- do.call(sum, lapply(nchar(v$value), sum))
      ProtLen <- do.call(sum, lapply(width(dt.seq), sum))
      
      protNum[i] <- length(unique(v$group))
      peptNum[i] <- nrow(v)
      seqCover[i] <- round(PepLen/ProtLen*100,2)
    }
    temp <- do.call(rbind.data.frame, Map('c', names(enzyme_list), protNum, peptNum, seqCover))
    result <- as.data.frame(temp)
    result[,c(2)] = as.numeric(as.character(result[,c(2)]))
    result[,c(3)] = as.numeric(as.character(result[,c(3)]))
    result[,c(4)] = as.numeric(as.character(result[,c(4)]))
    # result <- cbind(as.character(temp[,1]), as.numeric(as.character(temp[,2])), as.numeric(as.character(temp[,3])))
    names(result) <- c("Enzyme","Identifiable proteins","Identifiable peptides","Max. theoritical sequence coverage [%]")
    # print(result)
    remove_modal_progress()
    update_busy_bar(100)
    return(result)
  })
  
  #######################################
  ############MY FUNCTIONS###############
  createLink4seq <- function(unipID, seq, geneName, protName, start, end, identifiable) {
    link <- paste("<a href='#sequence-viewer' onclick='shinyjs.seque(\"",unipID,"\", \"",seq,"\", \"",geneName,"\", \"",protName,"\", \"",start,"\", \"",end,"\", \"",identifiable,"\");' class='btn btn-custom btn-sm'>Show sequence</a>", sep="")
    return(link)
  }
  
  digestedPep <- function(seq) { #as.numeric(miscleavage$choice)
    dpep <- cleave(seq, enzym = protease$choice, missedCleavages=0, custom = NULL, unique = FALSE)
    dpep <- as.data.frame(dpep)
    dpep$ID <- seq.int(nrow(dpep))
    return(dpep)
  }
  
  clRanges <- function(seq) {
    clRan <- cleavageRanges(seq, enzym = protease$choice, missedCleavages=0, custom = NULL)
    clRan <- as.data.frame(clRan)
    clRan$ID <- seq.int(nrow(clRan))
    return(clRan)
  }
  #####################################################################################################################
  #####################################################################################################################
  
  dt_peptides <- reactive({
    inFile_peptides <- input$peptides
    req(inFile_peptides)
    withProgress(message = 'Validating peptides, please wait...', {
      incProgress(1/1)
      tbl <- read.table(inFile_peptides$datapath, header = TRUE, fill = TRUE, quote = "", sep = "\t")
      required_columns <- c("id","Leading.razor.protein","Start.position","End.position","C.Count","Length")
      if(all(required_columns %in% names(tbl)) == F) {
        tbl <- NULL
        shinyjs::runjs("$('#peptides_progress').css('visibility','hidden')")
        js$swalErrorAlert()
      }
      else{
        tbl <- tbl[,c("id","Leading.razor.protein","Start.position","End.position","C.Count","Length")]
        ##replace '.' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "\\.", replacement = "_")
        ##replace '__' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "__", replacement = "_")
        ##replace tailing '_' symbol to '' in column names##
        names(tbl) <- gsub("_$", "", names(tbl))
        showNotification(paste("The file ", basename(inFile_peptides$name), " uploaded successfully!"), type="message", duration = 3)
      }
      return (tbl)
    })
  })
  
  
  dt_proteinGroups <- reactive({
    inFile_protein_groups <- input$proteinGroups
    req(inFile_protein_groups)
    withProgress(message = 'Validating protein groups, please wait...', {
      incProgress(1/1)
      tbl <- read.table(inFile_protein_groups$datapath, header = TRUE, fill = TRUE, quote = "", sep = "\t") 
      rep_int_corr_columnNames <- colnames(tbl[,grepl("Reporter.intensity.corrected.", colnames(tbl))])
      required_columns <- c(rep_int_corr_columnNames, "Mod..peptide.IDs", "Unique.peptides")
      if(all(required_columns %in% names(tbl)) == F) {
        tbl <- NULL
        shinyjs::runjs("$('#proteinGroups_progress').css('visibility','hidden')")
        js$swalErrorAlert()
      }
      else{
        tbl <- tbl[,c(rep_int_corr_columnNames, "Mod..peptide.IDs", "Unique.peptides")]
        #replace reporter to protein
        names(tbl) <- gsub(x = names(tbl), pattern = "Reporter.intensity.corrected.", replacement = "Protein.intensity.corrected.")
        ##replace '.' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "\\.", replacement = "_")
        ##replace '__' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "__", replacement = "_")
        ##replace tailing '_' symbol to '' in column names##
        names(tbl) <- gsub("_$", "", names(tbl))
        showNotification(paste("The file ", basename(inFile_protein_groups$name), " uploaded successfully!"), type="message", duration = 3)
      }
      return (tbl)
    })
  })
  
  
  dt_modSpecPeptides <- reactive({
    inFile_modSpecPept <- input$modSpecPeptides
    req(inFile_modSpecPept)
    withProgress(message = 'Validating modification specific peptides, please wait...', {
      incProgress(1/1)
      tbl <- read.table(inFile_modSpecPept$datapath, header = TRUE, fill = TRUE, quote = "", sep = "\t")
      rep_int_corr_columnNames <- colnames(tbl[,grepl("Reporter.intensity.corrected.", colnames(tbl))])
      required_columns <- c(rep_int_corr_columnNames,"Modifications","Raw.file","Reverse","Potential.contaminant","Unique..Groups.","Unique..Proteins.","Sequence","Proteins","Gene.Names","Protein.Names","Acetyl..Protein.N.term.","Carbamidomethyl..C.","Carbamidomethyl.Heavy..C..Var.","Oxidation..M.","id","Peptide.ID","Evidence.IDs","Carbamidomethyl..C..site.IDs","Carbamidomethyl.Heavy..C..Var..site.IDs","Mass","Retention.time","PEP","Score","Delta.score","Intensity","MS.MS.Count","Protein.group.IDs","Protein.Groups","Missed.cleavages","Charges")
      if(all(required_columns %in% names(tbl)) == F) {
        # print(setdiff(required_columns, names(tbl)))
        tbl <- NULL
        shinyjs::runjs("$('#modSpecPeptides_progress').css('visibility','hidden')")
        js$swalErrorAlert()
      }
      else{
        tbl <- tbl[,c(rep_int_corr_columnNames,"Modifications","Raw.file","Reverse","Potential.contaminant","Unique..Groups.","Unique..Proteins.","Sequence","Proteins","Gene.Names","Protein.Names","Acetyl..Protein.N.term.","Carbamidomethyl..C.","Carbamidomethyl.Heavy..C..Var.","Oxidation..M.","id","Peptide.ID","Evidence.IDs","Carbamidomethyl..C..site.IDs","Carbamidomethyl.Heavy..C..Var..site.IDs","Mass","Retention.time","PEP","Score","Delta.score","Intensity","MS.MS.Count","Protein.group.IDs","Protein.Groups","Missed.cleavages","Charges")]
        ##replace '.' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "\\.", replacement = "_")
        ##replace '__' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "__", replacement = "_")
        ##replace tailing '_' symbol to '' in column names##
        names(tbl) <- gsub("_$", "", names(tbl))
        showNotification(paste("The file ", basename(inFile_modSpecPept$name), " uploaded successfully!"), type="message", duration = 3)
      }
      return (tbl)
    })
  })
  
  
  dt_lightPeptides <- reactive({
    inFile_lightPeptides <- input$lightPeptides
    req(inFile_lightPeptides)
    withProgress(message = 'Validating light-labeled sites, please wait...', {
      incProgress(1/1)
      tbl <- read.table(inFile_lightPeptides$datapath, header = TRUE, fill = TRUE, quote = "", sep = "\t")
      required_columns <- c("id", "Proteins","Positions.within.proteins","Positions","Position","Position.in.peptide","Score.diff","Localization.prob")
      if(all(required_columns %in% names(tbl)) == F) {
        tbl <- NULL 
        shinyjs::runjs("$('#lightPeptides_progress').css('visibility','hidden')")
        js$swalErrorAlert()
      }
      else{
        tbl <- tbl[,c("id", "Proteins","Positions.within.proteins","Positions","Position","Position.in.peptide","Score.diff","Localization.prob")]
        names(tbl) <- paste(names(tbl), "L", sep = ".")
        ##replace '.' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "\\.", replacement = "_")
        ##replace '__' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "__", replacement = "_")
        ##replace tailing '_' symbol to '' in column names##
        names(tbl) <- gsub("_$", "", names(tbl))
        showNotification(paste("The file ", basename(inFile_lightPeptides$name), " uploaded successfully!"), type="message", duration = 3)
      }
      return (tbl)
    })
  })
  
  
  dt_heavyPeptides <- reactive({
    inFile_heavyPeptides <- input$heavyPeptides
    req(inFile_heavyPeptides)
    withProgress(message = 'Validating heavy-labeled sites, please wait...', {
      incProgress(1/1)
      tbl <- read.table(inFile_heavyPeptides$datapath, header = TRUE, fill = TRUE, quote = "", sep = "\t")
      required_columns <- c("id", "Proteins","Positions.within.proteins","Positions","Position","Position.in.peptide","Score.diff","Localization.prob")
      if(all(required_columns %in% names(tbl)) == F) {
        tbl <- NULL 
        shinyjs::runjs("$('#heavyPeptides_progress').css('visibility','hidden')")
        js$swalErrorAlert()
      }
      else{
        tbl <- tbl[,c("id", "Proteins","Positions.within.proteins","Positions","Position","Position.in.peptide","Score.diff","Localization.prob")]
        names(tbl) <- paste(names(tbl), "H", sep = ".")
        ##replace '.' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "\\.", replacement = "_")
        ##replace '__' symbol to '_' in column names##
        names(tbl) <- gsub(x = names(tbl), pattern = "__", replacement = "_")
        ##replace tailing '_' symbol to '' in column names##
        names(tbl) <- gsub("_$", "", names(tbl))
        showNotification(paste("The file ", basename(inFile_heavyPeptides$name), " uploaded successfully!"), type="message", duration = 3)
      }
      return (tbl)
    })
  })
  
  
  dt_main_dt <- reactive({
    if (!is.null(dt_peptides()) & !is.null(dt_proteinGroups()) & !is.null(dt_modSpecPeptides()) & !is.null(dt_lightPeptides()) & !is.null(dt_heavyPeptides())){
      withProgress(message = 'Calculation in progress. This step might take a while, please wait...', {
        incProgress(1/1)
        
        ##left outer join peptides.txt and ModSpecPeptides.txt##
        df2 <- left_join(dt_modSpecPeptides(), dt_peptides(), by = c("Peptide_ID" = "id"), na_matches = "never")
        
        ##remove those without Cysteine
        df2 <- df2[df2$C_Count > 0,]
        
        ##remove Reverse
        df2 <- df2[df2$Reverse != "+",]
        
        ##preprocess Protein Groups##
        ##keep unique peptides > 0##
        proteinGroups_data <- dt_proteinGroups()
        proteinGroups_data <- proteinGroups_data[proteinGroups_data$Unique_peptides > 0,]
        
        ##expand proteinGroups_data in order to join table with df2## 
        proteinGroups_data$Mod_peptide_IDs <- sub("^$", "-1", proteinGroups_data$Mod_peptide_IDs)
        x <- proteinGroups_data$Mod_peptide_IDs
        proteinGroups_data["number_of_mod_pep_IDs"] <- data.frame(matrix(unlist(lengths(regmatches(x, gregexpr(";", x))) + 1), nrow=length(lengths(regmatches(x, gregexpr(";", x)))), byrow=T))
        mod_pept_IDs <- unlist(strsplit(as.character(x), split=";"))
        proteinGroups_data <- proteinGroups_data[rep(row.names(proteinGroups_data), proteinGroups_data$number_of_mod_pep_IDs), 1:ncol(proteinGroups_data) - 1]
        proteinGroups_data["Mod_peptide_ID"] <- data.frame(matrix(as.numeric(mod_pept_IDs)))
        colNamesIntensCorr <- colnames(proteinGroups_data[,grepl("Protein_intensity_corrected_", colnames(proteinGroups_data))])
        
        df3 <- left_join(df2, proteinGroups_data[,c(colNamesIntensCorr, "Mod_peptide_ID")], by = c("id" = "Mod_peptide_ID"))
        
        #expand protein groups table again for heavy sites (this step should become dynamic)
        #fill in the blank cells with zero in order to expand properly
        df3$Carbamidomethyl_Heavy_C_Var_site_IDs <- sub("^$", "-1", df3$Carbamidomethyl_Heavy_C_Var_site_IDs)
        x <- df3$Carbamidomethyl_Heavy_C_Var_site_IDs
        df3["number_of_heavy_site_IDs"] <- data.frame(matrix(unlist(lengths(regmatches(x, gregexpr(";", x))) + 1), nrow=length(lengths(regmatches(x, gregexpr(";", x)))), byrow=T))
        heavy_site_IDs <- unlist(strsplit(as.character(x), split=";"))
        df3 <- df3[rep(row.names(df3), df3$number_of_heavy_site_IDs), 1:ncol(df3) - 1]
        df3["Heavy_site_ID"] <- data.frame(matrix(as.numeric(heavy_site_IDs)))
        
        df4 <- left_join(df3, dt_heavyPeptides(), by = c("Heavy_site_ID" = "id_H"), na_matches = "never")
        
        #expand protein groups table again for light sites (this step should become dynamic)
        #fill in the blank cells with zero in order to expand properly
        #names(LightFile_data) <- gsub(x = names(LightFile_data), pattern = "\\.", replacement = "_")
        df4$Carbamidomethyl_C_site_IDs <- sub("^$", "-1", df4$Carbamidomethyl_C_site_IDs)
        x <- df4$Carbamidomethyl_C_site_IDs
        df4["number_of_light_site_IDs"] <- data.frame(matrix(unlist(lengths(regmatches(x, gregexpr(";", x))) + 1), nrow=length(lengths(regmatches(x, gregexpr(";", x)))), byrow=T))
        light_site_IDs <- unlist(strsplit(as.character(x), split=";"))
        df4 <- df4[rep(row.names(df4), df4$number_of_light_site_IDs), 1:ncol(df4) - 1]
        df4["Light_site_ID"] <- data.frame(matrix(as.numeric(light_site_IDs)))
        
        df5 <- left_join(df4, dt_lightPeptides(), by = c("Light_site_ID" = "id_L"), na_matches = "never")
        
        df5 <- within(df5, Heavy_site_ID[Heavy_site_ID == '-1'] <- '')
        df5 <- within(df5, Light_site_ID[Light_site_ID == '-1'] <- '')
        df5 <- within(df5, Carbamidomethyl_Heavy_C_Var_site_IDs[Carbamidomethyl_Heavy_C_Var_site_IDs == '-1'] <- '')
        df5 <- within(df5, Carbamidomethyl_C_site_IDs[Carbamidomethyl_C_site_IDs == '-1'] <- '')
        
        df6 <- as.data.frame(df5) %>% 
          group_by(id) %>% 
          mutate_at(., vars(starts_with("Positions_"),starts_with("Position_"),ends_with("_site_ID")), function(x) ifelse(x=='', NA, paste(x, collapse=";"))  ) %>% 
          distinct(.)
        
        df6 <- sapply(df6, as.character)
        df6[is.na(df6)] <- "" 
        
        
        return(df6)
        
      })
    }
  })
  
  
  mainFunc <- function(exp_tbl) {
    
    df7 <- dt_main_dt()[apply(dt_main_dt()[,grep("Protein_intensity_corrected", colnames(dt_main_dt()))], 1, function(z) any(z==0)), ] 
    
    df8 <- dt_main_dt()[apply(dt_main_dt()[,grep("Protein_intensity_corrected", colnames(dt_main_dt()))], 1, function(z) any(z!=0)), ] 
    
    ix <- grep("Reporter_intensity_corrected", colnames(df8))
    iy <- grep("Protein_intensity_corrected", colnames(df8))
    df8a <- data.frame(matrix(unlist(as.numeric(df8[,ix]) / as.numeric(df8[,iy])), nrow=nrow(df8), byrow=F), stringsAsFactors=FALSE)
    colnames(df8a) <- paste("Normalised_intensity_corrected", seq(0, length = ncol(df8a)), sep = "_")
    df8 <- cbind(df8, df8a)
    
    #remove columns 'Reporter_intensity_corrected' and 'Protein_intensity_corrected'
    df8 <- select(df8, -contains("Reporter_intensity_corrected"))
    df8 <- select(df8, -contains("Protein_intensity_corrected"))
    
    ##divide each normalised intensity with the median per column and logarithmise##
    keep_names <- c(colnames(df8)[grep("Normalised_intensity_corrected",colnames(df8))])
    inds <- which(colnames(df8)%in%keep_names)
    df8[inds] <- df8[inds] / colMedians(data.matrix(df8[inds]), na.rm = TRUE)
    df8[inds] <- log(df8[inds] * exp(9),2)
    
    ##remove infinity values##
    #df8 <- df8[is.finite(rowSums(df8[,grep("_intensity_corrected_", colnames(df8))])),]
    ##replace inf values with NA##
    df8[sapply(df8, is.infinite)] <- NA
    
    df9 <- df8[(df8$Unique_Groups == 'no' & df8$Unique_Proteins == 'no') ,]
    df8 <- df8[((df8$Unique_Groups == 'yes' & df8$Unique_Proteins == 'yes') | (df8$Unique_Groups == 'yes' & df8$Unique_Proteins == 'no')) ,]
    
    
    if (exp_tbl == "1") {
      df <- df7
    }
    else if (exp_tbl == "2") {
      df <- df9
    }
    else if (exp_tbl == "3") {
      df <- df8
    }
    
    return (df)
    
}


  rSq <- function (x, y){
    df <- cor(x, y) ^ 2
    return(df)
  } 


   table_1 <- reactive ({
     mainFunc("1")
   })
   
   table_2 <- reactive ({
     mainFunc("2")
   })
   
   table_3 <- reactive ({
     mainFunc("3")
   })
  
   
     output$resultsmenuitem <- renderMenu({
       if (!is.null(dt_main_dt())){
         shinyalert(
           title = "Calculation completed",
           text = "The estimation of parameters was successful. <br /> Click on the 'Results' tab in the left panel to browse data or select 'Graphs' to check the quality of the uploaded data.",
           closeOnEsc = TRUE,
           closeOnClickOutside = TRUE,
           html = TRUE,
           type = "success",
           showConfirmButton = TRUE,
           showCancelButton = FALSE,
           confirmButtonText = "OK",
           confirmButtonCol = "#95a5a6",
           timer = 0,
           imageUrl = "",
           animation = TRUE
         )
         menuItem("Results", tabName = "resultsmenuitem", icon = icon("table"))
       }
     })
   
  
  
  output$resetsmenuitem <- renderMenu({
    if (!is.null(dt_main_dt())){
      menuItem("Reset session", tabName = "resetsmenuitem", icon = icon("redo"))
    }
  })
  
  
  output$qualitycontrolmenuitem <- renderMenu({
    if (!is.null(dt_main_dt())){
      menuItem("Graphs", tabName = "qualitycontrolmenuitem", icon = icon("chart-bar"))
    }
  })
  
  
  
  output$downloadBtn1 <- renderUI({
    #if(!is.null(DT::dataTableOutput("table_1"))) {
      downloadButton("downloadTable1", "Download")
    #}
  })
  
  output$downloadOptns1 <- renderUI({
    radioButtons("downloadOptions1","Specify output format", choiceNames= c("Standard", "Perseus"), choiceValues = c("std", "perseus"), inline = TRUE)
  })
  
  output$downloadBtn2 <- renderUI({
    #if(!is.null(DT::dataTableOutput("table_2"))) {
    downloadButton("downloadTable2", "Download")
    #}
  })
  
  output$downloadOptns2 <- renderUI({
    radioButtons("downloadOptions2","Specify output format", choiceNames= c("Standard", "Perseus"), choiceValues = c("std", "perseus"), inline = TRUE)
  })
  
  output$downloadBtn3 <- renderUI({
    #if(!is.null(DT::dataTableOutput("table_3"))) {
    downloadButton("downloadTable3", "Download")
    #}
  })
  
  output$downloadOptns3 <- renderUI({
    radioButtons("downloadOptions3","Specify output format", choiceNames= c("Standard", "Perseus"), choiceValues = c("std", "perseus"), inline = TRUE)
  })
  
###############################################
  #############################################
  output$identifiable_prot <- DT::renderDataTable({
    js$hidesequeNmolart()
    withProgress(message = 'Computation in progress. This step might take a while depending on the the number of proteins and configuration parameters.', {
      incProgress(1/1)
      identifiable_proteins()
    })
  }, escape = FALSE, selection = 'single',
  options = list(columnDefs = list(list(className = 'dt-center', targets = c(3,4,5,6))), searchHighlight = TRUE, hover = FALSE, scrollX = TRUE, dom = 'lfrtip', pageLength = 10, lengthMenu = c(10, 25, 50, 100)), rownames = FALSE)
  
  
  output$identifiable_pept <- DT::renderDataTable({
    js$hidesequeNmolart()
    withProgress(message = 'Computation in progress. This step might take a while depending on the number of proteins and configuration parameters.', {
      incProgress(1/1)
      identifiable_peptides()
      })
  }, escape = FALSE, selection = 'single',  
  options = list(columnDefs = list(list(className = 'dt-center', targets = c(3,4,5))), searchHighlight = TRUE, hover = FALSE, scrollX = TRUE, dom = 'lfrtip', pageLength = 10, lengthMenu = c(10, 25, 50, 100)), rownames = FALSE)
 
  output$download_identifiable_pept <- renderUI({
    downloadButton("download_iden_pept", "Download as TSV")
  })
  
  output$download_iden_pept <- downloadHandler(
    filename = function() {
      paste("Identifiable_peptides_", Sys.Date(), ".tsv", sep="")
    },
    content = function(file) {
        write.table(identifiable_peptides(), file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
  )
  
  output$frequent_pept <- DT::renderDataTable({
    js$hidesequeNmolart()
    withProgress(message = 'Computation in progress. This step might take a while depending on the number of proteins and configuration parameters.', {
      incProgress(1/1)
      frequent_peptides()
    })
  }, escape = FALSE, selection = 'single',
  options = list(columnDefs = list(list(className = 'dt-center', targets = c(1,2))), order = list(1, 'desc'),searchHighlight = TRUE, hover = FALSE, scrollX = TRUE, dom = 'lfrtip', pageLength = 10, lengthMenu = c(10, 25, 50, 100)), rownames = FALSE)
  
  
  output$non_identifiable_prot <- DT::renderDataTable({
    js$hidesequeNmolart()
    withProgress(message = 'Computation in progress. This step might take a while depending on the the number of proteins and configuration parameters.', {
      incProgress(1/1)
      non_identifiable_proteins()
    })
  }, escape = FALSE, selection = 'single',
  options = list(columnDefs = list(list(className = 'dt-center', targets = c(3,4,5,6))), searchHighlight = TRUE, hover = FALSE, scrollX = TRUE, dom = 'lfrtip', pageLength = 10, lengthMenu = c(10, 25, 50, 100)), rownames = FALSE)
  
  
  output$bulk_digestDT <- DT::renderDataTable({
    withProgress(message = 'Computation in progress. This step might take a while depending on the the number of proteins and configuration parameters.', {
      incProgress(1/1)
      bulk_digestion()
    })
  }, escape = FALSE, selection = 'single',
  options = list(autoWidth = FALSE, columnDefs = list(list(className = 'dt-center', targets = c(1,2,3)), list(width = '25%', targets = "_all")), order = list(3, 'desc'), searchHighlight = TRUE, hover = FALSE, scrollX = TRUE, dom = 'lfrtip', pageLength = 10, lengthMenu = c(10, 25, 50, 100)), rownames = FALSE)
  
  
  output$protGroups_digestion <- DT::renderDataTable({
    withProgress(message = 'Computation in progress. This step might take a while depending on the the number of proteins and configuration parameters.', {
      incProgress(1/1)
      protGroups_digest()
      # dt_out <- protGroups_digest()
      # 
      # temp_dt <- dt_out[,c("id", "Protein IDs", "Sequence coverage [%]", "Theoritical max coverage per protein ID [%]")]
      # # temp_dt$`Theoritical max. coverage [%]` <- sapply(strsplit(as.character(temp_dt$`Theoritical max coverage per protein ID [%]`),";"), function(x, na.strings = NA) if(all(is.na(as.numeric(x)))) NA_real_ else max(as.numeric(x), na.rm = TRUE) )
      # temp_dt$`Theoritical max. coverage [%]` <- sapply(strsplit(as.character(temp_dt$`Theoritical max coverage per protein ID [%]`),";"), function(x) {
      #   x1 <- suppressWarnings(as.numeric(x))
      #   if(all(is.na(x1))) NA_real_ else max(x1, na.rm = TRUE)
      # })
      # 
      # # write.table(temp_dt[order(as.numeric(temp_dt$id)),c("id", "Protein IDs", "Sequence coverage [%]", "Theoritical max coverage per protein ID [%]", "Theoritical max. coverage [%]")], "C:\\Users\\grigo\\Desktop\\SICyLIA_v.2.0\\myShinyApp\\tbl.txt", quote=F, sep = "\t", row.names=F, col.names=T)
      # temp_dt[,c("id", "Protein IDs", "Sequence coverage [%]", "Theoritical max coverage per protein ID [%]", "Theoritical max. coverage [%]")]

      })
  }, escape = FALSE, selection = 'single',
  options = list(columnDefs = list(list(className = 'dt-center', targets = c(0,2,3))), order = list(0, 'asc'), rowCallback = JS(rowCallback), searchHighlight = TRUE, hover = FALSE, scrollX = TRUE, dom = 'lfrtip', pageLength = 10, lengthMenu = c(10, 25, 50, 100)), rownames = FALSE)
  
  
  output$peptidesSeqCoverage <- renderPlot({
    withProgress(message = 'Generating figure, please wait...', {
      incProgress(1/1)
      if (input$switchButtonDigest == TRUE) {
        dt.seq <- fastafile_reader()
      }
      else if (input$switchButtonDigest == FALSE) {
        dt.seq <- proteinList_reader()
      }
      dg <- digestedPep(dt.seq)
      grouped_len <- data.frame(table(nchar(dg$value)))
      grouped_len$Lbl <- paste0(grouped_len$Var1,"")
      grouped_len$SeqCovPerPept <- as.numeric(as.vector(grouped_len$Var1)) * as.numeric(as.vector(grouped_len$Freq))
      grouped_len$Perc <- round(grouped_len$SeqCovPerPept/sum(grouped_len$SeqCovPerPept)*100,2)
      topGroupedLen <- head(grouped_len, 35)
      topGroupedLen %>% 
        arrange(Var1) %>%
        mutate(Lbl=factor(Lbl, levels=Lbl)) %>%
        ggplot(aes(x=Lbl, y=Perc)) + 
        geom_bar( fill=c(rep("#d12741",(min_peptide_length$choice-1)), rep("#22478a",(35-min_peptide_length$choice)+1)), color="black", stat="identity", alpha=0.5 ) +
        geom_text(aes(label=Perc), vjust=-0.5, size=2.7) +
        # theme(legend.position="none") +
        ggtitle("Contribution of cleaved peptides in the total sequence coverage of the uploaded proteins (per peptide length)") +
        xlab("Length of peptides") +
        ylab("Percentage [%]") +
        theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, family="serif"),
              axis.text.y = element_text(size = 14, vjust = 0, family="serif"),
              plot.title = element_text(size = 15, family="serif"))
    })
  })
  
  # colorScale <- reactive({
  #   
  #   return("grey40")
  # })
  
  
  output$protGroups_theorVSobserv <- renderPlot({
    req(input$protGroupsDigestion)
    ext <- tools::file_ext(input$protGroupsDigestion$datapath)
    validate(need(ext == "txt", "Please upload a txt file"))
    dt1 <- read.table(input$protGroupsDigestion$datapath, header = T, sep = "\t", fill = TRUE, quote = "", check.names = F, stringsAsFactors = F)
    required_columns <- c("id","Protein IDs", "Sequence coverage [%]")
    if(all(required_columns %in% names(dt1)) == T) {
    withProgress(message = 'Generating figure, please wait...', {
      incProgress(1/1)
      dt2 <- protGroups_digest()
      dt3 <- merge(dt2, dt1[,c("id", "Majority protein IDs", "Intensity", "Gene names", "Sequence length" ,"Sequence lengths" )], by = "id")
      dt4 <- dt3[,c("id", "Majority protein IDs", "Intensity", "Gene names", "Sequence length", "Sequence coverage [%]", "Theoritical max. coverage [%]")]
      dt4 <- dt4[order(dt4$`Sequence length`),]
      dt4$Residuals <- dt4$`Sequence coverage [%]` - dt4$`Theoritical max. coverage [%]`
      p <- ggplot(dt4, aes(x=as.numeric(`Sequence coverage [%]`),  
                           y=as.numeric(`Theoritical max. coverage [%]`))
      ) + 
        geom_point(color="grey40", size=2, shape=16, na.rm = TRUE) +
        # geom_point(aes(color=rank(`Intensity`)), size=2, shape=16, na.rm = TRUE) +
        # geom_point(aes(color=col4), size=2, shape=16, na.rm = TRUE) +
        # geom_text(aes(label=lbls4), size=2.3, hjust=-0.2, vjust=-0.2) +
        geom_abline(intercept = 0, slope = 1) + 
        geom_vline(aes(xintercept = 5), colour="#BB0000", linetype = "dashed") +
        # ggtitle('Comparison between theoritical maximum and observed sequence coverage per protein group (336 in total)') +
        # ggtitle('Protein groups with either IAA- or NEM-labeled peptides (proteinGroups: 45 IAA only, 21 NEM only out of total 336  ~ 20%)') +
        # ggtitle('Protein groups with IAA- & NEM-labeled site(s) in common (125 out of 130/[336] proteinGroups ~ 96%/[37%])') +
        xlab('Observed (MQ) sequence coverage [%]') + 
        ylab('Theoritical maximum sequence coverage [%]') + 
        # labs(color='Colored by\nHPA annotation') +
        # labs(color='Color ranking\nbased on\nIntensity\n(ascending)') +
        # labs(color="Legend") +
        # scale_color_viridis_c() +
        # scale_colour_manual(values = c( vir_cols[2], "grey80")) +
        theme_grey()
      p + theme(legend.position = c(0.9, 0.25))
    })
  }
  })
  
  
  output$peptidesDist <- renderPlot({
    withProgress(message = 'Generating figure, please wait...', {
      incProgress(1/1)
    if (input$switchButtonDigest == TRUE) {
      dt.seq <- fastafile_reader()
    }
    else if (input$switchButtonDigest == FALSE) {
      dt.seq <- proteinList_reader()
    }
    dg <- digestedPep(dt.seq)
    grouped_len <- data.frame(table(nchar(dg$value)))
    grouped_len$Lbl <- paste0(grouped_len$Var1,"")
    grouped_len$Perc <- round(grouped_len$Freq/sum(grouped_len$Freq)*100,2)
    topGroupedLen <- head(grouped_len, 35)
    topGroupedLen %>% 
      arrange(Var1) %>%
      mutate(Lbl=factor(Lbl, levels=Lbl)) %>%
      ggplot(aes(x=Lbl, y=Freq)) + 
      geom_bar( fill=c(rep("#d12741",(min_peptide_length$choice-1)), rep("#22478a",(35-min_peptide_length$choice)+1)), color="black", stat="identity", alpha=0.5 ) +
      geom_text(aes(label=paste0(Perc,"%")), vjust=-0.5, size=2.7) +
      # theme(legend.position="none") +
      ggtitle("Frequency of cleaved peptides (including non-detectable)") +
      xlab("Length of peptides") +
      ylab("Number of peptides") +
      theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, family="serif"),
            axis.text.y = element_text(size = 14, vjust = 0, family="serif"),
            plot.title = element_text(size = 15, family="serif"))
    })
  })
  ##############################################
  #############################################

  output$table_1 <- DT::renderDataTable({
    withProgress(message = 'Loading data, please wait...', {
      incProgress(1/1)
      table_1()
    })
  }, options = list(scrollX = TRUE, dom = 'lfrtip', pageLength = 10, lengthMenu = c(10, 25, 50, 100)), rownames = FALSE)
  
  
  output$downloadTable1 <- downloadHandler(
    filename = function() {
      paste("table1_", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      if (input$downloadOptions1 == "std") {
        #print(colnames(table_1()))
        write.table(table_1(), file, sep="\t", row.names=FALSE, col.names = TRUE, quote = FALSE)
      }
      else {
        #print(colnames(table_1()))
        write.table(add_row(table_1(), 
                            Modifications = "#!{Type}N", 
                            Raw_file = "T", 
                            .before = 1), file, sep="\t", row.names=FALSE, col.names = TRUE, quote = FALSE)
      }
    }
  )
  
  
  output$table_2 <- DT::renderDataTable({
    withProgress(message = 'Loading data, please wait...', {
      incProgress(1/1)
      table_2()
    })
  }, options = list(scrollX = TRUE, dom = 'lfrtip', pageLength = 10, lengthMenu = c(10, 25, 50, 100)), rownames = FALSE)
  
  
  output$downloadTable2 <- downloadHandler(
    filename = function() {
      paste("table2_", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      if (input$downloadOptions2 == "std") {
        write.table(table_2(), file, sep="\t", row.names=FALSE, col.names = TRUE, quote = FALSE)
      }
      else {
        write.table(add_row(table_2(), 
                            Modifications = "#!{Type}N", 
                            Raw_file = "T", 
                            .before = 1), file, sep="\t", row.names=FALSE, col.names = TRUE, quote = FALSE)
      }
    }
  )
  
  
  output$table_3 <- DT::renderDataTable({
    withProgress(message = 'Loading data, please wait...', {
      incProgress(1/1)
      table_3()
    })
  }, options = list(scrollX = TRUE, dom = 'lfrtip', pageLength = 10, lengthMenu = c(10, 25, 50, 100)), rownames = FALSE)
  
  
  output$downloadTable3 <- downloadHandler(
    filename = function() {
      paste("table3_", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      if (input$downloadOptions3 == "std") {
        write.table(table_3(), file, sep="\t", row.names=FALSE, col.names = TRUE, quote = FALSE)
      }
      else {
        write.table(add_row(table_3(), 
                            id = "#!{Type}C", 
                            Peptide_ID = "C", 
                            .before = 1), file, sep="\t", row.names=FALSE, col.names = TRUE, quote = FALSE)
      }
    }
  )
  

  observe({
    updateSelectInput(session = session, inputId = "rep_intensities_x", choices = names(table_3()[,grepl("Normalised_intensity_corrected_" , names(table_3()))]) )
  })

  observe({
    updateSelectInput(session = session, inputId = "rep_intensities_y", choices = names(table_3()[,grepl("Normalised_intensity_corrected_" , names(table_3()))]) )
  })


  output$plot <- renderPlot({
      if (!is.null(dt_main_dt())){
        withProgress(message = 'Generating figure, please wait...', {
          incProgress(1/1)
          plot(table_3()[,c(input$rep_intensities_x, input$rep_intensities_y)])
          abline(fit <- lm(table_3()[,c(input$rep_intensities_x)] ~ table_3()[,c(input$rep_intensities_y)], data=table_3()), col='red')
          legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=3)))
        })
      }
      else
        NULL
  })
  
  
  output$plot2 <- renderPlot({
    if (!is.null(dt_main_dt())){
      withProgress(message = 'Generating figure, please wait...', {
        incProgress(1/1)
        len <- length(grep(x = colnames(table_3()), pattern = "^Normalised_intensity_corrected_"))
        par(mfrow=c(len,len))
        par(mar=c(0,0,0,0)+0.1)
        #par(pty = "m")
        for (i in 0:(len-1)){
          for (j in 0:(len-1)){
            #print(paste0(i,"-",j))
             if (i==j) {
              plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab=NULL, ylab=NULL, main=NULL, yaxt='n',xaxt='n', ann=FALSE)
              legend("center", bty="n", cex=1.5, legend="1")
             }
            else if (i>j){
              plot(table_3()[,c(paste0("Normalised_intensity_corrected_",i), paste0("Normalised_intensity_corrected_",j))], xlab=NULL, ylab=NULL, main=NULL, yaxt='n',xaxt='n', ann=FALSE)
              abline(fit <- lm(table_3()[,c(paste0("Normalised_intensity_corrected_",i))] ~ table_3()[,c(paste0("Normalised_intensity_corrected_",j))], data=table_3()), col='red')
              legend("topleft", bty="n", cex=1.2, legend=format(summary(fit)$adj.r.squared, digits=3))
            }
            else {
              plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab=NULL, ylab=NULL, main=NULL, yaxt='n',xaxt='n', ann=FALSE)
              fit <- lm(table_3()[,c(paste0("Normalised_intensity_corrected_",i))] ~ table_3()[,c(paste0("Normalised_intensity_corrected_",j))], data=table_3())
              legend("center", bty="n", cex=1.5, legend=format(summary(fit)$adj.r.squared, digits=3))
              
            }
          }
        }
      })
    }
    else
      NULL
  })
  
 #observer to trigger shinyAlert when the question-mark (in the ui.R) is pressed  
  observeEvent(input$alertMessage, {
    alertData <- input$alertMessage
    req(alertData)
    shinyalert(
      title = alertData$title,
      text = alertData$message,
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = TRUE,
      type = alertData$type,
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#95a5a6",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
  })
  
  observe({
      #shinyjs::disable(selector = "[type=radio][value=11]")
      #shinyjs::runjs("$('[type=radio][value=11]').parent().addClass('disabled').css('opacity', 0.5)")
      shinyjs::disable(selector = "[type=radio][value=16]")
      shinyjs::runjs("$('[type=radio][value=16]').parent().addClass('disabled').css('opacity', 0.5)")
  })
 
  
  output$downloadDummyRec <- downloadHandler(
    filename <- function() {
      paste("output", ".zip", sep="")
    },
    content <- function(file) {
      file.copy("www/out.zip", file)
    },
    contentType = "application/zip"
  )
  
  

}
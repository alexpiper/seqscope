library(shiny)
library(ggplot2)
library(reshape2)
library(vegan)
library(dplyr)
library(phyloseq)
library(broom)
library(plotly)
library(tibble)
library(scales)
library(heatmaply)
library(markdown)
library(ranacapa)
library(speedyseq)

options(digits = 5, shiny.maxRequestSize = 10 * 1024 ^ 2)

server <- function(input, output)({
  
  # Setup and RenderUIs ---------------
  # RenderUI for which_variable_r, gets used in Panels 1, 3, 4, 5, 6
  output$which_variable_r <- renderUI({
    selectInput("var", "Select the variable", choices = heads())
  })
  output$which_variable_alphaDiv <- renderUI({
    selectInput("var_alpha", "Select the variable", choices = heads_alpha_Anova())
  })
  # RenderUIs for Panel 1
  output$seqtabSelect <- renderUI({
    req(input$mode)
    if (input$mode == "Real") {
      fileInput("in_seqtab", "Please select your sequence table.
                Note: this should be saved either as *.RDS or *.csv",
                accept = c(".rds", ".csv"))
    }
  })
  output$taxSelect <- renderUI({
    req(input$mode)
    if (input$mode == "Real") {
      fileInput("in_taxtab", "Please select your taxonomy table.
                Note: this should be saved either as *.RDS or *.csv",
                accept = c(".rds", ".csv"))
    }
  })
  output$metaSelect <- renderUI({
    req(input$mode)
    if (input$mode == "Real") {
      fileInput("in_metadata", "Please select the metadata file.
                Note: this should be saved as *.csv or *.txt",
                accept = c(".csv", ".txt"))
    }
  })

  
  # RenderUI for which_taxon_level, used for barplot and heatmap in Panels 7,8
  output$which_taxon_level <- renderUI({
    radioButtons("taxon_level",
                 "Pick the taxonomic level for making the plot",
                 choices = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
    
  })
  
  
  choices <- reactive({
    c("Select All", rownames(for_hm()))
  })
  output$select_species_heat <- renderUI({
    selectInput('which_taxa_heat', 'Select taxa to visualize', choices(),
                multiple=TRUE, selectize=FALSE, selected = "Select All")
  })

  # Read in data files, validate and make the physeq object -----
  
  seqtab <- reactive({
    if (input$mode == "Real") {
      if (grepl(input$in_seqtab$datapath, pattern = ".txt") |
          grepl(input$in_seqtab$datapath, pattern = ".tsv")) {
        read.table(input$in_seqtab$datapath, header = 1,
                   sep = "\t", stringsAsFactors = F,
                   quote = "", comment.char = "")
      }else if (grepl(input$in_seqtab$datapath, pattern = ".rds")) {
        readRDS(input$in_seqtab$datapath)
      }
    } else {
      readRDS("data/demo_seqtab.rds")
    }
  })
  
  taxtab <- reactive({
    if (input$mode == "Real") {
      if (grepl(input$in_taxtab$datapath, pattern = ".txt") |
          grepl(input$in_taxtab$datapath, pattern = ".tsv")) {
        read.table(input$in_taxtab$datapath, header = 1,
                   sep = "\t", stringsAsFactors = F,
                   quote = "", comment.char = "")
      }else if (grepl(input$in_taxtab$datapath, pattern = ".rds")) {
        readRDS(input$in_taxtab$datapath)
      }
    } else {
      readRDS("data/demo_taxtab.rds")
    }
  })
  
  samdf <- reactive({
    if (input$mode == "Real") {
      if (grepl(readLines(input$in_metadata$datapath, n = 1), pattern = "^#")) {
        phyloseq::import_qiime_sample_data(input$in_metadata$datapath) %>%
          as.matrix() %>%
          as.data.frame()
      } else {
        read.table(input$in_metadata$datapath,
                   header = 1, sep = "\t", stringsAsFactors = F,
                   quote = "", comment.char = "")
      }
    } else {
      read.csv("data/demo_samdf.csv") %>%
        dplyr::filter(!duplicated(SampleID)) %>%
        magrittr::set_rownames(.$SampleID) 
    }
  })
  
  #Validate input files - need to adapt validation for our particular onputs
#  
# output$fileStatus <- eventReactive(input$go, {
#   if (is.null(validate_input_files(taxtab(), samdf()))) {
#     paste("Data sucessfully loaded")
#   } else {
#     validate_input_files(taxtab(), samdf())
#   }
# })
  
  
  # Make physeq object 
  physeq <- eventReactive(input$go, {
    phyloseq(tax_table(taxtab()), sample_data(samdf()),
                 otu_table(seqtab(), taxa_are_rows = FALSE)) 

    
    # phy_tree(phytree(),)
  ##if(nrow(seqtab.nochim) > nrow(sample_data(ps))){warning("Warning: All samples not included in phyloseq object, check sample names match the sample metadata")}
  })
  
  heads <- reactive({
    base::colnames(samdf())
  })
  
  heads_numeric <- reactive({
    samdf() %>%
      dplyr::select_if(is.numeric) %>%
      base::colnames()
  })
  
  heads_alpha_Anova <- reactive({
    num_factors <- sapply(samdf(), function(col) length(unique(col)))
    heads()[num_factors > 2]
  })
  
  # Panel 2:  Print taxon table ---------
  
  output$print_taxon_table <- DT::renderDataTable({
    table <- taxtab() %>% 
      as.data.frame() %>%
      tibble::rownames_to_column(var="seq") %>%
      tidyr::unite(sum.taxonomy, 2:ncol(.), sep=";") %>%
      mutate(sum.taxonomy = stringr::str_replace(sum.taxonomy, pattern = "NA", replacement = "")) %>%
      right_join(seqtab() %>%
                   as.data.frame() %>%
                   t()  %>% 
                   as.data.frame() %>%
                   tibble::rownames_to_column(var="seq"),
                 by=c("seq")) %>%
      select(sum.taxonomy, everything())
    
    DT::datatable(table, options = list(scrollX = TRUE))
  })
  output$print_metadata_table <- DT::renderDataTable({
    DT::datatable(samdf(), options = list(scrollX = TRUE))
  }, options = list(pageLength = 5))
  
  
  # Panel 3: Filtering ---------- Need to change rarefied to relative abundance
  # Check if all samples have a non-NA value for the selected variable to plot by
  # If a sample has an NA for the selected variable, get rid of it from the
  # sample data and from the metadata and from the taxon table (the subset function does both)
  data_subset <- reactive({
    p2 <- physeq()
    sample_data(p2) <- physeq() %>%
      sample_data %>%
      subset(., !is.na(get(input$var)))
    p2
  })
 
 # Convert subsetted dataset to relative abundance
 data_subset_ra <- reactive({
   if (input$ra_method == "RA") {
     p2 <- data_subset() 
     newphyseq <- as(otu_table(p2), "matrix")[which(rowSums(as(otu_table(p2), "matrix")) > 0),] # remove empty samples
     newphyseq <- apply(newphyseq, 1, function(x) x/sum(x, na.rm=FALSE))
     otu_table(p2) <- otu_table(newphyseq, taxa_are_rows = TRUE)
     cat(file=stderr(), "Data transformed \n")
     p2
   } else {
     data_subset()
   }
 })
  
  
  # Panel 5: Taxonomy-by-site interactive barplot -------
  output$tax_bar <- renderPlotly({
    
    withProgress(message = 'Rendering taxonomy barplot', value = 0, {
      incProgress(0.5)
        physeqGlommed = speedyseq::tax_glom(data_subset_ra(), input$taxon_level, NArm = FALSE)
        cat(file=stderr(), "Glommed \n")
      plot_bar(physeqGlommed, fill = input$taxon_level) + theme_ranacapa() +
        theme(axis.text.x = element_text(angle = 45)) +
        theme(axis.title = element_blank())
      gp <- ggplotly() %>%
        layout(yaxis = list(title = "Abundance", titlefont = list(size = 16)),
               xaxis = list(title = "Sample", titlefont = list(size = 16)),
               margin = list(l = 70, b = 100))
      gp
    })
  })
  
  
  ## Panel 8: Heatmap of taxonomy by site ---------
  for_hm <- reactive({
    
    if (input$ra_method == "RA") {
      tt <-  otu_table(data_subset_ra())
      
    } else {
      tt <- t(otu_table(data_subset())) 
      
    }
    
    for_hm <- tt %>% 
    #for_hm <- t(otu_table(p2)) %>% 
      as.data.frame() %>%
      tibble::rownames_to_column(var="OTU") %>%
      left_join(speedyseq::psmelt(data_subset_ra()) %>% 
      #left_join(speedyseq::psmelt(p2) %>%
                  select("OTU", "loci", "Phylum", "Class", "Order", "Family", "Genus","Species") %>%
                  mutate_all(as.character),
                by="OTU") %>%
      unique()
    
    for_hm <- for_hm %>%
      mutate(loci = case_when(is.na(loci) ~ "unknown", TRUE ~ loci)) %>%
      mutate(Phylum = case_when(is.na(Phylum) ~ "unknown", TRUE ~ Phylum)) %>%
      mutate(Class = case_when(is.na(Class) ~ "unknown", TRUE ~ Class)) %>%
      mutate(Order = case_when(is.na(Order) ~ "unknown", TRUE ~ Order)) %>%
      mutate(Family = case_when(is.na(Family) ~ "unknown", TRUE ~ Family)) %>%
      mutate(Genus = case_when(is.na(Genus) ~ "unknown", TRUE ~ Genus)) %>%
      mutate(Species = case_when(is.na(Species) ~ "unknown", TRUE ~ Species)) 
    
    for_hm <- for_hm %>%
      group_by(get(input$taxon_level)) %>%
      #group_by(Species) %>%
      summarize_if(is.numeric, sum) %>%
      data.frame %>%
      #column_to_rownames("Species") #%>%
      column_to_rownames("get.input.taxon_level.")
    for_hm <- for_hm[which(rowSums(for_hm) > 0),]
    for_hm[for_hm == 0] <- NA
    for_hm
  })
  output$tax_heat <- renderPlotly({
    
    withProgress(message = 'Rendering taxonomy heatmap', value = 0, {
      incProgress(0.5)
      if("Select All" %in% input$which_taxa_heat){
        selected_taxa <- rownames(for_hm())
      } else{
        selected_taxa <- input$which_taxa_heat
      }
      for_hm <- for_hm()[selected_taxa,]
      heatmaply(for_hm, Rowv = F, Colv = F, hide_colorbar = F,
                grid_gap = 1, na.value = "white", key.title = "Number of \nSequences in \nSample")
    })
    
  })
  
  table_for_download <- reactive({
    taxcol <- reshape2::colsplit(taxtab()$sum.taxonomy, ";", paste0("V", 1:6))%>%
      mutate(V1 = paste0("p__", V1),
             V2 = paste0("c__", V2),
             V3 = paste0("o__", V3),
             V4 = paste0("f__", V4),
             V5 = paste0("g__", V5),
             V6 = paste0("s__", V6))  %>%
      mutate(taxonomy =  paste(V1, V2, V3, V4, V5, V6, sep = ";")) %>%
      select(-c(V1, V2, V3, V4, V5, V6))
    cbind(taxtab() %>% rename(taxID = sum.taxonomy), taxcol)
    
  })
  output$downloadTableForBiom <- downloadHandler(
    
    
    filename = function() {
      paste("taxonomy-for-biom.txt", sep = "")
    },
    content = function(file) {
      write.csv(table_for_download(), file, row.names = FALSE, quote = F)
    }
  )
  
  output$downloadPhyloseqObject <- downloadHandler(
    
    
    filename = function() {
      paste("phyloseq-object.Rds", sep = "")
    },
    content = function(file) {
      saveRDS(data_subset(), file)
    }
  )
  })

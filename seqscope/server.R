library(shiny)
library(ggplot2)
library(vegan)
library(dplyr)
library(phyloseq)
library(plotly)
library(tibble)
library(heatmaply)
library(markdown)
library(speedyseq)
library(RColorBrewer)

options(digits = 5, shiny.maxRequestSize = 10 * 1024 ^ 2)

server <- function(input, output)({
  

# Renderuis ---------------------------------------------------------------

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
  output$phyloSelect <- renderUI({
    req(input$mode)
    if (input$mode == "Real") {
      fileInput("in_phylo", "Please select the phylogeny file.
                Note: this should be saved as *.rds or *.txt",
                accept = c(".rds", ".txt"))
    }
  })

  
  # RenderUI for which_taxon_level, used for barplot and heatmap in Panels 7,8
  output$which_taxon_level <- renderUI({
    radioButtons("taxon_level",
                 "Select taxonomic level",
                 choices = c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"),
                 selected="Species")
    
  })
  
  
  choices <- reactive({
    c("Select All", rownames(for_hm()))
  })
  output$select_species_heat <- renderUI({
    selectInput('which_taxa_heat', 'Select taxa to visualize', choices(),
                multiple=TRUE, selectize=FALSE, selected = "Select All")
  })


# Read in data ------------------------------------------------------------


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
  
  phytree <- reactive({
    if (input$mode == "Real") {
      if (grepl(input$in_taxtab$datapath, pattern = ".txt")) {
        #read_phylo_txt
      } else {
        readRDS(input$in_phylo$datapath)
      }
    } else {
      readRDS("data/demo_phytree.rds")$tree
    }
  })
  

# Validate input files ----------------------------------------------------

# need to adapt validation for our particular inputs
#  
# output$fileStatus <- eventReactive(input$go, {
#   if (is.null(validate_input_files(taxtab(), samdf()))) {
#     paste("Data sucessfully loaded")
#   } else {
#     validate_input_files(taxtab(), samdf())
#   }
# })
  
  

# Create phyloseq object --------------------------------------------------

  physeq <- eventReactive(input$go, {
    phyloseq(tax_table(taxtab()), sample_data(samdf()),
                 otu_table(seqtab(), taxa_are_rows = FALSE), phy_tree(phytree())) 
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
  

# Print taxon table -------------------------------------------------------

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
  
  

# Filter and transform ----------------------------------------------------

  # Check if all samples have a non-NA value for the selected variable to plot by
  # If a sample has an NA for the selected variable, get rid of it from the
  # sample data and from the metadata and from the taxon table (the subset function does both)
  data_subset <- reactive({
    p2 <- physeq()
    sample_data(p2) <- physeq() %>%
      sample_data 
    taxa_names(p2) <- paste0("SV", seq(ntaxa(p2)),"-",tax_table(p2)[,7])
    p2
  })
 
 # Convert subsetted dataset to relative abundance
 data_subset_ra <- reactive({
   if (input$ra_method == "Relative Abundance") {
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
  
  

# Taxonomy barchart -------------------------------------------------------

   output$tax_bar <- renderPlotly({
    
    withProgress(message = 'Rendering taxonomy barplot', value = 0, {
      incProgress(0.5)
      
      if(!input$taxon_level =="OTU"){
        physeqGlommed = speedyseq::tax_glom(data_subset_ra(), input$taxon_level, NArm = FALSE)
        speedyseq::plot_bar(physeqGlommed, fill = input$taxon_level) + theme_bw() +
          theme(axis.text.x = element_text(angle = 45)) +
          theme(axis.title = element_blank()) 
      } else {
        physeqGlommed = data_subset_ra()
        speedyseq::plot_bar(physeqGlommed, fill="Species") + theme_bw() +
          theme(axis.text.x = element_text(angle = 45)) +
          theme(axis.title = element_blank()) 
        }
      
      
      gp <- ggplotly() %>%
        layout(yaxis = list(title = "Abundance", titlefont = list(size = 16)),
               xaxis = list(title = "Sample", titlefont = list(size = 16)),
               margin = list(l = 70, b = 100))
      gp
    })
  })
  


# Heatmap -----------------------------------------------------------------

  for_hm <- reactive({
    
    if (input$ra_method == "Relative Abundance") {
      tt <-  otu_table(data_subset_ra())
      
    } else {
      tt <- t(otu_table(data_subset())) 
      
    }
    
    for_hm <- tt %>% 
      as.data.frame() %>%
      tibble::rownames_to_column(var="OTU") %>%
      left_join(speedyseq::psmelt(data_subset_ra()) %>% 
                  select("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species") %>%
                  mutate_all(as.character),
                by="OTU") %>%
      unique()
    
    for_hm <- for_hm %>%
      mutate(Kingdom = case_when(is.na(Kingdom) ~ "unknown", TRUE ~ Kingdom)) %>%
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
      
      if (input$ra_method == "Relative Abundance") {
        heatmaply(for_hm, Rowv = F, Colv = F, hide_colorbar = F,
                  grid_gap = 1, na.value = "white", key.title = "Relative Abundance \n of Taxa in Sample")
      } else {
        heatmaply(for_hm, Rowv = F, Colv = F, hide_colorbar = F,
                  grid_gap = 1, na.value = "white", key.title = "Number of \nSequences in \nSample")
      }
    
    })
    
  })
  
  

# Phylogeny ---------------------------------------------------------------

  output$tax_phylo <- renderPlotly({
    
    withProgress(message = 'Rendering taxonomy barplot', value = 0, {
      incProgress(0.5)
      if(!input$taxon_level =="OTU"){
      physeqGlommed = speedyseq::tax_glom(data_subset_ra(), input$taxon_level, NArm = FALSE)
      } else (physeqGlommed = data_subset_ra())
      
      plot_tree(physeqGlommed, color="SampleID", label.tips="taxa_names", ladderize="left") + theme_bw() # need to rename taxa 
      gp <- ggplotly() #%>%
      #layout(yaxis = list(title = "Abundance", titlefont = list(size = 16)),
      #       xaxis = list(title = "Sample", titlefont = list(size = 16)),
      #       margin = list(l = 70, b = 100))
      gp
    })
  })
  
  

# Output table ------------------------------------------------------------

  
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

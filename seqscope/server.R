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
  

# Render uis ---------------------------------------------------------------

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
  
  


# UI Filter ------------------------------------------------------------------
 # UI subset_taxa expression cascade
 # filter_subset_taxa_expr
 output$filter_uix_subset_taxa_ranks <- renderUI({
   rankNames = list("NULL"="NULL")
   rankNames <- c(rankNames, as.list(rank_names(physeq(), errorIfNULL=FALSE)))
   rankNames <- c(rankNames, list(OTU="OTU"))
   return(
     selectInput("filter_rank", "Taxonomic Ranks", rankNames, "NULL", multiple = FALSE)
   )
 })
 output$filter_uix_subset_taxa_select <- renderUI({
   rank_list = list("NULL" = "NULL")
   if(!is.null(av(input$filter_rank))){
     # If a filter rank is specified, use it, and provide the multi-select widget for these rank classes
     if(input$filter_rank == "OTU"){
       rank_list <- c(rank_list, as.list(taxa_names(physeq())))
     } else {
       rank_list <- c(rank_list, as.list(get_taxa_unique(physeq(), input$filter_rank)))
     }
   }
   return(
     selectInput(inputId = "filter_rank_selection", label = "Select Taxa",
                 choices = rank_list, selected = "NULL", multiple = TRUE)
   )
 })
 

# UI Subset samples ----------------------------------------------------------
 # UI subset_samples expression cascade
 # filter_subset_samp_expr
 output$filter_uix_subset_sample_vars <- renderUI({
   sampVars = list("NULL"="NULL")
   sampVars <- c(sampVars, as.list(sample_variables(physeq(), errorIfNULL=FALSE)))
   sampVars <- c(sampVars, list(Sample="Sample"))
   return(
     selectInput("filter_samvars", "Sample Variables", sampVars, "NULL", multiple = FALSE)
   )
 })
 output$filter_uix_subset_sample_select <- renderUI({
   varLevels = list("NULL"="NULL")
   if(!is.null(av(input$filter_samvars))){
     if(input$filter_samvars == "Sample"){
       varLevels <- c(varLevels, as.list(sample_names(physeq())))
     } else {
       if(!is.null(sample_variables(physeq(), FALSE))){
         varvec = get_variable(physeq(), input$filter_samvars)
         if(plyr::is.discrete(varvec)){
           varLevels <- c(varLevels, as.list(unique(as(varvec, "character"))))
         }
       } 
     }
   }
   return(
     selectInput(inputId = "filter_samvars_selection", label = "Variable Classes",
                 choices = varLevels, selected = "NULL", multiple = TRUE)
   )  
 })
 

# Filtering ---------------------------------------------------------------
 # The main reactive data object. Returns a phyloseq-class instance.
 # This is considered the "filtered" data, used by all downstream panels,
 # And generally the input to any transformation options as well
 
 data_subset <- reactive({
   ps0 <- physeq() 
   if(input$actionb_filter == 0){
     # Don't execute filter if filter-button has never been clicked.
     if(inherits(ps0, "phyloseq")){
       return(ps0)
     } else {
       return(NULL)
     }
   }
   # Isolate all filter code so that button click is required for update
   isolate({
     if(inherits(ps0, "phyloseq")){
       message("ISOLATED")
       # Cascading selection filters
       if( !is.null(av(input$filter_rank_selection)) ){
         keepTaxa = NULL
         if(!is.null(tax_table(ps0, FALSE))){
           if(input$filter_rank == "OTU"){
             # OTU IDs directly
             keepTaxa = input$filter_rank_selection
           } else {
             TT = as(tax_table(ps0), "matrix")
             keepTaxa = TT[, input$filter_rank] %in% input$filter_rank_selection 
           }
           if(length(keepTaxa) > 1){
             ps0 <- prune_taxa(keepTaxa, ps0)
           } else {
             warning("Bad subset_taxa specification. ntaxa(ps0) one or fewer OTUs")
           }
         }
       }
       if( !is.null(av(input$filter_samvars_selection)) ){
         keepSamples = NULL
         if(!is.null(sample_data(ps0, FALSE))){
           if(input$filter_samvars == "Sample"){
             # Samples IDs directly
             keepSamples = input$filter_samvars_selection
           } else {
             varvec = as(get_variable(ps0, input$filter_samvars), "character")
             keepSamples = varvec %in% input$filter_samvars_selection 
           }
           if(length(keepSamples) > 1){
             ps0 <- prune_samples(keepSamples, ps0)
           } else {
             warning("Bad subset_taxa specification. ntaxa(ps0) one or fewer OTUs")
           }
         }
       }
       if( input$filter_taxa_sums_threshold > 0 ){
         # OTU sums filter
         ps0 <- prune_taxa({taxa_sums(ps0) > input$filter_taxa_sums_threshold}, ps0)
       }
       if( input$filter_sample_sums_threshold > 0 ){
         # Sample sums filtering
         ps0 <- prune_samples({sample_sums(ps0) > input$filter_sample_sums_threshold}, ps0)
       }
       if(inherits(input$filter_kOverA_sample_threshold, "numeric")){
         if(input$filter_kOverA_sample_threshold > 1){
           # kOverA OTU Filtering
           flist = genefilter::filterfun(
             genefilter::kOverA(input$filter_kOverA_sample_threshold,
                                input$filter_kOverA_count_threshold, na.rm=TRUE)
           )
           koatry = try(ps0 <- filter_taxa(ps0, flist, prune=TRUE), silent = TRUE)
           if(inherits(koatry, "try-error")){
             warning("kOverA parameters resulted in an error, kOverA filtering skipped.")
           }
         }
       }
       return(ps0)
     } else {
       return(NULL)
     }
   })
 })

  
# kOverA `k` Filter UI
maxSamples <- reactive({
  # Create logical indicated the samples to keep, or dummy logical if nonsense input
  if(inherits(physeq(), "phyloseq")){
    return(nsamples(physeq()))
  } else {
    # Dummy response.
    return(NULL)
  }
})
output$filter_ui_kOverA_k <- renderUI({
  numericInput("filter_kOverA_sample_threshold", "Across K samples",
                  min=0, max=maxSamples(), value=kovera_k, step=1)
})


# Filtering Histograms ----------------------------------------------------

sums_hist <- function(thesums=NULL, xlab="", ylab="", scale=TRUE){
  if(is.null(thesums)){
    p = qplot(0)
  } else {
    p = ggplot(data.frame(sums=thesums), aes(x=sums))
    p = p + geom_histogram()
    p = p + xlab(xlab) + ylab(ylab) 
    if(scale==TRUE){
      p = p + scale_x_log10(labels = scales::comma)
    }
  }
  return(p)
}

read_length_hist <- function(ps=NULL, xlab = "Sequence length", ylab = "Number of Reads (Counts)" ){
  if(is.null(ps)){
    p = qplot(0)
  } else {
  sums <- as.data.frame(taxa_sums(ps)) %>%
    magrittr::set_colnames("sums") %>%
    tibble::rownames_to_column("length") %>%
    dplyr::mutate(length = nchar(length)) %>%
    dplyr::group_by(length) %>%
    dplyr::summarise(sums = sum(sums))
  p = ggplot(sums, aes(x=length, y=sums))
  p = p + geom_bar(stat="identity")
  p = p + xlab(xlab) + ylab(ylab) 
  #p = p + scale_y_log2(labels = scales::comma)
  }
  return(p)
}

lib_size_hist <- reactive({
  xlab = "Number of Reads (Counts)"
  ylab = "Number of Libraries"
  return(sums_hist(sample_sums(physeq()), xlab, ylab))
})
otu_sum_hist <- reactive({
  xlab = "Number of Reads (Counts)"
  ylab = "Number of OTUs"
  return(sums_hist(taxa_sums(physeq()), xlab, ylab))    
})

length_reads_hist <- reactive({
  xlab = "Number of Reads (Counts)"
  ylab = "Sequence length"
  return(read_length_hist(physeq(), xlab, ylab))    
})
length_OTU_hist <- reactive({
  xlab = "Sequence length" 
  ylab = "Number of OTU's"
  return(sums_hist(nchar(colnames(get_taxa(physeq()))), xlab, ylab, scale=FALSE)) 
})

output$sample_variables <- renderText({return(
  paste0(sample_variables(physeq(), errorIfNULL=FALSE), collapse=", ")
)})
output$rank_names <- renderText({return(
  paste0(rank_names(physeq(), errorIfNULL=FALSE), collapse=", ")
)})

# Plot filter plots
output$filter_summary_plot <- renderPlot({
  plib0 = lib_size_hist() + ggtitle("Original Data")
  potu0 = otu_sum_hist() + ggtitle("Original Data")
  if(inherits(data_subset(), "phyloseq")){
    potu1 = sums_hist(taxa_sums(data_subset()), xlab = "Number of Reads (Counts)", ylab = "Number of OTUs" ) +
      ggtitle("Filtered Data")
    
    plib1 = sums_hist(sample_sums(data_subset()), xlab = "Number of Reads (Counts)", ylab = "Number of Libraries" ) +
      ggtitle("Filtered Data")

  } else {
    potu1 = plib1 = fail_gen()
  }
  gridExtra::grid.arrange(plib0, plib1, potu0,  potu1,  ncol=2) #, main="Histograms: Before and After Filtering")
  
})
output$length_summary_plot <- renderPlot({
  plreads0 = read_length_hist(ps = physeq()) + ggtitle("Original Data")
  plotu0 = length_OTU_hist() + ggtitle("Original Data") 
  if(inherits(data_subset(), "phyloseq")){
    plreads1 = read_length_hist(data_subset()) +
      ggtitle("Filtered Data")
    
    plotu1 = sums_hist(nchar(colnames(get_taxa(physeq()))), xlab = "Sequence length", ylab = "Number of OTU's", scale=FALSE) +
      ggtitle("Filtered Data")
    
  } else {
    plreads1 = plotu1 = fail_gen()
  }
  gridExtra::grid.arrange(plreads0, plreads1, plotu0, plotu1,  ncol=2) #, main="Histograms: Before and After Filtering")
  
})

# transform ----------------------------------------------------

# Check if all samples have a non-NA value for the selected variable to plot by
# If a sample has an NA for the selected variable, get rid of it from the
# sample data and from the metadata and from the taxon table (the subset function does both)
data_transformed <- reactive({
  p2 <- data_subset()
  sample_data(p2) <- data_subset() %>%
    sample_data 
  taxa_names(p2) <- paste0("SV", seq(ntaxa(p2)),"-",tax_table(p2)[,7])
  p2
})

# Convert subsetted dataset to relative abundance
data_transformed_ra <- reactive({
  if (input$ra_method == "Relative Abundance") {
    p2 <- data_transformed() 
    newphyseq <- as(otu_table(p2), "matrix")[which(rowSums(as(otu_table(p2), "matrix")) > 0),] # remove empty samples
    newphyseq <- apply(newphyseq, 1, function(x) x/sum(x, na.rm=FALSE))
    otu_table(p2) <- otu_table(newphyseq, taxa_are_rows = TRUE)
    cat(file=stderr(), "Data transformed \n")
    p2
  } else {
    data_transformed()
  }
})

 

# Taxonomy barchart -------------------------------------------------------

   output$tax_bar <- renderPlotly({
    
    withProgress(message = 'Rendering taxonomy barplot', value = 0, {
      incProgress(0.5)
      
      if(!input$taxon_level =="OTU"){
        physeqGlommed = speedyseq::tax_glom(data_transformed_ra(), input$taxon_level, NArm = FALSE)
        speedyseq::plot_bar(physeqGlommed, fill = input$taxon_level) + theme_bw() +
          theme(axis.text.x = element_text(angle = 45)) +
          theme(axis.title = element_blank()) 
      } else {
        physeqGlommed = data_transformed_ra()
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
      tt <-  otu_table(data_transformed_ra())
      
    } else {
      tt <- t(otu_table(data_transformed())) 
      
    }
    
    for_hm <- tt %>% 
      as.data.frame() %>%
      tibble::rownames_to_column(var="OTU") %>%
      left_join(speedyseq::psmelt(data_transformed_ra()) %>% 
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
      physeqGlommed = speedyseq::tax_glom(data_transformed_ra(), input$taxon_level, NArm = FALSE)
      } else (physeqGlommed = data_transformed_ra())
      
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
      saveRDS(data_transformed(), file)
    }
  )
  })

#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

shinyUI(pageWithSidebar(
                 
           headerPanel(
               list(HTML('<img src="imappests.jpg" width="240" height="70"/>'), "SeqScope"),
               windowTitle="SeqScope"

           
            # img(src='seqscope/agvic.jpg', align = "left"),"SeqScope"
            ),
           sidebarPanel(
             
             ## conditionalPanel() functions for selected tab
             conditionalPanel(condition="input.tabselected==1"),
             # For Panel 1, have an input option
             conditionalPanel(condition = "input.tabselected == 2",
                              h3("Run with demo or real dataset?"),
                              
                              radioButtons("mode", label = "",
                                           choices = c("Demo", "Real"), selected = "Demo"),
                              uiOutput("seqtabSelect"),
                              uiOutput("taxSelect"),
                              uiOutput("metaSelect"),
                              h3("Press the button below to process data"),
                              actionButton("go", "Process"),
                              textOutput("fileStatus")
             ),    
             # On panels 5 and 6 (barplot and heatmap), ask which taxonomic level they want to visualize to
             conditionalPanel(condition = "input.tabselected == 5 | input.tabselected == 6",
                              uiOutput("which_taxon_level")),
             conditionalPanel(condition = "input.tabselected == 6",
                              uiOutput("select_species_heat")),
             
             # On panels 5 and 6 (barplot and heatmap) ask user what variables should be visualised
             conditionalPanel(condition = "input.tabselected == 5 |
                              input.tabselected == 6",
                              uiOutput("which_variable_r")),
             
             # On panel 3, 5 and 6, ask if counts or RA shoudl be displayed
             conditionalPanel(condition = "input.tabselected == 3 | input.tabselected == 5 |
                              input.tabselected == 6",
                              radioButtons("ra_method", "Display samples as counts or relative abundances",
                                           choices = c("counts", "RA")),
                              uiOutput("rare_depth"))
           ),

      mainPanel(
        tabsetPanel(
          tabPanel("Welcome",value=1,
                     includeMarkdown("docs/welcome.md")),
          tabPanel("Data", value=2,
                   h2("Verify that the files below look as expected, and click on process to load all data"),
                   h4("Input taxonomy file"),
                   DT::dataTableOutput("print_taxon_table"),
                   h4("Input metadata file"),
                   DT::dataTableOutput("print_metadata_table"),
                   h4("Input QC file"),
                   DT::dataTableOutput("print_qc_table")
                   ),
          #tabPanel("Filter", value = 3),
          
          #   #tabPanel("Filter", value = 4),
              tabPanel("Barchart", value = 5,
                       plotlyOutput("tax_bar")),
              tabPanel("Heatmap", value = 6,
                       plotlyOutput("tax_heat", height = "750px", width = "1000px")),
          #   #tabPanel("Krona", value = 7),
          #   #tabPanel("Export", value = 8),
          id = "tabselected"
          )
  )
))

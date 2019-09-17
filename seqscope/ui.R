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
                 
           headerPanel("SeqScope"),
           sidebarPanel(
             
             ## conditionalPanel() functions for selected tab
             conditionalPanel(condition="input.tabselected==1"),
             # For Panel 1, have an input option
             conditionalPanel(condition = "input.tabselected == 2",
                              h3("Run with demo or real dataset?"),
                              
                              radioButtons("mode", label = "",
                                           choices = c("Demo", "Real"), selected = "Demo"),
                              uiOutput("biomSelect"),
                              uiOutput("metaSelect"),
                              h3("Press the button below to process data"),
                              actionButton("go", "Process"),
                              textOutput("fileStatus")
             ),    
             # On panels 5 and 6 (barplot and heatmap), ask which taxonomic level they want to visualize to
             conditionalPanel(condition = "input.tabselected == 5 | input.tabselected == 6",
                              uiOutput("which_taxon_level")),
             conditionalPanel(condition = "input.tabselected == 6",
                              uiOutput("select_species_heat"))
           ),

      mainPanel(
        tabsetPanel(
          tabPanel("Welcome",value=1,
                     includeMarkdown("docs/welcome.md")),
          tabPanel("Data", value=2,
                   h2("Please verify that the files below look as expected, and click on 'Run the app!' to get started!"),
                   h4("Input taxonomy file"),
                   DT::dataTableOutput("print_taxon_table"),
                   h4("Input metadata file"),
                   DT::dataTableOutput("print_metadata_table"),
                   h4("Input QC file"),
                   DT::dataTableOutput("print_qc_table")
                   ),
          #   #tabPanel("QC", value = 3),
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


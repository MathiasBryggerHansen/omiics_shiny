library(shiny)

ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
        ),
        mainPanel = mainPanel(
            dataTableOutput("gene_results_table")
        )
    )
)

server <- function(input, output) {
    gene_results_2 <- reactive(readRDS("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/DETable/data/gene_results_table.RDS"))
    output$gene_results_table <- renderDataTable({
        datatable(data = gene_results_2(),caption = "DE results including any extra datasets",filter = list(position = 'top'),escape = F, options = list(autoWidth = TRUE,scrollX = TRUE, columnDefs = list(list(width = '400px', targets = grep(colnames(gene_results_2()),pattern = "reactome|kegg|carta|wiki")))))
    })
}

shinyApp(ui = ui, server = server)

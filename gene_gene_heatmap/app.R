
#setwd("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/volcano/")

ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            textInput(inputId = "col_high", label = "Color of continuous high values", value = "red"),
            textInput(inputId = "col_low", label = "Color of continuous low values", value = "blue"),
        ),
        mainPanel(
            plotlyOutput("gene_gene",height = "500px")
        )
    )
)

server <- function(input, output) {
    gene_gene_cor <- reactive(readRDS("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/gene_gene_heatmap/data/gene_gene.RDS"))
    #pathway_dic <- reactive(readRDS("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/volcano/data/gene_dic_human.RDS"))
    output$gene_gene <- renderPlotly({
        #cor(data.frame(t(data())),method = "spearman")
        heatmaply(gene_gene_cor(),scale_fill_gradient_fun = scale_fill_gradient(low = input$col_low, high = input$col_high))
    })
}

shinyApp(ui = ui, server = server)

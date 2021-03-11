library(shiny)
#setwd("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/volcano/")

ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            textInput(inputId = "p", label = "p cutoff:",value = '10**-14'),
            numericInput(inputId = "fc", label = "FC cutoff:", value = 2),
            textInput(inputId = "experiment_id", label = "Use alternative dataset for volcano plot", value = ""),
            textInput(inputId = "volcano_col", label = "A column to annotate color in volcano plot",value = "gene_biotype"),
            textInput(inputId = "col_high", label = "Color of continuous high values", value = "red"),
            textInput(inputId = "col_low", label = "Color of continuous low values", value = "blue"),
            checkboxInput(inputId = "log_scale", label = "Scale color values"),
        ),
        mainPanel(
            plotlyOutput("volcano",height = "500px")
        )
    )
)

server <- function(input, output) {
    gene_results_filtered <- reactive(readRDS("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/volcano/data/volcano_data.RDS"))
    pathway_dic <- reactive(readRDS("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/volcano/data/gene_dic_human.RDS"))

    output$volcano <- renderPlotly({
        print((colnames(gene_results_filtered())))
        req(gene_results_filtered(), eval(parse(text = input$p))<0.005)
        volcano_plot(input = input, data = gene_results_filtered(), pathway_dic = pathway_dic())})
}

shinyApp(ui = ui, server = server)

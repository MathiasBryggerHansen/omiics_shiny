library(shiny)

ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
        ),
        mainPanel = mainPanel(
            plotlyOutput("pca",height = "500px")
        )
    )
)

server <- function(input, output) {
    pca_scores <- reactive(readRDS("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/pca/data/pca_scores.RDS"))
    pca_ann <- reactive(readRDS("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/pca/data/pca_ann.RDS"))
    #pathway_dic <- reactive(readRDS("C:/Users/Mathias/OneDrive/Documents/omiics/omiics_shiny/volcano/data/gene_dic_human.RDS"))
    output$pca <- renderPlotly({
        plot_ly(x=pca_scores()[,1], y=pca_scores()[,2], z=pca_scores()[,3], type = "scatter3d", mode="markers",name = pca_ann()[,2], color = pca_ann()[,1]) %>%
            layout(scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3")))
        })
}

shinyApp(ui = ui, server = server)

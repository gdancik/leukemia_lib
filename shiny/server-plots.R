library(ggplot2)

source("functions/queryMongoGeneratePlots.R")

plotMultiple <- function() {
  plots <- query_mongo(input$inputGene, "risk", datasets = "all", input$inputTest)
  plots_km <- query_mongo_km(input$inputGene, datasets = "all")
  forest_plot <- query_mongo_forest(input$inputGene, datasets = "all")
  
  output$riskPlots <- renderUI({
    if (all(is.na(plots))) {h1("No risk data to display.")}
    lapply(1:length(plots), function(i) {
      if (is.na(plots[[i]][1])) (return())
      id <- paste0("risk_", i)
      col <- column(width = 3, shinycssloaders::withSpinner(plotOutput(id)))
      output[[id]] <- renderPlot({
        plots[[i]]
      })
      col
      
    }
    )
  })
  
  output$survPlots <- renderUI({
    if (all(is.na(plots_km))) {h1("No survival data to display.")}
    lapply(1:length(plots_km), function(i) {
      if (is.na(plots_km[[i]][1])) (return())
      id <- paste0("surv_", i)
      col <- column(width = 3, shinycssloaders::withSpinner(plotOutput(id)))
      output[[id]] <- renderPlot({
        plots_km[[i]]
      })
      col
      
    }
    )
  })
  
  output$forestPlots <- renderUI({
    if (all(is.na(forest_plot))) {h1("No survival data to display.")}
    rendered_forest_plot <- shinycssloaders::withSpinner(plotOutput("forest_1"))
    output[["forest_1"]] <- renderPlot({
      forest_plot
    })
    rendered_forest_plot
  })
}
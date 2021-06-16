library(shiny)
library(shinythemes)

# Define UI for application

ui <- fluidPage(theme = shinytheme("cerulean"),
  titlePanel("Signaux des lymphocytes Th expliqués par les signaux des cellules dendritiques"),
  hr(),
  sidebarLayout(
    sidebarPanel(
      selectInput("lth", "Signal de lymphocytes Th : ", choices = colnames(data8)),
      selectInput("cd2", "Signal de cellules dendritiques : ", choices = colnames(data7))
    ),
    mainPanel(
      plotOutput("distPlot")
    )
  ),
  hr(),
  sidebarPanel(
    selectInput("lth2", "Signal de lymphocytes Th : ", choices = colnames(data8)),
    selectInput("cd", "Signal de cellules dendritiques : ", choices = colnames(data7)),
    sliderInput("n_group", label = "Nombre de classes : ", min = 2, max = 20, value = 6),
    tableOutput("tablecoef")
  ),
  mainPanel(
    textOutput("texte"),
    plotOutput("boxplot2")
  )
)

# Define server logic 
server <- function(input, output) {
  
  datacoef <- reactive({
    datata <- as.data.frame(beta1)
    datata
  })
  
  
  beta1 <- reactive({
    data4$ILscale <- scale(data4[,input$lth])
    data_sel <- data4[,c(input$lth, selectio)]
    mod_comp <- lm(paste(input$lth,"~."), data = data_sel)
    data4 <- data4[order(residuals(mod_comp)),]
    x_scale <- as.matrix(Diagonal(unlist(data4[,input$cd2]), n = length(data4$IL4)))
    mod_fused <- fusedlasso(data4$ILscale, x_scale, graph = gr, minlam = 5*(mean(unlist(data4[,input$lth]))))
    coef(mod_fused)$beta
  })
  
  
  
  
  
  output$distPlot <- renderPlot({
    #lt <- paste(input$lth)
    
    beta1() %>% as.data.frame() %>% 
      rowid_to_column() %>%
      gather(-rowid, key = "lambda", value = "coef") %>%
      ggplot(aes(x = log(as.numeric(lambda)), y = coef, color = as.character(rowid))) + geom_line() + theme(legend.position = "none")
  }, height = 400, width = 800)
  
  
  output$texte <- renderText({
    paste(input$cd, "," , input$n_group, " groupes")
  })
  
  
  beta2 <- reactive({
    data4$ILscale <- scale(data4[,input$lth2])
    data_sel <- data4[,c(input$lth2, selectio)]
    mod_comp <- lm(paste(input$lth2,"~."), data = data_sel)
    data4 <- data4[order(residuals(mod_comp)),]
    x_scale <- as.matrix(Diagonal(unlist(data4[,input$cd]), n = length(data4$IL4)))
    mod_fused <- fusedlasso(data4$ILscale, x_scale, graph = gr, minlam = 5*(mean(unlist(data4[,input$lth2]))))
    coef(mod_fused)$beta
  })

  
  output$boxplot2 <- renderPlot({
    data4$group <- as.character(as.numeric(as.factor(beta2()[,input$n_group])))

    ggboxplot(data4, x = "group", y = input$cd,
              color = "group",
              add = "jitter", shape = "group")
  }, height = 500, width = 800) 
  
  
  
  output$tablecoef <- renderTable({
    
    data4$ILscale <- scale(data4[,input$lth])
    data_sel <- data4[,c(input$lth, selectio)]
    mod_comp <- lm(paste(input$lth,"~."), data = data_sel)
    data4 <- data4[order(residuals(mod_comp)),]
    x_scale <- as.matrix(Diagonal(unlist(data4[,input$cd2]), n = length(data4$IL4)))
    mod_fused <- fusedlasso(data4$ILscale, x_scale, graph = gr, minlam = 5*(mean(unlist(data4[,input$lth]))))
    beta1 <- coef(mod_fused)$beta
    datacoef3 <- as.data.frame(format(unique(beta1[,input$n_group])*sd(unlist(data4[,input$lth])), scientific = TRUE))
    #datacoef3 <- as.data.frame(format(unique(beta1[,input$n_group]), scientific = TRUE))
    #datacoef3 <- t(datacoef3)
    datacoef3$Groupe <- as.character(as.numeric(as.factor(unique(beta1[,input$n_group]))))
    colnames(datacoef3) <- c("Coefficient", "Groupe")
    datacoef3 <- datacoef3[,c("Groupe", "Coefficient")]
    
    datacoef3
    #DT::datatable(datacoef3)
  })
  

  
}


# Run the application 
shinyApp(ui = ui, server = server)

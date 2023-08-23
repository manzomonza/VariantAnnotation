## HORAK scoring


# Read in cancerHotspots
# Read gnomad MAF

# generate scores from google sheets -- Horak scoring
#Libraries
library(shiny)
library(shinydashboard)
library(purrr)
library(dplyr)
library(shinyWidgets)
library(DT)
library(shinyBS)

tabledata1 <- data.frame(country = c("Australia","Australia","Australia","Australia","Australia"), Code =c("Rahil","Rahil","Rahil","Rahil","Rahil") ,cat = c(1,1,1,1,1))
un_sq <- unique(tabledata1$Code)
tempa <- list()

ui <- dashboardPage(
  dashboardHeader(title = "Data Updation with collapsibility functionality" , disable = FALSE),
  dashboardSidebar(),
  dashboardBody(

    fluidRow(
      box(
        title = "Table Structure", status = "primary", solidHeader = TRUE,
        collapsible = TRUE,width = 12,uiOutput("Works")
      )
    )
  )
)
server <- function(input,output,session){
  output$Works <- renderUI({
    lapply( 1:length(un_sq),function(i) {
      bsCollapse(id = "collapseExample",
                 bsCollapsePanel(un_sq[i],  dataTableOutput(paste("mytable", i , sep ="")), style = "primary")
      )

    }
    )})


  shinyValue = function(id, len) {
    unlist(lapply(seq_len(len), function(i) {
      value = input[[paste0(id, i)]]
      if (is.null(value)) NA else value
    }))
  }
  Rahil = function(FUN, len, id, ...) {
    inputs = character(len)
    for (i in seq_len(len)) {
      inputs[i] = as.character(FUN(paste0(id, i), label = NULL, ...))
    }
    inputs
  }
  th <- function(){

    for(i in 1:nrow(tabledata1)){

      tempa[[i]]=data.frame(tabledata1[i,],Answer=   Rahil(radioButtons , nrow(tabledata1[i,]),paste0("radio" , i),selected ="yes" ,
                                                           choices =c("yes", "no"), inline= T),OP =  tabledata1$cat[i])

    }
    return(do.call(rbind,tempa))


  }


  rea1 <- reactiveValues(df = th()
  )
  rea2<- reactive({rea1$df})

  observeEvent(input$radio11,{

    if(rea1$df$cat[1] == 1){
      rea1$df$cat[1] =2
      print(rea1$df$cat[1])
    }else{
      rea1$df$cat[1]=1
      print(rea1$df$cat[1])
    }
    rea2()


  })



  output$mytable1 = DT::renderDataTable({rea2()},selection='none',server = FALSE, escape = FALSE,class = 'cell-border stripe', options = list(columnDefs = list(list(width = '600px', targets = 2),list(visible=FALSE)), ordering=F,pageLength = 10000,   lengthMenu = c(5, 10, 20, 100, 1000, 10000) , dom ="t",
                                                                                                                                              preDrawCallback = JS('function() {
                                                                                                                                                                 Shiny.unbindAll(this.api().table().node()); }'),
                                                                                                                                              drawCallback = JS('function() {
                                                                                                                                                              Shiny.bindAll(this.api().table().node()); } ')))

}


shinyApp(ui, server)

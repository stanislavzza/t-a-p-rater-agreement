
data_load_ui <- function(id) {
  
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel( 
      fileInput(ns("dataFileUploader"), "Choose CSV File", accept = ".csv"),
      uiOutput(ns("dataType")),
      uiOutput(ns("subjectID")),
      uiOutput(ns("ratingNames")), # are ratings in their own columns?
      uiOutput(ns("ratingVals")), # values for the ratings
      uiOutput(ns("dataDoneButton"))
    ),
    mainPanel(
      fluidRow(
        column(12, tableOutput(ns("dataPanelHeader"))),
        column(12, textOutput(ns("dataTypeInstructions"))),
        column(12, tableOutput(ns("dataTypeExample"))),
        column(12, hr()),
        column(12, tableOutput(ns("varStats")))
      )
    )
  )
}

data_load_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    #INPUT Reactive expression to read the CSV file
    raw_data <- reactive({
      req(input$dataFileUploader)
      read_csv(input$dataFileUploader$datapath) |> 
        mutate(across(everything(), as.character))
    })
    
    #INPUT--------------------- Select the data format -----------------------
    output$dataType <- renderUI({ 
      req(raw_data())
      
      radio_choices <- c("Outcome Ratings in Columns"="wide","Rater Ratings in Columns"="raters","Outcome Names in a Column"="narrow")
      radioButtons(ns("dataType"), "Type of Data", radio_choices, selected = "wide", inline = FALSE, width = NULL)
    })
    
    #OUTPUT--------------------- Some instructions -----------------------
    output$dataTypeInstructions <- renderText({
      req(raw_data())
      return("Use this option if your data looks like the arrangement below.")
    })
    
    #OUTPUT--------------------- Display sample data type -----------------------
    output$dataTypeExample <- renderTable({
      req(raw_data(), input$dataType)
      #if(!is.null(processed_data())) return(NULL)
      
      # data type examples
      wide <- data.frame(SubjectID=c(111,111,112,112),Objective1=c(3,3,2,1),Objective2=c(1,4,5,2),Objective3=c(2,2,1,3))
      rater <- data.frame(Optional_SubjectID=c(111,112,113,114),Optional_Category=c(111,112,113,114),Rater1=c(3,3,2,1),Rater2=c(1,4,5,2),Rater3=c(2,2,1,3))
      narrow <- data.frame(SubjectID=c(111,111,112,112),ObjectiveName=c("Writing","Speaking","Writing","Speaking"),Ratings=c(1,4,5,2))
      
      if (input$dataType == "wide"){
        wide
      } else if(input$dataType == "narrow") {
        narrow
      } else {
        rater
      }
    },digits = 0)
    
    #INPUT--------------------- Select the Subject ID column-----------------------
    output$subjectID <- renderUI({ 
      req(raw_data())
      req(input$dataType)

      namelist <- as.list(sort(names(raw_data())))
      if(input$dataType =="raters") namelist <- c("[None]",namelist)
      selectInput(ns("subjectID"), "Subject ID",  namelist )
    })
    

    #INPUT--------------------- Select Ratings names column-----------------------
    output$ratingNames <- renderUI({ 
      req(raw_data(), input$dataType, input$subjectID)
      if(input$dataType == "wide") return()
      namelist <- as.list(sort(names(raw_data())))
      namelist <- namelist[namelist != input$subjectID]
      if(input$dataType == "raters") {
        namelist <- c("[None]",namelist)
      }
      selectInput(ns("ratingNames"), "Outcome Names",  namelist )
    })
    
    #INPUT--------------------- Select Ratings Values column-----------------------
    output$ratingVals <- renderUI({ 
      req(raw_data(), input$dataType, input$subjectID)
      if(input$dataType != "narrow") return()
      namelist <- as.list(sort(names(raw_data())))
      namelist <- namelist[namelist != input$subjectID]
      namelist <- namelist[namelist != input$ratingNames]
      selectInput(ns("ratingVals"), "Rating Values",  namelist )
    })
    
    #INPUT-----------------Button for finished loading data---------------
    output$dataDoneButton <- renderUI({
      req(raw_data(), input$dataType)
      actionButton(ns("dataDoneUploading"),"Done Uploading")
    })
    
    #PROCESS -- standardize the data for use based on inputs
    processed_data <- eventReactive(input$dataDoneUploading, {

      req(raw_data(), input$dataType) #, input$ratingNames, input$ratingVals)

      df <- raw_data()
      
      # need this in different places below
      fix_id <- function(df){
        if (input$subjectID == "[None]") {
          df$SubjectID__ <- 1:nrow(df)  # create an ID manually
        } else {
          df <- df |> rename(SubjectID__ = input$subjectID)
        }
        return(df)
      }
      
      # format for the chosen data input type
      if (input$dataType == "narrow"){
        
        df$ID__temp__ <- 1:nrow(df) # needed to make each row unique
        df <- spread(df, input$ratingNames,input$ratingVals)
        df$ID__temp__ <- NULL # get rid of the temp ID
        df <- fix_id(df) # standardize name of the subject ID column
        
        # cf http://stackoverflow.com/questions/21390141/specify-dplyr-column-names
      } else if(input$dataType == "raters"){
        
        df <- fix_id(df) # standardize name of the subject ID column
        
        # if there's a column for the dimension, leave it alone
        if(input$ratingNames == "[None]") {
          df <- gather(df, Rater, Rating, -SubjectID__)
        } else {
          df <- gather(df, Rater, Rating, -SubjectID__, -input$ratingNames) |> 
                spread(input$ratingNames,Rating)
        }
        
        df$Rater <- NULL # don't need this
        
      } else {
        df <- fix_id(df) # standardize name of the subject ID column
      }
      
      return(df)

    })
    
    output$varStats <- renderTable({
      req(processed_data())

      df <- processed_data()

      dftype <- c()
      dfN <- c()
      dfSummary <- c()
      dfNA <- c()
      dfUnique <- c()
      dfAction <- c()
      dfmatched <- c()
      
      dfnames <- as.character(sort(names(df)))
      for (t in dfnames) {
        dftype <- c(dftype,class(df[[t]]))
        
        dfN <- c(dfN, sum(!is.na(df[[t]])) )
        
        u <- length(unique(df[[t]][!is.null(df[[t]])]))
        dfUnique <- c(dfUnique, u )
        
        if (u < 11 && u > 1) { # only allow scores with up to 10 unique response values
          sum_tbl <- rowSums(table(df$SubjectID__,df[[t]]))
          m <- sum(sum_tbl[sum_tbl > 0] - 1)
          dfmatched <- c(dfmatched, m  )
        } else {
          m <- 0
          dfmatched  <-  c(dfmatched, NA )
        }
        
        # drop the column of data if it's useless. We require at least 5 shared ratings.
        if (t != "SubjectID__" && (u > 10 || m < 5)){
          df[[t]] <<- NULL
          dfAction <- c(dfAction,"Drop")
        } else {
          dfAction <- c(dfAction,"Keep")
        }
        
      }
      # return value
      data.frame(Variable = dfnames, 
                 Type = dftype, 
                 N = dfN, 
                 Unique = dfUnique, 
                 Matched = dfmatched, 
                 Action = dfAction)
      
    }, digits = 0)
    
    return(processed_data)
  })
}
ordinal_analysis_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel( 
      uiOutput(ns("targetVar")), 
      uiOutput(ns("action_tap"))
    ),
    
    mainPanel(
      h1("Ordinal parameter estimates"),
      plotOutput(ns("ordinal_distro"),height = 200),
      plotOutput(ns("ordinal_tap"),height = 600)
      
    )
  )
}

ordinal_analysis_server <- function(id, ratings) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rating_counts <- reactiveVal(NULL)
    ordinal_tap   <- reactiveVal(NULL)

    #INPUT--------------------- Select the Target Var to predict-----------------------
    output$targetVar <- renderUI({ 
      # turning off the newvars update for now--it's quite annoying
      #rvals$newVarExists # update the list of names if a new var is created
      req(ratings())
      namelist <- as.list(sort(names(ratings())))
      namelist <- namelist[namelist != "SubjectID__"]
      select <- isolate(input$targetVar) # save the old choice if there is one
      selectInput(ns("targetVar"), "Ratings",  choices = namelist , selected = select )
    })
    
    #INPUT -------------------- Button to compute the tap stats
    output$action_tap <- renderUI({
      req(ratings(), input$targetVar)
      actionButton(ns("action_tap"),"Compute")
    })
    
    #Process the data on action
    observeEvent(input$action_tap, {
      req(ratings(), input$targetVar)
      
      # ratings need to be filtered to the target var
      # set the reactive values to the output   
      tdf <- ratings() |> 
             select(SubjectID__, rating = input$targetVar) |> 
             na.omit()
     
      ordinal_tap(get_ordinal_tap(tdf))
      
    })
    
    # show the distribution of ratings for the chosen variable
    output$ordinal_distro <- renderPlot({
      
      req(ratings(), input$targetVar)

      # select the right column and count values
      tdf <- ratings() |> 
        select(SubjectID__, rating = input$targetVar) |> 
        na.omit() |> 
        count(rating)
      
      # prepare and plot a stacked bar chart
      tdf |> 
        mutate(inclass = cumsum(n),
               outclass = sum(n) - inclass,
               cutpoint = str_c(rating,"|", lead(rating)),
               p = scales::percent(inclass / (inclass + outclass))) |> 
        na.omit() |> 
        select(cutpoint, inclass, outclass, p) |> 
        gather(key = "class", value = "count", -cutpoint, -p) |>
        mutate(p = if_else(class == "inclass", p, ""),
               class = reorder(class, desc(class))) |> 
        ggplot(aes(x = cutpoint, y = count, fill = class, label = p)) +
        geom_col() +
        geom_text(position = position_stack(vjust = 0.5), color = "black") +
        theme_bw() +
        # set manual colors
        scale_fill_manual(values = c("inclass" = "#CCCCCC", "outclass" = "#888888")) +
        xlab("CutPoint") +
        theme(text=element_text(size=20))

    })
    

    output$ordinal_tap <- renderPlot({
      
      req(ordinal_tap())
      
      tap_stats <- ordinal_tap()
      
      tap_stats |> 
        select(CutPoint, a, t, p, degenerate, type) |> 
        # change the order to tap is on top and solid
        mutate(type = reorder(type, desc(type))) |> 
        gather(param, value, a, t, p) |> 
        ggplot(aes(x = CutPoint,
                   y = value,
                   shape = degenerate,
                   linetype = type,
                   group = type)) +
        geom_point(size = 3) +
        scale_color_manual(values = c("#222222", "#d7301f"), limits = c(FALSE, TRUE), drop = TRUE) +
        # scale_colour_brewer(palette = "Oranges", type = "qual") +
        scale_shape_manual(values=c(16, 1)) +
        geom_line() +
        #   geom_hline(yintercept = .5, linetype = "dashed") +
        #geom_hline(yintercept = .25, linetype = "dashed") +
        theme_bw()+
        facet_grid(param ~ .) +
        theme(text=element_text(size=20), strip.text.y = element_text(size = 20)) +
        ylab("Parameter estimate") +
        ggtitle("Inclass is x <= CutPoint") 
  
    })

  })
}
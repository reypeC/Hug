#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# Package names
packages <- c("shiny", "shinyFiles", "dendextend", "fields")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

colfunc<-paste0(colorRampPalette(c("white","royalblue","yellow","red",alpha=0.5))(10),"75")[-c(9,10)]

combined_result <- NULL


reactive_state <- reactiveValues(
  numberPlanes = NULL,
  planes = NULL,
  normalisedData = NULL,
  averageDist = 0.1,
  averageDistByPlane = NULL,
  select_folder_output = "~/",
  radius_planes_string_text = NULL
)

reactive_state_results<-reactiveValues(
  reverseNormalisedData = NULL,
  userNumericBounds = NULL,
  sourcesDetection = NULL,
  reverseSourcesDetection = NULL,
  statisticsByPlan = NULL,
  numberPlanes = NULL,
  numberPatternCooling = NULL,
  planes = NULL,
  normalisedData = NULL,
  rawData = NULL,
  centers = NULL,
  reverseCenters = NULL,
  proposedSources = NULL,
  reverseProposedSources = NULL,
  numberPatternDetection = NULL
  
)

ui <- navbarPage("Hug procedure",
                 
                 
                 ########################################################################################
                 ########################################################################################
                 ########################################################################################
                 
                 #                     Data Initialisation
                 
                 ########################################################################################
                 ########################################################################################
                 ########################################################################################
                 
                 tabPanel("Data Initialisation",
                          titlePanel("Data Initialisation"),
                          fluidPage(
                            titlePanel("Data Initialisation"),
                            fluidRow(
                              column(4,
                                     fileInput("file", "Choose .txt File", 
                                               accept = c("text", "text/comma-separated-values,text/plain", ".txt"))
                              ),
                              column(4,
                                     actionButton("normalise_btn", "Normalise Data")
                              ),
                              column(4,
                                     downloadButton("save_normalised_data", "Saved normalised data as .txt")
                              )
                            ),
                            h3("Lower and Upper Bounds for each Dimension"),
                            uiOutput("value_inputs"),       # Set bounds
                            
                            fluidRow(
                              column(6,
                                     h2("Raw Data Preview : "),
                                     tableOutput("filePreview")
                              ),
                              column(6,
                                     h2("Normalised Data Preview : "),
                                     tableOutput("normalisedData")
                              )
                            )
                            
                          )),
                 
                 
                 ########################################################################################
                 ########################################################################################
                 ########################################################################################
                 
                 #                     Region of Interest
                 
                 ########################################################################################
                 ########################################################################################
                 ########################################################################################
                 
                 
                 tabPanel("Region of Interest" ,
                          titlePanel("Region of Interest"),
                          fluidPage(
                            titlePanel("Definition of the Region of Interests"),
                            
                            sidebarLayout(
                              sidebarPanel(
                                br(), br(),
                                h3("Aggregated Centers from All Planes"),
                                tableOutput("combined_table"),
                                downloadButton("save_combined", "Saved Aggregated  centers as .txt"),
                                br(),
                                # actionButton("done", "Finish and Return Combined Data")
                              ),
                              mainPanel(
                                uiOutput("plots_ui")  # Dynamic UI for plots
                              )
                            )
                          )),
                 
                 
                 ########################################################################################
                 ########################################################################################
                 ########################################################################################
                 
                 #                     Parameters File
                 
                 ########################################################################################
                 ########################################################################################
                 ########################################################################################
                 
                 
                 tabPanel("Parameters File" ,
                          titlePanel("Parameters File"),
                          fluidPage(
                            titlePanel("Creation of the Parameters File"),
                            h1("Set seed :"),
                            fluidRow(
                              column(4,
                                     numericInput("val1", "seed 1", value = 4451)
                              ),
                              column(4,
                                     numericInput("val2", "seed 2", value = 2289)
                              ),
                              column(4,
                                     numericInput("val3", "seed 3", value = 5791)
                              )
                            ),
                            br(),
                            actionButton("randomise_btn", "Randomise seed"),
                            br(),br(),
                            h1("Set simulation parameters :"),
                            fluidRow(
                              column(4,
                                     numericInput("mH", "number of MH per iteration :", value = 50)
                              ),
                              column(4,
                                     numericInput("nIterSave", "number iterations between save :", value = 100)
                              ),
                              column(4,
                                     numericInput("nDetectedSource", "number of pattern for the detection :", value = 1000)
                              )
                            ),
                            br(),
                            h1("Set cooling scheme :"),
                            fluidRow(
                              column(4,
                                     numericInput("tI", "initial temperature :", value = 1000.0)
                              ),
                              column(4,
                                     numericInput("tF", "final temperature :", value = 0.0001)
                              ),
                              column(4,
                                     numericInput("c", "cooling coefficient :", value = 0.9999)
                              )
                            ),
                            br(),
                            h1("Set events :"),
                            fluidRow(
                              column(5,
                                     numericInput("pb", "probability of birth", value = 0.2)
                              ),
                              column(5,
                                     numericInput("pd", "probability of death", value = 0.2)
                              ),
                              column(5,
                                     numericInput("pcI", "probability of change in the region of interest", value = 0.2)
                              ),
                              column(5,
                                     numericInput("pcO", "probability of change out of the region of interest", value = 0.05)
                              ),
                              column(5,
                                     numericInput("rChange", "radius of change", value = 0.1)
                              )
                            ),
                            br(),
                            h1("names of the data files : "),
                            fluidPage(
                              fluidRow(
                                column(5,
                                       shinyFilesButton("data_file", "Choose .txt Data File", "Please select the data file", multiple = FALSE)
                                ),
                                column(5,
                                       textOutput("data_file_path_output")
                                )
                              ),
                              fluidRow(
                                column(5,
                                       shinyFilesButton("region_of_interest_file", "Choose .txt Region of Interest File", "Please select the region of interest file", multiple = FALSE)
                                ),
                                column(5,
                                       textOutput("region_of_interest_file_path_output")
                                )
                              )
                            ),
                            br(),
                            h1("Set value for the statistics :"),
                            fluidRow(
                              column(5,
                                     sliderInput("alpha", "percentage of data not explained on each plane :", min = 0, max = 1, value = 0.05,step = 0.001)
                              ),
                              column(5,
                                     numericInput("theta_a", "theta_a :", value = 1000.0)
                              )
                            ),
                            uiOutput("number_dimension"),
                            uiOutput("radius_interaction"),
                            
                            fluidRow(
                              column(3, h3("considered planes")),
                              column(9, textOutput("plane_string_text"))
                            ),
                            uiOutput("radius_plane_sliders"),
                          
                            h1("names of the saved files : "),
                            fluidPage(
                              column(4,
                                     shinyDirButton("folder", "Choose Folder for the Results", "Please select a folder")
                              ),
                              column(5,
                                     verbatimTextOutput("folder_path")
                              )
                            ),
                            textInput("filename_sources", "name of the file of saved sources generated during the cooling :", value = "sources"),
                            textInput("filename_sources_final", "name of the file of saved sources used for the final detection :", value = "sources_final"),
                            textInput("filename_statistics", "name of the file of saved statistics generated during the cooling :", value = "statistics"),
                            textInput("filename_statistics_final", "name of the file of saved statistics used for the final detection :", value = "statistics_final"),
                            textInput("filename_misc", "name of the miscelleanous file:", value = "misc"),
                            br(),
                            downloadButton("save_parameters", "Saved parameters file as .txt"),
                            br()
                            
                            
                          )),
                 
                 
                 ########################################################################################
                 ########################################################################################
                 ########################################################################################
                 
                 #                     Analyse Results Hug
                 
                 ########################################################################################
                 ########################################################################################
                 ########################################################################################
                 
                 
                 tabPanel("Analyse Results Hug",
                          titlePanel("Analyse Results Hug"),
                          fluidPage(
                            titlePanel("Import Results :"),
                            fluidRow(
                              column(5,
                                     shinyFilesButton("sources_final_file_results", "Choose .txt saved sources used for the final detection File", "Please select the saved sources used for the final detection file", multiple = FALSE)
                              ),
                              column(5,
                                     textOutput("file_path_sources_final_file_results_output")
                              )
                            ),
                            fluidRow(
                              column(5,
                                     shinyFilesButton("statistics_file_results", "Choose .txt saved statistics generated during the cooling", "Please select the saved statistics generated during the cooling file", multiple = FALSE)
                              ),
                              column(5,
                                     textOutput("file_path_statistics_file_results_output")
                              )
                            ),
                            fluidRow(
                              column(5,
                                     shinyFilesButton("statistics_final_file_results", "Choose .txt saved statistics used for the final detection File", "Please select the saved statistics used for the final detection file", multiple = FALSE)
                              ),
                              column(5,
                                     textOutput("file_path_statistics_final_file_results_output")
                              )
                            ),
                            br(),
                            actionButton("import_results_btn", "Import"),
                            br(),br(),
                            uiOutput("plot_stats_results"),
                            br(),
                            uiOutput("plot_stats_global_results"),
                            br()
                          ),
                          
                          fluidPage(
                            titlePanel("Import Data : "),
                            fluidRow(
                              column(5,
                                     fileInput("data_file_results", "Choose .txt used Data File",
                                               accept = c("text", "text/comma-separated-values,text/plain", ".txt"))
                                     
                              ),
                              column(5,
                                     fileInput("region_of_interest_file_results", "Choose .txt Region of Interest File",
                                               accept = c("text", "text/comma-separated-values,text/plain", ".txt"))
                              )
                            ),
                            fluidRow(
                              column(6,
                                     h3("Used Data : "),
                                     tableOutput("filePreview_import_data")
                                     
                              ),
                              column(6,
                                     h3("Centers of Region of Interest : "),
                                     tableOutput("filePreview_import_centers")
                              )
                            ),
                            
                            
                          ),
                          fluidPage(
                            br(),
                            br(),
                            checkboxInput("show_file_input", "Reverse the normalisation?", value = FALSE),
                            
                            conditionalPanel(
                              condition = "input.show_file_input == true",
                              fluidPage(
                                fluidRow(
                                  column(5,
                                         fileInput("data_file_results_raw", "Choose .txt Raw Data File", 
                                                   accept = c("text", "text/comma-separated-values,text/plain", ".txt"))
                                  ),
                                  column(5,
                                         actionButton("reverse_normalise_btn", "Reverse the normalisation of the Data")
                                  ),
                                ),
                                h3("Lower and Upper Bounds for each Dimension"),
                                uiOutput("value_inputs_results") ,      # Set bounds
                                tableOutput("user_values_results"),
                                
                                fluidRow(
                                  column(6,
                                         h3("Raw Data Preview : "),
                                         tableOutput("filePreview_results_raw")
                                  ),
                                  column(6,
                                         h3("Used Data after the Reverse Normalisation : "),
                                         tableOutput("filePreview_reverse_normalised_data")
                                  )
                                )
                                
                                
                                
                                
                              )
                            ),
                          ),
                          fluidPage(
                            br(),
                            uiOutput("number_sources_detected_slider"),
                            br(),br(),
                            uiOutput("plot_results_2d"),
                            br()
                          ),fluidPage(
                            br(),
                            uiOutput("plot_reverse_results_2d"),
                            br()
                          )
                 )
                 
                 
)

server <- function(input, output, session) {
  
  # When browser session ends, stop the Shiny app
  session$onSessionEnded(stopApp)
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  #                     Data Initialisation
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  # Reactive: Read uploaded data
  uploaded_data <- reactive({
    req(input$file)
    
    data=read.table(input$file$datapath, header = T)
    
    return(data)
  })
  
  
  # Render uploaded file
  output$filePreview <- renderTable({
    req(input$file)
    head(round(uploaded_data(),5))
  })
  
  
  # Generate side-by-side inputs for each numeric column
  output$value_inputs <- renderUI({
    req(input$file)
    df <- uploaded_data()
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    
    inputs <- lapply(numeric_cols, function(col) {
      fluidRow(
        column(12, strong(col)),
        column(4,
               numericInput(paste0(col, "_val1"), "Lower Bound", value = min(df[[col]],na.rm=T)-(max(df[[col]],na.rm=T)-min(df[[col]],na.rm=T)))),
        column(4,
               numericInput(paste0(col, "_val2"), "Upper Bound", value = max(df[[col]],na.rm=T)+(max(df[[col]],na.rm=T)-min(df[[col]],na.rm=T)))),
        br()
      )
    })
    do.call(tagList, inputs)
  })
  
  
  
  # Normalise data when button is clicked
  normalised_result <- eventReactive(input$normalise_btn, {
    df <- uploaded_data()
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    norm_df <- df
    
    for (col in numeric_cols) {
      val1 <- input[[paste0(col, "_val1")]]
      val2 <- input[[paste0(col, "_val2")]]
      if (!is.null(val1) && !is.null(val2) && val2 != val1) {
        norm_df[[col]] <- (df[[col]] - val1) / (val2 - val1)
      } else {
        norm_df[[col]] <- NA  # Invalid bounds
      }
    }
    
    
    # Update reactive values instead of global ones
    reactive_state$normalisedData <- norm_df
    numberPlanes=ncol(norm_df) * (ncol(norm_df) - 1) / 2
    reactive_state$numberPlanes <- numberPlanes
    planes=combn(1:ncol(norm_df), 2)
    reactive_state$planes <- planes
    
    reactive_state$averageDist <- mean(dist(na.omit(norm_df)))
    
    averageDistByPlane=rep(0.1,reactive_state$numberPlanes)
    for(i in 1:numberPlanes)
    {
      averageDistByPlane[i]=mean(dist(na.omit(norm_df[,planes[,i]])))
    }
    reactive_state$averageDistByPlane <- averageDistByPlane
    
    
    norm_df
  })
  
  # Render normalised data
  output$normalisedData <- renderTable({
    req(normalised_result())
    head(round(normalised_result(),5))
  })
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  #                     Region of Interest
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  # List to store reactiveValues for each plot
  values_list <- reactiveValues()
  
  # Generate dynamic UI for plots
  output$plots_ui <- renderUI({
    if (is.null(reactive_state$numberPlanes) || reactive_state$numberPlanes < 1) return(NULL)
    
    plot_uis <- lapply(1:reactive_state$numberPlanes, function(i) {
      ns <- paste0("plot_", i)
      wellPanel(
        h4(paste("Plot", i)),
        checkboxInput(paste0("include_", i), "Include in combined table", value = TRUE),
        actionButton(paste0("clear_", i), "Clear Points"),
        downloadButton(paste0("save_plot_", i), "Save Plot as .png"),
        plotOutput(paste0("plot_", i), click = paste0("click_", i)),
        # tableOutput(paste0("table_", i))
        DT::dataTableOutput(paste0("editable_table_", i))
      )
    })
    
    do.call(tagList, plot_uis)
  })
  
  # Initialize reactiveValues for each plot on number change
  observeEvent(reactive_state$numberPlanes, {
    
    # Initialize an empty data frame for each plot if it doesn't exist
    for (i in 1:reactive_state$numberPlanes) {
      if (is.null(values_list[[paste0("df_", i)]])) {
        values_list[[paste0("df_", i)]] <- data.frame(coord1 = numeric(0), coord2 = numeric(0),dim1 = numeric(0), dim2 = numeric(0))
      }
    }
  })
  
  observe({
    if (is.null(reactive_state$numberPlanes) || reactive_state$numberPlanes < 1) return()
    
    for (i in 1:reactive_state$numberPlanes) {
      local({
        ii <- i
        df_name <- paste0("df_", ii)
        clear_btn <- paste0("clear_", ii)
        click_input <- paste0("click_", ii)
        plot_output <- paste0("plot_", ii)
        table_output <- paste0("table_", ii)
        save_plot <- paste0("save_plot_", ii)
        
        
        # Add new points on click
        observeEvent(input[[click_input]], {
          new_point <- data.frame(coord1 = input[[click_input]]$x, coord2 = input[[click_input]]$y, dim1=reactive_state$planes[1,ii], dim2=reactive_state$planes[2,ii])
          values_list[[df_name]] <- rbind(values_list[[df_name]], new_point)
        })
        
        # Clear points
        observeEvent(input[[clear_btn]], {
          values_list[[df_name]] <-  data.frame(coord1 = numeric(0), coord2 = numeric(0),dim1 = numeric(0), dim2 = numeric(0))
        })
        
        # Render plot
        output[[plot_output]] <- renderPlot({
          par(pty="s")
          plot(reactive_state$normalisedData[,reactive_state$planes[,ii]], main = paste("Plane number ", ii),xlim=c(0,1),ylim=c(0,1),pch=19)
          points(values_list[[df_name]], col = "red", pch = 19)
          if(nrow(values_list[[df_name]])>0)
          {
            for(indice in 1:nrow(values_list[[df_name]]))
            {
              symbols(values_list[[df_name]][indice,], circles =reactive_state$averageDistByPlane[ii], add = TRUE, fg = "red",bg = rgb(1, 0, 0, alpha = 0.1),inches = FALSE)
            }
          }
          
          
        })
        
        # Render table
        # output[[table_output]] <- renderTable({
        #   values_list[[df_name]]
        # })
        output[[paste0("editable_table_", ii)]] <- DT::renderDT({
          DT::datatable(
            values_list[[df_name]][,c(1,2)],
            editable = TRUE,  # Enable editing globally
            options = list(dom = 'tp', pageLength = 5),
            rownames = FALSE
          ) 
        })
        
        observeEvent(input[[paste0("editable_table_", ii, "_cell_edit")]], {
          info <- input[[paste0("editable_table_", ii, "_cell_edit")]]
          i <- info$row
          j <- info$col+1
          v <- info$value
          
          df <- values_list[[df_name]]
          df[i, j] <- as.numeric(v)
          values_list[[df_name]] <- df
        })
        
        output[[paste0("save_plot_", ii)]] <- downloadHandler(
          filename = function() {
            paste0("plot_", ii, ".png")
          },
          content = function(file) {
            png(file, width = 1000,height = 1000)
            par(pty="s")
            plot(reactive_state$normalisedData[, reactive_state$planes[, ii]],
                 main = paste("Plane number ", ii), xlim = c(0, 1), ylim = c(0, 1), pch = 19)
            points(values_list[[df_name]], col = "red", pch = 19)
            
            if (nrow(values_list[[df_name]]) > 0) {
              for (indice in 1:nrow(values_list[[df_name]])) {
                symbols(values_list[[df_name]][indice, ],
                        circles = reactive_state$averageDistByPlane[ii],
                        add = TRUE,
                        fg = "red",
                        bg = rgb(1, 0, 0, alpha = 0.1),
                        inches = FALSE)
              }
            }
            dev.off()
          }
        )
      })
    }
  })
  
  
  
  # Reactive expression combining all points with plot ID
  combined_points <- reactive({
    
    if (is.null(reactive_state$numberPlanes) || reactive_state$numberPlanes < 1) return(NULL)
    
    all_points <- do.call(rbind, lapply(1:reactive_state$numberPlanes, function(i) {
      df <- values_list[[paste0("df_", i)]]
      include <- input[[paste0("include_", i)]]
      if (is.null(df) || nrow(df) > 0 && !isTRUE(include)) return(NULL)
      df
    }))
    
    if (is.null(all_points)) {
      return(data.frame())
    } else {
      rownames(all_points) <- NULL
      return(all_points)
    }
  })
  
  # Output combined table in UI below all plots
  output$combined_table <- renderTable({
    combined_points()
  })
  
  output$save_combined <- downloadHandler(
    filename = function() {
      paste0("aggregated_centers", ".txt")
    },
    content = function(file) {
      # write.table( combined_points(), file, row.names=F, col.names=T, append=F )
      write.table(combined_points(), file, sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE, append=F)
    }
  )
  
  output$save_normalised_data <- downloadHandler(
    filename = function() {
      paste0("normalised_data", ".txt")
    },
    content = function(file) {
      write.table(reactive_state$normalisedData, file, sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE, append=F)
      
    }
  )
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  #                     Parameters File
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  roots <- c(home = normalizePath("~"))  # top-level definition
  
  
  #Save path to the data file
  data_file_path <- reactiveVal(NULL)
  
  observe({
    shinyFileChoose(input, "data_file", roots = roots)
    
    if (!is.null(input$data_file)) {
      parsed <- parseFilePaths(roots, input$data_file)
      selected_path <- as.character(parsed$datapath[1])
      
      if (!is.na(selected_path) && selected_path != "") {
        data_file_path(selected_path)
        
        norm_df <- read.table(data_file_path(), header = T)
        
        
        # Update reactive values instead of global ones
        reactive_state$normalisedData <- norm_df
        numberPlanes=ncol(norm_df) * (ncol(norm_df) - 1) / 2
        reactive_state$numberPlanes <- numberPlanes
        planes=combn(1:ncol(norm_df), 2)
        reactive_state$planes <- planes
        
        reactive_state$averageDist <- mean(dist(na.omit(norm_df)))
        
        averageDistByPlane=rep(0.1,reactive_state$numberPlanes)
        for(i in 1:numberPlanes)
        {
          averageDistByPlane[i]=mean(dist(na.omit(norm_df[,planes[,i]])))
        }
        reactive_state$averageDistByPlane <- averageDistByPlane
        
      }
    }
  })
  
  output$data_file_path_output <- renderText({
    req(data_file_path())
    data_file_path()
  })
  
  #Save path to the region of interest file
  region_of_interest_file_path <- reactiveVal(NULL)
  
  observe({
    shinyFileChoose(input, "region_of_interest_file", roots = roots)
    
    if (!is.null(input$region_of_interest_file)) {
      parsed <- parseFilePaths(roots, input$region_of_interest_file)
      selected_path <- as.character(parsed$datapath[1])
      
      if (!is.na(selected_path) && selected_path != "") {
        region_of_interest_file_path(selected_path)
      }
    }
  })
  
  output$region_of_interest_file_path_output <- renderText({
    req(region_of_interest_file_path())
    region_of_interest_file_path()
  })
  
  
  
  output$save_parameters <- downloadHandler(
    filename = function() {
      paste0("parameters", ".txt")
    },
    content = function(file) {
      
      file_to_save_parameter <- file(file, open = "wt")
      
      writeLines("//seed_A", file_to_save_parameter)
      writeLines(as.character(input$val1), file_to_save_parameter)
      writeLines("//seed_B", file_to_save_parameter)
      writeLines(as.character(input$val2), file_to_save_parameter)
      writeLines("//seed_C", file_to_save_parameter)
      writeLines(as.character(input$val3), file_to_save_parameter)
      
      writeLines("//number_iteration", file_to_save_parameter)
      writeLines(as.character(input$nDetectedSource), file_to_save_parameter)
      writeLines("//number_application_of_MH_algorithm_per_iteration", file_to_save_parameter)
      writeLines(as.character(input$mH), file_to_save_parameter)
      writeLines("//number_of_iteration_between_save", file_to_save_parameter)
      writeLines(as.character(input$nIterSave), file_to_save_parameter)
      
      writeLines("//initial_temperature", file_to_save_parameter)
      writeLines(as.character(input$tI), file_to_save_parameter)
      writeLines("//minimal_temperature", file_to_save_parameter)
      writeLines(as.character(input$tF), file_to_save_parameter)
      writeLines("//cooling_coefficient_temperature", file_to_save_parameter)
      writeLines(as.character(input$c), file_to_save_parameter)
      
      writeLines("//probability_of_birth", file_to_save_parameter)
      writeLines(as.character(input$pb), file_to_save_parameter)
      writeLines("//probability_of_death", file_to_save_parameter)
      writeLines(as.character(input$pd), file_to_save_parameter)
      writeLines("//probability_change_in_source", file_to_save_parameter)
      writeLines(as.character(input$pcI), file_to_save_parameter)
      writeLines("//probability_change_out_source", file_to_save_parameter)
      writeLines(as.character(input$pcO), file_to_save_parameter)
      writeLines("//radius_for_the_proposal_in_change_event", file_to_save_parameter)
      writeLines(as.character(input$rChange), file_to_save_parameter)
      
      writeLines("//mean_of_theta_for_number_of_sources_notin_disk_truesources_statistic", file_to_save_parameter)
      writeLines(as.character(input$theta_a), file_to_save_parameter)
      writeLines("//percentage_of_data_not_explained", file_to_save_parameter)
      writeLines(as.character(input$alpha), file_to_save_parameter)
      
      writeLines("//radius_planes_:_radius_disk_region_of_interest", file_to_save_parameter)
      writeLines(as.character(reactive_state$radius_planes_string_text), file_to_save_parameter)
      writeLines("//radius_of_interaction", file_to_save_parameter)
      writeLines(as.character(input$r), file_to_save_parameter)
      
      writeLines("//number_of_dimensions", file_to_save_parameter)
      writeLines(as.character(input$n), file_to_save_parameter)
      writeLines("//considered_planes", file_to_save_parameter)
      writeLines(as.character(reactive_state$plane_string_text), file_to_save_parameter)
      
      writeLines("//name_of_the_data_file", file_to_save_parameter)
      writeLines(data_file_path(), file_to_save_parameter)
      writeLines("//name_of_the_region_of_interest_file", file_to_save_parameter)
      writeLines(region_of_interest_file_path(), file_to_save_parameter)
      
      
      writeLines("//name_of_the_saved_sources_file", file_to_save_parameter)
      writeLines(paste0(reactive_state$select_folder_output,"/",input$filename_sources,".txt"), file_to_save_parameter)
      writeLines("//name_of_the_final_saved_sources_file", file_to_save_parameter)
      writeLines(paste0(reactive_state$select_folder_output,"/",input$filename_sources_final,".txt"), file_to_save_parameter)
      writeLines("//name_of_the_saved_statistics_file", file_to_save_parameter)
      writeLines(paste0(reactive_state$select_folder_output,"/",input$filename_statistics,".txt"), file_to_save_parameter)
      writeLines("//name_of_the_final_saved_statistics_file", file_to_save_parameter)
      writeLines(paste0(reactive_state$select_folder_output,"/",input$filename_statistics_final,".txt"), file_to_save_parameter)
      writeLines("//name_of_the_miscelleanous_file", file_to_save_parameter)
      writeLines(paste0(reactive_state$select_folder_output,"/",input$filename_misc,".txt"), file_to_save_parameter)
      
      close(file_to_save_parameter)
      
      
      
    }
  )
  

  
  # Observe randomise button and assign random values
  observeEvent(input$randomise_btn, {
    updateNumericInput(session, "val1", value = sample(1:10000,1))
    updateNumericInput(session, "val2", value = sample(1:10000,1))
    updateNumericInput(session, "val3", value = sample(1:10000,1))
  })
  
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  shinyDirChoose(input, "folder", roots = volumes, session = session)
  
  folder_path <- reactiveVal(NULL)
  
  observeEvent(input$folder, {
    path <- parseDirPath(volumes, input$folder)
    folder_path(path)
  })
  
  output$folder_path <- renderText({
    req(folder_path())
    reactive_state$select_folder_output<-folder_path()
    folder_path()
  })
  
  # Output the current values as a table
  output$random_values <- renderTable({
    data.frame(
      Parameter = c("seed 1", "seed 2", "seed 3"),
      Value = c(input$val1, input$val2, input$val3)
    )
  })
  
  output$number_dimension <- renderUI({
    req(reactive_state$normalisedData)  # Make sure numberPlanes is available
    
    sliderInput("n", "number of dimension", min = 0, max = ncol(reactive_state$normalisedData), value = ncol(reactive_state$normalisedData),round=T,step = 1)
    
    
    
  })
  
  output$radius_interaction <- renderUI({
    req(reactive_state$normalisedData)  # Make sure numberPlanes is available
    
    sliderInput("r", "radius of interaction", min = 0, max = 1, value = reactive_state$averageDist,step=0.00001)
    
    
  })
  
  output$radius_plane_sliders <- renderUI({
    req(reactive_state$numberPlanes)
    
    sliders <- lapply(1:reactive_state$numberPlanes, function(i) {
      include <- input[[paste0("include_", i)]]
      if (is.null(reactive_state$numberPlanes) || !isTRUE(include)) return(NULL)
      sliderInput(
        inputId = paste0("slider_", i),
        label = paste("Radius of Disk for the Region of Interest of Plane", i),
        min = 0,
        max = 1,
        value = reactive_state$averageDistByPlane[i],
        step = 0.00001
      )
    })
    
    do.call(tagList, sliders)
  })
  
  observe({
    req(reactive_state$numberPlanes)
    
    for (i in 1:reactive_state$numberPlanes) {
      local({
        ii <- i
        slider_id <- paste0("slider_", ii)
        
        observeEvent(input[[slider_id]], {
          reactive_state$averageDistByPlane[ii] <- input[[slider_id]]
        }, ignoreInit = TRUE)
      })
    }
  })
  
  output$plane_string_text <- renderText({
    req(reactive_state$planes)
    included_planes <- which(sapply(1:reactive_state$numberPlanes, function(i) {
      isTRUE(input[[paste0("include_", i)]])
    }))
    
    plane_string <- paste(
      apply(reactive_state$planes[, included_planes, drop = FALSE], 2, function(pair) {
        paste0("(", pair[1], ",", pair[2], ")")
      }),
      collapse = ";"
    )
    
    radius_string <- paste0(reactive_state$averageDistByPlane[ included_planes, drop = FALSE],collapse=";")
    
    
    # Save the string to reactive_state for access elsewhere
    reactive_state$plane_string_text <- plane_string
    reactive_state$radius_planes_string_text <- radius_string
    
    plane_string
  })
  
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  #                     Analyse Results Hug
  
  ########################################################################################
  ########################################################################################
  ########################################################################################
  
  
  # Reactive: Read uploaded data
  observeEvent(input$data_file_results, {
    req(input$data_file_results)
    
    data <- read.table(input$data_file_results$datapath, header = TRUE)
    
    reactive_state_results$normalisedData <- data
  })
  
  # Reactive: Read uploaded region of interest
  observeEvent(input$region_of_interest_file_results, {
    req(input$region_of_interest_file_results)
    
    data <- read.table(input$region_of_interest_file_results$datapath, header = TRUE)
    
    reactive_state_results$centers <- data
  })
  
  
  
  # Render uploaded file
  output$filePreview_import_data <- renderTable({
    req(input$data_file_results)
    if(is.null(reactive_state_results$normalisedData)) return(NULL)
    round(reactive_state_results$normalisedData[1:min(c(5,nrow(reactive_state_results$normalisedData))),],5)
  })
  
  output$filePreview_import_centers <- renderTable({
    req(input$region_of_interest_file_results)
    reactive_state_results$centers
  })
  
  
  # Reactive: Read uploaded data
  observeEvent(input$data_file_results_raw, {
    req(input$data_file_results_raw)
    
    data <- read.table(input$data_file_results_raw$datapath, header = TRUE)
    
    reactive_state_results$rawData <- data
  })
  
  
  # Render uploaded file
  output$filePreview_results_raw <- renderTable({
    req(input$data_file_results_raw)
    if(is.null(reactive_state_results$rawData)) return(NULL)
    round(reactive_state_results$rawData[1:min(c(5,nrow(reactive_state_results$rawData))),],5)
  })
  
  # Create table from inputs
  output$user_values_results <- renderTable({
    req(input$data_file_results_raw)
    df <- reactive_state_results$rawData
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    
    result <- data.frame(
      Column = character(),
      Value1 = numeric(),
      Value2 = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (col in numeric_cols) {
      val1 <- input[[paste0(col, "_val1")]]
      val2 <- input[[paste0(col, "_val2")]]
      if (!is.null(val1) && !is.null(val2)) {
        result <- rbind(result, data.frame(Column = col, LowerBound = val1, UpperBound = val2))
      }
    }
    
  })
  
  
  # Generate side-by-side inputs for each numeric column
  output$value_inputs_results <- renderUI({
    req(input$data_file_results_raw)
    df <- reactive_state_results$rawData
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    
    inputs <- lapply(numeric_cols, function(col) {
      fluidRow(
        column(12, strong(col)),
        column(6,
               numericInput(paste0(col, "_val1"), "Lower Bound", value = min(df[[col]],na.rm=T)-(max(df[[col]],na.rm=T)-min(df[[col]],na.rm=T)))),
        column(6,
               numericInput(paste0(col, "_val2"), "Upper Bound", value = max(df[[col]],na.rm=T)+(max(df[[col]],na.rm=T)-min(df[[col]],na.rm=T)))),
        br()
      )
    })
    
    do.call(tagList, inputs)
  })
  
  roots <- c(home = normalizePath("~"))  # top-level definition
  
  
  #Save path to the sources used for the detection file
  file_path_sources_final_file_results <- reactiveVal(NULL)
  
  observe({
    shinyFileChoose(input, "sources_final_file_results", roots = roots)
    
    if (!is.null(input$sources_final_file_results)) {
      parsed <- parseFilePaths(roots, input$sources_final_file_results)
      selected_path <- as.character(parsed$datapath[1])
      
      if (!is.na(selected_path) && selected_path != "") {
        file_path_sources_final_file_results(selected_path)
      }
    }
  })
  
  output$file_path_sources_final_file_results_output <- renderText({
    req(file_path_sources_final_file_results())
    file_path_sources_final_file_results()
  })
  
  #Save path to the statistics during the cooling file
  file_path_statistics_file_results <- reactiveVal(NULL)
  
  observe({
    shinyFileChoose(input, "statistics_file_results", roots = roots)
    
    if (!is.null(input$statistics_file_results)) {
      parsed <- parseFilePaths(roots, input$statistics_file_results)
      selected_path <- as.character(parsed$datapath[1])
      
      if (!is.na(selected_path) && selected_path != "") {
        file_path_statistics_file_results(selected_path)
      }
    }
  })
  
  output$file_path_statistics_file_results_output <- renderText({
    req(file_path_statistics_file_results())
    file_path_statistics_file_results()
  })
  
  #Save path to the statistics used for the detection file
  file_path_statistics_final_file_results <- reactiveVal(NULL)
  
  observe({
    shinyFileChoose(input, "statistics_final_file_results", roots = roots)
    
    if (!is.null(input$statistics_final_file_results)) {
      parsed <- parseFilePaths(roots, input$statistics_final_file_results)
      selected_path <- as.character(parsed$datapath[1])
      
      if (!is.na(selected_path) && selected_path != "") {
        file_path_statistics_final_file_results(selected_path)
      }
    }
  })
  
  output$file_path_statistics_final_file_results_output <- renderText({
    req(file_path_statistics_final_file_results())
    file_path_statistics_final_file_results()
  })
  
  
  # Import results when button is clicked
  values_list_stats <- reactiveValues()
  
  observeEvent(input$import_results_btn, {
    req(file_path_sources_final_file_results(),
        file_path_statistics_file_results(),
        file_path_statistics_final_file_results())
    
    # Read input files
    sources <- read.table(file_path_sources_final_file_results(), header = TRUE)
    statisticsCooling <- read.table(file_path_statistics_file_results(), header = TRUE)
    statisticsDetection <- read.table(file_path_statistics_final_file_results(), header = FALSE)
    
    # Set column names to match
    colnames(statisticsDetection) <- colnames(statisticsCooling)
    
    # Combine and split by plan
    statisticsByPlan <- rbind(statisticsCooling, statisticsDetection)
    statisticsByPlan <- split(
      statisticsByPlan,
      factor(paste0("(", statisticsByPlan[,5], ";", statisticsByPlan[,6], ")"))
    )
    
    tableNumber=table(statisticsByPlan[[1]][,4])
    x_breaks <- seq(0, 1, length.out = 50)
    
    values_list_stats[["meanstat4"]] <- mean(statisticsByPlan[[1]][,4])
    values_list_stats[["medianstat4"]] <- median(statisticsByPlan[[1]][,4])
    values_list_stats[["modestat4"]] <- as.numeric(names(tableNumber)[which.max(tableNumber)])
   
    
    # Clustering sources
    clustering <- hclust(dist(sources), method = "ward.D2")
    values_list_stats[["clustering"]] <- clustering
    values_list_stats[["clusters"]] <- cutree(clustering, k = 20)
    values_list_stats[["cut_height"]] <- sort(clustering$height, decreasing = TRUE)[20]
    values_list_stats[["x_breaks"]] <- x_breaks
    
    # Compute the 2d density of the sources
    ndim=ncol(sources)
    numberPatternDetection=nrow(statisticsDetection)/length(statisticsByPlan)
    plan = 0
    for(dim1 in 1:(ndim-1))
    {
      for (dim2 in (dim1+1):ndim) {
        plan = plan+1
        x_bins <- cut(sources[, dim1], breaks = x_breaks, include.lowest = TRUE)
        y_bins <- cut(sources[, dim2], breaks = x_breaks, include.lowest = TRUE)
        
        # # Create a 2D table of counts
        grid_counts <- table(x_bins, y_bins)
        
        # # Convert table to matrix
        z <- as.matrix(grid_counts)
        
        # # Transpose to align with x/y axes
        z <- z/numberPatternDetection
        
        values_list_stats[[paste0("density2d_",plan)]] <- z
        
      }
    }
    values_list_stats[["numberPlanesTotal"]] <- plan
    
    # Store results in reactive state
    reactive_state_results$numberPatternDetection<-numberPatternDetection
    reactive_state_results$sourcesDetection <- sources
    reactive_state_results$statisticsByPlan <- statisticsByPlan
    reactive_state_results$numberPlanes <- length(statisticsByPlan)
    reactive_state_results$planes <- names(statisticsByPlan)
    reactive_state_results$numberPatternCooling <- nrow(statisticsCooling)/length(statisticsByPlan)
    
    
    
    
  })
  
  # Generate dynamic UI for statistics plots
  output$plot_stats_results <- renderUI({
    if (is.null(reactive_state_results$numberPlanes) || reactive_state_results$numberPlanes < 1) return(NULL)
    
    plot_uis_stats <- lapply(1:reactive_state_results$numberPlanes, function(i) {
      ns <- paste0("plot_stats_", i)
      wellPanel(
        h4(paste("Statistics for plane ", i)),
        downloadButton(paste0("save_plot_stats_", i), "Save Plot as .png"),
        plotOutput(paste0("plot_stats_", i))
      )
    })
    
    do.call(tagList, plot_uis_stats)
  })
  
  
  
  
  observe({
    if (is.null(reactive_state_results$numberPlanes) || reactive_state_results$numberPlanes < 1) return()
    
    for (i in 1:reactive_state_results$numberPlanes) {
      local({
        ii <- i
        plot_output <- paste0("plot_stats_", ii)
        save_plot <- paste0("save_plot_stats_", ii)
        
        statsplan=reactive_state_results$statisticsByPlan[[ii]]
        
        
        
        
        
        # Render plot
        output[[plot_output]] <- renderPlot({
          par(mfrow=c(1,2))
          par(mar=c(5,6,4,1)+.1)
          plot(statsplan[,1],type = "l",main = "",ylab = "n_e",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
          abline(v=reactive_state_results$numberPatternCooling,col="red",lwd=2,lty=2)
          plot(statsplan[,2],type = "l",main = "",ylab = "n_a",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
          abline(v=reactive_state_results$numberPatternCooling,col="red",lwd=2,lty=2)
        })
        
        
        output[[paste0("save_plot_stats_", ii)]] <- downloadHandler(
          filename = function() {
            paste0("plot_stats_", ii, ".png")
          },
          content = function(file) {
            png(file,width = 1000,height = 200)
            par(mfrow=c(1,2))
            par(mar=c(5,6,4,1)+.1)
            plot(statsplan[,1],type = "l",main = "",ylab = "n_e",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
            abline(v=reactive_state_results$numberPatternCooling,col="red",lwd=2,lty=2)
            plot(statsplan[,2],type = "l",main = "",ylab = "n_a",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
            abline(v=reactive_state_results$numberPatternCooling,col="red",lwd=2,lty=2)
            dev.off()
          }
        )
      })
    }
  })
  
  # Generate dynamic UI for global statistics plots
  output$plot_stats_global_results <- renderUI({
    if (is.null(reactive_state_results$numberPlanes) || reactive_state_results$numberPlanes < 1) return(NULL)
    
    wellPanel(
      h4("Plot for the Strauss statistic"),
      downloadButton("save_plot_strauss_stat", "Save Plot as .png"),
      plotOutput("plot_strauss_stat"),
      h4("Plot for the number of sources per pattern statistic"),
      downloadButton("save_plot_number_stat", "Save Plot as .png"),
      plotOutput("plot_number_stat"),
      h4("Dendrogramm of the sources used for the detection"),
      downloadButton("save_plot_dendrogramm", "Save Plot as .png"),
      plotOutput("plot_dendrogramm")
    )
    
    
  })
  
  
  
  # Render strauss plot
  output$plot_strauss_stat <- renderPlot({
    par(mfrow=c(1,1))
    par(mar=c(5,6,4,1)+.1)
    plot(reactive_state_results$statisticsByPlan[[1]][,3],type = "l",main = "",ylab = "n_r",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
    abline(v=reactive_state_results$numberPatternCooling,col="red",lwd=2,lty=2)
  })
  
  
  output$save_plot_strauss_stat <- downloadHandler(
    filename = function() {
      paste0("plot_stats_", 3, ".png")
    },
    content = function(file) {
      png(file,width = 1000,height = 200)
      par(mfrow=c(1,1))
      par(mar=c(5,6,4,1)+.1)
      plot(reactive_state_results$statisticsByPlan[[1]][,3],type = "l",main = "",ylab = "n_r",col=palette.colors()[6],cex.lab=2.5,cex.axis=2)
      abline(v=reactive_state_results$numberPatternCooling,col="red",lwd=2,lty=2)
      dev.off()
    }
  )
  
  
  # Render number plot
  output$plot_number_stat <- renderPlot({
    stat_values <- reactive_state_results$statisticsByPlan[[1]][,4]
    
    par(mfrow = c(1, 2))
    par(mar = c(5, 6, 4, 1) + 0.1)
    
    # Left Plot: Line Plot with Annotations
    plot(stat_values,
         type = "l",
         main = "",
         ylab = "n",
         col = palette.colors()[6],
         cex.lab = 2.5,
         cex.axis = 2)
    
    abline(h = values_list_stats$meanstat4, col = palette.colors()[8], lwd = 2)
    abline(h = values_list_stats$modestat4, col = palette.colors()[3], lwd = 2)
    abline(v = reactive_state_results$numberPatternCooling, col = "red", lwd = 2, lty = 2)
    abline(h = values_list_stats$medianstat4, col = palette.colors()[2], lwd = 2)
    
    # Right Plot: Histogram or Barplot based on unique levels
    
    if (length(levels(as.factor(stat_values))) == 1) {
      hist(table(stat_values) / nrow(stat_values),
           col = palette.colors()[6],
           xlab = "n",
           main = "",
           ylab = "Probability",
           cex.lab = 2,
           cex.axis = 2)
      
      legend("topright",
             legend = paste0(c("median = ", "mean = ", "mode = "),
                             round(c(values_list_stats$medianstat4,
                                     values_list_stats$meanstat4,
                                     values_list_stats$modestat4), 2)),
             cex = 1.5)
      
    } else {
      hist(stat_values,
           col = palette.colors()[6],
           prob = TRUE,
           nclass = 20,
           xlab = "n",
           main = "",
           cex.lab = 2,
           cex.axis = 2)
      
      abline(v = values_list_stats$meanstat4, col = palette.colors()[8], lwd = 2)
      abline(v = values_list_stats$modestat4, col = palette.colors()[3], lwd = 2)
      abline(v = values_list_stats$medianstat4, col = palette.colors()[2], lwd = 2)
      
      legend("topright",
             legend = paste0(c("median = ", "mean = ", "mode = "),
                             round(c(values_list_stats$medianstat4,
                                     values_list_stats$meanstat4,
                                     values_list_stats$modestat4), 2)),
             cex = 1.5)
    }
    
    
  })
  
  
  output$save_plot_number_stat <- downloadHandler(
    
    filename = function() {
      paste0("plot_stats_", 4, ".png")
    },
    content = function(file) {
      stat_values <- reactive_state_results$statisticsByPlan[[1]][,4]
      png(file,width = 1000,height = 300)
      
      par(mfrow = c(1, 2))
      par(mar = c(5, 6, 4, 1) + 0.1)
      
      # Left Plot: Line Plot with Annotations
      plot(stat_values,
           type = "l",
           main = "",
           ylab = "n",
           col = palette.colors()[6],
           cex.lab = 2.5,
           cex.axis = 2)
      
      abline(h = values_list_stats$meanstat4, col = palette.colors()[8], lwd = 2)
      abline(h = values_list_stats$modestat4, col = palette.colors()[3], lwd = 2)
      abline(v = reactive_state_results$numberPatternCooling, col = "red", lwd = 2, lty = 2)
      abline(h = values_list_stats$medianstat4, col = palette.colors()[2], lwd = 2)
      
      # Right Plot: Histogram or Barplot based on unique levels
      
      if (length(levels(as.factor(stat_values))) == 1) {
        hist(table(stat_values) / nrow(stat_values),
             col = palette.colors()[6],
             xlab = "n",
             main = "",
             ylab = "Probability",
             cex.lab = 2,
             cex.axis = 2)
        
        legend("topright",
               legend = paste0(c("median = ", "mean = ", "mode = "),
                               round(c(values_list_stats$medianstat4,
                                       values_list_stats$meanstat4,
                                       values_list_stats$modestat4), 2)),
               cex = 1.5)
        
      } else {
        hist(stat_values,
             col = palette.colors()[6],
             prob = TRUE,
             nclass = 20,
             xlab = "n",
             main = "",
             cex.lab = 2,
             cex.axis = 2)
        
        abline(v = values_list_stats$meanstat4, col = palette.colors()[8], lwd = 2)
        abline(v = values_list_stats$modestat4, col = palette.colors()[3], lwd = 2)
        abline(v = values_list_stats$medianstat4, col = palette.colors()[2], lwd = 2)
        
        legend("topright",
               legend = paste0(c("median = ", "mean = ", "mode = "),
                               round(c(values_list_stats$medianstat4,
                                       values_list_stats$meanstat4,
                                       values_list_stats$modestat4), 2)),
               cex = 1.5)
      }
      
      
      dev.off()
    }
  )
  # Render dendrogramm
  output$plot_dendrogramm <- renderPlot({
    
    plot(cut(as.dendrogram(values_list_stats$clustering),h=values_list_stats$cut_height)$upper,ylim=c(values_list_stats$cut_height,max(values_list_stats$clustering$height)),main="Dendrogram with a maximum of 20 clusters", leaflab = "none")
    
  })
  
  
  output$save_plot_dendrogramm <- downloadHandler(
    filename = function() {
      paste0("plot_dendrogramm_",".png")
    },
    content = function(file) {
      
      png(file,width = 1000,height = 1000)
      plot(cut(as.dendrogram(values_list_stats$clustering),h=values_list_stats$cut_height)$upper,ylim=c(values_list_stats$cut_height,max(values_list_stats$clustering$height)),main="Dendrogram with a maximum of 20 clusters", leaflab = "none")
      dev.off()
    }
  )
  
  # Reverse normalise data when button is clicked
  reverse_normalised_result <- eventReactive(input$reverse_normalise_btn, {
    req(input$data_file_results)
    df <- reactive_state_results$normalisedData
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    reverse_norm_df <- df
    reverseSourcesDetection <- reactive_state_results$sourcesDetection
    if(!is.null(reactive_state_results$centers))
    {
      reverseCenters <- reactive_state_results$centers
    }else{
      reverseCenters <- NULL
    }
    if(!is.null(reactive_state_results$proposedSources))
    {
      reverseProposedSources <- reactive_state_results$proposedSources
    }else{
      reverseProposedSources <- NULL
    }
    
    
    for (col in numeric_cols) {
      indice <- which(numeric_cols==col)[1]
      val1 <- input[[paste0(col, "_val1")]]
      val2 <- input[[paste0(col, "_val2")]]
      if (!is.null(val1) && !is.null(val2)) {
        reverse_norm_df[[col]] <- df[[col]] *(val2 - val1) + val1
        
        if (!is.null(reverseSourcesDetection[[col]])) {
          reverseSourcesDetection[[col]] <- reverseSourcesDetection[[col]] *(val2 - val1) + val1
        }
        if (!is.null(reverseProposedSources[[col]])) {
          reverseProposedSources[[col]] <- reverseProposedSources[[col]] *(val2 - val1) + val1
        }
        if(!is.null(reverseCenters)){
          reverseCenters[reverseCenters[,3]==indice,1] <- reverseCenters[reverseCenters[,3]==indice,1]*(val2 - val1) + val1
          reverseCenters[reverseCenters[,4]==indice,2] <- reverseCenters[reverseCenters[,4]==indice,2]*(val2 - val1) + val1
        }
        
      } else {
        reverse_norm_df[[col]] <- NA  # Invalid bounds
        
        if (!is.null(reverseSourcesDetection[[col]])) {
          reverseSourcesDetection[[col]] <- NA  # Invalid bounds
        }
        if (!is.null(reverseProposedSources[[col]])) {
          reverseProposedSources[[col]] <- NA  # Invalid bounds
        }
        if(!is.null(reverseCenters)){
          reverseCenters[reverseCenters[,3]==indice,1] <- NA
          reverseCenters[reverseCenters[,4]==indice,2] <- NA
        }
        
      }
      
      
    }
    
    # Update reactive values instead of global ones
    reactive_state_results$reverseNormalisedData <- reverse_norm_df
    reactive_state_results$reverseSourcesDetection <- reverseSourcesDetection
    reactive_state_results$reverseProposedSources <- reverseProposedSources
    reactive_state_results$reverseCenters <- reverseCenters
    
    numberColumn=length(numeric_cols)
    
    numberPatternDetection = reactive_state_results$numberPatternDetection
    plan = 0
    for(indice1 in 1:(numberColumn-1))
    {
      val1dim1 <- input[[paste0(numeric_cols[indice1], "_val1")]]
      val2dim1 <- input[[paste0(numeric_cols[indice1], "_val2")]]
      
      if (length(val1dim1) == 1 && length(val2dim1) == 1 && 
          !is.null(val1dim1) && !is.null(val2dim1) && 
          !is.na(val1dim1) && !is.na(val2dim1)) {
        x_breaks <- seq(val1dim1, val2dim1, length.out = 50)
        
        
        for (indice2 in (indice1+1):numberColumn) {
          
          plan = plan+1
          
          val1dim2 <- input[[paste0(numeric_cols[indice2], "_val1")]]
          val2dim2 <- input[[paste0(numeric_cols[indice2], "_val2")]]
          
          if (length(val1dim2) == 1 && length(val2dim2) == 1 && 
              !is.null(val1dim2) && !is.null(val2dim2) && 
              !is.na(val1dim2) && !is.na(val2dim2)) {
            y_breaks <- seq(val1dim2, val2dim2, length.out = 50)
            
            x_bins <- cut(reverseSourcesDetection[, indice1], breaks = x_breaks, include.lowest = TRUE)
            y_bins <- cut(reverseSourcesDetection[, indice2], breaks = y_breaks, include.lowest = TRUE)
            
            # # Create a 2D table of counts
            grid_counts <- table(x_bins, y_bins)
            
            # # Convert table to matrix
            z <- as.matrix(grid_counts)
            
            # # Transpose to align with x/y axes
            z <- z/numberPatternDetection
            
            values_list_stats[[paste0("density2d_reverse_",plan)]] <- z
          }
        }
      }
    }
    
    reverse_norm_df
  })
  
  # Render reverse normalised data
  output$filePreview_reverse_normalised_data <- renderTable({
    req(reverse_normalised_result())
    if(is.null(reactive_state_results$reverseNormalisedData)) return(NULL)
    round(reactive_state_results$reverseNormalisedData[1:min(c(5,nrow(reactive_state_results$reverseNormalisedData))),],5)
    
  })
  
  # Select the number of sources proposed
  output$number_sources_detected_slider <- renderUI({
    req(values_list_stats[["modestat4"]])
    if (is.null(values_list_stats[["modestat4"]]) ) return(NULL)
    
    ui_elements <- list(
      sliderInput(
        inputId = "number_sources_detected",
        label = "Number of sources detected",
        min = 0,
        max = 20,
        value = values_list_stats[["modestat4"]],
        step = 1
      ),
      br(), br(),
      fluidRow(
        column(6,
               h3("Normalised Proposed Sources :"),
               tableOutput("proposedSource"),
               downloadButton("save_proposed_source", "Saved Normalised Proposed Sources as .txt")
        ),
        if (!is.null(input$show_file_input) && input$show_file_input) {
          column(6,
                 h3("Reverse Normalised Proposed Sources :"),
                 tableOutput("reverseProposedSource"),
                 downloadButton("save_reverse_proposed_source", "Saved Reverse Normalised Proposed Sources as .txt")
          )
        }
      )
    )
    
    tagList(ui_elements)
    
    
    
    
  })
  
  output$save_proposed_source <- downloadHandler(
    filename = function() {
      paste0("proposed_source", ".txt")
    },
    content = function(file) {
      # write.table( combined_points(), file, row.names=F, col.names=T, append=F )
      write.table(reactive_state_results$proposedSources[,-1], file, sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE, append=F)
    }
  )
  output$save_reverse_proposed_source <- downloadHandler(
    filename = function() {
      paste0("reverse_proposed_source", ".txt")
    },
    content = function(file) {
      # write.table( combined_points(), file, row.names=F, col.names=T, append=F )
      write.table(reactive_state_results$reverseProposedSources[,-1], file, sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE, append=F)
    }
  )
  
  proposed_sources <- reactive({
    req(input$number_sources_detected)
    
    num_sources <- input$number_sources_detected
    
    # Make sure you have access to your data
    if (is.null(reactive_state_results$sourcesDetection)) return(NULL)
    
    sources <- reactive_state_results$sourcesDetection
    clusters <- cutree(values_list_stats[["clustering"]],k=num_sources)
    
    pSources <- aggregate(sources, by = list(Num = clusters), FUN = mean)
    reactive_state_results$proposedSources <- pSources
    
    if(!is.null(reactive_state_results$reverseProposedSources))
    {
      update_reverse_proposed_sources()
    }
    
    pSources
  })
  
  output$proposedSource <- renderTable({
    reactive_state_results$proposedSources
  })
  
  update_reverse_proposed_sources <-reactive({
    req(input$data_file_results)
    df <- reactive_state_results$normalisedData
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    reverseSourcesDetection <- reactive_state_results$sourcesDetection
    if(!is.null(reactive_state_results$proposedSources))
    {
      reverseProposedSources <- reactive_state_results$proposedSources
    }else{
      reverseProposedSources <- NULL
    }
    
    for (col in numeric_cols) {
      val1 <- input[[paste0(col, "_val1")]]
      val2 <- input[[paste0(col, "_val2")]]
      if (!is.null(val1) && !is.null(val2)) {
        if (!is.null(reverseProposedSources[[col]])) {
          reverseProposedSources[col] <- reverseProposedSources[[col]] *(val2 - val1) + val1
        }
        
      } else {
        if (!is.null(reverseProposedSources[[col]])) {
          reverseProposedSources[[col]] <- NA  # Invalid bounds
        }
        
        
      }
    }
    
    
    # Update reactive values instead of global ones
    reactive_state_results$reverseProposedSources <- reverseProposedSources
  })
  
  output$reverseProposedSource <- renderTable({
    reactive_state_results$reverseProposedSources
  })
  
  # Generate dynamic UI for 2d plots
  output$plot_results_2d <- renderUI({
    if (is.null(values_list_stats$numberPlanesTotal) || values_list_stats$numberPlanesTotal < 1) return(NULL)
    
    plot_uis_2d <- lapply(1:values_list_stats$numberPlanesTotal, function(i) {
      wellPanel(
        h4(paste("Projection on normalised plane ", i)),
        downloadButton(paste0("save_plot_2d_", i), "Save Plot as .png"),
        plotOutput(paste0("plot_2d_", i))
      )
    })
    
    do.call(tagList, plot_uis_2d)
  })
  
  observeEvent(input$number_sources_detected,{
    req(reactive_state_results$sourcesDetection)
    if (is.null(values_list_stats$numberPlanesTotal) || values_list_stats$numberPlanesTotal < 1) return()
    
    plan <- 0
    for (dim1 in 1:(ncol(reactive_state_results$sourcesDetection) - 1)) {
      for (dim2 in (dim1 + 1):ncol(reactive_state_results$sourcesDetection)) {
        plan <- plan + 1
        
        local({
          ii <- plan
          d1 <- dim1
          d2 <- dim2
          
          plot_id <- paste0("plot_2d_", ii)
          save_id <- paste0("save_plot_2d_", ii)
          density_data <- values_list_stats[[paste0("density2d_", ii)]]
          
          output[[plot_id]] <- renderPlot({
            req(values_list_stats$x_breaks, density_data)
            par(mar = c(5, 6, 4, 1) + 0.1)
            
            image.plot(
              x = values_list_stats$x_breaks, y = values_list_stats$x_breaks, z = density_data,  
              col = colfunc,
              cex.axis = 3, cex.lab = 3,
              xlab = colnames(reactive_state_results$sourcesDetection)[[d1]],
              ylab = colnames(reactive_state_results$sourcesDetection)[[d2]]
            )
            
            abline(h = seq(0, 1, length.out = 50), col = "gray")
            abline(v = seq(0, 1, length.out = 50), col = "gray")
            if(!is.null(reactive_state_results$normalisedData))
            {
              points(reactive_state_results$normalisedData[, c(d1, d2)], pch = 19, col = "black", cex = 1)
            }
            if(!is.null(proposed_sources()))
            {
              points(proposed_sources()[, c(d1+1, d2+1)], pch = 19, col = 'darkgreen', cex = 3)
            }
            
            if(!is.null(reactive_state_results$centers))
            {
              points(reactive_state_results$centers[(reactive_state_results$centers[,3]==d1) & (reactive_state_results$centers[,4]==d2), c(1,2)], pch = 17, col = "red", cex = 3)
            }
          })
          
          output[[save_id]] <- downloadHandler(
            filename = function() {
              paste0("plot_2d_plane_", ii, ".png")
            },
            content = function(file) {
              png(file, width = 800, height = 600)
              par(mar = c(5, 6, 4, 1) + 0.1)
              
              image.plot(
                x = values_list_stats$x_breaks, y = values_list_stats$x_breaks, z = density_data,  
                col = colfunc,
                cex.axis = 3, cex.lab = 3,
                xlab = colnames(reactive_state_results$sourcesDetection)[[d1]],
                ylab = colnames(reactive_state_results$sourcesDetection)[[d2]]
              )
              
              abline(h = seq(0, 1, length.out = 50), col = "gray")
              abline(v = seq(0, 1, length.out = 50), col = "gray")
              
              if(!is.null(reactive_state_results$normalisedData))
              {
                points(reactive_state_results$normalisedData[, c(d1, d2)], pch = 19, col = "black", cex = 1)
              }
              if(!is.null(proposed_sources()))
              {
                points(proposed_sources()[, c(d1+1, d2+1)], pch = 19, col = 'darkgreen', cex = 3)
              }
              if(!is.null(reactive_state_results$centers))
              {
                points(reactive_state_results$centers[(reactive_state_results$centers[,3]==d1) &(reactive_state_results$centers[,4]==d2), c(1,2)], pch = 17, col = "red", cex = 3)
              }
              dev.off()
            }
          )
        }) # local
      }
    }
  })
  
  # Generate dynamic UI for 2d plots after the reverse normalisation
  
  
  output$plot_reverse_results_2d <- renderUI({
    if (is.null(values_list_stats$numberPlanesTotal) || values_list_stats$numberPlanesTotal < 1) return(NULL)
    if(is.null(reactive_state_results$reverseProposedSources)) return(NULL)
    plot_uis_reverse_2d <- lapply(1:values_list_stats$numberPlanesTotal, function(i) {
      wellPanel(
        h4(paste("Projection on plane ", i)),
        downloadButton(paste0("save_plot_reverse_2d_", i), "Save Plot as .png"),
        plotOutput(paste0("plot_reverse_2d_", i))
      )
    })
    
    do.call(tagList, plot_uis_reverse_2d)
  })
  
  trigger_redraw <- reactiveVal(0)
  
  observeEvent(input$reverse_normalise_btn, {
    trigger_redraw(trigger_redraw() + 1)
  })
  
  observeEvent(input$number_sources_detected, {
    trigger_redraw(trigger_redraw() + 1)
  })
  
  observeEvent(trigger_redraw(),{
    req(reactive_state_results$reverseSourcesDetection)
    if (is.null(values_list_stats$numberPlanesTotal) || values_list_stats$numberPlanesTotal < 1) return()
    if(is.null(reactive_state_results$reverseProposedSources)) return(NULL)
    
    
    plan <- 0
    for (dim1 in 1:(ncol(reactive_state_results$reverseSourcesDetection) - 1)) {
      for (dim2 in (dim1 + 1):ncol(reactive_state_results$reverseSourcesDetection)) {
        plan <- plan + 1
        
        local({
          colnamesreverseSourcesDetection = colnames(reactive_state_results$reverseSourcesDetection)
          ii <- plan
          d1 <- dim1
          d2 <- dim2
          val1dim1 <- input[[paste0(colnamesreverseSourcesDetection[dim1], "_val1")]]
          val2dim1 <- input[[paste0(colnamesreverseSourcesDetection[dim1], "_val2")]]
          x_breaks <- seq(val1dim1, val2dim1, length.out = 50)
          val1dim2 <- input[[paste0(colnamesreverseSourcesDetection[dim2], "_val1")]]
          val2dim2 <- input[[paste0(colnamesreverseSourcesDetection[dim2], "_val2")]]
          y_breaks <- seq(val1dim2, val2dim2, length.out = 50)
          
          
          plot_id <- paste0("plot_reverse_2d_", ii)
          save_id <- paste0("save_plot_reverse_2d_", ii)
          density_data <- values_list_stats[[paste0("density2d_reverse_", ii)]]
          
          output[[plot_id]] <- renderPlot({
            req(values_list_stats$x_breaks, density_data)
            par(mar = c(5, 6, 4, 1) + 0.1)
            
            image.plot(
              x = x_breaks, y = y_breaks, z = density_data,  
              col = colfunc,
              cex.axis = 3, cex.lab = 3,
              xlab = colnames(reactive_state_results$reverseSourcesDetection)[[d1]],
              ylab = colnames(reactive_state_results$reverseSourcesDetection)[[d2]]
            )
            
            abline(v = x_breaks, col = "gray")
            abline(h = y_breaks, col = "gray")
            if(!is.null(reactive_state_results$rawData))
            {
              points(reactive_state_results$rawData[, c(d1, d2)], pch = 19, col = "black", cex = 1)
            }
            if(!is.null(reactive_state_results$reverseProposedSources))
            {
              points(reactive_state_results$reverseProposedSources[, c(d1+1, d2+1)], pch = 19, col = 'darkgreen', cex = 3)
            }
            
            if(!is.null(reactive_state_results$reverseCenters))
            {
              points(reactive_state_results$reverseCenters[(reactive_state_results$reverseCenters[,3]==d1) & (reactive_state_results$reverseCenters[,4]==d2), c(1,2)], pch = 17, col = "red", cex = 3)
            }
          })
          
          output[[save_id]] <- downloadHandler(
            filename = function() {
              paste0("plot_2d_plane_", ii, ".png")
            },
            content = function(file) {
              png(file, width = 800, height = 600)
              par(mar = c(5, 6, 4, 1) + 0.1)
              
              image.plot(
                x = x_breaks, y = y_breaks, z = density_data,  
                col = colfunc,
                cex.axis = 3, cex.lab = 3,
                xlab = colnames(reactive_state_results$reverseSourcesDetection)[[d1]],
                ylab = colnames(reactive_state_results$reverseSourcesDetection)[[d2]]
              )
              
              abline(v = x_breaks, col = "gray")
              abline(h = y_breaks, col = "gray")
              
              if(!is.null(reactive_state_results$rawData))
              {
                points(reactive_state_results$rawData[, c(d1, d2)], pch = 19, col = "black", cex = 1)
              }
              if(!is.null(reactive_state_results$reverseProposedSources))
              {
                points(reactive_state_results$reverseProposedSources[, c(d1+1, d2+1)], pch = 19, col = 'darkgreen', cex = 3)
              }
              if(!is.null(reactive_state_results$reverseCenters))
              {
                points(reactive_state_results$reverseCenters[(reactive_state_results$reverseCenters[,3]==d1) &(reactive_state_results$reverseCenters[,4]==d2), c(1,2)], pch = 17, col = "red", cex = 3)
              }
              dev.off()
            }
          )
        }) # local
      }
    }
  })
  
  
  
}

shinyApp(ui, server)







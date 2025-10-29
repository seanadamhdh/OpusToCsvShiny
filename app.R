

if(!require(shiny)){
  install.packages("shiny")
  require(shiny)
}
if(!require(shinyjs)){
  install.packages("shinyjs")
  reqiure(shinyjs)
}
if(!require(tidyverse)){
  install.packages("tidyverse")
  require(tidyverse)
}
if(!require(resemble)){
  install.packages("resemble")
  require(resemble)
}
if(!require(simplerspec)){
  devtools::install_github("https://github.com/philipp-baumann/simplerspec.git")
  require(simplerspec)}
if(!require(prospectr)){
  #install_packages("prospectr")
  devtools::install_github("https://github.com/l-ramirez-lopez/prospectr.git")
  require(prospectr)
}




OPUSraw_to_Preprocessed<-function(
    paths,
    sample_id,
    from=400,
    to=7500,
    id_string_location
){
  # read raw opus data
  raw_scans<-read_opus_univ(paths)
  
  
  # combine in readable dataset
  raw_spc<-gather_spc(raw_scans)
  
  
  # shortening some identifer col
  
  raw_spc$sample_id<-substr(sample_id,id_string_location[1],id_string_location[2])
  
  # formatting 1: from simplerspec to full unnested
  raw_spc%>%
    unnest(spc)->unnested_spc
  
  # formatting 2: prospectr friendly
  # averaging repetition spectra
  spc<-tibble(unnested_spc[c(3, # sample_id column
                             5:13864 #!!! spc colums... depends on your scan resolution and range
  )])%>%
    group_by(sample_id)%>%
    summarise_all(mean) #simple averageing of replicate scans
  
  spc0<-spc[2:13861] # !!! pulling spc columns
  spc_data<-tibble(spc[1],spc=spc0) #and reintroduicing them as a nested matix (prospectr-friendly)
  
  
  # metadata, ! using always metadata of first repetiton
  spc_data$metadata<-unnested_spc%>%group_by(sample_id)%>%summarise(metadata=first(metadata))
  
  
  # resampling 7500-400 @ 1 cm-1 resolution
  spc_data$spc_rs<-resample(spc_data$spc,
                            wav=as.numeric(names(spc_data$spc)),
                            new.wav=seq(7500, #start !!!
                                        400,  #stop  !!!
                                        -1    #step  !!!
                            ))
  
  spc_data$spc_sg <- savitzkyGolay(spc_data$spc_rs, m = 0, 
                                   p = 3, w = 21)
  spc_data$spc_sg_snv <- standardNormalVariate(spc_data$spc_sg)
  if (nrow(spc_data$spc_sg) > 1) {
    spc_data$spc_sg_bl <- baseline(spc_data$spc_sg, wav = as.numeric(colnames(spc_data$spc_sg)))
  }
  else {
    spc_data$spc_sg_bl <- baseline(spc_data$spc_sg, wav = as.numeric(colnames(spc_data$spc_sg))) %>% 
      as.list() %>% as_tibble %>% as.matrix()
  }
  spc_data$spc_sg_bl <- spc_data$spc_sg_bl[, -c(1, ncol(spc_data$spc_sg)), 
                                           drop = F]
  spc_data$spc_sg_bl_rs4 <- resample(spc_data$spc_sg_bl, wav = as.numeric(colnames(spc_data$spc_sg_bl)), 
                                     new.wav = seq(to - 14, from + 14, -4))
  spc_data$spc_sg_snv_rs4 <- resample(spc_data$spc_sg_snv, 
                                      wav = as.numeric(colnames(spc_data$spc_sg_snv)), 
                                      new.wav = seq(colnames(spc_data$spc_sg_snv) %>% 
                                                 as.numeric %>% max, colnames(spc_data$spc_sg_snv) %>% 
                                                 as.numeric %>% min, -4))
  spc_data$spc_rs4 = resample(spc_data$spc_rs, wav = as.numeric(colnames(spc_data$spc_rs)), 
                              new.wav = seq(max(as.numeric(colnames(spc_data$spc_rs))), 
                                            min(as.numeric(colnames(spc_data$spc_rs))), -4))
  spc_data$spc_sg_rs4 = resample(spc_data$spc_sg, wav
                                 = as.numeric(colnames(spc_data$spc_sg)), 
                                 new.wav = seq(max(as.numeric(colnames(spc_data$spc_sg))), 
                                               min(as.numeric(colnames(spc_data$spc_sg))), -4))
  spc_data$spc_sg1d = savitzkyGolay(spc_data$spc_rs, m = 1, 
                                    p = 3, w = 41)
  spc_data$spc_sg1d_rs4 = resample(spc_data$spc_sg1d, wav = as.numeric(colnames(spc_data$spc_sg1d)), 
                                   new.wav = seq(max(as.numeric(colnames(spc_data$spc_sg1d))), 
                                                 min(as.numeric(colnames(spc_data$spc_sg1d))), -4))
  spc_data$spc_sg2d = savitzkyGolay(spc_data$spc_rs, m = 2, 
                                    p = 3, w = 41)
  spc_data$spc_sg2d_rs4 = resample(spc_data$spc_sg2d, wav = as.numeric(colnames(spc_data$spc_sg1d)), 
                                   new.wav = seq(max(as.numeric(colnames(spc_data$spc_sg1d))), 
                                                 min(as.numeric(colnames(spc_data$spc_sg1d))), -4))
  return(spc_data)
}






ui <- fluidPage(
  
  # Application title
  titlePanel("Convert OPUS .0 to .csv"),
  
  
  sidebarLayout(
    
    
    sidebarPanel(
      
    selectInput(inputId = "set",
                  label="Preprocessing",
                  choices = list(
                    "raw"="spc",
                    "resampled 1 cm-1"="spc_rs",
                    "Savitzky-Golay"="spc_sg",
                     "Baseline-corrected"="spc_sg_bl",
                     "Baseline-corrected 4 cm-1"="spc_sg_bl_rs4",
                    # baseline corrected data produces this warning in mbl():
                    #' Warning in fit_and_predict(x = i_k_xr, y = i_k_yr, pred_method = method$method,  :
                    #' Variables with zero variance. Data will not be scaled
                    
                    "StandardNormalVariate"="spc_sg_snv",
                    "StandardNormalVariate 4 cm-1"="spc_sg_snv_rs4",
                    "First derivative"="spc_sg1d",
                    "Second derivative"="spc_sg2d",
                    "First derivative 4 cm-1"="spc_sg1d_rs4",
                    "Second derivative 4 cm-1"="spc_sg1d_rs4"
                  ),
                  selected ="spc_sg_snv_rs4",
                  multiple = F
      ),
      
      
    sliderInput(inputId = "wavenumber",
                label = "Wavenumber range",
                min = 400,
                max = 7500,
                step=1,
                value=c(400,7500)),
    
    
    sliderInput(inputId = "ID_string_location",
                label = "Sample identifier string position in filename",
                min = 1,
                max = 30,
                step=1,
                value=c(5,8)),
    
    fileInput(inputId = "new_spc",
              label = "Select OPUS files: Filenames must contain a consistantly placed unique sample idenifier 
              to aggregate extrenal replicates properly.
              Position can be selected above.",
              multiple = T
    ),
    
    withBusyIndicatorUI(
      actionButton("go","Run",
                   style="color: #fff; background-color: #22e",
                   class="btn-primary")
    ),

    
    uiOutput("download"),

      position="left"
    ), 
    
    
    mainPanel(
     #
      textOutput("action")

      
      
    )
  )
  
)

# SERVER ####
server <- function(input, output) {
  
  

    
    ## load and preprocess user spc ####
    load_spc<-eventReactive(input$go,{
      OPUSraw_to_Preprocessed(paths = input$new_spc$datapath,
                              sample_id = input$new_spc$name,id_string_location = input$ID_string_location,
                              from=input$wavenumber[1],
                              to=input$wavenumber[2])
    })
    
    #run MBL
    OUT<-eventReactive(input$go,{
      withBusyIndicatorServer("go",{
      withProgress(message="Reading spectra",detail = "This might take a while...",value=0.3,{
       #read user spc
      spectra<-load_spc()
      # init
      
      
      out<-list(spectra=tibble(spectra$sample_id,spectra[[input$set]]))
      setProgress(message="done",value=1)
      })
      out
    })
    })
    
  
    
   
    
    output$plot<-renderPlot({
      MBL<-MBL()
      plot(MBL$mod)
    })
    
    table<-reactive({
      OUT=OUT()
      OUT$spectra
    })
    ## main output ####
    output$table<-renderTable({
      
      table()
      
      # print(list("a"=1,"b"=2))
      #print(ncol(c(1:10))) 
      
    })
    
   output$downloadData <- downloadHandler(
      filename = function() {
        paste0(Sys.time(), ".csv")
      },
      content = function(file) {
        write.csv(table(),file, row.names = FALSE)
      })  
   
   
 output$download<-renderUI({
   req(input$go,table())
   downloadButton("downloadData","Download as csv")
   })    
 



 
 
 
 
 }

# Run the application 
shinyApp(ui = ui, server = server)




library(shiny)
library(shinydashboard)
library(RSQLite)
library(DBI)
library(factoextra)
library(ggfortify)
library("igraph")
library(corrplot)
library("AER") 
#install.packages("stargazer")
library(stargazer)
library("MASS")
library(ordinalNet)
library(readr)
library(readxl)
library(DT)
library(ggplot2)
library(gridExtra)
library(dplyr)
# We need the packages for Pie charts
if (!require("toxpiR"))install.packages("toxpiR")
# In case this simple installation does not function please try other options from here: https://github.com/ToxPi/toxpiR .
# There are at least 3 ways available to install this package.

# We need also the library
library(toxpiR)

#numbers are shown in full 0.000... not with "0.00...e"
options(scipen = 999)

shinyServer(function(input, output, session){
  
  #We create a small database
  helper_file <- "Database_Toxicology_Own.db"  
  
  #We open the connection
  toxi_con <- dbConnect(SQLite(), helper_file)  
  
  #We choose files
  observeEvent(input$press_here, {
  req(input$chemical_files)  
  
  #Files info
  for (i in seq_along(input$chemical_files$datapath)) {
  way <- input$chemical_files$datapath[i]
  endi <- tools::file_ext(input$chemical_files$name[i])  
      
  #Prepare to read file depending on extension
  if (endi %in% c("csv", "txt")) {
  dataset <- read_csv(way, show_col_types = FALSE)
  } 
  else if (endi %in% c("xls", "xlsx")) {
  dataset <- read_excel(way)
  }
  else {
  warning(paste("No supported:", endi, "Change file name and try again."))
  next
  }
      
      
  #Test table name
  tb <- gsub("[^a-zA-Z0-9_]", "_", as.character(tools::file_path_sans_ext(input$chemical_files$name[i])))
      
  #Try to print
  print(paste("Table name:", tb))
      
  #Checks
  if (nzchar(tb)) {
    
  #Check existence
  et <- dbListTables(toxi_con)
  if (tb %in% et) {
  print(paste("Table", tb, "is already in DB. We replace the old one with the new one."))
  } 
  else {
  print(paste("We add a new table:", tb))
  }
  
  #Improve names of columns
  colnames(dataset) <- make.names(colnames(dataset))
        
  #Show the improved names
  print(paste("Improved column names:", paste(colnames(dataset), collapse = ", ")))
        
  #Upload to DB and replace if needed
  dbWriteTable(toxi_con, tb, dataset, overwrite = TRUE)
  } 
  else {
  warning(paste("We cannot use this name:", tb, " Change name and try again."))
      }
  }
    
  #Message
  output$happening <- renderText("Files are now in our database!")
  })
  
  #connect to "Database_Toxicology_Own.db" and get chemical data
  funct_data_complete <- function(){
  #we try to manage errors
  chem_res <-tryCatch({
  toxi_con<- dbConnect(SQLite(), "Database_Toxicology_Own.db")
  dbListTables(toxi_con)
  #to delete:
  #dbExecute(toxi_con, "DROP TABLE README")
  #run query with some type of indexation
  dbExecute(toxi_con, "CREATE INDEX IF NOT EXISTS indexation_substance ON '1994_House_SupTab4' (Substance)")
  dbExecute(toxi_con, "CREATE INDEX IF NOT EXISTS indexation_p1_2_index ON '1994_House_SupTab4' (P1_2_index)")
  dbExecute(toxi_con, "CREATE INDEX IF NOT EXISTS indexation_p3_index ON '1994_House_SupTab4' (P3_index)")
  dbExecute(toxi_con, "CREATE INDEX IF NOT EXISTS indexation_p4_index ON '1994_House_SupTab4' (P4_index)")
  dbExecute(toxi_con, "CREATE INDEX IF NOT EXISTS indexation_p5_index ON '1994_House_SupTab4' (P5_index)")
  dbExecute(toxi_con, "CREATE INDEX IF NOT EXISTS indexation_p6_index ON '1994_House_SupTab4' (P6_index)")
  dbExecute(toxi_con, "CREATE INDEX IF NOT EXISTS indexation_p7_index ON '1994_House_SupTab4' (P7_index)")
  data_complete <- dbGetQuery(toxi_con, "select Substance, P1_2_index, P3_index, P4_index, P5_Index, P6_index, P7_index FROM  '1994_House_SupTab4'")
  #return message in case of need
  return(data_complete)}, chem_error = function(e){
    message("we identified a problem") 
    return(NULL)}
  )
  return(chem_res)
  }

  data_complete <- funct_data_complete()
  print(data_complete)
  #We are now preparing the data for K-Means
  
  data_complete.labels=data_complete$Substance
  table(data_complete.labels)
  
  data_complete_6lines <- data_complete[2:7]
  print(data_complete_6lines)
  
  # we will need to scale our data to prepare for applying the Euclidean distance
  data_complete_6lines_scale=scale(data_complete_6lines)
  print(data_complete_6lines_scale)
  
  data_complete_6lines <- dist(data_complete_6lines_scale)
  print(data_complete_6lines)
  
  # we need to calculate how many clusters we need by using a method called "elbow plot" 
  
  output$plotcluster <- renderPlot(fviz_nbclust(data_complete_6lines_scale, kmeans, method = "wss") + labs(subscribe = "Elbow method"))
  
  # now we apply K-Means
  
  output$plot1 <- renderPlot({
    km.out <- kmeans(data_complete_6lines_scale, centers = input$number_cluster, nstart=100)
    
    
  # show cluster plot
    
    km.clusters <- km.out$cluster
    rownames(data_complete_6lines_scale) <- data_complete$Substance
    
    fviz_cluster(list(data=data_complete_6lines_scale, cluster=km.clusters), show_labels=T, pointsize = 1.5, ellipse.type = "convex", geom = c("point", "text"), labelsize = 10, type = "phylogenic", ggtheme = theme_grey())
    
    
    
  })
  
  # Hierarchical clustering
  
  # Visualize dendrogram
  output$plot2 <- renderPlot({
    rownames(data_complete_6lines_scale) <- data_complete$Substance
    hc.cut <- hcut(data_complete_6lines_scale, k = 6,stand = T, hc_metric=input$distancemethod)
    fviz_dend(hc.cut, show_labels = TRUE, rect = TRUE, cex = 0.4, horiz = F, type = input$type, repel=T)
    
    # Visualize cluster
    #fviz_cluster(hc.cut, ellipse.type = "convex")
    
  })
  #PCA (Principal Component Analysis)
  
  output$plot3 <- renderPlot({ 
    #we open the connection to RSQLite
    toxi_con<- dbConnect(SQLite(), "Database_Toxicology_Own.db")
    # run our query
    data_complete_PCA <- dbGetQuery(toxi_con, "select Substance, Class, P3_7_index, P4_7_index, P5_7_index, P1_2_index, P1_index, P2_index, P3_index, P4_index, P5_Index, P6_index, P7_index FROM '1994_House_SupTab4'")
    #check output
    #print(data_complete_PCA)
    #prepare data for PCA
    data_frameworkPCA <- data_complete_PCA[,3:12]
    #check output
    #print(data_frameworkPCA)
    #scale data and apply function prcomp for PCA
    pca_results_1<- prcomp(data_frameworkPCA, scale. = TRUE)
    #calculate eigenvalues
    eigenvalues <- pca_results_1$sdev*pca_results_1$sdev
    eigenvalues
    #calculate eigenvectors
    eigenvectors <-pca_results_1$rotation
    eigenvectors
    #plot PCA in different ways
    #autoplot(pca_results_1)
    #autoplot(pca_results_1, data = data_complete_PCA , colour = 'Class', loadings = TRUE)
    autoplot(pca_results_1, data = data_complete_PCA, colour = 'Class',
             loadings = TRUE, loadings.colour = 'blue',
             loadings.label.repel = TRUE, loadings.label.size = 4, loadings.label.hjust=-1)
  })
  output$plot_5 <- renderPlot({ 
    #we open the connection to RSQLite
    toxi_con<- dbConnect(SQLite(), "Database_Toxicology_Own.db")
    
    # run our query
    data_complete_PCA <- dbGetQuery(toxi_con, "select Substance, Class, P3_7_index, P4_7_index, P5_7_index, P1_2_index, P1_index, P2_index, P3_index, P4_index, P5_Index, P6_index, P7_index FROM '1994_House_SupTab4'")
    
    
    #prepare data for PCA
    data_frameworkPCA <- data_complete_PCA[,3:12]
    
    #scale data and apply function prcomp for PCA
    pca_results_1<- prcomp(data_frameworkPCA, scale. = TRUE)
    fviz_pca_var(pca_results_1, col.var = "contrib", 
                 gradient.cols = c("yellow", "green", "red"),
                 ggtheme = theme_minimal(), repel=TRUE)
  })
  
  output$plot_4 <- renderPlot({
    #we open the connection to RSQLite
    toxi_con<- dbConnect(SQLite(), "Database_Toxicology_Own.db")
    # run our query
    data_complete_PCA <- dbGetQuery(toxi_con, "select Substance, Class, P3_7_index, P4_7_index, P5_7_index, P1_2_index, P1_index, P2_index, P3_index, P4_index, P5_Index, P6_index, P7_index FROM '1994_House_SupTab4'")
    #check output
    # print(data_complete_PCA)
    
    #prepare data for PCA
    data_frameworkPCA <- data_complete_PCA[,3:12]
    #check output
    #print(data_frameworkPCA)
    #scale data and apply function prcomp for PCA
    pca_results_1<- prcomp(data_frameworkPCA, scale. = TRUE)
    
    fviz_eig(pca_results_1, addlabels = T, label.size=40)
  })
  
  #TOCICITY CHARTS
  
  #Prepare to calculate Toxicity Indexes
  
  # Total number of samples that we have available from our initial data set for logit regression.
  toxi_con<- dbConnect(SQLite(), "Database_Toxicology_Own.db")
  data_complete_logit <- dbGetQuery(toxi_con, "select P1_index, P2_index, P3_index, P4_index, P5_index, P6_index, P7_index FROM '1994_House_SupTab4'")

  # We need to test the correlations
  output$plot_6 <- renderPlot({
    testcorr <- cor(data_complete_logit)
    round(testcorr, 2)
    corrplot(testcorr, type = "upper", order = "hclust", 
             tl.col = "red", tl.srt = 50)
    
  })
  
  output$plot_7 <- renderPlot({
    testcorr <- cor(data_complete_logit)
    round(testcorr, 2)
    corrplot(testcorr, method = 'number')
  })
  
  # Aromatic rings substances
  toxi_con<- dbConnect(SQLite(), "Database_Toxicology_Own.db")
  # Tables list
  dbListTables(toxi_con)
  
  arenes_data <- reactive({
    data_complete_arene <- dbGetQuery(toxi_con, "select * FROM '1994_House_SupTab2'")
    data_complete_arene <- data_complete_arene[order(data_complete_arene[,'Chemical.category'], decreasing = FALSE),]
    data_complete_arene <- data_complete_arene[,c("Chemical.Name","Chemical.category")]
    return(data_complete_arene)
  })
  
  
  output$table_arenes <- renderDT({
    req(arenes_data()) 
    datatable(arenes_data(),rownames = FALSE, options = list(pageLength = 5,searching = TRUE
      )
    )
  })
  
  #polr model
  data_complete_logit_reg_up<- reactive({
    data_complete_logit_reg <- dbGetQuery(toxi_con, "select P1_index, P2_index, P3_index, P4_index, P5_index, P6_index, P7_index, Substance, Class, mg, Per1, Per2, Per3, Per4, Per5, Per6, Per7, P3_7_index,P4_7_index,P5_7_index,P1_2_index,Primary_PAH FROM '1994_House_SupTab4'")
    # For this modelwe will have 4 slices
    number_of_chart_parts_model2 = 4
    print(number_of_chart_parts_model2)
    
    # For each part of the chart we can use 10 samples
    number_of_samples_per_part_model2 = 10
    #print(number_of_samples_per_part_model2)
    
    # Number of total samples that we need for this model
    total_samples_model2=number_of_chart_parts_model2*number_of_samples_per_part_model2
    #print(total_samples_model2)
    
    #We have 3 categories of toxicity in our models
    three_categories <- c(20,30,50)
    #print(three_categories)
    nc <- length(three_categories)
    #print(nc)
    
    #Calculate the total number of elements available
    
    total_number_rows <- nrow(data_complete_logit_reg)
    #print(total_number_rows)
    
    data_complete_logit_reg$eliminate <- rep(10, total_number_rows)
    # print(data_complete_logit_reg$eliminate)
    
    for (i in 1:nrow(data_complete_logit_reg)) {
      if(sum(data_complete_logit_reg$P1_index[i], data_complete_logit_reg$P2_index[i], data_complete_logit_reg$P3_index[i],
             data_complete_logit_reg$P4_index[i], data_complete_logit_reg$P5_index[i], data_complete_logit_reg$P6_index[i],
             data_complete_logit_reg$P7_index[i]) == 0){
        data_complete_logit_reg$eliminate[i] <- 20
      }
    }
    #print(data_complete_logit_reg$eliminate)
    
    #The next operation that we need to perform is to 
    #decide how many elements we need to have for each of the three categories.
    
    c1<-round(total_number_rows*three_categories[1]/100)
    print(c1)
    categories_vector <- c1
    print(categories_vector)
    
    c2 <- round(total_number_rows*three_categories[2]/100)
    print(c2)
    
    categories_vector <- c(c1, c2)
    print(categories_vector)
    
    c3 <- total_number_rows-c1-c2
    categories_vector <- c(c1, c2, c3)
    print(categories_vector)
    
    # For Model 2:
    
    categories_model_2 <- round(total_samples_model2*(categories_vector/total_number_rows))
    print(categories_model_2)
    
    if (sum(categories_model_2)>total_samples_model2) {categories_model_2[nc] <- categories_model_2[nc]-1}
    else if (sum(categories_model_2)>total_samples_model2) {categories_model_2[nc] <-categories_model_2[nc]+1}
    print(categories_model_2)
    
    #We have 3 categories of chemicals
    
    chemical_list1 <- list("195_UDAE", "102_CGO", "034_HFO", "058_HFO", "008_HFO", "089_UDAE", 
                           "134_HFO", "171_SRGO", "164_HFO", "182_VHGO", "172_SRGO", "106_CGO", 
                           "129_HFO", "130_CGO", "176_VHGO", "021_HFO", "135_TDAE", "096A_UDAE",
                           "078_HFO", "184_VHGO", "168_SRGO", "028_HFO", "097_HFO", "083_UATO",
                           "096B_UDAE", "070_CGO", "183_VHGOa", "064_UATO")
    
    
    chemical_list2<- list("104_kerosine","071_HFO","091_HFO","177_VHGO",       
                          "175_VHGO","050_HFO","031_HFO","080_HFO",        
                          "082_UATO","169_SRGO", "117_OLBO", "180_VHGO",     
                          "018_HFO", "197_HFO", "003_CGO", "166_HFO",      
                          "024_HFO", "187_1__OGO","178_VHGO", "193_OXASPH",
                          "012_CGO", "084_UATO","020_HFO", "025_HFO",       
                          "006_HFO", "170_SRGO", "017_HFO", "188_SRGO",     
                          "186_RAE", "192_BITUMENS", "075_OLBO", "079_HFO",        
                          "036_OLBO", "155_HFO", "173_OGO", "191_BITUMEN",   
                          "123_kerosine", "174_OGO", "179_VHGO", "185_RAE",      
                          "189_BITUMEN", "124_gasoline")
    
    chemical_list3<- list("069_TDAE", "146_gasoline", "139_OLBO", "118_OLBO", "081_OLBO", "149_OLBO",       
                          "131_HFO", "016_kerosine", "049_kerosine", "190_BITUMENS",   
                          "132_gasoline", "121_CGO", "073_OLBO", "152_slackwax",   
                          "041_CGO", "181_VHGO", "007_HFO", "090_FO",         
                          "160_FO", "092_OLBO", "076_gasoline", "061_P_Hwax",     
                          "085A_OLBO", "115_OLBO", "011_kerosine", "187_2__OGO",    
                          "067_P_Hwax", "085B_OLBO", "140_OLBO", "114_OLBO",     
                          "059_dieselfuel", "128_kerosine", "085C_OLBO", "147_OLBO",       
                          "046_gasoline", "127_slackwax", "074_OLBO", "148_OLBO",     
                          "138_OLBO", "060_OLBO", "119_OLBO", "108_OLBO",     
                          "141_gasoline", "151_OLBO", "009_gasoline", "063_P_Hwax",   
                          "085D_OLBO", "066_OLBO", "202_Paraffinwax", "112_OLBO",       
                          "087_kerosine", "072_OLBO", "109_FO", "113_OLBO",      
                          "086_kerosine", "062_P_Hwax", "196_petrolatum", "136_slackwax",  
                          "111_petrolatum", "153_OLBO", "056_gasoline", "137_OLBO",    
                          "120_gasoline", "093_petrolatum", "043_kerosine", "103_gasoline",  
                          "154_OLBO", "107_P_Hwax", "065_P_Hwax", "145_HRBO",       
                          "150_OLBO")
    
    # We eliminate elements with 0 for indexes
    
    data_for_reg <- data_complete_logit_reg[data_complete_logit_reg$eliminate != 20,]
    #print(nrow(data_for_reg))
    
    set.seed(200)
    datarandom1_m2 <- data_for_reg[ sample( which( data_for_reg$Substance %in% chemical_list1 ) , categories_model_2[1] ) , ]
    #print(datarandom1_m2)
    datarandom1_m2$category <- rep(1, nrow(datarandom1_m2))
    print(datarandom1_m2)
    
    datarandom2_m2 <- data_for_reg[ sample( which( data_for_reg$Substance %in% chemical_list2 ) , categories_model_2[2] ) , ]
    print(datarandom2_m2)
    datarandom2_m2$category <- rep(2, nrow(datarandom2_m2))
    print(datarandom2_m2)
    
    datarandom3_m2 <- data_for_reg[ sample( which( data_for_reg$Substance %in% chemical_list3 ) , categories_model_2[3] ) , ]
    print(datarandom3_m2)
    datarandom3_m2$category <- rep(3, nrow(datarandom3_m2))
    print(datarandom3_m2)
    
    data_model2_logit_0 <- rbind(datarandom1_m2,datarandom2_m2)
    data_model2_logit<- rbind(data_model2_logit_0, datarandom3_m2)
    print("Data for Model2:")
    print(data_model2_logit)
    print(nrow(data_model2_logit))
    
    # prepare data for polr and transform the dependent variables in a factor
    m2data <- data_model2_logit
    m2data$category <- as.factor(m2data$category)
    
    #order
    m2data$category <- ordered(m2data$category, levels= c(3, 2, 1))
    print(levels(m2data$category))
    print("md2data is:")
    print(m2data)
    return(m2data)
  })
  
  output$indep_index <- renderUI({
    req(data_complete_logit_reg_up())
    selectInput("independent_indexes", "Select Independent Variables", choices = names(data_complete_logit_reg_up()), multiple = TRUE)
  })
  
  polr_chem_model <- eventReactive(input$our_model, {
    req(data_complete_logit_reg_up())
    
    m2data <- data_complete_logit_reg_up()
    
    req(input$independent_indexes)
    
    print("independent indexes are:")
    print(input$independent_indexes)
    
    formula <- as.formula(paste("category", "~", paste(input$independent_indexes, collapse = " + ")))
    
    polr_chem_model <- polr(formula, data = m2data, Hess = TRUE)
  
    return(polr_chem_model)
  })
  
  output$chemical_summary <- renderPrint({
    req(polr_chem_model())
    summary(polr_chem_model())
  })
  
  output$chemical_stargazer <- renderUI({
    req(polr_chem_model())
    
    # Generate the stargazer output as HTML
    chemical_stargazer <- stargazer(polr_chem_model(), type = "html", out.header = FALSE)
    
    # Create an HTML output
    HTML(chemical_stargazer)
  })
  
  
  
  data_complete_logit_reg_up2<- reactive({
    data_complete_logit_reg <- dbGetQuery(toxi_con, "select P1_index, P2_index, P3_index, P4_index, P5_index, P6_index, P7_index, Substance, Class, mg, Per1, Per2, Per3, Per4, Per5, Per6, Per7, P3_7_index,P4_7_index,P5_7_index,P1_2_index,Primary_PAH FROM '1994_House_SupTab4'")
    #print(data_complete_logit_reg)
    
    req(input$user_slices)
    
    # For this model we will have 7 slices
    number_of_chart_parts_model1 = 7
    #print(number_of_chart_parts_model1)
    
    # For each part of the chart we can use 7 samples
    number_of_samples_per_part_model1 = as.numeric(input$user_slices)
    #print(number_of_samples_per_part_model1)
    
    # Number of total samples that we need for Model 1
    total_samples_model1=number_of_chart_parts_model1*number_of_samples_per_part_model1
    #print(total_samples_model1)
    
    #We have 3 categories of toxicity in our models 1 and 2
    three_categories <- c(20,30,50)
    #print(three_categories)
    nc <- length(three_categories)
    #print(nc)
    
    #Calculate the total number of elements available
    
    total_number_rows <- nrow(data_complete_logit_reg)
    #print(total_number_rows)
    
    #First, we need to choose the independent variables that we need
    indexes <- c('P3_7_index', 'P4_7_index', 'P5_7_index', 'P1_2_index', 'P1_index', 
                 'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index')
    #print(indexes)
    
    #From the initial dataset we need to eliminate the rows where all independent variables needed for the 
    #logit regression are 0.
    
    data_complete_logit_reg$eliminate <- rep(10, total_number_rows)
    # print(data_complete_logit_reg$eliminate)
    
    for (i in 1:nrow(data_complete_logit_reg)) {
      if(sum(data_complete_logit_reg$P1_index[i], data_complete_logit_reg$P2_index[i], data_complete_logit_reg$P3_index[i],
             data_complete_logit_reg$P4_index[i], data_complete_logit_reg$P5_index[i], data_complete_logit_reg$P6_index[i],
             data_complete_logit_reg$P7_index[i]) == 0){
        data_complete_logit_reg$eliminate[i] <- 20
      }
    }
    #print(data_complete_logit_reg$eliminate)
    
    #The next operation that we need to perform is to 
    #decide how many elements we need to have for each of the three categories.
    
    c1<-round(total_number_rows*three_categories[1]/100)
    #print(c1)
    categories_vector <- c1
    #print(categories_vector)
    
    c2 <- round(total_number_rows*three_categories[2]/100)
    #print(c2)
    
    categories_vector <- c(c1, c2)
    #print(categories_vector)
    
    c3 <- total_number_rows-c1-c2
    categories_vector <- c(c1, c2, c3)
    #print(categories_vector)
    
    # For Model 1:
    
    categories_model_1 <- round(total_samples_model1*(categories_vector/total_number_rows))
    #print(categories_model_1)
    
    if (sum(categories_model_1)>total_samples_model1) {categories_model_1[nc] <-categories_model_1[nc]-1}
    else if (sum(categories_model_1)>total_samples_model1) {categories_model_1[nc] <-categories_model_1[nc]+1}
    #print(categories_model_1)
    
    #We have 3 categories of chemicals
    
    chemical_list1 <- list("195_UDAE", "102_CGO", "034_HFO", "058_HFO", "008_HFO", "089_UDAE", 
                           "134_HFO", "171_SRGO", "164_HFO", "182_VHGO", "172_SRGO", "106_CGO", 
                           "129_HFO", "130_CGO", "176_VHGO", "021_HFO", "135_TDAE", "096A_UDAE",
                           "078_HFO", "184_VHGO", "168_SRGO", "028_HFO", "097_HFO", "083_UATO",
                           "096B_UDAE", "070_CGO", "183_VHGOa", "064_UATO")
    
    
    chemical_list2<- list("104_kerosine","071_HFO","091_HFO","177_VHGO",       
                          "175_VHGO","050_HFO","031_HFO","080_HFO",        
                          "082_UATO","169_SRGO", "117_OLBO", "180_VHGO",     
                          "018_HFO", "197_HFO", "003_CGO", "166_HFO",      
                          "024_HFO", "187_1__OGO","178_VHGO", "193_OXASPH",
                          "012_CGO", "084_UATO","020_HFO", "025_HFO",       
                          "006_HFO", "170_SRGO", "017_HFO", "188_SRGO",     
                          "186_RAE", "192_BITUMENS", "075_OLBO", "079_HFO",        
                          "036_OLBO", "155_HFO", "173_OGO", "191_BITUMEN",   
                          "123_kerosine", "174_OGO", "179_VHGO", "185_RAE",      
                          "189_BITUMEN", "124_gasoline")
    
    chemical_list3<- list("069_TDAE", "146_gasoline", "139_OLBO", "118_OLBO", "081_OLBO", "149_OLBO",       
                          "131_HFO", "016_kerosine", "049_kerosine", "190_BITUMENS",   
                          "132_gasoline", "121_CGO", "073_OLBO", "152_slackwax",   
                          "041_CGO", "181_VHGO", "007_HFO", "090_FO",         
                          "160_FO", "092_OLBO", "076_gasoline", "061_P_Hwax",     
                          "085A_OLBO", "115_OLBO", "011_kerosine", "187_2__OGO",    
                          "067_P_Hwax", "085B_OLBO", "140_OLBO", "114_OLBO",     
                          "059_dieselfuel", "128_kerosine", "085C_OLBO", "147_OLBO",       
                          "046_gasoline", "127_slackwax", "074_OLBO", "148_OLBO",     
                          "138_OLBO", "060_OLBO", "119_OLBO", "108_OLBO",     
                          "141_gasoline", "151_OLBO", "009_gasoline", "063_P_Hwax",   
                          "085D_OLBO", "066_OLBO", "202_Paraffinwax", "112_OLBO",       
                          "087_kerosine", "072_OLBO", "109_FO", "113_OLBO",      
                          "086_kerosine", "062_P_Hwax", "196_petrolatum", "136_slackwax",  
                          "111_petrolatum", "153_OLBO", "056_gasoline", "137_OLBO",    
                          "120_gasoline", "093_petrolatum", "043_kerosine", "103_gasoline",  
                          "154_OLBO", "107_P_Hwax", "065_P_Hwax", "145_HRBO",       
                          "150_OLBO")
    
    # We eliminate elements with 0 for indexes
    
    data_for_reg <- data_complete_logit_reg[data_complete_logit_reg$eliminate != 20,]
    #print(nrow(data_for_reg))
    
    set.seed(10)
    datarandom1 <- data_for_reg[ sample( which( data_for_reg$Substance %in% chemical_list1 ) , categories_model_1[1] ) , ]
    datarandom1$category <- rep(1, nrow(datarandom1))
    #print(datarandom1)
    
    datarandom2 <- data_for_reg[ sample( which( data_for_reg$Substance %in% chemical_list2 ) , categories_model_1[2] ) , ]
    datarandom2$category <- rep(2, nrow(datarandom2))
    #print(datarandom2)
    
    datarandom3 <- data_for_reg[ sample( which( data_for_reg$Substance %in% chemical_list3 ) , categories_model_1[3] ) , ]
    datarandom3$category <- rep(3, nrow(datarandom3))
    #print(datarandom3)
    
    data_model1_logit_0 <- rbind(datarandom1,datarandom2)
    data_model1_logit<- rbind(data_model1_logit_0, datarandom3)
    #print("Data for Model1:")
    #print(data_model1_logit)
    #print(nrow(data_model1_logit))
    
    #We keep only the variables needed for logit regression
    data_model1_logit_m <- data_model1_logit[,c('Substance','Class', 'P1_index', 
                                                'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index', 'category')]
    
    
    #print(data_model1_logit_m)
    
    indexes_logit <- c('P1_index', 'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index')
    length(indexes_logit)
    #print(length(indexes_logit))
    
    #apply standardization
    for(i in indexes_logit){
      data_model1_logit_m[,i] <- (data_model1_logit_m[,i]-min(data_model1_logit_m[,i]))/(max(data_model1_logit_m[,i])-min(data_model1_logit_m[,i]))
    }
    
    #print(data_model1_logit_m)
    
    return(data_model1_logit_m) 
  })
  
  
  ord_chem_model <- eventReactive(input$our_model2, {
    req(data_complete_logit_reg_up2())
    
    data_model1_logit_m <- data_complete_logit_reg_up2()
    
    #We have 3 categories of toxicity in our models
    three_categories <- c(20,30,50)
    #print(three_categories)
    nc <- length(three_categories)
    #print(nc)
    
    categories_ord <- as.character(1:nc)
    #print(categories_ord)
    
    fac_variable <- factor(data_model1_logit_m[,'category'], ordered = TRUE, levels = categories_ord)
    #print(fac_variable)
    
    indexes_logit <- c('P1_index', 'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index')
    length(indexes_logit)
    #print(length(indexes_logit))
    
    d_var <- rep(TRUE, length(indexes_logit))
    #print(d_var)
    
    m1 <- ordinalNet(as.matrix(data_model1_logit_m[, indexes_logit]),  fac_variable, positiveID = d_var, lambdaVals = 0, maxiterOut = 200)
    return(m1)
  })
  
  output$chemical_summary2 <- renderPrint({
    req(ord_chem_model())
    #We have 3 categories of toxicity in our models
    three_categories <- c(20,30,50)
    #print(three_categories)
    nc <- length(three_categories)
    #print(nc)
    indexes_logit <- c('P1_index', 'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index')
    length(indexes_logit)
    #print(length(indexes_logit))
    m1 <- ord_chem_model()
    
    print(m1)
    
    # Extract and normalize coefficients
    coefs_matrix <- coef(m1, matrix=TRUE)
    print(coefs_matrix)
    
    gewichteord <- m1$coefs[nc:(length(indexes_logit) + nc - 1)]
    print("Weights:")
    print(gewichteord)
    
    gewichteord <- gewichteord / sum(gewichteord)
    print("Weights after normalization:")
    print(gewichteord)
  })
  
  
  chemical_summary3_asfunction <- function(){
  
    req(ord_chem_model())
    #We have 3 categories of toxicity in our models
    three_categories <- c(20,30,50)
    #print(three_categories)
    nc <- length(three_categories)
    #print(nc)
    indexes_logit <- c('P1_index', 'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index')
    length(indexes_logit)
    #print(length(indexes_logit))
    m1 <- ord_chem_model()
    
    
    # Extract and normalize coefficients
    coefs_matrix <- coef(m1, matrix=TRUE)
   
    
    gewichteord <- m1$coefs[nc:(length(indexes_logit) + nc - 1)]
   
    
    gewichteord <- gewichteord / sum(gewichteord)
    
    #Toxicity scores will be calculated
    
    data_for_scores <- dbGetQuery(toxi_con, "select P1_index, P2_index, P3_index, P4_index, P5_index, P6_index, P7_index, P3_7_index, Substance, Class FROM '1994_House_SupTab4'")
    
    
    #We first standardize all indexes
    for(i in indexes_logit){
      data_for_scores[,i] <- (data_for_scores[,i]-min(data_for_scores[,i]))/(max(data_for_scores[,i])-min(data_for_scores[,i]))
    }
    
    
    toxicity_s <- rowSums(data_for_scores[,indexes_logit]*gewichteord[col(data_for_scores[,indexes_logit])])
    #print(toxicity_s)
    logit_result <- cbind(data_for_scores, toxicity_s)
    #print(logit_result)
    logit_result <- logit_result[order(logit_result[,'toxicity_s'], decreasing = TRUE),]
    #print(logit_result)
    return(logit_result)
    
  }    
    
    output$chemical_summary3 <- renderDT({
    logit_result <- chemical_summary3_asfunction()
    corS_P1 <- round(cor(logit_result$toxicity_s, logit_result$P1_index, method="spearman" ), 2)
    print("Corelation between Toxicity score and P1_index:")
    print(corS_P1)
    corS_P2 <- round(cor(logit_result$toxicity_s, logit_result$P2_index, method="spearman" ), 2)
    print("Corelation between Toxicity score and P2_index:")
    print(corS_P2)
    corS_P3 <- round(cor(logit_result$toxicity_s, logit_result$P3_index, method="spearman" ),2)
    print("Corelation between Toxicity score and P3_index:")
    print(corS_P3)
    corS_P4 <- round(cor(logit_result$toxicity_s, logit_result$P4_index, method="spearman" ),2)
    print("Corelation between Toxicity score and P4_index:")
    print(corS_P4)
    corS_P5 <- round(cor(logit_result$toxicity_s, logit_result$P5_index, method="spearman" ),2)
    print("Corelation between Toxicity score and P5_index:")
    print(corS_P5)
    corS_P6 <- round(cor(logit_result$toxicity_s, logit_result$P6_index, method="spearman" ),2)
    print("Corelation between Toxicity score and P6_index:")
    print(corS_P6)
    corS_P7 <- round(cor(logit_result$toxicity_s, logit_result$P7_index, method="spearman" ),2)
    print("Corelation between Toxicity score and P7_index:")
    print(corS_P7)
    corS_P8 <- round(cor(logit_result$toxicity_s, logit_result$P3_7_index, method="spearman" ),2)
    print("Corelation between Toxicity score and P3_7_index:")
    print(corS_P8)
    
    veccorr <- c(corS_P1, corS_P2, corS_P3, corS_P4, corS_P5, corS_P6, corS_P7, corS_P8)
    
    tdf <- data.frame(
      Variables = c("P1_index", "P2_index", "P3_index", "P4_index", "P5_index", "P6_index", "P7_index", "P3_7_index"),
      Correlation_with_Toxicity_Score = veccorr
    )

    datatable(tdf)
    
  })
  
  output$plot_corr1<- renderPlot({
    req(ord_chem_model())
    #We have 3 categories of toxicity in our models
    three_categories <- c(20,30,50)
    #print(three_categories)
    nc <- length(three_categories)
    #print(nc)
    indexes_logit <- c('P1_index', 'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index')
    length(indexes_logit)
    #print(length(indexes_logit))
    m1 <- ord_chem_model()
    
    
    # Extract and normalize coefficients
    coefs_matrix <- coef(m1, matrix=TRUE)
    
    
    gewichteord <- m1$coefs[nc:(length(indexes_logit) + nc - 1)]
    
    
    gewichteord <- gewichteord / sum(gewichteord)
    
    #Toxicity scores will be calculated
    
    data_for_scores <- dbGetQuery(toxi_con, "select P1_index, P2_index, P3_index, P4_index, P5_index, P6_index, P7_index, P3_7_index, Substance, Class FROM '1994_House_SupTab4'")
    
    
    #We first standardize all indexes
    for(i in indexes_logit){
      data_for_scores[,i] <- (data_for_scores[,i]-min(data_for_scores[,i]))/(max(data_for_scores[,i])-min(data_for_scores[,i]))
    }
    
    
    toxicity_s <- rowSums(data_for_scores[,indexes_logit]*gewichteord[col(data_for_scores[,indexes_logit])])
    #print(toxicity_s)
    logit_result <- cbind(data_for_scores, toxicity_s)
    #print(logit_result)
    logit_result <- logit_result[order(logit_result[,'toxicity_s'], decreasing = TRUE),]
    #print(logit_result)
    
    logit_result$Class <- as.factor(logit_result$Class)
    
    pchem1 <- ggplot(logit_result, aes(x = P1_index, y = toxicity_s, color = Class)) +
      geom_point() + geom_smooth(method = "lm", color = "grey", se=FALSE) +
      labs(x = "P1 Index", y = "Toxicity Score", color = "Class") +
      theme_minimal()  + theme(legend.text = element_text(size = 7),
                              legend.title = element_text(size = 9), 
                           legend.key.size = unit(0.3, "cm"))
    
    pchem2 <- ggplot(logit_result, aes(x = P2_index, y = toxicity_s, color = Class)) +
      geom_point() + geom_smooth(method = "lm", color = "grey", se=FALSE) +
      labs(x = "P2 Index", y = "Toxicity Score", color = "Class") +
      theme_minimal()+ theme(legend.text = element_text(size = 7), 
                             legend.title = element_text(size = 9), 
                             legend.key.size = unit(0.3, "cm"))
    
    pchem3 <- ggplot(logit_result, aes(x = P3_index, y = toxicity_s, color = Class)) +
      geom_point() + geom_smooth(method = "lm", color = "grey", se=FALSE) +
      labs(x = "P3 Index", y = "Toxicity Score", color = "Class") +
      theme_minimal()+ theme(legend.text = element_text(size = 7), 
                             legend.title = element_text(size = 9), 
                             legend.key.size = unit(0.3, "cm"))
    
    pchem4 <- ggplot(logit_result, aes(x = P4_index, y = toxicity_s, color = Class)) +
      geom_point() + geom_smooth(method = "lm", color = "grey", se=FALSE) +
      labs(x = "P4 Index", y = "Toxicity Score", color = "Class") +
      theme_minimal()+ theme(legend.text = element_text(size = 7), 
                             legend.title = element_text(size = 9), 
                             legend.key.size = unit(0.3, "cm"))
    
    
    library(gridExtra)
  grid.arrange(pchem1, pchem2, pchem3, pchem4, ncol = 2)
    
  })
  
  
  output$plot_corr2<- renderPlot({
    req(ord_chem_model())
    #We have 3 categories of toxicity in our models
    three_categories <- c(20,30,50)
    #print(three_categories)
    nc <- length(three_categories)
    #print(nc)
    indexes_logit <- c('P1_index', 'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index')
    length(indexes_logit)
    #print(length(indexes_logit))
    m1 <- ord_chem_model()
    
    
    # Extract and normalize coefficients
    coefs_matrix <- coef(m1, matrix=TRUE)
    
    
    gewichteord <- m1$coefs[nc:(length(indexes_logit) + nc - 1)]
    
    
    gewichteord <- gewichteord / sum(gewichteord)
    
    #Toxicity scores will be calculated
    
    data_for_scores <- dbGetQuery(toxi_con, "select P1_index, P2_index, P3_index, P4_index, P5_index, P6_index, P7_index, P3_7_index, Substance, Class FROM '1994_House_SupTab4'")
    
    
    #We first standardize all indexes
    for(i in indexes_logit){
      data_for_scores[,i] <- (data_for_scores[,i]-min(data_for_scores[,i]))/(max(data_for_scores[,i])-min(data_for_scores[,i]))
    }
    
    
    toxicity_s <- rowSums(data_for_scores[,indexes_logit]*gewichteord[col(data_for_scores[,indexes_logit])])
    #print(toxicity_s)
    logit_result <- cbind(data_for_scores, toxicity_s)
    #print(logit_result)
    logit_result <- logit_result[order(logit_result[,'toxicity_s'], decreasing = TRUE),]
    #print(logit_result)
    
    logit_result$Class <- as.factor(logit_result$Class)
    
    pchem5 <- ggplot(logit_result, aes(x = P5_index, y = toxicity_s, color = Class)) +
      geom_point() +
      geom_smooth(method = "lm", color = "grey", se=FALSE) +
      labs(x = "P5 Index", y = "Toxicity Score", color = "Class") +
      theme_minimal()  + theme(legend.text = element_text(size = 7), 
                               legend.title = element_text(size = 9), 
                               legend.key.size = unit(0.3, "cm"))
    
    pchem6 <- ggplot(logit_result, aes(x = P6_index, y = toxicity_s, color = Class)) +
      geom_point() +
      geom_smooth(method = "lm", color = "grey", se=FALSE) +
      labs(x = "P6 Index", y = "Toxicity Score", color = "Class") +
      theme_minimal()+ theme(legend.text = element_text(size = 7), 
                             legend.title = element_text(size = 9), 
                             legend.key.size = unit(0.3, "cm"))
    
    pchem7 <- ggplot(logit_result, aes(x = P7_index, y = toxicity_s, color = Class)) +
      geom_point() +
      geom_smooth(method = "lm", color = "grey", se=FALSE) +
      labs(x = "P7 Index", y = "Toxicity Score", color = "Class") +
      theme_minimal()+ theme(legend.text = element_text(size = 7), 
                             legend.title = element_text(size = 9), 
                             legend.key.size = unit(0.3, "cm"))
    pchem3_7 <- ggplot(logit_result, aes(x = P3_7_index, y = toxicity_s, color = Class)) +
      geom_point() +
      geom_smooth(method = "lm", color = "grey", se=FALSE) +
      labs(x = "P3_7 Index", y = "Toxicity Score", color = "Class") +
      theme_minimal()+ theme(legend.text = element_text(size = 7), 
                             legend.title = element_text(size = 9), 
                             legend.key.size = unit(0.3, "cm"))
    
  
    library(gridExtra)
    grid.arrange(pchem5, pchem6, pchem7, pchem3_7, ncol = 2)
    
  })
  
  # We write a function that helps us create a Boxplot for substances toxicity by class.
  
  wtoxplot_asfunction <- function() {
    req(ord_chem_model())
    #We have 3 categories of toxicity in our models
    three_categories <- c(20,30,50)
    #print(three_categories)
    nc <- length(three_categories)
    #print(nc)
    indexes_logit <- c('P1_index', 'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index')
    length(indexes_logit)
    #print(length(indexes_logit))
    m1 <- ord_chem_model()
    
    
    # Extract and normalize coefficients
    coefs_matrix <- coef(m1, matrix=TRUE)
    
    
    gewichteord <- m1$coefs[nc:(length(indexes_logit) + nc - 1)]
    
    
    gewichteord <- gewichteord / sum(gewichteord)
    
    #Toxicity scores will be calculated
    
    data_for_scores <- dbGetQuery(toxi_con, "select P1_index, P2_index, P3_index, P4_index, P5_index, P6_index, P7_index, P3_7_index, Substance, Class FROM '1994_House_SupTab4'")
    
    
    #We first standardize all indexes
    for(i in indexes_logit){
      data_for_scores[,i] <- (data_for_scores[,i]-min(data_for_scores[,i]))/(max(data_for_scores[,i])-min(data_for_scores[,i]))
    }
    
    
    toxicity_s <- rowSums(data_for_scores[,indexes_logit]*gewichteord[col(data_for_scores[,indexes_logit])])
    print("this is toxicity score:")
    print(toxicity_s)
    logit_result <- cbind(data_for_scores, toxicity_s)
    print(logit_result)
    logit_result <- logit_result[order(logit_result[,'toxicity_s'], decreasing = TRUE),]
    print(logit_result)
    
    logit_result$Class <- as.factor(logit_result$Class)
    return(logit_result)
  
  }
  
  #We use the function from above to render our chemical substances boxplot
    output$wtoxplot <- renderPlot({      
    logit_result <- wtoxplot_asfunction()
    ggplot(logit_result, aes(y = Class , x = toxicity_s, fill = Class)) + geom_boxplot()+
      scale_x_continuous(
        limits = c(0, 0.95),  
        breaks = seq(0, 0.95, by = 0.1) 
      ) +
      theme_minimal() +
      labs(title = "Our calculated Toxicity Index by Chemical Class", 
           y = "Chemical Class", 
           x = "Our calcualted Toxicity Index")
   })


  chemcalpie_asfunction <- function(our_chemical_class) {
    req(ord_chem_model())
    #We have 3 categories of toxicity in our models
    three_categories <- c(20,30,50)
    #print(three_categories)
    nc <- length(three_categories)
    #print(nc)
    indexes_logit <- c('P1_index', 'P2_index', 'P3_index', 'P4_index', 'P5_index', 'P6_index', 'P7_index')
    length(indexes_logit)
    #print(length(indexes_logit))
    m1 <- ord_chem_model()
    
    
    # Extract and normalize coefficients
    coefs_matrix <- coef(m1, matrix=TRUE)
    
    
    gewichteord <- m1$coefs[nc:(length(indexes_logit) + nc - 1)]
    
    
    gewichteord <- gewichteord / sum(gewichteord)
    
    #Toxicity scores will be calculated
    
    data_for_scores <- dbGetQuery(toxi_con, "select P1_index, P2_index, P3_index, P4_index, P5_index, P6_index, P7_index, P3_7_index, Substance, Class FROM '1994_House_SupTab4'")
    
    
    #We first standardize all indexes
    for(i in indexes_logit){
      data_for_scores[,i] <- (data_for_scores[,i]-min(data_for_scores[,i]))/(max(data_for_scores[,i])-min(data_for_scores[,i]))
    }
    
    # We start preparing the Pie Charts
    print("needed data")
    print(data_for_scores)
    data_for_scores <- data_for_scores %>% filter(Class == our_chemical_class)
    print(data_for_scores)
    print(gewichteord)
    
    chemical_parts <- TxpSliceList(PAC_1cycle = TxpSlice("P1_index"),
                                   PAC_2cycles = TxpSlice("P2_index"),
                                   PAC_3cycles = TxpSlice("P3_index"),
                                   PAC_4cycles = TxpSlice("P4_index"),
                                   PAC_5cycles = TxpSlice("P5_index"),
                                   PAC_6cycles = TxpSlice("P6_index"),
                                   PAC_7cycles = TxpSlice("P7_index"))
    print(chemical_parts)
    chemical_tr<- TxpTransFuncList(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
    print(chemical_tr)
    print(gewichteord)
    chem_m <- TxpModel(txpSlices = chemical_parts, txpWeights = gewichteord, txpTransFuncs = chemical_tr)
    print(chem_m)
    
    chemr <- txpCalculateScores(model = chem_m,
                                    input = data_for_scores,
                                    id.var = 'Substance')
    print(chemr)
    return(chemr)
  }
    
  output$chemcalpie <- renderPlot({
    chemr <- chemcalpie_asfunction(input$chem_name)
    txpSliceScores(chemr)
    txpWeights(chemr)
    txpMissing(chemr)
    library("ggplot2")
    
    
    plot(chemr, package = "gg")
    
    }, width = 1200, height = 400)
  

# TESTING  
# We are now testing our dashboard. The code from below is very time consuming as it runs
# the dashboard 100 times for each function to calculate the average time.
  
# Calculation of the average time needed to run the function for chemical pies (chemcalpie_asfunction("UDAE")). 
  #trials_chem_pie <- numeric(100)
  #observe({
  #for (k in 1:100){
  # trials_chem_pie[k] <- system.time({
  #  chemcalpie_asfunction("UDAE")
  # })["elapsed"]
  #}
  
  #trial_avg <- mean(trials_chem_pie)
  #print(trial_avg)
  #print(summary(trials_chem_pie))
  #})
  
# Calculation of the average time needed to run the function for chemical substances boxplot.   
  #trials_wtoxplot_asfunction <- numeric(100)
  #observe({
  #for (m in 1:100){
  #trials_wtoxplot_asfunction[m] <- system.time({
  #  wtoxplot_asfunction()
  #})["elapsed"]
  #}
  
  #trial_avg2 <- mean(trials_wtoxplot_asfunction)
  #print(trial_avg2)
  #print(summary(trials_wtoxplot_asfunction))
  #})

# Calculation of the average time needed to run the function chemical_summary3_asfunction(). 
#trials_chemical_summary3_asfunction <- numeric(100)
#observe({
#for (k in 1:100){
# trials_chemical_summary3_asfunction[k] <- system.time({
#   chemical_summary3_asfunction()
# })["elapsed"]
#}

#trial_avg <- mean(trials_chemical_summary3_asfunction)
#print(trial_avg)
#print("Summary values:")
#print(summary(trials_chemical_summary3_asfunction))
#})

#END
})

# Shiny-App-Toxicology

#1. This web app uses the database file called Database_Toxicology_Own.db. All the needed data to run the app is already in this database.

#2. This app has been created in R Studio with version R 4.4.3 and you will need either this R version or a more recent one.

#3. You will need to install the ToxPi R package. The app does this automatically, but in case it does not work for you please install it from here (3 options): https://github.com/ToxPi/toxpiR.

#Literature Sources:
  
#  [1] J. S. House, F. A. Grimm, W. D. Klaren, A. Dalzell, S. Kuchi, S.-D. Zhang,
#K. Lenz, P. J. Boogaard, H. B. Ketelslegers, T. W. Gant, and et al., “Grouping
#of uvcb substances with new approach methodologies (nams) data,” ALTEX -
#  Alternatives to animal experimentation, vol. 38, p. 123–137, Jan. 2021.

#[2] F. Sewell, C. Alexander-White, S. Brescia, R. A. Currie, R. Roberts, C. Roper,
#C. Vickers, C. Westmoreland, and I. Kimber, “New approach methodologies
#(nams): identifying and overcoming hurdles to accelerated adoption,” Toxicol-
#  ogy Research, vol. 13, p. tfae044, 03 2024.

#[3] J. F. Fleming, D. L. Filer, D. T. Lloyd, P. Thunga, S. W. Marvel, A. A.
#Motsinger-Reif, and D. M. Reif, “toxpir: Create toxpi prioritization mod-
#  els,” 2024d. R package, version 1.3.0, https://cran.r-project.org/web/packages/toxpiR/index.html.

#[4] N. Kleinstreuer and T. Hartung, “Artificial intelligence (ai)—it’s the end of
#the tox as we know it (and i feel fine),” Archives of Toxicology, vol. 98, no. 3,
#pp. 735–754, 2024.

#[5] A. Lai, A. M. Clark, B. I. Escher, M. Fernandez, L. R. McEwen, Z. Tian,
#Z. Wang, and E. L. Schymanski, “The next frontier of environmental unknowns:
#  substances of unknown or variable composition, complex reaction products, or
#biological materials (uvcbs),” Environmental Science & Technology, vol. 56,
#no. 12, pp. 7448–7466, 2022.

#[6] T. Bl¨ummel, J. Rehn, C. Mereu, F. Graf, F. Bazing, C. Kneuer, A. Sonnen-
#  burg, P. Wittkowski, F. Padberg, K. Bech, et al., “Exploring the use of artificial
#intelligence (ai) for extracting and integrating data obtained through new ap-
#  proach methodologies (nams) for chemical risk assessment,” EFSA Supporting
#Publications, vol. 21, no. 1, p. 8567E, 2024.

#[7] J. F. Fleming, J. S. House, J. R. Chappel, A. A. Motsinger-Reif, and D. M.
#Reif, “Guided optimization of toxpi model weights using a semi-automated
#approach,” Computational Toxicology, vol. 29, p. 100294, 2024a.

#[8] H. Kamp, N. A. Kocabas, F. Faulhammer, N. Synhaeve, E. Rushton, B. Flick,
#V. Giri, S. Sperber, L. Higgins, M. Penman, et al., “Utility of in vivo
#metabolomics to support read-across for uvcb substances under reach,” Archives
#of Toxicology, vol. 98, no. 3, pp. 755–768, 2024.

#[9] S. Heux, T. Fuchs, J. Buhmann, N. Zamboni, and U. Sauer, “A high-throughput
#metabolomics method to predict high concentration cytotoxicity of drugs from
#low concentration profiles,” Metabolomics, vol. 8, pp. 433–443, 06 2011.

#[10] C. B. Pestana, J. W. Firman, and M. T. Cronin, “Incorporating lines of ev-
#  idence from new approach methodologies (nams) to reduce uncertainties in a
#category based read-across: a case study for repeated dose toxicity,” Regulatory
#Toxicology and Pharmacology, vol. 120, p. 104855, 2021.

#[11] S. W. Marvel, K. To, F. A. Grimm, F. A. Wright, I. Rusyn, and D. M. Reif,
#“Toxpi graphical user interface 2.0: Dynamic exploration, visualization, and
#sharing of integrated data models,” BMC bioinformatics, vol. 19, pp. 1–7, 2018.

#[12] D. M. Reif, M. T. Martin, S. W. Tan, K. A. Houck, R. S. Judson, A. M.
#Richard, T. B. Knudsen, D. J. Dix, and R. J. Kavlock, “Endocrine profiling and
#prioritization of environmental chemicals using toxcast data,” Environmental
#health perspectives, vol. 118, no. 12, pp. 1714–1720, 2010.

#[13] “Appendix c toxicological priority index (toxpi).” National Research Council.
#A Framework to Guide Selection of Chemical Alternatives. Washington, DC:
#  The National Academies Press, 2014. https://nap.nationalacademies.org/read/18872/chapter/19, accessed on 2025-04-23.

#[14] J. Fleming, D. Filer, D. Lloyd, P. Thunga, S. W. Marvel, A. A. Motsinger-
#  Reif, and D. M. Reif, “Toxpi r introduction,” 2024b. https://cran.rstudio.com/web/packages/toxpiR/vignettes/introduction.html#:~:text=What%20is%20ToxPi%3F,i.e.%20multimodal%20or%20multiscale%20information, accessed on 2024-11-21.

#[15] J. Fleming, S. W. Marvel, S. Supak, A. A. Motsinger-Reif, and D. M. Reif,
#“Toxpi* gis toolkit: creating, viewing, and sharing integrative visualizations
#for geospatial data using arcgis,” Journal of Exposure Science & Environmental
#Epidemiology, vol. 32, no. 6, pp. 900–907, 2022.

#[16] D. M. Reif, M. Sypa, E. F. Lock, F. A. Wright, A. Wilson, T. Cathey, R. R.
#Judson, and I. Rusyn, “Toxpi gui: an interactive visualization tool for trans-
#  parent integration of data from diverse sources of evidence,” Bioinformatics,
#vol. 29, no. 3, pp. 402–403, 2013.

#[17] J. G. Speight, Environmental organic chemistry for engineers. Butterworth-
#  Heinemann, 2016.

#[18] S. Katoch, S. S. Chauhan, and V. Kumar, “A review on genetic algorithm: past,
#present, and future,” Multimedia Tools and Applications, vol. 80, pp. 8091–8126,
#Feb 2021.

#[19] “Polycyclic aromatic compounds.” National Toxicology Program, U.S. Depart-
#  ment of Health and Human Services, 2014. https://ntp.niehs.nih.gov/whatwestudy/topics/pacs, accessed on 2024-11-21.

#[20] O. C. Ifegwu and C. Anyakora, “Chapter six - polycyclic aromatic hydrocar-
#  bons: Part i. exposure,” in Advances in Clinical Chemistry (G. S. Makowski,
#                                                              ed.), vol. 72, pp. 277–304, Elsevier, 2015.

#[21] A. Schatten, M. Demolsky, D. Winkler, S. Biffl, E. Gostischa-Franta, and
#T. ¨Ostreicher, Best Practice Software-Engineering. Verlag Heidelberg, 2010.

#[22] H. Wickham and K. M¨uller, “Dbi: R database interface,” 2024. R package, ver-
#  sion 1.2.3, https://cran.r-project.org/web/packages/DBI/index.html.

#[23] K. M¨uller and H. Wickham, “Rsqlite: Sqlite interface for r,” 2024. R pack-
#  age, version 2.3.9, https://cran.r-project.org/web/packages/RSQLite/index.html.

#[24] H. Wickham, J. Bryan, M. Kalicinski, K. Valery, C. Leitienne, B. Colbert,
#D. Hoerl, and E. Miller, readxl: Read Excel Files, 2025. R package version
#1.4.5, https://cran.r-project.org/web/packages/readxl/index.html.

#[25] S. Karnouskos, R. Sinha, P. Leit˜ao, L. Ribeiro, and T. I. Strasser, “The appli-
#  cability of iso/iec 25023 measures to the integration of agents and automation
#systems,” in IECON 2018 - 44th Annual Conference of the IEEE Industrial
#Electronics Society, pp. 2927–2934, 2018.

#[26] C. Abras, D. Maloney-Krichmar, and J. Preece, “User-centered design,” User-
#  Centered Design, pp. 445–456, 01 2004.

#[27] International Organization for Standardization (ISO) and International Elec-
#  trotechnical Commission (IEC), ISO/IEC 25023:2016 - Systems and software
#engineering — Systems and software Quality Requirements and Evaluation
#(SQuaRE) — Measurement of system and software product quality, 2016. On-
#  line at: https://www.iso.org/standard/35747.html.

#[28] T. Galli, F. Chiclana, and F. Siewe, “Software product quality models, devel-
#  opments, trends, and evaluation,” SN Computer Science, vol. 1, 05 2020.

#[29] M. J. Wurm, P. J. Rathouz, and B. M. Hanlon, “Regularized ordinal regression
#and the ordinalNet R package,” Journal of Statistical Software, vol. 99, no. 6,
#pp. 1–42, 2021.

#[30] J. G. Speight, “Chapter 2 - organic chemistry,” in Environmental Or-
#  ganic Chemistry for Engineers (J. G. Speight, ed.), pp. 43–86, Butterworth-
#  Heinemann, 2017.

#[31] “The r foundation, the comprehensive r archive network,” 2025. R version 4.4.3,
#https://cran.r-project.org/bin/windows/base/.

#[32] W. Chang, J. Cheng, J. Allaire, C. Sievert, B. Schloerke, Y. Xie, J. Allen,
#J. McPherson, A. Dipert, and B. Borges, “shiny: Web application framework
#for r,” 2024. R version 1.10.0, https://cran.r-project.org/web/packages/shiny/index.html.

#[33] W. Chang and B. B. Ribeiro, shinydashboard: Create Dashboards with ’Shiny’,
#2021. R package version 0.7.2, https://cran.r-project.org/web/packages/shinydashboard/index.html.

#[34] C. M. Bishop and N. M. Nasrabadi, Pattern recognition and machine learning,
#vol. 4. Springer, 2006.

#[35] A. Kassambara and F. Mundt, factoextra: Extract and Visualize the Results of
#Multivariate Data Analyses, 2020. R version 1.0.7, https://cran.r-project.org/web/packages/factoextra/index.html.

#[36] M. M. Deza and E. Deza, Encyclopedia of Distances. Springer, 1st ed., 2009.

#[37] D. Jurafsky and J. H. Martin, Speech and Language Processing: An Introduc-
#  tion to Natural Language Processing, Computational Linguistics, and Speech
#Recognition with Language Models. published at: https://web.stanford.edu/~jurafsky/slp3/, 3rd ed., 2025. Online manuscript released January 12, 2025.

#[38] B. Ripley, B. Venables, D. M. Bates, K. Hornik, A. Gebhardt, and
#D. Firth, MASS: Support Functions and Datasets for Venables and Ripley’s
#MASS, 2025. R package version 7.3-64. https://cran.r-project.org/web/packages/MASS/index.html.

#[39] B. Ripley, B. Venables, D. M. Bates, K. Hornik, A. Gebhardt, and D. Firth,
#“Mass: Support functions and datasets for venables and ripley’s mass,” 2025. R
#package, version 7.3-64, https://cran.r-project.org/web/packages/MASS/index.html.

#[40] H. Wickham, W. Chang, L. Henry, T. L. Pedersen, K. Takahashi, C. Wilke,
#K. Woo, H. Yutani, D. Dunnington, and T. van den Brand, ggplot2: Create
#Elegant Data Visualisations Using the Grammar of Graphics, 2024. R pack-
#  age version 3.5.1, https://cran.r-project.org/web/packages/ggplot2/index.html.

#[41] Y. Xie, J. Cheng, X. Tan, J. Allaire, M. Girlich, G. F. Ellis, J. Rauh, S. Lim-
#  ited, B. Reavis, L. Gersen, B. Szopka, A. Pickering, W. Holmes, M. Marttila,
#A. Quintero, and S. Laurent, DT: A Wrapper of the JavaScript Library ’DataT-
#  ables’, 2024. R version 0.33, https://cran.r-project.org/web/packages/DT/index.html.

#[42] B. Auguie and A. Antonov, gridExtra: Miscellaneous Functions for ”Grid”
#Graphics, 2017. R package version 2.3, https://cran.r-project.org/web/packages/gridExtra/index.html.

#[43] S. Marvel and D. M. Reif, “Toxpi™ user manual v2.3,” 2018. https://toxpi.org/dist/ToxPi%20User%20Manual.pdf, accessed on 2024-11-21.

#[44] J. Fleming, D. Filer, D. Lloyd, P. Thunga, S. W. Marvel, A. A. Motsinger-Reif,
#and D. M. Reif, “Manual package ‘toxpir’,” 2024c. https://cran.r-project.org/web/packages/toxpiR/toxpiR.pdf, accessed on 2024-11-21.

#[45] A. Kassambara and F. Mundt, “factoextra : Extract and visualize the results
#of multivariate data analysesg,” 2020. https://cran.r-project.org/web/packages/factoextra/readme/README.html.

#[46] J. R. Groff, P. N. Weinberg, and A. J. Oppel, SQL: The Complete Reference,
#Third Edition. McGraw-Hill, 2010.

#[47] S. Parth, “Understanding k-means and hierarchical clus-
#  tering in r,” 2020. https://medium.com/@shreytparth/understanding-k-means-and-hierarchical-clustering-in-r-9899342cc4fb,
#accessed on 2024-11-21.

#[48] Z. Keita, “Introduction to principal component analysis (pca),” 2023. https://www.datacamp.com/tutorial/pca-analysis-r.

#[49] A. Kassambara, “Pca - principal component analysis essen-
#  tials,” 2017. https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/.

#[50] “R shiny tutorial — r shiny dashboard —enabling menu items for their re-
#  spective pages — r programming.” Data Science Tutorial, 2017. https://www.youtube.com/watch?v=fUXBL5bk20M.

#[51] “Shiny web app tutorial — how to upload a file in shiny app — r programming
#tutorial.” Data Science Tutorial, 2020. https://www.youtube.com/watch?v=A6VYSCB0TZM&t=199s.

#[52] “R shiny basics and user interface.” Statistical Learning Group, 2020. https://www.youtube.com/watch?v=6mJaw5pLtso.

#[53] Y. Dai and Y. Xiang, “Connecting sqlite database and shiny app for
#business intelligence,” 2020. https://shanghai.hosting.nyu.edu/data/r/case-3-sql-shiny.html, accessed on 2024-11-21.
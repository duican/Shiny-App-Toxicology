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

shinyUI(
  dashboardPage(
    skin = "green",
    dashboardHeader(title="Toxicology Analyis"),
    
    dashboardSidebar(
      sidebarMenu(
      menuItem("K-Means", tabName = "kmeans", icon = icon("home")),
      menuItem("Hierarchical Clustering", tabName = "hierarchical", icon = icon("home")),
      menuItem("Principal Component Analysis", tabName = "pca", icon = icon("home")),
      menuItem("Toxicity Charts", tabName = "toxpi", icon = icon("home")),
      menuItem("Upload Data", tabName = "upload", icon = icon("calendar"))
    )),
    
    dashboardBody(
      tabItems(
        tabItem(tabName = "kmeans",h4("K-Means"),
               fluidRow(column(12,
               tabBox(id="tabchart1",width=12,
                      tabPanel("K-Means results",numericInput('number_cluster', 'Please choose the prefered number of clusters:', 4, min = 1, max = 10), plotOutput('plot1')),
                      tabPanel("Elbow Plot", plotOutput('plotcluster'))
                      )
               )),
        ),
        tabItem(tabName = "hierarchical", h4("Hierarchical Clustering"),
                fluidRow(column(
                  height = 12,
                  width = 12
                ),
                tabBox(width='100%', height='100%',id="tabchart2",
                   splitLayout(cellWidths = c("50%", "50%"),radioButtons("distancemethod", "Which distance method do you prefer?",
                                                                              choices = c("manhattan","euclidian", "minkowski")), radioButtons("type", "Which shape do you prefer for the chart?",
                                                                                                                                                                     choices = c("rectangle", "circular", "phylogenic")))),
                
                #fluidRow(column(12,selectInput("distancemethod", "Which distance method do you prefer?",
                                               # choices = c("manhattan","euclidian", "maximum", "canberra", "minkowski")
                                              #  ))),
                #fluidRow(column(12,selectInput("type", "Which shape do you prefer for the chart?",
                                              # choices = c("rectangle", "circular", "phylogenic")
               # ))),
                plotOutput('plot2'))
               ),
        tabItem(tabName = "toxpi", h4("Toxicity Charts"),
                fluidRow(  column(
                  height = 12,
                  width = 12
                ),
                tabBox(width='100%', height='100%',id="tabchart4",
                       tabPanel("Features Correlations", 
                                fluidRow( column(6,plotOutput('plot_6')),
                                         column(6,plotOutput('plot_7'))
                                )),
                       tabPanel("Aromatic rings examples", 
                                fluidRow( column(6,DTOutput("table_arenes"))
                                ) 
                       ),
                       tabPanel("Logistic Regressions",  
                         uiOutput("indep_index"),
                         actionButton("our_model", "Run polr"),
                         verbatimTextOutput("chemical_summary"),
                        uiOutput("chemical_stargazer")
                         ),
                       tabPanel("OrdinalNET",  
                                fluidPage(
                                  selectInput("user_slices", "Select the number of slices:", choices = c(7, 8, 9, 10, 11)),
                                  textOutput("received_slices")
                                ),
                                actionButton("our_model2", "Run ord reg"),
                               verbatimTextOutput("chemical_summary2")
                       ),
                       tabPanel("Correlation values",  
                                fluidPage(
                                  selectInput("user_slices", "Select the number of slices:", choices = c(7, 8, 9, 10, 11)),
                                  textOutput("received_slices")
                                ),
                                actionButton("our_model2", "Show results"),
                                fluidRow( column(6,DTOutput("chemical_summary3"))),
                                plotOutput('wtoxplot')
                       ),
                       tabPanel("Correlation Plots",  
                                
                                actionButton("our_model2", "Show charts:"),
                                plotOutput("plot_corr1"),
                                plotOutput("plot_corr2")
                       ),
                       tabPanel("Pie Plots",  
                                radioButtons("chem_name", "Which chemical do you prefer?",
                                             choices = c("UDAE", "CGO", "HFO", "UATO","BIT", "VHGO", "NAPTHA", "OGO", "RAE", "SRGO", "TDAE", "P.LAT","BO","WAX","FO", "KER"), inline = TRUE),
                                actionButton("our_model2", "Show charts:"),
                                plotOutput("chemcalpie")
                       ),
                       
                )
                )
                ),
        tabItem(tabName = "pca", h4("Principal Component Analysis"),
                fluidRow(  column(
                  height = 12,
                  width = 12
                ),
                  tabBox(width='100%', height='100%',id="tabchart3",
                         tabPanel("PCA results", plotOutput('plot3')),
                         tabPanel("Percentage of explained variance", plotOutput('plot_4')),
                         tabPanel("PCA variables", plotOutput('plot_5')))
                )
        ),
        tabItem(tabName = "upload",
                fluidPage(
                  titlePanel("Upload your chemical data files"),
                  fileInput("chemical_files", "Select your needed files", 
                            multiple = TRUE, 
                            accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".xls", ".xlsx")),
                  actionButton("press_here", "Press here to upload chemicals"), 
                  verbatimTextOutput("happening")
                )
                
        )  
        )
      )
    )
)

#Sources:
#1. ToxPi in R installation: https://github.com/ToxPi/toxpiR
#2. R shiny basics and user interface.” Statistical Learning Group, 2020. https://www.youtube.com/watch?v=6mJaw5pLtso.
#3. “Shiny web app tutorial — how to upload a file in shiny app — r programming tutorial.” Data Science Tutorial, 2020. https://www.youtube.com/watch?v=A6VYSCB0TZM&t=199s.
#4. “R shiny tutorial — r shiny dashboard —enabling menu items for their respective pages — r programming.” Data Science Tutorial, 2017. https://www.youtube.com/watch?v=fUXBL5bk20M.
#5. Z. Keita, “Introduction to principal component analysis (pca),” 2023. https://www.datacamp.com/tutorial/pca-analysis-r.
#6. A. Kassambara, “Pca - principal component analysis essentials,” 2017. https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/.
#7.  J. Fleming, D. Filer, D. Lloyd, P. Thunga, S. W. Marvel, A. A. Motsinger-Reif, and D. M. Reif, “Toxpi r introduction,” 2024b. https://cran.rstudio.com/web/packages/toxpiR/vignettes/introduction.html#:~:text=What%20is%20ToxPi%3F,i.e.%20multimodal%20or%20multiscale%20information, accessed on 2024-11-21.
#8. S. Parth, “Understanding k-means and hierarchical clustering in r,” 2020. https://medium.com/@shreytparth/understanding-k-means-and-hierarchical-clustering-in-r-9899342cc4fb, accessed on 2024-11-21.
#9. Y. Dai and Y. Xiang, “Connecting sqlite database and shiny app for business intelligence,” 2020. https://shanghai.hosting.nyu.edu/data/r/case-3-sql-shiny.html, accessed on 2024-11-21.
#10. J. Fleming, D. Filer, D. Lloyd, P. Thunga, S. W. Marvel, A. A. Motsinger-Reif,and D. M. Reif, “Manual package ‘toxpir’,” 2024c. https://cran.r-project.org/web/packages/toxpiR/toxpiR.pdf, accessed on 2024-11-21.
#11. S. Marvel and D. M. Reif, “Toxpi™ user manual v2.3,” 2018. https://toxpi.org/dist/ToxPi%20User%20Manual.pdf, accessed on 2024-11-21.
#12. H. Wickham and K. M¨uller, “Dbi: R database interface,” 2024. R package, version 1.2.3, https://cran.r-project.org/web/packages/DBI/index.html.
#13. K. M¨uller and H. Wickham, “Rsqlite: Sqlite interface for r,” 2024. R package, version 2.3.9, https://cran.r-project.org/web/packages/RSQLite/index.html.
#14. H. Wickham, J. Bryan, M. Kalicinski, K. Valery, C. Leitienne, B. Colbert, D. Hoerl, and E. Miller, readxl: Read Excel Files, 2025. R package version 1.4.5, https://cran.r-project.org/web/packages/readxl/index.html.
#15. “The r foundation, the comprehensive r archive network,” 2025. R version 4.4.3,https://cran.r-project.org/bin/windows/base/.
#16. W. Chang, J. Cheng, J. Allaire, C. Sievert, B. Schloerke, Y. Xie, J. Allen,J. McPherson, A. Dipert, and B. Borges, “shiny: Web application framework for r,” 2024. R version 1.10.0, https://cran.r-project.org/web/packages/shiny/index.html.
#17. W. Chang and B. B. Ribeiro, shinydashboard: Create Dashboards with ’Shiny’, 2021. R package version 0.7.2, https://cran.r-project.org/web/packages/shinydashboard/index.html.
#18. C. M. Bishop and N. M. Nasrabadi, Pattern recognition and machine learning, vol. 4. Springer, 2006.
#19. A. Kassambara and F. Mundt, factoextra: Extract and Visualize the Results of Multivariate Data Analyses, 2020. R version 1.0.7, https://cran.r-project.org/web/packages/factoextra/index.html.
#20. B. Ripley, B. Venables, D. M. Bates, K. Hornik, A. Gebhardt, and D. Firth, MASS: Support Functions and Datasets for Venables and Ripley’s MASS, 2025. R package version 7.3-64. https://cran.r-project.org/web/packages/MASS/index.html.
#21. H. Wickham, W. Chang, L. Henry, T. L. Pedersen, K. Takahashi, C. Wilke, K. Woo, H. Yutani, D. Dunnington, and T. van den Brand, ggplot2: Create Elegant Data Visualisations Using the Grammar of Graphics, 2024. R package version 3.5.1, https://cran.r-project.org/web/packages/ggplot2/index.html.
#22. Y. Xie, J. Cheng, X. Tan, J. Allaire, M. Girlich, G. F. Ellis, J. Rauh, S. Limited, B. Reavis, L. Gersen, B. Szopka, A. Pickering, W. Holmes, M. Marttila, A. Quintero, and S. Laurent, DT: A Wrapper of the JavaScript Library ’DataTables’, 2024. R version 0.33, https://cran.r-project.org/web/packages/DT/index.html.
#23. B. Auguie and A. Antonov, gridExtra: Miscellaneous Functions for ”Grid” Graphics, 2017. R package version 2.3, https://cran.r-project.org/web/packages/gridExtra/index.html.
#24. J. G. Speight, “Chapter 2 - organic chemistry,” in Environmental Organic Chemistry for Engineers (J. G. Speight, ed.), pp. 43–86, Butterworth-Heinemann, 2017.
#25. A. Kassambara and F. Mundt, “factoextra : Extract and visualize the results of multivariate data analysesg,” 2020. https://cran.r-project.org/web/packages/factoextra/readme/README.html.
#26. “Appendix c toxicological priority index (toxpi).” National Research Council. A Framework to Guide Selection of Chemical Alternatives. Washington, DC: TheNational Academies Press, 2014. https://doi.org/10.17226/18872, accessed on 2024-11-21.
#27. J. F. Fleming, D. L. Filer, D. T. Lloyd, P. Thunga, S. W. Marvel, A. A. Motsinger-Reif, and D. M. Reif, “toxpir: Create toxpi prioritization models,” 2024d. R package, version 1.3.0, https://cran.r-project.org/web/packages/toxpiR/index.html.
#28. Polycyclic aromatic compounds.” National Toxicology Program, U.S. Department of Health and Human Services, 2014. https://ntp.niehs.nih.gov/whatwestudy/topics/pacs, accessed on 2024-11-21.
#29. J. S. House, F. A. Grimm, W. D. Klaren, A. Dalzell, S. Kuchi, S.-D. Zhang, K. Lenz, P. J. Boogaard, H. B. Ketelslegers, T. W. Gant, and et al., “Grouping of uvcb substances with new approach methodologies (nams) data,” ALTEX -Alternatives to animal experimentation, vol. 38, p. 123–137, Jan. 2021.
#30. F. Sewell, C. Alexander-White, S. Brescia, R. A. Currie, R. Roberts, C. Roper, C. Vickers, C. Westmoreland, and I. Kimber, “New approach methodologies (nams): identifying and overcoming hurdles to accelerated adoption,” Toxicology Research, vol. 13, p. tfae044, 03 2024.
#31. N. Kleinstreuer and T. Hartung, “Artificial intelligence (ai)—it’s the end of the tox as we know it (and i feel fine),” Archives of Toxicology, vol. 98, no. 3, pp. 735–754, 2024.
#32. A. Lai, A. M. Clark, B. I. Escher, M. Fernandez, L. R. McEwen, Z. Tian,Z. Wang, and E. L. Schymanski, “The next frontier of environmental unknowns: substances of unknown or variable composition, complex reaction products, or biological materials (uvcbs),” Environmental Science & Technology, vol. 56, no. 12, pp. 7448–7466, 2022.
#33. J. F. Fleming, J. S. House, J. R. Chappel, A. A. Motsinger-Reif, and D. M.Reif, “Guided optimization of toxpi model weights using a semi-automated approach,” Computational Toxicology, vol. 29, p. 100294, 2024a.
#34. H. Kamp, N. A. Kocabas, F. Faulhammer, N. Synhaeve, E. Rushton, B. Flick, V. Giri, S. Sperber, L. Higgins, M. Penman, et al., “Utility of in vivo metabolomics to support read-across for uvcb substances under reach,” Archives of Toxicology, vol. 98, no. 3, pp. 755–768, 2024.
#35. S. Heux, T. Fuchs, J. Buhmann, N. Zamboni, and U. Sauer, “A high-throughput metabolomics method to predict high concentration cytotoxicity of drugs from low concentration profiles,” Metabolomics, vol. 8, pp. 433–443, 06 2011.
#36. C. B. Pestana, J. W. Firman, and M. T. Cronin, “Incorporating lines of evidence from new approach methodologies (nams) to reduce uncertainties in a category based read-across: a case study for repeated dose toxicity,” Regulatory Toxicology and Pharmacology, vol. 120, p. 104855, 2021.
#37. S. W. Marvel, K. To, F. A. Grimm, F. A. Wright, I. Rusyn, and D. M. Reif, “Toxpi graphical user interface 2.0: Dynamic exploration, visualization, and sharing of integrated data models,” BMC bioinformatics, vol. 19, pp. 1–7, 2018.
#38. J. Fleming, S. W. Marvel, S. Supak, A. A. Motsinger-Reif, and D. M. Reif, “Toxpi* gis toolkit: creating, viewing, and sharing integrative visualizations for geospatial data using arcgis,” Journal of Exposure Science & Environmental Epidemiology, vol. 32, no. 6, pp. 900–907, 2022.
#39. D. M. Reif, M. Sypa, E. F. Lock, F. A. Wright, A. Wilson, T. Cathey, R. R. Judson, and I. Rusyn, “Toxpi gui: an interactive visualization tool for transparent integration of data from diverse sources of evidence,” Bioinformatics, vol. 29, no. 3, pp. 402–403, 2013.
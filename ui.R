source("./Rsource/SwitchButton.R")
source("./Rsource/customMenuSubItem.R")
tagList(
dashboardPage(skin = "black", title="Theseus suite - A collection of software tools for the proteomics community",
  dashboardHeader(disable = F
                  ,title=div(id="title-div", a(img(src="Logo_Theseus_suite_small.png"), href="#", onclick="shinyjs.reset();") #javascript:location.reload();
                              )
                  ),
  dashboardSidebar(
    sidebarMenu(
    id = "tabs",
    menuItem("Welcome", tabName = "welcome", icon = icon("home")),
    menuItem("About", tabName = "about", icon = icon("users")),
    menuItem("Documentation", tabName = "howToUse", icon = icon("book-open")),
    # menuItem("Protein cleaver", tabName = "insilico_digestion", icon = icon("cut"), badgeLabel = "beta", badgeColor = "blue"),
    menuItem("Theseus toolbox", icon = icon("laptop-code"), 
             customMenuSubItem("Protein Cleaver", tabName = "insilico_digestion"),
             customMenuSubItem("SICyLIA TMT", tabName = "datauploadTMT", badgeLabel = "beta", badgeColor = "blue")
             # menuSubItem("SICyLIA LFQ", tabName = "datauploadLFQ")
             ),
    menuItem("Databases", icon = icon("database"),
             customMenuSubItem("C-omics Explorer", href = "http://83.212.98.39/ComicsExplorer/Default", badgeLabel = "beta", badgeColor = "blue")
             ),
    menuItem("Downloads", icon = icon("download"), 
             menuSubItem("Sample data", tabName = "sampleDownload"),
             menuSubItem("Publications", tabName = "manuscriptDownload")),
    menuItem("Disclaimer", tabName = "disclaimer", icon = icon("exclamation-circle")),
    menuItem("How to cite", tabName = "howToCite", icon = icon("feather-alt")),
    hr(),
    menuItemOutput("resultsmenuitem"),
    menuItemOutput("qualitycontrolmenuitem"),
    menuItemOutput("resetsmenuitem")
    #menuItemOutput("graphsmenuitem"),
    #,withSpinner(hidden(menuItemOutput("hidden")))
    #menuItem("Reset session", tabName = "reset", icon = icon("redo"))
    )
    # ,div(id="cruklogodiv", 
    #     h6("Supported by: The Proteomics group"),
    #     img(src="cruk-beatson-logo.jpg", height="15%", width="15%")
    #     )
    ,tags$footer(div(class="footerdiv", "Supported by:", br(), "The ", tags$u(a(href="http://www.beatson.gla.ac.uk/Advanced-Technologies/proteomics.html", target="_blank", "Proteomics Facility"))), id="cruklogodiv", align = "center", a(href="http://www.beatson.gla.ac.uk/", target="_target", img(src="cruk-beatson-logo.jpg", height="80%", width="80%")))
  ),
  dashboardBody(
    fluidPage(theme = shinytheme("yeti"),
    busy_start_up(loader = spin_epic("orbit", color = "#FFF"),
                text = "Initialization...",
                timeout = 1500,
                color = "#FFF",
                background = "#222d32"
              ),
     use_busy_bar(color="#22478a", centered = TRUE, height = "4px"),
     useShinyjs(debug=F),
     # useShinyalert(),
     # extendShinyjs(script = "www/jsCode.js", functions = c("seque", "hidesequeNmolart", "reset", "swalErrorAlert", "swalAlert", "enableTab", "disableTab")),  ### LINUX
     extendShinyjs(script = "jsCode.js", functions = c("seque", "hidesequeNmolart", "reset", "swalErrorAlert", "swalAlert", "enableTab", "disableTab")),  ####### WINDOWS
      tags$head(tags$style("#plot{height:70vh !important;}"),
                tags$style("#plot2{height:83vh !important;}"), 
                tags$script(src = "handlebars.js"),
                # tags$script(src = "sweetalert.min.js"),
                tags$script(src = "sweetalert2.min.js"),
                tags$script(src = "sequence-viewer.min.js"),
                tags$script(src = "molart.js"),
                tags$link(rel = "stylesheet", type = "text/css", href = "custom.css?201117"),
                tags$link(rel = "stylesheet", type = "text/css", href = "button.css?201117"),
                tags$link(rel="shortcut icon", href="Logo_Theseus_suite_ico16x16.png")),
      fluidRow(
        column(1),
        column(10,
      tabItems(
        tabItem(tabName = "welcome", 

                fluidRow(
                  column(12, align="center",
                  # br(),
                  img(src="Logo_Theseus_suite_full.png"),
                  br(),br(),hr()
                  )
                ),
                fluidRow(
                  column(12, align="justify",
                          br(),
                      h4("The ",em("Theseus Suite")," is a collection of auxiliary tools for downstream analysis and visualization of ",a(em("MaxQuant"), href="https://www.maxquant.org/", target="_blank")," output proteomics data in the R Shiny web interface. It provides an automated pipeline for the analysis of amino acid modifications and an integrated framework for the systematic comparison of theoritical and observed proteolitic cleavage sites.")
                      )
                  )
        ),
        
        tabItem(tabName = "howToUse", 
                
                fluidRow(
                  column(12, align="left",
                         h4("How to use Theseus suite"),
                         hr(),
                         br()
                  )
                )
        ),
        
        tabItem(tabName = "howToCite", 
                
                fluidRow(
                  column(12, align="left",
                         h4("How to cite Theseus suite"),
                         hr(),
                         br()
                  )
                )
        ),
        
        tabItem(tabName = "about", 
                
                fluidRow(
                  column(12, align="left",
                         h4("What is Theseus suite"),
                         hr(),
                         br()
                  )
                )
        ),
        
        tabItem(tabName = "disclaimer", 
                
                fluidRow(
                  column(12, align="left",
                         h4("Disclaimer of liability"),
                         hr(),
                         p("The software is provided \"as is\", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and non-infringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software. The authors make no warranties about the accuracy, reliability, completeness, or timeliness of the material, services, software, text, graphics, and links. Description of, or references to, products or publications does not imply endorsement of that product or publication. By using this site you agree to assume all risks associated with your use or transfer of any and all information contained on this site."),
                         br(),br(),
                         h4("Privacy policy"),
                         hr(),
                         p("We collect information that your browser sends whenever you visit the \"Theseus Suite\" website (\"Log Data\"). This Log Data may include information such as your computer's Internet Protocol (\"IP\") address, browser type, browser version, the pages of our website that you visit, the time and date of your visit, the time spent on those pages and other statistics. In addition, we may use third party services such as Google Analytics that collect, monitor and analyze website traffic, in order to track usage statistics and to identify possible operational problems. This information is not used to identify individuals or organizations, and is never shared with third parties. Cookies may be used by the search pages in order to remember your analysis settings. Some cookies persist after you exit the browser, but they are never used for either identification or tracking purposes. Alternatively, you may want to use the standalone version which runs \"Theseus Suite\" locally on your computer.")
                  )
                )
        ),
        
        tabItem(tabName = "sampleDownload", 
                
                fluidRow(
                  column(12, align="left",
                         # h4("Download sample data"),
                         h4("Sample data can be found under this section. Follow our step-by-step quick quide to learn how to use Theseus suite and analyze your own proteomics data."),
                         hr(),
                         br(),
                         
                         tags$table(class="table table-striped",
                                    tags$tr(
                                      tags$th("Description"),
                                      tags$th("Action")
                                    ),
                                    tags$tr(
                                      tags$td("TMT 10plex sample dataset analyzed in MaxQuant"),
                                      tags$td(downloadButton("downloadDummyRec", label = "Download", class = "btn btn-xs"))
                                    )
                         )
                         
                  )
                )
        ),
        
        tabItem(tabName = "insilico_digestion",
                fluidRow(
                  column(3, align="center",
                         img(src="Logo_Protein_cleaver_ico.png")
                  ),
                  column(9, align="center",
                         h4("In silico prediction of protease-induced cleavage sites in protein sequences")
                  )
                ),
                hr(),
                br(),
                tabsetPanel(id="digestion_tabSet",
                  tabPanel(title="Configuration panel", id="digestion_configuration_tab", value="digestion_configuration_tab",
                br(),br(),
                p(tags$u("Step 1:"), " Configure digestion parameters"),
                fluidRow(
                  column(3, align="center",
                           selectInput(
                           inputId="protease",
                           label = "Enzyme",
                           choices = enzyme_list,
                           selected = "trypsin",
                           multiple = FALSE,
                           selectize = FALSE
                         )
                  ),
                  # column(3, align="center",
                  #        sliderInput(inputId="peptide_length", label="Min/max peptide length", min = 1, 
                  #                    max = 35, value = c(7, 25))
                  # ),
                  column(3, align="center",
                         numericInput(inputId="peptide_length", label="Min peptide length", value=7, step=1, width="80%" 
                         )
                  ),
                  column(3, align="center",
                         numericInput(inputId="mol_weight", label="Max peptide mass [Da]", value=4600, step=10, width="80%"
                         )
                  ),
                  column(3, align="center",
                         selectInput(
                           inputId="no_of_misceavages",
                           label = "Allowed miscleavages",
                           choices = c("No miscleavage"="0",
                                      "Up to 1"="1",
                                      "Up to 2"="2"),
                           selected = "0",
                           multiple = FALSE,
                           selectize = FALSE
                         )
                  )
                  
                ),
                br(),br(),
                p(tags$u("Step 2:"), " Upload a fasta file or submit a list of UniProt identifiers"),
                fluidRow(
                  column(3, align="center",
                         switchButton(inputId = "switchButtonDigest",
                                      label = "Select upload method", 
                                      value = TRUE, col = "GB", type = "OO"),
                         br(),br(),
                         actionBttn(inputId="uploadExampleDT",
                                    label = "Example dataset",
                                    icon = NULL,
                                    style = "minimal",
                                    color = "primary",
                                    size = "md",
                                    block = FALSE,
                                    no_outline = TRUE)
                  ),
                  column(6, align="center",
                         fileInput(inputId="fastaFile", 
                                   span("Protein fasta file",
                                        id="fastafileSpan"
                                   ),
                                   multiple=FALSE,
                                   accept = c(
                                     'text/plain',
                                     '.fasta'
                                   )
                         ),
                         textAreaInput(inputId="proteinListTextInput", 
                                   label="Protein list", 
                                   value = "", 
                                   width = "100%", 
                                   height = "200px",
                                   placeholder = "Paste a list of UniProt accession IDs or try the example data set",
                                   resize = "vertical")
                     ),
                  column(3, align="center",
                         selectInput(
                           inputId="fastaformat",
                           label = "Parsing rule",
                           choices = c("UniProt identifier"="uniprot"),
                           selected = "uniprot",
                           multiple = FALSE,
                           selectize = FALSE
                         ),
                         br(),br(),
                         actionBttn(inputId="submitProtList",
                                    label = "Submit",
                                    icon = NULL,
                                    style = "pill",
                                    color = "primary",
                                    size = "md",
                                    block = FALSE,
                                    no_outline = TRUE)
                  )
                )
                ),
                tabPanel(title="Context browser", id="digestion_results_tab", value="digestion_results_tab",
                   br(),
                   div(id="prot-viewer-div",
                       # box(width = 12, title = "Selected parameters", solidHeader = TRUE, collapsible = FALSE,
                           fluidRow(
                            column(3, align="center", h6(textOutput("selected_enzyme"))),
                            column(3, align="center", h6(textOutput("minimum_pep_len"))),
                            column(3, align="center", h6(textOutput("maximum_pep_mass"))),
                            column(3, align="center", h6(textOutput("missed_cleavages")))
                      ),
                      # ),
                      br(),
                      fluidRow(id="fluidRow-sequence-viewer",
                        column(12, align="left",
                              div(id="sequence-viewer"),
                              
                        )
                      ),
                      br(),
                      fluidRow(id="fluidRow-molart-viewer",
                        column(12, align="left",
                              div(id="molart-viewer")
                        ) 
                      )
                    ),
                  
                  br(),
                  div(id="identifiable-prot-pept-datatable",
                    verticalTabsetPanel(id="vertical_results_tabSet", contentWidth = 10, color = "#22478a",
                      verticalTabPanel(id="vertical_entries_tab", value="vertical_entries_tab",
                        title = "Data & Reports", box_height = "100px;", icon = icon("table"),
                        tabsetPanel(
                          tabPanel("Identifiable peptides", icon = icon("table"), br(), DT::dataTableOutput("identifiable_pept"),
                                   fluidRow(column(12, align="right", uiOutput("download_identifiable_pept"))), br() ),
                          tabPanel("Identifiable proteins", icon = icon("table"), br(), DT::dataTableOutput("identifiable_prot"), br() ),
                          tabPanel("Non-detectable proteins", icon = icon("table"), br(), DT::dataTableOutput("non_identifiable_prot"), br() ),
                          # tabPanel("Common identifiable peptides", icon = icon("table"), br(), DT::dataTableOutput("frequent_pept"), br() ),
                          tabPanel("Cleaved peptides frequency", icon = icon("chart-bar"),
                                   br(), 
                                   fluidRow(
                                     column(12,
                                            align="center", 
                                            plotOutput("peptidesDist")
                                     )
                                   ) , br() 
                          ),
                          tabPanel("Cleaved peptides coverage", icon = icon("chart-bar"),
                                   br(), 
                                   fluidRow(
                                     column(12,
                                            align="center", 
                                            plotOutput("peptidesSeqCoverage")
                                     )
                                   ) , br() 
                          )
                        )
                      ),
                      verticalTabPanel(id="vertical_seq_coverage_estim", value="vertical_seq_coverage_estim",
                                       title = "Coverage comparison", box_height = "100px;", icon = icon("calculator"),
                                       tabsetPanel(
                                         tabPanel("Data panel", icon = icon("table"),
                                                  br(), 
                                                  fluidRow(
                                                    column(12,
                                                           align="center",
                                                           fileInput(inputId="protGroupsDigestion", 
                                                                     span("ProteinGroups.txt file",
                                                                          id="fastafileSpan"
                                                                     ),
                                                                     multiple=FALSE,
                                                                     accept = c(
                                                                       'text/plain',
                                                                       '.txt'
                                                                     )
                                                           ),
                                                           br(), 
                                                           DT::dataTableOutput("protGroups_digestion")
                                                    )
                                                  ), br()
                                         ),
                                         tabPanel("Scatter plot", icon = icon("chart-bar"),
                                                  br(),
                                                  fluidRow(
                                                    column(12,
                                                           align="center",
                                                           plotOutput("protGroups_theorVSobserv")
                                                    )
                                                  ), br()
                                           
                                         )
                                       )
                      ),
                      verticalTabPanel(id="vertical_bulk_digest_tab", value="vertical_bulk_digest_tab",
                        title = "Bulk digestion", box_height = "100px;", icon = icon("cut"),
                        tabsetPanel(
                          tabPanel("Summary", icon = icon("table"),
                                 br(), 
                                 fluidRow(
                                   column(12,
                                          align="center",
                                          DT::dataTableOutput("bulk_digestDT")
                                   )
                                 ), br()
                        )
                      )
                      )
                      
                      ) # verticalTabsetPanel
                    ) # div
                  ) # tabPanel
                ), # tabsetPanel
                br()
          
        ),#end of tabItem
        
          tabItem(tabName = "datauploadTMT",
                  fluidRow(
                    column(3, align="center",
                           img(src="Logo_SICyLIA-TMT_ico.png")
                    ),
                    column(9, align="center",
                           h4("A computational workflow for the analysis of cysteine oxidation from Tandem Mass Tag (TMT) proteomics datasets of various complexity.")
                    )
                  ),
                hr(),
                br(),
                p(tags$u("Step 1:"), " Select the multiplexing capacity of your experiment."),
                #radioButtons("plex_radiobtn","Plex:", c("10-plex" = "10", "11-plex" = "11", "16-plex" = "16"), inline=T, selected = "10"),
                
                fluidRow(
                  column(12, align="center",
                    awesomeRadio(inputId = "Id049", 
                             label = NULL,
                             choices = c("10-plex" = "10", "11-plex" = "11", "16-plex" = "16"),
                             inline = T,
                             checkbox = T,
                             selected = "10"),
                
                  )),
                hr(),
                br(),

                p(tags$u("Step 2:"), " Upload the required MaxQuant output files. Upon successful upload and validation the analysis will run automatically."),
                #br(),
                
                fluidRow(
                  column(6, align="center",
                
                fileInput('peptides', 
                          span("Peptides  ",
                               id="peptidesSpan",
                               tags$a(
                                 tags$i(class='fa fa-question-circle'),
                                 href = "#",
                                 onclick = "shinyjs.swalAlert('<h3><strong>Required fields</strong></h3>', '<h4>The following columns must exist in the uploaded <strong>peptides.txt</strong> file:</h4> <br /> <h5><strong>id</strong>  <br />  <strong>Leading razor protein</strong>  <br />  <strong>Start position</strong>  <br />  <strong>End position</strong>  <br />  <strong>C Count</strong>  <br />  <strong>Length</strong></h5> <br />', 'info'); return false;")
                          ),
                          multiple=FALSE,
                          accept = c(
                            'text/csv',
                            'text/comma-separated-values',
                            'text/tab-separated-values',
                            'text/plain',
                            '.csv',
                            '.tsv'
                          )
                ),

                br(),
                
                fileInput('proteinGroups', 
                          span("Protein groups  ",
                               id="proteinGroupsSpan",
                               tags$a(
                                 tags$i(class='fa fa-question-circle'),
                                 href = "#",
                                 onclick = "shinyjs.swalAlert('<h3><strong>Required fields</strong></h3>', '<h4>The following columns must exist in the uploaded <strong>proteinGroups.txt</strong> file:</h4> <br /> <h5><strong>Reporter intensity corrected 0</strong>  <br /> <strong>Reporter intensity corrected 1</strong> <br /> <strong>Reporter intensity corrected 2</strong> <br /> <strong>Reporter intensity corrected 3</strong> <br /> <strong>Reporter intensity corrected 4</strong> <br /> <strong>Reporter intensity corrected 5</strong> <br /> <strong>Reporter intensity corrected 6</strong> <br /> <strong>Reporter intensity corrected 7</strong> <br /> <strong>Reporter intensity corrected 8</strong> <br /> <strong>Reporter intensity corrected 9</strong> <br /> <strong>Mod. peptide IDs</strong> <br /> <strong>Unique peptides</strong></h5> <br />', 'info'); return false;")
                          ),
                          multiple=FALSE,
                          accept = c(
                            'text/csv',
                            'text/comma-separated-values',
                            'text/tab-separated-values',
                            'text/plain',
                            '.csv',
                            '.tsv'
                          )
                ),
                
                br(),
                
                fileInput('modSpecPeptides', 
                          span("Modification specific peptides  ", 
                               id="modSpecPeptidesSpan",
                               tags$a(
                                 tags$i(class='fa fa-question-circle'),
                                 href = "#",
                                 onclick = "shinyjs.swalAlert('<h3><strong>Required fields</strong></h3>', '<h4>The following columns must exist in the uploaded <strong>modificationSpecificPeptides.txt</strong> file:</h4> <br /> <h5><strong>Reporter intensity corrected 0</strong>  <br /> <strong>Reporter intensity corrected 1</strong> <br /> <strong>Reporter intensity corrected 2</strong> <br /> <strong>Reporter intensity corrected 3</strong> <br /> <strong>Reporter intensity corrected 4</strong> <br /> <strong>Reporter intensity corrected 5</strong> <br /> <strong>Reporter intensity corrected 6</strong> <br /> <strong>Reporter intensity corrected 7</strong> <br /> <strong>Reporter intensity corrected 8</strong> <br /> <strong>Reporter intensity corrected 9</strong> <br /> <strong>Modifications</strong> <br /> <strong>Raw file</strong> <br /> <strong>Reverse</strong> <br /> <strong>Potential contaminant</strong> <br /> <strong>Unique  Groups</strong> <br /> <strong>Unique  Proteins</strong> <br /> <strong>Sequence</strong> <br /> <strong>Proteins</strong> <br /> <strong>Gene Names</strong> <br /> <strong>Protein Names</strong> <br /> <strong>Acetyl  Protein N term</strong> <br /> <strong>Carbamidomethyl  C</strong> <br /> <strong>Carbamidomethyl Heavy  C  Var</strong> <br /> <strong>Oxidation  M</strong> <br /> <strong>id</strong> <br /> <strong>Peptide ID</strong> <br /> <strong>Evidence IDs</strong> <br /> <strong>Carbamidomethyl  C  site IDs</strong> <br /> <strong>Carbamidomethyl Heavy  C  Var  site IDs</strong> <br /> <strong>Mass</strong> <br /> <strong>Retention time</strong> <br /> <strong> PEP </strong> <br /> <strong> Score </strong> <br /> <strong> Delta score </strong> <br /> <strong> Intensity </strong> <br /> <strong> MS MS Count </strong> <br /> <strong> Protein group IDs </strong> <br /> <strong> Protein Groups </strong> <br /> <strong> Missed cleavages </strong> <br /> <strong> Charges </strong></h5> <br />', 'info'); return false;")
                               ),
                          multiple=FALSE,
                          accept = c(
                            'text/csv',
                            'text/comma-separated-values',
                            'text/tab-separated-values',
                            'text/plain',
                            '.csv',
                            '.tsv'
                          )
                )
                
                  
                )
                
                ,column(6, align="center",
                        
                         
                        fileInput('lightPeptides', 
                                  span("Light-labeled sites  ",
                                       id="lightPeptidesSpan",
                                       tags$a(
                                         tags$i(class='fa fa-question-circle'),
                                         href = "#",
                                         onclick = "shinyjs.swalAlert('<h3><strong>Required fields</strong></h3>', '<h4>The following columns are required in the <strong>Carbamidomethyl (C)Sites.txt</strong> file:</h4> <br />  <h5><strong>id</strong>  <br />  <strong>Proteins</strong>  <br />  <strong>Positions within proteins</strong>  <br />  <strong>Positions</strong>  <br />  <strong>Position</strong>  <br />  <strong>Position in peptide</strong>  <br />  <strong>Score diff</strong>  <br />  <strong>Localization prob</strong></h5> <br />', 'info'); return false;")
                                  ),
                                   multiple=FALSE,
                                   accept = c(
                                     'text/csv',
                                     'text/comma-separated-values',
                                     'text/tab-separated-values',
                                     'text/plain',
                                     '.csv',
                                     '.tsv'
                                   )
                         ),
                        
                        br(),
                        
                        fileInput('heavyPeptides', 
                                  span("Heavy-labeled sites  ",
                                       id="heavyPeptidesSpan",
                                       tags$a(
                                         tags$i(class='fa fa-question-circle'),
                                         href = "#",
                                         onclick = "shinyjs.swalAlert('<h3><strong>Required fields</strong></h3>', '<h4>The following columns are required in the <strong>Carbamidomethyl Heavy (C)Sites</strong> file:</h4> <br />  <h5><strong>id</strong>  <br />  <strong>Proteins</strong>  <br />  <strong>Positions within proteins</strong>  <br />  <strong>Positions</strong>  <br />  <strong>Position</strong>  <br />  <strong>Position in peptide</strong>  <br />  <strong>Score diff</strong>  <br />  <strong>Localization prob</strong></h5> <br />', 'info'); return false;")
                                  ),
                                  multiple=FALSE,
                                  accept = c(
                                    'text/csv',
                                    'text/comma-separated-values',
                                    'text/tab-separated-values',
                                    'text/plain',
                                    '.csv',
                                    '.tsv'
                                  )
                        )
                         
                )
            )

        ),
        
        # tabItem(tabName = "datauploadLFQ",
        #         # h4("SICyLIA LFQ"),
        #         # hr(),
        #         br(),br(),
        #         fluidRow(
        #           column(12, align="center",
        #                 img(src="coming_soon.png") 
        #           )
        #         )
        #     
        # ),

        tabItem(tabName = "resultsmenuitem",
                tabsetPanel(
                    
                    tabPanel("Unquantified proteins", br(), DT::dataTableOutput("table_1"), 
                             fluidRow(
                               column(2, align="center",
                                  uiOutput("downloadBtn1")
                               ),
                               column(10, align="left",
                                  uiOutput("downloadOptns1")
                               )
                             )
                             ),
                    tabPanel("Non-unique peptides", br(), DT::dataTableOutput("table_2"), 
                             fluidRow(
                               column(2, align="center",
                                      uiOutput("downloadBtn2")
                               ),
                               column(10, align="left",
                                      uiOutput("downloadOptns2")
                               )
                             )),
                    tabPanel("Final data", br(), DT::dataTableOutput("table_3"), 
                             fluidRow(
                               column(2, align="center",
                                      uiOutput("downloadBtn3")
                               ),
                               column(10, align="left",
                                      uiOutput("downloadOptns3")
                               )
                             ))
                )
        ),
        
        tabItem(tabName = "qualitycontrolmenuitem",
                 tabsetPanel(
                   tabPanel("QC plot 1", br(),
                            fluidRow(
                              column(6, align="center",
                                     selectInput(inputId = "rep_intensities_x", label = "Select Reporter Intenisty", choices = NULL)
                              ),
                              column(6, align="center",
                                     selectInput(inputId = "rep_intensities_y", label = "Select Background Intenisty", choices = NULL)
                              )
                            ),
                            br(),
                            fluidRow(
                              #column(1),
                              column(12, align="center", plotOutput("plot"))
                              #column(1)
                            )
                            # tableOutput("rSquare")
                            
                   ),
                   tabPanel("QC plot 2", br(),
                            fluidRow(
                              #column(1),
                              column(12, align="center", plotOutput("plot2"))
                              #column(1)
                            )
                            # tableOutput("rSquare")
                            
                   )
                 )
        )
        
        
        
        
      )
      ),

      column(1)
    )
    )
  )

) #end dashboardPage
)

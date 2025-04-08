#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
#library(plotly)
library(DT)
library(patchwork)
library(latex2exp)
#library(RColorBrewer)
library(shinybrowser)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # detect window size #
  tags$head(tags$script('
                                  var dimension = [0, 0];
                                  $(document).on("shiny:connected", function(e) {
                                      dimension[0] = window.innerWidth;
                                      dimension[1] = window.innerHeight;
                                      Shiny.onInputChange("dimension", dimension);
                                  });
                                  $(window).resize(function(e) {
                                      dimension[0] = window.innerWidth;
                                      dimension[1] = window.innerHeight;
                                      Shiny.onInputChange("dimension", dimension);
                                  });
                              ')),
  
  # Application title
  titlePanel("Power calculation for cross-sectional stepped-wedge cluster randomized trials with time-to-event endpoints"),
  
  sidebarLayout(
    sidebarPanel(
      # cluster or power #
      radioButtons("n_power",
                   "Output display:",
                   c("Power"="power",
                     "Number of clusters (n)"="n")
      ),
      # balanced design or upload own #
      conditionalPanel(
        condition="input.n_power == 'power'",
        radioButtons("balanced",
                     "Design type:",
                     c("Balanced"="balanced",
                       "Unbalanced (upload your own design)"="unbalanced")
        ),
        helpText("Balanced designs refer to randomization schemes where an equal number of clusters are randomized to each treatment sequence.", style="margin-top:-1em; margin-bottom:1em;"),
        
        
        ## upload own design ##
        conditionalPanel(
          condition = "input.balanced == 'unbalanced'",
          fileInput("file1", "Upload a design matrix:", accept=c('text/plain', '.csv')),
          helpText("The file must be a comma separated .csv file consisting of 0s (control) and 1s (treatment), with a column for each time period and a row for each cluster. There may not be any missing cluster-periods. Do not include row or column names.", style="margin-top:-1em; margin-bottom:1em;")
        ),
        
        tags$u(h3("Design constraints")),
        conditionalPanel(
          condition = "input.balanced == 'balanced'",
          # clusters #
          numericInput("n",
                       "Number of clusters (n):",
                       min=2,
                       max=1000,
                       value=2),
          helpText("Total number of clusters to be randomized. For balanced designs, n should be a multiple of number of sequences (J-1).", style="margin-top:-1em; margin-bottom:1em;"),
        ) #end balanced conditional
      ), #end power conditional
      
      conditionalPanel(
        condition = "input.n_power == 'n'",
        tags$u(h3("Design constraints"))
      ),
      # cluster size #
      numericInput("m",
                   "Cluster-period size (m):",
                   min=2,
                   max=100000,
                   value=2),
      helpText("Cluster-period size is the number of participants recruited per cluster in one time period. Power calculations assume cluster-period size is equal across clusters and time.", style="margin-top:-1em; margin-bottom:1em;"),
      
      # time periods #
      conditionalPanel(
        condition = "input.balanced == 'balanced' | input.n_power == 'n'",
        numericInput("J",
                     "Number of time periods (J):",
                     min = 3,
                     max = 50,
                     value = 3),
        ),
      conditionalPanel(
        condition = "input.balanced == 'balanced' & input.n_power == 'power'",
        helpText("The number of time periods will be 1 larger than the number of sequences or 'steps' clusters can be randomized to. For balanced designs, n should be a multiple of number of sequences (J-1).", style="margin-top:-1em; margin-bottom:1em;")
      ),
      conditionalPanel(
        condition = "input.balanced == 'unbalanced' | input.n_power == 'n'",
        helpText("The number of time periods will be 1 larger than the number of sequences or 'steps' clusters can be randomized to.", style="margin-top:-1em; margin-bottom:1em;")
      ),
      
      conditionalPanel(
        condition="input.n_power == 'n'",
        # clusters #
        numericInput("power",
                     "Power:",
                     min=0,
                     max=1,
                     value=0.8,
                     step=0.01)
      ),
      # treatment effect #
      numericInput("beta",
                   "Treatment effect size",
                   min=0,
                   max=999999999,
                   value=1),
      helpText("Treatment effect is measured as a log hazard ratio.", style="margin-top:-1em; margin-bottom:1em;"),
      # correlations or ICC? #
      # radioButtons("icc_tau",
      #              "Specify within- and between-period Kendall's tau or ICCs?",
      #              c("Kendall's tau" = "tau",
      #                "ICC" = "icc")),
      
      # correlations #
      tags$u(h3("Correlation structure")),
      # conditionalPanel(
      #   condition = "input.icc_tau == 'tau'",
      numericInput("tau_w",
                   withMathJax("Within-period Kendall's tau (\\(\\tau_w\\))"),
                   min=0,
                   max=1,
                   value=0.1,
                   step=0.01),
      numericInput("tau_b",
                   withMathJax("Between-period Kendall's tau (\\(\\tau_b\\))"),
                   min=0,
                   max=1,
                   value=0.1,
                   step=0.01),
      helpText("Kendall's tau is a type of rank correlation. The within-period Kendall's tau specifies how related survival times within the same cluster and period are. The within-period Kendall's tau specifies how related survival times within the same cluster but in different periods are.", style="margin-top:-1em; margin-bottom:1em;"),
      # ),
      # icc #
      # conditionalPanel(
      #   condition = "input.icc_tau == 'icc'",
      #   numericInput("icc_w",
      #                withMathJax("Within-period ICC (\\(\\rho_w\\))"),
      #                min=0,
      #                max=1,
      #                value=0.1,
      #                step=0.01),
      #   numericInput("icc_b",
      #                withMathJax("Between-period ICC (\\(\\rho_b\\))"),
      #                min=0,
      #                max=1,
      #                value=0.1,
      #                step=0.01)
      # ),
      
      tags$u(h3("Censoring & event rate constraints")),
      # admin censoring #
      numericInput("admin_censor",
                   "Administrative censoring (proportion)",
                   min=0,
                   max=1,
                   value=0,
                   step=0.01),
      
      # baseline hazard #
      radioButtons("constant_baseline",
                   "Baseline hazard:",
                   c("Constant"="constant",
                     "Change by constant over time"="change_t")
      ),
      helpText(withMathJax("Constant baseline hazard: \\(\\lambda_{0j}(t) = \\lambda_0(t)\\); Baseline hazard changing by constant \\(C\\) over time: \\(\\lambda_{0j}(t) = \\lambda_0(t) + C(j-1)\\)"), style="margin-top:-1em; margin-bottom:1em;"),
      
      conditionalPanel(
        condition="input.constant_baseline == 'change_t'",
        numericInput("baseline_change",
                     "Baseline hazard change constant:",
                     min=-99999,
                     max=99999,
                     value=0)
      ),
      # numericInput("other_censor",
      #              "Other censoring (proportion)",
      #              min=0,
      #              max=1,
      #              value=0,
      #              step=0.01),
      conditionalPanel(
        condition = "input.n_power == 'power'",
        tags$u(h3("Type I error & degrees of freedom"))
      ),
      
      conditionalPanel(
        condition = "input.n_power == 'n'",
        tags$u(h3("Type I error"))
      ),
      
      # sig level #
      numericInput(inputId="sig",
                   label="Significance level",
                   0.05, min=0.000001, max=1, step=0.001),
      
      # degrees of freedom #
      conditionalPanel(
        condition = "input.n_power == 'power'",
        radioButtons("DoF",
                     "Degrees of Freedom:",
                     c("(n-1)"="n1",
                       "(n-2)"="n2")
                     )
      ),

      ## loading message when the app is calculating ##
      tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 8%;
               left: 33.5%;
               width: 50%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 150%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          ")),
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                       tags$div("Loading...",id="loadmessage")),
      
      actionButton("submit", "Update View")
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Results",
                  span( textOutput("text_context"), 
                        tags$ul(
                          tags$li( textOutput("text_t") ),
                          tags$li( textOutput("text_smScore") ),
                          tags$li( textOutput("text_tangScore") )
                        ),
                        textOutput("text_ICCs"),
                        
                        style="font-size:20px" )
                  ),
        
        tabPanel("Design Matrix",
                 tableOutput("design_matrix"),
                 conditionalPanel(
                   condition="input.n_power == 'n'",
                   print("*Calculations are made assuming total number of clusters calculated in 'Results' tab are evenly distributed to each of the above sequences.")
                 ),
                 
      ),
      tabPanel("References and Resources", 
               h3("References and Resources"),
               HTML("<p style = 'font-size:13pt;'>This application has been written by Mary Ryan Baumann (University of Wisconsin - Madison USA), Fan Li (Yale School of Public Health USA), and with input from Monica Taljaard (University of Ottawa - Ottawa Hospital Research Institute CA). Please email <a href='mailto:mary.ryan@wisc.edu'>mary.ryan@wisc.edu</a> if you need to report errors or would like to submit comments or feedback.</p>"),
               p("The code repository for this application can be found at: ", tags$a(href="https://github.com/maryryan/survivalSWCRT", "https://github.com/maryryan/survivalSWCRT"), style = "font-size:13pt;"),
               p("A preprint detailing the methods this application uses can be found at: ", tags$a(href="https://doi.org/10.48550/arXiv.2312.13097", "https://doi.org/10.48550/arXiv.2312.13097"), style = "font-size:13pt;")
      )
        
    )#end tablesetPanel
  )#end mainPanel
  )#end sibebar layout
))

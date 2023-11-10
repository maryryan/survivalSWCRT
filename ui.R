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
  titlePanel("Power calculation for stepped-wedge cluster randomized trials with time-to-event endpoints"),
  
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
                       "Unbalanced"="upload")
        ),
        
        ## upload own design ##
        conditionalPanel(
          condition = "input.balanced == 'upload'",
          fileInput("file1", "Upload a design matrix:", accept=c('text/plain', '.csv')),
          helpText("The file must be a comma separated .csv file consisting of 0s (control) and 1s (treatment), with a column for each time period and a row for each cluster. There may not be any missing cluster-periods. Do not include row or column names.", style="margin-top:-0.5em; margin-bottom:1em;")
        ),
        
        conditionalPanel(
          condition = "input.balanced == 'balanced'",
          # time periods #
          numericInput("J",
                       "Number of time periods:",
                       min = 3,
                       max = 50,
                       value = 3),
          # clusters #
          numericInput("n",
                       "Number of clusters:",
                       min=2,
                       max=1000,
                       value=2)
        ) #end balanced conditional
      ), #end power conditional
      
      conditionalPanel(
        condition = "input.n_power == 'n'",
        # time periods #
        numericInput("J",
                     "Number of time periods:",
                     min = 3,
                     max = 50,
                     value = 3)
      ),
      
      # cluster size #
      numericInput("m",
                   "Cluster-period size:",
                   min=2,
                   max=100000,
                   value=2),
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
      # correlations or ICC? #
      # radioButtons("icc_tau",
      #              "Specify within- and between-period Kendall's tau or ICCs?",
      #              c("Kendall's tau" = "tau",
      #                "ICC" = "icc")),
      # correlations #
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
      
      # baseline hazard #
      radioButtons("constant_baseline",
                   "Baseline hazard:",
                   c("Constant"="constant",
                     "Change by constant over time"="change_t")
      ),
      conditionalPanel(
        condition="input.constant_baseline == 'change_t'",
        numericInput("baseline_change",
                     "Baseline hazard change constant:",
                     min=-99999,
                     max=99999,
                     value=0)
      ),
      # censoring #
      # numericInput("admin_censor",
      #              "Administrative censoring (proportion)",
      #              min=0,
      #              max=1,
      #              value=0,
      #              step=0.01),
      # numericInput("other_censor",
      #              "Other censoring (proportion)",
      #              min=0,
      #              max=1,
      #              value=0,
      #              step=0.01),
      # sig level #
      numericInput(inputId="sig",
                   label="Significance level",
                   0.05, min=0.000001, max=1, step=0.001)
      
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
      tabPanel("Design Matrix",#br(),
               tableOutput("design_matrix")
               )
    )#end tablesetPanel
  )#end mainPanel
  )#end sibebar layout
))

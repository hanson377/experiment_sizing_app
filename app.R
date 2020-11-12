## in this file, I generate a shiny app that allows users to visualize two proportions with a beta distribution
## additionally, the distribution of the % difference between the two proportions from 1m simulations is visualized

library(shiny)
library(dplyr)
library(knitr)
library(kableExtra)
library(ggplot2)

calc_stats <- function(data) {

 trials1 <- data$n_control
 alpha1 <- round(data$n_control*data$baseline,digits=0)
 beta1 <- trials1-alpha1

 trials2 <- data$n_variant
 alpha2 <- round(data$n_variant*data$variant,digits=0)
 beta2 <- trials2-alpha2

 sim_trials <- 1000000

 sample1 <- rbeta(sim_trials,alpha1,beta1)
 sample2 <- rbeta(sim_trials,alpha2,beta2)

 prob <- sum(sample2 >= sample1)/sim_trials

 ## calculate quantiles for differences from simulations
 sample <- data.frame(control = sample1, variant = sample2)
 sample$diff <- (sample$variant/sample$control)-1

 median_diff <- quantile(sample$diff,.5)
 lower_diff <- quantile(sample$diff,.01)
 upper_diff <- quantile(sample$diff,.99)

 ##
 summary_stats <- data.frame(prob,median_diff,lower_diff,upper_diff)

 return(summary_stats)
}


# Define UI ----
ui <- fluidPage(
  titlePanel("Bayesian Comparision of Proportions"),
  sidebarLayout(
    sidebarPanel(
      h1('Experiment Sizing Inputs'),
      ##
      h3('Proportion Comparision Inputs'),
      sliderInput("baseline", "Baseline Proportion:",min = 0, max = 1, value = 0.45),
      numericInput("min_delta", "Minimum Relative Lift:",min = 0, max = 1, value = 0.0),
      numericInput("max_delta", "Maximum Relative Lift:",min = 0, max = 1, value = 0.1),
      numericInput("incremental_differences", "Increments Between Min and Max",min = .001, max = .1, value = 0.02),
      numericInput("min_effect_desized", "Minimum Effect Desired",min = .001, max = .1, value = 0.01),
      ##
      h3('Sample Size Inputs'),
      numericInput("n_per_day", "Volume per Day:", 10000, min = 1, max = 100000),
      sliderInput("sample_split", "Sample Size Split:",min = 0.01, max = .99, value = 0.5),
      ##
      h3('Runtime Inputs'),
      sliderInput("max_days", "Maxium Runtime",min = 1, max = 28, value = 14),
      sliderInput("day_increments", "Day Increments",min = 1, max = 7, value = 7)
    ),

    mainPanel(
      h1('Summary Views of Time to Run'),
      ##splitLayout(cellWidths = c("30%", "30%", "40%"), plotOutput("thing"), plotOutput("thing"), plotOutput("thing")),
      plotOutput("thing"),
      tableOutput("table_summary")
    )
  )
)

# Define server logic ----
server <- function(input, output,session) {

  base <- reactive({

    ## ## create skeleton frame for visuals
    skeleton <- data.frame(baseline=input$baseline,delta=seq(input$min_delta,input$max_delta,input$incremental_differences),binary=1)
    skeleton$variant <- skeleton$baseline*(1+skeleton$delta)

    ## create sample size skeleton to cross join to skeleton
    skeleton2 <- data.frame(days = seq(1,input$max_days,input$day_increments),binary=1)
    skeleton2$days <- skeleton2$days-1+input$day_increments
    skeleton2$n_total <- skeleton2$days*input$n_per_day
    skeleton2$n_control <- round(skeleton2$n_total*input$sample_split,digits=0)
    skeleton2$n_variant <- skeleton2$n_total-skeleton2$n_control
    ## cross join!
    base <- skeleton %>% inner_join(skeleton2,by='binary')
    base <- base %>% select(baseline,variant,delta,days,n_total,n_control,n_variant) %>% arrange(days,delta) %>% mutate(row = row_number())

    ## run simulations to calculate what we would know given each of these hypothetical deltas and sample sizes
    unique_rows <- unique(base$row)

    results <- NA

    for (i in unique_rows) {

    base_filtered <- subset(base,row==i)
    stats <- calc_stats(base_filtered)

    temp <- data.frame(base_filtered,stats)
    results <- rbind(results,temp)

    }
    results <- subset(results,is.na(prob) == FALSE)

  })

  summary_table <- reactive({

    ## ## create skeleton frame for visuals
    skeleton <- data.frame(baseline=input$baseline,delta=seq(input$min_delta,input$max_delta,input$incremental_differences),binary=1)
    skeleton$variant <- skeleton$baseline*(1+skeleton$delta)

    ## create sample size skeleton to cross join to skeleton
    skeleton2 <- data.frame(days = seq(1,input$max_days,input$day_increments),binary=1)
    skeleton2$days <- skeleton2$days-1+input$day_increments
    skeleton2$n_total <- skeleton2$days*input$n_per_day
    skeleton2$n_control <- round(skeleton2$n_total*input$sample_split,digits=0)
    skeleton2$n_variant <- skeleton2$n_total-skeleton2$n_control
    ## cross join!
    base <- skeleton %>% inner_join(skeleton2,by='binary')
    base <- base %>% select(baseline,variant,delta,days,n_total,n_control,n_variant) %>% arrange(days,delta) %>% mutate(row = row_number())

    unique_rows <- unique(base$row)

    results <- NA

    for (i in unique_rows) {

    base_filtered <- subset(base,row==i)
    stats <- calc_stats(base_filtered)

    temp <- data.frame(base_filtered,stats)
    results <- rbind(results,temp)

    }
    results <- subset(results,is.na(prob) == FALSE)

    ## calculate observations closest to 99% mark
    final <- results %>% mutate(diff = abs(prob-.99)) %>% group_by(factor(days)) %>% mutate(diff_ranker = rank(diff)) %>% filter(diff_ranker == 1)

  })


  output$thing<-renderPlot({
    ggplot(base(),aes(x=delta,y=prob,colour=factor(days))) + geom_line() + theme(legend.position='bottom',legend.title=element_blank()) + scale_y_continuous(label=scales::percent) + xlab('Observed Delta') + ylab('Probability Variant > Control') + geom_hline(yintercept=.99,colour='red',linetype='dashed') + scale_x_continuous(label=scales::percent)
  })

  output$table_summary<-renderTable({
    ##kable(summary_data(),format = 'html', col.names = c('lower diff','upper diff'),caption = 'summary of two proportions')
    summary_table()
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)

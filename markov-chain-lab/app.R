library(shiny)

options(shiny.maxRequestSize = 10 * 1024^2)

examples <- list(
  "All states recurrent (finite irreducible)" = list(
    P = matrix(c(.55,.45,0, .20,.55,.25, 0,.35,.65), 3, byrow=TRUE),
    note = "Every state communicates with every other state. A finite irreducible chain has only positive recurrent states."
  ),
  "One transient state" = list(
    P = matrix(c(.35,.65,0, 0,.65,.35, 0,.25,.75), 3, byrow=TRUE),
    note = "State 1 can leave for the closed class {2,3}, but can never be reached from it. State 1 is transient; states 2 and 3 are recurrent."
  ),
  "Two closed classes and one transient state" = list(
    P = matrix(c(.40,.30,.30, 0,1,0, 0,0,1), 3, byrow=TRUE),
    note = "States 2 and 3 are absorbing recurrent states. State 1 eventually enters one of them and is transient."
  ),
  "Periodic but recurrent" = list(
    P = matrix(c(0,1,0, 0,0,1, 1,0,0), 3, byrow=TRUE),
    note = "The chain returns to each state every three steps. Periodicity does not imply transience."
  )
)

classify_states <- function(P) {
  n <- nrow(P); reach <- P > 1e-12; diag(reach) <- TRUE
  for (k in seq_len(n)) for (i in seq_len(n)) for (j in seq_len(n))
    reach[i,j] <- reach[i,j] || (reach[i,k] && reach[k,j])
  recurrent <- logical(n)
  for (i in seq_len(n)) {
    cls <- which(reach[i,] & reach[,i])
    recurrent[i] <- all(rowSums(P[cls, setdiff(seq_len(n), cls), drop=FALSE]) < 1e-10)
  }
  ifelse(recurrent, "Recurrent", "Transient")
}

simulate_chain <- function(P, start, steps) {
  x <- integer(steps + 1); x[1] <- start
  for (t in seq_len(steps)) x[t+1] <- sample.int(nrow(P), 1, prob=P[x[t],])
  x
}

first_return_curve <- function(P, state, horizon) {
  n <- nrow(P); Q <- P; Q[,state] <- 0
  v <- numeric(n); v[state] <- 1; out <- numeric(horizon)
  for (k in seq_len(horizon)) {
    hit <- sum(v * P[,state]); out[k] <- hit + if (k > 1) out[k-1] else 0
    v <- as.numeric(v %*% Q)
  }
  pmin(out, 1)
}

ui <- fluidPage(
  tags$head(tags$style(HTML("body{background:#f7f8fc}.well{background:white;border:0;box-shadow:0 2px 12px #dde1eb}.hero{padding:22px 26px;background:linear-gradient(120deg,#172554,#2563eb);color:white;border-radius:12px;margin:18px 0}.hero h2{margin-top:0}.value{font-size:28px;font-weight:700;color:#1d4ed8}.small-note{color:#596579}.recurrent{color:#137333;font-weight:700}.transient{color:#b42318;font-weight:700}.tab-content{padding-top:18px}"))),
  div(class="hero", h2("Return or escape?"),
      p("An interactive laboratory for recurrent and transient Markov chains")),
  tabsetPanel(
    tabPanel("Finite chains",
      sidebarLayout(
        sidebarPanel(width=4,
          selectInput("example", "Choose a chain", choices=names(examples)),
          uiOutput("matrix_ui"),
          sliderInput("steps", "Length of sample path", 10, 300, 80),
          selectInput("start", "Starting state", 1:3),
          actionButton("simulate", "Run a new path", class="btn-primary"),
          hr(), p(class="small-note", "A state is recurrent when the probability of ever returning to it is 1. One short path cannot prove recurrence or transience.")),
        mainPanel(width=8,
          wellPanel(h4("What the structure says"), uiOutput("classification"), textOutput("example_note")),
          plotOutput("path_plot", height=250),
          fluidRow(
            column(6, plotOutput("return_plot", height=260)),
            column(6, tableOutput("powers_table"))
          )
        )
      )
    ),
    tabPanel("Biased random walk",
      sidebarLayout(
        sidebarPanel(width=4,
          sliderInput("p_right", "Probability p of moving right", 0.05, 0.95, 0.50, step=.05),
          sliderInput("rw_steps", "Steps per walk", 20, 1000, 200, step=20),
          sliderInput("walks", "Number of simulated walks", 100, 10000, 2000, step=100),
          actionButton("run_rw", "Simulate", class="btn-primary"),
          hr(), p(class="small-note", "For the walk on all integers, p = 1/2 is recurrent. Any drift, however small, makes the starting state transient.")),
        mainPanel(width=8,
          fluidRow(
            column(4, wellPanel(h4("Classification"), uiOutput("rw_class"))),
            column(4, wellPanel(h4("Ever-return probability"), div(class="value", textOutput("rw_theory")))),
            column(4, wellPanel(h4("Returned by horizon"), div(class="value", textOutput("rw_empirical"))))
          ),
          plotOutput("rw_paths", height=300),
          plotOutput("rw_curve", height=280)
        )
      )
    ),
    tabPanel("Key ideas",
      wellPanel(
        h3("The distinction in one line"),
        withMathJax(p("State \\(i\\) is recurrent if \\(f_{ii}=\\Pr_i(T_i^+<\\infty)=1\\), and transient if \\(f_{ii}<1\\).")),
        h4("What to notice"),
        tags$ul(
          tags$li("A recurrent state need not return quickly; recurrence concerns eventual return."),
          tags$li("A transient state can return many times on a particular path; it still has a positive probability of never returning."),
          tags$li("In a finite chain, every closed communicating class is recurrent."),
          tags$li("Periodicity controls when returns are possible, not whether eventual return has probability one."),
          tags$li("For the biased random walk, the return probability is 1 at p = 1/2 and 2 min(p, 1-p) otherwise."))
      )
    )
  )
)

server <- function(input, output, session) {
  P <- reactive(examples[[input$example]]$P)
  output$matrix_ui <- renderUI({
    m <- P()
    header_cells <- lapply(seq_len(ncol(m)), function(j) {
      tags$th(paste0("to ", j))
    })
    body_rows <- lapply(seq_len(nrow(m)), function(i) {
      probability_cells <- lapply(m[i, ], function(z) {
        tags$td(sprintf("%.2f", z))
      })
      tags$tr(tags$th(paste0("from ", i)), probability_cells)
    })

    tagList(
      h4("Transition matrix"),
      tags$table(class="table table-bordered table-condensed",
        tags$thead(tags$tr(tags$th(""), header_cells)),
        tags$tbody(body_rows)
      )
    )
  })
  output$classification <- renderUI({
    z <- classify_states(P())
    tags$p(lapply(seq_along(z), function(i) tags$span(class=tolower(z[i]), paste0("State ",i,": ",z[i], if(i<length(z)) "   |   " else ""))))
  })
  output$example_note <- renderText(examples[[input$example]]$note)

  path <- eventReactive(list(input$simulate, input$example, input$start, input$steps), {
    simulate_chain(P(), as.integer(input$start), input$steps)
  }, ignoreInit=FALSE)
  output$path_plot <- renderPlot({
    x <- path(); plot(0:(length(x)-1), x, type="s", lwd=2, col="#2563eb", yaxt="n", xlab="Time", ylab="State", main="One realization")
    axis(2, 1:nrow(P())); abline(h=1:nrow(P()), col="#d9deea", lty=3)
  })
  output$return_plot <- renderPlot({
    h <- max(20, input$steps); f <- first_return_curve(P(), as.integer(input$start), h)
    plot(seq_len(h), f, type="l", lwd=3, col="#7c3aed", ylim=c(0,1.02), xlab="Horizon n", ylab="P(return by n)", main=paste("Return to state",input$start))
    abline(h=1, lty=2, col="#64748b")
  })
  output$powers_table <- renderTable({
    m <- diag(nrow(P())); ks <- c(1,2,5,10,25,50); ans <- data.frame(n=ks)
    vals <- numeric(length(ks)); at <- 1
    for(k in seq_len(max(ks))){m <- m %*% P(); if(k %in% ks){vals[at] <- m[as.integer(input$start),as.integer(input$start)]; at <- at+1}}
    ans[[paste0("P^n[",input$start,",",input$start,"]")]] <- round(vals,4); ans
  }, caption="Probability of being back at the starting state at time n")

  rw <- eventReactive(list(input$run_rw, input$p_right, input$rw_steps, input$walks), {
    inc <- matrix(ifelse(runif(input$walks*input$rw_steps)<input$p_right,1,-1), input$walks)
    pos <- t(apply(inc,1,cumsum)); returned <- apply(pos==0,1,any)
    list(pos=pos, returned=returned)
  }, ignoreInit=FALSE)
  output$rw_class <- renderUI({
    rec <- abs(input$p_right-.5)<1e-10
    tags$p(class=if(rec) "recurrent" else "transient", if(rec) "Recurrent" else "Transient")
  })
  output$rw_theory <- renderText(sprintf("%.3f", if(abs(input$p_right-.5)<1e-10) 1 else 2*min(input$p_right,1-input$p_right)))
  output$rw_empirical <- renderText(sprintf("%.3f", mean(rw()$returned)))
  output$rw_paths <- renderPlot({
    z <- rw()$pos; k <- min(20,nrow(z)); matplot(0:ncol(z), t(cbind(0,z[seq_len(k),,drop=FALSE])), type="l", lty=1, col=adjustcolor("#2563eb",.25), xlab="Time", ylab="Position", main="First 20 simulated walks"); abline(h=0,lwd=2)
  })
  output$rw_curve <- renderPlot({
    z <- rw()$pos; cum <- apply(z==0,1,cumsum)>0; curve <- colMeans(cum)
    theory <- if(abs(input$p_right-.5)<1e-10) 1 else 2*min(input$p_right,1-input$p_right)
    plot(curve,type="l",lwd=3,col="#db2777",ylim=c(0,1.02),xlab="Simulation horizon",ylab="Fraction returned",main="Finite simulation approaches the eventual-return probability"); abline(h=theory,lty=2,lwd=2,col="#172554"); legend("bottomright",c("Simulated by horizon","Theoretical eventual return"),col=c("#db2777","#172554"),lty=c(1,2),lwd=2,bty="n")
  })
}

shinyApp(ui, server)

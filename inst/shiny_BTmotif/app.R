library(shiny)
library(shinyjs)
library(shinydashboard)
library(h2o)
library(ggplot2)
library(ggseqlogo)

NT <- c("A", "C", "G", "T")

shinyApp(
  #
  # Define UI -------
  #
  ui=dashboardPage(
    dashboardHeader(title="BTmotif: Extract sequence determinants", titleWidth=500),
    dashboardSidebar(disable=TRUE),
    dashboardBody(
      fluidRow(
        shinyjs::useShinyjs(),
        column(width=4,
          box(width=NULL, title="Training data", status="warning",
            textAreaInput("targets", label="Paste here the BI blunt rate analysis", height="200px",
                          value=paste("# Columns are <tab> or <comma> separated.",
                                      "# Columns are expected to be in this order: protospacer_sequence | blunt_rate",
                                      "# Columns are not named.",
                                      "# Lines starting with '#' are ignored", sep="\n")
            ),
            fluidRow(column(12, actionLink("example", "example"), align="right")),
            hr(),
            numericInput('ntrees', 'Number of trees', 1000, step=100),
            numericInput('nfolds', 'Number of folds for CV', 5),
            checkboxInput("scaleMotif", "Scale motif with importance", TRUE),
            hr(),
            fluidRow(
              column(6, actionButton("go", "Go!", icon=icon("play")), align="left"),
              column(6, shinyjs::disabled(downloadButton("download_model", "Download model")), align="right")
            )
          )
        ),
        column(width=8,
          box(width=NULL, title="Sequence determinants", status="warning",
            div(style='overflow-x: scroll', plotOutput("variableImportance")),
            div(style='overflow-x: scroll', plotOutput("motif"))
          ),
          box(width=NULL, title="Model Performance", status="warning",
            div(style='overflow-x: scroll', plotOutput("OEplot"))
          )#{{Impressum-placeholder}}
        )
      )
    )
  ),
  #
  # Define server logic -----
  #
  server=function(input, output, session) {
    #
    # init stuff and support functions
    #
    withProgress(message="Initializing H2O...", value=0, {
      # Require h2o 3.36.1.2 for the breaktag predictor
      # https://stackoverflow.com/questions/64801593/get-specific-h2o-version-in-r
      # install.packages("h2o", type="source", repos="http://h2o-release.s3.amazonaws.com/h2o/rel-zumbo/2/R")
      h2o.init()
    })
    onStop(function() h2o.shutdown(prompt=FALSE))

    # one-hot encode the sequences
    onehot <- function(x) {    # onehot a sequence of nt
      onehot1 <- function(x) {
        m <- setNames(integer(4), NT)
        m[x] <- 1
        m
      }
      unname(do.call(c, lapply(unlist(strsplit(x, "")), onehot1)))
    }
    onehotV <- Vectorize(onehot, SIMPLIFY=FALSE)

    #
    # reactive components -----
    #
    targets <- reactive({
      input$go
      isolate({
        tryCatch({
          read.table(stringsAsFactors=FALSE, strip.white=TRUE, fill=TRUE, text=gsub("[\t|,]+", "\t", input$targets))
        }, error=function(e) NULL)
      })
    })

    model <- reactive({
      req(!is.null(targets()))

      withProgress(message="Validating input targets...", value=0, {
        validate(
          need(ncol(targets()) == 2, "Columns are expected to be in this order: protospacer_sequence | blunt_rate"),
          need(is(targets()[[1]], "character"), "First column must contain a valid protospacer sequence"),
          need(all(nchar(targets()[[1]]) == nchar(targets()[[1]][1])), "Protospacers must be all of the same length"),
          need(is(targets()[[2]], "numeric"), "Second column must contain valid numbers")
        )
      })

      withProgress(message="One-hot encoding protospacer sequences...", value=0, {
        tryCatch({
          onehot <- as.list(as.data.frame(onehotV(toupper(targets()[[1]]))))
        }, error=function(e) NULL)
      })

      withProgress(message="Training model, this may take a while. Check the R console for progress...", value=0, {
        tryCatch({
          df <- as.data.frame(t(as.data.frame(onehot)))
          colnames(df) <- paste("p", paste(rep(1:nchar(targets()[[1]][1]), each=4), NT, sep="_"), sep="_")
          df$blunt_rate <- targets()[[2]]
          response   <- "blunt_rate"
          predictors <- setdiff(colnames(df), response)

          # Train blunt/staggered regression model
          df.h2o <- as.h2o(df)
          h2o.xgboost(x=predictors,
                      y=response,
                      training_frame=df.h2o,
                      ntrees=isolate(input$ntrees),
                      max_depth=0,
                      nfolds=isolate(input$nfolds),
                      keep_cross_validation_models=TRUE,
                      keep_cross_validation_predictions=TRUE,
                      keep_cross_validation_fold_assignment=TRUE,
                      booster="dart",
                      normalize_type="tree",
                      seed=1234)
        }, error=function(e) NULL)
      })
    })

    observeEvent(input$example, {
      x <- read.delim(system.file("extdata/BTmotif.txt", package="breakinspectoR"), head=FALSE, sep="\t")
      x <- paste(apply(x, 1, paste, collapse=","), collapse="\n")
      updateTextInput(session, "targets", value=x)
    })

    observe({
      req(model(), cancelOutput=TRUE)
      #updateActionButton(session, "download_model", disabled=FALSE)
      shinyjs::enable("download_model")
    })

    #
    # UI ------
    #
    output$variableImportance <- renderPlot({
      req(model(), cancelOutput=TRUE)

      varimp <- as.data.frame(h2o.varimp(model()))
      varimp$pos  <- as.numeric(sub("^p_(\\d+)_([ACTG])", "\\1", as.character(varimp$variable)))
      varimp$to   <- sub("^p_(\\d+)_([ACTG])", "\\2", as.character(varimp$variable))

      ggplot(varimp, aes(x=pos, y=relative_importance, color=to)) +
        geom_line() +
        labs(title="XGBoost blunt regressor", subtitle="Position/Nucleotide importance", y="relative importance", x="protospacer position") +
        scale_color_manual("protospacer is", values=palette()) +
        scale_x_continuous(breaks=1:nchar(targets()[[1]][1])) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              panel.grid.minor=element_blank(),
              legend.position="bottom")
    })

    output$motif <- renderPlot({
      req(model(), cancelOutput=TRUE)

      varimp <- as.data.frame(h2o.varimp(model()))
      varimp$pos  <- as.numeric(sub("^p_(\\d+)_([ACTG])", "\\1", as.character(varimp$variable)))
      varimp$to   <- sub("^p_(\\d+)_([ACTG])", "\\2", as.character(varimp$variable))

      y <- as.data.frame(t(as.data.frame(strsplit(targets()[[1]], ""))))
      coefs <- apply(y, 2, function(y) tryCatch({
        coef(lm(scale(targets()[[2]], center=TRUE, scale=FALSE) ~ 0+y))
      }, error=function(e) rep(0, 4)))
      rownames(coefs) <- NT

      if(input$scaleMotif) {
        i <- sapply(1:nchar(targets()[[1]][1]), function(x) with(subset(varimp, pos == x), tapply(scaled_importance, to, sum)))
        coefs <- coefs * i
      }

      ggplot() +
        ggseqlogo::geom_logo(coefs, method="custom", seq_type="dna") +
        labs(title="Local explanation of the observed blunt rate",
             y="scaled regression coefficient",
             x="protospacer position") +
        ggseqlogo::theme_logo()
    })

    output$OEplot <- renderPlot({
      req(model(), cancelOutput=TRUE)
      
      cvpreds <- h2o.getFrame(model()@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
      x <- data.frame(expected=as.data.frame(cvpreds)$predict,
                      observed=targets()[[2]])

      opar <- par(mfrow=c(1, 2))
      smoothScatter(x$expected ~ x$observed, nbin=512, nrpoints=0,
                    main=paste0("Prediction performance (on cross-valitaded data)\nR=", round(cor(x$observed, x$expected), 2)),
                    xlab="observed blunt rate", ylab="predicted blunt rate",
                    colramp=colorRampPalette(c("black", "blue", "cyan", "yellow", "red")))
      abline(0, 1, col="blue", lty=1, lwd=2)
      abline(lm(x$expected ~ x$observed), col="blue", lty=2, lwd=2)

      hist(x$observed - x$expected, breaks=100, main="Distribution of residuals")
      par(opar)
    })

    output$download_model <- downloadHandler("model.h2o", content=function(file) {
      req(model(), cancelOutput=TRUE)
      h2o.saveModel(model(), force=TRUE, export_cross_validation_predictions=TRUE, filename=file)
    })
  }
)

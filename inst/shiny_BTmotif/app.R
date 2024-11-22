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
              footer="Right click and 'Save Image as' to get a vectorized version for your publication",
            div(style='overflow-x: scroll', imageOutput("variableImportance")),
            div(style='overflow-x: scroll', imageOutput("motif"))
          ),
          box(width=NULL, title="Model Performance", status="warning",
              footer="Right click and 'Save Image as' to get a vectorized version for your publication",
            div(style='overflow-x: scroll', imageOutput("OEplot"))
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
    # simpler version which doesn't take into account protospacer-sgRNA mismatches as in the NatBiotech paper
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
    variableImportanceSVG <- reactive({
      req(model(), cancelOutput=TRUE)

      varimp <- as.data.frame(h2o.varimp(model()))
      varimp$pos  <- as.numeric(sub("^p_(\\d+)_([ACTG])", "\\1", as.character(varimp$variable)))
      varimp$to   <- sub("^p_(\\d+)_([ACTG])", "\\2", as.character(varimp$variable))

      p <- ggplot(varimp, aes(x=pos, y=relative_importance, color=to)) +
        geom_line() +
        labs(title="XGBoost blunt regressor", subtitle="Position/Nucleotide importance", y="relative importance", x="protospacer position") +
        scale_color_manual("protospacer is", values=palette()) +
        scale_x_continuous(breaks=1:nchar(targets()[[1]][1])) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              panel.grid.minor=element_blank(),
              legend.position="bottom")

      file <- htmltools::capturePlot(p, tempfile(fileext=".svg"), grDevices::svg, width=8, heigh=4)
    })
    output$variableImportance <- renderImage(list(src=variableImportanceSVG()), deleteFile=TRUE)

    motifSVG <- reactive({
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

      p <- ggplot() +
        ggseqlogo::geom_logo(coefs, method="custom", seq_type="dna") +
        labs(title="Local explanation of the observed blunt rate",
             y="scaled regression coefficient",
             x="protospacer position") +
        ggseqlogo::theme_logo()

      file <- htmltools::capturePlot(p, tempfile(fileext=".svg"), grDevices::svg, width=8, heigh=4)
    })
    output$motif <- renderImage(list(src=motifSVG()), deleteFile=TRUE)

    OEplotSVG <- reactive({
      req(model(), cancelOutput=TRUE)
      
      cvpreds <- h2o.getFrame(model()@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
      x <- data.frame(expected=as.data.frame(cvpreds)$predict,
                      observed=targets()[[2]])

      p1 <- ggplot(x, aes(x=observed, y=expected)) +
        geom_density_2d_filled(show.legend=FALSE) +
        geom_smooth(method="lm", color="blue", lty=2) +
        geom_abline(slope=1, intercept=0, color="blue", lty=1) +
        labs(title=paste0("Prediction performance (on cross-valitaded data)\nR=", round(cor(x$observed, x$expected), 2)),
            x="observed blunt rate", y="predicted blunt rate") +
        theme_bw() +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

      p2 <- ggplot(x, aes(x=expected)) +
        geom_histogram() +
        labs(title="Distribution of residuals") +
        theme_bw()

      file <- htmltools::capturePlot(cowplot::plot_grid(p1, p2, cols=2),
                                     tempfile(fileext=".svg"), grDevices::svg, width=8, heigh=4)
    })
    output$OEplot <- renderImage(list(src=OEplotSVG()), deleteFile=TRUE)

    output$download_model <- downloadHandler("model.h2o", content=function(file) {
      req(model(), cancelOutput=TRUE)
      h2o.saveModel(model(), force=TRUE, export_cross_validation_predictions=TRUE, filename=file)
    })
  }
)

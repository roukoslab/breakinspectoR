library(shiny)
library(shinydashboard)
library(DT)
library(h2o)

options(shiny.maxRequestSize = 100 * 1024^2)   # allow uploading models up to 100MB
NT <- c("A", "C", "G", "T")

shinyApp(
  #
  # Define UI -------
  #
  ui=dashboardPage(
    dashboardHeader(title="Blunt rate prediction (Chakrabarti's)", titleWidth=300),
    dashboardSidebar(disable=TRUE),
    dashboardBody(
      fluidRow(
        column(width=4,
          box(width=NULL, title="gRNAs", status="warning",
            textAreaInput("gRNAs", label="Paste here your gRNAs", height="200px",
                          value=paste("# Columns are <tab> or <comma> separated.",
                                      "# First column is mandatory and must contain gRNA target sequences.",
                                      "# If using the default model used in the Nat Biotech paper, sequences must be at least 10nt.",
                                      "# If using a custom model, sequences must be the same length as in the training set.",
                                      "# Build custom models using the companion breakinspectoR::shiny_BTmotif() app.",
                                      "# Other columns are optional.",
                                      "# Columns are not named.",
                                      "# Lines starting with '#' are ignored", sep="\n")
            ),
            fluidRow(column(12, actionLink("example", "example"), align="right")),
            hr(),
            fluidRow(
              column(6, actionButton("predict", "Predict!", icon=icon("play")), align="left"),
              column(6, fileInput("uploadModel", "Upload a custom model", accept=".h2o"), align="right")
            )
          )
        ),
        column(width=8,
          box(width=NULL, title="Predicted table", status="warning",
            div(style='overflow-x: scroll', DTOutput("predTable")),
            downloadButton("downloadPredTable", "Download predicted table")
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

    # load models
    model <- reactiveValues()
    withProgress(message="Loading default spCas9 XBoost regressor model trained with Hiplex1 data...", value=0, {
      # defaults to the XGBoost regressor trained with HiPlex1 data that we produced with the Nat Biotechnology paper
      model$h2o     <- h2o.loadModel(system.file("extdata/bluntPred.hiplex1.h2o", package="breakinspectoR"))
      model$default <- TRUE
      model$default_columns <- isolate(model$h2o@model$names)
    })

    observe({
      req(input$uploadModel, cancelOutput=TRUE)
      withProgress(message="Loading custom model...", value=0, {
        model$h2o <- h2o.loadModel(input$uploadModel$datapath)
        model$default <- length(isolate(model$h2o@model$names)) == length(isolate(model$default_columns)) &&
                         all(isolate(model$h2o@model$names) == isolate(model$default_columns))
      })
    })

    # one-hot encode the sequences as in the Nat Biotech paper
    onehotV <- Vectorize(function(x, y) {    # onehot a sequence of nt
      onehot1 <- function(x, y) { # onehot one single nt
        if(y == "N") {
          y <- x  # nt "N" in `y` always match `x` (eg, NGG protospacer, N will always match)
        }
        if(x == "N") {
          x <- y  # same for `x`
        }
        m <- matrix(0, nrow=4, ncol=4, dimnames=list(observed=NT, expected=NT))
        m[x, y] <- 1
        as.vector(m)
      }

      unname(do.call(c, Map(unlist(strsplit(x, "")), unlist(strsplit(y, "")), f=onehot1)))
    })

    # one-hot encode the sequences as in shiny_BTmotif, for custom models
    # the main difference is that here we use the whole protospacer and we don't include mismatches between protospacer-sgRNA
    onehot2V <- Vectorize(function(x) {    # onehot a sequence of nt
      onehot1 <- function(x) {
        m <- setNames(integer(4), NT)
        m[x] <- 1
        m
      }
      unname(do.call(c, lapply(unlist(strsplit(x, "")), onehot1)))
    }, SIMPLIFY=FALSE)

    #
    # reactive components -----
    #
    gRNAs <- reactive({
      input$predict
      isolate({
        tryCatch({
          read.table(stringsAsFactors=FALSE, strip.white=TRUE, fill=TRUE, text=gsub("[\t|,]+", "\t", input$gRNAs))
        }, error=function(e) NULL)
      })
    })

    predgRNAs <- reactive({
      req(!is.null(gRNAs()))
      withProgress(message="Predicting blunt rates, this may take a while...", value=0, {
        tryCatch({
          # encode
          # support either the model trained with Hiplex data for the Nat Biotech paper included by default in this app (last 10 nt)
          # and support custom models trained with the shiny_BTmotif() app included in this paper
          if(isolate(model$default)) {
            validate(need(all(nchar(gRNAs()[[1]]) >= 10), "Length of predicted sequences must be at least 10nt for the default Hiplex 1 model"))
            guide <- toupper(substr(gRNAs()[[1]], nchar(gRNAs()[[1]]) - 10 + 1, nchar(gRNAs()[[1]])))
            onehot <- as.list(as.data.frame(onehotV(guide, guide)))
            df <- as.data.frame(t(as.data.frame(onehot)))
            rownames(df) <- NULL
            colnames(df) <- paste("p", paste(rep(1:10, each=16), paste(NT, rep(NT, each=4), sep="_instead_of_"), sep="_"), sep="_")
          } else {
            n_training_vars <- (length(isolate(model$h2o@model$names)) - 1) / 4  # length of sequences in the training set
            validate(need(all(nchar(gRNAs()[[1]]) >= n_training_vars), paste("Length of predicted sequences must be at least", n_training_vars, "nt for this custom model")))
            guide <- toupper(gRNAs()[[1]])                   # use the whole sequence
            onehot <- as.list(as.data.frame(onehot2V(guide)))
            df <- as.data.frame(t(as.data.frame(onehot)))
            rownames(df) <- NULL
            colnames(df) <- paste("p", paste(rep(1:nchar(gRNAs()[[1]][1]), each=4), NT, sep="_"), sep="_")
          }

          # predict the outcome for each gRNA
          df.h2o <- as.h2o(df)
          blunt_rate <- tryCatch({
            h2o.predict(isolate(model$h2o), df.h2o)
          }, error=function(e) {
            NULL
          })
          validate(need(!is.null(blunt_rate), "Prediction failed. Ensure the length of the predicted sequences matches the length of the training sequences"))
          x <- cbind(round(as.data.frame(blunt_rate)$predict, digits=2), gRNAs())
          colnames(x)[1] <- "predicted blunt rate"
          x
        })
      })
    })

    observeEvent(input$example, {
      x <- read.delim(system.file("extdata/bluntPred.txt", package="breakinspectoR"), head=FALSE, sep="\t")
      x <- paste(apply(x, 1, paste, collapse=","), collapse="\n")
      updateTextInput(session, "gRNAs", value=x)
    })

    #
    # UI ------
    #
    output$predTable <- renderDT({
      req(predgRNAs(), cancelOutput=TRUE)
      datatable(predgRNAs(), rownames=FALSE, selection="none", options=list(pageLength=5))
    })

    output$downloadPredTable <- downloadHandler(
      filename="predTable.csv",
      content=function(f) {
        if(!is.null(predgRNAs())) {
          write.csv(predgRNAs(), file=f, row.names=FALSE)
        }
      }
    )
  }
)

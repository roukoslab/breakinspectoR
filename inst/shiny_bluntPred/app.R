library(shiny)
library(shinydashboard)
library(DT)
library(h2o)

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
                                             "# First column is mandatory and must contain gRNA target sequences of at least 10nt.",
                                             "# Other columns are optional.",
                                             "# Columns are not named.",
                                             "# Lines starting with '#' are ignored", sep="\n")
                   ),
                   hr(),
                   fluidRow(column(6, actionLink("predict", "Predict!")),
                            column(6, actionLink("example", "example"), align="right"),
                   ),
               )
        ),
        column(width=8,
               box(width=NULL, title="Predicted table", status="warning",
                   div(style='overflow-x: scroll', DTOutput("predTable")),
                   downloadLink("downloadPredTable", "Download predicted table")
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

    withProgress(message="Loading XBoost regressor model...", value=0, {
      # defaults to the XGBoost regressor trained with HiPlex1 data that we produced with the Nat Biotechnology paper
      blunt_regressor.xgb <- h2o.loadModel(system.file("extdata/bluntPred.hiplex1.h2o", package="breakinspectoR"))
    })

    # one-hot encode the sequences
    onehot <- function(x, y) {    # onehot a sequence of nt
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
    }
    onehotV <- Vectorize(onehot)

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
          guide10 <- toupper(substr(gRNAs()[[1]], 11, 20))
          onehot10 <- as.list(as.data.frame(onehotV(guide10, guide10)))

          # predict the outcome for each gRNA
          df <- as.data.frame(t(as.data.frame(onehot10)))
          rownames(df) <- NULL
          colnames(df) <- paste("p", paste(rep(1:10, each=16), paste(NT, rep(NT, each=4), sep="_instead_of_"), sep="_"), sep="_")
          df.h2o <- as.h2o(df)
          blunt_rate <- h2o.predict(blunt_regressor.xgb, df.h2o)
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

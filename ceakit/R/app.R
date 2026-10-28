#' Shiny app for any ceakit project
#'
#' Upload or edit every input table, change settings, run each analysis and
#' explore results (summary table, key metrics, cost-effectiveness, value of
#' information, project figures) and download them.
#' @param dir Project folder.
#' @return A `shiny.appobj`. In the project's `app.R` use `ceakit::cea_app(".")`.
#' @export
cea_app <- function(dir = ".") {
  for (p in c("shiny", "bslib", "DT")) if (!requireNamespace(p, quietly = TRUE))
    stop("install the ", p, " package to use the app", call. = FALSE)
  proj <- load_project(dir)
  `%||%` <- function(a, b) if (is.null(a)) b else a
  std <- c("wtp", "wtp_grid", "n_sim", "seed", "owsa_pct", "evppi", "evppi_groups", "population")
  extra <- Filter(function(x) is.numeric(x) && length(x) == 1, proj$settings[setdiff(names(proj$settings), std)])
  tables0 <- c(list(parameters = proj$params_raw), proj$data[vapply(proj$data, is.data.frame, logical(1))],
               if (!is.null(proj$counts)) list(counts = proj$counts))
  labels_of <- function(P) stats::setNames(dplyr::coalesce(P$label, P$parameter), P$parameter)

  results_ui <- function(id, a) {
    ns <- shiny::NS(id)
    shiny::tagList(
      shiny::uiOutput(ns("status")),
      bslib::navset_card_tab(
        bslib::nav_panel("Summary table", DT::DTOutput(ns("wide"))),
        bslib::nav_panel("Key results", shiny::plotOutput(ns("metrics"), height = "640px")),
        bslib::nav_panel("Cost-effectiveness", DT::DTOutput(ns("inc")),
          bslib::layout_columns(shiny::plotOutput(ns("plane"), height = "440px"),
                                shiny::plotOutput(ns("ceac"), height = "440px")),
          shiny::plotOutput(ns("tornado"), height = "520px")),
        bslib::nav_panel("Value of information",
          bslib::layout_columns(shiny::plotOutput(ns("evpi"), height = "440px"),
                                shiny::plotOutput(ns("evppi"), height = "440px"))),
        if (length(a$plots)) bslib::nav_panel("Model figures",
          shiny::selectInput(ns("fig"), NULL, names(a$plots)), shiny::plotOutput(ns("fig_out"), height = "680px")),
        bslib::nav_panel("Download", shiny::p("Excel workbook, or a zip of all tables and figures."),
          shiny::downloadButton(ns("xlsx"), "Results (Excel)"), " ",
          shiny::downloadButton(ns("zip"), "Tables + figures (zip)"))))
  }

  results_server <- function(id, res, data, labels) shiny::moduleServer(id, function(input, output, session) {
    output$status <- shiny::renderUI({
      r <- res()
      if (is.null(r)) shiny::div(class = "alert alert-info", "Not run yet: check inputs, then press Run in the sidebar.")
      else shiny::div(class = "alert alert-light", sprintf("Last run %s | PSA runs: %s | WTP %s",
        format(r$time, "%H:%M:%S"), if (is.null(r$psa)) 0 else dim(r$psa)[3], fmt(r$settings$wtp)))
    })
    output$wide <- DT::renderDT({ shiny::req(res())
      DT::datatable(res()$wide, rownames = FALSE, options = list(pageLength = 100, dom = "ft", scrollX = TRUE)) })
    output$metrics <- shiny::renderPlot(res = 100, { shiny::req(res()); plot_metrics(res()) })
    output$inc <- DT::renderDT({ shiny::req(res())
      d <- res()$incremental
      d$cost <- scales::dollar(d$cost, 0.01); d$effect <- round(d$effect, 5)
      d$icer <- ifelse(is.na(d$icer), "-", scales::dollar(d$icer, 1))
      DT::datatable(d, rownames = FALSE, options = list(dom = "t"), caption = "Incremental analysis") })
    output$plane <- shiny::renderPlot(res = 100, { shiny::req(res()); plot_ce_plane(res()) })
    output$ceac <- shiny::renderPlot(res = 100, { shiny::req(res()); plot_ceac(res()) })
    output$tornado <- shiny::renderPlot(res = 100, { shiny::req(res()); plot_tornado(res(), labels()) })
    output$evpi <- shiny::renderPlot(res = 100, { shiny::req(res()); plot_evpi(res()) })
    output$evppi <- shiny::renderPlot(res = 100, { shiny::req(res()); plot_evppi(res(), labels()) })
    output$fig_out <- shiny::renderPlot(res = 100, {
      shiny::req(res(), input$fig)
      p <- res()$analysis$plots[[input$fig]](res(), res()$p, data())
      if (inherits(p, "ggplot")) print(p)
    })
    output$xlsx <- shiny::downloadHandler(
      filename = function() paste0("results_", res()$analysis$name, ".xlsx"),
      content = function(file) { r <- res()
        openxlsx::write.xlsx(Filter(Negate(is.null), list(Results = r$wide, Incremental = r$incremental,
          Long = r$long, EVPPI = r$evppi, OWSA = r$owsa$table, Parameters = r$params)), file) })
    output$zip <- shiny::downloadHandler(
      filename = function() paste0("results_", res()$analysis$name, ".zip"),
      content = function(file) {
        d <- file.path(tempdir(), paste0("cea_", as.integer(Sys.time())))
        shiny::withProgress(message = "Writing files", value = 0.5, write_cea(res(), d, data(), labels()))
        old <- setwd(d); on.exit(setwd(old)); utils::zip(file, list.files(".", recursive = TRUE), flags = "-q")
      })
  })

  table_panel <- function(id) bslib::nav_panel(id,
    bslib::layout_columns(fill = FALSE, col_widths = c(5, 7),
      shiny::fileInput(paste0("up_", id), NULL, accept = c(".csv", ".xlsx", ".xls"), buttonLabel = "Upload CSV/Excel"),
      shiny::div(shiny::downloadButton(paste0("dl_", id), "Download table"), " ",
                 shiny::actionButton(paste0("reset_", id), "Reset", class = "btn-outline-secondary"))),
    DT::DTOutput(paste0("tbl_", id)))

  panels <- c(
    list(bslib::nav_panel("Inputs", shiny::uiOutput("validation"),
      do.call(bslib::navset_card_tab, unname(lapply(names(tables0), table_panel))))),
    unname(lapply(names(proj$analyses), function(n)
      bslib::nav_panel(proj$analyses[[n]]$title, results_ui(n, proj$analyses[[n]])))),
    list(bslib::nav_panel("About", bslib::card(bslib::card_body(
      if (file.exists(file.path(proj$dir, "README.md")) && requireNamespace("markdown", quietly = TRUE))
        shiny::includeMarkdown(file.path(proj$dir, "README.md")) else shiny::p("ceakit project"))))))
  sidebar <- bslib::sidebar(width = 300,
    shiny::numericInput("wtp", "Willingness to pay", proj$settings$wtp, min = 0, step = 1000),
    shiny::numericInput("n_sim", "PSA simulations (0 = none)", min(proj$settings$n_sim, 500), min = 0, step = 100),
    shiny::numericInput("seed", "Seed", proj$settings$seed),
    shiny::numericInput("owsa_pct", "One-way SA range (+/-)", proj$settings$owsa_pct, min = 0.01, max = 1, step = 0.05),
    shiny::checkboxInput("evppi", "Compute EVPPI", isTRUE(proj$settings$evppi)),
    lapply(names(extra), function(n) shiny::numericInput(paste0("x_", n), n, extra[[n]])),
    shiny::uiOutput("badge"),
    shiny::selectInput("which", "Analysis", stats::setNames(names(proj$analyses),
                       vapply(proj$analyses, `[[`, "", "title"))),
    shiny::actionButton("run", "Run", class = "btn-primary w-100"))
  ui <- do.call(bslib::page_navbar, c(list(
    id = "main", title = proj$env$project_title %||% "Decision model",
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"), fillable = FALSE,
    header = shiny::tags$style(shiny::HTML(".alert{margin-bottom:.75rem}")), sidebar = sidebar), panels))

  server <- function(input, output, session) {
    tabs <- do.call(shiny::reactiveValues, tables0)
    for (id in names(tables0)) local({
      id <- id
      tidy <- if (id == "parameters") tidy_parameters else identity
      output[[paste0("tbl_", id)]] <- DT::renderDT(DT::datatable(tabs[[id]], editable = "cell", rownames = FALSE,
        options = list(pageLength = 25, scrollX = TRUE)), server = TRUE)
      shiny::observeEvent(input[[paste0("tbl_", id, "_cell_edit")]], {
        tabs[[id]] <- tidy(DT::editData(tabs[[id]], input[[paste0("tbl_", id, "_cell_edit")]], rownames = FALSE)) })
      shiny::observeEvent(input[[paste0("up_", id)]], {
        f <- input[[paste0("up_", id)]]
        d <- tryCatch(tidy(read_table(f$datapath, f$name)), error = function(e) {
          shiny::showNotification(conditionMessage(e), type = "error"); NULL })
        if (!is.null(d)) tabs[[id]] <- d })
      shiny::observeEvent(input[[paste0("reset_", id)]], tabs[[id]] <- tables0[[id]])
      output[[paste0("dl_", id)]] <- shiny::downloadHandler(filename = function() paste0(id, ".csv"),
        content = function(file) readr::write_csv(tabs[[id]], file, na = ""))
    })
    settings <- shiny::reactive({
      s <- proj$settings
      s$wtp <- input$wtp; s$wtp_grid <- seq(0, 3 * input$wtp, length.out = 61)
      s$n_sim <- input$n_sim; s$seed <- input$seed; s$owsa_pct <- input$owsa_pct; s$evppi <- input$evppi
      for (n in names(extra)) s[[n]] <- input[[paste0("x_", n)]]
      s
    })
    params <- shiny::reactive(tryCatch(apply_counts(tabs$parameters, if ("counts" %in% names(tables0)) tabs$counts),
                                       error = function(e) structure(tabs$parameters, err = conditionMessage(e))))
    data <- shiny::reactive({ d <- proj$data; for (n in intersect(names(d), names(tables0))) d[[n]] <- tabs[[n]]; d })
    check <- shiny::reactive({
      v <- check_project(proj, params(), data(), settings())
      if (!is.null(attr(params(), "err"))) v$errors <- c(attr(params(), "err"), v$errors)
      v })
    output$validation <- shiny::renderUI({ v <- check(); shiny::tagList(
      if (length(v$errors)) shiny::div(class = "alert alert-danger", shiny::strong("Fix before running:"),
                                       shiny::tags$ul(lapply(v$errors, shiny::tags$li))),
      if (length(v$warnings)) shiny::div(class = "alert alert-warning", shiny::tags$details(
        shiny::tags$summary(sub(":.*", "", v$warnings[1]), " (click to list)"), paste(v$warnings, collapse = " "))),
      if (!length(v$errors) && !length(v$warnings)) shiny::div(class = "alert alert-success", "All inputs valid.")) })
    output$badge <- shiny::renderUI({ v <- check()
      if (length(v$errors)) shiny::span(class = "badge bg-danger mb-2", sprintf("%d input error(s)", length(v$errors)))
      else if (length(v$warnings)) shiny::span(class = "badge bg-warning text-dark mb-2", "Valid (DUMMY values remain)")
      else shiny::span(class = "badge bg-success mb-2", "Inputs valid") })
    results <- do.call(shiny::reactiveValues, stats::setNames(vector("list", length(proj$analyses)), names(proj$analyses)))
    shiny::observeEvent(input$run, {
      if (length(check()$errors)) { shiny::showNotification("Fix input errors first (Inputs tab).", type = "error"); return() }
      a <- proj$analyses[[input$which]]
      r <- shiny::withProgress(message = paste("Running", a$title), value = 0.02, tryCatch(
        run_cea(a, params(), data(), settings(), progress = function(x) shiny::setProgress(x)),
        error = function(e) { shiny::showNotification(conditionMessage(e), type = "error", duration = NULL); NULL }))
      if (!is.null(r)) { results[[input$which]] <- r; bslib::nav_select("main", a$title) }
    })
    for (n in names(proj$analyses)) local({ n <- n
      results_server(n, shiny::reactive(results[[n]]), data, shiny::reactive(labels_of(params()))) })
  }
  shiny::shinyApp(ui, server)
}

#' Launch the app for a project
#' @param dir Project folder.
#' @param ... Passed to [shiny::runApp()].
#' @export
run_app <- function(dir = ".", ...) shiny::runApp(cea_app(dir), ...)

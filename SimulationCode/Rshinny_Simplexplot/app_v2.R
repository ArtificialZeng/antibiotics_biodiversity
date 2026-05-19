# Three-Strain / Two-Antibiotic Frequency Dynamics (Replicator Map)
# ---------------------------------------------------------------
# Discrete-time map on the simplex (N=3, M=2) with a ternary trajectory plot
# and a flow/stream-density plot for the whole state space.

suppressPackageStartupMessages({
  library(shiny)
  library(ggtern)
  library(dplyr)
  library(raster)
  library(rasterVis)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
})

# -------------------------
# Phenotype matrices (edit as needed)
# -------------------------
S <- matrix(1, nrow = 3, ncol = 2, byrow = TRUE,
            dimnames = list(paste0("strain", 1:3), paste0("ab", 1:2)))

Pmat <- matrix(c(1,0,
                 0,0,
                 0,1), nrow = 3, byrow = TRUE,
               dimnames = list(paste0("strain", 1:3), paste0("ab", 1:2)))

Dmat <- matrix(c(0,0,
                 1,0,
                 0,0), nrow = 3, byrow = TRUE,
               dimnames = list(paste0("strain", 1:3), paste0("ab", 1:2)))

Rmat <- matrix(c(0,0,
                 0,1,
                 1,0), nrow = 3, byrow = TRUE,
               dimnames = list(paste0("strain", 1:3), paste0("ab", 1:2)))

# -------------
# UI
# -------------
ui <- fluidPage(
  titlePanel("Three-Strain / Two-Antibiotic Frequency Dynamics"),
  sidebarLayout(
    sidebarPanel(
      h4("Growth & interaction strength"),
      fluidRow(
        column(4, sliderInput("g1", HTML("g<sub>1</sub>"), min = 0, max = 2, value = 1, step = 0.01)),
        column(4, sliderInput("g2", HTML("g<sub>2</sub>"), min = 0, max = 2, value = 1, step = 0.01)),
        column(4, sliderInput("g3", HTML("g<sub>3</sub>"), min = 0, max = 2, value = 1, step = 0.01))
      ),
      fluidRow(
        column(6, sliderInput("KD", HTML("K<sub>D</sub> (degradation radius)"), min = 0, max = 30, value = 10, step = 0.1)),
        column(6, sliderInput("KP", HTML("K<sub>P</sub> (killing radius)"),    min = 0, max = 30, value = 15, step = 0.1))
      ),
      tags$hr(),
      h4("Metabolic costs (per antibiotic)"),
      fluidRow(
        column(4, sliderInput("cres1", HTML("c<sub>res</sub>[1]"), min = 0, max = 0.5, value = 0.05, step = 0.005)),
        column(4, sliderInput("cdeg1", HTML("c<sub>deg</sub>[1]"), min = 0, max = 0.5, value = 0.105, step = 0.005)),
        column(4, sliderInput("cprod1", HTML("c<sub>prod</sub>[1]"), min = 0, max = 0.5, value = 0.15, step = 0.005))
      ),
      fluidRow(
        column(4, sliderInput("cres2", HTML("c<sub>res</sub>[2]"), min = 0, max = 0.5, value = 0.05, step = 0.005)),
        column(4, sliderInput("cdeg2", HTML("c<sub>deg</sub>[2]"), min = 0, max = 0.5, value = 0.105, step = 0.005)),
        column(4, sliderInput("cprod2", HTML("c<sub>prod</sub>[2]"), min = 0, max = 0.5, value = 0.15, step = 0.005))
      ),
      tags$hr(),
      h4("Initial frequencies"),
      sliderInput("x1", HTML("X<sub>1</sub>"), min = 0, max = 1, value = 1/3, step = 0.01),
      sliderInput("x2", HTML("X<sub>2</sub>"), min = 0, max = 1, value = 1/3, step = 0.01),
      div("X3 is computed as 1 - X1 - X2 and must be ≥ 0."),
      tags$hr(),
      sliderInput("T", "Iterations (discrete steps)", min = 1, max = 400, value = 150, step = 1),
      actionButton("reset", "Reset defaults")
    ),
    mainPanel(
      navbarPage(title = "Outputs",
                 tabPanel(title = "Trajectory",
                          fluidPage(
                            plotOutput("ternplot", height = "650px"),
                            br(), htmlOutput("x3txt"), br(),
                            h4("Phenotype matrices (fixed in this demo; edit in code)"),
                            tableOutput("matrices")
                          )
                 ),
                 tabPanel(title = "Flow / Streamdensity",
                          fluidPage(
                            fluidRow(
                              column(8,
                                     wellPanel(
                                       sliderInput("gridN", "Grid resolution (per side)", min = 40, max = 600, value = 200, step = 10),
                                       checkboxInput("normvel", "Normalize arrow lengths", TRUE),
                                       sliderInput("arrows", "Number of arrows", min = 200, max = 5000, value = 1200, step = 50)
                                     ),
                                     plotOutput("flowplot", height = "650px")
                              ),
                              column(4,
                                     helpText("Vector field of the discrete map: arrows show ΔX = X' − X on the simplex."),
                                     helpText("Colors encode relative speed (|ΔX|).")
                              )
                            )
                          )
                 )
      )
    )
  )
)

# -------------
# Helpers
# -------------
clip01 <- function(x) pmin(pmax(x, 0), 1)

# One-step update of the discrete replicator map
step_map <- function(X, params) {
  g    <- params$g
  KD   <- params$KD
  KP   <- params$KP
  cres <- params$cres
  cdeg <- params$cdeg
  cprod<- params$cprod
  
  r <- numeric(3)
  for (i in 1:3) {
    r[i] <- g[i] - sum(Pmat[i,] * cprod) - sum(Dmat[i,] * cdeg) - sum(Rmat[i,] * cres)
  }
  
  Aj <- as.numeric(t(X) %*% Dmat)
  Bj <- as.numeric(t(X) %*% Pmat)
  
  pij <- matrix(0, nrow = 3, ncol = 2)
  for (i in 1:3) for (j in 1:2)
    pij[i, j] <- S[i, j] * exp(-KD * Aj[j]) * (1 - exp(-KP * Bj[j]))
  
  fi <- r * (1 - pij[,1]) * (1 - pij[,2])
  fi[fi < 0] <- 0
  
  denom <- sum(X * fi)
  if (!is.finite(denom) || denom <= 1e-12) return(X)
  as.numeric(X * fi / denom)
}

# Map-based vector field on the triangle using raster + rasterVis
compute_flow_field <- function(params, Ngrid = 250, normalize = TRUE) {
  u <- v <- raster(xmn=0, xmx=1, ymn=0, ymx=1, ncol=Ngrid, nrow=Ngrid)
  x <- init(u, v='x'); y <- init(u, v='y')
  
  x1 <- 0;   y1 <- 0
  x2 <- 1;   y2 <- 0
  x3 <- 0.5; y3 <- sqrt(3)/2
  
  denom <- (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
  lam1 <- ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / denom   # X1
  lam3 <- ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / denom   # X3
  lam2 <- 1 - lam1 - lam3                                      # X2
  
  maskTri <- (lam1 >= 0) & (lam2 >= 0) & (lam3 >= 0)
  
  X1v <- values(lam1); X2v <- values(lam2); X3v <- values(lam3)
  inside <- which(values(maskTri))
  U <- V <- rep(NA_real_, length(X1v))
  
  if (length(inside) > 0) {
    for (idx in inside) {
      X <- c(X1v[idx], X2v[idx], X3v[idx])
      Xn <- step_map(X, params)
      dX <- Xn - X
      u1 <- dX[1]; v1 <- dX[3]
      U[idx] <- (x1 - x3)*u1 + (x2 - x3)*v1
      V[idx] <- (y1 - y3)*u1 + (y2 - y3)*v1
    }
  }
  
  field <- stack(u, v)
  names(field) <- c('u', 'v')
  values(field[[1]]) <- U
  values(field[[2]]) <- V
  
  if (normalize) {
    mag <- sqrt(U^2 + V^2)
    maxmag <- suppressWarnings(max(mag, na.rm = TRUE))
    if (is.finite(maxmag) && maxmag > 0) {
      values(field[[1]]) <- U / maxmag
      values(field[[2]]) <- V / maxmag
    }
  }
  field
}

# -------------
# Server
# -------------
server <- function(input, output, session) {
  # Reset button
  observeEvent(input$reset, {
    updateSliderInput(session, "g1", value = 1)
    updateSliderInput(session, "g2", value = 1)
    updateSliderInput(session, "g3", value = 1)
    updateSliderInput(session, "KD", value = 10)
    updateSliderInput(session, "KP", value = 15)
    updateSliderInput(session, "cres1", value = 0.05)
    updateSliderInput(session, "cres2", value = 0.05)
    updateSliderInput(session, "cdeg1", value = 0.105)
    updateSliderInput(session, "cdeg2", value = 0.105)
    updateSliderInput(session, "cprod1", value = 0.15)
    updateSliderInput(session, "cprod2", value = 0.15)
    updateSliderInput(session, "x1", value = 1/3)
    updateSliderInput(session, "x2", value = 1/3)
    updateSliderInput(session, "T",  value = 150)
  })
  
  params <- reactive({
    list(
      g     = c(input$g1, input$g2, input$g3),
      KD    = input$KD,
      KP    = input$KP,
      cres  = c(input$cres1, input$cres2),
      cdeg  = c(input$cdeg1, input$cdeg2),
      cprod = c(input$cprod1, input$cprod2)
    )
  })
  
  x3 <- reactive({
    x3val <- 1 - input$x1 - input$x2
    if (x3val < 0) NA_real_ else x3val
  })
  
  output$x3txt <- renderUI({
    if (is.na(x3())) {
      HTML("<b>Warning:</b> X1 + X2 > 1. Please reduce one of them so that X3 ≥ 0.")
    } else {
      HTML(paste0("Computed X<sub>3</sub> = ", sprintf("%.3f", x3())))
    }
  })
  
  # Trajectory under the discrete map
  trajectory <- reactive({
    req(!is.na(x3()))
    Tn <- input$T
    X0 <- c(input$x1, input$x2, x3())
    X0 <- pmax(X0, 0); X0 <- X0 / sum(X0)
    P  <- params()
    
    out <- matrix(NA_real_, nrow = Tn + 1, ncol = 3)
    out[1, ] <- X0
    for (t in 1:Tn) out[t + 1, ] <- step_map(out[t, ], P)
    df <- as.data.frame(out)
    names(df) <- c("X1", "X2", "X3")
    df$step <- 0:Tn
    df
  })
  
  output$ternplot <- renderPlot({
    df <- trajectory()
    Tn <- nrow(df) - 1
    
    ggplot(df, aes(x = X1, y = X2, z = X3)) +
      coord_tern() +
      geom_path(size = 1, alpha = 0.85) +
      geom_point(data = df[1,], size = 3, shape = 21, fill = NA, stroke = 1.2) +
      geom_point(data = df[nrow(df),], size = 3) +
      labs(T = "Strain 1", L = "Strain 2", R = "Strain 3") +
      theme_bw() +
      theme(
        tern.panel.grid.major = element_line(linetype = "dotted"),
        plot.title = element_text(face = "bold", size = 16)
      ) +
      labs(title = "Trajectory on the simplex", subtitle = paste0("Iterations: ", Tn))
  })
  
  # Removed stray ggtern block that caused errors
  
  output$matrices <- renderTable({
    list(S = S, P = Pmat, D = Dmat, R = Rmat) %>%
      lapply(as.data.frame) %>% bind_rows(.id = "Matrix")
  }, rownames = TRUE)
  
  # Flow / streamdensity tab
  output$flowplot <- renderPlot({
    P <- params()
    field <- compute_flow_field(P, Ngrid = input$gridN, normalize = input$normvel)
    
    my_colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100)
    color_breaks <- seq(0, 1, length.out = 100)
    
    vp <- vectorplot(field[[c(1,2)]],
                     ylab = NULL, xlab = NULL,
                     isField = 'dXY',
                     lwd.arrows = 1.2,
                     scaleSlope = FALSE,
                     arrow.length = 0.08,
                     col.arrows = '#303030',
                     interpolate = FALSE,
                     maxpixels = 1e5,
                     narrows = input$arrows,
                     col.regions = my_colors,
                     at = color_breaks,
                     colorkey = FALSE,
                     par.settings = modifyList(PuOrTheme(), list(
                       axis.line = list(col = "transparent")
                     )),
                     scales = list(alternating = 0),
                     xlim = c(-0.05, 1.05),
                     ylim = c(-0.05, 1.05),
                     aspX = 0.2,
                     aspY = 0.2,
                     alpha = 0.7)
    
    original_panel <- vp$panel
    vp$panel <- function(...) {
      original_panel(...)
      panel.lines(x = c(0, 0.5, 1, 0),
                  y = c(0, sqrt(3)/2, 0, 0),
                  col = "black", lwd = 2.5)
      grid.text("Strain 1", x = unit(0, "native"), y = unit(-0.05, "native"),
                just = c("center","top"), gp = gpar(cex = 1.1))
      grid.text("Strain 2", x = unit(0.5, "native"), y = unit(sqrt(3)/2 + 0.04, "native"),
                just = c("center","bottom"), gp = gpar(cex = 1.1))
      grid.text("Strain 3", x = unit(1, "native"), y = unit(-0.05, "native"),
                just = c("center","top"), gp = gpar(cex = 1.1))
    }
    
    grid.newpage()
    vp_grob <- grid.grabExpr(print(vp))
    
    ckey <- draw.colorkey(
      key = list(
        col = my_colors,
        at = color_breaks,
        labels = list(cex = 0.8),
        height = 0.3, width = 1, space = "right"
      ), draw = FALSE)
    title_grob <- textGrob("Relative speed |ΔX|", rot = 90, gp = gpar(cex = 0.9))
    legend_grob <- arrangeGrob(title_grob, ckey, ncol = 2, widths = unit(c(1.5, 5), "lines"))
    
    grid.arrange(vp_grob, legend_grob, ncol = 2,
                 widths = unit.c(unit(1, "null"), unit(6, "lines")))
  })
}

shinyApp(ui, server)

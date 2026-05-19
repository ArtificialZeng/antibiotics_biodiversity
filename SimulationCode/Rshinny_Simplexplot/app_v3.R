# Minimal Shiny simplex with parameter sliders (ggplot2 only)
# -----------------------------------------------------------
# Draws a ternary simplex with a vector field of the discrete-time map
# for 3 strains and 2 antibiotics. Users can tune g, K_D, K_P, and costs c_*.

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(grid)
})

# -------------------------
# Fixed phenotype matrices (edit if needed)
# -------------------------
# S <- matrix(1, nrow = 3, ncol = 2, byrow = TRUE)           # susceptibility
   S <- matrix(c(0,0,
                 1,0,
                 0,1), nrow = 3, byrow = TRUE)             # producers
Pmat <- matrix(c(1,0,
                 0,1,
                 0,0), nrow = 3, byrow = TRUE)             # producers
Dmat <- matrix(c(0,1,
                 0,0,
                 0,0), nrow = 3, byrow = TRUE)             # degraders
Rmat <- matrix(c(0,0,
                 0,0,
                 0,0), nrow = 3, byrow = TRUE)             # intrinsic resistance

# Triangle vertices (Cartesian)
V1 <- c(0, 0)            # X1 = 1
V2 <- c(1, 0)            # X2 = 1
V3 <- c(0.5, sqrt(3)/2)  # X3 = 1

# Barycentric -> XY
bary_to_xy <- function(l1, l2, l3) {
  x <- l1*V1[1] + l2*V2[1] + l3*V3[1]
  y <- l1*V1[2] + l2*V2[2] + l3*V3[2]
  c(x, y)
}

# Discrete-time map X' = (X * f) / <X,f>
step_map <- function(X, pars) {
  g     <- pars$g
  KD    <- pars$KD
  KP    <- pars$KP
  cres  <- pars$cres
  cdeg  <- pars$cdeg
  cprod <- pars$cprod
  
  r <- numeric(3)
  for (i in 1:3) {
    r[i] <- g[i] - sum(Pmat[i,] * cprod) - sum(Dmat[i,] * cdeg) - sum(Rmat[i,] * cres)
  }
  Aj <- as.numeric(t(X) %*% Dmat)
  Bj <- as.numeric(t(X) %*% Pmat)
  
  pij <- matrix(0, nrow = 3, ncol = 2)
  for (i in 1:3) for (j in 1:2)
    pij[i,j] <- S[i,j] * exp(-KD * Aj[j]) * (1 - exp(-KP * Bj[j]))
  
  fi <- r * (1 - pij[,1]) * (1 - pij[,2])
  fi[fi < 0] <- 0
  denom <- sum(X * fi)
  if (!is.finite(denom) || denom <= 1e-12) return(X)
  as.numeric(X * fi / denom)
}

# Grid of barycentric points inside the triangle
make_bary_grid <- function(n = 28) {
  pts <- list(); k <- 1
  for (i in 0:n) {
    for (j in 0:(n - i)) {
      l1 <- i / n; l2 <- j / n; l3 <- 1 - l1 - l2
      pts[[k]] <- c(l1, l2, l3); k <- k + 1
    }
  }
  do.call(rbind, pts)
}

# Build arrow segments for current parameters
compute_arrows <- function(pars, n = 28, arrow_scale = 0.18) {
  B <- make_bary_grid(n)
  XY <- t(apply(B, 1, function(b) bary_to_xy(b[1], b[2], b[3])))
  Xn <- t(apply(B, 1, function(b) step_map(b, pars)))
  dX <- Xn - B
  # Map differential in barycentric to Cartesian: dXY = (V1-V3)*dl1 + (V2-V3)*dl2
  dL1 <- dX[,1]; dL2 <- dX[,2]
  u <- (V1[1] - V3[1]) * dL1 + (V2[1] - V3[1]) * dL2
  v <- (V1[2] - V3[2]) * dL1 + (V2[2] - V3[2]) * dL2
  mag <- sqrt(u^2 + v^2); maxmag <- max(mag, na.rm = TRUE)
  if (!is.finite(maxmag) || maxmag == 0) maxmag <- 1
  s <- arrow_scale / maxmag
  data.frame(x = XY[,1], y = XY[,2], xend = XY[,1] + s*u, yend = XY[,2] + s*v)
}

# -------------------------
# UI with sliders for g, KD/KP, and costs per antibiotic
# -------------------------
ui <- fluidPage(
  titlePanel("Simplex flow (tunable parameters)"),
  sidebarLayout(
    sidebarPanel(
      h4("Growth rates"),
      sliderInput("g1", "g1", min = 0, max = 2, value = 1, step = 0.01),
      sliderInput("g2", "g2", min = 0, max = 2, value = 1, step = 0.01),
      sliderInput("g3", "g3", min = 0, max = 2, value = 1, step = 0.01),
      tags$hr(),
      h4("Interaction scales"),
      sliderInput("KD", "K_D (degradation radius)", min = 0, max = 30, value = 10, step = 0.1),
      sliderInput("KP", "K_P (killing radius)",    min = 0, max = 30, value = 15, step = 0.1),
      tags$hr(),
      h4("Metabolic costs (per antibiotic)"),
      fluidRow(
        column(4, sliderInput("cres1",  "c_res[1]",  min = 0, max = 0.5, value = 0.05,  step = 0.005)),
        column(4, sliderInput("cdeg1",  "c_deg[1]",  min = 0, max = 0.5, value = 0.105, step = 0.005)),
        column(4, sliderInput("cprod1", "c_prod[1]", min = 0, max = 0.5, value = 0.15,  step = 0.005))
      ),
      fluidRow(
        column(4, sliderInput("cres2",  "c_res[2]",  min = 0, max = 0.5, value = 0.05,  step = 0.005)),
        column(4, sliderInput("cdeg2",  "c_deg[2]",  min = 0, max = 0.5, value = 0.105, step = 0.005)),
        column(4, sliderInput("cprod2", "c_prod[2]", min = 0, max = 0.5, value = 0.15,  step = 0.005))
      ),
      tags$hr(),
      actionButton("reset", "Reset defaults")
    ),
    mainPanel(
      plotOutput("simplex", height = "640px")
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  observeEvent(input$reset, {
    updateSliderInput(session, "g1", value = 1)
    updateSliderInput(session, "g2", value = 1)
    updateSliderInput(session, "g3", value = 1)
    updateSliderInput(session, "KD", value = 10)
    updateSliderInput(session, "KP", value = 15)
    updateSliderInput(session, "cres1", value = 0.05)
    updateSliderInput(session, "cdeg1", value = 0.105)
    updateSliderInput(session, "cprod1", value = 0.15)
    updateSliderInput(session, "cres2", value = 0.05)
    updateSliderInput(session, "cdeg2", value = 0.105)
    updateSliderInput(session, "cprod2", value = 0.15)
  })
  
  pars <- reactive({
    list(
      g     = c(input$g1, input$g2, input$g3),
      KD    = input$KD,
      KP    = input$KP,
      cres  = c(input$cres1, input$cres2),
      cdeg  = c(input$cdeg1, input$cdeg2),
      cprod = c(input$cprod1, input$cprod2)
    )
  })
  
  output$simplex <- renderPlot({
    segs <- compute_arrows(pars(), n = 28, arrow_scale = 0.18)
    
    tri_df <- data.frame(
      x = c(V1[1], V3[1], V2[1], V1[1]),
      y = c(V1[2], V3[2], V2[2], V1[2])
    )
    
    ggplot() +
      geom_polygon(data = tri_df, aes(x, y), fill = NA, color = "black", size = 1) +
      geom_segment(data = segs, aes(x = x, y = y, xend = xend, yend = yend),
                   size = 0.4,
                   arrow = arrow(length = unit(3, "pt"), type = "closed"),
                   alpha = 0.85) +
      annotate("text", x = V1[1], y = V1[2] - 0.04, label = "Strain 1") +
      annotate("text", x = V3[1], y = V3[2] + 0.04, label = "Strain 2") +
      annotate("text", x = V2[1], y = V2[2] - 0.04, label = "Strain 3") +
      coord_equal(xlim = c(-0.05, 1.05), ylim = c(-0.05, sqrt(3)/2 + 0.08), expand = FALSE) +
      theme_void()
  })
}

shinyApp(ui, server)

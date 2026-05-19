# Shiny simplex with parameter sliders + user trajectory (ggplot2 only)
# -------------------------------------------------------------------
# - Draws a ternary simplex with a vector field (fixed resolution)
# - Lets user set g, K_D, K_P, c_* via sliders
# - Lets user set an initial frequency (X1,X2,X3) and number of iterations
#   and plots the discrete-time trajectory on the simplex

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(grid)
})

# -------------------------
# Phenotype matrices (as requested)
# -------------------------
S <- matrix(c(0,0,
              1,0,
              0,1), nrow = 3, byrow = TRUE)           # susceptibility
Pmat <- matrix(c(1,0,
                 0,1,
                 0,0), nrow = 3, byrow = TRUE)         # producers
Dmat <- matrix(c(0,1,
                 0,0,
                 0,0), nrow = 3, byrow = TRUE)         # degraders
Rmat <- matrix(c(0,0,
                 0,0,
                 0,0), nrow = 3, byrow = TRUE)         # intrinsic resistance

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
  B  <- make_bary_grid(n)
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

# Continuous background density via regular XY grid inside triangle
compute_density_raster <- function(pars, Nx = 220, Ny = 190) {
  # Regular grid covering the triangle's bounding box
  xs <- seq(0, 1, length.out = Nx)
  ys <- seq(0, sqrt(3)/2, length.out = Ny)
  grid <- expand.grid(x = xs, y = ys)
  
  # Barycentric coords for each grid point
  x1 <- V1[1]; y1 <- V1[2]
  x2 <- V2[1]; y2 <- V2[2]
  x3 <- V3[1]; y3 <- V3[2]
  denom <- (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
  l1 <- ((y2 - y3)*(grid$x - x3) + (x3 - x2)*(grid$y - y3)) / denom
  l3 <- ((y3 - y1)*(grid$x - x3) + (x1 - x3)*(grid$y - y3)) / denom
  l2 <- 1 - l1 - l3
  
  inside <- (l1 >= 0) & (l2 >= 0) & (l3 >= 0)
  grid <- grid[inside, , drop = FALSE]
  l1 <- l1[inside]; l2 <- l2[inside]; l3 <- l3[inside]
  
  # Speed at each inside point
  step_speed <- function(a,b,c) {
    X  <- c(a,b,c)
    Xn <- step_map(X, pars)
    dX <- Xn - X
    dL1 <- dX[1]; dL2 <- dX[2]
    u <- (V1[1] - V3[1]) * dL1 + (V2[1] - V3[1]) * dL2
    v <- (V1[2] - V3[2]) * dL1 + (V2[2] - V3[2]) * dL2
    sqrt(u*u + v*v)
  }
  mags <- mapply(step_speed, l1, l2, l3)
  r <- range(mags, finite = TRUE)
  mags <- (mags - r[1]) / (r[2] - r[1] + 1e-12)
  data.frame(x = grid$x, y = grid$y, mag = mags)
}

# Compute trajectory from X0 for T steps
compute_trajectory <- function(pars, X0, steps = 150) {
  X <- as.numeric(X0)
  # ensure valid simplex point (clip negatives, renorm)
  X[X < 0] <- 0
  s <- sum(X); if (!is.finite(s) || s == 0) X <- c(1,0,0) else X <- X / s
  out <- matrix(NA_real_, nrow = steps + 1, ncol = 3)
  out[1, ] <- X
  for (t in 1:steps) {
    X <- step_map(X, pars)
    out[t + 1, ] <- X
  }
  XY <- t(apply(out, 1, function(b) bary_to_xy(b[1], b[2], b[3])))
  data.frame(x = XY[,1], y = XY[,2], step = 0:steps)
}

# -------------------------
# UI with sliders for g, KD/KP, costs, and trajectory inputs
# -------------------------
ui <- fluidPage(
  titlePanel("Three-Strain / Two-Antibiotic Frequency Dynamics"),
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
      h4("Trajectory initial condition & time"),
      numericInput("x1", "Initial X1", value = 0.33, min = 0, max = 1, step = 0.01),
      numericInput("x2", "Initial X2", value = 0.33, min = 0, max = 1, step = 0.01),
      # numericInput("x3", "Initial X3", value = 0.34, min = 0, max = 1, step = 0.01),
      numericInput("steps", "Iterations (T)", value = 150, min = 1, max = 1000, step = 1),
      helpText("Initial values are automatically normalized to sum to 1."),
      tags$hr(),
      actionButton("reset", "Reset defaults"),
      tags$hr(),
      h4("Color scale"),
      selectInput("pal", "ColorBrewer palette (colorblind-friendly)",
                  choices = c("PuBuGn","Blues","Purples","Greys","Oranges","YlGnBu"),
                  selected = "YlGnBu")
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
    updateNumericInput(session, "x1", value = 0.33)
    updateNumericInput(session, "x2", value = 0.33)
    # updateNumericInput(session, "x3", value = 0.34)
    updateNumericInput(session, "steps", value = 150)
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
    P <- pars()
    segs <- compute_arrows(P, n = 28, arrow_scale = 0.18)
    
    # Trajectory path from user-specified initial condition
    X0 <- c(input$x1, input$x2, 1-input$x1-input$x2)
    traj <- compute_trajectory(P, X0, steps = input$steps)
    
    # mark start (open circle) and end (closed gray)
    start_pt <- traj[1, , drop = FALSE]
    end_pt   <- traj[nrow(traj), , drop = FALSE]
    
    # Continuous speed background
    dens <- compute_density_raster(P, Nx = 200, Ny = 200)
    
    tri_df <- data.frame(
      x = c(V1[1], V3[1], V2[1], V1[1]),
      y = c(V1[2], V3[2], V2[2], V1[2])
    )
    
    ggplot() +
      # background color density (smooth)
      geom_raster(data = dens, aes(x = x, y = y, fill = mag)) +
      scale_fill_distiller(palette = input$pal, type = "seq", direction = 1, name = "speed") +
      # scale_fill_viridis_c(option = "C", name = "speed") +
      # triangle outline
      geom_polygon(data = tri_df, aes(x, y), fill = NA, color = "black", size = 1) +
      # flow field arrows
      geom_segment(data = segs, aes(x = x, y = y, xend = xend, yend = yend),
                   size = 0.4,
                   arrow = grid::arrow(length = grid::unit(3, "pt"), type = "closed"),
                   alpha = 0.75) +
      # trajectory path
      geom_path(data = traj, aes(x = x, y = y), linewidth = 1.1, alpha = 0.95) +
      # start: open circle; end: closed gray circle
      geom_point(data = start_pt, aes(x = x, y = y),
                 shape = 21, size = 3.8, stroke = 1.1, fill = NA, color = "black") +
      geom_point(data = end_pt, aes(x = x, y = y),
                 shape = 21, size = 3.8, stroke = 0.5, fill = "gray40", color = "gray40") +
      annotate("text", x = V1[1], y = V1[2] - 0.04, label = "Strain 1") +
      annotate("text", x = V3[1], y = V3[2] + 0.04, label = "Strain 2") +
      annotate("text", x = V2[1], y = V2[2] - 0.04, label = "Strain 3") +
      coord_equal(xlim = c(-0.05, 1.05), ylim = c(-0.05, sqrt(3)/2 + 0.08), expand = FALSE) +
      theme_void() +
      scale_color_viridis_c(option = "C", name = "speed")
  })
}

shinyApp(ui, server)

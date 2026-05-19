# Three-Strain / Two-Antibiotic Frequency Dynamics (Replicator Map)
# ---------------------------------------------------------------
# This Shiny app visualizes the discrete-time frequency dynamics of N=3 strains
# under M=2 antibiotics on a ternary (simplex-3) diagram. The map is
#   X_{t+1,i} = X_{t,i} * f_i(X_t) / sum_k X_{t,k} f_k(X_t)
# with
#   f_i = r_i * prod_{j=1..2} [ 1 - S_{ij} * exp(-K_D * A_j) * (1 - exp(-K_P * B_j)) ]
#   r_i = g_i - sum_j ( P_{ij}*c_prod_j + D_{ij}*c_deg_j + R_{ij}*c_res_j )
#   A_j = sum_k D_{kj} X_k,  B_j = sum_k P_{kj} X_k.
# Matrices S, P, D, R are set below (editable in code). Sliders control g, K_D, K_P, and cost vectors c_*.

suppressPackageStartupMessages({
  library(shiny)
  library(ggtern)
  library(dplyr)
})

# -------------------------
# Phenotype matrices (edit)
# -------------------------
# 3 strains x 2 antibiotics matrices.
# S: susceptibility coefficients (default all 1 = fully susceptible)
S <- matrix(1, nrow = 3, ncol = 2, byrow = TRUE,
            dimnames = list(paste0("strain", 1:3), paste0("ab", 1:2)))

# P: producer phenotype (1 if produces antibiotic j, else 0)
# (demo) strain1 produces ab1; strain3 produces ab2
Pmat <- matrix(c(1,0,
                 0,0,
                 0,1), nrow = 3, byrow = TRUE,
               dimnames = list(paste0("strain", 1:3), paste0("ab", 1:2)))

# D: degrader phenotype (1 if degrades antibiotic j, else 0)
# (demo) strain2 degrades ab1
Dmat <- matrix(c(0,0,
                 1,0,
                 0,0), nrow = 3, byrow = TRUE,
               dimnames = list(paste0("strain", 1:3), paste0("ab", 1:2)))

# R: intrinsic resistance phenotype (1 if resistant to antibiotic j, else 0)
# (demo) strain2 resistant to ab2; strain3 resistant to ab1
Rmat <- matrix(c(0,0,
                 0,1,
                 1,0), nrow = 3, byrow = TRUE,
               dimnames = list(paste0("strain", 1:3), paste0("ab", 1:2)))

# -------------
# Shiny UI
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
      div("X3 is computed as 1 - X1 - X2 and must be \u2265 0."),
      tags$hr(),
      sliderInput("T", "Iterations (discrete steps)", min = 1, max = 400, value = 150, step = 1),
      actionButton("reset", "Reset defaults")
    ),
    mainPanel(
      fluidRow(
        column(12,
               plotOutput("ternplot", height = "650px"),
               br(),
               htmlOutput("x3txt"),
               br(),
               h4("Phenotype matrices (fixed in this demo; edit in code)"),
               tableOutput("matrices")
        )
      )
    )
  )
)

# -------------
# Helper functions
# -------------
clip01 <- function(x) pmin(pmax(x, 0), 1)

# One-step update of the discrete replicator map
step_map <- function(X, params) {
  # Unpack
  g    <- params$g
  KD   <- params$KD
  KP   <- params$KP
  cres <- params$cres
  cdeg <- params$cdeg
  cprod<- params$cprod
  
  # Growth rates corrected for metabolic costs
  r <- numeric(3)
  for (i in 1:3) {
    r[i] <- g[i] - sum(Pmat[i,] * cprod) - sum(Dmat[i,] * cdeg) - sum(Rmat[i,] * cres)
  }
  
  # Global totals for each antibiotic
  Aj <- as.numeric(t(X) %*% Dmat)  # length 2: A1, A2
  Bj <- as.numeric(t(X) %*% Pmat)  # length 2: B1, B2
  
  # Kill probabilities for each (i,j)
  pij <- matrix(0, nrow = 3, ncol = 2)
  for (i in 1:3) {
    for (j in 1:2) {
      pij[i, j] <- S[i, j] * exp(-KD * Aj[j]) * (1 - exp(-KP * Bj[j]))
    }
  }
  
  # Fitnesses f_i (M=2 => product over j=1..2)
  fi <- r * (1 - pij[,1]) * (1 - pij[,2])
  fi[fi < 0] <- 0  # guard against negative fitness
  
  denom <- sum(X * fi)
  if (!is.finite(denom) || denom <= 1e-12) return(X)
  
  Xnext <- X * fi / denom
  as.numeric(Xnext)
}

# -------------
# Shiny server
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
      HTML("<b>Warning:</b> X1 + X2 > 1. Please reduce one of them so that X3 \u2265 0.")
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
    for (t in 1:Tn) {
      out[t + 1, ] <- step_map(out[t, ], P)
    }
    df <- as.data.frame(out)
    names(df) <- c("X1", "X2", "X3")
    df$step <- 0:Tn
    df
  })
  
  output$ternplot <- renderPlot({
    df <- trajectory()
    Tn <- nrow(df) - 1
    
    ggtern(df, aes(X1, X2, X3)) +
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
  
  output$matrices <- renderTable({
    list(
      S = S,
      P = Pmat,
      D = Dmat,
      R = Rmat
    ) |> lapply(as.data.frame) |> bind_rows(.id = "Matrix")
  }, rownames = TRUE)
}

shinyApp(ui, server)

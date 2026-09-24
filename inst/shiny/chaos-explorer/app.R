# chaosrcpp explorer: a small Shiny front end for the package.
# Launch with chaosrcpp::run_chaos_explorer().

library(shiny)
library(chaosrcpp)

ui <- navbarPage(
  "chaosrcpp explorer",
  theme = NULL,

  # ---------------------------------------------------------------- Lorenz
  tabPanel(
    "Lorenz",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        sliderInput("sigma", HTML("&sigma;"), min = 1, max = 30, value = 10, step = 0.5),
        sliderInput("rho", HTML("&rho;"), min = 0.5, max = 60, value = 28, step = 0.5),
        sliderInput("beta", HTML("&beta;"), min = 0.5, max = 6, value = 8 / 3, step = 0.05),
        sliderInput("tmax", "t max", min = 5, max = 60, value = 30, step = 1),
        radioButtons("lorenz_view", "View",
                     c("3D attractor (fp64)" = "attractor",
                       "fp32 vs fp64 from the same start" = "divergence")),
        helpText("Both runs use the same RK4 code and step size; only the",
                 "floating-point type differs. Watch the two curves part ways.")
      ),
      mainPanel(
        width = 9,
        conditionalPanel("input.lorenz_view == 'attractor'",
                         plotly::plotlyOutput("lorenz3d", height = "620px")),
        conditionalPanel("input.lorenz_view == 'divergence'",
                         plotly::plotlyOutput("lorenz_div3d", height = "400px"),
                         plotOutput("lorenz_div", height = "260px"))
      )
    )
  ),

  # ---------------------------------------------------------- Logistic map
  tabPanel(
    "Logistic map",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        sliderInput("r", "r", min = 2.5, max = 4, value = 3.7, step = 0.001),
        sliderInput("x0", HTML("x<sub>0</sub>"), min = 0.01, max = 0.99, value = 0.1, step = 0.01),
        sliderInput("nsteps", "iterations shown", min = 10, max = 200, value = 60, step = 10),
        helpText("The vertical line on the bifurcation diagram marks the current r.",
                 "The Lyapunov exponent below it is positive exactly where the map is chaotic.")
      ),
      mainPanel(
        width = 9,
        fluidRow(
          column(6, plotOutput("cobweb", height = "360px")),
          column(6, plotOutput("orbit", height = "360px"))
        ),
        plotOutput("bifurcation", height = "380px")
      )
    )
  ),

  # ---------------------------------------------------------- Standard map
  tabPanel(
    "Standard map",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        sliderInput("K", "K (kick strength)", min = 0, max = 3, value = 0.97, step = 0.01),
        sliderInput("norbits", "orbits", min = 10, max = 150, value = 60, step = 10),
        sliderInput("niter", "iterates per orbit", min = 100, max = 3000, value = 800, step = 100),
        helpText("K = 0 is integrable (horizontal lines). Around K = 0.97 the last",
                 "invariant torus breaks and the chaotic sea connects.")
      ),
      mainPanel(width = 9, plotOutput("standard", height = "640px"))
    )
  )
)

server <- function(input, output, session) {

  # ---- Lorenz ------------------------------------------------------------
  lorenz_fp64 <- reactive({
    lorenz(sigma = input$sigma, rho = input$rho, beta = input$beta,
           dt = 0.01, n = round(input$tmax / 0.01), precision = "fp64")
  })
  lorenz_fp32 <- reactive({
    lorenz(sigma = input$sigma, rho = input$rho, beta = input$beta,
           dt = 0.01, n = round(input$tmax / 0.01), precision = "fp32")
  })

  output$lorenz3d <- plotly::renderPlotly({
    plot_attractor_3d(lorenz_fp64())
  })

  output$lorenz_div3d <- plotly::renderPlotly({
    plot_attractor_3d(fp32 = lorenz_fp32(), fp64 = lorenz_fp64())
  })

  output$lorenz_div <- renderPlot({
    a <- lorenz_fp32(); b <- lorenz_fp64()
    d <- sqrt((a$x - b$x)^2 + (a$y - b$y)^2 + (a$z - b$z)^2)
    ok <- d > 0
    par(mar = c(4, 4.5, 1, 1))
    plot(a$t[ok], d[ok], type = "l", log = "y", col = "#D55E00", lwd = 1.5,
         xlab = "t", ylab = "|fp32 - fp64|", xlim = c(0, input$tmax))
    abline(h = 1, lty = 3, col = "grey40")
    lines(a$t, pmin(2^-24 * exp(0.9056 * a$t), 100), lty = 2, col = "grey30")
    legend("bottomright", bty = "n",
           legend = c("distance between the two runs", "2^-24 * exp(0.906 t)"),
           col = c("#D55E00", "grey30"), lty = c(1, 2), lwd = c(1.5, 1))
  })

  # ---- Logistic map ------------------------------------------------------
  bif <- logistic_bifurcation(seq(2.5, 4, length.out = 1200), n_transient = 300, n_keep = 120)
  lyap <- logistic_lyapunov(seq(2.5, 4, length.out = 1200), n_transient = 300, n = 1500)

  output$cobweb <- renderPlot({
    cw <- logistic_cobweb(input$r, input$x0, input$nsteps)
    xs <- seq(0, 1, length.out = 400)
    par(mar = c(4, 4.5, 2.5, 1), pty = "s")
    plot(xs, input$r * xs * (1 - xs), type = "l", lwd = 2, col = "#1f2a44",
         xlim = c(0, 1), ylim = c(0, 1), xlab = expression(x[n]), ylab = expression(x[n + 1]),
         main = sprintf("Cobweb, r = %.3f", input$r))
    abline(0, 1, col = "grey55")
    s <- cw$segments
    cols <- hcl.colors(input$nsteps, "Plasma", rev = TRUE)
    segments(s$x, s$y, s$xend, s$yend, col = cols[s$step], lwd = 0.9)
  })

  output$orbit <- renderPlot({
    x <- logistic_orbit(input$r, input$x0, input$nsteps)
    par(mar = c(4, 4.5, 2.5, 1))
    plot(seq_along(x) - 1, x, type = "b", pch = 16, cex = 0.7, col = "#0072B2",
         ylim = c(0, 1), xlab = "n", ylab = expression(x[n]), main = "Orbit")
  })

  output$bifurcation <- renderPlot({
    lam <- lyap$lambda[which.min(abs(lyap$r - input$r))]
    layout(matrix(1:2, 2), heights = c(2, 1))
    par(mar = c(0.5, 4.5, 2.5, 1))
    plot(bif$r, bif$x, pch = ".", col = adjustcolor("#1f2a44", 0.4), xaxt = "n",
         xlab = "", ylab = "x", xaxs = "i",
         main = sprintf("Bifurcation diagram, lambda(r = %.3f) = %.3f", input$r, lam))
    abline(v = input$r, col = "#D55E00", lwd = 2)
    par(mar = c(4, 4.5, 0.5, 1))
    plot(lyap$r, lyap$lambda, type = "l", lwd = 0.8, col = "#1f2a44", xaxs = "i",
         ylim = c(-3, 0.8), xlab = "r", ylab = expression(lambda))
    polygon(c(lyap$r, rev(lyap$r)), c(pmax(lyap$lambda, 0), rep(0, nrow(lyap))),
            col = adjustcolor("#D55E00", 0.3), border = NA)
    abline(h = 0, col = "grey40")
    abline(v = input$r, col = "#D55E00", lwd = 2)
  })

  # ---- Standard map ------------------------------------------------------
  output$standard <- renderPlot({
    sm <- standard_map(K = input$K, n_orbits = input$norbits, n = input$niter, seed = 1)
    par(bg = "white")
    plot(sm)
  })
}

shinyApp(ui, server)

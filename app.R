# ==============================================================================
# App Shiny: atípicos, alta palanca y observaciones influyentes
#            con dos variables explicativas  y = b0 + b1 x1 + b2 x2 + e
#
# Ejecutar:  install.packages("shiny")   # solo la primera vez
#            shiny::runApp("ruta/a/la/carpeta/de/app.R")
#
# Solo usa shiny y gráficos base de R.
# ==============================================================================

library(shiny)

# ---- Parámetros "verdaderos" del modelo ---------------------------------------
b0r <- 2; b1r <- 1; b2r <- 0.8
pred <- function(a, b) b0r + b1r * a + b2r * b
p <- 3

# ---- Recorrido automático (igual que en el GIF) -------------------------------
ini <- c(5, 5); fin <- c(22, 19); A <- 9
ta <- seq(0, 1, length.out = 40)
sb <- seq(0, 1, length.out = 70)
tc <- seq(0, 1, length.out = 30)
recorrido <- rbind(
  data.frame(x1 = ini[1] + ta * (fin[1] - ini[1]), x2 = ini[2] + ta * (fin[2] - ini[2]), d = 0, fase = 1),
  data.frame(x1 = fin[1], x2 = fin[2], d = A * sin(2 * pi * sb), fase = 2),
  data.frame(x1 = fin[1] + tc * (ini[1] - fin[1]), x2 = fin[2] + tc * (ini[2] - fin[2]), d = 0, fase = 3)
)
textos_fase <- c(
  "Fase 1 · El punto se aleja en (x₁, x₂) siguiendo la tendencia: alta palanca, poca atipicidad.",
  "Fase 2 · Lejos en (x₁, x₂) y fuera de la tendencia: palanca + atipicidad = influencia.",
  "Fase 3 · Vuelve al centro siguiendo la tendencia: la influencia desaparece."
)

# ---- Utilidades ---------------------------------------------------------------
unit <- function(v) v / sqrt(sum(v^2))

col_cook <- function(D, U, activos) {
  supera <- c(F = D > U[["F"]], uno = D > U[["uno"]], cuatro = D > U[["cuatro"]])
  k <- sum(supera[activos])
  if (k == 0) "#27ae60" else if (k < length(activos)) "#f39c12" else "#c0392b"
}

llave <- function(x0, x1, y, alto, ...) {
  w <- x1 - x0
  t <- seq(0, 1, length.out = 120)
  f <- 1 / (1 + exp(-18 * t)) + 1 / (1 + exp(-18 * (t - 1))) - 0.5
  lines(x0 + c(t * w / 2, w - rev(t) * w / 2), y - alto * c(f, rev(f)), ...)
  c(x0 + w / 2, y - alto)
}

# ---- Ajuste del modelo con el punto móvil -------------------------------------
ajustar <- function(base, x1n, x2n, yn) {
  x1 <- c(base$x1, x1n); x2 <- c(base$x2, x2n); y <- c(base$y, yn)
  m  <- lm(y ~ x1 + x2)
  n  <- length(y)
  U  <- c(F = qf(0.5, p, n - p), uno = 1, cuatro = 4 / (n - p - 1))
  list(x1 = x1, x2 = x2, y = y, n = n, m = m,
       beta = coef(m), beta_sin = coef(lm(base$y ~ base$x1 + base$x2)),
       yhat = unname(fitted(m)), e = unname(residuals(m)), h = unname(hatvalues(m)),
       r = unname(rstandard(m)), D = unname(cooks.distance(m)), U = U,
       sigma = summary(m)$sigma)
}

# ==============================================================================
# Gráficos
# ==============================================================================
g_variables <- function(a, theta, phi, mostrar_sin, colp) {
  n <- a$n
  g1 <- range(c(-1, 26, a$x1)); g2 <- range(c(-1, 26, a$x2))
  zr <- range(c(a$y, pred(g1, g2) + c(-12, 12)))
  plano <- function(b) outer(g1, g2, function(u, v) b[1] + b[2] * u + b[3] * v)
  par(mar = c(1, 1, 2.5, 1))
  pm <- persp(g1, g2, plano(a$beta), zlim = zr, theta = theta, phi = phi,
              col = adjustcolor("#e74c3c", 0.15), border = "#e74c3c", lwd = 1.6,
              xlab = "x₁", ylab = "x₂", zlab = "y", expand = 0.85)
  if (mostrar_sin) {
    zs <- plano(a$beta_sin)
    q <- trans3d(c(g1[1], g1[2], g1[2], g1[1], g1[1]), c(g2[1], g2[1], g2[2], g2[2], g2[1]),
                 c(zs[1, 1], zs[2, 1], zs[2, 2], zs[1, 2], zs[1, 1]), pm)
    polygon(q$x, q$y, col = adjustcolor("#7f8c8d", 0.2), border = "#555555", lty = 2, lwd = 1.4)
  }
  pb <- trans3d(a$x1[-n], a$x2[-n], a$y[-n], pm)
  points(pb$x, pb$y, pch = 21, bg = "#5b8db8", col = "#2c3e50", cex = 1.1)
  s1 <- trans3d(a$x1[n], a$x2[n], c(a$yhat[n], a$y[n]), pm)
  lines(s1$x, s1$y, col = colp, lwd = 1.6)
  points(s1$x[2], s1$y[2], pch = 23, bg = colp, cex = 2.4, lwd = 1.5)
  title(main = "Espacio de las variables (x₁, x₂, y)", cex.main = 1, font.main = 2, adj = 0)
  mtext(if (mostrar_sin) "rojo: plano con el punto · gris: sin el punto" else "rojo: plano ajustado con el punto",
        side = 3, line = 0.1, adj = 0, cex = 0.8, col = "#555555")
}

g_observaciones <- function(a, theta, phi, k_e = 1) {
  n <- a$n
  x1c <- a$x1 - mean(a$x1); x2c <- a$x2 - mean(a$x2); yc <- a$y - mean(a$y)
  u1 <- unit(x1c); u2 <- unit(x2c - sum(x2c * u1) * u1)
  ne <- sqrt(sum(a$e^2))
  u3 <- if (ne > 1e-9) sign(a$e[n] + 1e-12) * a$e / ne else rep(0, n)
  B  <- rbind(u1, u2, u3)
  Y  <- as.vector(B %*% yc); Yh <- as.vector(B %*% (a$yhat - mean(a$y)))
  E  <- as.vector(B %*% a$e)
  E[3] <- k_e * E[3]; Y[3] <- k_e * Y[3]     # ampliación opcional del eje de e
  X1 <- as.vector(B %*% x1c); X2 <- as.vector(B %*% x2c)
  L  <- max(10, sqrt(sum(Y^2)) * 1.15)
  sc <- 0.75 * L / max(sqrt(sum(X1^2)), sqrt(sum(X2^2)))

  par(mar = c(1, 1, 2.5, 1))
  q <- c(-0.3 * L, L)
  pm <- persp(q, q, matrix(0, 2, 2), zlim = c(-0.4 * L, 0.4 * L), theta = theta, phi = phi,
              col = adjustcolor("#5b8db8", 0.18), border = adjustcolor("#5b8db8", 0.6),
              xlab = "", ylab = "", zlab = "", expand = 0.85)
  P <- function(v) trans3d(v[1], v[2], v[3], pm)
  fl <- function(d, h, ...) { a1 <- P(d); b1 <- P(h); arrows(a1$x, a1$y, b1$x, b1$y, length = 0.08, angle = 20, ...) }
  et <- function(v, lab, ...) { pp <- P(v); text(pp$x, pp$y, lab, ...) }
  o <- c(0, 0, 0)
  lr <- trans3d(c(0, 0), c(0, 0), c(-0.4, 0.4) * L, pm)
  lines(lr$x, lr$y, lty = 2, lwd = 1.5, col = "#27ae60")
  fl(o, X1 * sc, col = "#2c3e50", lwd = 1.8)
  fl(o, X2 * sc, col = "#2c3e50", lwd = 1.8)
  fl(o, Y,  col = "black",   lwd = 2.4)
  fl(o, Yh, col = "#e74c3c", lwd = 2.4)
  if (ne > 0.01 * L) {
    fl(Yh, Y, col = "#27ae60", lwd = 2.4)
    fl(o, E, col = adjustcolor("#27ae60", 0.45), lwd = 1.3)
    t <- 0.05 * L; d1 <- unit(-Yh); d2 <- c(0, 0, sign(E[3]))
    ar <- rbind(Yh + t * d1, Yh + t * (d1 + d2), Yh + t * d2)
    pa <- trans3d(ar[, 1], ar[, 2], ar[, 3], pm)
    lines(pa$x, pa$y, col = "#333333", lwd = 1.1)
  }
  et(X1 * sc * 1.1, "x₁", col = "#2c3e50")
  et(X2 * sc * 1.1, "x₂", col = "#2c3e50")
  et(Y + c(0, 0, 0.06 * L * sign(E[3] + 1e-12)), "y", font = 2, cex = 1.2)
  et(Yh * 1.05 - c(0, 0, 0.07 * L), "ŷ = Hy", col = "#e74c3c")
  if (ne > 0.1 * L) et(Yh + 0.5 * E + c(0.1 * L, 0, 0), "e = (I − H)y", col = "#1e8449")
  et(c(0.85 * L, -0.3 * L, 0), "col(X)", col = "#2e6da4", font = 2)
  et(c(0, 0, 0.45 * L), "col(I − H)", col = "#1e8449", cex = 0.9)
  title(main = "Espacio de las observaciones (ℝⁿ)", cex.main = 1, font.main = 2, adj = 0)
  mtext(if (k_e == 1) "vista exacta en span{x₁, x₂, y} (centradas) · ŷ ⊥ e"
        else sprintf("span{x₁, x₂, y} (centradas) · eje de e ampliado ×%g · ŷ ⊥ e", k_e),
        side = 3, line = 0.1, adj = 0, cex = 0.8, col = "#555555")
}

g_cook <- function(a, activos, log_y) {
  n <- a$n; U <- a$U
  cols <- sapply(a$D, col_cook, U = U, activos = activos)
  yl <- if (log_y) c(1e-3, max(2, max(a$D) * 1.5)) else c(0, max(1.2, max(a$D) * 1.1))
  par(mar = c(4, 4.5, 2.5, 1))
  Dp <- if (log_y) pmax(a$D, 1e-3) else a$D
  plot(seq_len(n), Dp, type = "h", lwd = 2, col = cols, ylim = yl, log = if (log_y) "y" else "",
       bty = "l", las = 1, xlab = "observación i", ylab = expression(D[i]))
  points(seq_len(n), Dp, pch = 19, col = cols)
  points(n, Dp[n], pch = 23, bg = cols[n], cex = 1.8)
  lt <- c(F = 1, uno = 2, cuatro = 3)
  for (k in activos) abline(h = U[k], lty = lt[k], col = "#7f8c8d", lwd = 1.3)
  leg <- c(F = sprintf("C1: F0.5(p, n−p) = %.2f", U["F"]), uno = "C2: 1",
           cuatro = sprintf("C3: 4/(n−p−1) = %.2f", U["cuatro"]))
  if (length(activos)) legend("topleft", leg[activos], lty = lt[activos], col = "#7f8c8d",
                              bty = "n", cex = 0.85)
  title(main = sprintf("Distancia de Cook de cada observación (la %d es la móvil)", n),
        adj = 0, cex.main = 1)
}

g_leverage <- function(a, activos) {
  n <- a$n; U <- a$U
  hmax <- max(0.95, max(a$h) * 1.05)
  rmax <- max(4, abs(a$r) * 1.1)
  par(mar = c(4, 4.5, 2.5, 1))
  plot(a$h, a$r, xlim = c(0, hmax), ylim = c(-rmax, rmax), pch = 19, col = "#5b8db8",
       bty = "l", las = 1, xlab = expression("Leverage  " * h[ii]),
       ylab = expression("Residuo estudentizado  " * r[i] / (hat(sigma) * sqrt(1 - h[ii]))))
  abline(h = c(-2, 2), col = "grey80", lty = 2)
  abline(v = 2 * p / n, col = "grey80", lty = 2)
  hh <- seq(0.001, hmax, length.out = 300)
  lt <- c(F = 1, uno = 2, cuatro = 3)
  for (k in activos) {
    rr <- sqrt(p * U[k] * (1 - hh) / hh)   # curva de nivel D = U
    lines(hh, rr, lty = lt[k], col = "#c0392b"); lines(hh, -rr, lty = lt[k], col = "#c0392b")
  }
  colp <- col_cook(a$D[n], U, activos)
  points(a$h[n], a$r[n], pch = 23, bg = colp, cex = 2.2, lwd = 1.5)
  text(2 * p / n, -rmax, " 2p/n", adj = c(0, 0), cex = 0.8, col = "grey40")
  title(main = "Palanca vs. atipicidad (curvas: Dᵢ = umbral)", adj = 0, cex.main = 1)
}

g_formula <- function(a, hist) {
  n <- a$n; hn <- a$h[n]; rn <- a$r[n]
  layout(matrix(1:2, 1), widths = c(1, 1.25))
  par(mar = c(0.5, 0.5, 2.5, 0.5))
  plot.new()
  title(main = "Distancia de Cook (influencia de la obs. i)", adj = 0, cex.main = 1)
  cex <- 1.5; x <- 0.02; y <- 0.66; lims <- list()
  partes <- list(quote(D[i] == phantom(0)),
                 quote(bgroup("(", frac(r[i], hat(sigma) * sqrt(1 - h[ii])), ")")^2),
                 quote(phantom(0) * frac(h[ii], 1 - h[ii])),
                 quote(phantom(0) * frac(1, p)))
  for (k in seq_along(partes)) {
    e <- as.expression(partes[[k]]); w <- strwidth(e, cex = cex)
    text(x, y, e, adj = c(0, 0.5), cex = cex); lims[[k]] <- c(x, x + w); x <- x + w
  }
  yb <- y - 0.17
  m1 <- llave(lims[[2]][1] + 0.01, lims[[2]][2] - 0.005, yb, 0.05, col = "#c0392b", lwd = 1.4)
  text(m1[1], m1[2] - 0.1, sprintf("Atipicidad\n%.2f", rn^2), col = "#c0392b", cex = 0.85)
  m2 <- llave(lims[[3]][1] + 0.02, lims[[3]][2] - 0.005, yb, 0.05, col = "#2e6da4", lwd = 1.4)
  text(m2[1], m2[2] - 0.1, sprintf("leverage\n%.2f", hn / (1 - hn)), col = "#2e6da4", cex = 0.85)
  text(0.02, 0.04, sprintf("Punto móvil:  D = %.2f × %.2f × 1/%d = %.3f",
                           rn^2, hn / (1 - hn), p, a$D[n]), adj = 0, cex = 0.95, font = 2)

  par(mar = c(3, 5.5, 2.5, 1))
  if (length(hist) < 2) {
    plot.new(); title(main = "Historial de Dᵢ del punto móvil", adj = 0, cex.main = 1)
    text(0.5, 0.5, "Mueve el punto para ver el historial", col = "grey50")
  } else {
    hp <- pmax(hist, 1e-3)
    plot(seq_along(hp), hp, type = "l", lwd = 2, col = "#e74c3c", log = "y",
         ylim = c(1e-3, max(2, max(hp) * 1.5)), bty = "l", las = 1, xlab = "", ylab = "", cex.axis = 0.85)
    mtext("D (escala log)", side = 2, line = 4, cex = 0.85)
    abline(h = a$U, lty = 1:3, col = "#7f8c8d")
    points(length(hp), hp[length(hp)], pch = 21, bg = "#e74c3c", cex = 1.4)
    title(main = "Historial de Dᵢ del punto móvil", adj = 0, cex.main = 1)
  }
}

tabla_resumen <- function(a, activos) {
  n <- a$n; y <- a$y
  SCT <- sum((y - mean(y))^2); SCR <- sum((a$yhat - mean(y))^2); SCE <- sum(a$e^2)
  hn <- a$h[n]; rn <- a$r[n]; Dn <- a$D[n]; U <- a$U
  ok <- function(c) if (c) "<b style='color:#c0392b'>SÍ</b>" else "no"
  crit <- c(
    F      = sprintf("<tr><td>C1: D &gt; F<sub>0.5</sub>(p, n−p) = %.3f</td><td>%s</td></tr>", U["F"], ok(Dn > U["F"])),
    uno    = sprintf("<tr><td>C2: D &gt; 1</td><td>%s</td></tr>", ok(Dn > 1)),
    cuatro = sprintf("<tr><td>C3: D &gt; 4/(n−p−1) = %.3f</td><td>%s</td></tr>", U["cuatro"], ok(Dn > U["cuatro"])))
  infl <- if (length(activos) && any(c(F = Dn > U[["F"]], uno = Dn > 1, cuatro = Dn > U[["cuatro"]])[activos]))
    "<span style='color:#c0392b;font-weight:bold'>INFLUYENTE según al menos un criterio activo</span>"
  else "<span style='color:#27ae60;font-weight:bold'>No influyente según los criterios activos</span>"
  HTML(paste0(
    "<table class='table table-condensed' style='font-size:13px'>",
    "<tr><th colspan=3>Análisis de varianza (con el punto móvil)</th></tr>",
    "<tr><td>Fuente</td><td>SC</td><td>gl</td></tr>",
    sprintf("<tr><td>Regresión (SCR = Σ(ŷ<sub>i</sub> − ȳ)²)</td><td>%.2f</td><td>%d</td></tr>", SCR, p - 1),
    sprintf("<tr><td>Error / residual (SCE = Σe<sub>i</sub>²)</td><td>%.2f</td><td>%d</td></tr>", SCE, n - p),
    sprintf("<tr><td>Total (SCT = Σ(y<sub>i</sub> − ȳ)²)</td><td>%.2f</td><td>%d</td></tr>", SCT, n - 1),
    sprintf("<tr><td>R² = SCR/SCT</td><td>%.3f</td><td></td></tr>", SCR / SCT),
    sprintf("<tr><td>σ̂² = SCE/(n−p)</td><td>%.3f</td><td></td></tr>", SCE / (n - p)),
    "</table>",
    "<table class='table table-condensed' style='font-size:13px'>",
    sprintf("<tr><th colspan=2>Punto móvil (obs. %d), n = %d, p = %d</th></tr>", n, n, p),
    sprintf("<tr><td>Leverage h<sub>ii</sub> (corte 2p/n = %.2f)</td><td>%.3f %s</td></tr>",
            2 * p / n, hn, if (hn > 2 * p / n) "<b>ALTO</b>" else ""),
    sprintf("<tr><td>Residuo estudentizado r<sub>i</sub>/(σ̂√(1−h<sub>ii</sub>))</td><td>%+.2f %s</td></tr>",
            rn, if (abs(rn) > 2) "<b>ATÍPICO</b>" else ""),
    sprintf("<tr><td>Distancia de Cook D<sub>i</sub></td><td><b>%.3f</b></td></tr>", Dn),
    paste(crit[activos], collapse = ""),
    sprintf("<tr><td colspan=2>%s</td></tr>", infl),
    "</table>",
    "<table class='table table-condensed' style='font-size:13px'>",
    "<tr><th>Coeficientes</th><th>b₀</th><th>b₁</th><th>b₂</th></tr>",
    sprintf("<tr><td>Sin el punto</td><td>%+.3f</td><td>%+.3f</td><td>%+.3f</td></tr>",
            a$beta_sin[1], a$beta_sin[2], a$beta_sin[3]),
    sprintf("<tr><td>Con el punto</td><td>%+.3f</td><td>%+.3f</td><td>%+.3f</td></tr>",
            a$beta[1], a$beta[2], a$beta[3]),
    sprintf("<tr><td>Verdaderos</td><td>%+.3f</td><td>%+.3f</td><td>%+.3f</td></tr>", b0r, b1r, b2r),
    "</table>"))
}

# ==============================================================================
# Interfaz
# ==============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-size: 13px; }
    .well { padding: 12px; }
    .fase { background:#fdf2e9; border-left:4px solid #e67e22; padding:8px 12px; margin-bottom:10px; }
    h4 { margin-top: 14px; }
    /* Botón de reproducción del slider: grande y centrado */
    .slider-animate-container { text-align: center !important; margin-top: 14px; }
    .slider-animate-button { display: block; opacity: 1 !important; text-decoration: none !important; }
    .btn-anim { width: 100%; font-size: 17px; font-weight: bold; padding: 10px 0; }
  "))),
  titlePanel("Atípicos, alta palanca y observaciones influyentes (dos variables explicativas)"),
  sidebarLayout(
    sidebarPanel(width = 3,
      h4("Punto móvil"),
      sliderInput("x1n", "x₁ del punto", min = -5, max = 30, value = 5, step = 0.25),
      sliderInput("x2n", "x₂ del punto", min = -5, max = 30, value = 5, step = 0.25),
      sliderInput("dn", "Desviación vertical respecto a la tendencia verdadera",
                  min = -20, max = 20, value = 0, step = 0.25),
      helpText("x₁ y x₂ controlan la palanca; la desviación vertical, la atipicidad."),
      h4("Escenarios rápidos"),
      fluidRow(
        column(6, actionButton("esc_centro", "Centro", width = "100%")),
        column(6, actionButton("esc_lev", "Solo palanca", width = "100%"))),
      div(style = "height:6px"),
      fluidRow(
        column(6, actionButton("esc_atip", "Solo atípico", width = "100%")),
        column(6, actionButton("esc_infl", "Influyente", width = "100%"))),
      h4("Recorrido automático"),
      checkboxInput("auto", "Usar el recorrido del GIF (se activa solo al reproducir o mover el paso)", FALSE),
      sliderInput("paso", "Paso", min = 1, max = nrow(recorrido), value = 1, step = 1,
                  animate = animationOptions(
                    interval = 150, loop = TRUE,
                    playButton  = tags$span(class = "btn btn-success btn-lg btn-anim",
                                            icon("play"), " Reproducir recorrido"),
                    pauseButton = tags$span(class = "btn btn-warning btn-lg btn-anim",
                                            icon("pause"), " Pausar"))),
      h4("Datos base"),
      sliderInput("n0", "Número de observaciones fijas", min = 8, max = 80, value = 19, step = 1),
      sliderInput("sd", "Desviación estándar del ruido σ", min = 0.2, max = 4, value = 1, step = 0.1),
      numericInput("semilla", "Semilla", value = 7, min = 1, step = 1),
      h4("Criterios de influencia"),
      checkboxGroupInput("crit", NULL,
        choices = c("C1: D > F0.5(p, n−p)" = "F", "C2: D > 1" = "uno", "C3: D > 4/(n−p−1)" = "cuatro"),
        selected = c("F", "uno", "cuatro")),
      checkboxInput("logy", "Escala logarítmica en el gráfico de Cook", TRUE),
      h4("Vista 3D"),
      sliderInput("theta", "Giro horizontal", min = -180, max = 180, value = 30, step = 5),
      sliderInput("phi", "Inclinación", min = 0, max = 60, value = 20, step = 2),
      checkboxInput("sin", "Mostrar plano sin el punto", TRUE),
      sliderInput("k_e", "Ampliar el eje de e = (I − H)y en el espacio de las observaciones",
                  min = 1, max = 5, value = 1, step = 0.5),
      helpText("Con 1 la vista conserva longitudes; con más de 1 se exagera ‖e‖ para verlo mejor."),
      actionButton("borrar", "Borrar historial de D", width = "100%")
    ),
    mainPanel(width = 9,
      uiOutput("fase"),
      fluidRow(
        column(6, plotOutput("p_var", height = "400px")),
        column(6, plotOutput("p_obs", height = "400px"))),
      fluidRow(
        column(6, plotOutput("p_cook", height = "300px")),
        column(6, plotOutput("p_lev", height = "300px"))),
      fluidRow(
        column(7, plotOutput("p_formula", height = "260px")),
        column(5, uiOutput("tabla")))
    )
  )
)

# ==============================================================================
# Servidor
# ==============================================================================
server <- function(input, output, session) {

  base <- reactive({
    set.seed(input$semilla)
    n0 <- input$n0
    x1 <- runif(n0, 0, 10); x2 <- runif(n0, 0, 10)
    list(x1 = x1, x2 = x2, y = pred(x1, x2) + rnorm(n0, 0, input$sd))
  })

  punto <- reactive({
    if (isTRUE(input$auto)) {
      r <- recorrido[input$paso, ]
      list(x1 = r$x1, x2 = r$x2, d = r$d, fase = r$fase)
    } else list(x1 = input$x1n, x2 = input$x2n, d = input$dn, fase = NA)
  })

  ajuste <- reactive({
    pt <- punto()
    ajustar(base(), pt$x1, pt$x2, pred(pt$x1, pt$x2) + pt$d)
  })

  activos <- reactive(if (is.null(input$crit)) character(0) else input$crit)

  # Escenarios rápidos
  fijar <- function(a, b, d) {
    updateCheckboxInput(session, "auto", value = FALSE)
    updateSliderInput(session, "x1n", value = a)
    updateSliderInput(session, "x2n", value = b)
    updateSliderInput(session, "dn", value = d)
  }
  observeEvent(input$esc_centro, fijar(5, 5, 0))
  observeEvent(input$esc_lev,    fijar(24, 22, 0))
  observeEvent(input$esc_atip,   fijar(5, 5, 8))
  observeEvent(input$esc_infl,   fijar(24, 22, -9))

  # Al reproducir o mover el paso se activa el recorrido automático;
  # al mover a mano el punto se vuelve al modo manual.
  observeEvent(input$paso, updateCheckboxInput(session, "auto", value = TRUE), ignoreInit = TRUE)
  observeEvent(list(input$x1n, input$x2n, input$dn),
               updateCheckboxInput(session, "auto", value = FALSE), ignoreInit = TRUE)

  # Historial de la distancia de Cook del punto móvil
  hist <- reactiveVal(numeric(0))
  observeEvent(ajuste(), {
    a <- ajuste()
    hist(tail(c(hist(), a$D[a$n]), 200))
  })
  observeEvent(input$borrar, hist(numeric(0)))

  output$fase <- renderUI({
    pt <- punto()
    txt <- if (!is.na(pt$fase)) textos_fase[pt$fase] else
      "Modo manual · mueve el punto con los controles o usa un escenario rápido."
    div(class = "fase", txt)
  })

  colp <- reactive({ a <- ajuste(); col_cook(a$D[a$n], a$U, activos()) })

  output$p_var     <- renderPlot(g_variables(ajuste(), input$theta, input$phi, input$sin, colp()))
  output$p_obs     <- renderPlot(g_observaciones(ajuste(), input$theta, input$phi, input$k_e))
  output$p_cook    <- renderPlot(g_cook(ajuste(), activos(), input$logy))
  output$p_lev     <- renderPlot(g_leverage(ajuste(), activos()))
  output$p_formula <- renderPlot(g_formula(ajuste(), hist()))
  output$tabla     <- renderUI(tabla_resumen(ajuste(), activos()))
}

shinyApp(ui, server)

# Regenerates the README figures. Run from the package root:
#   Rscript tools/make-figures.R
# Needs ggplot2, patchwork and Rmpfr. Takes about two minutes, most of it in
# the Rmpfr reference run of the divergence figure.

library(chaosrcpp)
library(ggplot2)
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
fig <- function(name) file.path("man/figures", name)

# 1. The butterfly effect of rounding error -----------------------------------
message("divergence (mpfr reference, this is the slow one) ...")
d <- precision_divergence(t_max = 60, reference = "mpfr", mpfr_bits = 128)
print(d)
p <- plot_divergence(d)
ggsave(fig("divergence.png"), p, width = 9, height = 6.5, dpi = 130, bg = "white")

# 2. Bifurcation diagram with the Lyapunov exponent ---------------------------
message("bifurcation ...")
r <- seq(2.8, 4, length.out = 2500)
bif <- logistic_bifurcation(r, n_transient = 600, n_keep = 250)
ly <- logistic_lyapunov(r, n_transient = 600, n = 3000)
p <- plot_bifurcation(bif, ly, point_alpha = 0.06)
ggsave(fig("bifurcation-lyapunov.png"), p, width = 9, height = 6.5, dpi = 130, bg = "white")

# 3. Lorenz attractor and its recurrence plot ---------------------------------
message("recurrence plot ...")
tr <- lorenz(n = 4000, dt = 0.02)
png(fig("lorenz-recurrence.png"), width = 1500, height = 720, res = 130)
layout(matrix(1:2, 1), widths = c(1.15, 1))
cols <- hcl.colors(nrow(tr), "Viridis")
plot_attractor(tr, vars = c("x", "z"), col = cols, lwd = 0.8)
plot_recurrence(tr, rr = 0.04, main = "Recurrence plot (RR = 4%)")
dev.off()

# 4. Chirikov standard map ----------------------------------------------------
message("standard map ...")
sm <- standard_map(K = 0.971635, n_orbits = 90, n = 1500, seed = 11)
p <- plot_standard_map(sm, size = 0.12)
ggsave(fig("standard-map.png"), p, width = 7.5, height = 7.8, dpi = 130, bg = "white")

# 5. Henon and Ikeda attractors -----------------------------------------------
message("henon + ikeda ...")
png(fig("henon-ikeda.png"), width = 1500, height = 640, res = 130)
layout(matrix(1:2, 1))
plot(henon(n = 300000), col = adjustcolor("#1f2a44", 0.35))
plot(ikeda(n = 400000), col = adjustcolor("#7b1e3a", 0.25))
dev.off()

# 6. Koopman view: Lorenz DMD spectrum, logistic sweep, HAVOK forcing ----------
message("dmd + havok ...")
dl <- dmd(lorenz(n = 8000, dt = 0.01), delays = 60, rank = 40)
rs <- seq(3.4, 4, length.out = 241)
ly <- logistic_lyapunov(rs)$lambda
e1 <- vapply(rs, function(r) {
  x <- logistic_orbit(r, x0 = 0.2, n = 1500)[-(1:500)]
  dmd_features(x, delays = 16, center = TRUE)$dmd_one_step_nrmse
}, numeric(1))
x0 <- unlist(utils::tail(lorenz(n = 5000, dt = 0.002), 1)[c("x", "y", "z")])  # start on the attractor
h <- havok(lorenz(x0 = x0, n = 60000, dt = 0.002), delays = 50, rank = 11)
png(fig("dmd-koopman.png"), width = 1800, height = 620, res = 130)
layout(matrix(1:3, 1), widths = c(1, 1.25, 1.6))
par(mar = c(4, 4, 2.5, 1))
th <- seq(0, 2 * pi, length.out = 361)
plot(cos(th), sin(th), type = "l", col = "grey60", asp = 1, xlab = "Re", ylab = "Im",
     main = "Lorenz: Hankel-DMD eigenvalues")
points(Re(dl$eigenvalues), Im(dl$eigenvalues), pch = 16,
       cex = 0.6 + 3 * sqrt(dl$energy / max(dl$energy)), col = adjustcolor("#1f2a44", 0.7))
par(mar = c(4, 4, 2.5, 4))
plot(rs, ly, type = "l", xlab = "r", ylab = "Lyapunov exponent",
     main = sprintf("Logistic map (Spearman %.2f)", cor(ly, e1, method = "spearman")))
abline(h = 0, col = "grey75")
par(new = TRUE)
plot(rs, e1, type = "l", col = "#c0392b", axes = FALSE, xlab = "", ylab = "")
axis(4, col.axis = "#c0392b")
mtext("DMD one-step error", 4, line = 2.5, col = "#c0392b", cex = 0.8)
par(mar = c(4, 4, 2.5, 1))
keep <- h$t < 60
act <- abs(h$forcing) > 2 * sd(h$forcing)
plot(h$t[keep], h$x[keep], type = "l", col = "grey45", xlab = "t", ylab = "x",
     main = "HAVOK: active forcing (red) on Lorenz x")
points(h$t[keep & act], h$x[keep & act], pch = 16, cex = 0.35, col = "#c0392b")
dev.off()

# 7. Benchmark ----------------------------------------------------------------
message("benchmark ...")
b <- benchmark_lorenz(n = 1e6)
print(b)
saveRDS(b, "tools/benchmark.rds")
message("done")

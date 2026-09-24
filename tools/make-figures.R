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

# 6. Benchmark ----------------------------------------------------------------
message("benchmark ...")
b <- benchmark_lorenz(n = 1e6)
print(b)
saveRDS(b, "tools/benchmark.rds")
message("done")

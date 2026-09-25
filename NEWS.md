# chaosrcpp 0.2.0

* New `dmd()`: exact Dynamic Mode Decomposition with optional delay
  coordinates (Hankel DMD), with `print()`, `plot()`, `fitted()` and
  `predict()` methods. Accepts vectors, state matrices, `chaos_trajectory`
  and `chaos_map` objects.
* New `dmd_features()`: nine DMD summary measures as a one-row data frame,
  the counterpart of `rqa_features()`.
* New `havok()`: Hankel Alternative View Of Koopman (Brunton et al. 2017),
  a linear model of eigen-time-delay coordinates with intermittent forcing.
* New vignette "Linear lenses on chaos: DMD, Hankel DMD and HAVOK".

# chaosrcpp 0.1.0

* First release.

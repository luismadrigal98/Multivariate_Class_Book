## ============================================================================
##  MS_LJMR :: 00_utils.R
##  Shared utilities for the Multivariate Data Analysis course code base
##  (Soberon / Jimenez / Cobos originals, revised & integrated with the ML
##   course by Luis J. Madrigal-Roca).
##
##  Design goals
##  ------------
##  * No hard-coded setwd() and no x11(): scripts run on Windows, macOS, Linux
##    and in batch (Rscript) sessions.
##  * Every data set is obtained through get_data(): it loads the real course
##    file when it is present on disk, and otherwise falls back to a
##    reproducible, seeded simulation with the SAME column structure, so the
##    scripts always run and figures always render. Drop the real CSV/RDS into
##    MS_LJMR/data/ (or set options(msljmr.data_dir=...)) and the loaders will
##    pick it up automatically -- no code change needed.
##  * open_dev() replaces x11(): interactive on screen, PNG in batch mode.
##
##  Author: Luis J. Madrigal-Roca   |   License: GPL-2 (matches course code)
## ============================================================================

## ---- paths -----------------------------------------------------------------
.msljmr_root <- function() {
  ## directory that holds this file's parent (…/MS_LJMR)
  getOption("msljmr.root", default = normalizePath(file.path(
    dirname(dirname(sys.frame(1)$ofile %||% ".")), "."), mustWork = FALSE))
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

data_dir <- function() {
  getOption("msljmr.data_dir",
            default = file.path(getOption("msljmr.root", "."), "data"))
}

fig_dir <- function() {
  d <- getOption("msljmr.fig_dir",
                 default = file.path(getOption("msljmr.root", "."), "figs"))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

## ---- graphics device -------------------------------------------------------
##  Where a figure goes is decided by one option, `msljmr.device`:
##
##    "auto"   (default) screen when the session is interactive, PNG under
##             Rscript / R CMD BATCH
##    "screen" always draw on the CURRENT device. Inside RStudio that is the
##             Plot pane, so figures pile up in the plot history and you can
##             page back and forth through them while teaching
##    "png"    always write figs/<name>.png
##
##  Set it for a whole session with options(msljmr.device = "png"), or per
##  figure with with_fig(..., device = "png").
##
##  Screen mode deliberately opens nothing: the plotting call itself brings up
##  the default device, which is RStudioGD inside RStudio and X11/quartz/windows
##  elsewhere. An earlier version called dev.new(noRStudioGD = TRUE), which is
##  exactly the flag that throws a separate window in front of the class.
##
##  Two more options, both off by default:
##    msljmr.save_figs  TRUE  -> in screen mode ALSO copy each figure to figs/
##    msljmr.pause      TRUE  -> wait for <Enter> between figures, for live demos

in_rstudio <- function()
  identical(Sys.getenv("RSTUDIO"), "1") || identical(.Platform$GUI, "RStudio")

.fig_mode <- function(device = NULL) {
  mode <- device %||% getOption("msljmr.device", "auto")
  if (identical(mode, "auto")) mode <- if (interactive()) "screen" else "png"
  if (!mode %in% c("screen", "png"))
    stop("msljmr.device must be \"auto\", \"screen\" or \"png\"", call. = FALSE)
  mode
}

##  open_dev("name", w, h): starts a figure and returns the mode it used, so the
##  caller knows whether it owns a device that must be closed.
open_dev <- function(name = "plot", width = 7, height = 6, res = 150,
                     device = NULL) {
  mode <- .fig_mode(device)
  if (mode == "png") {
    grDevices::png(file.path(fig_dir(), paste0(name, ".png")),
                   width = width, height = height, units = "in", res = res)
  } else if (!in_rstudio() && grDevices::dev.cur() == 1L) {
    ## plain console with no device yet: give it a window of the right shape
    grDevices::dev.new(width = width, height = height)
  }
  invisible(mode)
}

## convenience wrapper: with_fig("scree", { plot(...) })
##
##  with_fig() also saves and restores par(). That matters in screen mode, where
##  every figure shares one device: without it a par(mfrow = c(1, 2)) set for one
##  figure would silently split the next one too. Note that an on.exit() written
##  *inside* the { } block does NOT do this -- the block is a promise, so its
##  on.exit registers against the caller's frame and never fires. Restore par
##  here, once, rather than in each script.
with_fig <- function(name, expr, width = 7, height = 6, res = 150,
                     device = NULL) {
  mode <- open_dev(name, width, height, res, device)
  old  <- if (mode == "screen" && grDevices::dev.cur() != 1L)
            graphics::par(no.readonly = TRUE) else NULL

  on.exit({
    if (mode == "png") grDevices::dev.off()
    else if (!is.null(old)) try(graphics::par(old), silent = TRUE)
  }, add = TRUE)

  force(expr)

  if (mode == "screen") {
    if (isTRUE(getOption("msljmr.save_figs", FALSE))) {
      grDevices::dev.copy(grDevices::png,
                          filename = file.path(fig_dir(), paste0(name, ".png")),
                          width = width, height = height, units = "in", res = res)
      grDevices::dev.off()
    }
    if (isTRUE(getOption("msljmr.pause", FALSE)) && interactive())
      invisible(readline(paste0("[", name, "]  <Enter> for the next figure ")))
  }
  invisible(NULL)
}

## ---- package helper --------------------------------------------------------
##  need("vegan","MASS"): load packages, telling the user how to install any
##  that are missing instead of crashing mid-script.
need <- function(...) {
  pkgs <- c(...)
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    message("Missing packages: ", paste(miss, collapse = ", "),
            "\n  install.packages(c(",
            paste(sprintf('\"%s\"', miss), collapse = ", "), "))")
  }
  invisible(lapply(setdiff(pkgs, miss), function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))))
}

## ---- optional-dependency guard ---------------------------------------------
##  Several classroom sections need a package that may not be installed, is no
##  longer on CRAN (maptools, BiplotGUI), needs a live network (WorldClim), or
##  opens an interactive window (rgl). Those sections are wrapped in
##
##      if (has_pkg("rgl")) { ... } else skip_note("interactive 3-D", "rgl")
##
##  so the demo runs in class on a machine that has the package, and merely
##  prints a note everywhere else instead of aborting halfway through a session.
has_pkg <- function(...)
  all(vapply(c(...), requireNamespace, logical(1), quietly = TRUE))

skip_note <- function(what, pkgs) {
  message("  [skipped] ", what, " -- needs ", paste(pkgs, collapse = ", "),
          ":  install.packages(c(",
          paste(sprintf('"%s"', pkgs), collapse = ", "), "))")
  invisible(FALSE)
}

## ---- evaluation helper -----------------------------------------------------
##  auc(p, y): area under the ROC curve from predicted scores `p` and binary
##  outcomes `y`, computed from the rank-sum identity AUC = (W - n1(n1+1)/2) /
##  (n1 n0) so no extra package is needed. Shared by the two ENM sessions.
auc <- function(p, y) {
  r  <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ============================================================================
##  get_data(): the real-file-or-fallback loader
## ============================================================================
##  get_data("iris")                  -> Anderson's iris (built-in)
##  get_data("butterflies")           -> Quintana Roo wing-pattern counts
##  get_data("biodiv")                -> biodiversity-of-countries table
##  get_data("taxon")                 -> Crawley 4-taxa morphometrics
##  get_data("crawley")               -> speciesCrawley3 mixed-type table
##  get_data("europe")                -> 21 European cities (lon/lat)
##  get_data("countries_live")        -> ordinal liveability rankings
##  get_data("leukemia")              -> gene-expression subset (genes x samples)
##  get_data("limenitis")             -> grayscale wing image matrix
##  get_data("doubs")                 -> list(fish, env, xy) river gradient
## ----------------------------------------------------------------------------
get_data <- function(name, seed = 1998, quiet = FALSE) {
  name <- tolower(name)
  fmap <- list(
    butterflies    = "ButterfliesQRoo2.csv",
    biodiv         = "BiodivCountries.csv",
    taxon          = "taxon.csv",
    crawley        = "speciesCrawley3.csv",
    europe         = "CitiesEurope.csv",
    countries_live = "CountriesToLive.csv",
    leukemia       = "leukemiaExpressionSubset.rds",
    limenitis      = "Limenitis_archippus.csv",
    neotoma        = "NeotomaMorphoEnvir.csv",
    hanta          = "hanta_virtual.csv",
    biodiv_pc      = "BiodiversityCountriesPCValues.csv",
    biodiv_gv      = "BiodiversityCountriesBiGv.csv",
    ssa_factor     = "BiodiversityCountriesSSAFactanal2.csv",
    insatisf       = "Insatisf2.csv"
  )

  ## Datasets whose real file is a headerless numeric grid: no header row and no
  ## row names, every line pure data. read.csv()'s default header = TRUE would
  ## silently promote the first row of pixels to column names and hand back a
  ## data frame, which then fails in anything requiring a true matrix (norm(),
  ## svd() on a subset, image()). Read these with header = FALSE and coerce.
  matrix_files <- c("limenitis")
  ## 1) built-ins
  if (name == "iris")   { utils::data("iris",   envir = environment()); return(iris) }
  if (name == "mtcars") { utils::data("mtcars", envir = environment()); return(mtcars) }

  ## 2) real course file on disk?
  f <- fmap[[name]]
  if (!is.null(f)) {
    path <- file.path(data_dir(), f)
    if (file.exists(path)) {
      if (!quiet) message("get_data('", name, "'): using real file ", path)
      if (grepl("\\.rds$", f)) return(readRDS(path))
      if (name %in% matrix_files)
        return(as.matrix(utils::read.csv(path, header = FALSE)))
      return(utils::read.csv(path, stringsAsFactors = TRUE))
    }
  }

  ## 2a) GENUINE in-package datasets the class actually used: prefer the real
  ##     data when the package is installed on the user's machine.
  if (name == "doubs" && requireNamespace("ade4", quietly = TRUE)) {
    if (!quiet) message("get_data('doubs'): using the REAL ade4::doubs dataset")
    e <- new.env(); utils::data("doubs", package = "ade4", envir = e)
    return(list(fish = e$doubs$fish, env = e$doubs$env, xy = e$doubs$xy))
  }

  ## 2b) frozen DEMO file (data/demo/) -- shared byte-for-byte with the Python
  ##     verifier, so deterministic results (eigenvalues, inertias, loadings)
  ##     reproduce exactly. Generated by verify/freeze_and_reference.py.
  demo_dir  <- file.path(data_dir(), "demo")
  demo_read <- function(fn) utils::read.csv(file.path(demo_dir, fn), stringsAsFactors = TRUE)
  if (name == "doubs" && file.exists(file.path(demo_dir, "doubs_fish.csv"))) {
    if (!quiet) message("get_data('doubs'): using frozen demo data")
    return(list(fish = demo_read("doubs_fish.csv"), env = demo_read("doubs_env.csv")))
  }
  if (name == "limenitis" && file.exists(file.path(demo_dir, "limenitis.csv"))) {
    if (!quiet) message("get_data('limenitis'): using frozen demo data")
    return(as.matrix(utils::read.csv(file.path(demo_dir, "limenitis.csv"), header = FALSE)))
  }
  if (name == "leukemia" && file.exists(file.path(demo_dir, "leukemia.csv"))) {
    if (!quiet) message("get_data('leukemia'): using frozen demo data")
    m <- utils::read.csv(file.path(demo_dir, "leukemia.csv"), row.names = 1, check.names = FALSE)
    return(as.matrix(m))
  }
  demo1 <- c(taxon = "taxon.csv", crawley = "crawley.csv", hanta = "hanta.csv",
             pam = "pam.csv")
  if (!is.null(demo1[[name]]) && file.exists(file.path(demo_dir, demo1[[name]]))) {
    if (!quiet) message("get_data('", name, "'): using frozen demo data")
    return(demo_read(demo1[[name]]))
  }

  ## 3) reproducible fallback (only if no real or demo file is present)
  if (!quiet) message("get_data('", name,
                      "'): real file not found -- using seeded simulation.")
  set.seed(seed)
  sim <- switch(name,
    butterflies    = .sim_butterflies(),
    biodiv         = .sim_biodiv(),
    taxon          = .sim_taxon(),
    crawley        = .sim_crawley(),
    europe         = .sim_europe(),
    countries_live = .sim_countries_live(),
    leukemia       = .sim_leukemia(),
    limenitis      = .sim_limenitis(),
    doubs          = .sim_doubs(),
    neotoma        = .sim_neotoma(),
    hanta          = .sim_hanta(),
    ssa_factor     = .sim_ssa_factor(),
    biodiv_gv      = .sim_biodiv_gv(),
    insatisf       = .sim_insatisf(),
    pam            = .sim_pam(),
    stop("Unknown dataset: ", name))
  attr(sim, "simulated") <- TRUE
  sim
}

## ---- fallback generators (structure-faithful, seeded) ----------------------
##
##  A fallback is only useful if a script written against it also runs against
##  the real file, so each generator should emit the REAL file's column names.
##  Verified aligned: taxon, crawley, biodiv, hanta, limenitis, leukemia.
##  All generators now emit the REAL files' column names, so a script written
##  against a simulation runs unchanged on the classroom data. (Earlier editions
##  left butterflies, europe, countries_live and neotoma divergent because no
##  script used them; the restored teaching sections do, so they were aligned.)

## Quintana Roo butterflies: species x (taxonomy + wing pattern + 7 sites along
## a successional gradient). Columns match ButterfliesQRoo2.csv exactly.
.sim_butterflies <- function() {
  patterns <- c("Sand","BlackOrange","ShrubOrange","DarkCarpet","OpenContrast",
                "BlackRed","SilentBand","Bark","CanopyOrange","Adelpha",
                "CloseReflect","Tiger")
  sites <- c("HD","SD","GA","YA","MA","OA","PF")           # succession gradient
  ns    <- 128
  genus <- paste0("Gen", sprintf("%02d", sample(30, ns, replace = TRUE)))
  sp    <- paste0("sp", sprintf("%03d", seq_len(ns)))
  opt   <- runif(ns, 1, 7)                     # latent successional optimum
  M <- sapply(seq_along(sites), function(k)
    rpois(ns, lambda = 6 * exp(-0.5 * (k - opt)^2)))
  colnames(M) <- sites
  data.frame(Genus = genus, Sp = sp, Name = paste(genus, sp),
             Pattern = sample(patterns, ns, replace = TRUE), M,
             stringsAsFactors = TRUE)
}

## Biodiversity of countries: RegionCode plus the three variable blocks the
## clustering session uses -- diversity, governance and scientific capacity.
## Column names mirror BiodivCountries.csv exactly (62 columns in the real file;
## the ones any script touches are reproduced here) so a script written against
## the simulation runs unchanged on the classroom data.
.sim_biodiv <- function() {
  regions <- c("APC","CAR","EECA","LAM","MENA","NAM","SSA","WEU")
  n    <- 186
  reg  <- sample(regions, n, replace = TRUE)
  base <- match(reg, regions)
  rich <- function(mult) round(abs(rnorm(n, 50 * mult + 8 * base, 25)))
  dens <- function(mult) round(abs(rnorm(n, 5 * mult + 0.8 * base, 3)), 3)
  gov  <- function(sd = 1) round(rnorm(n, 0, sd), 3)

  out <- data.frame(
    ID = seq_len(n),
    Country    = paste0("C", sprintf("%03d", seq_len(n))),
    RegionCode = factor(reg),
    ## --- diversity block ---
    AmphRich = rich(1.0), Rept_rich = rich(1.2),
    BirdRich = rich(3.0), MamsRich  = rich(1.5),
    CountCentPlantDiv = rpois(n, 2 + base),
    CountEcoregions   = rpois(n, 5 + base),
    CountHotspots     = rpois(n, 1 + base / 2),
    DensAmphRich = dens(1.0), DensRept_rich = dens(1.2),
    DensBirdRich = dens(3.0), DensMamsRich  = dens(1.5),
    ## --- governance block (World Bank WGI: mean and sd per indicator) ---
    VA_Mean = gov(), VA_SD = abs(gov(.3)),
    PS_Mean = gov(), PS_SD = abs(gov(.3)),
    GE_Mean = gov(), GE_SD = abs(gov(.3)),
    RL_Mean = gov(), RL_SD = abs(gov(.3)),
    CC_Mean = gov(), CC_SD = abs(gov(.3)),
    ## --- scientific-capacity block ---
    PapersCapita        = round(rlnorm(n, -1 + base / 6, 1), 3),
    GBIF_OwnToTotal     = round(runif(n), 3),
    CountHerbaria       = rpois(n, 3 + base),
    log10NmbSpc         = round(rnorm(n, 4 + base / 8, .6), 3),
    InternetPenetration = round(runif(n, 5, 95), 2),
    GDRE_Mean10Years    = round(abs(rnorm(n, .5 + base / 20, .4)), 3),
    logGDP2010          = round(rnorm(n, 10 + base / 5, 1), 3),
    logGDP2011          = round(rnorm(n, 10 + base / 5, 1), 3),
    stringsAsFactors = TRUE)

  ## the real file is not complete; imitate that so NA handling gets exercised
  for (j in c("AmphRich","Rept_rich","BirdRich","MamsRich","DensAmphRich",
              "DensRept_rich","DensBirdRich","DensMamsRich","VA_Mean",
              "PapersCapita","logGDP2010"))
    out[sample(n, 3), j] <- NA
  out
}

## Sub-Saharan Africa capacity/governance table used by the factor-analysis
## session: 42 countries x 17 numeric indicators (plus a leading label column X).
.sim_ssa_factor <- function() {
  n <- 42
  f <- function(k) matrix(rnorm(n * k), n, k)
  ## three latent factors -- capacity, biodiversity density, governance --
  ## so factanal() has something real to find
  L <- f(3)
  mk <- function(w, sd) round(as.vector(L %*% w) + rnorm(n, 0, sd), 3)
  data.frame(
    X = paste0("SSA", sprintf("%02d", seq_len(n))),
    PapersCapita        = mk(c(1.0, 0, .2), .4),
    GBIF_Total          = mk(c(0.9, 0, .1), .5),
    CountHerbaria       = mk(c(0.8, 0, .1), .5),
    CountGEFGrants      = mk(c(0.7, .2, .2), .5),
    InternetPenetration = mk(c(0.8, 0, .4), .5),
    CountEcoregions     = mk(c(0.1, .8, 0), .5),
    DensAmphRich        = mk(c(0, 1.0, 0), .4),
    DensRept_rich       = mk(c(0, 0.9, 0), .4),
    DensBirdRich        = mk(c(0, 1.0, 0), .4),
    DensMamsRich        = mk(c(0, 0.9, 0), .4),
    VA_Mean             = mk(c(0, 0, 1.0), .3),
    PS_Mean             = mk(c(0, 0, 0.9), .3),
    GE_Mean             = mk(c(.3, 0, 0.9), .3),
    RL_Mean             = mk(c(.2, 0, 1.0), .3),
    CC_Mean             = mk(c(.1, 0, 0.9), .3),
    logGDP2010          = mk(c(.9, 0, .3), .4),
    logGDP2011          = mk(c(.9, 0, .3), .4),
    stringsAsFactors = TRUE)
}

## US-states social indicators (the "Insatisf" table): 50 states x 9 numerics.
.sim_insatisf <- function() {
  n <- 50
  L <- matrix(rnorm(n * 2), n, 2)
  mk <- function(w, mu, sd) round(mu + as.vector(L %*% w) * sd, 3)
  data.frame(
    State = state.name[seq_len(n)],
    Murder         = mk(c(1.0, .1),  7.8,  4.3),
    Assault        = mk(c(1.0, .1), 170,   83),
    UrbanPop       = mk(c(.2,  .9),  65,   14),
    Rape           = mk(c(.8,  .3),  21,    9),
    Unemploym      = mk(c(.6, -.4),  .06,  .02),
    AllStudents    = mk(c(-.2, .7),  .70,  .09),
    Disabilities   = mk(c(.5, -.2),  .30,  .06),
    LimitedEnglish = mk(c(.1,  .8),  .36,  .12),
    Poor           = mk(c(.9, -.2),  .40,  .10),
    stringsAsFactors = TRUE)
}

## Biodiversity + governance, complete cases only (BiodiversityCountriesBiGv.csv):
## the two-table set-up the canonical-correlation session needs.
.sim_biodiv_gv <- function() {
  regions <- c("APC","CAR","EECA","LAM","MENA","NAM","SSA","WEU")
  n   <- 130
  reg <- sample(regions, n, replace = TRUE)
  k   <- match(reg, regions)
  lat <- rnorm(n)                      # one latent axis both blocks respond to
  data.frame(
    RegionCode = factor(reg),
    ISO3V10    = paste0("C", sprintf("%03d", seq_len(n))),
    AmphRich   = round(abs(rnorm(n,  60 + 8 * k + 20 * lat, 25))),
    Rept_rich  = round(abs(rnorm(n,  70 + 8 * k + 18 * lat, 25))),
    BirdRich   = round(abs(rnorm(n, 200 + 9 * k + 40 * lat, 60))),
    MamsRich   = round(abs(rnorm(n, 120 + 7 * k + 30 * lat, 40))),
    VA_Mean    = round(rnorm(n, -0.2 * lat, 1), 3),
    PS_Mean    = round(rnorm(n, -0.2 * lat, 1), 3),
    GE_Mean    = round(rnorm(n, -0.3 * lat, 1), 3),
    RL_Mean    = round(rnorm(n, -0.3 * lat, 1), 3),
    stringsAsFactors = TRUE)
}

## Crawley "taxon": 4 taxa x 30, seven morphometric variables (Gaussian blobs)
.sim_taxon <- function() {
  g <- 4; per <- 30; p <- 7
  centers <- matrix(rnorm(g * p, 0, 3), g, p)
  X <- do.call(rbind, lapply(seq_len(g), function(k)
    matrix(rnorm(per * p), per, p) + matrix(centers[k, ], per, p, byrow = TRUE)))
  vn <- c("Petals","Internode","Sepal","Bract","Petiole","Leaf","Fruit")
  colnames(X) <- vn
  data.frame(Taxon = factor(rep(c("I","II","III","IV"), each = per)), X)
}

## Crawley "speciesCrawley3": counts + continuous env + soil factor
.sim_crawley <- function() {
  n <- 90
  Tmp    <- runif(n, -2, 12)
  Precip <- rgamma(n, 2, scale = 60)
  pH     <- runif(n, 3.5, 8.5)
  Soil   <- factor(sample(c("clay","loam","sand"), n, TRUE))
  mu_b    <- exp(1.6 + 0.04 * Tmp - 0.06 * (pH - 6)^2)   # positive mean
  Biomass <- rgamma(n, shape = 3, scale = mu_b / 3)       # right-skewed (non-normal)
  Species <- rpois(n, lambda = exp(1.2 + 0.06 * Tmp + 0.15 * (pH - 4)))
  data.frame(Species, Biomass, Tmp, Precip, pH, Soil)
}

## 21 European cities (real public coordinates); columns City, Long, Lat as in
## CitiesEurope.csv.
.sim_europe <- function() {
  data.frame(
    City = c("Athens","Barcelona","Berlin","Brussels","Dublin","Geneva",
             "Helsinki","Lisbon","London","Madrid","Milan","Munich","Oslo",
             "Paris","Prague","Rome","Stockholm","Vienna","Warsaw","Zurich","Hamburg"),
    Long = c(23.73,2.17,13.40,4.35,-6.26,6.14,24.94,-9.14,-0.13,-3.70,9.19,11.58,
             10.75,2.35,14.42,12.50,18.07,16.37,21.01,8.54,9.99),
    Lat  = c(37.98,41.39,52.52,50.85,53.35,46.20,60.17,38.72,51.51,40.42,45.46,
             48.14,59.91,48.86,50.08,41.90,59.33,48.21,52.23,47.38,53.55),
    stringsAsFactors = TRUE)
}

## Ordinal liveability rankings: 13 countries ranked on six criteria. Columns
## match CountriesToLive.csv (living, climate, food, security, hospitality,
## infrastructure), each a permutation of 1..13.
.sim_countries_live <- function() {
  countries <- c("Norway","Switzerland","Canada","Germany","Japan","France",
                 "Spain","Portugal","Mexico","Brazil","India","Kenya","Vietnam")
  crit <- c("living","climate","food","security","hospitality","infrastructure")
  M <- sapply(crit, function(.) sample(seq_along(countries)))
  data.frame(Country = countries, M, stringsAsFactors = TRUE)
}

## Leukemia-like expression: genes x samples, 3 subtypes with block signatures
.sim_leukemia <- function() {
  types <- rep(c("ALL","AML","CLL"), each = 20)          # 60 samples
  g <- 300
  base <- matrix(rnorm(g * length(types), 6, 1), g, length(types))
  for (t in unique(types)) {
    idx <- which(types == t)
    sig <- sample(g, 40)
    base[sig, idx] <- base[sig, idx] + rnorm(1, 3, 0.3)
  }
  colnames(base) <- paste0(types, ".", ave(seq_along(types), types, FUN = seq_along))
  rownames(base) <- paste0("g", seq_len(g))
  base
}

## Synthetic grayscale "wing" image matrix (rank structure for SVD/NMF)
.sim_limenitis <- function(nr = 69, nc = 100) {
  x <- seq(-3, 3, length.out = nc); y <- seq(-2, 2, length.out = nr)
  G <- outer(y, x, function(yy, xx)
    exp(-(xx^2)/3 - (yy^2)/1.2) * (1 + 0.6 * cos(3 * xx)) )
  G <- G + 0.15 * outer(y, x, function(yy, xx) sin(2 * xx) * cos(2 * yy))
  G <- (G - min(G)) / (max(G) - min(G))
  round(G, 4)
}

## Doubs-like river gradient: fish counts, env, spatial coords
.sim_doubs <- function() {
  n <- 30; s <- 27
  grad <- seq(0, 1, length.out = n)             # up- to down-stream
  opt  <- seq(0, 1, length.out = s)
  fish <- sapply(seq_len(s), function(k)
    rpois(n, lambda = 8 * exp(-((grad - opt[k])^2) / 0.02)))
  colnames(fish) <- paste0("Sp", seq_len(s))
  env <- data.frame(
    alt = 900 - 850 * grad + rnorm(n, 0, 10),
    slo = exp(4 - 4 * grad + rnorm(n, 0, .3)),
    flo = 10 + 400 * grad + rnorm(n, 0, 20),
    pH  = 8 - grad + rnorm(n, 0, .1),
    oxy = 12 - 4 * grad + rnorm(n, 0, .5),
    nit = 1 + 3 * grad + rnorm(n, 0, .2))
  xy <- data.frame(x = cumsum(runif(n, .5, 1.5)), y = cumsum(rnorm(n, 0, .4)))
  list(fish = as.data.frame(fish), env = env, xy = xy)
}

## Neotoma-like table: 615 specimens, morphometrics + environment. Columns match
## NeotomaMorphoEnvir.csv: ID, sp, long, lat, No..cat, Code, d1-d14 (dental),
## v1-v28, m1-m14, L1-L20 (landmarks) and BIO1-BIO19 (bioclim).
.sim_neotoma <- function() {
  n  <- 615
  sp <- factor(sample(paste0("N", 1:15), n, TRUE))
  k  <- as.integer(sp)
  blk <- function(pre, m, mu, sd) {
    X <- sapply(seq_len(m), function(j) rnorm(n, mu + 0.2 * k, sd))
    colnames(X) <- paste0(pre, seq_len(m)); X
  }
  bio <- sapply(1:19, function(j) rnorm(n, 100 + 5 * k, 30))
  colnames(bio) <- paste0("BIO", 1:19)
  data.frame(ID = seq_len(n), sp = sp,
             long = runif(n, -115, -70), lat = runif(n, 20, 45),
             `No..cat` = seq_len(n), Code = paste0("C", k),
             blk("d", 14, 10, 2), blk("v", 28, 8, 1.5),
             blk("m", 14, 12, 2), blk("L", 20, 0, 1), bio,
             check.names = FALSE, stringsAsFactors = TRUE)
}

## Hantavirus-like presence/absence for ENM (binomial on two bioclim predictors)
.sim_hanta <- function() {
  n <- 400
  bio_1  <- runif(n, -5, 30); bio_12 <- runif(n, 100, 2500)
  eta <- -2 + 0.15 * bio_1 - 0.0015 * bio_12 - 0.004 * (bio_1 - 15)^2
  Sp  <- rbinom(n, 1, 1 / (1 + exp(-eta)))
  data.frame(Sp = Sp, bio_1 = bio_1, bio_12 = bio_12)
}

## Presence-absence matrix: sites x species along a latitudinal gradient
.sim_pam <- function() {
  nsite <- 80; nsp <- 120
  lat <- seq(0, 1, length.out = nsite)
  opt <- runif(nsp); br <- runif(nsp, 0.05, 0.25)
  prob <- exp(-outer(lat, opt, `-`)^2 / (2 * rep(br, each = nsite)^2))
  pam <- matrix(as.integer(runif(nsite * nsp) < prob), nsite, nsp)
  colnames(pam) <- paste0("sp", seq_len(nsp))
  data.frame(lat = lat, pam)
}

message("MS_LJMR utils loaded. get_data(<name>), open_dev(), with_fig(), need().")

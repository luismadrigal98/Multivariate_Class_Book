# Reproducible environments for MS_LJMR

Two conda environments live here:

| File | Use it when |
|---|---|
| `environment.yml` | Default. CPU-only TensorFlow. Works on Linux, macOS (Intel + Apple Silicon) and Windows. |
| `environment-gpu.yml` | Linux with an NVIDIA card and a CUDA 12-capable driver. |

Both cover the whole repository in one env: the multivariate book's R stack, the
`UnitCART`/`UnitCARET` machine-learning stack, and — the reason this exists — the
R `keras` + Python TensorFlow pair that `UnitNN` and Track B of `R/21_NN_DL.R`
need.

## Why this is needed

`R/21_NN_DL.R` has two tracks. Track A uses `nnet` and always runs. Track B is
the real deep-learning path, and it is guarded:

```r
if (requireNamespace("keras", quietly = TRUE)) { ... } else { message("keras not installed — Track B is a template only.") }
```

Without a Python TensorFlow underneath R's `keras`, that guard always takes the
`else` branch, so the deep-learning half of the chapter never executes. Installing
`keras` from CRAN alone does **not** fix it — the R package is a binding, and it
needs a matching Python interpreter with TensorFlow in it. Conda is the least
painful way to get both halves pinned together, which is exactly what these files do.

## Setup

```bash
# 1. create (mamba is far faster than conda for an env this size)
mamba env create -f envs/environment.yml     # or: conda env create -f ...
conda activate msljmr

# 2. point reticulate at THIS env's python, not a system one
R -e 'reticulate::use_condaenv("msljmr", required = TRUE); reticulate::py_config()'

# 3. verify the R -> Python -> TensorFlow chain end to end
R -e 'library(tensorflow); tf$constant("hello"); tf$config$list_physical_devices()'

# 4. two CRAN-only comparators used in UnitCARET/src/main.class.R have no
#    conda-forge feedstock; add them inside the activated env:
R -e 'install.packages(c("hda","HiDimDA"), repos="https://cloud.r-project.org")'
```

Step 3 is the one that matters. If it prints a tensor and a device list, Track B
will run. If it errors, `reticulate` is still bound to the wrong Python — rerun
step 2 and restart R.

## Version pinning, and why

`r-keras`, `r-tensorflow` and `tensorflow` are all pinned to the **2.16** line.
This is deliberate:

- All of `UnitNN/src/create_model_*.R` uses the Keras 2 API — `keras_model_sequential()`,
  `layer_conv_2d()`, `%>%` chaining into `compile()`. Keras 3 reorganised much of
  that surface and changed default behaviours (notably around multi-output models,
  which the CRISPR network relies on: five named heads off one shared trunk).
- The R binding and the Python engine must agree. A 2.16 R package against a
  Keras 3 Python install fails at import, usually with an opaque reticulate error.

If you deliberately want Keras 3, bump all three together and expect to port the
model constructors.

`r-base` is pinned to 4.4 to match `UnitCARET/R_session_info.txt` (R 4.4.0), the
version the recorded results were produced under.

## GPU notes

`UnitNN/main.R` sources `src/Sec_unable_gpu_processing.R` as its first action,
which hard-disables CUDA:

```r
Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")
tf$config$set_visible_devices(list(), 'GPU')
```

That was the right call for portability — but it means the GPU env alone changes
nothing. Comment out that `source()` line to actually use the card.

Everything in `UnitNN` is CPU-tractable anyway: images are downsampled to
120×160, the largest model (CNN2) is ~9.3 M parameters, and the transfer-learning
model freezes ResNet-50's 23.6 M weights and trains only 525 829 of them.

## Status of these files

The two `.yml` files were authored against this repository's actual dependency
surface — every package listed is imported somewhere in `R/`, `UnitCART/`,
`UnitCARET/` or `UnitNN/`. They have **not** been solved or built in the sandbox
where they were written: that machine's network allowlist blocks
`repo.anaconda.com` and `conda.anaconda.org`, so `conda env create` cannot run
there. `r-keras`'s presence and 2.16.1 version on conda-forge were confirmed
through the Anaconda API. Expect to run `mamba env create` once on your own
machine and adjust if the solver objects to a pin.

## If you would rather not use conda

The R packages other than `keras`/`tensorflow` are all available from CRAN, and
on Debian/Ubuntu most are packaged as `r-cran-*`. The deep-learning pair is the
only genuinely awkward dependency, and `keras::install_keras()` will build it a
private virtualenv if you prefer that to a full conda env:

```r
install.packages("keras")
keras::install_keras(method = "virtualenv", version = "2.16")
```

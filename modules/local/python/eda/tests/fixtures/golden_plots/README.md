# EDA Golden Plots

PNG and SVG files saved by `equivalence.nf.test` on its first successful run.
They capture the visual output of the current EDA implementation so you can
side-by-side compare them against T18b (polars-bio rewrite) plots.

## How they are populated

The test's `then {}` block copies all `plots/*.png` and `plots/*.svg` files here
automatically the first time it runs successfully (when the directory contains no
PNGs).  Subsequent runs (including T18b) do NOT overwrite them.

## Usage in T18b

After T18b is implemented, run the equivalence test and compare the new plots
against the files here.  Any visual difference is an intentional behavior change
that should be documented in the T18b commit message.

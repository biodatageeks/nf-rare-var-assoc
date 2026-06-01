# `calc_f_outliers` test fixtures

Two hand-rolled plink2-style `.het` files for `CALCULATE_F_OUTLIERS` unit tests.
Both fixtures share the same 22 samples and the same `F` values; they differ only
in whether the header carries an `FID` column, which exercises the column-selection
branch at [`main.nf:44-47`](../../main.nf#L44-L47).

## Sample layout (both fixtures)

| Sample | F | Role |
|---|---|---|
| `SAMPLE_01` .. `SAMPLE_20` | `0.01` | Normal — drawn from a degenerate "tight" distribution |
| `SPIKE_HIGH` | `0.5` | Engineered outlier on the upper tail |
| `SPIKE_LOW` | `-0.5` | Engineered outlier on the lower tail |

`O(HOM)`, `E(HOM)`, and `OBS_CT` are filler values that the script ignores; only
the `F` column is read.

## Why these spike values

With N=22, normals all at 0.01, spikes at +/-0.5, `numpy.std` (ddof=0) gives:

```
mean F  = 0.009091
std F   = 0.150783
2-sigma bounds  = (-0.2925, 0.3107)   -> both spikes outside, all normals inside
10-sigma bounds = (-1.4987, 1.5169)   -> both spikes inside (no outliers)
```

So `range_stds=2` yields exactly 2 outliers (the two spikes), and `range_stds=10`
yields zero. These two regimes are what the tests assert.

## Files

- `test_het.het` — full plink2 header `#FID IID O(HOM) E(HOM) OBS_CT F`; exercises
  the `FID in data.columns` branch (output rows are `FID\tIID`).
- `test_het_no_fid.het` — header without `FID` (`#IID O(HOM) E(HOM) OBS_CT F`);
  exercises the fallback branch (output rows are `IID\tIID`).

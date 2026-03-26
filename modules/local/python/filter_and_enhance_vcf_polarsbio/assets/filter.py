import argparse
from pathlib import Path
import polars as pl
import polars_bio as pb


parser = argparse.ArgumentParser()
parser.add_argument("--input-vcf-path", help="Input VCF file")
parser.add_argument("--output-vcf-path", help="Output VCF file")
parser.add_argument("--samples-file", help="File containing sample names")
parser.add_argument("--qual-min", type=int, default=25, help="Minimum quality score")
parser.add_argument("--avg-gq-min", type=int, default=25, help="Minimum average genotype quality")
parser.add_argument("--avg-dp-min", type=int, default=25, help="Minimum average depth")
parser.add_argument("--avg-dp-max", type=int, default=200, help="Maximum average depth")
parser.add_argument("--sample-gq-min", type=int, default=20, help="Minimum sample genotype quality")
parser.add_argument("--sample-dp-min", type=int, default=20, help="Minimum sample depth")
parser.add_argument("--sample-dp-max", type=int, default=250, help="Maximum sample depth")
parser.add_argument("--calc-ds-min-gq", type=int, default=1, help="Minimum genotype quality for DS calculation")

args = parser.parse_args()


samples = [s.strip() for s in Path(args.samples_file).read_text().splitlines() if s.strip()]


pb.register_vcf(
    args.input_vcf_path,
    name="vcf_table",
    # info_fields=[],
    format_fields=["GT", "DP", "GQ", "PL"],
    samples=samples,
)

schema_df = pb.sql("DESCRIBE vcf_table").collect()
name_col = next((c for c in ["column_name", "field_name", "name"] if c in schema_df.columns), None)
if name_col is None:
    raise ValueError(f"Could not infer schema column-name field from DESCRIBE output columns: {schema_df.columns}")

input_columns = schema_df[name_col].to_list()
passthrough_cols = [c for c in input_columns if c not in ["genotypes", "AN", "AC", "AF"]]


def quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


select_passthrough_sql = ", ".join("i." + quote_ident(col) for col in passthrough_cols)


FILTER_SQL = f"""\
WITH site_filtered AS (
  SELECT *
  FROM vcf_table
  WHERE qual IS NOT NULL
    AND qual = qual              -- filters out NaN
    AND qual >= {args.qual_min}
    AND list_avg(genotypes."GQ") >= {args.avg_gq_min}
    AND list_avg(genotypes."DP") >= {args.avg_dp_min}
    AND list_avg(genotypes."DP") <= {args.avg_dp_max}
),
indexed AS (
  SELECT
    *,
    ROW_NUMBER() OVER () AS variant_idx
  FROM site_filtered
),
samples_unnested AS (
  SELECT
    variant_idx,
    UNNEST(genotypes."GT") AS gt,
    UNNEST(genotypes."GQ") AS gq,
    UNNEST(genotypes."DP") AS dp,
    UNNEST(genotypes."PL") AS pl
  FROM indexed
),
sample_processed AS (
  SELECT
    *,
    (gt IN ('0', '1', '2')) AS gt_haploid,
    (gq IS NOT NULL AND dp IS NOT NULL
      AND gq >= {args.sample_gq_min}
      AND dp >= {args.sample_dp_min}
      AND dp <= {args.sample_dp_max}) AS is_good,
    (gt IN ('0/0', '0|0', '0')) AS is_hom_ref,
    COALESCE(pl[1], 0) AS pl0,
    COALESCE(pl[2], 0) AS pl1,
    COALESCE(pl[3], 0) AS pl2,
    (pl IS NOT NULL AND array_length(pl) = 3) AS pl_valid,
    (gq IS NOT NULL AND gq >= {args.calc_ds_min_gq}) AS gq_sufficient
  FROM samples_unnested
),
pl_corrected AS (
  SELECT
    *,
    (is_hom_ref AND pl_valid AND pl0 = 0 AND pl1 = 0 AND pl2 = 0 AND gq_sufficient) AS needs_correction,
    POWER(10.0, -COALESCE(gq, 0) / 10.0) AS p_wrong
  FROM sample_processed
),
pl_calc AS (
  SELECT
    *,
    (-1.0 + SQRT(1.0 + 4.0 * p_wrong)) / 2.0 AS x
  FROM pl_corrected
),
final_pl AS (
  SELECT
    *,
    CASE WHEN needs_correction THEN
      CAST(LEAST(255, GREATEST(0, ROUND(-10.0 * LOG10(GREATEST(1e-25, 1.0 - x - x*x))))) AS INTEGER)
    ELSE pl0 END AS final_pl0,
    CASE WHEN needs_correction THEN
      CAST(LEAST(255, GREATEST(0, ROUND(-10.0 * LOG10(GREATEST(1e-25, x))))) AS INTEGER)
    ELSE pl1 END AS final_pl1,
    CASE WHEN needs_correction THEN
      CAST(LEAST(255, GREATEST(0, ROUND(-10.0 * LOG10(GREATEST(1e-25, x*x))))) AS INTEGER)
    ELSE pl2 END AS final_pl2,
    CASE WHEN needs_correction THEN
      [
        CAST(LEAST(255, GREATEST(0, ROUND(-10.0 * LOG10(GREATEST(1e-25, 1.0 - x - x*x))))) AS INTEGER),
        CAST(LEAST(255, GREATEST(0, ROUND(-10.0 * LOG10(GREATEST(1e-25, x))))) AS INTEGER),
        CAST(LEAST(255, GREATEST(0, ROUND(-10.0 * LOG10(GREATEST(1e-25, x*x))))) AS INTEGER)
      ]
    ELSE pl
    END AS pl_final
  FROM pl_calc
),
with_ds AS (
  SELECT
    *,
    CASE WHEN final_pl0 < 255 THEN POWER(10.0, -final_pl0/10.0) ELSE 0.0 END AS l_ref,
    CASE WHEN final_pl1 < 255 THEN POWER(10.0, -final_pl1/10.0) ELSE 0.0 END AS l_het,
    CASE WHEN final_pl2 < 255 THEN POWER(10.0, -final_pl2/10.0) ELSE 0.0 END AS l_alt
  FROM final_pl
),
ds_calc AS (
  SELECT
    *,
    CASE
      WHEN is_good THEN gt
      ELSE CASE WHEN gt_haploid THEN '.' ELSE './.' END
    END AS gt_final,
    CASE
      WHEN pl_valid AND gq_sufficient AND (l_ref + l_het + l_alt) > 0 THEN
        CASE
          WHEN ((l_het + 2.0 * l_alt) / (l_ref + l_het + l_alt)) > 0
               AND ((l_het + 2.0 * l_alt) / (l_ref + l_het + l_alt)) < 0.0001
          THEN 0.0
          ELSE (l_het + 2.0 * l_alt) / (l_ref + l_het + l_alt)
        END
      ELSE NULL
    END AS ds
  FROM with_ds
),
genotypes_aggregated AS (
  SELECT
    variant_idx,
    STRUCT(
      array_agg(gt_final) AS GT,
      array_agg(gq) AS GQ,
      array_agg(dp) AS DP,
      array_agg(pl_final) AS PL,
      array_agg(CAST(ds AS FLOAT)) AS DS
    ) AS genotypes
  FROM ds_calc
  GROUP BY variant_idx
)
SELECT {select_passthrough_sql},
  vcf_an(g.genotypes."GT") AS "AN",
  vcf_ac(g.genotypes."GT", i.alt) AS "AC",
  vcf_af(g.genotypes."GT", i.alt) AS "AF",
  g.genotypes
FROM indexed i join genotypes_aggregated g ON i.variant_idx = g.variant_idx
"""

# other aggreagation apart from array_agg currently not supported by FusedArrayTransform
#
#   SUM(
#    CASE
#      WHEN gt_final IN ('0/0', '0|0', '0') THEN 0
#      WHEN gt_final IN ('0/1', '1/0', '0|1', '1|0', '1') THEN 1
#      WHEN gt_final IN ('1/1', '1|1') THEN 2
#      ELSE 0
#    END
#  )::INTEGER AS "AC",
#  SUM(
#    CASE
#      WHEN gt_final IN ('0/0', '0|0', '0', '0/1', '1/0', '0|1', '1|0', '1', '1/1', '1|1') THEN 2
#      ELSE 0
#    END
#  )::INTEGER AS "AN",

# another option:
#  list_sum(array_agg(
#    CASE
#      WHEN gt_final IN ('0/0', '0|0', '0') THEN 0
#      WHEN gt_final IN ('0/1', '1/0', '0|1', '1|0', '1') THEN 1
#      WHEN gt_final IN ('1/1', '1|1') THEN 2
#      ELSE 0
#    END)) AS "AC"
# but list_sum would have to be added to datafusion-bio-formats



lf = pb.sql("EXPLAIN select * FROM (" + FILTER_SQL + ")")

optim_used = False
with pl.Config(fmt_str_lengths=None, tbl_rows=-1, tbl_cols=-1, fmt_table_cell_list_len=-1, tbl_width_chars=20000):
    logical_plan = lf.collect()["plan"][0]
    physical_plan = lf.collect()["plan"][1]
    print(f"sql result:")
    print("")
    print("  logical plan:")
    print(logical_plan)
    print("")
    print("  physical plan:")
    print(physical_plan)
    if "FusedArrayTransform" in physical_plan:
        optim_used = True

print("")
print(f"optim_used = {optim_used}")

if not optim_used:
    print("Error: for some reason FusedArrayTransform optimization is not used. This optimization is mandatory for this query. Interrupting..")
    exit(1)

lf = pb.sql(FILTER_SQL)
pb.sink_vcf(lf, args.output_vcf_path)
print(f"Wrote: {args.output_vcf_path}")

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
passthrough_cols = [c for c in input_columns if c != "genotypes"]


def quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


select_passthrough_sql = ", ".join(f"vcf_table.{quote_ident(col)}" for col in passthrough_cols)

FILTER_SQL = f"""\
SELECT {select_passthrough_sql}, \
  named_struct(\
    'GT', vcf_set_gts(vcf_table."genotypes"."GT", \
      list_and(list_and(list_gte(vcf_table."genotypes"."GQ", CAST({args.sample_gq_min} AS INT)), list_gte(vcf_table."genotypes"."DP", CAST({args.sample_dp_min} AS INT))), \
               list_lte(vcf_table."genotypes"."DP", CAST({args.sample_dp_max} AS INT))), './.'), \
    'GQ', vcf_table."genotypes"."GQ", \
    'DP', vcf_table."genotypes"."DP", \
    'PL', vcf_table."genotypes"."PL"\
  ) AS genotypes \
FROM vcf_table \
WHERE vcf_table."qual" >= {args.qual_min}
  AND list_avg(vcf_table."genotypes"."GQ") >= {args.avg_gq_min}
  AND list_avg(vcf_table."genotypes"."DP") >= {args.avg_dp_min}
  AND list_avg(vcf_table."genotypes"."DP") <= {args.avg_dp_max}"""

lf = pb.sql("EXPLAIN select * FROM (" + FILTER_SQL + ")")

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

lf = pb.sql(FILTER_SQL)
pb.sink_vcf(lf, args.output_vcf_path)
print(f"Wrote: {args.output_vcf_path}")

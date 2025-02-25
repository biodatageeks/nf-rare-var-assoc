import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_sep", help="input file path", type=str, required=True)
parser.add_argument("--output", help="output file path", type=str, required=True)
parser.add_argument("--output_sep", help="output file path", type=str, required=True)
parser.add_argument("--inputs", nargs="+", help="input files paths", required=True)
args = parser.parse_args()

dfs = []
for input_path in args.inputs:
    df = pd.read_csv(input_path, sep=args.input_sep)
    dfs.append(df)

pd.concat(dfs).to_csv(args.output, sep=args.output_sep)

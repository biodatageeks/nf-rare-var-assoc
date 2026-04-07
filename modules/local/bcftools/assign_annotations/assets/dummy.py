#!/usr/bin/env python3

if __name__ == "__main__":
    # create dummy expected output files
    with open("dummy.annotations", "w") as f:
        # f.write("key\tSymbol\tConsequence\n")  # no header
        f.write("chr1_1000_A_G\tGENE1\tmissense_variant\n")
    with open("dummy.setlist", "w") as f:
        # f.write("symbol\tchrom\tpos\tvariants\n")  # no header
        f.write("GENE1\tchr1\t1000\tchr1_1000_A_G\n")

    print(f"Ok")

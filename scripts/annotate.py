# annotate integration sites with features from one or more gff files
# add two columns for each gff - the closest feature to each integration site and 
# the distance to that feature

import argparse
import tempfile

import pandas as pd
from pybedtools import BedTool

def main():

    # parse command line arguments
    args = parse_args()

    # read input file
    df = pd.read_csv(args.input, sep="\t", header=0)
    
    # generate BedTool
    ints = integration_bedtool(df, args.virus)

    # get closest feature from each gff
    for gff in args.gff:

        # get closest feature
        df = get_closest_feature(ints, gff, df)

    # write output
    df.to_csv(args.output, sep="\t", index=False)


def parse_args():

    parser = argparse.ArgumentParser(description="Annotate integration sites with features from one or more gff files")
    parser.add_argument("-i", "--input", help="input file of integration sites", required=True)
    parser.add_argument("-g", "--gff", help="gff file(s) to annotate with", nargs="+")
    parser.add_argument("-v", "--virus", help="Use viral coordinates instead of host coordinates", action="store_true")
    parser.add_argument("-o", "--output", help="output file name", required = True)
    args = parser.parse_args()

    return args

def get_closest_feature(ints, gff, df):

    gt = BedTool(gff).sort()

    # get closest feature
    closest = ints.closest(gt, d=True, t="all")

    # combine closest feature and distance if there is more than one match
    # and add to dataframe
    closest_df = (pd.read_csv(closest.fn, sep="\t", header=None)
                    .rename(columns={3: "ID", 12: "Feature", 13: "Distance"})
                    .groupby("ID")
                    .agg(Feature = ("Feature", lambda x: ','.join(str(i) for i in x)), 
                         Distance = ("Distance", lambda x: ','.join(str(i) for i in x)))
                    .rename(columns={"Feature": gff + "_Feature", "Distance": gff + "_Distance"})
                    .reset_index()
                  )

    # join with original dataframe
    df = df.merge(closest_df, how="left", left_on="SiteID", right_on="ID")
    print(df)

    return df



def integration_bedtool(df, virus=False):

    # get columns
    if virus:
        bed = zip(df["Virus"], df["VirusStart"], df["VirusStop"], df["SiteID"])
    else:
        bed = zip(df["Chr"], df["IntStart"], df["IntStop"], df["SiteID"])

    # return BedTool
    return BedTool(list(bed))


if __name__ == "__main__":
    main()
import os
import tempfile

from pandas.core.frame import DataFrame
import pytest
import pandas as pd
from pybedtools import BedTool
from scripts.annotate import get_closest_feature, parse_args, integration_bedtool

@pytest.fixture
def ints():
    ints = [("chr1", 100, 200, "AAV2", 10, 20, "Site1"), 
            ("chr1", 300, 400, "AAV2", 10, 20, "Site2"), 
            ("chr2", 500, 600, "AAV2", 10, 20, "Site3")]
    return ints

@pytest.fixture
def gff_bedtool():
    gff = [("chr1", "HAVANA", "gene", 80, 90, ".", "+", ".", "Gene1"),
           ("chr1", "HAVANA", "gene", 300, 400, ".", "+", ".", "Gene2"),]
    return BedTool(gff)

@pytest.fixture
def ints_bedtool(ints: list):
    return BedTool(ints)

@pytest.fixture
def df(ints: list):
    return pd.DataFrame(ints, columns=["Chr", "IntStart", "IntStop", 
                                       "Virus", "VirusStart", "VirusStop", 
                                       "SiteID"])

@pytest.fixture
def df_empty():
    return pd.DataFrame(columns=["Chr", "IntStart", "IntStop", 
                                       "Virus", "VirusStart", "VirusStop", 
                                       "SiteID"])

#### get_closest_feature ####

def test_get_closest_feature(df: DataFrame, ints_bedtool: BedTool, gff_bedtool: BedTool):
    # expected output
    expected_output = df.copy()
    fn = os.path.basename(gff_bedtool.fn)
    expected_output[f"{fn}_Feature"] = ["Gene1", "Gene2", "."]
    expected_output[f"{fn}_Distance"] = [11, 0, -1]

    # test function
    output = get_closest_feature(ints_bedtool, gff_bedtool, df)

    # check output
    assert all(output.columns == expected_output.columns)
    assert all(output == expected_output)

def test_get_closest_feature_empty_gff(df: DataFrame, ints_bedtool: BedTool):
    
    # empty bedtool
    gff = BedTool([])

    # expected output
    expected_output = df.copy()
    fn = os.path.basename(gff.fn)
    expected_output[f"{fn}_Feature"] = [".", ".", "."]
    expected_output[f"{fn}_Distance"] = [-1, -1, -1]

    # test function
    output = get_closest_feature(ints_bedtool, gff, df)

    # check output
    assert all(output.columns == expected_output.columns)
    assert all(output == expected_output)

def test_get_closest_feature_empty_ints(df_empty: DataFrame, gff_bedtool: BedTool):
    
    # empty 
    ints_bedtool = BedTool([])

    # expected output
    expected_output = df_empty.copy()
    fn = os.path.basename(gff_bedtool.fn)
    expected_output[f"{fn}_Feature"] = []
    expected_output[f"{fn}_Distance"] = []

    # test function
    output = get_closest_feature(ints_bedtool, gff_bedtool, df_empty)

    # check output
    assert all(output.columns == expected_output.columns)
    assert all(output == expected_output)

#### parse_args ####
    
def test_parse_args_1():
    # Test that required arguments are present
    with pytest.raises(SystemExit):
        parse_args([])
    with pytest.raises(SystemExit):
        parse_args(["-o", "output.txt"])
    with pytest.raises(SystemExit):
        parse_args(["-i", "input.txt"])
    
def test_parse_args_2():
    # Test that optional arguments are correctly parsed
    args = parse_args(["-i", "input.txt", "-g", "file1.gff", "file2.gff", "-v", "-o", "output.txt"])
    assert args.input == "input.txt"
    assert args.gff == ["file1.gff", "file2.gff"]
    assert args.virus == True
    assert args.output == "output.txt"
    
def test_parse_args_3():
    args = parse_args(["-i", "input.txt", "-o", "output.txt"])
    assert args.input == "input.txt"
    assert args.gff == None
    assert args.virus == False
    assert args.output == "output.txt"
    
def test_parse_args_4():
    args = parse_args(["-i", "input.txt", "-g", "file1.gff", "-o", "output.txt"])
    assert args.input == "input.txt"
    assert args.gff == ["file1.gff"]
    assert args.virus == False
    assert args.output == "output.txt"

#### integration_bedtool ####

def test_integration_bedtool_host(df: DataFrame, ints: list):
    # Test with virus=False
    expected_output = BedTool([i[0:3]+i[6:7] for i in ints])
    assert integration_bedtool(df, virus=False) == expected_output

def test_integration_bedtool_virus(df: DataFrame, ints: list):
    # Test with virus=False
    expected_output = BedTool([i[3:6]+i[6:7] for i in ints])
    assert integration_bedtool(df, virus=True) == expected_output

def test_integration_bedtool_empty(df: DataFrame, ints: list):
    # Test with virus=False
    expected_output = BedTool([i[0:3]+i[6:7] for i in ints])
    assert integration_bedtool(df, virus=False) == expected_output
import pytest

def test_merge_parquet_and_convert():
    """
    Test if merging parquet files and converting to ROOT format with the helper scripts works
    """
    import os

    if not os.path.exists("./tests/samples/parquet_files/merged/"):
        os.makedirs("./tests/samples/parquet_files/merged/")
    os.system("cp ./tests/test_cat.json ./higgs_dna/category.json")
    x = os.system("python ./scripts/postprocessing/merge_parquet.py --source ./tests/samples/parquet_files/singles/ --target ./tests/samples/parquet_files/merged/ --cats test_cat.json --skip-normalisation")  
    x += os.system("python ./scripts/postprocessing/convert_parquet_to_root.py ./tests/samples/parquet_files/merged/merged.parquet ./tests/samples/parquet_files/merged/merged.root mc --cats category.json") 

    assert x==0

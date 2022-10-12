import pytest

def test_merge_parquet_and_convert():
    """
    Test if merging parquet files and converting to ROOT format with the helper scripts works
    """
    import os

    if not os.path.exists("./tests/samples/parquet_files/merged/"):
        os.makedirs("./tests/samples/parquet_files/merged/")
    x = os.system("python ./scripts/merge_parquet.py ./tests/samples/parquet_files/singles/ ./tests/samples/parquet_files/merged/")  
    x += os.system("python ./scripts/convert_parquet_to_root.py ./tests/samples/parquet_files/merged/merged.parquet ./tests/samples/parquet_files/merged/merged.root") 

    assert x==0
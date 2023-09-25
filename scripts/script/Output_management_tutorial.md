# Output management scripts instructions

## Disclamer

This tutorial makes a few assumptions that show the weak points of the routine, if someone feels like proposing/making improvements they're more than welcome:
* The output to be analysied are supposed to be in the format that is given by default by HiggsDNA, i. e.:
```
        out_dir -> sample_1 -> .parquet_1
                |                 ...
                |              .parquet_n
                -> sample_2 -> .parquet_1
                |                 ...
                |              .parquet_n
                ...
                |
                -> sample_n -> .parquet_1
                |                 ...
                |              .parquet_n
```
* The script distinguish between MC and data samples by looking for `M125` in the file name (not ideal since not everybody uses standard model Higgs mass in the samples).
* The sample should contain the year to which they correspond in the name at the end, i. e.: `ggh_M125_2017`
* In the main script, `prepare_output_file.py`, are specified two dictionaries that contain a map between the name of the sample and the process that contains, and the systematics variations that are contained in the samples as separate collections (not simply as weight collections), you have to update these properly before using the script.
```
process_dict = {
    "ggh_M125_2017": "ggh",
    "tth_M125_2017": "tth",
    "vbf_M125_2017": "vbf",
    "vh_M125_2017": "vh",
}

var_dict = {
    "NOMINAL": "nominal",
    "FNUFUp": "FNUF_up",
    "FNUFDown": "FNUF_down",
    "ShowerShapeUp": "ShowerShape_up",
    "ShowerShapeDown": "ShowerShape_down",
}
```

## Merging step
During this step the main script calls `merge_parquet.py` multiple times. The starting point is the output of HiggsDNA, i.e. `out_dir -> sample_1/n` in the disclamer graph. These directory **must** contain only `.parquet` files that have to be merged. 
The script will create a new directory called `merged` under `out_dir`, if this directory already exists it will throw an error and exit.
When converting the data (in my case they were split per era, `DoubleEG_B_2017`, `DoubleEG_C_2017` etc.) the script will put them in a new directory `Data_2017` and then merge again the output in a `.parquet` called `allData_2017.parquet`.
During this step the events are also split into categories according to the boundaries defined in the `cat_dict` in the main file.

## Root step 

During this step the script calls multiple times the script `convert_parquet_to_root.py`. The arguments to pass to the script, for instance if you want the systematic variation included in the output `root tree` are specified when calling `prepare_output_file.py` using `--args "--do_syst"`.
As before the script creates a new called `root` under `out_dir`, if this directory already exists it will throw an error and exit. In the script there is a dictionary called `outfiles` that contains the name of the output root file that will be created according to the process tipe, if the wf is run using the main script this correspond to the proces containd in `process_dict`.

## Workspace step

During this step the main script uses multiple time the `Flashgg_FinalFit`, it moves to the directory defined in the `--final_fit` option (improvable) and uses the Tree2WS script there on the content of the `root` directory previously created. The output is stored in `out_dir/root/smaple_name/ws/`.

## Commands

The workflow is meant to be run in one go using the `prepare_output_file.py` script, it can be also split in different steps or run with the single auxiliary files but it can result a bit cumbersome.

To run everything starting from the output of HiggsDNA with categories and systematic variatrion one can use 
```
python3 prepare_output_file.py --input [path to output dir] --merge --root --ws --syst --cats --args "--do_syst"
```
and everithing should run smoothly, it does for me at least (modulo the things contained in the disclamer).
Some options can be removed. If you want to use `--syst` and `--root` you should also add `--args "--do_syst"`.

The complete list of options for the main file is here:

* "--merge", "Do merging of the .parquet files"
* "--root", "Do root conversion step"
* "--ws", "Do root to workspace conversion step"
* "--ws_config", "configuration file for Tree2WS, as it is now it must be stored in Tree2WS directory in FinalFit",
* "--final_fit", "FlashggFinalFit path" # the default is just for me, it should be changed but I don't see a way to make this generally vali* "--syst",
    dest="syst", "Do systematics variation treatment"
* "--cats", ="Split into categories",
* "--args", "additional options for root converter: --do_syst, --notag",
* "--verbose", "verbose lefer for the logger: INFO (default), DEBUG",

The merging step can also be run separately using:

```
python3 merge_parquet.py --source [path to the directory containing .paruets] --target [target directory path] --cats [cat_dict]
```

the script works also without the `--cats` option, it creates a dummy selection of Pt > -1 and call the category `NOTAG`.

Same for the root step:
```
python3 convert_parquet_to_root.py [/path/to/merged.parquet] [path to output file containing also the filename] mc (or data depending what you're doing) --process [process name (should match one of the outfiles dict entries)] --do_syst --cats [cat_dict] --vars [variation.json]
```

`--do_syst` is not mandatory, but if it's there also the dictionary containing the variations must be specified with the `--var` option. As before the script works also without the `--cats` option, it creates a dummy selection of Pt > -1 and call the category `NOTAG`.

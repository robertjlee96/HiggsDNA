Workspace preparation for FinalFit interface
============================================

Standard Procedure
------------------

The standard way to get HiggsDNA Ntuples and transform them in FinalFit friendly output is to use the ``prepare_output_file.py`` script, provided and maintained in the ``script`` repository.
The script will perform multiple steps:
* Merge all the ``.parquet`` files and categorise the events, obtaining one file for each category of each sample.
* Convert the ``merged.parquet`` into ``ROOT`` trees.
* Convert the ``ROOT`` trees into FinalFit compatible ``RooWorkspace``s.

All the steps can be performed in one go with a command more or less like this::

        python3 prepare_output_file.py --input [path to output dir] --merge --root --ws --syst --cats --args "--do_syst"

or the single steps can be performed by running the auxiliary files (``merge_parquet.py``, ``convert_parquet_to_root.py``, ``Tree2WS``) separately.
A complete set of options for the main script is listed below.

Merging step
------------
During this step the main script calls ``merge_parquet.py`` multiple times. The starting point is the output of HiggsDNA, i.e. ``out_dir/sample_n/``. These directory **must** contain only ``.parquet`` files that have to be merged. 
The script will create a new directory called ``merged`` under ``out_dir``, if this directory already exists it will throw an error and exit.
When converting the data (in my case they were split per era, ``Data_B_2017``, ``Data_C_2017`` etc.) the script will put them in a new directory ``Data_2017`` and then merge again the output in a ``.parquet`` called ``allData_2017.parquet``.
During this step the events are also split into categories according to the boundaries defined in the ``cat_dict`` in the main file. An example of such dictionary is presented here::

        if opt.cats:
        cat_dict = {
            "best_resolution": {
                "cat_filter": [
                    ("sigma_m_over_m_decorr", "<", 0.005),
                    ("lead_mvaID", ">", 0.43),
                    ("sublead_mvaID", ">", 0.43),
                ]
            },
            "medium_resolution": {
                "cat_filter": [
                    ("sigma_m_over_m_decorr", ">", 0.005),
                    ("sigma_m_over_m_decorr", "<", 0.008),
                    ("lead_mvaID", ">", 0.43),
                    ("sublead_mvaID", ">", 0.43),
                ]
            },
            "worst_resolution": {
                "cat_filter": [
                    ("sigma_m_over_m_decorr", ">", 0.008),
                    ("lead_mvaID", ">", 0.43),
                    ("sublead_mvaID", ">", 0.43),
                ]
            },
        }

if you don't provide the dictionary to the script all the events will be put in a single file labelled as ``UNTAGGED``.

During the merging step MC samples can also be normalised to the ``efficiency x acceptance`` value as required later on by FinalFits, this step can be skipped using the tag ``--skip-normalisation``.

Root step 
---------

During this step the script calls multiple times the script ``convert_parquet_to_root.py``. The arguments to pass to the script, for instance if you want the systematic variation included in the output ``ROOT tree`` are specified when calling ``prepare_output_file.py`` using ``--args "--do_syst"``.
As before the script creates a new called ``root`` under ``out_dir``, if this directory already exists it will throw an error and exit. In the script there is a dictionary called ``outfiles`` that contains the name of the output root file that will be created according to the process tipe, if the wf is run using the main script this correspond to the proces containd in ``process_dict``.

Workspace step
--------------

During this step the main script uses multiple time the ``Flashgg_FinalFit``, it moves to the directory defined in the ``--final_fit`` option (improvable) and uses the ``Tree2WS`` script there on the content of the ``root`` directory previously created. The output is stored in ``out_dir/root/smaple_name/ws/``.

Commands
--------

The workflow is meant to be run in one go using the ``prepare_output_file.py`` script, it can be also split in different steps or run with the single auxiliary files but it can result a bit cumbersome.

To run everything starting from the output of HiggsDNA with categories and systematic variatrion one can use::

        python3 prepare_output_file.py --input [path to output dir] --merge --root --ws --syst --cats --args "--do_syst"

and everithing should run smoothly, it does for me at least (I've not tried the scripts in a while so thing may have to be adjusted in this document).
Some options can be removed. If you want to use ``--syst`` and ``--root`` you should also add ``--args "--do_syst"``.

The complete list of options for the main file is here:

    * ``--merge``, "Do merging of the .parquet files"
    * ``--root``, "Do root conversion step"
    * ``--ws``, "Do root to workspace conversion step"
    * ``--ws_config``, "configuration file for Tree2WS, as it is now it must be stored in Tree2WS directory in FinalFit",
    * ``--final_fit``, "FlashggFinalFit path" # the default is just for me, it should be changed but I don't see a way to make this generally valid
    * ``--syst``, "Do systematics variation treatment"
    * ``--cats``, ="Split into categories",
    * ``--args``, "additional options for root converter: --do_syst, --notag",
    * ``--skip-normalisation``, "Independent of file type, skip normalisation step",
    * ``--verbose``, "verbose lefer for the logger: INFO (default), DEBUG",

The merging step can also be run separately using::

        python3 merge_parquet.py --source [path to the directory containing .paruets] --target [target directory path] --cats [cat_dict]

the script works also without the ``--cats`` option, it creates a dummy selection of ``Pt > -1`` and call the category ``UNTAGGED``.

Same for the root step::

        python3 convert_parquet_to_root.py [/path/to/merged.parquet] [path to output file containing also the filename] mc (or data depending what you're doing) --process [process name (should match one of the outfiles dict entries)] --do_syst --cats [cat_dict] --vars [variation.json]

``--do_syst`` is not mandatory, but if it's there also the dictionary containing the variations must be specified with the ``--var`` option. As before the script works also without the ``--cats`` option.



from typing import List, Optional

import awkward
import pandas
import os
import pathlib
import shutil
import pyarrow.parquet as pq


def diphoton_list_to_pandas(self, diphotons: awkward.Array) -> pandas.DataFrame:
    """
    Convert diphoton array to pandas dataframe.
    By default the observables related to each item of the diphoton pair are
    stored preceded by its prefix (e.g. 'lead', 'sublead').
    The observables related to the diphoton pair are stored with no prefix.
    To change the behavior, you can redefine the `diphoton_list_to_pandas` method in the
    derived class.
    """
    output = pandas.DataFrame()
    for field in awkward.fields(diphotons):
        prefix = self.prefixes.get(field, "")
        if len(prefix) > 0:
            for subfield in awkward.fields(diphotons[field]):
                if subfield != "__systematics__":
                    output[f"{prefix}_{subfield}"] = awkward.to_numpy(
                        diphotons[field][subfield]
                    )
        else:
            output[field] = awkward.to_numpy(diphotons[field])
    return output


def dump_pandas(
    self,
    pddf: pandas.DataFrame,
    fname: str,
    location: str,
    subdirs: Optional[List[str]] = None,
) -> None:
    """
    Dump a pandas dataframe to disk at location/'/'.join(subdirs)/fname.
    """
    subdirs = subdirs or []
    xrd_prefix = "root://"
    pfx_len = len(xrd_prefix)
    xrootd = False
    if xrd_prefix in location:
        try:
            import XRootD  # type: ignore
            import XRootD.client  # type: ignore

            xrootd = True
        except ImportError as err:
            raise ImportError(
                "Install XRootD python bindings with: conda install -c conda-forge xroot"
            ) from err
    local_file = (
        os.path.abspath(os.path.join(".", fname))
        if xrootd
        else os.path.join(".", fname)
    )
    merged_subdirs = "/".join(subdirs) if xrootd else os.path.sep.join(subdirs)
    destination = (
        location + merged_subdirs + f"/{fname}"
        if xrootd
        else os.path.join(location, os.path.join(merged_subdirs, fname))
    )
    pddf.to_parquet(local_file)
    if xrootd:
        copyproc = XRootD.client.CopyProcess()
        copyproc.add_job(local_file, destination)
        copyproc.prepare()
        copyproc.run()
        client = XRootD.client.FileSystem(
            location[: location[pfx_len:].find("/") + pfx_len]
        )
        status = client.locate(
            destination[destination[pfx_len:].find("/") + pfx_len + 1 :],
            XRootD.client.flags.OpenFlags.READ,
        )
        assert status[0].ok
        del client
        del copyproc
    else:
        dirname = os.path.dirname(destination)
        if not os.path.exists(dirname):
            pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
        shutil.copy(local_file, destination)
        assert os.path.isfile(destination)
    pathlib.Path(local_file).unlink()


def diphoton_ak_array(self, diphotons: awkward.Array) -> awkward.Array:
    """
    Adjust the prefix.
    By default the observables related to each item of the diphoton pair are
    stored preceded by its prefix (e.g. 'lead', 'sublead').
    The observables related to the diphoton pair are stored with no prefix.
    """
    output = {}
    for field in awkward.fields(diphotons):
        prefix = self.prefixes.get(field, "")
        if len(prefix) > 0:
            for subfield in awkward.fields(diphotons[field]):
                if subfield != "__systematics__":
                    output[f"{prefix}_{subfield}"] = diphotons[field][subfield]
        else:
            output[field] = diphotons[field]
    return awkward.Array(output)


def dump_ak_array(
    self,
    akarr: awkward.Array,
    fname: str,
    location: str,
    metadata: None,
    subdirs: Optional[List[str]] = None,
) -> None:
    """
    Dump an awkward array to disk at location/'/'.join(subdirs)/fname.
    """
    subdirs = subdirs or []
    xrd_prefix = "root://"
    pfx_len = len(xrd_prefix)
    xrootd = False
    if xrd_prefix in location:
        try:
            import XRootD  # type: ignore
            import XRootD.client  # type: ignore

            xrootd = True
        except ImportError as err:
            raise ImportError(
                "Install XRootD python bindings with: conda install -c conda-forge xroot"
            ) from err
    local_file = (
        os.path.abspath(os.path.join(".", fname))
        if xrootd
        else os.path.join(".", fname)
    )
    merged_subdirs = "/".join(subdirs) if xrootd else os.path.sep.join(subdirs)
    destination = (
        location + merged_subdirs + f"/{fname}"
        if xrootd
        else os.path.join(location, os.path.join(merged_subdirs, fname))
    )

    pa_table = awkward.to_arrow_table(akarr)
    # If metadata is not None then write to pyarrow table
    if metadata:
        merged_metadata = {**metadata, **(pa_table.schema.metadata or {})}
        pa_table = pa_table.replace_schema_metadata(merged_metadata)

    # Write pyarrow table to parquet file
    pq.write_table(pa_table, local_file)

    if xrootd:
        copyproc = XRootD.client.CopyProcess()
        copyproc.add_job(local_file, destination)
        copyproc.prepare()
        copyproc.run()
        client = XRootD.client.FileSystem(
            location[: location[pfx_len:].find("/") + pfx_len]
        )
        status = client.locate(
            destination[destination[pfx_len:].find("/") + pfx_len + 1 :],
            XRootD.client.flags.OpenFlags.READ,
        )
        assert status[0].ok
        del client
        del copyproc
    else:
        dirname = os.path.dirname(destination)
        pathlib.Path(dirname).mkdir(parents=True, exist_ok=True)
        shutil.copy(local_file, destination)
        assert os.path.isfile(destination)
    pathlib.Path(local_file).unlink()

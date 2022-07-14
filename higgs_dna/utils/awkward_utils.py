import awkward
import vector
import numpy

import logging

logger = logging.getLogger(__name__)


def missing_fields(array, fields):
    """
    Check if every field of fields is present in a given array.
    Fields of records (e.g. Photon.pt) should be passed as tuples: ("Photon", "pt")
    :param array: array to be checked
    :type array: awkward array
    :param fields: list of fields
    :type fields: list of str or tuple
    :return: list of missing fields in array
    :rtype: list
    """

    missing_fields = []

    for field in fields:
        if isinstance(field, str):
            if field not in array.fields:
                missing_fields.append(field)
        elif isinstance(field, list):
            field = tuple(field)
        elif isinstance(field, tuple):
            sub_array = array
            for sub_field in field:
                if sub_field not in sub_array.fields:
                    missing_fields.append(field)
                    break
                sub_array = sub_array[sub_field]

        else:
            message = (
                "Each entry in the <fields> argument should be either a str or tuple, not %s which is the type of entry %s"
                % (str(type(field)), str(field))
            )
            logger.exception(message)
            raise TypeError(message)

    return missing_fields


def construct_jagged_array(offsets, contents):
    """
    Use numpy arrays of the offsets and contents to create a jagged awkward array.
    See: https://awkward-array.readthedocs.io/en/latest/ak.layout.ListOffsetArray.html

    :param offsets: specifies the starting and stopping index of each list
    :type offsets: numpy.array
    :param contents: underlying data for all lists
    :type contents: numpy.array
    :return: awkward array with unequal-length lists
    :rtype: awkward.layout.ListOffsetArray
    """
    array = awkward.Array(
        awkward.layout.ListOffsetArray64(
            awkward.layout.Index64(offsets), awkward.layout.NumpyArray(contents)
        )
    )
    return array


def add_field(events, name, data, overwrite=False):
    """
    Add a field or record to an awkward array.
    Checks if the field is already present and returns the existing field if so,
    unless the overwrite option is selected (be careful using this option!).

    :param events: base events array which you want to add a field to
    :type events: awkward.highlevel.Array
    :param name: name of the field or record you want to add
    :type name: str or tuple
    :param data: either a dictionary of subfields : awkward array (for creating a nested record) or a single awkward array (for creating a single field)
    :type data: dict or awkward.highlevel.Array or numpy.ndarray
    :param overwrite: whether to overwrite this field in events (only applicable if it already exists)
    :type overwrite: bool
    :return: events array with field or record added
    :rtype: awkward.highlevel.Array
    """

    already_present = len(missing_fields(events, [name])) == 0
    if already_present and not overwrite:
        logger.warning(
            "[awkward_utils.py : add_field] You tried to write the field %s, but it is already present and overwrite option was not selected. Not overwriting existing field."
            % name
        )
        return events[name]

    if isinstance(data, awkward.highlevel.Array) or isinstance(data, numpy.ndarray):
        events[name] = data
        return events[name]
    elif isinstance(data, dict):
        return create_record(events, name, data, overwrite)
    else:
        message = (
            "[awkward_utils.py : add_field] argument <data> should be either an awkward.highlevel.Array (in the case of adding a single field) or a dictionary (in the case of creating a record), not %s as you have passed."
            % (str(type(data)))
        )
        logger.exception(message)
        raise TypeError(message)


def add_object_fields(
    events, name, objects, n_objects, dummy_value=-999.0, fields="all", overwrite=False
):
    """
    For a collection of jagged-length objects (e.g. jets or leptons),
    add fixed-length flat fields to the events array, storing information for each of the
    first n_objects objects in each event, and filling missing events with dummy_value.

    For example add_object_fields(events, "jet", selected_jets, 4) will add every field belonging to the
    selected_jets record as individual entries for each of the first four jets:
    events.jet_1_pt
    events.jet_1_eta
    ...
    events.jet_4_phi
    where events with less than 4 jets will receive values of the dummy_value, -999.

    :param events: base events array which you want to add a field to
    :type events: awkward.highlevel.Array
    :param name: prefix to give the flattended arrays made from the original jagged record
    :type name: str
    :param objects: jagged array or record
    :type objects: awkward.highlevel.Array
    :param n_objects: number of objects to store for each event
    :type n_objects: int
    :param dummy_value: dummy value to give for events with less than n_objects objects
    :type dummy_value: float, defaults to -999
    :param fields: which fields in the <objects> record to add to the events array
    :type fields: str, list, defaults to "all"
    :param overwrite: whether to overwrite this field in events (only applicable if it already exists)
    :type overwrite: bool
    """

    padded_objects = awkward.pad_none(objects, n_objects, clip=True)
    if isinstance(fields, str):
        if fields == "all":
            fields = objects.fields
    elif not isinstance(fields, list):
        message = (
            "[awkward_utils.py : add_object_fields] argument <fields> should either be a string 'all' to save all fields in the original record, or a list of fields which is a subset of the fields in the original record, not '%s' as you have passed."
            % (str(type(fields)))
        )
        logger.exception(message)
        raise TypeError(message)

    for field in fields:
        for i in range(n_objects):
            add_field(
                events=events,
                name="%s_%d_%s" % (name, i + 1, field),
                data=awkward.fill_none(padded_objects[field][:, i], dummy_value),
                overwrite=overwrite,
            )


def create_record(events, name, data, overwrite):
    """
    Adds a nested record to events array.
    If there are any fields containing "p4", this is assumed to be an array of <vector> objects and the pt/eta/phi/mass will be added as top-level fields in the record for easier access.

    :param name: name of the record to be added
    :type name: str or tuple
    :param data: dictionary of subfields : awkward array
    :type data: dict
    :param overwrite: whether to overwrite this field in events (only applicable if it already exists)
    :type overwrite: bool
    :return: events array with record added
    :rtype: awkward.highlevel.Array
    """

    # If there is a field named "p4", save its properties as additional fields for easier access
    additional_fields = {}
    for key, array in data.items():
        if "p4" in key:
            logger.debug(
                "[awkward_utils.py : create_record] Found a field %s in your data which looks like a four vector. For convenience, we will make pt/eta/phi/mass accessible as e.g. events.%s.pt (in addition to events.%s.%s.pt)"
                % (key, name, name, key)
            )
            additional_fields[key.replace("p4", "pt")] = array.pt
            additional_fields[key.replace("p4", "eta")] = array.eta
            additional_fields[key.replace("p4", "phi")] = array.phi
            additional_fields[key.replace("p4", "mass")] = array.mass

    for key, array in additional_fields.items():
        data[key] = array

    events[name] = awkward.zip(data)

    return events[name]


def create_four_vectors(events, offsets, contents):
    """
    Zip four vectors as record in awkard array.
    Record these as `Momentum4D` vector objects so all typical
    four vector properties and methods can be called on them.
    The contents must be given as a numpy array with 4 indicies in the last
    dimension, in the order: pt, eta, phi, mass

    :param events: awkward array of events
    :type events: awkward array
    :param offsets: offsets describing the object locations
    :type offsets: numpy array
    :param contents: pt, eta, phi, mass of the object to be zipped
    :type contents: numpy array
    :return: awkward array of the four momentums
    :rtype: awkward array of Momentum4D vector
    """

    # First convert from awkward.layout.Index & awkward.layout.Array -> to awkward array
    objects = {}
    for idx, field in enumerate(["pt", "eta", "phi", "mass"]):
        objects[field] = construct_jagged_array(offsets, contents[:, idx])

    # Second, convert to awkward array of Momentum4D vector objects
    objects_p4 = vector.awk(
        {
            "pt": objects["pt"],
            "eta": objects["eta"],
            "phi": objects["phi"],
            "mass": objects["mass"],
        },
        with_name="Momentum4D",
    )

    return objects_p4

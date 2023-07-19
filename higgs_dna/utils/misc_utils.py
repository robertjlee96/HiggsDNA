import awkward as ak


def choose_jet(jets_variable, n, fill_value):
    """
    this helper function is used to create flat jets from a jagged collection,
    parameters:
    * jet_variable: (ak array) selected variable from the jet collection
    * n: (int) nth jet to be selected
    * fill_value: (float) value with wich to fill the padded none.
    """
    leading_jets_variable = jets_variable[
        ak.local_index(jets_variable) == n
    ]
    leading_jets_variable = ak.pad_none(
        leading_jets_variable, 1
    )
    leading_jets_variable = ak.flatten(
        ak.fill_none(leading_jets_variable, fill_value)
    )
    return leading_jets_variable

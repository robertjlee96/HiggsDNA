import os
import awkward
import pandas
import numpy as np
import higgs_dna.tools.decorrelator as decorr


def decorrelate_mass_resolution(events: awkward.Array):

    # reading the CDFs files
    var,dVar,ref = "sigma_m_over_m","mass",125.0   # Varable to be decorrelated, decorrelated w.r.t and the reference bin
    dummyDf = pandas.DataFrame({'{}'.format(var): [0], '{}'.format(dVar): [0]})
    decl = decorr.decorrelator(dummyDf, var, dVar, np.linspace(100., 180., 161))
    decl.loadCdfs(os.path.dirname(__file__) + '/Smeared_Diphoton_CDFs.pkl.gz')

    # setting up the decorrelator
    df = pandas.DataFrame()
    df["sigma_m_over_m"] = events.sigma_m_over_m_Smeared.to_numpy()  # Reading directly the smeared sigma_m_over_m
    df["mass"] = events.mass.to_numpy()
    df["weight"] = events.weight.to_numpy()

    decl.df = df.loc[:, [var, dVar]]
    decl.df.reset_index(inplace=True)

    # options.ref is the mass bin (125.)
    df['{}_decorr'.format(var)] = decl.doDecorr(ref)

    # performing the decorrelation
    events["sigma_m_over_m_decorr"] = decl.doDecorr(ref)

    # returning the array with the decorrelated mass resolution
    return events["sigma_m_over_m_decorr"]

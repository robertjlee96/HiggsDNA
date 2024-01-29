import os
import awkward
import pandas
import numpy as np
import higgs_dna.tools.decorrelator as decorr


def decorrelate_mass_resolution(events: awkward.Array, type: str, year):

    # type = "nominal","smeared","corr","corr_smeared"

    # reading the CDFs files
    var,dVar,ref = "sigma_m_over_m","mass",125.0   # Varable to be decorrelated, decorrelated w.r.t and the reference bin
    dummyDf = pandas.DataFrame({'{}'.format(var): [0], '{}'.format(dVar): [0]})
    decl = decorr.decorrelator(dummyDf, var, dVar, np.linspace(100., 180., 161))

    # setting up the decorrelator
    df = pandas.DataFrame()

    if (type == "nominal"):

        if (year == "2022postEE"):
            decl.loadCdfs(os.path.dirname(__file__) + '/decorrelation_CDFs/postEE/nominal_sigma_m_postEE_CDFs.pkl.gz')
        else:
            decl.loadCdfs(os.path.dirname(__file__) + '/decorrelation_CDFs/preEE/nominal_sigma_m_preEE_CDFs.pkl.gz')
        df["sigma_m_over_m"] = events.sigma_m_over_m.to_numpy()

    elif (type == "smeared"):

        if (year == "2022postEE"):
            decl.loadCdfs(os.path.dirname(__file__) + '/decorrelation_CDFs/postEE/sigma_m_smeared_postEE_CDFs.pkl.gz')
        else:
            decl.loadCdfs(os.path.dirname(__file__) + '/decorrelation_CDFs/preEE/sigma_m_smeared_preEE_CDFs.pkl.gz')
        df["sigma_m_over_m"] = events.sigma_m_over_m_Smeared.to_numpy()

    elif (type == "corr"):

        if (year == "2022postEE"):
            decl.loadCdfs(os.path.dirname(__file__) + '/decorrelation_CDFs/postEE/sigma_m_corr_postEE_CDFs.pkl.gz')
        else:
            decl.loadCdfs(os.path.dirname(__file__) + '/decorrelation_CDFs/preEE/sigma_m_corr_preEE_CDFs.pkl.gz')
        df["sigma_m_over_m"] = events.sigma_m_over_m_corr.to_numpy()

    elif (type == "corr_smeared"):

        if (year == "2022postEE"):
            decl.loadCdfs(os.path.dirname(__file__) + '/decorrelation_CDFs/postEE/sigma_m_smeared_corr_postEE_CDFs.pkl.gz')
        else:
            decl.loadCdfs(os.path.dirname(__file__) + '/decorrelation_CDFs/preEE/sigma_m_smeared_corr_preEE_CDFs.pkl.gz')
        df["sigma_m_over_m"] = events.sigma_m_over_m_Smeared_corr.to_numpy()

    else:
        print("Specify a valid type: nominal,smeared,corr,corr_smeared")
        exit()

    # Reading directly the smeared sigma_m_over_m
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

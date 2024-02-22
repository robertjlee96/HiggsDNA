import decorrelator as decorr
import awkward as ak
import argparse
import os
import pandas as pd
import numpy as np
import glob


def printProgressBar(iteration,total,prefix='',suffix='',decimals=1,length=100,fill=chr(9608),printEnd="\r"):

    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)
    if iteration == total:
        print()


def diphoton_ak_array(diphotons: ak.Array) -> ak.Array:

    output = {}
    for field in ak.fields(diphotons):
        output[field] = diphotons[field]
    return ak.Array(output)


def getArrayBranchName(branchname, fieldname, index):
    if index != ():
        return '{}{}'.format(branchname, index[0])
    return '{}'.format(branchname)


def main(options):

    if options.infile == options.outFile:
        raise RuntimeError('Outfile will be recreated, cannot be the same as infile')

    if os.path.exists(options.outFile):
        print("WARNING: outfile exists.")

    dummyDf = pd.DataFrame({'{}'.format(options.var): [0], '{}'.format(options.dVar): [0]})
    decl = decorr.decorrelator(dummyDf, options.var, options.dVar, np.linspace(100., 180., 161))
    decl.loadCdfs(options.cdfFile)

    files = glob.glob(str(options.infile) + "*.parquet")
    data = [pd.read_parquet(f) for f in files]
    # data = pd.read_parquet("/net/scratch_cms3a/daumann/massresdecorrhiggsdna/big_bkg/Diphoton.parquet")
    events = pd.concat(data,ignore_index=True)

    df = pd.DataFrame()
    df["sigma_m_over_m"] = events.sigma_m_over_m_Smeared.to_numpy()
    df["mass"] = events.mass.to_numpy()
    df["weight"] = events.weight.to_numpy()

    print("var, dVar:", options.var, options.dVar)
    decl.df = df.loc[:, [options.var, options.dVar]]
    decl.df.reset_index(inplace=True)

    df['{}_decorr'.format(options.var)] = decl.doDecorr(options.ref)

    if 'sigmaMoM_decorr' in df.columns:
        df['sigmaMoM_decorrOld'] = df['sigmaMoM_decorr']

    if options.var == 'sigmarv':
        df['sigmaMoM_decorr'] = df['sigmarv_decorr']

    if options.var == 'sigmaRV':
        df['sigmaMoM_decorr'] = df['sigmaRV_decorr']

    events["sigma_m_over_m_decorr"] = decl.doDecorr(options.ref)

    events.to_parquet(options.outFile)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    requiredArgs = parser.add_argument_group('Required Arguements')
    requiredArgs.add_argument('-t','--tree', nargs='+', required=True)
    requiredArgs.add_argument('-i', '--infile', action='store', type=str, required=True)
    requiredArgs.add_argument('-c','--cdfFile', action='store', type=str, required=True)
    requiredArgs.add_argument('-v','--var', action='store', type=str, required=True)
    requiredArgs.add_argument('-d','--dVar', action='store', type=str, required=True)
    requiredArgs.add_argument('-o','--outFile', action='store', type=str, required=True)
    optArgs = parser.add_argument_group('Optional Arguments')
    optArgs.add_argument('-r', '--ref', action='store', type=float, default=125.)
    optArgs.add_argument('--columns', nargs='+')
    optArgs.add_argument('--nomColumns', nargs='+')
    optArgs.add_argument('--vecColumns', nargs='+')
    options = parser.parse_args()
    main(options)

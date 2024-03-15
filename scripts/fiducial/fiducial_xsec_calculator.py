import argparse
import awkward as ak
import numpy as np

available_processes = ['GluGluH', 'VBFH', 'VH', 'ttH', 'all']
available_years = ['2022']
available_eras = ['preEE', 'postEE', 'all']

# Setup command-line argument parsing
parser = argparse.ArgumentParser(description = "Calculate the inclusive fiducial cross section of pp->H(yy)+X process(es) based on processed samples without detector-level selections.")
parser.add_argument('path', type = str, help = "Path to the top-level folder containing the different directories.")
parser.add_argument('--process', type = str, choices = available_processes, default = 'GluGluH', help = "Please specify the process(es) for which you want to calculate the inclusive fiducial xsec.")
parser.add_argument('--year', type = str, choices = available_years, default = '2022', help = 'Please specify the desired year if you want to combine samples from multiple eras.')
parser.add_argument('--era', type = str, choices = available_eras, default = 'postEE', help = 'Please specify the era(s) that you want to run over. You do not have to specify all as the eras (corresponding to data-taking conditions) should not change the fiducial cross section, but for consistency, a lumi-weighting of the samples from different eras is implemented.')
parser.add_argument('--doLumiWeighting', action = 'store_true', help = 'Specify this flag if you want to weight the respective eras with corresponding integrated luminosity during data taking. If not specified, instead do inverse variance weighting by MC stat as the proper statistical combination procedure for point estimators.')

args = parser.parse_args()

path_folder = args.path # Use the specified folder path
# Pepare the processes array appropriately
if args.process == 'all':
    processes = available_processes
    processes.remove('all')
else:
    processes = [args.process]
year = args.year
# Pepare the eras array appropriately
if args.era == 'all':
    eras = available_eras
    eras.remove('all')
else:
    eras = [args.era]
do_lumi_weighting = args.doLumiWeighting
if not do_lumi_weighting:
    variances = {}

# See also the following pages (note numbers always in picobarn)
# 13: for 125
# 13p6: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG136TeVxsec_extrap, for 125.38
# 14: for 125
XS_map = {'13':   {'GluGluH': 48.58, 'VBFH': 3.782, 'VH': 2.2569, 'ttH': 0.5071}, 
         '13p6': {'GluGluH': 51.96, 'VBFH': 4.067, 'VH': 2.3781, 'ttH': 0.5638}, 
         '14':   {'GluGluH': 54.67, 'VBFH': 4.278, 'VH': 2.4991, 'ttH': 0.6137},}

# c.f. https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#DATA_AN2
lumi_map = {'2022': {'preEE': 7.98, 'postEE': 26.67}}
lumi_map = lumi_map[year]

# Extract the total lumi based on the specified era(s)
lumi_list = [lumi_map[era] for era in eras]
total_lumi = np.sum(lumi_list)

# This depends on how you named your samples in HiggsDNA
processMap = {'GluGluH': 'GluGluHtoGG_M-125',
              'VBFH': 'VBFHtoGG_M-125',
              'VH': 'VHtoGG_M-125',
              'ttH': 'ttHtoGG_M-125',}

BR = 0.2270/100 # SM value for mH close to 125: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR

# Already saving some code snippets 
#differential_var = 'HTXS_Higgs_pt' # HTXS_Higgs_y and HTXS_njets30
#binning = np.array([0., 1., 2., 3., 4., np.infty]) # for njets
#differential_flag = arr[differential_var] > 15
#sumwIn = ak.sum(arr.weight[(inFiducialFlag & differential_flag)])

fid_xsecs = []
for era in eras:
    print(f'INFO: Now extracting numbers for era: {era} ...')
    if not do_lumi_weighting: variances_tmp = []
    for process in processes:
        print(f'INFO: Now extracting fraction of in-fiducial events for process {process} ...')
        # Extract the events
        process_string = processMap[process]
        arr = ak.from_parquet(path_folder + process_string + '_' + era + '/nominal')
        # Calculating the relevant fractions
        inFiducialFlag = arr.fiducialGeometricTagger_20 == 21 # Only for this type of tagger right now, can be customised in the future
        sumwAll = ak.sum(arr.weight)

        sumwIn = ak.sum(arr.weight[(inFiducialFlag)])
        in_frac = sumwIn/sumwAll

        print(f"INFO: Fraction of in-fiducial events: {in_frac} ...")

        if do_lumi_weighting:
            rel_fraction = lumi_map[era] / total_lumi
        else:
            sumw2 = ak.sum(arr.weight[(inFiducialFlag)]**2)
            rel_fraction = 1/sumw2 # inverse variance
            variances_tmp.append(sumw2)

        fid_xsecs.append(in_frac * XS_map['13p6'][process] * 1000 * BR * rel_fraction) # also converting to femtobarn on the way
    
    if not do_lumi_weighting:
        variances_tmp = np.asarray(variances_tmp)
        variances[era] = np.sqrt(np.sum(variances_tmp))

if do_lumi_weighting:
    result = np.sum(fid_xsecs)
else:
    variances_list = np.asarray([variances[era] for era in eras])
    result = np.sum(fid_xsecs)/np.sum(1/variances_list)

print(f"The inclusive fiducial cross section is given by: {result} fb")
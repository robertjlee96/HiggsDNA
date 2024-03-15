import argparse
import awkward as ak
import numpy as np

available_processes = ['GluGluH', 'VBFH', 'VH', 'ttH', 'all']
available_years = ['2022']
available_eras = ['preEE', 'postEE', 'all']

# Setup command-line argument parsing
parser = argparse.ArgumentParser(description = "Calculate the inclusive fiducial cross section of pp->H(yy)+X process(es) based on processed samples without detector-level selections. ")
parser.add_argument('path', type = str, help = "Path to the top-level folder containing the different directories.")
parser.add_argument('--process', type = str, choices = available_processes, default = 'GluGluH', help = "Please specify the process(es) for which you want to calculate the inclusive fiducial xsec.")
parser.add_argument('--year', type = str, choices = available_years, default = '2022', help = 'Please specify the desired year if you want to combine samples from multiple eras.')
parser.add_argument('--era', type = str, choices = available_eras, default = 'postEE', help = "Please specify the era(s) that you want to run over. If you specify 'all', an inverse variance weighting is performed to increase the precision.")

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

# See also the following pages (note numbers always in picobarn)
# 13: for 125
# 13p6: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG136TeVxsec_extrap, for 125.38
# 14: for 125
XS_map = {'13':   {'GluGluH': 48.58, 'VBFH': 3.782, 'VH': 2.2569, 'ttH': 0.5071}, 
         '13p6': {'GluGluH': 51.96, 'VBFH': 4.067, 'VH': 2.3781, 'ttH': 0.5638}, 
         '14':   {'GluGluH': 54.67, 'VBFH': 4.278, 'VH': 2.4991, 'ttH': 0.6137},}

# c.f. https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#DATA_AN2
#lumi_map = {'2022': {'preEE': 7.98, 'postEE': 26.67}}

# This depends on how you named your samples in HiggsDNA
processMap = {'GluGluH': 'GluGluHtoGG_M-125',
              'VBFH': 'VBFHtoGG_M-125',
              'VH': 'VHtoGG_M-125',
              'ttH': 'ttHtoGG_M-125',}

BR = 0.2270/100 # SM value for mH close to 125: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR

# Already saving some code snippets for differentials
#differential_var = 'HTXS_Higgs_pt' # HTXS_Higgs_y and HTXS_njets30
#binning = np.array([0., 1., 2., 3., 4., np.infty]) # for njets
#differential_flag = arr[differential_var] > 15
#sumwIn = ak.sum(arr.weight[(inFiducialFlag & differential_flag)])

fid_xsecs_per_process = {}
for process in processes:
    print(f'INFO: Now extracting fraction of in-fiducial events for process {process} ...')
    fid_xsecs_per_era = {} # These numbers are already weighted with the inverse of the MC stat variance
    sumw2_tmp = []
    for era in eras:
        print(f'INFO: Now extracting numbers for era: {era} ...')
        # Extract the events
        process_string = processMap[process]
        arr = ak.from_parquet(path_folder + process_string + '_' + era + '/nominal')
        # Calculating the relevant fractions
        inFiducialFlag = arr.fiducialGeometricTagger_20 == 21 # Only for this type of tagger right now, can be customised in the future
        sumwAll = ak.sum(arr.weight)

        sumwIn = ak.sum(arr.weight[(inFiducialFlag)])
        in_frac = sumwIn/sumwAll

        print(f"INFO: Fraction of in-fiducial events: {in_frac} ...")

        sumw2 = ak.sum(arr.weight[(inFiducialFlag)]**2) # This is the MC stat variance
        sumw2_tmp.append(sumw2)

        fid_xsecs_per_era[era] = in_frac * XS_map['13p6'][process] * 1000 * BR * 1/sumw2 # also converting to femtobarn on the way
    
    sumw2_tmp = np.asarray(sumw2_tmp)
    result = np.sum(np.asarray([fid_xsecs_per_era[era] for era in eras]))
    fid_xsecs_per_process[process] = result / np.sum(1/sumw2_tmp)

final_fid_xsec = np.sum(np.asarray([fid_xsecs_per_process[process] for process in processes]))

print(f"The inclusive fiducial cross section is given by: {final_fid_xsec} fb")
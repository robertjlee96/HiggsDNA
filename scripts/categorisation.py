"""
Script for optimization of event categories based on boundaries on decorrelated mass resolution estimator and photon MVA ID cut.
The optimization is done by fitting signal and background histograms with appropriate functions, determining S/sqrt(B), and comparing the sensitivities achieved with different category partitionings.
This code can be run on 1 to 4 categories in mass resolution.

Usage:
    1. Run HiggsDNA on Diphoton, GJet and Hgg signal samples.
    2. Set variable `base_path` to point to the output of HiggsDNA.
    3. Define resolution and MVA ID boundaries in 'resboundaries' and 'MVAboundaries'. NB: the script will take a long time to run over 4 categories, so it is helpful to first run it on a coarse grid and then on a fine grid around the optimum.
    4. Define over how many resolution categories to run, e.g. cats = [1, 3] to test one vs. 3 resolution categories.
    5. Run the script.
"""

import numpy as np
import hist
import os
import glob
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import quad
import pyarrow.parquet as pq
import matplotlib.pyplot as plt
import mplhep as hep
import time
hep.style.use("CMS")


def exponential(x, a, b, c):
    return a * np.exp(- x / b) + c


def double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
    g1 = A1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2))
    g2 = A2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2))
    return g1 + g2


def get_central_interval(y_vals, dx, area=0.68):
    sorted_indices = np.argsort(y_vals)[::-1]  # Sort in descending order
    cumulative_area = 0
    central_interval = []

    for idx in sorted_indices:
        cumulative_area += y_vals[idx] * dx
        central_interval.append(x_steps[idx])
        if cumulative_area >= area:
            break
    central_interval = [min(central_interval), max(central_interval)]
    return central_interval


def plot_category(hist_sig, hists_bkg, x_vals, y_vals_sig, y_vals_bkg, params_sig, ax, hist_unc_bkg=None):

    hep.histplot(
        hists_bkg,
        label=["GJet", "Diphoton"],
        histtype="fill",
        linewidth=3,
        stack=True,
        alpha=0.5,
        color=["tab:gray","tab:blue"],
        ax=ax
    )

    if hist_unc_bkg is not None:

        errors_bkg = np.sqrt(hist_unc_bkg.to_numpy()[0])
        lower_bound = hists_bkg[0].to_numpy()[0] + hists_bkg[1].to_numpy()[0] - errors_bkg
        upper_bound = hists_bkg[0].to_numpy()[0] + hists_bkg[1].to_numpy()[0] + errors_bkg
        _, bin_edges = hists_bkg[0].to_numpy()
        bin_centers = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(bin_edges) - 1)]
        # Plot hatched region
        ax.fill_between(
            bin_centers,
            lower_bound,
            upper_bound,
            hatch='XXX',
            facecolor="none",
            edgecolor="tab:gray",
            linewidth=0
        )
    ax.plot(x_vals, y_vals_bkg, label="Background fit", color="black", linewidth=3, linestyle="--")
    hep.histplot(hist_sig, label=r"ggH+VBF+VH+ttH $(\times 10)$", histtype="fill", linewidth=3, color="tab:orange", ax=ax)
    ax.plot(x_vals, y_vals_sig, label=r"Signal fit $(\times 10)$", color="tab:red", linewidth=3, linestyle="--", alpha=0.7)
    ax.text(0.05, 0.73, r"Gauss 1: $\mu={:.2f}\,$GeV, $\sigma={:.2f}\,$GeV".format(params_sig[1], params_sig[2]), fontsize=18, transform=ax.transAxes)
    ax.text(0.05, 0.68, r"Gauss 2: $\mu={:.2f}\,$GeV, $\sigma={:.2f}\,$GeV".format(params_sig[4], params_sig[5]), fontsize=18, transform=ax.transAxes)
    ax.set_xlabel('Invariant diphoton mass [GeV]')
    ax.set_ylabel(r'Events / GeV')
    ax.set_ylim(0., 1.5 * ax.get_ylim()[1])
    ax.legend(ncol=2)
    ax.set_xlim(100.,180.)
    hep.cms.label(data=True, ax=ax, loc=0, label="Simulation Work in Progress", com=13.6, lumi=20.7, fontsize=22)


# values from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG136TeVxsec_extrap and XSDB
dict_xSecs = {
    # in fb
    "ggH": 52.23e3 * 0.00227,
    "VBF": 4.078e3 * 0.00227,
    "VH": 2.4009e3 * 0.00227,
    "ttH": 0.5700e3 * 0.00227,
    "Diphoton": 89.14e3,
    "GJetPT20to40": 242.5e3,
    "GJetPT40": 919.1e3,
}
lumi = 20.7  # /fb, 21.7 is recorded 2022F+G, 20.7 is GoldenJSON, https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis

base_path = "/path/to/my/samples/processed/with/HiggsDNA/"

dict_paths_signal = {
    "ggH": base_path + "GluGluHToGG_postEE_M125_2022/nominal/",
    "VBF": base_path + "VBFtoGG_postEE_M125_2022/nominal/",
    "VH": base_path + "VHtoGG_postEE_M125_2022/nominal/",
    "ttH": base_path + "ttHtoGG_postEE_M125_2022/nominal/"
}
dict_paths_bkg = {
    "Diphoton": base_path + "Diphoton/nominal/",
    "GJetPT20to40": base_path + "GJetPt20to40/nominal/",
    "GJetPT40": base_path + "GJetPt40toInf/nominal/"
}

signal_events = []
for process in dict_paths_signal.keys():
    files_signal = glob.glob(dict_paths_signal[process] + "/*.parquet")
    data_signal = [pd.read_parquet(f) for f in files_signal]
    events = pd.concat(data_signal,ignore_index=True)
    events = events[
        (events.mass > 100) & (events.mass < 180)
    ]
    sum_genw_beforesel = 0
    for file in files_signal:
        sum_genw_beforesel += float(pq.read_table(file).schema.metadata[b'sum_genw_presel'])
    events["weight"] *= (lumi * dict_xSecs[process] / sum_genw_beforesel)
    signal_events.append(events)
signal_events = pd.concat(signal_events,ignore_index=True)
signal_events["min_mvaID"] = np.min([signal_events.lead_mvaID.values, signal_events.sublead_mvaID.values], axis=0)

bkg_events = []
for process in dict_paths_bkg.keys():
    files_bkg = glob.glob(dict_paths_bkg[process] + "/*.parquet")
    data_bkg = [pd.read_parquet(f) for f in files_bkg]
    events = pd.concat(data_bkg, ignore_index=True)
    events = events[
        (events.mass > 100) & (events.mass < 180)
    ]
    sum_genw_beforesel = 0
    for file in files_bkg:
        sum_genw_beforesel += float(pq.read_table(file).schema.metadata[b'sum_genw_presel'])
    events["weight"] *= (lumi * dict_xSecs[process] / sum_genw_beforesel)
    events["process"] = process
    if "GJet" in process:
        print("INFO: applying overlap removal for sample", process)
        events = events[events.lead_genPartFlav + events.sublead_genPartFlav != 2]
    bkg_events.append(events)
bkg_events = pd.concat(bkg_events,ignore_index=True)
bkg_events["min_mvaID"] = np.min([bkg_events.lead_mvaID.values, bkg_events.sublead_mvaID.values], axis=0)

relevant_columns_sig = ['mass', 'weight', 'sigma_m_over_m_decorr', 'min_mvaID']
signal_events = signal_events[relevant_columns_sig]
signal_events["weight"] = 10 * signal_events["weight"].values
relevant_columns_bkg = ['mass', 'weight', 'sigma_m_over_m_decorr', 'min_mvaID', 'process']
bkg_events = bkg_events[relevant_columns_bkg]

resboundaries = np.linspace(0.005, 0.025, 21)
MVAboundaries = np.linspace(0., 0.3, 31)
relative_metrics = []
cats = [1, 2, 3, 4]  # MEOW

for n_cat in cats:

    print(f"\n INFO: starting calculation for {n_cat} categories.\n")

    if n_cat == 1:
        path_plots_1_cat = "./Plots_1cat/"
        if not os.path.exists(path_plots_1_cat):
            os.makedirs(path_plots_1_cat)
        metrics = []

        for i, MVAbound in enumerate(MVAboundaries):
            # Filter events based on MVAbound
            signal_events_cat0 = signal_events[signal_events['min_mvaID'] > MVAbound]
            bkg_events_cat0 = bkg_events[bkg_events['min_mvaID'] > MVAbound]

            ##### background fit
            hist_bkg = hist.Hist(hist.axis.Regular(160, 100, 180))
            hist_bkg.fill(bkg_events_cat0.mass.values, weight=bkg_events_cat0.weight.values)
            histo, bin_edges = hist_bkg.to_numpy()
            bin_centers = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(histo))]
            params_bkg, covariance = curve_fit(exponential, bin_centers, histo, p0=([10000, 40, 0]))

            ##### signal fit
            hist_sig = hist.Hist(hist.axis.Regular(160, 100, 180))
            hist_sig.fill(signal_events_cat0.mass.values, weight=signal_events_cat0.weight.values)
            histo, bin_edges = hist_sig.to_numpy()
            initial_guess = [5000, 124, 5, 5000, 125, 1]
            params_sig, covariance = curve_fit(double_gaussian, bin_centers, histo, p0=initial_guess)

            x_steps = np.linspace(bin_edges[0], bin_edges[-1], 10000)
            y_vals = double_gaussian(x_steps, *params_sig)

            # Normalize to form a PDF
            dx = x_steps[1] - x_steps[0]
            y_vals /= np.sum(y_vals * dx)
            central_interval = get_central_interval(y_vals, dx, area=0.68)
            signal_integral, _ = quad(double_gaussian, central_interval[0], central_interval[1], args=tuple(params_sig))
            background_integral, _ = quad(exponential, central_interval[0], central_interval[1], args=tuple(params_bkg))

            # Generating values for plotting
            y_vals_bkg = exponential(x_steps, params_bkg[0], params_bkg[1], params_bkg[2])
            y_vals_sig = double_gaussian(x_steps, params_sig[0], params_sig[1], params_sig[2], params_sig[3], params_sig[4], params_sig[5])

            metric = (signal_integral / np.sqrt(background_integral))
            metrics.append(metric)

            # Plotting
            hist_diphoton = hist.Hist(hist.axis.Regular(80, 100, 180))
            hist_diphoton.fill(bkg_events_cat0.mass[bkg_events_cat0.process == "Diphoton"].values, weight=bkg_events_cat0.weight[bkg_events_cat0.process == "Diphoton"].values)
            hist_GJet = hist.Hist(hist.axis.Regular(80, 100, 180))
            hist_GJet.fill(bkg_events_cat0.mass[bkg_events_cat0.process != "Diphoton"].values, weight=bkg_events_cat0.weight[bkg_events_cat0.process != "Diphoton"].values)
            # hist of squared weights for error bars of weighted hists. I am not aware of an easier method to get correct error bars with the hist package....
            hist_bkg_unc = hist.Hist(hist.axis.Regular(80, 100, 180))
            hist_bkg_unc.fill(bkg_events_cat0.mass.values, weight=bkg_events_cat0.weight.values**2)

            # put to figure
            fig, ax = plt.subplots(1, 1, figsize=(10 * n_cat, 10))
            plot_category(hist_sig[::2j], [hist_GJet, hist_diphoton], x_steps, y_vals_sig * 2, y_vals_bkg * 2, params_sig, ax=ax, hist_unc_bkg=hist_bkg_unc)
            fig.suptitle("MVA ID bound: {:.3f}, abs. metric: {:.4f}".format(MVAbound, metric), fontsize=24)
            fig.tight_layout()
            path = path_plots_1_cat
            if not os.path.exists(path):
                os.makedirs(path)
            fig.savefig(path + f"fit_MVA_{i}.pdf")
            plt.close()

        metric_1cat = np.max(np.array(metrics))
        relative_metrics.append(metric_1cat)
        print("Metric for 1 category:", metric_1cat)

    if 1 not in cats:
        metric_1cat = 33.055  # just for later comparison

    if n_cat == 2:
        t0 = time.time()
        metrics = []
        for i_res, resbound in enumerate(resboundaries):
            print("Resolution bound", i_res)
            _metrics = []
            params = []  # for later plotting
            for i, MVAbound in enumerate(MVAboundaries):

                events_sig_cat0 = signal_events[(signal_events['sigma_m_over_m_decorr'] < resbound) & (signal_events['min_mvaID'] > MVAbound)]
                events_sig_cat1 = signal_events[(signal_events['sigma_m_over_m_decorr'] >= resbound) & (signal_events['min_mvaID'] > MVAbound)]
                events_bkg_cat0 = bkg_events[(bkg_events['sigma_m_over_m_decorr'] < resbound) & (bkg_events['min_mvaID'] > MVAbound)]
                events_bkg_cat1 = bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= resbound) & (bkg_events['min_mvaID'] > MVAbound)]

                ##### background fit
                hist_bkg_cat0 = hist.Hist(hist.axis.Regular(160, 100, 180))
                hist_bkg_cat0.fill(events_bkg_cat0.mass.values, weight=events_bkg_cat0.weight.values)
                histo, bin_edges = hist_bkg_cat0.to_numpy()
                bin_centers = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(histo))]
                params_bkg_cat0, _ = curve_fit(exponential, bin_centers, histo, p0=([10000, 40, 0]))

                hist_bkg_cat1 = hist.Hist(hist.axis.Regular(160, 100, 180))
                hist_bkg_cat1.fill(events_bkg_cat1.mass.values, weight=events_bkg_cat1.weight.values)
                histo, bin_edges = hist_bkg_cat1.to_numpy()
                bin_centers = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(histo))]
                params_bkg_cat1, _ = curve_fit(exponential, bin_centers, histo, p0=([10000, 40, 0]))

                ##### signal fit
                hist_sig_cat0 = hist.Hist(hist.axis.Regular(160, 100, 180))
                hist_sig_cat0.fill(events_sig_cat0.mass.values, weight=events_sig_cat0.weight.values)
                histo, bin_edges = hist_sig_cat0.to_numpy()
                initial_guess = [5000, 124, 5, 5000, 125, 1]
                params_sig_cat0, _ = curve_fit(double_gaussian, bin_centers, histo, p0=initial_guess)

                hist_sig_cat1 = hist.Hist(hist.axis.Regular(160, 100, 180))
                hist_sig_cat1.fill(events_sig_cat1.mass.values, weight=events_sig_cat1.weight.values)
                histo, bin_edges = hist_sig_cat1.to_numpy()
                initial_guess = [5000, 124, 5, 5000, 125, 1]
                params_sig_cat1, _ = curve_fit(double_gaussian, bin_centers, histo, p0=initial_guess)

                params.append((params_sig_cat0, params_bkg_cat0, params_sig_cat1, params_bkg_cat1))

                # Generate a dense array of x values
                x_steps = np.linspace(bin_edges[0], bin_edges[-1], 10000)

                # Evaluate the fitted double Gaussian function at each x value
                y_vals_cat0 = double_gaussian(x_steps, *params_sig_cat0)
                y_vals_cat1 = double_gaussian(x_steps, *params_sig_cat1)

                # Normalize to form a PDF
                dx = x_steps[1] - x_steps[0]
                y_vals_cat0 /= np.sum(y_vals_cat0 * dx)
                y_vals_cat1 /= np.sum(y_vals_cat1 * dx)

                # Find the central interval
                central_interval_cat0 = get_central_interval(y_vals_cat0, dx, area=0.68)
                central_interval_cat1 = get_central_interval(y_vals_cat1, dx, area=0.68)

                signal_integral_cat0, _ = quad(double_gaussian, central_interval_cat0[0], central_interval_cat0[1], args=tuple(params_sig_cat0))
                signal_integral_cat1, _ = quad(double_gaussian, central_interval_cat1[0], central_interval_cat1[1], args=tuple(params_sig_cat1))

                background_integral_cat0, _ = quad(exponential, central_interval_cat0[0], central_interval_cat0[1], args=tuple(params_bkg_cat0))
                background_integral_cat1, _ = quad(exponential, central_interval_cat1[0], central_interval_cat1[1], args=tuple(params_bkg_cat1))

                metric = np.sqrt((signal_integral_cat0 / np.sqrt(background_integral_cat0))**2 + (signal_integral_cat1 / np.sqrt(background_integral_cat1))**2)
                _metrics.append(metric)

            # Plot best MVA category (if above threshold)
            if np.max(_metrics) / metric_1cat > 1.08:
                max_MVAbound = MVAboundaries[np.argmax(_metrics)]

                # select arrays for best MVA bound again
                events_sig_cat0 = signal_events[(signal_events['sigma_m_over_m_decorr'] < resbound) & (signal_events['min_mvaID'] > max_MVAbound)]
                events_sig_cat1 = signal_events[(signal_events['sigma_m_over_m_decorr'] >= resbound) & (signal_events['min_mvaID'] > max_MVAbound)]
                events_bkg_cat0 = bkg_events[(bkg_events['sigma_m_over_m_decorr'] < resbound) & (bkg_events['min_mvaID'] > max_MVAbound)]
                events_bkg_cat1 = bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= resbound) & (bkg_events['min_mvaID'] > max_MVAbound)]

                # histogramming
                hists_signal = [hist.Hist(hist.axis.Regular(80, 100, 180)) for _ in range(2)]
                hists_signal[0].fill(events_sig_cat0.mass.values, weight=events_sig_cat0.weight.values)
                hists_signal[1].fill(events_sig_cat1.mass.values, weight=events_sig_cat1.weight.values)

                # Assuming the process distinction is available in 'process' column in events dataframe
                hists_diphoton = [hist.Hist(hist.axis.Regular(80, 100, 180)) for _ in range(2)]
                hists_diphoton[0].fill(events_bkg_cat0.mass[events_bkg_cat0.process == "Diphoton"].values, weight=events_bkg_cat0.weight[events_bkg_cat0.process == "Diphoton"].values)
                hists_diphoton[1].fill(events_bkg_cat1.mass[events_bkg_cat1.process == "Diphoton"].values, weight=events_bkg_cat1.weight[events_bkg_cat1.process == "Diphoton"].values)

                hists_GJet = [hist.Hist(hist.axis.Regular(80, 100, 180)) for _ in range(2)]
                hists_GJet[0].fill(events_bkg_cat0.mass[events_bkg_cat0.process != "Diphoton"].values, weight=events_bkg_cat0.weight[events_bkg_cat0.process != "Diphoton"].values)
                hists_GJet[1].fill(events_bkg_cat1.mass[events_bkg_cat1.process != "Diphoton"].values, weight=events_bkg_cat1.weight[events_bkg_cat1.process != "Diphoton"].values)

                # extract fitted values
                params_sig_cat0, params_bkg_cat0, params_sig_cat1, params_bkg_cat1 = params[np.argmax(_metrics)]

                # generate x_steps based on the previously defined range
                x_steps = np.linspace(100, 180, 10000)

                y_vals_bkg_cat0 = exponential(x_steps, *params_bkg_cat0)
                y_vals_sig_cat0 = double_gaussian(x_steps, *params_sig_cat0)
                y_vals_bkg_cat1 = exponential(x_steps, *params_bkg_cat1)
                y_vals_sig_cat1 = double_gaussian(x_steps, *params_sig_cat1)

                # put to figure
                fig, axes = plt.subplots(1, n_cat, figsize=(10 * n_cat, 10))
                plot_category(hists_signal[0], [hists_GJet[0], hists_diphoton[0]], x_steps, y_vals_sig_cat0 * 2, y_vals_bkg_cat0 * 2, params_sig_cat0, ax=axes[0])
                plot_category(hists_signal[1], [hists_GJet[1], hists_diphoton[1]], x_steps, y_vals_sig_cat1 * 2, y_vals_bkg_cat1 * 2, params_sig_cat1, ax=axes[1])
                fig.suptitle("MVA ID bound: {:.3f}, resolution bound: {:.3f}, metric: {:.4f}".format(max_MVAbound, resbound, np.max(_metrics) / metric_1cat), fontsize=24)
                fig.tight_layout()
                path = "./Plots_2cats/"
                if not os.path.exists(path):
                    os.makedirs(path)
                fig.savefig(path + f"ResBound_{str(i_res)}.pdf")
                plt.close()

            metrics.append(np.max(_metrics))

        relative_metrics.append(np.max(metrics))

        fig, ax = plt.subplots()
        plt.plot(resboundaries, np.array(metrics) / metric_1cat)
        plt.xlabel("Boundary in sigma m / m")
        plt.ylabel("Relative sensitivity")
        plt.savefig("./sensitivity_2_cats.pdf")
        print("Metric 2 categories:", np.max(metrics))
        t1 = time.time()
        print(f"time for the whole 2 cat scan: {t1-t0:.2f}s.")

    if n_cat == 3:
        t0 = time.time()
        metrics = np.zeros((len(resboundaries), len(resboundaries)))
        metrics_df = pd.DataFrame(columns=["resBound1", "resBound2", "MVA_bound", "relativeMetric"])

        for i, bound1 in enumerate(resboundaries):
            print("\n bound:", i)
            if bound1 > 0.012:
                print("WARNING: skipping cats with bound 1 > 0.012 for time reasons. Make sure that this fits your values!")  # these are worse anyways with the values we have here

            for j, bound2 in enumerate(resboundaries):
                _metrics = []
                params = []
                for iMVA, MVAbound in enumerate(MVAboundaries):

                    events_sig_cats = [
                        signal_events[(signal_events['sigma_m_over_m_decorr'] < bound1) & (signal_events['min_mvaID'] > MVAbound)],
                        signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound1) & (signal_events['sigma_m_over_m_decorr'] < bound2) & (signal_events['min_mvaID'] > MVAbound)],
                        signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound2) & (signal_events['min_mvaID'] > MVAbound)]
                    ]

                    events_bkg_cats = [
                        bkg_events[(bkg_events['sigma_m_over_m_decorr'] < bound1) & (bkg_events['min_mvaID'] > MVAbound)],
                        bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound1) & (bkg_events['sigma_m_over_m_decorr'] < bound2) & (bkg_events['min_mvaID'] > MVAbound)],
                        bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound2) & (bkg_events['min_mvaID'] > MVAbound)]
                    ]

                    _params = []  # for later plotting
                    metric_val = 0

                    for k in range(3):
                        hist_bkg = hist.Hist(hist.axis.Regular(160, 100, 180))
                        hist_bkg.fill(events_bkg_cats[k].mass.values, weight=events_bkg_cats[k].weight.values)
                        histo, bin_edges = hist_bkg.to_numpy()
                        bin_centers = [(bin_edges[_l] + bin_edges[_l + 1]) / 2 for _l in range(len(histo))]
                        params_bkg, _ = curve_fit(exponential, bin_centers, histo, p0=([10000, 40, 0]))

                        hist_sig = hist.Hist(hist.axis.Regular(160, 100, 180))
                        hist_sig.fill(events_sig_cats[k].mass.values, weight=events_sig_cats[k].weight.values)
                        histo, bin_edges = hist_sig.to_numpy()
                        initial_guess = [5000, 124, 5, 5000, 125, 1]

                        try:
                            params_sig, _ = curve_fit(double_gaussian, bin_centers, histo, p0=initial_guess)
                        except RuntimeError:
                            print(f"WARNING: signal fit in cat. with res. bounds [{bound1:.3f},{bound2:.3f}] and MVA bound {MVAbound:.3f} failed.")
                            params_sig = [0., 124, 999, 0., 125, 999]

                        _params.append(params_sig)
                        _params.append(params_bkg)

                        x_steps = np.linspace(bin_edges[0], bin_edges[-1], 10000)
                        y_vals = double_gaussian(x_steps, *params_sig)
                        dx = x_steps[1] - x_steps[0]
                        y_vals /= np.sum(y_vals * dx)

                        central_interval = get_central_interval(y_vals, dx, area=0.68)
                        signal_integral, _ = quad(double_gaussian, central_interval[0], central_interval[1], args=tuple(params_sig))
                        background_integral, _ = quad(exponential, central_interval[0], central_interval[1], args=tuple(params_bkg))

                        metric_val += (signal_integral / np.sqrt(background_integral)) ** 2

                    metric = np.sqrt(metric_val)
                    _metrics.append(metric)
                    params.append(_params)

                metrics[i, j] = np.max(_metrics)

                # Plot best MVA category for the current resolution category
                if np.max(_metrics) / metric_1cat > 1.10:
                    metrics_df.loc[len(metrics_df)] = [bound1, bound2, MVAboundaries[np.argmax(_metrics)], np.max(_metrics) / metric_1cat]
                    max_MVAbound = MVAboundaries[np.argmax(_metrics)]

                    # Select arrays for best MVA bound again
                    events_sig_cats = [
                        signal_events[(signal_events['sigma_m_over_m_decorr'] < bound1) & (signal_events['min_mvaID'] > max_MVAbound)],
                        signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound1) & (signal_events['sigma_m_over_m_decorr'] < bound2) & (signal_events['min_mvaID'] > max_MVAbound)],
                        signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound2) & (signal_events['min_mvaID'] > max_MVAbound)]
                    ]

                    events_bkg_cats = [
                        bkg_events[(bkg_events['sigma_m_over_m_decorr'] < bound1) & (bkg_events['min_mvaID'] > max_MVAbound)],
                        bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound1) & (bkg_events['sigma_m_over_m_decorr'] < bound2) & (bkg_events['min_mvaID'] > max_MVAbound)],
                        bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound2) & (bkg_events['min_mvaID'] > max_MVAbound)]
                    ]

                    # Histogramming
                    hists_signal = [hist.Hist(hist.axis.Regular(80, 100, 180)) for _ in range(3)]
                    hists_diphoton = [hist.Hist(hist.axis.Regular(80, 100, 180)) for _ in range(3)]
                    hists_GJet = [hist.Hist(hist.axis.Regular(80, 100, 180)) for _ in range(3)]
                    # hist of squared weights for error bars of weighted hists. I am not aware of an easier method to get correct error bars with the hist package....
                    hists_bkg_unc = [hist.Hist(hist.axis.Regular(80, 100, 180)) for _ in range(3)]
                    for i_, (event_sig, event_bkg) in enumerate(zip(events_sig_cats, events_bkg_cats)):
                        hists_signal[i_].fill(event_sig.mass.values, weight=event_sig.weight.values)
                        hists_diphoton[i_].fill(event_bkg.mass[event_bkg.process == "Diphoton"].values, weight=event_bkg.weight[event_bkg.process == "Diphoton"].values)
                        hists_GJet[i_].fill(event_bkg.mass[event_bkg.process != "Diphoton"].values, weight=event_bkg.weight[event_bkg.process != "Diphoton"].values)
                        hists_bkg_unc[i_].fill(event_bkg.mass.values, weight=event_bkg.weight.values**2)

                    # Extract fitted values
                    params_list = params[np.argmax(_metrics)]
                    y_vals_bkg = [exponential(x_steps, *params_list[_i * 2 + 1]) for _i in range(3)]
                    y_vals_sig = [double_gaussian(x_steps, *params_list[_i * 2]) for _i in range(3)]

                    # Put to figure
                    fig, axes = plt.subplots(1, 3, figsize=(30, 10))
                    for _i in range(3):
                        plot_category(hists_signal[_i], [hists_GJet[_i], hists_diphoton[_i]], x_steps, y_vals_sig[_i] * 2, y_vals_bkg[_i] * 2, params_list[_i * 2], ax=axes[_i], hist_unc_bkg=hists_bkg_unc[_i])
                    fig.suptitle("MVA ID bound: {:.3f}, resolution bounds: [{:.3f},{:.3f}], metric: {:.4f}".format(max_MVAbound, bound1, bound2, np.max(_metrics) / metric_1cat), fontsize=24)
                    fig.tight_layout()
                    path = "./Plots_3cats/bound1_{:.3f}/".format(bound1)
                    if not os.path.exists(path):
                        os.makedirs(path)
                    fig.savefig(path + "bound2_{}.pdf".format(j))
                    plt.close()

        relative_metrics.append(np.max(metrics))
        t1 = time.time()
        print(f"time for the whole 3 cat scan: {t1-t0:.2f}s.")
        metrics_df.to_csv("./metrics_3cats.csv", float_format='%.5f')

        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, 8))
        width = (resboundaries[-1] - resboundaries[0]) / (len(resboundaries) - 1)
        c = ax.imshow(metrics.T / metric_1cat, cmap='viridis', vmin=1, origin="lower", extent=[resboundaries[0], resboundaries[-1] + width, resboundaries[0], resboundaries[-1] + width])

        # Set labels and title
        ax.set_xlabel('Boundary 1 Value')
        ax.set_ylabel('Boundary 2 Value')
        ax.set_title('Sensitivity 3 cats / 1 cat')

        # Display colorbar
        cbar = fig.colorbar(c, ax=ax)
        cbar.set_label('Metric Value')

        num_rows, num_cols = metrics.T.shape
        for i in range(num_rows):
            for j in range(num_cols):
                value = metrics.T[i, j] / metric_1cat
                x = resboundaries[j] + width / 2
                y = resboundaries[i] + width / 2
                if value == 0:
                    continue
                fontsize = 14 if len(resboundaries) < 15 else 10
                ax.text(x, y, f"{value:.2f}", va='center', ha='center', color='white', fontsize=fontsize)
                print(x, y, f"{value:.2f}")

        plt.tight_layout()
        plt.savefig("./sensitivity_3_cats.pdf")

        path_plots_4_cats = "./Plots_4cats/"
        if not os.path.exists(path_plots_4_cats):
            os.makedirs(path_plots_4_cats)

        metrics = np.zeros((len(resboundaries), len(resboundaries), len(resboundaries)))
        metrics_df = pd.DataFrame(columns=["resBound1", "resBound2", "resBound3", "MVA_bound", "relativeMetric"])
        t0 = time.time()

    if n_cat == 4:
        path_plots_4_cats = "./Plots_4cats/"
        if not os.path.exists(path_plots_4_cats):
            os.makedirs(path_plots_4_cats)

        t0 = time.time()
        metrics = np.zeros((len(resboundaries), len(resboundaries), len(resboundaries)))
        metrics_df = pd.DataFrame(columns=["resBound1", "resBound2", "resBound3", "MVA_bound", "relativeMetric"])

        for i, bound1 in enumerate(resboundaries):
            print("\n bound 1:", i)
            if bound1 > 0.012:
                print("WARNING: skipping cats with bound 1 > 0.012 for time reasons. Make sure that this fits your values!")  # these are worse anyways with the values we have here
            for j, bound2 in enumerate(resboundaries):
                print("\t bound 2:", j)
                for k, bound3 in enumerate(resboundaries):
                    if bound2 <= bound1 or bound3 <= bound2:  # Ensuring increasing order of resboundaries
                        continue
                    _metrics = []
                    params = []
                    for iMVA, MVAbound in enumerate(MVAboundaries):
                        events_sig_cats = [
                            signal_events[(signal_events['sigma_m_over_m_decorr'] < bound1) & (signal_events['min_mvaID'] > MVAbound)],
                            signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound1) & (signal_events['sigma_m_over_m_decorr'] < bound2) & (signal_events['min_mvaID'] > MVAbound)],
                            signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound2) & (signal_events['sigma_m_over_m_decorr'] < bound3) & (signal_events['min_mvaID'] > MVAbound)],
                            signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound3) & (signal_events['min_mvaID'] > MVAbound)]
                        ]

                        events_bkg_cats = [
                            bkg_events[(bkg_events['sigma_m_over_m_decorr'] < bound1) & (bkg_events['min_mvaID'] > MVAbound)],
                            bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound1) & (bkg_events['sigma_m_over_m_decorr'] < bound2) & (bkg_events['min_mvaID'] > MVAbound)],
                            bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound2) & (bkg_events['sigma_m_over_m_decorr'] < bound3) & (bkg_events['min_mvaID'] > MVAbound)],
                            bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound3) & (bkg_events['min_mvaID'] > MVAbound)]
                        ]

                        metric_val = 0
                        _params = []  # for later plotting
                        for l_ in range(4):
                            hist_bkg = hist.Hist(hist.axis.Regular(160, 100, 180))
                            hist_bkg.fill(events_bkg_cats[l_].mass.values, weight=events_bkg_cats[l_].weight.values)
                            histo, bin_edges = hist_bkg.to_numpy()
                            bin_centers = [(bin_edges[_l] + bin_edges[_l + 1]) / 2 for _l in range(len(histo))]
                            params_bkg, _ = curve_fit(exponential, bin_centers, histo, p0=([10000, 40, 0]))

                            hist_sig = hist.Hist(hist.axis.Regular(160, 100, 180))
                            hist_sig.fill(events_sig_cats[l_].mass.values, weight=events_sig_cats[l_].weight.values)
                            histo, bin_edges = hist_sig.to_numpy()
                            initial_guess = [5000, 124, 5, 5000, 125, 1]

                            try:
                                params_sig, _ = curve_fit(double_gaussian, bin_centers, histo, p0=initial_guess)
                            except RuntimeError:
                                print(f"WARNING: signal fit in cat. with res. bounds [{bound1:.3f},{bound2:.3f},{bound3:.3f}] and MVA bound {MVAbound:.3f} failed.")
                                params_sig = [0., 124, 999, 0., 125, 999]

                            _params.append(params_sig)
                            _params.append(params_bkg)

                            x_steps = np.linspace(bin_edges[0], bin_edges[-1], 10000)
                            y_vals = double_gaussian(x_steps, *params_sig)
                            dx = x_steps[1] - x_steps[0]
                            y_vals /= np.sum(y_vals * dx)

                            central_interval = get_central_interval(y_vals, dx, area=0.68)
                            signal_integral, _ = quad(double_gaussian, central_interval[0], central_interval[1], args=tuple(params_sig))
                            background_integral, _ = quad(exponential, central_interval[0], central_interval[1], args=tuple(params_bkg))

                            metric_val += (signal_integral / np.sqrt(background_integral)) ** 2

                        metric = np.sqrt(metric_val)
                        _metrics.append(metric)
                        params.append(_params)

                    metrics[i, j, k] = np.max(_metrics)
                    print(np.max(_metrics))

                    # Plot best MVA category for the current resolution category
                    if np.max(_metrics) / metric_1cat > 1.12:
                        metrics_df.loc[len(metrics_df)] = [bound1, bound2, bound3, MVAboundaries[np.argmax(_metrics)], np.max(_metrics) / metric_1cat]
                        max_MVAbound = MVAboundaries[np.argmax(_metrics)]

                        # Select arrays for best MVA bound again
                        events_sig_cats = [
                            signal_events[(signal_events['sigma_m_over_m_decorr'] < bound1) & (signal_events['min_mvaID'] > max_MVAbound)],
                            signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound1) & (signal_events['sigma_m_over_m_decorr'] < bound2) & (signal_events['min_mvaID'] > max_MVAbound)],
                            signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound2) & (signal_events['sigma_m_over_m_decorr'] < bound3) & (signal_events['min_mvaID'] > max_MVAbound)],
                            signal_events[(signal_events['sigma_m_over_m_decorr'] >= bound3) & (signal_events['min_mvaID'] > max_MVAbound)]
                        ]

                        events_bkg_cats = [
                            bkg_events[(bkg_events['sigma_m_over_m_decorr'] < bound1) & (bkg_events['min_mvaID'] > max_MVAbound)],
                            bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound1) & (bkg_events['sigma_m_over_m_decorr'] < bound2) & (bkg_events['min_mvaID'] > max_MVAbound)],
                            bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound2) & (bkg_events['sigma_m_over_m_decorr'] < bound3) & (bkg_events['min_mvaID'] > max_MVAbound)],
                            bkg_events[(bkg_events['sigma_m_over_m_decorr'] >= bound3) & (bkg_events['min_mvaID'] > max_MVAbound)]
                        ]

                        # Histogramming
                        hists_signal = [hist.Hist(hist.axis.Regular(160, 100, 180)) for _ in range(4)]
                        hists_diphoton = [hist.Hist(hist.axis.Regular(160, 100, 180)) for _ in range(4)]
                        hists_GJet = [hist.Hist(hist.axis.Regular(160, 100, 180)) for _ in range(4)]
                        for _i, (event_sig, event_bkg) in enumerate(zip(events_sig_cats, events_bkg_cats)):
                            hists_signal[_i].fill(event_sig.mass.values, weight=event_sig.weight.values)
                            hists_diphoton[_i].fill(event_bkg.mass[event_bkg.process == "Diphoton"].values, weight=event_bkg.weight[event_bkg.process == "Diphoton"].values)
                            hists_GJet[_i].fill(event_bkg.mass[event_bkg.process != "Diphoton"].values, weight=event_bkg.weight[event_bkg.process != "Diphoton"].values)

                        # Extract fitted values
                        params_list = params[np.argmax(_metrics)]
                        y_vals_bkg = [exponential(x_steps, *params_list[i_ * 2 + 1]) for i_ in range(4)]
                        y_vals_sig = [double_gaussian(x_steps, *params_list[i_ * 2]) for i_ in range(4)]

                        # Put to figure
                        fig, axes = plt.subplots(1, 4, figsize=(40, 10))
                        for _i in range(4):
                            plot_category(hists_signal[_i], [hists_GJet[_i], hists_diphoton[_i]], x_steps, y_vals_sig[_i], y_vals_bkg[_i], params_list[i * 2], ax=axes[_i])

                        fig.suptitle("MVA ID bound: {:.3f}, resolution bounds: [{:.3f},{:.3f},{:.3f}], metric: {:.4f}".format(max_MVAbound, bound1, bound2, bound3, np.max(_metrics) / metric_1cat), fontsize=24)
                        fig.tight_layout()
                        path = "./Plots_4cats/bound1_{:.3f}/".format(bound1)
                        if not os.path.exists(path):
                            os.makedirs(path)
                        fig.savefig(path + "bound2_{}_bound3_{}.pdf".format(j,k))
                        plt.close()

            ### plot slices of 3D tensor as heatmap
            metric_ = metrics[i]
            fig, ax = plt.subplots(figsize=(10, 8))
            width = (resboundaries[-1] - resboundaries[0]) / (len(resboundaries) - 1)
            c = ax.imshow(metric_.T / metric_1cat, cmap='viridis', vmin=1, origin="lower", extent=[resboundaries[0], resboundaries[-1] + width, resboundaries[0], resboundaries[-1] + width])
            ax.set_xlabel('Boundary 2 Value')
            ax.set_ylabel('Boundary 3 Value')
            ax.set_title(f'Sensitivity 4 cats / 1 cat, bound 1: {resboundaries[i]:.3f}')
            cbar = fig.colorbar(c, ax=ax)
            cbar.set_label('Metric Value')

            num_rows, num_cols = metric_.T.shape
            for i_ in range(num_rows):
                for j_ in range(num_cols):
                    value = metric_.T[i_, j_] / metric_1cat
                    x = resboundaries[j_] + width / 2
                    y = resboundaries[i_] + width / 2
                    if value == 0:
                        continue
                    fontsize = 14 if len(resboundaries) < 15 else 10
                    ax.text(x, y, f"{value:.2f}", va='center', ha='center', color='white', fontsize=fontsize)
            plt.tight_layout()
            plt.savefig(path_plots_4_cats + f"sensitivity_{i}.pdf")
            plt.clf()

        relative_metrics.append(np.max(metrics))
        metrics_df.to_csv("./metrics_4cats.csv", float_format='%.5f')
        t1 = time.time()
        print(f"time for the whole 4 cat scan: {t1-t0:.2f}s.")

relative_metrics = np.array(relative_metrics) / metric_1cat

labels = ["1 cat.", "2 cat.", "3 cat.", "4 cat."]
categories = [1, 2, 3, 4]
plt.figure(figsize=(10, 6))
plt.scatter(categories, relative_metrics, color='blue', s=100)
for i, label in enumerate(labels):
    plt.annotate(label, (categories[i], relative_metrics[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=14)
plt.xlabel("Number of categories")
plt.ylabel("Relative significance")
plt.xticks(categories)  # Set the x-ticks to be the category numbers
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.ylim(0.9, 1.3)
plt.tight_layout()
plt.savefig("relative_significance.pdf")

from TTbarAnalyzer import TTbarAnalyzer
from Plotter import Plotter
from collections import OrderedDict
from Fitter import Fitter
import matplotlib.pyplot as plt

if __name__ == "__main__":
    """
    Main analysis script. Here you run the analysis and evaluate the results.
    """

    #run_all = True to run both data and Monte Carlo simulation.
    #run_all = False to only run the Monte Carlo simulation. Use this as you default for the design and optimization of the analysis.
    run_all = True

    # List of datasets to be analyzed
    if run_all:
        datasets = OrderedDict([('Data', 'data.root'),
                                ('QCD', 'qcd.root'),
                                ('Diboson', 'diboson.root'),
                                ('DY+jets', 'dy.root'),
                                ('single top', 'single_top.root'),
                                ('TTbar', 'ttbar.root'),
                                ('W+jets', 'wjets.root'),
        ]
        )
    else:
        datasets = OrderedDict([('QCD', 'qcd.root'),
                                ('Diboson', 'diboson.root'),
                                ('DY+jets', 'dy.root'),
                                ('single top', 'single_top.root'),
                                ('TTbar', 'ttbar.root'),
                                ('W+jets', 'wjets.root'),
        ]
        )

    # Options for the event builder
    event_options = {'JEC': 'nominal', # Jet Energy corrections: change to "up" or "down" to evaluate the systematic uncertainties
                     'muon_isolation': 0.1 # muon isolation, you can leave this at the default value
                     }

    analyzers = OrderedDict()

    # Analyze all datasets:
    for name, file_name in datasets.items():
        analyzer = TTbarAnalyzer(name, file_name, event_options) # create an Analyzer for each dataset
        analyzer.run() # run the Analyzer
        analyzers[name] = analyzer # store the results


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Exercise 2: Measurement of the ttbar production cross section
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    L = 50  # integrated luminosity in pb^-1
    del_L = 0.05 * L  # uncertainty of the integrated luminosity in pb^-1
    xs_lit = 173  # literature value of the ttbar production cross section in pb

    
    def eff_pur_xs(n_datasel, n_bkgMC, n_sigsel, n_sigtot, n_MCsel, L, del_L):
        """
        Calculate the efficiency and purity of the selection and its uncertainty.
        Then calculate the ttbar production cross section and its uncertainty.
        
        n_data: number of selected data events
        n_bkg: number of selected background events
        n_sigsel: number of selected signal events
        n_sigtot: total number of signal events
        L: integrated luminosity in pb^-1
        del_L: uncertainty of the integrated luminosity in pb^-1
        
        returns: (cross section, uncertainty)
        """
        e = n_sigsel / n_sigtot  # efficiency
        del_e = ((n_sigsel / n_sigtot**2) + (n_sigsel**2 / n_sigtot**3))**0.5

        p = n_sigsel / n_MCsel
        del_p = ((n_sigsel / n_MCsel**2) + (n_sigsel**2 / n_MCsel**3))**0.5

        xs = (n_datasel - n_bkgMC) / (L * e)  # cross section in pb
        del_xs = ((n_datasel**0.5 / (e * L))**2
                + (n_bkgMC**0.5 / (e * L))**2
                + ((n_bkgMC - n_datasel) * del_e / (L * e**2))**2
                + ((n_bkgMC - n_datasel) * del_L / (L**2 * e))**2 )**0.5
        
        rel_del_xs = del_xs * 100 / xs  # relative uncertainty

        return (e, del_e, p, del_p, xs, del_xs, rel_del_xs)
    

    n_ttbar_total = analyzers['TTbar'].n_total  # get total number of ttbar events
    n_ttbar_trigger = analyzers['TTbar'].n_trigger  # get number of ttbar events passing the trigger
    n_ttbar_n_jets = analyzers['TTbar'].n_n_jets  # get number of ttbar events passing the n_jets selection
    n_ttbar_n_muon = analyzers['TTbar'].n_n_muon  # get number of ttbar events passing the n_muon selection
    n_ttbar_n_b_jets = analyzers['TTbar'].n_n_b_jets  # get number of ttbar events passing the n_b_jets selection
    n_ttbar_met = analyzers['TTbar'].n_met  # get number of ttbar events passing the met selection

    n_data_total = analyzers['Data'].n_total  # get total number of data events
    n_data_trigger = analyzers['Data'].n_trigger  # get number of data events passing the trigger
    n_data_n_jets = analyzers['Data'].n_n_jets  # get number of data events passing the n_jets selection
    n_data_n_muon = analyzers['Data'].n_n_muon  # get number of data events passing the n_muon selection
    n_data_n_b_jets = analyzers['Data'].n_n_b_jets  # get number of data events passing the n_b_jets selection
    n_data_met = analyzers['Data'].n_met  # get number of data events passing the met selection

    n_bkg_total = sum([an.n_total for key, an in analyzers.items() if key != 'Data' and key != 'TTbar'])  # get total number of background events
    n_bkg_trigger = sum([an.n_trigger for key, an in analyzers.items() if key != 'Data' and key != 'TTbar'])  # get number of background events passing the trigger
    n_bkg_n_jets = sum([an.n_n_jets for key, an in analyzers.items() if key != 'Data' and key != 'TTbar'])  # get number of background events passing the n_jets selection
    n_bkg_n_muon = sum([an.n_n_muon for key, an in analyzers.items() if key != 'Data' and key != 'TTbar'])  # get number of background events passing the n_muon selection
    n_bkg_n_b_jets = sum([an.n_n_b_jets for key, an in analyzers.items() if key != 'Data' and key != 'TTbar'])  # get number of background events passing the n_b_jets selection
    n_bkg_met = sum([an.n_met for key, an in analyzers.items() if key != 'Data' and key != 'TTbar'])  # get number of background events passing the met selection

    n_mc_total = sum([an.n_total for key, an in analyzers.items() if key != 'Data'])  # get total number of MC events
    n_mc_trigger = sum([an.n_trigger for key, an in analyzers.items() if key != 'Data'])  # get number of MC events passing the trigger
    n_mc_n_jets = sum([an.n_n_jets for key, an in analyzers.items() if key != 'Data'])  # get number of MC events passing the n_jets selection
    n_mc_n_muon = sum([an.n_n_muon for key, an in analyzers.items() if key != 'Data'])  # get number of MC events passing the n_muon selection
    n_mc_n_b_jets = sum([an.n_n_b_jets for key, an in analyzers.items() if key != 'Data'])  # get number of MC events passing the n_b_jets selection
    n_mc_met = sum([an.n_met for key, an in analyzers.items() if key != 'Data'])  # get number of MC events passing the met selection

    n_sel_total_total = sum([an.n_total for key, an in analyzers.items()])  # get total number of selected events
    n_sel_total_trigger = sum([an.n_trigger for key, an in analyzers.items()])  # get total number of selected events passing the trigger
    n_sel_total_n_jets = sum([an.n_n_jets for key, an in analyzers.items()])  # get total number of selected events passing the n_jets selection
    n_sel_total_n_muon = sum([an.n_n_muon for key, an in analyzers.items()])  # get total number of selected events passing the n_muon selection
    n_sel_total_n_b_jets = sum([an.n_n_b_jets for key, an in analyzers.items()])  # get total number of selected events passing the n_b_jets selection
    n_sel_total_met = sum([an.n_met for key, an in analyzers.items()])  # get total number of selected events passing the met selection

    e_total, del_e_total, p_total, del_p_total, xs_total, del_xs_total, rel_del_xs_total = eff_pur_xs(
        n_data_total,
        n_bkg_total,
        n_ttbar_total,
        n_ttbar_total,
        n_mc_total,
        L,
        del_L
        )

    e_trigger, del_e_trigger, p_trigger, del_p_trigger, xs_trigger, del_xs_trigger, rel_del_xs_trigger = eff_pur_xs(
        n_data_trigger,
        n_bkg_trigger,
        n_ttbar_trigger,
        n_ttbar_total,
        n_mc_trigger,
        L,
        del_L
        )

    e_n_jets, del_e_n_jets, p_n_jets, del_p_n_jets, xs_n_jets, del_xs_n_jets, rel_del_xs_n_jets = eff_pur_xs(
        n_data_n_jets,
        n_bkg_n_jets,
        n_ttbar_n_jets,
        n_ttbar_total,
        n_mc_n_jets,
        L,
        del_L
        )

    e_n_muon, del_e_n_muon, p_n_muon, del_p_n_muon, xs_n_muon, del_xs_n_muon, rel_del_xs_n_muon = eff_pur_xs(
        n_data_n_muon,
        n_bkg_n_muon,
        n_ttbar_n_muon,
        n_ttbar_total,
        n_mc_n_muon,
        L,
        del_L
        )

    e_n_b_jets, del_e_n_b_jets, p_n_b_jets, del_p_n_b_jets, xs_n_b_jets, del_xs_n_b_jets, rel_del_xs_n_b_jets = eff_pur_xs(
        n_data_n_b_jets,
        n_bkg_n_b_jets,
        n_ttbar_n_b_jets,
        n_ttbar_total,
        n_mc_n_b_jets,
        L,
        del_L
        )

    e_met, del_e_met, p_met, del_p_met, xs_met, del_xs_met, rel_del_xs_met = eff_pur_xs(
        n_data_met,
        n_bkg_met,
        n_ttbar_met,
        n_ttbar_total,
        n_mc_met,
        L,
        del_L
        )

    print("________________________________________________________")
    print("Efficiency and purity of the selection:")
    print("________________________________________________________")
    print(f"Total: {e_total:.3f} +- {del_e_total:.3f}, {p_total:.3f} +- {del_p_total:.3f}")
    print(f"Trigger: {e_trigger:.3f} +- {del_e_trigger:.3f}, {p_trigger:.3f} +- {del_p_trigger:.3f}")
    print(f"Number of jets: {e_n_jets:.3f} +- {del_e_n_jets:.3f}, {p_n_jets:.3f} +- {del_p_n_jets:.3f}")
    print(f"Number of muons: {e_n_muon:.3f} +- {del_e_n_muon:.3f}, {p_n_muon:.3f} +- {del_p_n_muon:.3f}")
    print(f"Number of b jets: {e_n_b_jets:.3f} +- {del_e_n_b_jets:.3f}, {p_n_b_jets:.3f} +- {del_p_n_b_jets:.3f}")
    print(f"MET: {e_met:.3f} +- {del_e_met:.3f}, {p_met:.3f} +- {del_p_met:.3f}")
    print("________________________________________________________")
    print("Cross section of ttbar production after each selection step:")
    print("________________________________________________________")
    print(f"Total: {xs_total:.3f} +- {del_xs_total:.3f} pb, ({rel_del_xs_total:.3f}%)")
    print(f"Trigger: {xs_trigger:.3f} +- {del_xs_trigger:.3f} pb, ({rel_del_xs_trigger:.3f}%)")
    print(f"Number of jets: {xs_n_jets:.3f} +- {del_xs_n_jets:.3f} pb, ({rel_del_xs_n_jets:.3f}%)")
    print(f"Number of muons: {xs_n_muon:.3f} +- {del_xs_n_muon:.3f} pb, ({rel_del_xs_n_muon:.3f}%)")
    print(f"Number of b jets: {xs_n_b_jets:.3f} +- {del_xs_n_b_jets:.3f} pb, ({rel_del_xs_n_b_jets:.3f}%)")
    print(f"MET: {xs_met:.3f} +- {del_xs_met:.3f} pb, ({rel_del_xs_met:.3f}%)")
    print("________________________________________________________")

    # plot efficiency and purity in the same plot
    e_p, ax_e_p = plt.subplots()
    ax_e_p.errorbar([1, 2, 3, 4, 5], [e_trigger, e_n_jets, e_n_muon, e_n_b_jets, e_met], yerr=[del_e_trigger, del_e_n_jets, del_e_n_muon, del_e_n_b_jets, del_e_met], fmt='o', label='efficiency')
    ax_e_p.errorbar([0, 1, 2, 3, 4, 5], [p_total, p_trigger, p_n_jets, p_n_muon, p_n_b_jets, p_met], yerr=[del_p_total, del_p_trigger, del_p_n_jets, del_p_n_muon, del_p_n_b_jets, del_p_met], fmt='o', label='purity')
    ax_e_p.set_xticks([0, 1, 2, 3, 4, 5])
    ax_e_p.set_xticklabels(['Total', 'Trigger', 'Number of jets', 'Number of muons', 'Number of b jets', 'MET'], rotation=25)
    ax_e_p.legend()
    ax_e_p.set_title('Efficiency and purity')
    plt.tight_layout()
    plt.savefig('eff_pur.pdf')

    # plot cross section after each selection step
    xs, ax_xs = plt.subplots()
    ax_xs.errorbar([0, 1, 2, 3, 4, 5], [xs_total, xs_trigger, xs_n_jets, xs_n_muon, xs_n_b_jets, xs_met], yerr=[del_xs_total, del_xs_trigger, del_xs_n_jets, del_xs_n_muon, del_xs_n_b_jets, del_xs_met], fmt='o')
    ax_xs.set_xticks([0, 1, 2, 3, 4, 5])
    ax_xs.set_xticklabels(['Total', 'Trigger', 'Number of jets', 'Number of muons', 'Number of b jets', 'MET'], rotation=25)
    ax_xs.set_title('Cross section after each selection step')
    plt.tight_layout()
    plt.savefig('xs.pdf')

    # plot cross section without trigger selection step
    xs_no_trigger, ax_xs_no_trigger = plt.subplots()
    ax_xs_no_trigger.errorbar([0, 1, 2, 3, 4], [xs_total, xs_n_jets, xs_n_muon, xs_n_b_jets, xs_met], yerr=[del_xs_total, del_xs_n_jets, del_xs_n_muon, del_xs_n_b_jets, del_xs_met], fmt='o')
    ax_xs_no_trigger.set_xticks([0, 1, 2, 3, 4])
    ax_xs_no_trigger.set_xticklabels(['Total', 'Number of jets', 'Number of muons', 'Number of b jets', 'MET'], rotation=25)
    ax_xs_no_trigger.set_title('Cross section after each selection step without trigger selection')
    plt.tight_layout()
    plt.savefig('xs_no_trigger.pdf')

    def dif_to_lit(xs, del_xs):
        dif = xs - xs_lit
        sigma = dif / del_xs
        if sigma < 0:
            sigma = -sigma
        return dif, sigma


    print("Difference to literature value and significance:")
    print("________________________________________________________")
    print(f"Total: {dif_to_lit(xs_total, del_xs_total)}")
    print(f"Trigger: {dif_to_lit(xs_trigger, del_xs_trigger)}")
    print(f"Number of jets: {dif_to_lit(xs_n_jets, del_xs_n_jets)}")
    print(f"Number of muons: {dif_to_lit(xs_n_muon, del_xs_n_muon)}")
    print(f"Number of b jets: {dif_to_lit(xs_n_b_jets, del_xs_n_b_jets)}")
    print(f"MET: {dif_to_lit(xs_met, del_xs_met)}")
    print("________________________________________________________")

    # Plot all histograms filled in the Analysis
    plotter = Plotter(analyzers)
    plotter.process()

    # Run the fit of the top mass distribution
    fitter = Fitter(analyzers)
    fitter.fit(130., 210.)
    # fitter.fit(x,y)
    # (x,y) = fit range

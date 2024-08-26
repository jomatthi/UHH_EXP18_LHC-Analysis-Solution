from TTbarAnalyzer import TTbarAnalyzer
from Plotter import Plotter
from collections import OrderedDict
from Fitter import Fitter

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

    # You can access variables from TTbarAnalyzer like this:
    n_background_total = sum([an.n_total for key, an in analyzers.items() if not (key == 'Data' or key == 'TTbar')]) # Sum total number of all events, except data and ttbar

    n_ttbar_total = analyzers['TTbar'].n_total # get total number of ttbar events
    n_ttbar_trigger = analyzers['TTbar'].n_trigger # get number of ttbar events passing the trigger
    n_ttbar_n_muon = analyzers['TTbar'].n_n_muon # get number of ttbar events passing the n_muon selection
    n_ttbar_n_jets = analyzers['TTbar'].n_n_jets # get number of ttbar events passing the n_jets selection
    n_ttbar_n_b_jets = analyzers['TTbar'].n_n_b_jets # get number of ttbar events passing the n_b_jets selection
    n_ttbar_met = analyzers['TTbar'].n_met # get number of ttbar events passing the met selection

    n_data_total = analyzers['Data'].n_total # get total number of data events
    n_data_trigger = analyzers['Data'].n_trigger # get number of data events passing the trigger
    n_data_n_muon = analyzers['Data'].n_n_muon # get number of data events passing the n_muon selection
    n_data_n_jets = analyzers['Data'].n_n_jets # get number of data events passing the n_jets selection
    n_data_n_b_jets = analyzers['Data'].n_n_b_jets # get number of data events passing the n_b_jets selection
    n_data_met = analyzers['Data'].n_met # get number of data events passing the met selection

    bkg_processes = ['QCD', 'Diboson', 'DY+jets', 'single top', 'W+jets']

    n_bkg_total = sum([an.n_total for key, an in analyzers.items() if key in bkg_processes])
    n_bkg_trigger = sum([an.n_trigger for key, an in analyzers.items() if key in bkg_processes])
    n_bkg_n_muon = sum([an.n_n_muon for key, an in analyzers.items() if key in bkg_processes])
    n_bkg_n_jets = sum([an.n_n_jets for key, an in analyzers.items() if key in bkg_processes])
    n_bkg_n_b_jets = sum([an.n_n_b_jets for key, an in analyzers.items() if key in bkg_processes])
    n_bkg_met = sum([an.n_met for key, an in analyzers.items() if key in bkg_processes])

    n_mc_total = sum([an.n_total for key, an in analyzers.items() if key != 'Data'])
    n_mc_trigger = sum([an.n_trigger for key, an in analyzers.items() if key != 'Data'])
    n_mc_n_muon = sum([an.n_n_muon for key, an in analyzers.items() if key != 'Data'])
    n_mc_n_jets = sum([an.n_n_jets for key, an in analyzers.items() if key != 'Data'])
    n_mc_n_b_jets = sum([an.n_n_b_jets for key, an in analyzers.items() if key != 'Data'])
    n_mc_met = sum([an.n_met for key, an in analyzers.items() if key != 'Data'])

    n_sel_total_met = sum([an.n_met for key, an in analyzers.items()])
    n_sel_total_n_b_jets = sum([an.n_n_b_jets for key, an in analyzers.items()])

    L = 50  # integrated luminosity in pb^-1
    del_L = 0.05 * L  # uncertainty of the integrated luminosity in pb^-1

    
    def calculate_xs(n_data, n_bkg, n_sigsel, n_sigtot, L, del_L):
        """
        Calculate the cross section and its uncertainty
        
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

        xs = (n_data - n_bkg) / (L * e)  # cross section in pb
        del_xs = ((n_data**0.5 / (e * L))**2
                + (n_bkg**0.5 / (e * L))**2
                + ((n_bkg - n_data) * del_e / (L * e**2))**2
                + ((n_bkg - n_data) * del_L / (L**2 * e))**2 )**0.5
        
        rel_del_xs = del_xs / xs  # relative uncertainty

        return (xs, del_xs, rel_del_xs)


    print("________________________________________________________")
    print("Total Number of ttbar events: {0}".format(n_ttbar_total)) # write to the console
    print("Total Number of background events: {0}".format(n_background_total))
    print("Total Number of data events: {0}".format(n_data_total))
    print("Number of data events passing the n_b_jets selection: {0}".format(n_data_n_b_jets))
    print(f"Number of ttbar events passing the n_b_jets selection: {n_ttbar_n_b_jets}")
    print(f"Number of background events passing the n_b_jets selection: {n_bkg_n_b_jets}")
    print(f"Number of data events passing the met selection: {n_data_met}")
    print(f"Number of ttbar events passing the met selection: {n_ttbar_met}")
    print(f"Number of background events passing the met selection: {n_bkg_met}")

    print(f"Number of selected events (with MET cut): {n_sel_total_met}" if not n_sel_total_met == n_data_met + n_ttbar_met + n_bkg_met else " ")
    print(f"Number of selected events (w/o MET cut): {n_sel_total_n_b_jets}" if not n_sel_total_n_b_jets == n_data_n_b_jets + n_ttbar_n_b_jets + n_bkg_n_b_jets else " ")

    print(f"Efficiency after no selection: {n_ttbar_total/n_ttbar_total}")
    print(f"Efficiency after trigger selection: {n_ttbar_trigger/n_ttbar_total}")
    print(f"Efficiency after n_jets selection: {n_ttbar_n_jets/n_ttbar_total}")
    print(f"Efficiency after n_muon selection: {n_ttbar_n_muon/n_ttbar_total}")
    print(f"Efficiency after n_b_jets selection: {n_ttbar_n_b_jets/n_ttbar_total}")
    print(f"Efficiency after met selection: {n_ttbar_met/n_ttbar_total}")

    print("________________________________________________________")

    print(f"Purity after no selection: {n_ttbar_total/n_mc_total}")
    print(f"Purity after trigger selection: {n_ttbar_trigger/n_mc_trigger}")
    print(f"Purity after n_jets selection: {n_ttbar_n_jets/n_mc_n_jets}")
    print(f"Purity after n_muon selection: {n_ttbar_n_muon/n_mc_n_muon}")
    print(f"Purity after n_b_jets selection: {n_ttbar_n_b_jets/n_mc_n_b_jets}")
    print(f"Purity after met selection: {n_ttbar_met/n_mc_met}")

    print("________________________________________________________")

    print(f"Cross section after MET selection: {calculate_xs(n_data_met, n_bkg_met, n_ttbar_met, n_ttbar_total, L, del_L)}")
    print(f"Cross section after n_b_jets selection: {calculate_xs(n_data_n_b_jets, n_bkg_n_b_jets, n_ttbar_n_b_jets, n_ttbar_total, L, del_L)}")
    print(f"Cross section after n_jets selection: {calculate_xs(n_data_n_jets, n_bkg_n_jets, n_ttbar_n_jets, n_ttbar_total, L, del_L)}")
    print(f"Cross section after n_muon selection: {calculate_xs(n_data_n_muon, n_bkg_n_muon, n_ttbar_n_muon, n_ttbar_total, L, del_L)}")
    print(f"Cross section after trigger selection: {calculate_xs(n_data_trigger, n_bkg_trigger, n_ttbar_trigger, n_ttbar_total, L, del_L)}")
    print(f"Cross section after no selection: {calculate_xs(n_data_total, n_bkg_total, n_ttbar_total, n_ttbar_total, L, del_L)}")

    # Plot all histograms filled in the Analysis
    plotter = Plotter(analyzers)
    plotter.process()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Exercise 3: Reconstruction of the top quark mass
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Run the fit of the top mass distribution
    fitter = Fitter(analyzers)
    fitter.fit(130., 210.)
    # fitter.fit(x,y)
    # (x,y) = fit range

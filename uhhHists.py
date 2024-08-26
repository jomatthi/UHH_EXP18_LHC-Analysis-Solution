from Histograms import Histograms
from ROOT import TH1F
from collections import OrderedDict

class DefaultHistograms(Histograms):
    """
    A default set of Histograms for the ttbar analysis.
    """
    def __init__(self, name):
        self.hists = OrderedDict([('muons_number', TH1F('muons_number','N_{#mu}', 11,-0.5,10.5)),
                      ('muon1_pt',     TH1F('muon1_pt','p_{T,#mu} [GeV]', 60,0,300)),
                      ('muon1_eta',    TH1F('muon1_eta','#eta_{#mu}', 50,-5.0,5.0)),
                      ('muon1_phi',    TH1F('muon1_phi','#phi_{#mu}', 40,-3.2,3.2)),
                      ('jets_number' , TH1F('jets_number', 'N_{jets}', 11, -0.5, 10.5)),
                      ('jet1_pt',      TH1F('jet1_pt','p_{T,jet} [GeV]', 60,0,300)),
                      ('jet1_eta',     TH1F('jet1_eta','#eta_{jet}', 50,-5.0,5.0)),
                      ('jet1_phi',     TH1F('jet1_phi','#phi_{jet}', 40,-3.2,3.2)),
                      ('jet2_pt',      TH1F('jet2_pt','p_{T,jet} [GeV]', 60,0,300)),
                      ('jet2_eta',     TH1F('jet2_eta','#eta_{jet}', 50,-5.0,5.0)),
                      ('jet2_phi',     TH1F('jet2_phi','#phi_{jet}', 40,-3.2,3.2)),
                      ('jet3_pt',      TH1F('jet3_pt','p_{T,jet} [GeV]', 60,0,300)),
                      ('jet3_eta',     TH1F('jet3_eta','#eta_{jet}', 50,-5.0,5.0)),
                      ('jet3_phi',     TH1F('jet3_phi','#phi_{jet}', 40,-3.2,3.2)),
                      ('met_pt',       TH1F('met_pt','p_{T,miss}', 60,0,300)),
                      ('met_phi',      TH1F('met_phi','#phi_{MET}', 40,-3.2,3.2)),
                      ('bjets_number' , TH1F('bjets_number', 'N_{b-jets}', 11, -0.5, 10.5)),
                      ('bjet1_pt',      TH1F('bjet1_pt','p_{T,b-jet} [GeV]', 60,0,300)),
                      ('bjet1_eta',     TH1F('bjet1_eta','#eta_{b-jet}', 50,-5.0,5.0)),
                      ('bjet1_phi',     TH1F('bjet1_phi','#phi_{b-jet}', 40,-3.2,3.2)),
                      ('bjet2_pt',      TH1F('bjet2_pt','p_{T,b-jet} [GeV]', 60,0,300)),
                      ('bjet2_eta',     TH1F('bjet2_eta','#eta_{b-jet}', 50,-5.0,5.0)),
                      ('bjet2_phi',     TH1F('bjet2_phi','#phi_{b-jet}', 40,-3.2,3.2)),
                      ])
        ## DO NOT TOUCH THIS PART ##
        name = name + "_default"
        super(DefaultHistograms, self).__init__(name)

    def fill(self, event):
        """
        Here the histograms are filled.
        """

        event_weight = event.weight
        
        # fill muon hists
        self.hists['muons_number'].Fill(event.n_muons(), event_weight)
        if event.n_muons() >= 1:
            muon = event.muons[0]
            self.hists['muon1_pt'].Fill(muon.pt(), event_weight)
            self.hists['muon1_eta'].Fill(muon.eta(), event_weight)
            self.hists['muon1_phi'].Fill(muon.phi(), event_weight)
        
        # fill jet hists
        self.hists['jets_number'].Fill(event.n_jets(), event_weight)
        i_jet = 0
        for jet in event.jets:
            i_jet += 1
            if i_jet == 1:
                self.hists['jet1_pt'].Fill(jet.pt(), event_weight)
                self.hists['jet1_eta'].Fill(jet.eta(), event_weight)
                self.hists['jet1_phi'].Fill(jet.phi(), event_weight)
            elif i_jet == 2:
                self.hists['jet2_pt'].Fill(jet.pt(), event_weight)
                self.hists['jet2_eta'].Fill(jet.eta(), event_weight)
                self.hists['jet2_phi'].Fill(jet.phi(), event_weight)
            elif i_jet == 3:
                self.hists['jet3_pt'].Fill(jet.pt(), event_weight)
                self.hists['jet3_eta'].Fill(jet.eta(), event_weight)
                self.hists['jet3_phi'].Fill(jet.phi(), event_weight)
            elif i_jet > 3:
                break

        # fill met hists
        self.hists['met_pt'].Fill(event.met.pt(), event_weight)
        self.hists['met_phi'].Fill(event.met.phi(), event_weight)

        # fill b-jet hists
        self.hists['bjets_number'].Fill(event.n_b_jets(), event_weight)
        i_bjet = 0
        for b_jet in event.b_jets:
            i_bjet += 1
            if i_bjet == 1:
                self.hists['bjet1_pt'].Fill(b_jet.pt(), event_weight)
                self.hists['bjet1_eta'].Fill(b_jet.eta(), event_weight)
                self.hists['bjet1_phi'].Fill(b_jet.phi(), event_weight)
            elif i_bjet == 2:
                self.hists['bjet2_pt'].Fill(b_jet.pt(), event_weight)
                self.hists['bjet2_eta'].Fill(b_jet.eta(), event_weight)
                self.hists['bjet2_phi'].Fill(b_jet.phi(), event_weight)
            elif i_bjet > 2:
                break
        

        # other histograms go here

class TopMassHist(Histograms):
    def __init__(self, name):
        self.hists = OrderedDict([('top_mass', TH1F('top_mass',';M_{T}', 120,0,300)),
                                  ])
        ## DO NOT TOUCH THIS PART ##
        name = name + "_default"
        super(TopMassHist, self).__init__(name)

    def fill(self, event):
        """
        Here the histograms are filled.
        """

        event_weight = event.weight

        #fill top hists
        if event.top_mass > 0.0:
            self.hists['top_mass'].Fill(event.top_mass, event_weight)

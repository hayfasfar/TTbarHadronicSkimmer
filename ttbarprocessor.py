#!/usr/bin/env python 
# coding: utf-8

from coffea import processor, nanoevents
from coffea import util
from coffea.btag_tools import BTagScaleFactor
from coffea.nanoevents.methods import candidate
from coffea.nanoevents.methods import vector
from coffea.jetmet_tools import JetResolutionScaleFactor
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory
from coffea.lookup_tools import extractor
from coffea.analysis_tools import Weights, PackedSelection
from collections import defaultdict
import sys
import os, psutil
import copy
import scipy.stats as ss
import numpy as np
import itertools
import pandas as pd
from numpy.random import RandomState
import random
import correctionlib
import hist
import json
import logging
import psutil
import time

import awkward as ak

# for dask, `from python.corrections import` does not work
sys.path.append(os.getcwd()+'/python/')

from corrections import (
    GetFlavorEfficiency,
    HEMCleaning,
    HEMVeto,
    GetL1PreFiringWeight,
    GetJECUncertainties,
    GetPDFWeights,
    GetPUSF,
    GetQ2weights,
    getLumiMaskRun2,
    getMETFilter,
    pTReweighting,
)
from btagCorrections import btagCorrections
from functions import getRapidity




# logging
# logfile = 'coffea_' + str(int(time.time())) + '.log'
# print(logfile)
# logging.basicConfig(filename=logfile, level=logging.DEBUG)
logger = logging.getLogger('__main__')
logger.setLevel(logging.DEBUG)


#ak.behavior.update(candidate.behavior)
ak.behavior.update(vector.behavior)


# --- Define 'Manual bins' to use for mistag plots for aesthetic purposes--- #
manual_bins = [400, 500, 600, 800, 1000, 1500, 2000, 3000, 7000, 10000]




def get_memory_usage():
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    memory_usage_bytes = memory_info.rss
    memory_usage_mb = memory_usage_bytes / (1024 * 1024)

    return memory_usage_mb



def update(events, collections):
    # https://github.com/nsmith-/boostedhiggs/blob/master/boostedhiggs/hbbprocessor.py
    """Return a shallow copy of events array with some collections swapped out"""
    out = events
#     logger.debug('update:%s:%s', time.time(), collections)
    
    for name, value in collections.items():
        out = ak.with_field(out, value, name)

    return out


"""Package to perform the data-driven mistag-rate-based ttbar hadronic analysis. """
class TTbarResProcessor(processor.ProcessorABC):
    def __init__(self,
                 htCut=1400.,
                 ak8PtMin=400.,
                 minMSD=105.,
                 maxMSD=210.,
                 tau32Cut=0.65,
                 bdisc=0.5847,
                 deepAK8Cut='medium',
                 useDeepAK8=True,
                 useDeepCSV=True,
                 iov='2016',
                 bkgEst=False,
                 noSyst=False,
                 blinding=False,
                 systematics = ['nominal', 'pileup'],
                 anacats = ['2t0bcen'],
                 #rpf_params = {'params':[1.0], 'errors':[0.0]},
                ):
                 
        self.iov = iov
        self.htCut = htCut
        self.minMSD = minMSD
        self.maxMSD = maxMSD
        self.tau32Cut = tau32Cut
        self.ak8PtMin = ak8PtMin
        self.bdisc = bdisc
        self.deepAK8Cut = deepAK8Cut
        self.useDeepAK8 = useDeepAK8
        self.useDeepCSV = useDeepCSV
        self.means_stddevs = defaultdict()
        self.bkgEst = bkgEst
        self.noSyst = noSyst
        self.systematics = systematics
        self.blinding = blinding
        #self.rpf_params = rpf_params        
        
#         self.transfer_function = np.load('plots/save.npy')

#         deepak8cuts = {
#             'loose':{ # 1%
#                 '2016APV': 0.486, 
#                 '2016': 0.475,
#                 '2017': 0.487,
#                 '2018': 0.477,
#             },
#             'medium':{ # 0.5%
#                 '2016APV': 0.677, 
#                 '2016': 0.666,
#                 '2017': 0.673,
#                 '2018': 0.669,
#             },
#             'tight': { # 0.1%
#                 '2016APV': 0.902, 
#                 '2016': 0.897,
#                 '2017': 0.898,
#                 '2018': 0.900,
#             } 
#         }
    
        # from https://twiki.cern.ch/twiki/bin/view/CMS/DeepAK8Tagging2018WPsSFs#2016_Data
        deepak8cuts = {
            'loose':{ # 1%
                '2016APV': 0.435, 
                '2016':    0.435,
                '2017':    0.344,
                '2018':    0.470,
            },
            'medium':{ # 0.5%
                '2016APV': 0.632, 
                '2016':    0.632,
                '2017':    0.554,
                '2018':    0.685,
            },
            'tight': { # 0.1%
                '2016APV': 0.889, 
                '2016':    0.889,
                '2017':    0.863,
                '2018':    0.920,
            } 
        }
    
        
        
        self.weights = {}
    
        
        
        self.deepAK8disc = deepak8cuts[deepAK8Cut][self.iov]
        
        if deepAK8Cut == 'tight':
#             self.deepAK8low = deepak8cuts['loose'][self.iov]
            self.deepAK8low = deepak8cuts['medium'][self.iov]
        elif deepAK8Cut == 'medium':
             self.deepAK8low = deepak8cuts['loose'][self.iov]
#            self.deepAK8low = deepak8cuts['medium'][self.iov]
        else:
            self.deepAK8low = 0.2
        
        
        
        
        
        # analysis categories #
        self.anacats = anacats
        self.label_dict = {i: label for i, label in enumerate(self.anacats)}
        self.label_to_int_dict = {label: i for i, label in enumerate(self.anacats)}

        
        # systematics
        syst_category_strings = ['nominal']
        if not self.noSyst:
            for s in self.systematics:
                if (s != 'nominal'):
                    
                    if ('hem' in s):
                        syst_category_strings.append(s)
                    else:
                        syst_category_strings.append(s+'Down')
                        syst_category_strings.append(s+'Up')
        
#         syst_category_strings = ['nominal', 'test1', 'test2', 'test3', 'test4']
        
        # axes
        dataset_axis     = hist.axis.StrCategory([], growth=True, name="dataset", label="Primary Dataset")
        syst_axis        = hist.axis.StrCategory(syst_category_strings, name="systematic")
        ttbarmass_axis   = hist.axis.Regular(50, 800, 8000, name="ttbarmass", label=r"$m_{t\bar{t}}$ [GeV]")
        jetmass_axis     = hist.axis.Regular(50, 0, 500, name="jetmass", label=r"Jet $m$ [GeV]")
        jetmsd_axis      = hist.axis.Regular(20, 0, 500, name="jetmass", label=r"Jet $m_{SD}$ [GeV]")
        ttbarmass2D_axis = hist.axis.Regular(60, 800, 6800, name="ttbarmass", label=r"$m_{t\bar{t}}$ [GeV]")
        jetmass2D_axis   = hist.axis.Regular(100, 0, 500, name="jetmass", label=r"Jet $m_{SD}$ [GeV]")
        ht_axis        = hist.axis.Regular(40, 400, 4400, name="ht", label=r"$H_T$ [GeV]")
        cats_axis        = hist.axis.IntCategory(range(len(self.anacats)), name="anacat", label="Analysis Category")
        manual_axis      = hist.axis.Variable(manual_bins, name="jetp", label=r"Jet Momentum [GeV]")
        nsub_axis        = hist.axis.Regular(10, 0, 1, name="nsub", label=r"$\tau_{3} / \tau_{2}$")
        njet_axis        = hist.axis.Regular(10,1,11, name="njet", label=r"$N_{FatJet}$")

        
        # output
        self.histo_dict = {

            
            # histograms
            'ttbarmass'  : hist.Hist(syst_axis, cats_axis, ttbarmass2D_axis, storage="weight", name="Counts"),
            'mtt_unwgt'  : hist.Hist(syst_axis, cats_axis, ttbarmass2D_axis, storage="weight", name="Counts"),
            'numerator'  : hist.Hist(cats_axis, manual_axis, storage="weight", name="Counts"),
            'denominator': hist.Hist(cats_axis, manual_axis, storage="weight", name="Counts"),
            'jetmass' : hist.Hist(syst_axis, cats_axis, jetmass2D_axis, storage="weight", name="Counts"),
            'jetmsd' : hist.Hist(syst_axis, cats_axis, jetmsd_axis, storage="weight", name="Counts"),

            # second leading jet
            'jetmass1' : hist.Hist(syst_axis, cats_axis, jetmass2D_axis, storage="weight", name="Counts"),
            'jetmsd1' : hist.Hist(syst_axis, cats_axis, jetmsd_axis, storage="weight", name="Counts"),

            
            'ht'  : hist.Hist(syst_axis, cats_axis, ht_axis, storage="weight", name="Counts"),
            
            
            'mtt_vs_mt' : hist.Hist(syst_axis, cats_axis, jetmass2D_axis, ttbarmass2D_axis, storage="weight", name="Counts"),
            
            
            
            # accumulators
            'cutflow': processor.defaultdict_accumulator(int),
            'weights': processor.defaultdict_accumulator(float),
            'systematics': processor.defaultdict_accumulator(float),
        }
        
      

        
    @property
    def accumulator(self):
        return self._accumulator
    
    
    
    def process(self, events):
                
        
        logger.debug('memory:%s:start %s preprocessor: %s', time.time(), events.metadata['dataset'], get_memory_usage())

                
        # reference for return processor.accumulate
        # https://github.com/nsmith-/boostedhiggs/blob/master/boostedhiggs/hbbprocessor.py
        
        nEvents = len(events.event)

        
        # Remove events with large weights
        if "QCD" in events.metadata['dataset']: # and ('2017' not in self.iov): 
            events = events[ events.Generator.binvar > 400 ]
        
            if events.metadata['dataset'] not in self.means_stddevs : 
                average = np.average( events.genWeight )
                stddev = np.std( events.genWeight )
                self.means_stddevs[events.metadata['dataset']] = (average, stddev)            
            average,stddev = self.means_stddevs[events.metadata['dataset']]
            vals = (events.genWeight - average ) / stddev
            events = events[(np.abs(vals) < 2)]

        
        isData = ('JetHT' in events.metadata['dataset']) or ('SingleMu' in events.metadata['dataset'])
        
        noCorrections = (not 'jes' in self.systematics and not 'jer' in self.systematics)

        if noCorrections or self.noSyst or isData:
            return self.process_analysis(events, 'nominal', nEvents)
        
        
     #   if isData:
     #       
     #       return processor.accumulate([
     #           self.process_analysis(events, 'nominal', nEvents),
     #           self.process_analysis(events, 'hemVeto', nEvents)
     #       ]) 
        
        
        FatJets = events.FatJet
        GenJets = events.GenJet
        Jets = events.Jet
        Met = events.MET
        #print("dictionary of MET ", events.MET.fields)
        
                
        
        FatJets["p4"] = ak.with_name(FatJets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        GenJets["p4"] = ak.with_name(GenJets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        Jets["p4"]    = ak.with_name(Jets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")

        
        FatJets["p4"] = ak.with_name(FatJets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        GenJets["p4"] = ak.with_name(GenJets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        Jets["p4"]    = ak.with_name(Jets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")

        FatJets["matched_gen_0p2"] = FatJets.p4.nearest(GenJets.p4, threshold=0.2)
        FatJets["pt_gen"] = ak.values_astype(ak.fill_none(FatJets.matched_gen_0p2.pt, 0), np.float32)

        Jets["matched_gen_0p2"] = Jets.p4.nearest(GenJets.p4, threshold=0.2)
        Jets["pt_gen"] = ak.values_astype(ak.fill_none(Jets.matched_gen_0p2.pt, 0), np.float32)


        corrected_fatjets = GetJECUncertainties(FatJets, events, self.iov, R='AK8', isData=isData)
        corrected_jets = GetJECUncertainties(Jets, events, self.iov, R='AK4', isData=isData)
        
        
        if (len(corrected_jets.pt[0]) > 1) and  (len(corrected_fatjets.pt[0]) > 1) :
        
            logger.debug('JEC:%s:JES up, nom, down:%s:%s:%s', time.time(), 
                         corrected_fatjets.JES_jes.up.pt[0][0],
                         corrected_fatjets.pt[0][0],
                         corrected_fatjets.JES_jes.down.pt[0][0])

            logger.debug('JEC:%s:JER up, nom, down:%s:%s:%s', time.time(), 
                         corrected_fatjets.JER.up.pt[0][0],
                         corrected_fatjets.pt[0][0],
                         corrected_fatjets.JER.down.pt[0][0])

            logger.debug('JEC:%s:JES up, nom, down AK4:%s:%s:%s', time.time(), 
                         corrected_jets.JES_jes.up.pt[0][0],
                         corrected_jets.pt[0][0],
                         corrected_jets.JES_jes.down.pt[0][0])

            logger.debug('JEC:%s:JER up, nom, down AK4:%s:%s:%s', time.time(), 
                         corrected_jets.JER.up.pt[0][0],
                         corrected_jets.pt[0][0],
                         corrected_jets.JER.down.pt[0][0])



        
        
        if 'jes' in self.systematics:
            corrections = [
                ({"Jet": corrected_jets, "FatJet": corrected_fatjets}, 'nominal'),
                ({"Jet": corrected_jets.JES_jes.up, "FatJet": corrected_fatjets.JES_jes.up}, "jesUp"),
                ({"Jet": corrected_jets.JES_jes.down, "FatJet": corrected_fatjets.JES_jes.down}, "jesDown"),
            ]
        if 'jer' in self.systematics:
            corrections.extend([
                ({"Jet": corrected_jets.JER.up, "FatJet": corrected_fatjets.JER.up}, "jerUp"),
                ({"Jet": corrected_jets.JER.down, "FatJet": corrected_fatjets.JER.down}, "jerDown"),
            ])
            
            
        
            
#         if ('hem' in self.systematics) and ('2018' in self.iov):
            
#             if 'hemVeto' in self.systematics:
                
            
#                 corrections.extend([
#                     ({"Jet": corrected_jets, "FatJet": corrected_fatjets}, "hemVeto"),
#                     ])

            
#             corrected_jets_hem    = HEMCleaning(corrected_jets)
#             corrected_fatjets_hem = HEMCleaning(corrected_fatjets)

#             corrections.extend([
#             ({"Jet": corrected_jets_hem, "FatJet": corrected_fatjets_hem}, "hem"),
#             ])
            
#             corrections.extend([
#             ({"Jet": Jets, "FatJet": FatJets}, "hem"),
#             ])

                

        
        
#         # get nominal output
#         output_total = self.process_analysis(update(events, corrections[0][0]), 'nominal', nEvents)
        
# #         logger.debug('output:%s:nominal:%s:%s', time.time(), output_total['cutflow'], output_total['systematics'])
# #         logger.debug('output:%s:nominal:%s', time.time(), output_total['weights'])
        
#         # loop through corrections
#         outputs = {}
#         for collections, name in corrections[1:]:
#             process_output = self.process_analysis(update(events, collections), name, nEvents)
#             outputs[name] = process_output
            
# #             logger.debug('output:%s:%s:%s', time.time(), name, process_output['weights'])


#         # combine outputs
#         for name, output_correction in outputs.items():
#             for key in output_total.keys():

#                 if 'hist' in str(type(output_total[key])):
#                     if 'systematic' in list(output_total[key].axes.name):
#                         output_total[key] += output_correction[key]

#                 elif 'accumulator' in str(type(output_total[key])):
#                     if key != 'cutflow':
#                         output_total[key][name] = process_output[key][name]


        
        # loop through corrections
        outputs = []
        for collections, name in corrections:
            outputs.append(self.process_analysis(update(events, collections), name, nEvents))


        return outputs[0]
     


    def process_analysis(self, events, correction, nEvents):
        
        dataset = events.metadata['dataset']
        filename = events.metadata['filename']
        
        logger.debug('memory:%s: start processor %s:%s', time.time(), correction, get_memory_usage())

                
        isNominal = (correction=='nominal')
        isData = ('JetHT' in dataset) or ('SingleMu' in dataset)

        
        
        if (self.iov == '2018'):
            
            if isData:
                                    
                events = events[HEMVeto(events.Jet, events.FatJet, events.run)]


            else:
                events = events[HEMVeto(events.Jet, events.FatJet, events.run)]
                

        
                
        output = self.histo_dict 
        
        if isNominal:
            output['cutflow']['all events 1'] += nEvents
        
        
        # lumi mask #
        if (isData):
            
            lumi_mask = np.array(getLumiMaskRun2(self.iov)(events.run, events.luminosityBlock), dtype=bool)
            events = events[lumi_mask]
            del lumi_mask

        
        
        # blinding #
        if self.blinding:
            if isData: #and (('2017' in self.iov) or ('2018' in self.iov)):
                events = events[::10]
         
            
        
        # event selection #
        selection = PackedSelection()

        # trigger cut #

        triggernames = { 

        "2016APV": ["PFHT900", "AK8PFJet450"],
        "2016" : ["PFHT900", "AK8PFJet450"],
        "2017" : ["PFHT1050", "AK8PFJet500"],
        "2018" : ["PFHT1050", "AK8PFJet500"],

        }

        try:
            selection.add('trigger', (events.HLT[triggernames[self.iov][0]] | events.HLT[triggernames[self.iov][1]]) )
        except:
            selection.add('trigger', (events.HLT[triggernames[self.iov][0]]))


        # objects #
        
            
        FatJets = events.FatJet
        SubJets = events.SubJet
        Jets    = events.Jet
        Met = events.MET


       
        FatJets["p4"] = ak.with_name(FatJets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        SubJets["p4"] = ak.with_name(SubJets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        Jets["p4"]    = ak.with_name(Jets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")

        if not isData:
            GenJets = events.GenJet
            GenJets["p4"] = ak.with_name(GenJets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
                    
        
        logger.debug('memory:%s: get nanoAOD objects %s:%s', time.time(), correction, get_memory_usage())

        
        
        
        # ---- Get event weights from dataset ---- #

        # if blinding + trigger results in too few events
        if (len(events) < 10): return output
        

                
        
        if isData:
            evtweights = np.ones(len(events))
        else:
            if "LHEWeight_originalXWGTUP" not in events.fields: 
                evtweights = events.genWeight
#                 print('genWeight', evtweights) 
            else: 
                evtweights = events.LHEWeight_originalXWGTUP
#                 print('LHEWeight_originalXWGTUP', evtweights) 
#                 print('events.genWeight', events.genWeight) 


        if isNominal:
            output['cutflow']['all events'] += len(FatJets)
            output['cutflow']['sumw'] += np.sum(evtweights)
            output['cutflow']['sumw2'] += np.sum(evtweights**2)

            
        
            
        
        # ---- event selection and object selection ---- #


        jetht = ak.sum(Jets[(Jets.pt>30) & (np.abs(Jets.eta)<3.0)].pt, axis=1)
        


        # ht cut #
        selection.add('htCut',
            ak.sum(Jets[(Jets.pt > 30) & (np.abs(Jets.eta) < 3.0)].pt, axis=1) > self.htCut
        )

        
        # met filters #
        selection.add('metfilter', getMETFilter(self.iov, events))
                
        # jet id #
        selection.add('jetid', ak.any((FatJets.jetId > 1), axis=1))
        FatJets = FatJets[FatJets.jetId > 1]
                
        # jet kinematics # 
        jetkincut = (FatJets.pt > self.ak8PtMin) & (np.abs(getRapidity(FatJets.p4)) < 2.4)
        
        selection.add('jetkincut', ak.any(jetkincut, axis=1))
        FatJets = FatJets[jetkincut]
        del jetkincut

  
        
        # at least 2 ak8 jets #
        selection.add('twoFatJets', (ak.num(FatJets) > 1))
        

        # event cuts #
        
        testCut = (ak.sum(Jets[(Jets.pt>30) & (np.abs(Jets.eta)<3.0)].pt, axis=1) > self.htCut) & \
                  (ak.any((FatJets.jetId > 1), axis=1)) & \
                  (ak.any((FatJets.pt > self.ak8PtMin) & (np.abs(getRapidity(FatJets.p4)) < 2.4))) & \
                  ((ak.num(FatJets) > 1))
                
        # save cutflow
        if isNominal:
            cuts = []
            for cut in selection.names:
                cuts.append(cut)
                output['cutflow'][cut] += len(FatJets[selection.all(*cuts)])
            del cuts
            
        jetpt_nocut = FatJets[selection.all('twoFatJets','jetid')].pt[:,0]
        jety_nocut = getRapidity(FatJets[selection.all('twoFatJets','jetid')].p4)[:,0]
        jeteta_nocut = FatJets[selection.all('twoFatJets','jetid')].eta[:,0]
        jetphi_nocut = FatJets[selection.all('twoFatJets','jetid')].eta[:,0]
        
        


        eventCut = selection.all(*selection.names)
                            
        FatJets = FatJets[eventCut]
        SubJets = SubJets[eventCut]
        Jets    = Jets[eventCut]
        Met = Met[eventCut]
        evtweights = evtweights[eventCut]
        events = events[eventCut]
        
        # if event cut results in few events
        if (len(events) < 10): return output

        if not isData: GenJets = GenJets[eventCut]
            


    
        logger.debug('JEC:%s:ttbar cand JES:%s:%s', time.time(), FatJets.pt, correction)    





        # sort jets by pt to select two leading jets
        FatJet_pt_argsort = ak.argsort(FatJets.pt, ascending=False) 
        SortedFatJets = FatJets[FatJet_pt_argsort]

        # higher deepak8 discriminator will be used for jet in mt of mt vs mtt distribution
        jet0 = ak.where(SortedFatJets[:,0].deepTagMD_TvsQCD > SortedFatJets[:,1].deepTagMD_TvsQCD,
                            SortedFatJets[:,0],
                            SortedFatJets[:,1]
                                   )
            
        jet1 = ak.where(SortedFatJets[:,0].deepTagMD_TvsQCD > SortedFatJets[:,1].deepTagMD_TvsQCD,
                            SortedFatJets[:,1],
                            SortedFatJets[:,0]
                                   )

        mcut_s0 = ((self.minMSD < jet0.msoftdrop) & (jet0.msoftdrop < self.maxMSD) )
        mcut_s1 = ((self.minMSD < jet1.msoftdrop) & (jet1.msoftdrop < self.maxMSD) )



        logger.debug('SortedFatJets:%s:FatJets.pt:%s:%s', time.time(), FatJets.pt, correction)
        logger.debug('SortedFatJets:%s:SortedFatJets.pt:%s:%s', time.time(), SortedFatJets.pt, correction)
        logger.debug('SortedFatJets:%s:SortedFatJets.deepTagMD_TvsQCD:%s:%s', time.time(), SortedFatJets.deepTagMD_TvsQCD, correction)
        logger.debug('SortedFatJets:%s:jet0.pt:%s:%s', time.time(), jet0.pt, correction)
        logger.debug('SortedFatJets:%s:jet1.pt:%s:%s', time.time(), jet1.pt, correction)


        del FatJet_pt_argsort, SortedFatJets

        # signal = pass region for 2DAlphabet
        # both jets pass deepak8 tagger
        ttag_s0 = (jet0.deepTagMD_TvsQCD > self.deepAK8disc)
        ttag_s1 = (jet1.deepTagMD_TvsQCD > self.deepAK8disc) & (mcut_s1)

            
        # antitag = fail region for 2DAlphabet
        # leading (in deepak8 disc) jet passes deepak8 tagger
        # subleading (in deepak8 disc) jet fails deepak8 tagger         
        antitag_disc = ((jet1.deepTagMD_TvsQCD < self.deepAK8disc) & (jet1.deepTagMD_TvsQCD > self.deepAK8low))
        antitag = (antitag_disc) & (ttag_s0) & (mcut_s1)


            

                

        
        # ---- Apply Delta Phi Cut for Back to Back Topology ---- #
        dPhiCut = (np.abs(jet0.p4.delta_phi(jet1.p4)) > 2.1)

        
        
        # ttbar candidates have 2 subjets #
        hasSubjets0 = ((jet0.subJetIdx1 > -1) & (jet0.subJetIdx2 > -1))
        hasSubjets1 = ((jet1.subJetIdx1 > -1) & (jet1.subJetIdx2 > -1))
        GoodSubjets = ((hasSubjets0) & (hasSubjets1))
        
        # apply ttbar event cuts #
        if isNominal:
            output['cutflow']['dPhiCut'] += len(FatJets[(dPhiCut)])
            output['cutflow']['Good Subjets'] += len(FatJets[(dPhiCut & GoodSubjets)])

        ttbarcandCuts = (dPhiCut & GoodSubjets)
        antitag = antitag[ttbarcandCuts]
        ttag_s0 = ttag_s0[ttbarcandCuts]
        ttag_s1 = ttag_s1[ttbarcandCuts]
        jet0 = jet0[ttbarcandCuts]
        jet1 = jet1[ttbarcandCuts]
        FatJets = FatJets[ttbarcandCuts]
        Jets = Jets[ttbarcandCuts]
        SubJets = SubJets[ttbarcandCuts]
        events = events[ttbarcandCuts]
        evtweights = evtweights[ttbarcandCuts]
        
        if not isData: GenJets = GenJets[ttbarcandCuts]
        del dPhiCut, ttbarcandCuts, hasSubjets0, hasSubjets1, GoodSubjets
                              
        logger.debug('memory:%s: apply event cuts %s:%s', time.time(), correction, get_memory_usage())

        
        # ttbarmass
        ttbarmass = (jet0.p4 + jet1.p4).mass
        
        # subjets
        SubJet00 = ak.flatten(SubJets[ak.unflatten(jet0.subJetIdx1, np.ones(len(FatJets), dtype='i'))])
        SubJet01 = ak.flatten(SubJets[ak.unflatten(jet0.subJetIdx2, np.ones(len(FatJets), dtype='i'))])
        SubJet10 = ak.flatten(SubJets[ak.unflatten(jet1.subJetIdx1, np.ones(len(FatJets), dtype='i'))])
        SubJet11 = ak.flatten(SubJets[ak.unflatten(jet1.subJetIdx2, np.ones(len(FatJets), dtype='i'))])

        # print('\n----subjets----')

        # print('SubJet00.pt', SubJet00.pt, events.event[0], correction)
        # print('SubJet01.pt', SubJet01.pt, events.event[0], correction)
        # print('SubJet10.pt', SubJet10.pt, events.event[0], correction)
        # print('SubJet11.pt', SubJet11.pt, events.event[0], correction)
        # print('SubJet00.btagDeepB', SubJet00.btagDeepB, events.event[0], correction)
        # print('SubJet01.btagDeepB', SubJet01.btagDeepB, events.event[0], correction)
        # print('SubJet10.btagDeepB', SubJet10.btagDeepB, events.event[0], correction)
        # print('SubJet11.btagDeepB', SubJet11.btagDeepB, events.event[0], correction)
        
        # print('----subjets----\n')
    
        
        

        # discriminators for plotting
        tau32_s0 = np.where(jet0.tau2>0,jet0.tau3/jet0.tau2, 0 )
        tau32_s1 = np.where(jet1.tau2>0,jet1.tau3/jet1.tau2, 0 )

        taucut_s0 = (tau32_s0 < self.tau32Cut)
        taucut_s1 = (tau32_s1 < self.tau32Cut)
        
        bdisc_s0 = np.maximum(SubJet00.btagDeepB , SubJet01.btagDeepB)
        bdisc_s1 = np.maximum(SubJet10.btagDeepB , SubJet11.btagDeepB)
        
        
        tdisc_s0 = jet0.deepTagMD_TvsQCD
        tdisc_s1 = jet1.deepTagMD_TvsQCD
        
        
        
        # ---- Define "Top Tag" Regions ---- #
        antitag_probe = np.logical_and(antitag, ttag_s1) # Found an antitag and ttagged probe pair for mistag rate (AT&Pt)
        pretag =  (ttag_s0)                    # Only jet0 (pret)
        ttag0 =   ((~ttag_s0) & (~ttag_s1))    # No tops tagged (0t)
        ttag1 =   (ttag_s0 ^ ttag_s1)          # Exclusively one top tagged (1t)
        ttagI =   (ttag_s0 | ttag_s1)          # At least one top tagged ('I' for 'inclusive' tagger; >=1t; 1t+2t)
        ttag2 =   (ttag_s0 & ttag_s1)          # Both jets top tagged (2t)
        Alltags = (ttag0 | ttagI)              # Either no tag or at least one tag (0t+1t+2t)
        

        # rapidity #
        rapidity = getRapidity(jet0.p4) - getRapidity(jet1.p4)
        cen = (np.abs(rapidity) < 1.0)
        fwd = (~cen)





    
        # rapidity, btag and top tag categories
        regs = {'cen': cen, 'fwd': fwd}
        #btags = {'0b': btag0, '1b':btag1, '2b':btag2}
        ttags = {
            "at":antitag, # 2Dalphabet fail region
            "2t":ttag2, # 2Dalphabet pass region
            # "AT&Pt": antitag_probe, 
            #  "pret":pretag, 
            #  "0t":ttag0, 
            #  "1t":ttag1, 
            #  ">=1t":ttagI, 
            #  ">=0t":Alltags
                }
        
        
        
        # get all analysis category masks
        categories = { t[0]+y[0] : (t[1]&y[1])  for t,y in itertools.product( ttags.items(), 
                                                                        regs.items())
            }
        
        # use subset of analysis category masks from ttbaranalysis.py
        labels_and_categories = {label:categories[label] for label in self.anacats}
    
    
    
        logger.debug('memory:%s: get analysis categories %s:%s', time.time(), correction, get_memory_usage())





        
        
        self.weights[correction] = Weights(len(evtweights))
        
        self.weights[correction].add('genWeight', evtweights)
                        
        
        # kinematics variables for plotting
        jetpt = jet0.p4.pt
        jeteta = jet0.p4.eta
        jety = getRapidity(jet0.p4)
        jetphi = jet0.p4.phi
        jetmass = jet0.p4.mass
        jetp = jet0.p4.p

        jetpt1 = jet1.p4.pt
        jeteta1 = jet1.p4.eta
        jety1 = getRapidity(jet1.p4)
        jetphi1 = jet1.p4.phi
        jetmass1 = jet1.p4.mass
        jetp1 = jet1.p4.p
        
        # plot same jetmass as pre-tagged, anti-tagged jet
        jetmsd = jet0.msoftdrop
        jetmsd1 = jet1.msoftdrop

        
        # values for mistag rate calculation #
        numerator = np.where(antitag_probe, jet1.p4.p, -1)
        denominator = np.where(antitag, jet1.p4.p, -1)
        
        # pt reweighting #
        if ('TTbar' in dataset):
             ttbar_wgt = pTReweighting(jet0.pt, jet1.pt)
             self.weights[correction].add('ptReweighting', ttbar_wgt)
                 
        if (not self.noSyst) and (not isData):
                    
            if 'pileup' in self.systematics:
                
                puNom, puUp, puDown = GetPUSF(events, self.iov)
                self.weights[correction].add("pileup", 
                    weight=puNom, 
                    weightUp=puUp, 
                    weightDown=puDown,
                           )
                
            logger.debug('memory:%s: pileup systematics %s:%s', time.time(), correction, get_memory_usage())

            if ('prefiring' in self.systematics) and ("L1PreFiringWeight" in events.fields):
                if ('2016' in self.iov) or ('2017' in self.iov):
                
                    prefiringNom, prefiringUp, prefiringDown = GetL1PreFiringWeight(events)
                    self.weights[correction].add("prefiring", 
                        weight=prefiringNom, 
                        weightUp=prefiringUp, 
                        weightDown=prefiringDown,
                               )
                    
            logger.debug('memory:%s: prefiring systematics %s:%s', time.time(), correction, get_memory_usage())
                    
            if 'pdf' in self.systematics:
                
                pdfNom, pdfUp, pdfDown  = GetPDFWeights(events)
                self.weights[correction].add("pdf", 
                    weight=pdfNom, 
                    weightUp=pdfUp, 
                    weightDown=pdfDown,
                           )    
                
            logger.debug('memory:%s: pdf systematics %s:%s', time.time(), correction, get_memory_usage())
            
            if 'q2' in self.systematics:
                
                q2Nom, q2Up, q2Down = GetQ2weights(events)
                
                self.weights[correction].add("q2", 
                    weight=q2Nom, 
                    weightUp=q2Up, 
                    weightDown=q2Down,
                           )   
                
                
            logger.debug('memory:%s: q2 systematics %s:%s', time.time(), correction, get_memory_usage())

                


            if 'ttag_pt1' in self.systematics:

                ptbins = [0,400,480,600]


  

                ttag_scale_factors = {
                                        '2016':{
                                            'tight':  {'nominal': [0.92,1.01,0.84,1.00], 'up': [1.04,1.18,0.9,1.07], 'down': [0.82,0.84,0.78,0.94]},
                                            'medium': {'nominal': [0.89,1.02,0.93,1.00], 'up': [0.97,1.07,0.97,1.05], 'down': [0.81,0.97,0.89,0.95]},
                                            'loose':  {'nominal': [0.98,0.95,0.97,1.00], 'up': [0.98,0.99,1.0,1.04], 'down': [0.91,0.91,0.94,0.96]}  
                                        },
                                        '2016APV':{
                                            'tight':  {'nominal': [0.92,1.01,0.84,1.00], 'up': [1.04,1.18,0.9,1.07], 'down': [0.82,0.84,0.78,0.94]},
                                            'medium': {'nominal': [0.89,1.02,0.93,1.00], 'up': [0.97,1.07,0.97,1.05], 'down': [0.81,0.97,0.89,0.95]},
                                            'loose':  {'nominal': [0.98,0.95,0.97,1.00], 'up': [0.98,0.99,1.0,1.04], 'down': [0.91,0.91,0.94,0.96]}  
                                        },
                                        '2017':{
                                            'tight':  {'nominal': [0.88,0.9,0.95,0.97], 'up': [0.96,0.95,1.0,1.03], 'down': [0.8,0.85,0.9,0.91]},
                                            'medium': {'nominal': [0.95,1.0,0.98,0.98], 'up': [1.01,1.04,1.02,1.02], 'down': [0.89,0.96,0.94,0.94]},
                                            'loose':  {'nominal': [0.95,0.98,0.97,0.97], 'up': [1.0,1.01,1.0,1.01], 'down': [0.9,0.95,0.94,0.93]}  
                                        },
                                        '2018':{
                                            'tight':  {'nominal': [0.81,0.93,0.96,0.93], 'up': [0.88,0.98,1.02,1.01], 'down': [0.74,0.88,0.92,0.85]},
                                            'medium': {'nominal': [0.90,0.97,0.98,0.95], 'up': [0.95,1.0,1.01,0.98], 'down': [0.85,0.94,0.95,0.92]},
                                            'loose':  {'nominal': [0.96,1.00,0.98,0.99], 'up': [1.0,1.03,1.0,1.02], 'down': [0.92,0.97,0.96,0.96]}
                                        }   
                                    }
                
                nomsf = np.array(ttag_scale_factors[self.iov][self.deepAK8Cut]['nominal'])
                upsf = np.array(ttag_scale_factors[self.iov][self.deepAK8Cut]['up'])
                downsf = np.array(ttag_scale_factors[self.iov][self.deepAK8Cut]['down'])
                
                nomsf_fail = np.array(ttag_scale_factors[self.iov]['loose']['nominal'])
                upsf_fail = np.array(ttag_scale_factors[self.iov]['loose']['up'])
                downsf_fail = np.array(ttag_scale_factors[self.iov]['loose']['down'])

                jet0_ptbins = np.digitize(ak.to_numpy(jetpt), ptbins) - 1
                jet1_ptbins = np.digitize(ak.to_numpy(jetpt1), ptbins) - 1


                ttagSFNom = np.ones(len(events))
                ttagSFUp = np.ones(len(events))
                ttagSFDown = np.ones(len(events))

    


                ttagSFNom_1 = np.where(jet0_ptbins==1, nomsf[1], 1.0) * np.where((jet1_ptbins==1 & ttag2), nomsf[1], 1.0) * np.where((jet1_ptbins==1 & antitag), nomsf_fail[1], 1.0)
                ttagSFUp_1 = np.where(jet0_ptbins==1, upsf[1], 1.0) * np.where((jet1_ptbins==1 & ttag2), upsf[1], 1.0) * np.where((jet1_ptbins==1 & antitag), upsf_fail[1], 1.0)
                ttagSFDown_1 = np.where(jet0_ptbins==1, downsf[1], 1.0) * np.where((jet1_ptbins==1 & ttag2), downsf[1], 1.0) * np.where((jet1_ptbins==1 & antitag), downsf_fail[1], 1.0)

                ttagSFNom_2 = np.where(jet0_ptbins==2, nomsf[2], 1.0) * np.where((jet1_ptbins==2 & ttag2), nomsf[2], 1.0) * np.where((jet1_ptbins==2 & antitag), nomsf_fail[2], 1.0)
                ttagSFUp_2 = np.where(jet0_ptbins==2, upsf[2], 1.0) * np.where((jet1_ptbins==2 & ttag2), upsf[2], 1.0) * np.where((jet1_ptbins==2 & antitag), upsf_fail[2], 1.0)
                ttagSFDown_2 = np.where(jet0_ptbins==2, downsf[2], 1.0) * np.where((jet1_ptbins==2 & ttag2), downsf[2], 1.0) * np.where((jet1_ptbins==2 & antitag), downsf_fail[2], 1.0)

                ttagSFNom_3 = np.where(jet0_ptbins==3, nomsf[3], 1.0) * np.where((jet1_ptbins==3 & ttag2), nomsf[3], 1.0) * np.where((jet1_ptbins==3 & antitag), nomsf_fail[3], 1.0)
                ttagSFUp_3 = np.where(jet0_ptbins==3, upsf[3], 1.0) * np.where((jet1_ptbins==3 & ttag2), upsf[3], 1.0) * np.where((jet1_ptbins==3 & antitag), upsf_fail[3], 1.0)
                ttagSFDown_3 = np.where(jet0_ptbins==3, downsf[3], 1.0) * np.where((jet1_ptbins==3 & ttag2), downsf[3], 1.0) * np.where((jet1_ptbins==3 & antitag), downsf_fail[3], 1.0)


              
            
                self.weights[correction].add("ttag_pt1", 
                    weight=ttagSFNom_1, 
                    weightUp=ttagSFUp_1, 
                    weightDown=ttagSFDown_1,
                           )

                self.weights[correction].add("ttag_pt2", 
                    weight=ttagSFNom_2, 
                    weightUp=ttagSFUp_2, 
                    weightDown=ttagSFDown_2,
                           )
                
                self.weights[correction].add("ttag_pt3", 
                    weight=ttagSFNom_3, 
                    weightUp=ttagSFUp_3, 
                    weightDown=ttagSFDown_3,
                           )
            
             
                            
                            
        logger.debug('JEC:%s:histogram JES:%s:%s', time.time(), FatJets.pt, correction)
        
        
        for i, [ilabel,icat] in enumerate(labels_and_categories.items()):
    
            # final cutflow per analysis category
            
            if isNominal:
                output['cutflow'][ilabel] += len(events.event[icat])
                                                           
            
            output['jetmass'].fill(
                                   systematic=correction,
                                   anacat = i,
                                   jetmass = jetmass[icat],
                                   weight = self.weights[correction].weight()[icat],
                                  )
            output['jetmsd'].fill(
                                   systematic=correction,
                                   anacat = i,
                                   jetmass = jetmsd[icat],
                                   weight = self.weights[correction].weight()[icat],
                                  )
            output['jetmass1'].fill(
                                   systematic=correction,
                                   anacat = i,
                                   jetmass = jetmass1[icat],
                                   weight = self.weights[correction].weight()[icat],
                                  )
            output['jetmsd1'].fill(
                                   systematic=correction,
                                   anacat = i,
                                   jetmass = jetmsd1[icat],
                                   weight = self.weights[correction].weight()[icat],
                                  )

            
            
            output['mtt_vs_mt'].fill(
                                     systematic=correction,
                                     anacat = i,
                                     jetmass = jetmsd[icat],
                                     ttbarmass = ttbarmass[icat],
                                     weight = self.weights[correction].weight()[icat],
                                    )

            
            output['ttbarmass'].fill(systematic=correction,
                                         anacat = i,
                                         ttbarmass = ttbarmass[icat],
                                         weight = self.weights[correction].weight()[icat],
                                        )
            
            output['mtt_unwgt'].fill(systematic=correction,
                                         anacat = i,
                                         ttbarmass = ttbarmass[icat],
                                         weight = np.ones_like(self.weights[correction].weight()[icat]),
                                        )
            
            
            
            
            
            
            # save weights
            
            output['weights'][correction] += np.sum(self.weights[correction].weight())
            output['systematics'][correction] += len(events.event[icat])

            if isNominal:  

                for syst in self.weights[correction].variations:
                    
                    
                    output['weights'][syst] += np.sum(self.weights[correction].weight(syst))
                    output['systematics'][syst] += len(events.event[icat])
                    
                    output['jetmass'].fill(
                                   systematic=syst,
                                   anacat = i,
                                   jetmass = jetmass[icat],
                                   weight = self.weights[correction].weight(syst)[icat],
                                  )
                    output['jetmsd'].fill(
                                           systematic=syst,
                                           anacat = i,
                                           jetmass = jetmsd[icat],
                                           weight = self.weights[correction].weight(syst)[icat],
                                          )
                    output['jetmass1'].fill(
                                   systematic=syst,
                                   anacat = i,
                                   jetmass = jetmass1[icat],
                                   weight = self.weights[correction].weight(syst)[icat],
                                  )
                    output['jetmsd1'].fill(
                                           systematic=syst,
                                           anacat = i,
                                           jetmass = jetmsd1[icat],
                                           weight = self.weights[correction].weight(syst)[icat],
                                          )
                    

                    output['ttbarmass'].fill(systematic=syst,
                                         anacat = i,
                                         ttbarmass = ttbarmass[icat],
                                         weight = self.weights[correction].weight(syst)[icat],
                                        )
                    output['mtt_unwgt'].fill(systematic=syst,
                                         anacat = i,
                                         ttbarmass = ttbarmass[icat],
                                         weight = np.ones_like(self.weights[correction].weight(syst)[icat]),
                                        )

                    output['mtt_vs_mt'].fill(systematic=syst,
                                         anacat = i,
                                         ttbarmass = ttbarmass[icat],
                                         jetmass = jetmsd[icat],
                                         weight = self.weights[correction].weight(syst)[icat],
                                        )




                    
                    
                    
        logger.debug('memory:%s: fill histograms %s:%s', time.time(), correction, get_memory_usage())
                    
                    
                    
                    


        

        return output

    def postprocess(self, accumulator):
        logger.debug('memory:%s: finish processor:%s', time.time(), get_memory_usage())
        return accumulator
        
        
        

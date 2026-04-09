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
from coffea.analysis_tools import PackedSelection
from collections import defaultdict
import sys
import os, psutil
import copy
import scipy.stats as ss
import numpy as np
import pandas as pd
from numpy.random import RandomState
import random
import correctionlib
import json
import logging
import psutil
import time

import awkward as ak

# for dask, `from python.corrections import` does not work
sys.path.append(os.getcwd()+'/python/')

from corrections import (
    GetFlavorEfficiency,
    getLumiMask,
    getMETFilter,
)
from btagCorrections import btagCorrections
from functions import getRapidity
from categories import build_analysis_categories
from jets import Run3JetManager
from hists import build_output_histograms
from weights import Run3WeightManager
from truthstudy import truthstudy_counts, build_gen_top_match_info, build_top_aligned_genjetak8_match_info



# logfile = 'coffea_' + str(int(time.time())) + '.log'
# print(logfile)
# logging.basicConfig(filename=logfile, level=logging.DEBUG)
logger = logging.getLogger('__main__')
logger.setLevel(logging.DEBUG)


#ak.behavior.update(candidate.behavior)
ak.behavior.update(vector.behavior)


# ...existing code...
def get_memory_usage(human_readable=True, precision=2):
    process = psutil.Process(os.getpid())
    memory_usage_bytes = process.memory_info().rss

    if not human_readable:
        # return MB as before if numeric value is needed
        return memory_usage_bytes / (1024 * 1024)

    units = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    size = float(memory_usage_bytes)
    unit_index = 0
    while size >= 1024 and unit_index < len(units) - 1:
        size /= 1024.0
        unit_index += 1

    return f"{size:.{precision}f} {units[unit_index]}"


class Logger:
    DEBUG = 10
    INFO = 20

    def __init__(self, mode='debug'):
        self.level = self.DEBUG if mode == 'debug' else self.INFO

    def debug(self, msg, *args):
        if self.level <= self.DEBUG:
            print('[DEBUG]', msg % args if args else msg)

    def info(self, msg, *args):
        if self.level <= self.INFO:
            print('[INFO]', msg % args if args else msg)


def update(events, collections):
    # https://github.com/nsmith-/boostedhiggs/blob/master/boostedhiggs/hbbprocessor.py
    """Return a shallow copy of events array with some collections swapped out"""
    out = events
#     logger.debug('update:%s:%s', time.time(), collections)
    
    for name, value in collections.items():
        out = ak.with_field(out, value, name)

    return out


"""Skimmer Class to apply event selections and store needed variables."""
class TTbarResProcessor(processor.ProcessorABC):
    def __init__(self,
                 htCut=900.,
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
                 systematics = ['nominal', 'pileup', 'pdf', 'q2', "ttag_pt1"],
                 anacats = ['2t0bcen'],
                 debug = False
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
        self.debug = debug  
        #self.rpf_params = rpf_params        
    
        # from https://twiki.cern.ch/twiki/bin/view/CMS/DeepAK8Tagging2018WPsSFs#2016_Data
        deepak8cuts = {
            'loose':{ # 1%
                '2022': 0.435, 
                '2023':    0.435,
                '2024':    0.344,
                '2025':    0.470,
            },
            'medium':{ # 0.5%
                '2022': 0.632, 
                '2023':    0.632,
                '2024':    0.554,
                '2025':    0.685,
            },
            'tight': { # 0.1%
                '2016APV': 0.889, 
                '2016':    0.889,
                '2017':    0.863,
                '2018':    0.920,
            } 
        }
    
        
        self.logger = Logger(mode='debug' if debug else 'info')

        self.weights = {}
        self.jet_manager = Run3JetManager(
            iov=self.iov,
            systematics=self.systematics,
            no_syst=self.noSyst,
            ak8_pt_min=self.ak8PtMin,
            ht_cut=self.htCut,
        )
        self.weight_manager = Run3WeightManager(
            iov=self.iov,
            systematics=self.systematics,
            no_syst=self.noSyst,
            deepak8_cut=self.deepAK8Cut,
        )
    
        
        
        self.deepAK8disc = deepak8cuts[deepAK8Cut][self.iov]
        
        if deepAK8Cut == 'tight':
            self.deepAK8low = deepak8cuts['medium'][self.iov]
        elif deepAK8Cut == 'medium':
             self.deepAK8low = deepak8cuts['loose'][self.iov]
        else:
            self.deepAK8low = 0.2
        
        
        
        
        
        # analysis categories #
        self.anacats = anacats
        self.label_dict = {i: label for i, label in enumerate(self.anacats)}
        self.label_to_int_dict = {label: i for i, label in enumerate(self.anacats)}

        
        # output histograms/accumulators (organized in a ROOT-like grouped spec in hists_run3.py)
        self.histo_dict = build_output_histograms(
            anacats=self.anacats,
            systematics=self.systematics,
            no_syst=self.noSyst,
        )
        
      

        
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

        isData = ('data' in events.metadata['dataset']) or ('SingleMu' in events.metadata['dataset'])
        corrections = self.jet_manager.build_corrections(events, isData)
        
        if corrections is None:
            return self.process_analysis(events, 'nominal', nEvents)


        # loop through corrections
        outputs = []
        for collections, name in corrections:
            outputs.append(self.process_analysis(update(events, collections), name, nEvents))


        return processor.accumulate(outputs)
     


    def process_analysis(self, events, correction, nEvents):
        
        dataset = events.metadata['dataset']
        filename = events.metadata['filename']
        
        self.logger.debug('start processor: correction=%s, memory=%s', correction, get_memory_usage())
        

                
        isNominal = (correction=='nominal')
        isData = ('data' in dataset) or ('SingleMu' in dataset)
    
        output = self.histo_dict 
        
        if isNominal:
            output['cutflow']['all events 1'] += nEvents
        
        
        # lumi mask #
        if (isData):
            
            lumi_mask = np.array(getLumiMask(self.iov)(events.run, events.luminosityBlock), dtype=bool)
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

        "2022": ["PFHT1050"],
        "2023" :   ["PFHT1050"],
        "2024" :   ["PFHT1050"],
        "2025" :   ["PFHT1050"],

        }
        '''
        if "HLT" in events.fields:
         print("HLT paths in this file:")
         for hlt_name in events.HLT.fields:
            print(hlt_name)
        '''
        try:
            selection.add('trigger', (events.HLT[triggernames[self.iov][0]] | events.HLT[triggernames[self.iov][1]]) )
        except:
            selection.add('trigger', (events.HLT[triggernames[self.iov][0]]))


        # objects #
        FatJets, SubJets, Jets, GenJets, GenJetAK8, SubGenJetAK8 = self.jet_manager.prepare_analysis_objects(events, isData)
        #Met = events.MET
        run = events.run.to_numpy()
        lumi = events.luminosityBlock.to_numpy()
        evt = events.event.to_numpy()
        '''
        print("\n--- FatJet variables ---")
        for var in events.FatJet.fields:
            print("FatJet_" + var)

        print("\n--- SubJet variables ---")
        for var in events.SubJet.fields:
            print("SubJet_" + var)
        '''
        logger.debug('memory:%s: get nanoAOD objects %s:%s', time.time(), correction, get_memory_usage())

        
        
        
        # ---- Get event weights from dataset ---- #

        # if blinding + trigger results in too few events
        if len(events) < 10: return output
        

                
        
        if isData:
            evtweights = np.ones(len(events))
        else:
            if "LHEWeight_originalXWGTUP" not in events.fields: 
                evtweights = events.genWeight
            else: 
                evtweights = events.LHEWeight_originalXWGTUP


        if isNominal:
            output['cutflow']['all events'] += len(FatJets)
            output['cutflow']['sumw'] += np.sum(evtweights)
            output['cutflow']['sumw2'] += np.sum(evtweights**2)

            
        
            
        
        # ---- event selection and object selection ---- #


        FatJets, jet_masks = self.jet_manager.baseline_masks(events, FatJets, Jets)

        # ht cut #
        selection.add('htCut', jet_masks['htCut'])

        
        # met filters #
        selection.add('metfilter', getMETFilter(self.iov, events))
                
        # jet id #
        #selection.add('jetid', jet_masks['jetid'])
                
        # jet kinematics # 
        selection.add('jetkincut', jet_masks['jetkincut'])

  
        
        # at least 2 ak8 jets #
        selection.add('twoFatJets', jet_masks['twoFatJets'])
        
        # event cuts #
        ''' 
        # save cutflow
        if isNominal:
            cuts = []
            for cut in selection.names:
                cuts.append(cut)
                output['cutflow'][cut] += len(FatJets[selection.all(*cuts)])
            del cuts
        '''    
        eventCut = selection.all(*selection.names)
                            
        FatJets = FatJets[eventCut]
        SubJets = SubJets[eventCut]
        Jets    = Jets[eventCut]
        #Met = Met[eventCut]
        evtweights = evtweights[eventCut]
        events = events[eventCut]
        run = run[eventCut]
        lumi = lumi[eventCut]
        evt = evt[eventCut]
        # if event cut results in few events
        logger.debug(f"Length of event {len(events)}")
        if (len(events) < 10): return output

        if not isData:
            GenJets = GenJets[eventCut]
            if GenJetAK8 is not None:
                GenJetAK8 = GenJetAK8[eventCut]
            if SubGenJetAK8 is not None:
                SubGenJetAK8 = SubGenJetAK8[eventCut]
            
        ##Add GenTruth study
        if isNominal and (not isData):
            truth_counts = truthstudy_counts(
                genparts=events.GenPart,
                fatjets=FatJets,
                subJets=SubJets,
                jets=Jets,
                dr_ak8=0.8,
                dr_ak4=1.2,
            )
            output["truthstudy"]["n_hadtop"] += truth_counts["n_hadtop"]
            output["truthstudy"]["n_hadtop_ak8"] += truth_counts["n_hadtop_ak8"]
            output["truthstudy"]["n_hadtop_ak8_ak4"] += truth_counts["n_hadtop_ak8_ak4"]

        
  

    
        #logger.debug('JEC:%s:ttbar cand JES:%s:%s', time.time(), FatJets.pt, correction)    

        # sort jets by pt to select two leading jets
       
        FatJet_pt_argsort = ak.argsort(FatJets.pt, ascending=False) 
        SortedFatJets = FatJets[FatJet_pt_argsort]
        
        # higher deepak8 discriminator will be used for jet in mt of mt vs mtt distribution
        if (self.iov == '2023'):
            jet0 = ak.where(SortedFatJets[:,0].particleNet_XttVsQCD > SortedFatJets[:,1].particleNet_XttVsQCD,
                            SortedFatJets[:,0],
                            SortedFatJets[:,1]
                                   )
            
            jet1 = ak.where(SortedFatJets[:,0].particleNet_XttVsQCD > SortedFatJets[:,1].particleNet_XttVsQCD,
                            SortedFatJets[:,1],
                            SortedFatJets[:,0]
                                )
        elif (self.iov == '2024'):

              jet0 = ak.where(SortedFatJets[:,0].globalParT3_TopbWqq > SortedFatJets[:,1].globalParT3_TopbWqq,
                            SortedFatJets[:,0],
                            SortedFatJets[:,1]
               )
            
              jet1 = ak.where(SortedFatJets[:,0].globalParT3_TopbWqq > SortedFatJets[:,1].globalParT3_TopbWqq,
                            SortedFatJets[:,1],
                            SortedFatJets[:,0]

              )
        mcut_s0 = ((self.minMSD < jet0.msoftdrop) & (jet0.msoftdrop < self.maxMSD) )
        mcut_s1 = ((self.minMSD < jet1.msoftdrop) & (jet1.msoftdrop < self.maxMSD) )
        del FatJet_pt_argsort, SortedFatJets


        '''
        logger.debug('SortedFatJets:%s:FatJets.pt:%s:%s', time.time(), FatJets.pt, correction)
        logger.debug('SortedFatJets:%s:SortedFatJets.pt:%s:%s', time.time(), SortedFatJets.pt, correction)
        logger.debug('SortedFatJets:%s:SortedFatJets.deepTagMD_TvsQCD:%s:%s', time.time(), SortedFatJets.particleNet_XttVsQCD, correction)
        logger.debug('SortedFatJets:%s:jet0.pt:%s:%s', time.time(), jet0.pt, correction)
        logger.debug('SortedFatJets:%s:jet1.pt:%s:%s', time.time(), jet1.pt, correction)
        '''

        #del FatJet_pt_argsort, SortedFatJets

        # signal = pass region for 2DAlphabet
        # both jets pass deepak8 tagger
        if (self.iov == '2023'):
        
            ttag_s0 = (jet0.particleNet_XttVsQCD > self.deepAK8disc)
            ttag_s1 = (jet1.particleNet_XttVsQCD > self.deepAK8disc) & (mcut_s1)
            ttag_s0_1 = (jet0.particleNet_XttVsQCD > self.deepAK8disc)
            ttag_s1_1 = (jet1.particleNet_XttVsQCD > self.deepAK8disc) & (mcut_s1)
            
            # antitag = fail region for 2DAlphabet
            # leading (in deepak8 disc) jet passes deepak8 tagger
            # subleading (in deepak8 disc) jet fails deepak8 tagger         
            antitag_disc = ((jet1.particleNet_XttVsQCD < self.deepAK8disc) & (jet1.particleNet_XttVsQCD > self.deepAK8low))
            
        elif (self.iov == '2024'):
             
            ttag_s0 = (jet0.globalParT3_TopbWqq > self.deepAK8disc)
            ttag_s1 = (jet1.globalParT3_TopbWqq > self.deepAK8disc) & (mcut_s1)
            ttag_s0_1 = (jet0.globalParT3_TopbWqq > self.deepAK8disc)
            ttag_s1_1 = (jet1.globalParT3_TopbWqq > self.deepAK8disc) & (mcut_s1)
            
            # antitag = fail region for 2DAlphabet
            # leading (in deepak8 disc) jet passes deepak8 tagger
            # subleading (in deepak8 disc) jet fails deepak8 tagger         
            antitag_disc = ((jet1.globalParT3_TopbWqq < self.deepAK8disc) & (jet1.globalParT3_TopbWqq > self.deepAK8low))
        
            antitag = (antitag_disc) & (ttag_s0) & (mcut_s1)


            

                

        
        # ---- Apply Delta Phi Cut for Back to Back Topology ---- #
        dPhiCut = (np.abs(jet0.p4.delta_phi(jet1.p4)) > 2.1)

        
        
        # ttbar candidates have 2 subjets #
        hasSubjets0 = ((jet0.subJetIdx1 > -1) & (jet0.subJetIdx2 > -1))
        hasSubjets1 = ((jet1.subJetIdx1 > -1) & (jet1.subJetIdx2 > -1))
        GoodSubjets = ((hasSubjets0) & (hasSubjets1))
        '''       
        # apply ttbar event cuts #
        if isNominal:
            output['cutflow']['dPhiCut'] += len(FatJets[(dPhiCut)])
            output['cutflow']['Good Subjets'] += len(FatJets[(dPhiCut & GoodSubjets)])
        ''' 

        
        

        ttbarcandCuts = (dPhiCut & GoodSubjets)

        #signal_region_cuts = (dPhiCut & GoodSubjets & ttag_s0 & ttag_s1 & mcut_s0 i)
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
        run = run[ttbarcandCuts & ttag_s0_1 & ttag_s1_1]
        lumi = lumi[ttbarcandCuts  &  ttag_s0_1 & ttag_s1_1]
        evt = evt[ttbarcandCuts & ttag_s0_1 & ttag_s1_1]

        

        if isNominal:
          before = len(output["event_list"]["run"])
 
          output["event_list"]["run"]   += list(run)
          output["event_list"]["lumi"]   += list(lumi)
          output["event_list"]["event"]   += list(evt)
          if self.debug:
              after = len(output["event_list"]["run"])
              print(after - before,  len(run), "   if different then is wrong")
              df = pd.DataFrame(output["event_list"])
              print("rows:", len(df))
              print("duplicates:", df.duplicated(["run","lumi","event"]).sum())
       

        if not isData:
            GenJets = GenJets[ttbarcandCuts]
            if GenJetAK8 is not None:
                GenJetAK8 = GenJetAK8[ttbarcandCuts]
            if SubGenJetAK8 is not None:
                SubGenJetAK8 = SubGenJetAK8[ttbarcandCuts]
        del dPhiCut, ttbarcandCuts, hasSubjets0, hasSubjets1, GoodSubjets
                              
        logger.debug('memory:%s: apply event cuts %s:%s', time.time(), correction, get_memory_usage())

        # Find if there is a third AK8 jet and find the delta_R between the third jet and the two leading jets
        third_jet_mask = ak.num(FatJets) > 2
        jet2 = FatJets[third_jet_mask][:,2]
        dR_jet0_jet2 = jet0[third_jet_mask].p4.delta_r(jet2.p4)
        dR_jet1_jet2 = jet1[third_jet_mask].p4.delta_r(jet2.p4)
        ttbarmass = (jet0.p4 + jet1.p4).mass 


        # ttbarmass
        #print ("ttbarmass" , ttbarmass)


        
        

        # discriminators for plotting
        tau32_s0 = np.where(jet0.tau2>0,jet0.tau3/jet0.tau2, 0 )
        tau32_s1 = np.where(jet1.tau2>0,jet1.tau3/jet1.tau2, 0 )

        taucut_s0 = (tau32_s0 < self.tau32Cut)
        taucut_s1 = (tau32_s1 < self.tau32Cut)
        
        #bdisc_s0 = np.maximum(SubJet00.btagDeepB , SubJet01.btagDeepB)
        #bdisc_s1 = np.maximum(SubJet10.btagDeepB , SubJet11.btagDeepB)
        
        if (self.iov == '2023'):
            tdisc_s0 = jet0.particleNet_XttVsQCD
            tdisc_s1 = jet1.particleNet_XttVsQCD
        elif (self.iov == '2024'):
             tdisc_s0 = jet0.globalParT3_TopbWqq
             tdisc_s1 = jet1.globalParT3_TopbWqq
        
        
        rapidity = getRapidity(jet0.p4) - getRapidity(jet1.p4)
        jet0_abs_eta = np.abs(jet0.eta)
        jet1_abs_eta = np.abs(jet1.eta)

        jet0_ak4_pairs = ak.cartesian({"ak8": ak.singletons(jet0), "ak4": Jets}, axis=1, nested=True)
        jet0_ak4_dr = jet0_ak4_pairs["ak8"].p4.delta_r(jet0_ak4_pairs["ak4"].p4)
        jet0_has_nearby_ak4 = ak.to_numpy(
            ak.flatten(ak.any((jet0_ak4_dr > 0.4) & (jet0_ak4_dr < 0.8), axis=2), axis=1)
        )
        del jet0_ak4_pairs, jet0_ak4_dr

        jet1_ak4_pairs = ak.cartesian({"ak8": ak.singletons(jet1), "ak4": Jets}, axis=1, nested=True)
        jet1_ak4_dr = jet1_ak4_pairs["ak8"].p4.delta_r(jet1_ak4_pairs["ak4"].p4)
        jet1_has_nearby_ak4 = ak.to_numpy(
            ak.flatten(ak.any((jet1_ak4_dr > 0.4) & (jet1_ak4_dr < 0.8), axis=2), axis=1)
        )
        del jet1_ak4_pairs, jet1_ak4_dr

        jet0_ak8_pairs = ak.cartesian({"ak8": ak.singletons(jet0), "fat": FatJets}, axis=1, nested=True)
        jet0_ak8_dr = jet0_ak8_pairs["ak8"].p4.delta_r(jet0_ak8_pairs["fat"].p4)
        jet0_has_nearby_ak8 = ak.to_numpy(
            ak.flatten(
                ak.any(
                    (jet0_ak8_dr < 0.8) & (jet0_ak8_dr > 1e-6),
                    axis=2,
                ),
                axis=1,
            )
        )
        del jet0_ak8_pairs, jet0_ak8_dr

        jet1_ak8_pairs = ak.cartesian({"ak8": ak.singletons(jet1), "fat": FatJets}, axis=1, nested=True)
        jet1_ak8_dr = jet1_ak8_pairs["ak8"].p4.delta_r(jet1_ak8_pairs["fat"].p4)
        jet1_has_nearby_ak8 = ak.to_numpy(
            ak.flatten(
                ak.any(
                    (jet1_ak8_dr < 0.8) & (jet1_ak8_dr > 1e-6),
                    axis=2,
                ),
                axis=1,
            )
        )
        del jet1_ak8_pairs, jet1_ak8_dr
        jet0_nearby_label = np.where(
            jet0_has_nearby_ak4,
            "ak4_nearby",
            np.where(jet0_has_nearby_ak8, "ak8_nearby", "no_jet_nearby"),
        )
        jet1_nearby_label = np.where(
            jet1_has_nearby_ak4,
            "ak4_nearby",
            np.where(jet1_has_nearby_ak8, "ak8_nearby", "no_jet_nearby"),
        )

        gen_top_match_info = None
        genjetak8_match_info = None
        if isNominal and (not isData):
            gen_top_match_info = build_gen_top_match_info(
                genparts=events.GenPart,
                jet0=jet0,
                jet1=jet1,
                dr_match=0.8,
            )
            genjetak8_match_info = build_top_aligned_genjetak8_match_info(
                genparts=events.GenPart,
                genjetak8=GenJetAK8,
                subgenjetak8=SubGenJetAK8,
                jet0=jet0,
                jet1=jet1,
                dr_match=0.8,
                top_align_dr=0.8,
            )
            output["truthstudy"]["n_selected_allhad"] += int(ak.sum(gen_top_match_info["event_mask"]))
            output["truthstudy"]["n_selected_jet0_matched"] += int(ak.sum(gen_top_match_info["jet0_is_matched"]))
            output["truthstudy"]["n_selected_jet1_matched"] += int(ak.sum(gen_top_match_info["jet1_is_matched"]))
            output["truthstudy"]["n_selected_both_matched"] += int(ak.sum(gen_top_match_info["both_jets_matched"]))
            if genjetak8_match_info is not None:
                output["truthstudy"]["n_selected_jet0_genak8_matched"] += int(ak.sum(genjetak8_match_info["jet0_is_matched"]))
                output["truthstudy"]["n_selected_jet1_genak8_matched"] += int(ak.sum(genjetak8_match_info["jet1_is_matched"]))
                output["truthstudy"]["n_selected_both_genak8_matched"] += int(ak.sum(genjetak8_match_info["both_jets_matched"]))

        labels_and_categories = build_analysis_categories(
            antitag=antitag,
            ttag_s0=ttag_s0,
            ttag_s1=ttag_s1,
            rapidity=rapidity,
            anacats=self.anacats,
        )
    
    
    
        #logger.debug('memory:%s: get analysis categories %s:%s', time.time(), correction, get_memory_usage())

        
        antitag_probe = np.logical_and(antitag, ttag_s1)
        self.weights[correction] = self.weight_manager.build_weights(
            dataset=dataset,
            events=events,
            evtweights=evtweights,
            is_data=isData,
            jet0=jet0,
            jet1=jet1,
            ttag2=(ttag_s0 & ttag_s1),
            antitag=antitag,
        )
                        
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
        
        #logger.debug('JEC:%s:histogram JES:%s:%s', time.time(), FatJets.pt, correction)
        #self.logger.debug("Labels and categories: %s", labels_and_categories.items())
        
        for i, [ilabel,icat] in enumerate(labels_and_categories.items()):
    
            # final cutflow per analysis category
            
            if isNominal:
                output['cutflow'][ilabel] += len(events.event[icat])

                dR_min_jet2 = ak.where(
                    dR_jet0_jet2 < dR_jet1_jet2,
                    dR_jet0_jet2,
                    dR_jet1_jet2,
                )

                output['dR_min_jet2'].fill(
                    systematic=correction,
                    dr=dR_min_jet2,
                    ttbarmass=ttbarmass[third_jet_mask],
                    anacat = i,
                    weight=self.weights[correction].weight()[third_jet_mask],
                )

            output['jetmsd'].fill(
                                   systematic=correction,
                                   anacat = i,
                                   jetmsd = jetmsd[icat],
                                   weight = self.weights[correction].weight()[icat],
                                  )
            output['ttbarmass'].fill(systematic=correction,
                                         anacat = i,
                                         ttbarmass = ttbarmass[icat],
                                         weight = self.weights[correction].weight()[icat],
                                        )

            if gen_top_match_info is not None:
                truth_cat_mask = icat[gen_top_match_info["event_mask"]]
                truth_event_weights = self.weights[correction].weight()[gen_top_match_info["event_mask"]]
                truth_weights = truth_event_weights[truth_cat_mask]

                output["gen_mt"].fill(
                    systematic=correction,
                    anacat=i,
                    gentopmass=gen_top_match_info["gen_top0"].mass[truth_cat_mask],
                    weight=truth_weights,
                )
                output["gen_mt"].fill(
                    systematic=correction,
                    anacat=i,
                    gentopmass=gen_top_match_info["gen_top1"].mass[truth_cat_mask],
                    weight=truth_weights,
                )
                output["gen_mttbar"].fill(
                    systematic=correction,
                    anacat=i,
                    ttbarmass=gen_top_match_info["top_pair_mass"][truth_cat_mask],
                    weight=truth_weights,
                )
                output["jet0_gen_dr"].fill(
                    systematic=correction,
                    anacat=i,
                    dr=gen_top_match_info["jet0_dr"][truth_cat_mask],
                    weight=truth_weights,
                )
                output["jet1_gen_dr"].fill(
                    systematic=correction,
                    anacat=i,
                    dr=gen_top_match_info["jet1_dr"][truth_cat_mask],
                    weight=truth_weights,
                )
            if genjetak8_match_info is not None:
                genak8_truth_cat_mask = icat[genjetak8_match_info["event_mask"]]
                genak8_event_weights = self.weights[correction].weight()[genjetak8_match_info["event_mask"]]
                genak8_truth_weights = genak8_event_weights[genak8_truth_cat_mask]

                # Fill response for each top-aligned AK8: jet0->jetmsd and jet1->jetmsd1.
                output["gen_jetmsd_reco_jetmsd"].fill(
                    systematic=correction,
                    anacat=i,
                    abs_eta=jet0_abs_eta[genjetak8_match_info["event_mask"]][genak8_truth_cat_mask],
                    jet_nearby=jet0_nearby_label[genjetak8_match_info["event_mask"]][genak8_truth_cat_mask],
                    genjetmass=genjetak8_match_info["jet0_genjet"].mass[genak8_truth_cat_mask],
                    jetmsd=jetmsd[genjetak8_match_info["event_mask"]][genak8_truth_cat_mask],
                    weight=genak8_truth_weights,
                )
                output["gen_jetmsd_reco_jetmsd"].fill(
                    systematic=correction,
                    anacat=i,
                    abs_eta=jet1_abs_eta[genjetak8_match_info["event_mask"]][genak8_truth_cat_mask],
                    jet_nearby=jet1_nearby_label[genjetak8_match_info["event_mask"]][genak8_truth_cat_mask],
                    genjetmass=genjetak8_match_info["jet1_genjet"].mass[genak8_truth_cat_mask],
                    jetmsd=jetmsd1[genjetak8_match_info["event_mask"]][genak8_truth_cat_mask],
                    weight=genak8_truth_weights,
                )

                jet0_genak8_massres = (
                    jetmsd[genjetak8_match_info["event_mask"]] - genjetak8_match_info["jet0_genjet"].mass
                ) / genjetak8_match_info["jet0_genjet"].mass
                jet1_genak8_massres = (
                    jetmsd1[genjetak8_match_info["event_mask"]] - genjetak8_match_info["jet1_genjet"].mass
                ) / genjetak8_match_info["jet1_genjet"].mass
                jet0_genak8_truth_cat_mask = genak8_truth_cat_mask & genjetak8_match_info["jet0_is_matched"]
                jet1_genak8_truth_cat_mask = genak8_truth_cat_mask & genjetak8_match_info["jet1_is_matched"]

                output["jet_mass_resolution"].fill(
                    systematic=correction,
                    anacat=i,
                    abs_eta=jet0_abs_eta[genjetak8_match_info["event_mask"]][jet0_genak8_truth_cat_mask],
                    jet_nearby=jet0_nearby_label[genjetak8_match_info["event_mask"]][jet0_genak8_truth_cat_mask],
                    massres=jet0_genak8_massres[jet0_genak8_truth_cat_mask],
                    weight=genak8_event_weights[jet0_genak8_truth_cat_mask],
                )
                output["jet_mass_resolution"].fill(
                    systematic=correction,
                    anacat=i,
                    abs_eta=jet1_abs_eta[genjetak8_match_info["event_mask"]][jet1_genak8_truth_cat_mask],
                    jet_nearby=jet1_nearby_label[genjetak8_match_info["event_mask"]][jet1_genak8_truth_cat_mask],
                    massres=jet1_genak8_massres[jet1_genak8_truth_cat_mask],
                    weight=genak8_event_weights[jet1_genak8_truth_cat_mask],
                )
               
            '''
            output['jetdy'].fill(
                                  systematic=correction,
                                  anacat = i,
                                  jetdy = rapidity[icat],
                                  weight = self.weights[correction].weight()[icat],
                                  )   
            
            output['jetmass'].fill(
                                   systematic=correction,
                                   anacat = i,
                                   jetmass = jetmass[icat],
                                   weight = self.weights[correction].weight()[icat],
                                  )
            output['jetmsd'].fill(
                                   systematic=correction,
                                   anacat = i,
                                   jetmsd = jetmsd[icat],
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
            
            
            
            
            
            '''
            # save weights
            
            output['weights'][correction] += np.sum(self.weights[correction].weight())
            output['systematics'][correction] += len(events.event[icat])

            if isNominal:  


                for syst in self.weights[correction].variations: # Filling for non jet systematics here
                    output['jetmsd'].fill(
                                systematic=syst,
                                anacat = i,
                                jetmsd = jetmsd[icat],
                                weight = self.weights[correction].weight(syst)[icat],
                                )
                    output['ttbarmass'].fill(
                                        systematic=syst,
                                        anacat = i,
                                        ttbarmass = ttbarmass[icat],
                                        weight = self.weights[correction].weight(syst)[icat],
                                        )
                    
                    '''
                    output['weights'][syst] += np.sum(self.weights[correction].weight(syst))
                    output['systematics'][syst] += len(events.event[icat])
                    output['jetdy'].fill(
                                  systematic=syst,
                                  anacat = i,
                                  jetdy = rapidity[icat],
                                  weight = self.weights[correction].weight(syst)[icat],
                                  )                    
                    
                    output['jetmass'].fill(
                                   systematic=syst,
                                   anacat = i,
                                   jetmass = jetmass[icat],
                                   weight = self.weights[correction].weight(syst)[icat],
                                  )
                    output['jetmsd'].fill(
                                           systematic=syst,
                                           anacat = i,
                                           jetmsd = jetmsd[icat],
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




                    '''
                    
                    
        logger.debug('memory:%s: fill histograms %s:%s', time.time(), correction, get_memory_usage())
        del self.weights[correction]

        return output

    def postprocess(self, accumulator):
        logger.debug('memory:%s: finish processor:%s', time.time(), get_memory_usage())
        return accumulator
        

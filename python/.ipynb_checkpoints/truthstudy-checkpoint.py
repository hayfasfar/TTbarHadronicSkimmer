import awkward as ak
import numpy as np


def _has_flag(status_flags, bit):
    #print("status flag is: ", status_flags)
    #print("bit is : ", bit, "  1 << bit is: ", 1 << bit)
    #print("return statement is: ",  (status_flags & (1 << bit)) != 0)
    return (status_flags & (1 << bit)) != 0


def get_hadronic_tops(genparts):
    # GenPart statusFlags bit 13 = isLastCopy
    last_copy = _has_flag(genparts.statusFlags, 13)
    is_top = (abs(genparts.pdgId) == 6) & last_copy

    
    ## if you want to debug uncomment the following lines: 
    #print("last copy is: ", last_copy)
    #print("pdgID is : ", abs(genparts.pdgId))
    #print("is top? " , is_top) 

    # local indices per event
    idx = ak.local_index(genparts, axis=1) ## indexing all gen particles per event. Instead of for loop in C++
    top_idx = idx[is_top]  ## indexing tops per event.

    
    ## safety statement in case of 0 tops in the event.
    if ak.all(ak.num(top_idx, axis=1) == 0):
        return genparts[is_top]

    # children of each top (shape: [events, nTop, nGenPart])
    child_mask = genparts.genPartIdxMother[:, None, :] == top_idx[:, :, None]
    
    
    #print ("idx is : ", idx)
    #print("top_idx is: ", top_idx)
    #print ("ak.num(top_idx, axis=1) : ", ak.num(top_idx, axis=1))
    #print ("ak.all(ak.num(top_idx, axis=1) == 0) " , ak.all(ak.num(top_idx, axis=1) == 0) )
    #print("child_mask any True (event 0):", ak.any(child_mask[0]))
    #print("child_mask sum (event 0):", ak.sum(child_mask[0]))

    # W child per top: pick index of any W (or -1 if none)
    is_w = child_mask & (abs(genparts.pdgId)[:, None, :] == 24 )
    
    #print("is it W anywhere (event 0):", ak.any(is_w[0]))
    
    w_idx = ak.max(ak.where(is_w, idx[:, None, :], -1), axis=2)

    # W daughters (shape: [events, nTop, nGenPart])
    w_daughter_mask = genparts.genPartIdxMother[:, None, :] == w_idx[:, :, None]


    # W decays has intermediate states,  need to find the last-copy of it before decaying to quarks.
    
    w_last = _has_flag(genparts.statusFlags, 13) & (abs(genparts.pdgId) == 24)
    w_chain_mask = genparts.genPartIdxMother[:, None, :] == w_idx[:, :, None]
    w_last_mask = w_chain_mask & w_last[:, None, :]
    w_last_idx = ak.max(ak.where(w_last_mask, idx[:, None, :], -1), axis=2)

    # W daughters from last-copy W
    w_daughter_mask = genparts.genPartIdxMother[:, None, :] == w_last_idx[:, :, None]
    
    
    is_quark = w_daughter_mask & (abs(genparts.pdgId)[:, None, :] >= 1) & (abs(genparts.pdgId)[:, None, :] <= 5)
    n_quarks = ak.sum(is_quark, axis=2)

    hadronic = n_quarks >= 2
    # hadronic is [events, nTop] bool
    n_had = ak.num(hadronic[hadronic], axis=1)  # number of hadronic tops per event

    n_allhad = ak.sum(n_had == 2)
    n_semilep = ak.sum(n_had == 1)
    n_dilep = ak.sum(n_had == 0)

    print("all-had:", int(n_allhad))
    print("semi-lep:", int(n_semilep))
    print("di-lep:", int(n_dilep))
    tops = genparts[is_top]
    return tops[hadronic]


def truthstudy_counts(genparts, fatjets, subJets, jets, dr_ak8=0.8, dr_ak4=1.2):
    tops = get_hadronic_tops(genparts)
    mask = (ak.num(tops, axis=1) == 2)

    # keep only all‑hadronic events
    tops = tops[mask]
    fatjets = fatjets[mask]
    jets = jets[mask]
    subJets=subJets[mask]
    
    if ak.sum(mask) == 0:
     return {"n_hadtop": 0, "n_hadtop_ak8": 0, "n_hadtop_ak8_ak4": 0}
    
    #all_tops = genparts[abs(genparts.pdgId) == 6]
    #print(ak.num(all_tops, axis=1))

    #build p4 if not present
    print("tops fileds are ", tops.fields)
    
    if "p4" not in tops.fields:
        tops["p4"] = ak.with_name(tops[["pt", "eta", "phi", "mass"]], "PtEtaPhiMLorentzVector")
    if "p4" not in fatjets.fields:
        fatjets["p4"] = ak.with_name(fatjets[["pt", "eta", "phi", "mass"]], "PtEtaPhiMLorentzVector")
    if "p4" not in jets.fields:
        jets["p4"] = ak.with_name(jets[["pt", "eta", "phi", "mass"]], "PtEtaPhiMLorentzVector")
    
    #Match top to AK8 using cartesian pairing (Awkward-1 safe)
    top_fat_pairs = ak.cartesian({"top": tops, "fat": fatjets}, axis=1, nested=True)
    dr_top_fat = top_fat_pairs["top"].p4.delta_r(top_fat_pairs["fat"].p4)  # [evt, nTop, nFat]
    min_dr_fat = ak.fill_none(ak.min(dr_top_fat, axis=2), 999)
    has_ak8 = min_dr_fat < dr_ak8

    # AK4 near top
    top_ak4_pairs = ak.cartesian({"top": tops, "jet": jets}, axis=1, nested=True)
    dr_top_ak4 = top_ak4_pairs["top"].p4.delta_r(top_ak4_pairs["jet"].p4)  # [evt, nTop, nJet]
    ak4_near_top = dr_top_ak4 < dr_ak4

    #dr between all AK8 and AK4
    fat_ak4_pairs = ak.cartesian({"fat": fatjets, "jet": jets}, axis=1, nested=True)
    dr_fat_ak4_all = fat_ak4_pairs["fat"].p4.delta_r(fat_ak4_pairs["jet"].p4)  # [evt, nFat, nJet]

    #pick the AK8 closest to the top, then measure AK4 distance to that AK8
    closest_fat_idx = ak.argmin(dr_top_fat, axis=2)  # [evt, nTop]
    #Awkward-1: use advanced indexing instead of ak.take
    dr_fat_ak4 = dr_fat_ak4_all[closest_fat_idx]
    dr_fat_ak4 = ak.fill_none(dr_fat_ak4, 999)

    ak4_outside_ak8 = dr_fat_ak4 > dr_ak8
    has_ak4_outside = ak.any(ak4_near_top & ak4_outside_ak8, axis=2)
    
    btag_wp = 0.5  # à adapter selon CMS
    is_btag = (jets.btagDeepFlavB > btag_wp)
    has_btag_outside = ak.any(
    ak4_near_top & ak4_outside_ak8 & is_btag[:, None, :], axis=2)  # -> [evt, nTop]
    
    print("number of TRUE AK8  jets:  ", int(ak.sum(ak.ones_like(has_ak8))), 
           "number of MATCHED AK8 jets: ", int(ak.sum(has_ak8)), 
           "number of AK4 close to AK8: ",  int(ak.sum(has_ak8 & has_ak4_outside)),
            "number of AK8 + AK4(btag):", int(ak.sum(has_ak8 & has_btag_outside)))

    evt_has_btag = ak.any(has_btag_outside, axis=1)   # [evt]

    # ---- top tag only for matched AK8 in extra‑AK4 events ----
    matched_fat = fatjets[closest_fat_idx]           # [evt, nTop]
    ttag_wp = 0.4
    is_tagged = (matched_fat.globalParT3_TopbWqq > ttag_wp) & has_ak8
    
    # keep only events with extra AK4 btag
    #is_tagged_btag = is_tagged[evt_has_btag]
    #print("matched+tagged in extra‑AK4 events:", int(ak.sum(is_tagged_btag)))
    matched_tagged_with_btag = is_tagged & has_ak4_outside #has_btag_outside
    print("matched+tagged WITH extra AK4 btag:", int(ak.sum(matched_tagged_with_btag)))
    
    # ---- subjet counts only in extra‑AK4 events ----
    fat_sel = fatjets[evt_has_btag]
    evt_has_2fat = ak.num(fat_sel, axis=1) >= 2
    fat_sel = fat_sel[evt_has_2fat]
    
    order = ak.argsort(fat_sel.pt, ascending=False)
    fat2 = fat_sel[order][:, :2]
    
    nsub_ak8 = ak.values_astype(fat2.subJetIdx1 >= 0, np.int32) + ak.values_astype(fat2.subJetIdx2 >= 0, np.int32)
    
    print("nsubjets in AK8(lead, sublead) per event:", nsub_ak8[:10])
    print("mean nsubjets per AK8:", ak.mean(ak.flatten(nsub_ak8)))


    
    return {
        "n_hadtop": int(ak.sum(ak.ones_like(has_ak8))),
        "n_hadtop_ak8": int(ak.sum(has_ak8)),
        "n_hadtop_ak8_ak4": int(ak.sum(has_ak8 & has_ak4_outside)),
    }

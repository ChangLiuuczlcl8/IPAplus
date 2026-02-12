import os
import pickle
import numpy as np
import pandas as pd
from ipaPy2 import ipa
from ipaplus import util


# load data
db = pd.read_csv('Dataset/Database/IPA_MS1.csv')
DBMS2 = pd.read_csv('Dataset/Database/IPA_MS2.csv')
adducts = pd.read_csv('Dataset/Database/adducts.csv')
dfMS1 = pd.read_csv("Dataset/MTBLS2207/ms1.csv")
dfMS2 = pd.read_csv(os.path.join("Dataset/MTBLS2207/ms2.csv"))
ionisation = 1

# preprocessed results from single-peak annotation models
with open("Dataset/MTBLS2207/annotations_sirius.pkl","rb") as fp:
    annotations_sirius = pickle.load(fp)
with open("Dataset/MTBLS2207/annotations_metfrag.pkl", "rb") as fp:
    annotations_metfrag = pickle.load(fp)
with open("Dataset/MTBLS2207/annotations_msfinder.pkl", "rb") as fp:
    annotations_msfinder = pickle.load(fp)

# add new compounds into IPAplus database
# add sirius prediction into IPAplus database
keys_sirius = list(annotations_sirius.keys())
keys_metfrag = list(annotations_metfrag.keys())
keys_msfinder = list(annotations_msfinder.keys())
r_sirius = 0.5
r_metfrag = 0.7
r_msfinder = 0.5
for f1 in keys_sirius:
    sirius_results = annotations_sirius[f1]
    sirius_results = sirius_results.iloc[:min(10, len(sirius_results.index)), :]
    for i in range(len(sirius_results.index)):
        if sirius_results.iloc[i, -1].startswith("sirius"):
            idx = (db["inchi"].str.contains(sirius_results.iloc[i, 3], na=False, regex=False) & db["id"].str.startswith("peak_" + str(f1)))
            if sum(idx) == 0:
                db.loc[len(db.index)] = ["peak_" + str(f1) + "_" + sirius_results.iloc[i, -1]] + [sirius_results.iloc[i, 0]] + [
                                        sirius_results.iloc[i, 1]] + [sirius_results.iloc[i, 3]] + [sirius_results.iloc[i, 5]] + [
                                        np.nan] + [sirius_results.iloc[i, 2]] + [db.iloc[0, 7]] + [
                                        np.nan] + [0.5] + [np.nan] + [np.nan]
            else:
                db.loc[idx, "id"] = db.loc[idx, "id"] + "_" + sirius_results.iloc[i, -1]
                adduct_types = db.loc[idx, "adductsPos"].str.split(";").iloc[0]
                if sirius_results.iloc[i, 2] not in adduct_types:
                    db.loc[db["id"] == sirius_results.iloc[i, -1], "adductsPos"] = db.loc[db["id"] == sirius_results.iloc[i, -1], "adductsPos"] + ";" + sirius_results.iloc[i, 2]
        else:
            adduct_types = db.loc[db["id"] == sirius_results.iloc[i, -1], "adductsPos"].str.split(";").iloc[0]
            if sirius_results.iloc[i, 0] not in adduct_types:
                db.loc[db["id"] == sirius_results.iloc[i, -1], "adductsPos"] = db.loc[db["id"] == sirius_results.iloc[i, -1], "adductsPos"] + ";" + sirius_results.iloc[i, 2]
# add metfrag prediction into IPAplus database
for f1 in keys_metfrag:
    metfrag_results = annotations_metfrag[f1]
    metfrag_results = metfrag_results.iloc[:min(10, len(metfrag_results.index)), :]
    for i in range(len(metfrag_results.index)):
        if metfrag_results.iloc[i, -1].startswith("metfrag"):
            idx = (db["inchi"].str.contains(metfrag_results.iloc[i, 2], na=False, regex=False) & db["id"].str.startswith("peak_" + str(f1)))
            if sum(idx) == 0:
                db.loc[len(db.index)] = ["peak_" + str(f1) + "_" + metfrag_results.iloc[i, -1]] + [metfrag_results.iloc[i, 7]] + [
                    metfrag_results.iloc[i, 4]] + [metfrag_results.iloc[i, 2]] + [metfrag_results.iloc[i, 5]] + [
                                          np.nan] + ["M+H"] + [db.iloc[0, 7]] + [np.nan] + [0.5] + [
                                          np.nan] + [np.nan]
            else:
                db.loc[idx, "id"] = db.loc[idx, "id"] + "_" + metfrag_results.iloc[i, -1]
                adduct_types = db.loc[idx, "adductsPos"].str.split(";").iloc[0]
                if "M+H" not in adduct_types:
                    db.loc[idx, "adductsPos"] = db.loc[idx, "adductsPos"] + ";M+H"
# add msfinder prediction into IPAplus database
for f1 in keys_msfinder:
    msfinder_results = annotations_msfinder[f1]
    msfinder_results = msfinder_results.iloc[:min(10, len(msfinder_results.index)), :]
    for i in range(len(msfinder_results.index)):
        if msfinder_results.iloc[i, -1].startswith("msfinder"):
            idx = (db["inchi"].str.contains(msfinder_results.iloc[i, 6], na=False, regex=False) & db["id"].str.startswith("peak_" + str(f1)))
            if sum(idx) == 0:
                db.loc[len(db.index)] = ["peak_" + str(f1) + "_" + msfinder_results.iloc[i, -1]] + [msfinder_results.iloc[i, 1]] + [
                                        msfinder_results.iloc[i, 3]] + [msfinder_results.iloc[i, 6]] + [msfinder_results.iloc[i, 5]] + [
                                        np.nan] + [msfinder_results.iloc[i, 0]] + [db.iloc[0, 7]] + [np.nan] + [0.5] + [
                                        np.nan] + [np.nan]
            else:
                db.loc[idx, "id"] = db.loc[idx, "id"] + "_" + msfinder_results.iloc[i, -1]
                adduct_types = db.loc[idx, "adductsPos"].str.split(";").iloc[0]
                if msfinder_results.iloc[i, 0] not in adduct_types:
                    db.loc[idx, "adductsPos"] = db.loc[idx, "adductsPos"] + ";" + msfinder_results.iloc[i, 0]
        else:
            adduct_types = db.loc[db["id"] == msfinder_results.iloc[i, -1], "adductsPos"].str.split(";").iloc[0]
            if msfinder_results.iloc[i, 0] not in adduct_types:
                db.loc[db["id"] == msfinder_results.iloc[i, -1], "adductsPos"] = db.loc[db["id"] == msfinder_results.iloc[i, -1], "adductsPos"] + ";" + msfinder_results.iloc[i, 0]
db.to_csv("Dataset/Database/IPA_MS1_updated.csv", index=False)

# run IPAplus integration
annotations_ipaplus = ipa.simpleIPA(df=dfMS1, dfMS2=dfMS2, ionisation=ionisation, DB=db, DBMS2=DBMS2, noits=0, adductsAll=adducts, ppm=5, ncores=2, isodiff=1.003355)
for f1 in annotations_ipaplus:
    annotation_ipaplus_ = annotations_ipaplus[f1]
    annotation_ipaplus_ = annotation_ipaplus_[~(annotation_ipaplus_.id.str.startswith("peak") & (~annotation_ipaplus_.id.str.startswith("peak_" + str(f1))))]
    annotations_ipaplus[f1] = annotation_ipaplus_
keys_ipaplus = list(annotations_ipaplus.keys())
for f1 in keys_ipaplus:
    annotation_ipaplus_ = annotations_ipaplus[f1]
    annotation_ipaplus_ = annotation_ipaplus_.assign(inchi=None)
    for m in range(len(annotation_ipaplus_.index)):
        if annotation_ipaplus_.iloc[m, 0] != "Unknown":
            annotation_ipaplus_.iloc[m, -1] = db["inchi"][db["id"] == annotation_ipaplus_.iloc[m, 0]].values[0]
    prior_prob = annotation_ipaplus_["post"].values
    annotation_ipaplus_ = annotation_ipaplus_.assign(post_original=prior_prob)
    n = len(prior_prob)
    # integrate sirius results
    if f1 in keys_sirius:
        sirius_results = annotations_sirius[f1]
        sirius_results = sirius_results.iloc[:10]
        cond_sirius = np.zeros(n)
        for cp in range(len(sirius_results.index)):
            if sirius_results.iloc[cp, -1].startswith("sirius"):
                new_candidates_idx = annotation_ipaplus_.id.apply(lambda x: util.split_every_n_special(x, "_", 2)).apply(
                    lambda x: "peak_" + str(f1) in x and "sirius_" + str(cp + 1) in x)
                cond_sirius[new_candidates_idx] += sirius_results.iloc[cp, 7]
            else:
                cond_sirius[annotation_ipaplus_['id'] == sirius_results.iloc[cp, -1]] += sirius_results.iloc[cp, 7]
    # integrate metfrag results
    if f1 in keys_metfrag:
        metfrag_results = annotations_metfrag[f1]
        metfrag_results = metfrag_results.iloc[:10]
        cond_metfrag = np.zeros(n)
        for cp in range(len(metfrag_results.index)):
            if metfrag_results.iloc[cp, -1].startswith("metfrag"):
                new_candidates_idx = annotation_ipaplus_.id.apply(lambda x: util.split_every_n_special(x, "_", 2)).apply(
                    lambda x: "peak_" + str(f1) in x and metfrag_results.iloc[cp, -1] in x)
                cond_metfrag[new_candidates_idx] += metfrag_results.iloc[cp, 1]
            else:
                cond_metfrag[annotation_ipaplus_['id'] == metfrag_results.iloc[cp, -1]] += metfrag_results.iloc[cp, 1]
    # integrate msfinder results
    if f1 in keys_msfinder:
        msfinder_results = annotations_msfinder[f1]
        msfinder_results = msfinder_results.iloc[:10]
        cond_msfinder = np.zeros(n)
        for cp in range(len(msfinder_results.index)):
            if msfinder_results.iloc[cp, -1].startswith("msfinder"):
                new_candidates_idx = annotation_ipaplus_.id.apply(lambda x: util.split_every_n_special(x, "_", 2)).apply(
                    lambda x: "peak_" + str(f1) in x and msfinder_results.iloc[cp, -1] in x)
                cond_msfinder[new_candidates_idx] += msfinder_results.iloc[cp, 2]
            else:
                cond_msfinder[annotation_ipaplus_['id'] == msfinder_results.iloc[cp, -1]] += msfinder_results.iloc[cp, 2]
    if f1 in keys_sirius:
        prior_prob = util.cal_posterior_3(prior_prob, cond_sirius, r_sirius)
    if f1 in keys_metfrag:
        prior_prob = util.cal_posterior_3(prior_prob, cond_metfrag, r_metfrag)
    if f1 in keys_msfinder:
        prior_prob = util.cal_posterior_3(prior_prob, cond_msfinder, r_msfinder)
    annotation_ipaplus_["post"] = prior_prob
    annotations_ipaplus[f1] = annotation_ipaplus_
ipa.Gibbs_sampler_add(dfMS1, annotations_ipaplus, noits=30000, delta_add=0.1)
with open("Dataset/MTBLS2207/output/annotations_ipaplus.pkl", "wb") as fp:
    pickle.dump(annotations_ipaplus, fp)


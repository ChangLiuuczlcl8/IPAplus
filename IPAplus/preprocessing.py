import os
import pymzml
import pickle
import requests
import numpy as np
import pandas as pd


# load data
db = pd.read_csv('/Users/user/Documents/dataset/ipapy2/DB/IPA_MS1.csv')
dfMS1 = pd.read_csv("/Users/user/Documents/dataset/MTBLS2207/ipa_input/ms1/M3T-Std_Yeast_pos_DDA_3mz.csv")

# mzml data preprocessing
df = pd.read_csv("/Users/user/Documents/dataset/MTBLS2207/mzmatch/M3T-Std_Yeast_pos_DDA_3mz.csv")
group_label = 0
group_labels = [0]
for i in range(1, df.shape[0]):
    if abs(df.iloc[i, 3] - df.iloc[i - 1, 3]) > 0.002:
        group_label += 1
    group_labels.append(group_label)
df["cluster"] = group_labels
df = df.sort_values(["cluster", "RT", "AVGMZ"], ascending=[True, True, False]).reset_index(drop=True)
df["ids"] = range(1, df.shape[0]+1)
df = df[["ids", "cluster", "AVGMZ", "RT", "MAXINTENSITY"]] #use max intensity or average intensity
df = df.rename(columns={"cluster": "rel.ids", "AVGMZ": "mzs", "RT": "RTs", "MAXINTENSITY": "Int"})
df.to_csv("/Users/user/Documents/dataset/MTBLS2207/ipa_input/ms1/M3T-Std_Yeast_pos_DDA_3mz.csv", index=False)
mzml_file = "/Users/user/Documents/dataset/MTBLS2207/mzml/M3T-Std_Yeast_pos_DDA_3mz.mzML"
run = pymzml.run.Reader(mzml_file)
ms2_spectra = pd.DataFrame(columns=["id", "spectrum", "eV"])
for n, spec in enumerate(run):
    if spec.ms_level == 2:
        ms2 = []
        precursur_mz = spec.selected_precursors[0]["mz"]
        precursur_rt = spec.scan_time[0] * 60
        spec_id = df["ids"][df["mzs"].between(precursur_mz -0.01, precursur_mz + 0.01) & df["RTs"].between(precursur_rt - 30, precursur_rt + 30)].values
        if len(spec_id) > 0:
            if len(spec_id) > 1:
                precursor_candidates = df.iloc[spec_id - 1].copy()
                precursor_candidates["RTs"] = abs(precursor_candidates["RTs"] - precursur_rt)
                precursor_candidates = precursor_candidates.sort_values(by=['RTs'])
                ms2.append(precursor_candidates.iloc[0, 0])
            else:
                ms2.append(spec_id[0])
            ms2_spectrum = spec.centroidedPeaks
            ms2_spectrum = pd.DataFrame(ms2_spectrum, columns=["mzs", "int"])
            ms2_spectrum["merged"] = ms2_spectrum["mzs"].round(4).astype(str) + ":" + ms2_spectrum["int"].round(4).astype(str)
            ms2.append(ms2_spectrum["merged"].str.cat(sep=' '))
            ms2.append(int(spec['collision energy']))
            ms2_spectra.loc[len(ms2_spectra)] = ms2
ms2_spectra = ms2_spectra.sort_values(by=['id'])
ms2_spectra.to_csv("/Users/user/Documents/dataset/MTBLS2207/ipa_input/ms2/M3T-Std_Yeast_pos_DDA_3mz.csv", index=False)

# add inchi to all candidates in IPA database
for i in range(len(db.index)):
    if pd.isnull(db.iloc[i, 3]):
        name__ = db.iloc[i, 0]
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/InChI/TXT" % name__
        response = requests.get(url)
        inchi = response.text.strip().split('\n')[0]
        if inchi.startswith("Status:"):
            name__ = db.iloc[i, 1]
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/InChI/TXT" % name__
            response = requests.get(url)
            inchi = response.text.strip().split('\n')[0]
            if inchi.startswith("Status:"):
                inchi = np.nan
        db.iloc[i, 3] = inchi
db.to_csv('/Users/user/Documents/dataset/ipapy2/DB/IPA_MS1.csv', index=False)

# SIRIUS results
sirius_results = pd.read_csv("/Users/user/Documents/dataset/MTBLS2207/sirius/compound_identifications_all.csv")
annotations_sirius = {}
for feature in set(sirius_results["mappingFeatureId"]):
    candidate_structures = sirius_results.loc[sirius_results["mappingFeatureId"] == feature]
    candidate_structures = candidate_structures.iloc[:min(len(candidate_structures.index), 20)]
    candidate_mz = candidate_structures.iloc[0, 19]
    candidate_rt = candidate_structures.iloc[0, 20]
    candidate_id = dfMS1["ids"][dfMS1["mzs"].between(candidate_mz - 0.01, candidate_mz + 0.01) & dfMS1["RTs"].between(candidate_rt - 0.6, candidate_rt + 0.6)].values
    if len(candidate_id) == 0:
        continue
    annotation_sirius = candidate_structures[["name", "molecularFormula", "adduct", "InChI", "InChIkey2D", "smiles", "CSI:FingerIDScore"]]
    annotation_sirius = annotation_sirius.assign(conditional=-1 / annotation_sirius["CSI:FingerIDScore"].copy())
    annotation_sirius.loc[:, "conditional"] = annotation_sirius["conditional"].copy() / sum(annotation_sirius["conditional"])
    annotation_sirius = annotation_sirius.assign(id=["sirius_" + str(i + 1) for i in range(len(annotation_sirius.index))])
    annotation_sirius.InChI = annotation_sirius.InChI.str.split("/").str[:4].str.join("/")
    for i in range(len(annotation_sirius.index) - 1, 0, -1):
        inchi = annotation_sirius.iloc[i, 3]
        if inchi in annotation_sirius.InChI.values[:i]:
            annotation_sirius = annotation_sirius.drop([i])
    annotation_sirius = annotation_sirius.reset_index(drop=True)
    for i in range(len(annotation_sirius.index)):
        if sum(db["inchi"].str.contains(annotation_sirius.iloc[i, 3], na=False, regex=False)) != 0:
            annotation_sirius.iloc[i, -1] = db["id"][db["inchi"].str.contains(annotation_sirius.iloc[i, 3], na=False, regex=False)].iloc[0]
    annotation_sirius.adduct = annotation_sirius.adduct.str.split(pat='\[|\]').str[1].str.replace(" ", "")
    annotations_sirius[candidate_id[0]] = annotation_sirius
annotations_sirius = dict(sorted(annotations_sirius.items()))
with open("/Users/user/Documents/dataset/MTBLS2207/results/sirius/annotations_sirius.pkl", "wb") as fp:
    pickle.dump(annotations_sirius, fp)

# MetFrag data preprocessing
os.mkdir("/Users/user/Documents/dataset/MTBLS2207/metfrag_input/M3T-Std_Yeast_pos_DDA_3mz")
os.mkdir("/Users/user/Documents/dataset/MTBLS2207/metfrag/M3T-Std_Yeast_pos_DDA_3mz")
ms1_spec = pd.read_csv("/Users/user/Documents/dataset/MTBLS2207/ipa_input/ms1/M3T-Std_Yeast_pos_DDA_3mz.csv")
mzml_file = '/Users/user/Documents/dataset/MTBLS2207/mzml/M3T-Std_Yeast_pos_DDA_3mz.mzML'
run = pymzml.run.Reader(mzml_file)
ms2_spectra = pd.DataFrame(columns=["id", "spectrum", "precursor", "RT", "max_int"])
for n, spec in enumerate(run):
    if spec.ms_level == 2:
        precursur_mz = spec.selected_precursors[0]["mz"]
        precursur_rt = spec.scan_time[0] * 60
        spec_id = ms1_spec[["ids", "mzs", "RTs"]][
            ms1_spec["mzs"].between(precursur_mz - 0.01, precursur_mz + 0.01) & ms1_spec["RTs"].between(
                precursur_rt - 30, precursur_rt + 30)].copy()
        if len(spec_id.index) > 0:
            ms2_spectrum = spec.centroidedPeaks
            ms2_spectrum = pd.DataFrame(ms2_spectrum, columns=["mzs", "int"])
            if len(ms2_spectrum.index) == 0:
                continue
            if len(spec_id.index) > 1:
                spec_id["RTs"] = abs(spec_id["RTs"] - precursur_rt)
                spec_id = spec_id.sort_values(by=['RTs'])
                ms2_spectra.loc[len(ms2_spectra), "id"] = spec_id.iloc[0, 0]
                ms2_spectra.iloc[len(ms2_spectra) - 1, 2] = spec_id.iloc[0, 1]
                ms2_spectra.iloc[len(ms2_spectra) - 1, 3] = precursur_rt
            else:
                ms2_spectra.loc[len(ms2_spectra), "id"] = spec_id.iloc[0, 0]
                ms2_spectra.iloc[len(ms2_spectra) - 1, 2] = spec_id.iloc[0, 1]
                ms2_spectra.iloc[len(ms2_spectra) - 1, 3] = precursur_rt
            ms2_spectrum["merged"] = ms2_spectrum["mzs"].round(5).astype(str) + " " + ms2_spectrum[
                "int"].round(5).astype(str)
            ms2_spectrum_ = ms2_spectrum["merged"].tolist()
            ms2_spectra.iloc[len(ms2_spectra) - 1, 1] = ms2_spectrum_
            ms2_spectra.iloc[len(ms2_spectra) - 1, 4] = max(ms2_spectrum["int"])
for i in set(ms2_spectra.id):
    ms2_spectra_ = ms2_spectra.loc[ms2_spectra.id == i]
    max_int = max(ms2_spectra_["max_int"])
    for j in ms2_spectra_.index:
        if ms2_spectra_.loc[j, "max_int"] < max_int / 100:
            ms2_spectra = ms2_spectra.drop(index=j)
ms2_spectra = ms2_spectra.sort_values(by=['id', 'RT'])
ms2_spectra = ms2_spectra.reset_index(drop=True)
for d in set(ms2_spectra.iloc[:, 0]):
    os.mkdir("/Users/user/Documents/dataset/MTBLS2207/metfrag_input/pos/" + f1[:-4] + "/" + str(d))
    os.mkdir("/Users/user/Documents/dataset/MTBLS2207/metfrag/" + f1[:-4] + "/" + str(d))
    ms2_spectra_ = ms2_spectra.loc[ms2_spectra['id'] == d]
    for k in range(len(ms2_spectra_.index)):
        pklist = ms2_spectra_.iloc[k, 1]
        pl_d = os.path.join("/Users/user/Documents/dataset/MTBLS2207/metfrag_input/pos", f1[:-4], str(d))
        with open(pl_d + '/peaklist_' + str(k + 1) + '.txt', 'w') as f2:
            for line in pklist:
                f2.write(f"{line}\n")
        with open(pl_d + '/parameters_' + str(k + 1) + '.cfg', 'w') as f3:
            parameters_fix = ["MetFragDatabaseType = PubChem\n",
                              "PrecursorIonMode = 1\n",
                              "IsPositiveIonMode = true\n",
                              "DatabaseSearchRelativeMassDeviation = 5.0\n",
                              "FragmentPeakMatchRelativeMassDeviation = 5.0\n",
                              "FragmentPeakMatchAbsoluteMassDeviation = 0.001\n",
                              "MetFragScoreTypes = FragmenterScore\n",
                              "MetFragScoreWeights = 1.0\n",
                              "MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter,IsotopeFilter\n",
                              "MaximumTreeDepth = 2\n",
                              "NumberThreads = 2\n",
                              "UseSmiles = true\n",
                              "MetFragCandidateWriter = CSV\n"
                              "IonizedPrecursorMass = " + str(ms2_spectra_.iloc[k, 2]) + "\n",
                              "PeakListPath = " + pl_d + "/peaklist_" + str(k + 1) + ".txt\n",
                              "ResultsPath = /Users/user/Documents/dataset/MTBLS2207/metfrag/" + f1[:-4] + "/" + str(ms2_spectra_.iloc[k, 0]) + "\n",
                              "SampleName = MetFragCL_Candidates_" + str(k + 1) + "\n"]
            f3.writelines(parameters_fix)

# MetFrag results
annotations_metfrag = {}
d_names = os.listdir("/Users/user/Documents/dataset/MTBLS2207/metfrag/M3T-Std_Yeast_pos_DDA_3mz")
for d1 in d_names:
    if not os.path.isdir("/Users/user/Documents/dataset/MTBLS2207/metfrag/M3T-Std_Yeast_pos_DDA_3mz/" + d1):
        continue
    root_ = "/Users/user/Documents/dataset/MTBLS2207/metfrag/M3T-Std_Yeast_pos_DDA_3mz/" + d1
    f_names_ = os.listdir(root_)
    metfrag_results = pd.DataFrame(columns=["FragmenterScore", "Score", "InChI", "Identifier", "MolecularFormula", "SMILES", "InChIKey", "IUPACName"])
    for f1 in f_names_:
        if f1.endswith(".csv"):
            metfrag_result = pd.read_csv(root_ + "/" + f1)
            if len(metfrag_result.index) == 0:
                continue
            metfrag_result = metfrag_result.iloc[:min(len(metfrag_result.index), 100), [13, 0, 6, 8, 14, 2, 3, 12]]
            metfrag_results = pd.merge(metfrag_result, metfrag_results, how="outer", on=["Identifier", "InChI", "MolecularFormula", "SMILES", "InChIKey", "IUPACName"])
            metfrag_results.iloc[:, 0] = metfrag_results.iloc[:, 0].astype(float).fillna(0) + metfrag_results.iloc[:, 8].astype(float).fillna(0)
            metfrag_results.iloc[:, 1] = metfrag_results.iloc[:, 1].astype(float).fillna(0) + metfrag_results.iloc[:, 9].astype(float).fillna(0)
            metfrag_results = metfrag_results.iloc[:, :-2]
    metfrag_results.iloc[:, 0] = metfrag_results.iloc[:, 0]/len(f_names_)
    metfrag_results.iloc[:, 1] = metfrag_results.iloc[:, 1]/len(f_names_)
    metfrag_results = metfrag_results.sort_values([metfrag_results.columns[0], metfrag_results.columns[1]], ascending=[False, False]).reset_index(drop=True)
    metfrag_results = metfrag_results.assign(id=["metfrag_" + str(i + 1) for i in range(len(metfrag_results.index))])
    metfrag_results.to_csv(root_ + "/MetFragCL_Candidates.csv", index=False)
    if len(metfrag_results.index) == 0:
        continue
    metfrag_results = metfrag_results.iloc[:min(len(metfrag_results.index), 20)]
    metfrag_results.InChI = metfrag_results.InChI.str.split("/").str[:4].str.join("/")
    for i in range(len(metfrag_results.index) - 1, 0, -1):
        inchi = metfrag_results.iloc[i, 2]
        if inchi in metfrag_results.InChI.values[:i]:
            metfrag_results = metfrag_results.drop([i])
    metfrag_results = metfrag_results.reset_index(drop=True)
    for i in range(len(metfrag_results.index)):
        if sum(db["inchi"].str.contains(metfrag_results.iloc[i, 2], na=False, regex=False)) != 0:
            metfrag_results.iloc[i, -1] = db["id"][db["inchi"].str.contains(metfrag_results.iloc[i, 2], na=False, regex=False)].iloc[0]
    annotations_metfrag[int(d1)] = metfrag_results
annotations_metfrag = dict(sorted(annotations_metfrag.items()))
with open("/Users/user/Documents/dataset/MTBLS2207/results/metfrag/annotations_metfrag.pkl", "wb") as fp:
    pickle.dump(annotations_metfrag, fp)

# ms-finder data preprocessing
ms1_spec = pd.read_csv("/Users/user/Documents/dataset/MTBLS2207/ipa_input/ms1/M3T-Std_Yeast_pos_DDA_3mz.csv")
for root, d_names, f_names in os.walk("/Users/user/Documents/dataset/MTBLS2207/metfrag_input/pos/M3T-Std_Yeast_pos_DDA_3mz"):
    for d1 in d_names:
        ms1list = ms1_spec.loc[ms1_spec["rel.ids"] == ms1_spec.iloc[int(d1)-1, 1]]
        ms1list["mz_int"] = ms1list["mzs"].astype(str) + " " + ms1list["Int"].astype(str)
        ms1list_ = '\n'.join(ms1list["mz_int"].values)
        files = os.listdir("/Users/user/Documents/dataset/MTBLS2207/metfrag_input/M3T-Std_Yeast_pos_DDA_3mz/" + d1)
        for f1 in files:
            if f1.endswith(".txt"):
                with open("/Users/user/Documents/dataset/MTBLS2207/metfrag_input/M3T-Std_Yeast_pos_DDA_3mz/" + d1 + "/" + f1) as pl:
                    ms2list = pl.read()
                k = f1[9:-4]
                parameters_fix = ["NAME: yeast_pos_3_peak_" + d1 + "_" + k + "\n",
                                  "PRECURSORMZ: " + str(ms1_spec.iloc[int(d1)-1, 2]) + "\n",
                                  "PRECURSORTYPE: [M+H]+ \n",
                                  "IONMODE: Positive \n",
                                  "MSTYPE: MS1 \n",
                                  "Num Peaks: " + str(len(ms1list.index)) + "\n",
                                  ms1list_ + "\n",
                                  "MSTYPE: MS2 \n",
                                  "Num Peaks: " + str(len(ms2list.splitlines())) + "\n",
                                  ms2list]
                with open("/Users/user/Documents/dataset/MTBLS2207/msfinder_input/pos/M3T-Std_Yeast_pos_DDA_3mz/yeast_pos_3_peak_" + d1 + "_" + k + ".mat", 'w') as f3:
                    f3.writelines(parameters_fix)

# ms-finder results
annotations_msfinder = {}
msfinder = pd.read_csv('/Users/user/Documents/dataset/MTBLS2207/msfinder/M3T-Std_Yeast_pos_DDA_3mz.txt', sep="\t")
msfinder = msfinder.iloc[:, [2, 6, 7, 8, 10, 12, 13]]
msfinder = msfinder.loc[msfinder["Formula"] != "Spectral DB search"]
msfinder = msfinder.loc[msfinder["Total score"] > 0]
msfinder = msfinder.reset_index(drop=True)
msfinder["peak"] = msfinder["Title"].str.split("_", expand=True).iloc[:, 5]
peaks = list(set(msfinder["peak"]))
peaks = sorted(peaks)
for peak in peaks:
    msfinder_peak = msfinder.loc[msfinder["peak"] == peak]
    msfinder_peak = msfinder_peak.iloc[:, :-1]
    msfinder_specs = pd.DataFrame(columns=["Total score", "Structure", "Formula", "SMILES", "InChIKey", "Precursor type"])
    specs = list(set(msfinder_peak["Title"]))
    for spec in specs:
        msfinder_spec = msfinder_peak.loc[msfinder_peak["Title"] == spec]
        msfinder_spec = msfinder_spec.iloc[:, 1:]
        msfinder_specs = pd.merge(msfinder_spec, msfinder_specs, how="outer", on=["InChIKey", "SMILES", "Structure", "Formula", "Precursor type"])
        msfinder_specs.iloc[:, -1] = msfinder_specs.iloc[:, -1].astype(float).fillna(0)
        msfinder_specs.iloc[:, 2] = msfinder_specs.iloc[:, 2].astype(float).fillna(0)
        msfinder_specs.iloc[:, 2] = msfinder_specs.iloc[:, 2] + msfinder_specs.iloc[:, -1].astype(float)
        msfinder_specs = msfinder_specs.iloc[:, :-1]
    msfinder_specs = msfinder_specs.sort_values(by=[msfinder_specs.columns[2]], ascending=[False]).reset_index(drop=True)
    msfinder_specs = msfinder_specs.iloc[:min(10, len(msfinder_specs.index)), :]
    msfinder_specs.iloc[:, 2] = msfinder_specs.iloc[:, 2]/len(specs)
    msfinder_specs["InChI"] = None
    for i in range(len(msfinder_specs.index)):
        inchikey = msfinder_specs.iloc[i, 4].split("-")[0]
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/property/InChI/TXT" % inchikey
        response = requests.get(url)
        inchi = response.text.strip().split('\n')[0]
        if inchi.startswith("Status:"):
            smiles = msfinder_specs.iloc[i, 5]
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/property/InChI/TXT" % smiles
            response = requests.get(url)
            inchi = response.text.strip().split('\n')[0]
            if inchi.startswith("Status:") or inchi == '':
                name = msfinder_specs.iloc[0, 1]
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/InChI/TXT" % name
                response = requests.get(url)
                inchi = response.text.strip().split('\n')[0]
                if inchi.startswith("Status:"):
                    inchi = np.nan
        msfinder_specs.iloc[i, 6] = inchi
    msfinder_specs.InChI = msfinder_specs.InChI.str.split("/").str[:4].str.join("/")
    for i in range(len(msfinder_specs.index) - 1, 0, -1):
        inchi = msfinder_specs.iloc[i, 6]
        if inchi in msfinder_specs.InChI.values[:i]:
            msfinder_specs = msfinder_specs.drop([i])
    msfinder_specs = msfinder_specs.reset_index(drop=True)
    msfinder_specs = msfinder_specs.assign(id=["msfinder_" + str(i + 1) for i in range(len(msfinder_specs.index))])
    for i in range(len(msfinder_specs.index)):
        if len(db["id"].loc[db["inchi"] == msfinder_specs.iloc[i, 6]]) != 0:
            msfinder_specs.iloc[i, -1] = db["id"].loc[db["inchi"] == msfinder_specs.iloc[i, 6]].iloc[0]
    msfinder_specs = msfinder_specs.rename(columns={"Precursor type": "adduct"})
    msfinder_specs.adduct = msfinder_specs.adduct.str.split(pat='\[|\]').str[1]
    annotations_msfinder[int(peak)] = msfinder_specs
annotations_msfinder = dict(sorted(annotations_msfinder.items()))
for peak in list(annotations_msfinder.keys()):
    annotation_msfinder = annotations_msfinder[peak]
    annotation_msfinder = annotation_msfinder.loc[~annotation_msfinder.InChI.isnull()]
    annotation_msfinder = annotation_msfinder.reset_index(drop=True)
    annotations_msfinder[peak] = annotation_msfinder
    if len(annotation_msfinder.index) == 0:
        del annotations_msfinder[peak]
with open("/Users/user/Documents/dataset/MTBLS2207/results/msfinder/annotations_msfinder.pkl", "wb") as fp:
    pickle.dump(annotations_msfinder, fp)


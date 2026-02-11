# IPAplus

### Motivation: 
Metabolomics is established in biomedical research, aiming to characterise the composition and functions of metabolites within biological systems. Mass spectrometry methods play a crucial role. Methods for annotating a single spectral peak have seen significant development, however spectral annotation of mixtures remains challenging. Probabilistic frameworks like Integrated Probabilistic Annotation ([IPA](https://github.com/francescodc87/ipaPy2)) leverage database scoring and contextual information to improve annotations, but are hindered by limited native databases and assumptions. Single-peak annotation tools can provide complementary, additional method-specific identifications. 
### Results: 
IPA is a Bayesian tool for estimating the posterior probabilities of the presence of metabolites in a particular sample given the MS data. Here, we extend the previous generic IPA pipeline using a systematic and rigorous integration of multiple single-peak annotation methods. We optimised different likelihood functions to achieve a Bayesian explainable integration from different single-peak annotation methods, testing three state-of-the-art single-peak identification tools, SIRIUS, MetFrag and MS-FINDER. The novel integrated IPAplus pipeline adds new compounds to the database and leads to more confident and plausible annotations, as well as a general framework for incorporating further sources of evidence and unifying the diverse toolset of metabolome annotation. 

![plot](./IPAplus/plot.png)


### Example Datasets:
[MTBLS2207] (https://www.ebi.ac.uk/metabolights/editor/MTBLS2207/descriptors) 
[MSV000089047] (https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=e3ec6a27ae3f419aa3c2fdc7a72af71b#%7B%22table_sort_history%22%3A%22main.collection_asc%22%7D)

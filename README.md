# 3D-QSAR Model for Prediction of hERG Channel Blockade

This repository contains a three-dimensional quantitative structure–activity relationship (3D-QSAR) model developed to predict hERG (KCNH2) potassium channel blockade, a major cardiotoxicity liability in small-molecule drug discovery.

The model is based on a Random Forest classifier trained on molecular descriptors computed from RDKit-generated three-dimensional conformers.

---

## Dataset

Bioactivity data were obtained from the ChEMBL database for compounds annotated against the hERG (KCNH2) target.

- **ChEMBL Target ID:** CHEMBL240  
- **Processed dataset:** `hERG_activity.csv`  
- **Activity types retained:** IC₅₀, Kᵢ  
- **Units:** nM  

Compounds were labeled using a binary classification scheme based on pChEMBL values:

- **Blocker (1):** pChEMBL ≥ 5.5  
- **Non-blocker (0):** pChEMBL < 5.5  

---

## Molecular Preparation

SMILES strings were converted to RDKit molecule objects and explicit hydrogens were added prior to three-dimensional embedding.

Three-dimensional conformers were generated using the ETKDGv3 algorithm, followed by geometry optimization with the UFF force field. Molecules for which conformer generation or optimization failed were excluded from further analysis.

---

## Descriptor Calculation

A set of physicochemical and three-dimensional shape descriptors was computed for each optimized conformer using RDKit.

The descriptor set includes:

- Molecular weight  
- LogP  
- Topological polar surface area (TPSA)  
- Hydrogen bond donors and acceptors  
- Number of rotatable bonds  
- Radius of gyration  
- Asphericity  
- Spherocity index  
- Principal moments of inertia (PMI1, PMI2, PMI3)  

The resulting descriptor matrix is provided as `CHEMBL240_3D_descriptors.csv`.

---

## Model Development

- **Algorithm:** Random Forest classifier  
- **Number of trees:** 800  
- **Maximum tree depth:** 14  
- **Minimum samples per leaf:** 8  
- **Train–test split:** 80:20 (stratified)  
- **Random seed:** 42  

The model was trained using the descriptor matrix as input and the binary hERG activity labels as output.

---

## Model Performance

Predictive performance was evaluated on a held-out test set using the area under the receiver operating characteristic curve (ROC–AUC).

- **ROC–AUC:** 0.8139  

---

## Inference Procedure

For prediction on new compounds, a stochastic conformer-based inference strategy is employed. For a given SMILES string:

1. Fifteen independent 3D conformers are generated using ETKDGv3  
2. Each conformer is optimized using UFF  
3. Descriptors are computed for each conformer  
4. The Random Forest model predicts a hERG blocker probability for each conformer  

The final reported prediction is the arithmetic mean of the conformer-level probabilities. Minimum, maximum, and standard deviation values are also reported to reflect conformational variability.

---

## Repository Contents

```
.
├── model.py                         # Data processing, training, and inference code
├── hERG_activity.csv                # Curated bioactivity dataset
├── CHEMBL240_3D_descriptors.csv     # 3D descriptor matrix
├── CHEMBL240_3D_molecules.pkl       # RDKit molecule objects with conformers
├── CHEMBL240_hERG_with_3D.pkl       # Integrated dataset with 3D structures
├── herg_rf_model_v1.joblib          # Trained Random Forest model
├── herg_rf_features_v1.joblib       # Descriptor set used for inference
└── README.md
```

---

## Limitations

- Predictions depend on successful 3D conformer generation  
- Conformational sampling is limited to a finite number of embeddings  
- Model applicability is restricted to chemical space represented in the training dataset  

---

## Disclaimer

This model is intended for research and exploratory purposes only and should not be used as a substitute for experimental evaluation or regulatory decision-making.

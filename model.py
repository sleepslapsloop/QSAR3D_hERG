import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, Descriptors3D
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
import joblib
import numpy as np
from sklearn.metrics import pairwise_distances

def smiles_to_mol(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        return mol
    except:
        return None

def embed_3d(mol):
    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        AllChem.EmbedMolecule(mol, params)
        AllChem.UFFOptimizeMolecule(mol)
        return mol
    except:
        return None
    
def compute_descriptors(mol):
    if mol is None:
        return None
    if mol.GetNumConformers() == 0:
        return None

    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotB": Descriptors.NumRotatableBonds(mol),
        "RadiusGyr": Descriptors3D.RadiusOfGyration(mol),
        "Asphericity": Descriptors3D.Asphericity(mol),
        "Spherocity": Descriptors3D.SpherocityIndex(mol),
        "PMI1": Descriptors3D.PMI1(mol),
        "PMI2": Descriptors3D.PMI2(mol),
        "PMI3": Descriptors3D.PMI3(mol),
    }

def loadModel(model_path, features_path):
    model = joblib.load(model_path)
    features = joblib.load(features_path)
    return model, features

def predict_herg_blocker(smiles, model, features):
    probabilities = []
    for i in range(15):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.UFFOptimizeMolecule(mol)

        x = compute_descriptors(mol)
        X_new = pd.DataFrame([x])[features]

        prob = model.predict_proba(X_new)[0, 1]
        probabilities.append(prob)
        #print("Predicted hERG blocker probability:", prob)
    
    avg_prob = np.mean(probabilities)
    prob_max = np.max(probabilities)
    prob_min = np.min(probabilities)
    prob_std = np.std(probabilities)
    print("Minimum predicted probability:", prob_min)
    print("Average predicted hERG blocker probability:", avg_prob)
    print("Maximum predicted probability:", prob_max)
    print("Standard deviation of predictions:", prob_std)

"""
df = pd.read_csv("hERG_activity.csv", sep=";")

df = df[
    df["Smiles"].notna() &
    df["pChEMBL Value"].notna()
].copy()

df = df[df["Standard Type"].isin(["IC50", "Ki"])].copy()
df = df[df["Standard Units"] == "nM"].copy()
df = df[df["pChEMBL Value"].notna()].copy()

df = df.rename(columns={
    "Smiles": "SMILES",
    "pChEMBL Value": "pIC50"
})

df = df[["SMILES", "pIC50"]]

df["Mol"] = df["SMILES"].apply(smiles_to_mol)

print("Before:", len(df))
df = df[df["Mol"].notna()].copy()
print("After:", len(df))

df["Mol3D"] = df["Mol"].apply(embed_3d)

print("3D success:", df["Mol3D"].notna().sum()) #Extremely long computation

df.to_pickle("CHEMBL240_3D_molecules.pkl")
df = df.drop(columns=["Mol"])
df.to_pickle("CHEMBL240_hERG_with_3D.pkl") #no mol, only mol3d

df = df.dropna(subset=["Mol3D"]).copy()

X = pd.DataFrame(df["Mol3D"].apply(compute_descriptors).tolist())
y = (df["pIC50"] >= 5.5).astype(int)

X = X.reset_index(drop=True)
df = df.reset_index(drop=True)

y = (df["pIC50"] >= 5.5).astype(int)
X.to_csv("CHEMBL240_3D_descriptors.csv", index=False)

X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    stratify=y,
    random_state=42
)

model = RandomForestClassifier(
    n_estimators=800,
    max_depth=14,
    min_samples_leaf=8,
    n_jobs=-1,
    random_state=42
)

model.fit(X_train, y_train)

y_pred = model.predict_proba(X_test)[:, 1]
print("ROC AUC:", roc_auc_score(y_test, y_pred))

pd.Series(
    model.feature_importances_,
    index=X.columns
).sort_values(ascending=False)

joblib.dump(model, "herg_rf_model_v1.joblib")
joblib.dump(list(X.columns), "herg_rf_features_v1.joblib")
"""

model, features = loadModel("./herg_rf_model_v1.joblib", "./herg_rf_features_v1.joblib")

smiles = "Fc1ccc(cc1)Cn2c5ccccc5nc2NC4CCN(CCc3ccc(OC)cc3)CC4"  # astemizole

predict_herg_blocker(smiles, model, features)

"""
# Distance from new compound to training set centroid
centroid = X_train.mean().values
dist = np.linalg.norm(X_new.values - centroid)

# Normalize (rough heuristic)
max_dist = np.percentile(
    np.linalg.norm(X_train.values - centroid, axis=1), 95
)

confidence_ad = max(0, 1 - dist / max_dist)

print("AD confidence score:", confidence_ad)
"""
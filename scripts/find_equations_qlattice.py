# %%
import feyn
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

# %%
data = pd.read_table("../results/sym_reg/genes_for_reg.tsv")

UNWANTED_COLUMNS_FOR_TRAINING = ["run", "phenotype"]
data = data[data.columns.difference(UNWANTED_COLUMNS_FOR_TRAINING)].dropna()

# Let's record the categorical data types in our dataset (note features will be treated as numerical by default).
stypes = {}
for f in data.columns:
    if data[f].dtype == "object":
        stypes[f] = "c"

# %%
# Set random seed for reproducibility
random_seed = 1024

# Define the target variable
target = "phenotype_reg"

# Split
train, test = train_test_split(
    data, test_size=0.33, stratify=data[target], random_state=random_seed
)
# %%
ql = feyn.QLattice(random_seed=random_seed)

# %%
models = ql.auto_run(
    data=train,
    output_name=target,
    kind="classification",
    stypes=stypes,
    n_epochs=20,
)
# %%
best = models[0]
best.plot(train, test)
# %%
best.sympify()
# %%
best.plot_signal(train)

# %%
best.save("../results/sym_reg/qlattice_model_all_genes.json")

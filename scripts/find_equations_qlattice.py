import feyn
import pandas as pd
from sklearn.model_selection import train_test_split


def get_train_test_types(dataset, random_seed=1024, target="phenotype_reg"):
    UNWANTED_COLUMNS_FOR_TRAINING = ["run", "phenotype"]
    data = dataset[dataset.columns.difference(UNWANTED_COLUMNS_FOR_TRAINING)].dropna()

    # Let's record the categorical data types in our dataset (note features will be treated as numerical by default).
    stypes = {}
    for f in data.columns:
        if data[f].dtype == "object":
            stypes[f] = "c"

    # Split
    train, test = train_test_split(
        data, test_size=0.33, stratify=data[target], random_state=random_seed
    )

    return train, test, stypes


def get_best_model(
    training_data, stypes, priors, target="phenotype_reg", epochs=20, random_seed=1024
):

    ql = feyn.QLattice(random_seed=random_seed)
    ql.update_priors(priors)

    models = ql.auto_run(
        data=training_data,
        output_name=target,
        kind="classification",
        stypes=stypes,
        n_epochs=epochs,
    )

    best = models[0]

    return best


def save_model(model, train, test, filename):
    model.plot(train, test, filename=f"../results/sym_reg/{filename}_summary.html")
    model.plot_signal(train, filename=f"../results/sym_reg/{filename}_signal.svg")
    model.save(f"../results/sym_reg/{filename}_model.json")


def run_models():

    three_gene_data = pd.read_table("../results/sym_reg/selected_genes_for_reg.tsv")
    full_data = pd.read_table("../results/sym_reg/genes_for_reg.tsv")

    threeg_train, threeg_test, threeg_types = get_train_test_types(three_gene_data)
    full_train, full_test, full_types = get_train_test_types(full_data)

    threeg_priors = feyn.tools.estimate_priors(threeg_train, "phenotype_reg", floor=0.1)
    full_priors = feyn.tools.estimate_priors(full_train, "phenotype_reg", floor=0.1)

    threeg_model = get_best_model(threeg_train, threeg_types, threeg_priors)
    full_model = get_best_model(full_train, full_types, full_priors)

    save_model(threeg_model, threeg_train, threeg_test, "three_gene")
    save_model(full_model, full_train, full_test, "all_genes")


if __name__ == "__main__":
    run_models()

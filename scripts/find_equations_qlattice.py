# %%
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


def get_models(
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

    return models


def save_model(model, train, test, filename):
    model.plot(train, test, filename=f"../results/sym_reg/{filename}_summary.html")
    model.plot_signal(train, filename=f"../results/sym_reg/{filename}_signal.svg")
    model.save(f"../results/sym_reg/{filename}_model.json")


# From https://github.com/abzu-ai/QLattice-clinical-omics/blob/main/notebooks/functions.py
def get_models_table(models, train, test, model_name):
    model_list = []
    auc_list_train = []
    auc_list_test = []
    bic_list = []
    accuracy_train = []
    accuracy_test = []
    feat_list = []
    function_list = []
    loss_list = []
    i = 0
    for x in models:
        model_list.append(str(i))
        auc_list_train.append(str(x.roc_auc_score(train).round(2)))
        auc_list_test.append(str(x.roc_auc_score(test).round(2)))
        accuracy_train.append(str(x.accuracy_score(train).round(2)))
        accuracy_test.append(str(x.accuracy_score(test).round(2)))
        bic_list.append(str(x.bic.round(2)))
        feat_list.append(len(x.features))
        function_list.append(
            x.sympify(symbolic_lr=False, symbolic_cat=True, include_weights=False)
        )
        loss_list.append(x.loss_value)
        i += 1
    df = pd.DataFrame(
        list(
            zip(
                model_list,
                auc_list_train,
                auc_list_test,
                accuracy_train,
                accuracy_test,
                bic_list,
                feat_list,
                function_list,
                loss_list,
            )
        ),
        columns=[
            "Model",
            "AUC Train",
            "AUC Test",
            "Accuracy Train",
            "Accuracy Test",
            "BIC",
            "N. Features",
            "Functional form",
            "Loss",
        ],
    )

    df.to_csv(f"../results/sym_reg/{model_name}_table.csv", index=False)


def run_models():

    two_gene_data = pd.read_table("../results/sym_reg/selected_genes_for_reg.tsv")
    full_data = pd.read_table("../results/sym_reg/genes_for_reg.tsv")

    twog_train, twog_test, twog_types = get_train_test_types(two_gene_data)
    full_train, full_test, full_types = get_train_test_types(full_data)

    twog_priors = feyn.tools.estimate_priors(twog_train, "phenotype_reg", floor=0.1)
    full_priors = feyn.tools.estimate_priors(full_train, "phenotype_reg", floor=0.1)

    twog_models = get_models(twog_train, twog_types, twog_priors)
    full_models = get_models(full_train, full_types, full_priors)

    get_models_table(twog_models, twog_train, twog_test, "two_gene")
    get_models_table(full_models, full_train, full_test, "all_genes")

    save_model(twog_models[0], twog_train, twog_test, "two_gene")
    save_model(full_models[0], full_train, full_test, "all_genes")


if __name__ == "__main__":
    run_models()

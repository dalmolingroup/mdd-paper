# %%
import pandas as pd
import pysr

# %%
df = pd.read_table("../results/sym_reg/selected_genes_for_reg.tsv").dropna()

BINARY_OPERATORS = [
    "+",
    "-",
    "*",
    "/",
    "^",
]

UNARY_OPERATORS = [
    "square",
    "cube",
    "exp",
    "log",
]

PREVENT_NESTED_UNARY_OPERATORS = {
    unary_operator: {unary_operator: 0} for unary_operator in UNARY_OPERATORS
}

UNWANTED_COLUMNS_FOR_TRAINING = [
    "run",
    "phenotype",
    "phenotype_reg",
    "gender",
    "region",
]

# %%
regression_model = pysr.PySRRegressor(
    binary_operators=BINARY_OPERATORS,
    unary_operators=UNARY_OPERATORS,
    nested_constraints=PREVENT_NESTED_UNARY_OPERATORS,
    maxsize=8,  # default=20
    maxdepth=6,  # default=None
    complexity_of_constants=4,  # default=1
    precision=64,  # default=32
    equation_file=f"equations.csv",
    warm_start=False,
    random_state=1024,
)

# %%
regression_model.fit(
    X=df[df.columns.difference(UNWANTED_COLUMNS_FOR_TRAINING)],
    y=df["phenotype_reg"],
)


# %%
# Gambiarra
# Penalizes predictions when the sign differs
CUSTOM_SIGN_LOSS = "loss(x, y) = (x < 0) != (y < 0)"

regression_model_with_custom_loss = pysr.PySRRegressor(
    binary_operators=BINARY_OPERATORS,
    unary_operators=UNARY_OPERATORS,
    nested_constraints=PREVENT_NESTED_UNARY_OPERATORS,
    maxsize=8,  # default=20
    maxdepth=6,  # default=None
    complexity_of_constants=4,  # default=1
    precision=64,  # default=32
    equation_file=f"equations2.csv",
    warm_start=False,
    loss=CUSTOM_SIGN_LOSS,
)

regression_model_with_custom_loss.fit(
    X=df[df.columns.difference(UNWANTED_COLUMNS_FOR_TRAINING)],
    y=df["phenotype_reg"],
)

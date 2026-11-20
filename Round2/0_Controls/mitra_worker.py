"""Mitra (AutoGluon) in its own environment (Code/.venv-mitra: AutoGluon needs pandas < 3 and torch 2.10, which the
Round1 environment does not have). Reads Xtr, ytr, Xte and the task from an .npz file, writes predictions to an
.npz file. Called by common.py; not a step script.

  .venv-mitra/bin/python Round2/0_Controls/mitra_worker.py in.npz out.npz
"""

import sys
import tempfile

import numpy as np
import pandas as pd
from autogluon.tabular import TabularPredictor

d = np.load(sys.argv[1], allow_pickle=False)
task = str(d["task"])
cols = [f"f{i}" for i in range(d["Xtr"].shape[1])]
tr = pd.DataFrame(d["Xtr"], columns=cols).assign(y=d["ytr"])
te = pd.DataFrame(d["Xte"], columns=cols)
with tempfile.TemporaryDirectory() as tmp:
    p = TabularPredictor(label="y", problem_type="binary" if task == "classification" else "regression", path=tmp, verbosity=0)
    p.fit(tr, hyperparameters={"MITRA": {"fine_tune": False}}, ag_args_fit={"num_gpus": 0})
    out = p.predict_proba(te)[1].to_numpy() if task == "classification" else p.predict(te).to_numpy()
np.savez(sys.argv[2], pred=out)

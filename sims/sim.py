import geonomics as gnx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.random.seed(42)

def make_unif_array(n):
    """Makes a square array of ones, size n x n cells."""
    array = np.ones((n,n))
    return array

mod = gnx.make_model(parameters="GNX_params.py", verbose=True)
mod.run()
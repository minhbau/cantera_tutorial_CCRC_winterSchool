#! simple test to ensure that cantera importeds and runs properly
import numpy as np
import cantera as ct

# check install path
print(ct.__path__)
# check cantera version
print(ct.__version__)

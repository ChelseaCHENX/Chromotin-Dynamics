import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas
import bioframe
import cooltools
import cooler
from math import ceil, floor

def resol2bin(resol='50kb'):
    return int(resol.strip('kb')) * 1000
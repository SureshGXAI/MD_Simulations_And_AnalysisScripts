#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm

fig = plt.subplots(figsize=(12, 6))

hbond = pd.read_csv('hbond_complex.dat', sep='[ ]', header=None, names=["Time", "Distance"])


data = norm.rvs(hbond['Distance'])
 
ax = sns.distplot(data, bins=50, kde=True, color='green', hist_kws={"linewidth": 15,'alpha':1})

plt.xlabel('Hbond Distance ($\AA$)', fontsize=15)
plt.ylabel('Frequency', fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0, 6)
plt.ylim(0, 1)
plt.savefig('hbond.jpg')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from gimmemotifs.maelstrom import MaelstromResult
import matplotlib.pyplot as plt
import sys

folder = sys.argv[1]
mr = MaelstromResult(folder)
mr.plot_heatmap(threshold=3, indirect=False)
plt.savefig(folder + "/heatmap.pdf", dpi=300)

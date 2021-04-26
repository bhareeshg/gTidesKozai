import matplotlib
matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import rebound
from matplotlib.ticker import FormatStrFormatter

print(sys.argv)
datbs=pd.read_csv("test.txt");

print(datbs);

cmap = plt.cm.get_cmap('RdYlBu');
fig, axs = plt.subplots(2,2);

axs[0][0].plot(datbs.time,(1-datbs.e1),linewidth=0.6);
axs[0][1].plot(datbs.time,np.multiply(datbs.a1,1-datbs.e1),linewidth=0.6,label="q");
axs[0][1].plot(datbs.time,np.multiply(datbs.a1,1),linewidth=0.6,label="a");

axs[0][1].legend();
axs[1][0].plot(datbs.time,datbs.i1*180/np.pi,linewidth=0.6);
axs[1][1].plot(datbs.time,datbs.s1,linewidth=0.6,label="Planet");
axs[1][1].plot(datbs.time,datbs.s0,linewidth=0.6,label="Star");
axs[1][1].legend();

axs[0][0].set_ylabel(r'$e$')
axs[1][0].set_ylabel(r'$i^\circ$')
axs[0][1].set_ylabel(r'$a[AU]$')
axs[1][1].set_ylabel(r'$spin(days)$')

axs[1][0].set_xlabel(r'$time$')
axs[1][1].set_xlabel(r'$time$')

fig.subplots_adjust(wspace=.3)
axs[0][0].set_xlim(1e6,3.2e9);
axs[1][0].set_xlim(1e6,3.2e9);
axs[0][1].set_xlim(1e6,3.2e9);
axs[1][1].set_xlim(1e6,3.2e9);

axs[0][0].set_xscale("log");
axs[0][1].set_xscale("log");
axs[1][0].set_xscale("log");
axs[1][1].set_xscale("log");

axs[0][0].set_yscale("log");
axs[0][1].set_yscale("log");
#axs[1][0].set_yscale("log");
axs[1][1].set_yscale("log");


sp=[axs[0][0], axs[0][1]]
plt.setp([a.get_xticklabels() for a in sp], visible=False)

filename="comp_int.pdf";
plt.savefig(filename)


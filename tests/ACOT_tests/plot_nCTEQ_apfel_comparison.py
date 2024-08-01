import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

path = "/home/users/p_riss01/codes/apfelxx_ACOT/tests/ACOT_tests"
file_nCTEQ = path + "/nCTEQ_results/NC_F2_SACOT.csv"
file_apfel = path + "/APFELxx_results/NC_F2_SACOT.csv"

colors = ["orange","red","royalblue","darkgreen"]

keys = ["NCF2"]

#prepare data
nCTEQ = pd.read_csv(file_nCTEQ)
apfel = pd.read_csv(file_apfel)

Q2 = nCTEQ[nCTEQ["x"]==nCTEQ["x"].iloc[0]]["Q2"].to_numpy()
x = nCTEQ[apfel["Q2"]==Q2[0]]["x"].to_numpy()

#plot each structure function
for key in keys:
  fig,ax = plt.subplots(2,1)
  tot = ax[0]
  ratio = ax[1]
  for i,Q2i in enumerate(Q2):
    val_nCTEQ = nCTEQ[nCTEQ["Q2"]==Q2i][key]
    val_apfel = apfel[apfel["Q2"]==Q2i][key]
    tot.plot(x,val_nCTEQ,color=colors[i],label='Q={Q:.2f}'.format(Q=np.sqrt(Q2i)))
    tot.plot(x,val_apfel,ls='dashed',color=colors[i])
    ratio.plot(x,val_nCTEQ/val_apfel,color=colors[i])
    # print("Q:",np.sqrt(Q2i))
    # print(F2_nCTEQ/F2_apfel)

  tot.legend()
  tot.set_yscale('log')
  #symmetric ylim around 1
  yabs_max = abs(max(ratio.get_ylim(), key=abs)-1)
  print('yabs_max:',yabs_max)
  #clip plotting range as decimal precision of the output is set to 1e-10
  yabs_max = max(yabs_max*1.2,1e-9) 
  ratio.set_ylim(ymin=1-yabs_max, ymax=1+yabs_max)

  ratio.axhline(1,color='black',zorder=-1)

  for i in range(2):
    ax[i].set_xscale('log')
    ax[i].set_xlim(min(x),max(x))

  tot.set_title(key)
  # plt.suptitle(key)
  plt.show()  

#collect as PDF
#TODO
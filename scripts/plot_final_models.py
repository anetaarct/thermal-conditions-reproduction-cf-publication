from pathlib import Path
import sys
P=Path(__file__).resolve().parent.parent/'figures'

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
C=P.parent/'tables/final_models_20261019'; F=P.parent/'tables/primary_ZIyear_ML_20261021'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
blue='#0072B2';orange='#D55E00'
def curve(ax,d,label,xlabel,col):
    ax.plot(d.x,d.estimate,color=col,lw=1.7)
    ax.fill_between(d.x,d.CI_low,d.CI_high,color=col,alpha=.18,lw=0)
    ax.set(xlabel=xlabel,ylabel=label)
for kind,name,xlabel in [('temperature','Fig3_COM_M1_temperature','Mean spring temperature (°C)'),('year','Fig4_COM_M1_time','Year'),('date','Fig5_COM_M1_date','Relative laying date (days)')]:
    fig,ax=plt.subplots(2,2,figsize=(9.5,7.4),layout='constrained')
    for j,comp in enumerate(['mean','sd']):
        curve(ax[0,j],pd.read_csv(C/f'COM_{kind}_{comp}.csv'),['Mean clutch size (eggs)','Clutch conditional SD (eggs)'][j],xlabel,[blue,orange][j])
    for j,comp in enumerate(['success','phi']):
        d=pd.read_csv(F/f'ZIyear_{kind}_{comp}.csv')
        if kind=='temperature':d.x+=5.098595980711215
        label=['Expected egg-to-fledging success','Component precision (phi)'][j]
        if kind=='date':
            for age,col in [('YOUNG',blue),('OLD',orange)]:curve(ax[1,j],d[d.age_class==age],label,xlabel,col)
            ax[1,j].legend(ax[1,j].lines,['Yearling','Older'],frameon=False)
        else:curve(ax[1,j],d,label,xlabel,[blue,orange][j])
    for a,l in zip(ax.flat,'ABCD'):a.set_title(l,loc='left',fontweight='bold')
    fig.savefig(P/(name+'.pdf'));fig.savefig(P/(name+'.png'),dpi=150);plt.close(fig)
fig,ax=plt.subplots(2,2,figsize=(8,6),layout='constrained')
for a,d,label,l in zip(ax.flat,[pd.read_csv(C/'COM_age_mean.csv'),pd.read_csv(C/'COM_age_sd.csv'),pd.read_csv(F/'ZIyear_age_success.csv'),pd.read_csv(F/'ZIyear_age_phi.csv')],['Mean clutch size (eggs)','Clutch conditional SD (eggs)','Expected egg-to-fledging success','Component precision (phi)'],'ABCD'):
    for i,age in enumerate(['YOUNG','OLD']):
        z=d[d.age_class==age].iloc[0];a.errorbar(i,z.estimate,yerr=[[z.estimate-z.CI_low],[z.CI_high-z.estimate]],fmt='o',color=[blue,orange][i],capsize=4)
    a.set(xticks=[0,1],xticklabels=['Yearling','Older'],ylabel=label,xlim=(-.6,1.6));a.set_title(l,loc='left',fontweight='bold')
fig.savefig(P/'Fig6_COM_M1_age.pdf');fig.savefig(P/'Fig6_COM_M1_age.png',dpi=150)
print('Four primary ML COM-Poisson / M1 figures generated from archived prediction tables.')

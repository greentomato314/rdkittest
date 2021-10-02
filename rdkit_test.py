import pandas as pd
import os
import collections

from rdkit import Chem
from rdkit.Chem import Draw

max_n = 10
dic_name_mon = {'PSS':'C=Cc1ccc(S(=O)(=O)[O-])cc1','PAA':'C=CC(=O)O','PAAm':'C=CC(N)=O'} #ビニルモノマー
names = []
ns = []
smiles = []


for name,monomer in dic_name_mon.items():
    if not os.path.isdir(name):
        os.makedirs(name)
    sokusa = monomer[3:]
    for n in range(1,max_n+1):
        d = collections.deque(['C'])
        for i in range(n):
            d.appendleft('CC(')
            d.append(')'+sokusa)
        polymer = ''.join(d)
        rd = Chem.MolFromSmiles(polymer)
        # output
        names.append(name)
        ns.append(n)
        smiles.append(Chem.MolToSmiles(rd))
        Draw.MolToFile(rd,'{}/{}_polymer{}.png'.format(name,name,n),size=(300, 300))
    
df = pd.DataFrame({'names':names,'n':ns,'smiles':smiles})
df.to_csv('data.csv',index=False)
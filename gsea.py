import pandas as pd
import numpy as np

# gene set
g=pd.read_table('h.all.v5.1.symbols.gmt',index_col=0).loc['HALLMARK_GLYCOLYSIS'][1:].values

# gseapy sort high to low
d=pd.read_table('example.rnk',names=['lfc']).sort_values('lfc',ascending=0).reset_index()
d.columns='gene lfc'.split()

# overlap
d['hit']=(d['gene']).isin(g)

# Enrichment Score, data (index, columns: gene, lfc, hit), gene_set 
def ES(d,g):
	# norm
	hit= 1./ d[d.hit==1].lfc.abs().sum()
	es= ( d.hit*d.lfc.abs()*hit - (1.-d.hit)/len(d[d.hit==0]) ).cumsum()
	# es
	if np.abs(es.max()) > np.abs(es.min()):
		es1=es.max()
	else:
		es1=es.min()
	return es1
es= ES(d,g)

# random es, 1k
esr=[]

# overlap
o=len(d[d.hit==1])

# Permulation test
for i in range(1000):
	sam=np.random.choice(len(d),o)
	d['hit']=(d.index).isin(sam)
	esr.append(ES(d,g))

esr=pd.DataFrame(esr)

# Normalized ES
if es>0:
	nes = es/esr[esr>=0].mean()
else:
	nes = es/esr[esr<0].mean()

import pandas as pd

# glycolysis STAD c/s, ES 0.3
# gene set
s='HALLMARK_GLYCOLYSIS'
#s='HALLMARK_INFLAMMATORY_RESPONSE'
g=pd.read_table('/media/c/c/purity/gsea/db/h.all.v5.1.symbols.gmt',index_col=0).loc[s][1:].values

# gseapy sort high to low
d=pd.read_table('STAD.rnk',names=['lfc']).sort_values('lfc',ascending=0).reset_index()
d.columns='gene lfc'.split()

# overlap
d['hit']=(d['gene']).isin(g)

# data (index, columns: gene, lfc, hit), gene_set 
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

for i in range(1000):
	sam=np.random.choice(len(d),o)

	d['hit']=(d.index).isin(sam)

	esr.append(ES(d,g))

esr=pd.DataFrame(esr)

if es>0:
	nes = es/esr[esr>=0].mean()
else:
	nes = es/esr[esr<0].mean()

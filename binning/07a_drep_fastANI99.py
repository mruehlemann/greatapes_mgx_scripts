import pandas as pd

import drep
import drep.d_cluster.external
import drep.d_cluster.utils
import drep.d_cluster.compare_utils


wd=drep.WorkDirectory.WorkDirectory(".")
Ndb = wd.get_db('Ndb')
id="primary_cluster"

Cdb = pd.DataFrame()
c2ret = {}
for name, ndb in Ndb.groupby(id):
	cdb, cluster_ret = drep.d_cluster.cluster_utils.genome_hierarchical_clustering(ndb, cluster=name, comp_method="fastANI", cluster_method="average")
	cdb[id] = name
	Cdb = pd.concat([Cdb,cdb], ignore_index=True)
	c2ret[name] = cluster_ret

wd.store_db(Cdb, "Cdb_99")

Cdb97 = pd.DataFrame()
c2ret97 = {}
for name, ndb in Ndb.groupby(id):
	cdb97, cluster_ret97 = drep.d_cluster.cluster_utils.genome_hierarchical_clustering(ndb, cluster=name, comp_method="fastANI", cluster_method="average", S_ani=0.97)
	cdb97[id] = name
	Cdb97 = pd.concat([Cdb97,cdb97], ignore_index=True)
	c2ret97[name] = cluster_ret97

wd.store_db(Cdb97, "Cdb_97")

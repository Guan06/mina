# MINA
### by guan@mpipz.mpg.de
## addTable

#### Arguments
* __minaObj__ an object of the class __mina__.
* __countTable__ a feature table (e.g. OTU table).

#### Values
* __minaObj$countTable__.

## addDesign

#### Arguments
* __minaObj__ an object of the class __mina__.
* __design__ an object of the class __designTable__.

#### Values
* __minaObj$designTable__.

## normalizeTable

#### Arguments
* __minaObj__ an object of the class __mina__. Must contain a count table.
* __method__ normalization method to be used.
* __depth__ depth to perform rarefection if __method="rarefy"__.

#### Values
* __minaObj$normTable__.

## ComDis

#### Arguments
* __minaObj__ an object of the class __mina__. Must contain a normalized feature
table.
* __method__ distance or dissimilarity method to be used.

#### Values
* __minaObj$dis__.

## dmr

#### Arguments
* __minaObj__ an object of the class __mina__. Must contain a distance or
dissimilarity table.
* __method__ dimmensionality reduction method to be used.
* __k__ number of dimensions to use.

#### Values

* __minaObj$dmr$points__ coordinates of samples in the projection.
* __minaObj$dmr$variance__ variance of each axis (if __method="pcoa"__).

## adjacency

#### Arguments
* __minaObj__ an object of the class __mina__. Must contain a normalized feature
table.
* __method__ correlation method to be used.

#### Values

* __minaObj$Adj__. Adjacency matrix between all features in the normalized
table.

## netCluster

#### Arguments
* __minaObj__ an object of the class __mina__. Must contain a adjacency matrix.

#### Values

* an object of the class __netCls__
* __minaObj$netCls$feature__
* __minaObj$netCls$cluster__

## multAdj

#### Arguments
* __countTable__ normalized feature table.
* __method__ correlation method to be used.
* __n__ size of the subsampled matrices.
* __m__ number of subsamples.

#### Values
* __multAObj__ an object of the class __multA__ which contains a list of
adjacency matrices from each subsampled feature table.

## multNet

#### Arguments
* __minaObj__ an object of the class __mina__. Must contain a feature table and
a corresponding design.
* __method__ correlation method to be used.
* __n__ size of the subsampled matrices.
* __m__ number of subsamples.
* __env__ name of the environmental variable to be used for network
reconstruction (must be found in the design table).

#### Values
* __minaObj$envMultA__ list of object of the class __multA__ which contains
subsampled adjacency matrices for every value of __env__.

## permNet

#### Arguments
* __minaObj__ an object of the class __mina__. Must contain a feature table and
a corresponding design.
* __method__ correlation method to be used.
* __n__ size of the permutated matrices.
* __m__ number of permutations.
* __env__ name of the environmental variable to be used for network
reconstruction (must be found in the design table).

#### Values
* __minaObj$envPermA__ list of object of the class __permA__ which contains
permutated adjacency matrices for every pair of values in __env__.

## distNet

#### Arguments
* __A__, __B__ adjacency matrices of two networks to be compared. Should have
the same dimensions and matching row names (node IDs).
* __method__ distance method to be employed.

#### Values

* __dist__ distance between networks represented in __A__ and __B__.

## netDistSig

#### Arguments
* __envMultAObj__ an object of the class __envMultA__ which contains a list of
bootstrapped adjacency matrices for a set of environmental factors.
* __envPermAObj__ an object of the class __envPermA__ containing a list of
permutated adjacency matrices for each pairwise combination of environmental
factors.
* __method__ distance method to be employed.

#### Values
* __pvals__ list of P-values for each pair of environmental factors tested.
* __netDists__ list of network distances for each pair of environmental factors
tested.

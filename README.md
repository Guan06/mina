# MINA
### by guan@mpipz.mpg.de

## addTab

#### Arguments
* __mina__ an object of the class __mina__.
* __countTable__ a feature table (e.g. OTU table).

#### Values
* __mina$Tab__.

## addDes

#### Arguments
* __mina__ an object of the class __mina__.
* __design__ an object of the class __designTable__.

#### Values
* __mina$desTab__.

## normTab

#### Arguments
* __mina__ an object of the class __mina__. Must contain a count table.
* __method__ normalization method to be used.
* __depth__ depth to perform rarefection if __method="rarefy"__.

#### Values
* __mina$normTab__.

## comDis

#### Arguments
* __mina__ an object of the class __mina__. Must contain a normalized feature
table.
* __method__ distance or dissimilarity method to be used.

#### Values
* __mina$dis__.

## dmr

#### Arguments
* __mina__ an object of the class __mina__. Must contain a distance or
dissimilarity table.
* __method__ dimmensionality reduction method to be used.
* __k__ number of dimensions to use.

#### Values

* __mina$dmr$points__ coordinates of samples in the projection.
* __mina$dmr$variance__ variance of each axis (if __method="pcoa"__).

## adj

#### Arguments
* __mina__ an object of the class __mina__. Must contain a normalized feature
table.
* __method__ correlation method to be used.

#### Values

* __mina$Adj__. Adjacency matrix between all features in the normalized
table.

## netCl

#### Arguments
* __mina__ an object of the class __mina__. Must contain a adjacency matrix.

#### Values

* an object of the class __netCls__
* __mina$netCls$feature__
* __mina$netCls$cluster__

## multiNet

#### Arguments
* __mina__ an object of the class __mina__. Must contain a feature table and
a corresponding design.
* __method__ correlation method to be used.
* __n__ size of the subsampled matrices.
* __m__ number of subsamples.
* __env__ name of the environmental variable to be used for network
reconstruction (must be found in the design table). When env==NA,
construct network based on subset of all samples.

#### Values
* __mina$envMultA__ list of object of the class __multA__ which contains
subsampled adjacency matrices for every value of __env__.

## permNet

#### Arguments
* __mina__ an object of the class __mina__. Must contain a feature table and
a corresponding design.
* __method__ correlation method to be used.
* __n__ size of the permutated matrices.
* __m__ number of permutations.
* __env__ name of the environmental variable to be used for network
reconstruction (must be found in the design table).

#### Values
* __mina$envPermA__ list of object of the class __permA__ which contains
permutated adjacency matrices for every pair of values in __env__.

## netDis

#### Arguments
* __A__, __B__ adjacency matrices of two networks to be compared. Should have
the same dimensions and matching row names (node IDs).
* __method__ distance method to be employed.

#### Values

* __dist__ distance between networks represented in __A__ and __B__.

## netSig

#### Arguments
* __envMultA__ an object of the class __envMultA__ which contains a list of
bootstrapped adjacency matrices for a set of environmental factors.
* __envPermA__ an object of the class __envPermA__ containing a list of
permutated adjacency matrices for each pairwise combination of environmental
factors.
* __method__ distance method to be employed.

#### Values
* __pvals__ list of P-values for each pair of environmental factors tested.
* __netDists__ list of network distances for each pair of environmental factors
tested.

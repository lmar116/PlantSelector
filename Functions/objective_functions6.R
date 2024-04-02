## pollinator abundance function
abundance.ulb <- function(m) sum(rowSums(ULB[,m]))
abundance.natag <- function(m) sum(rowSums(NATAG[,m]))
abundance.combi <- function(m) sum(rowSums(COMBI[,m]))
abundance.user <- function(m) sum(rowSums(USERdf[,m]))

## pollinator richness function
richness.ulb <- function(m) sum(rowSums(ULB[,m]>0)>0)
richness.natag <- function(m) sum(rowSums(NATAG[,m]>0)>0)
richness.combi <- function(m) sum(rowSums(COMBI[,m]>0)>0)
richness.user <- function(m) sum(rowSums(USERdf[,m]>0)>0)

## combination functions
abundance.ulb.richness.ulb <- function(m)
  abundance.ulb(m)*richness.ulb(m)
abundance.natag.richness.natag <- function(m)
  abundance.natag(m)*richness.natag(m)
abundance.combi.richness.combi <- function(m)
  abundance.combi(m)*richness.combi(m)
abundance.user.richness.user <- function(m) 
  abundance.user(m)*richness.user(m)

## Xtree function for the calculation of FD
## 17/1/03. Written by Jens Schumacher

Xtree <- function(h)
## evaluate species branch matrix (sensu Petchey & Gaston 2002) from a dendrogram
## The object returned by Xtree() is a list containing a species by branch matrix (H1) and a branch length
## vector (h2)
## see http://onlinelibrary.wiley.com/doi/10.1046/j.1461-0248.2002.00339.x/abstract

## tested for results of hclust and agnes
## hclust - hierarchical clustering
## agnes - agglomerative clustering

## used components:
## merge - history of cluster merging
## height - actual heights at merging
## order - permutation to achieve nice output (needed only for agnes)
{
    
    species.names <- h$labels
    
    
    H1 <- matrix(0, length(h$order), 2 * length(h$order) - 2)
    l <- vector("numeric", 2 * length(h$order) - 2)
    for(i in 1:(length(h$order) - 1)) {
        # evaluate branch lengths
        #
        if(h$merge[i, 1] < 0) {
             l[2 * i - 1] <- h$height[order(h$height)[i]]
            H1[ - h$merge[i, 1], 2 * i - 1] <- 1
        }
        else {
            l[2 * i - 1] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 1]]]
            H1[, 2 * i - 1] <- H1[, 2 * h$merge[i, 1] - 1] + H1[
            , 2 * h$merge[i, 1]]
        }
        if(h$merge[i, 2] < 0) {
            l[2 * i] <- h$height[order(h$height)[i]]
            H1[ - h$merge[i, 2], 2 * i] <- 1
        }
        else {
            l[2 * i] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 2]]]
            H1[, 2 * i] <- H1[, 2 * h$merge[i, 2] - 1] + H1[, 2 *
            h$merge[i, 2]]
        }
    }
    dimnames(H1) <- list(species.names,NULL)
    list(h2.prime=l, H1=H1)
    ## l contains the length of all the tiny branches
    ## H1: each row represents one species, each column represents one branch
    ##     1 indicates that a branch is part of the pathway from species to top of the dendrogram
    ##     0 otherwise
}


## compute score for model m
score <- function(f.list, m) {
  g <- function(f.list, m) sapply(f.list, function(y) y(m))
  prod(g(f.list, m))
}

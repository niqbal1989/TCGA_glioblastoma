
MicroRNAfilter <- function(MicroRNAid, index){#use index to return list of relevant MicroRNAs
    cheese <- list()
    for(i in 1:length(index)){
        cheese[i] <- MicroRNAid[index[i]]
    }
    return(cheese)
}
subtypecounter <- function(x) {#count number per each subtype
    one<-0
    two <- 0
    three <- 0
    four <- 0
    for(i in 1:length(x)) {
        if (x[i]==1) one <- one+1
        if (x[i]==2) two <- two+1
        if (x[i]==3) three <- three+1
        if (x[i]==4) four <- four+1
    }
    return(c(one,two,three,four))
}
#this function returns a subtype list which also definitively has microRNA expression data for each sample
MicroRNASamplefilter <- function(MicroRNAsamplelist, subtypesamplelist){
    w <- list()
    x <- list()
    y <- list()
    z <- list()
    f <- 1
    for (i in 1:(length(MicroRNAsamplelist))){
        x[i] <- substr(MicroRNAsamplelist[[i]]$barcode[1],9,12)
    }
    for (i in 1:(length(subtypesamplelist))){
        y[i] <- substr(subtypesamplelist[i],9,12)
    }
    z <- intersect(x,y)
    for (i in 1:length(subtypesamplelist)){
        if (z[f] == substr(subtypesamplelist[i],9,12)){
            w[f] <- subtypesamplelist[i]
            f <- f+1
        }
    }
    return(w)

}
subtypesorter <- function(x,number) { #use clusterindex for x
    a <- 1
    cheese <- list()
    for (i in 1:length(x)) {
        if (x[i] == number) {
            cheese[a] <- names(x[i])
            a <- a+1
        }
    }
   return(cheese)
}
positivecounter <- function(x){
    counter <- 0
    for (i in 1:length(x[,3])) {
        if (x[i,3] > 0) counter <- counter+1
    }
    return(counter)
}


#positive sillhouette takes the a silhouette class of type sil of data and the clusterindex aka the consensus class vector for k clusters. The output is the clusterindex of only positive values
positivesilhouette <- function(x, index){
    g<-1
    y<-integer()
    samples <- list()
    for(i in 1:length(index)){
        if (x[i,3]>0){
           y[g]<-index[i]
           samples<- c(samples,names(index[i]))
           g <- g+1
        }
    }
    names(y) <- samples
    return(y)
}
#the following function takes list of samples in each subtypes, x, and the original list of genes and creates a dataframe where the rows are the gene expression while the columns are the samples in each subtype x
subtypedfmRNA <- function(x, genes){
    a <- 1
    b <- length(rownames(genes))
    z <- matrix(nrow=b)
    for (i in 1:length(colnames(genes))){
        if (colnames(genes)[i] == x[a]){
            a <- a+1
            z <- cbind(z, genes[,i])
        }
    }
    z<-z[,-1]
    colnames(z) <- x
    rownames(z) <- rownames(genes)
    return(z)
}
#this function creates a data matrix where the columes are sample barcodes and the rows are human microRNAs
subtypedfmiRNA <- function(miRNAlist, data, samplelist){#miRNAlist is list of microRNAs that we need to use. data is microRNAdata that was read in. samplelist is list of samples per subtype
    a <- 1
    b <- length(miRNAlist)
    c <- length(samplelist)
    d <- 1
    e <- length(data)
    f <- 1
    z <- matrix(nrow=b,ncol=c)
    rownames(z) <- miRNAlist
    colnames(z) <- samplelist
    for (i in 1:e){
        d <- 1
        if (substr(data[[i]]$barcode[1],9,12) == substr(samplelist[a],9,12)){
            for (j in 1:length(data[[1]]$miRNA.id)){
                if (miRNAlist[d]== data[[i]]$miRNA.id[j]){
                    z[d,a] <- data[[i]]$value[j]

                    d <- d+1
                }
            }
            a <- a+1
        }
    }
    return(z)
}
targetlist <- function(names, n){
    maturemirna <- vector()
    targetsymbol <- vector()
    if (n==1){#for validated targets
        for(i in 1:length(names)){
            test <- get.multimir(mirna=names[i])
            print(i)
            if (length(test)==0){
                next
            }else{
                maturemirna <- c(maturemirna, as.character(test$validated$mature_mirna_id))
                targetsymbol <- c(targetsymbol, as.character(test$validated$target_symbol))
            }
        }
    }
    if (n==2){
        for (i in 1:length(names)){
            test <- get.multimir(org='hsa', target= names[i], table='predicted')
            print(i)
            if (length(test)==0){
                next
            }else{
                maturemirna <- c(maturemirna, as.character(test$predicted$mature_mirna_id))
                targetsymbol <- c(targetsymbol, as.character(test$predicted$target_symbol))
            }
        }
    }
    stuff <- cbind(maturemirna, targetsymbol)
    return(stuff)
}
miRtargetbasefilter <- function(mirna,tardatabase,genes){#mirna is list of human microRNAs, tardatabase is the target database, genes is the list of gene names
    b <- intersect(mirna,tardatabase[,1])#find common miRNAs between target database and list of microRNAs
    f <- intersect(genes,tardatabase[,2])#find common genes in target database and list of genes
    database <- list()#create empty databases
    for (i in 1:length(b)){#for loop for length of number of common miRNAs
        print(i)
        print("l")
        z <- data.matrix(tardatabase[which(tardatabase[,1]==b[i]),])#create data frame from original database that is only those microRNAs from mirna
        if(length(z)==0)
            next
        if(length(z)==2){
            if (any(which(f == z[2]))){
            database$miRNA <- c(database$miRNA, as.character(z[1]))
            database$Targetgene <- c(database$Targetgene,as.character(z[2]))
        }
        }else{
            for (j in 1:length(z[,2])){#for loop for length of data matrix z
                if (any(which(f == z[j,2]))){#find which genes, common from gene and target database are in the z data matrix. When they are common, specifically, when there is a target interaction that has elements from both genes list and mirna list, stick entry into new database
                    database$miRNA <- c(database$miRNA, as.character(z[j,1]))
                    database$Targetgene <- c(database$Targetgene,as.character(z[j,2]))
                }
            }
    }
    }
    database<-data.frame(lapply(database,function(x)factor(unlist(x))))
    return(database)
}
miRtargetbasefilter2 <- function(mirna,tardatabase,genes){#mirna is list of human microRNAs, tardatabase is the target database, genes is the list of gene names. This is the original function from thesis.r. The new one kept returning a list/database of 0 so used this one instead and it worked
    b <- intersect(mirna,tardatabase[,1])#find common miRNAs between target database and list of microRNAs
    f <- intersect(genes,tardatabase[,2])#find common genes in target database and list of genes
    database <- list()#create empty databases
    for (i in 1:length(b)){#for loop for length of number of common miRNAs
        z <- tardatabase[which(tardatabase[,1]==b[i]),]#create data frame from original database that is only those microRNAs from mirna
        for (j in 1:length(z[,2])){#for loop for length of data matrix z
            if (any(which(f == z[j,2]))){#find which genes, common from gene and target database are in the z data matrix. When they are common, specifically, when there is a target interaction that has elements from both genes list and mirna list, stick entry into new database
                database$miRNA <- c(database$miRNA, as.character(z[j,1]))
                database$Targetgene <- c(database$Targetgene,as.character(z[j,2]))
            }
        }
    }

    return(database)
}
Networkformation2 <- function(miRNAexpr,targetdatabase){
    Network <- list()
    x <- rownames(miRNAexpr)
    x <- intersect(x,targetdatabase$miRNA)
    for (i in 1:length(x)){
        y <- which(targetdatabase$miRNA == x[i])#y is list of index values
        y <- targetdatabase[y,1:2]#y is data frame that has x[i] in column and has target genes in third column
        for (j in 1:length(y[,1])){
            z <- which(targetdatabase$Targetgene == y[j,2])
            z <- targetdatabase[z,1:2]#z is data frame that has all the same gene symbols in the third column and the miRNA's that target them in the first column
            for (k in 1:length(z[,1])){
                a <- which(Network$primarymiR== as.character(z[k,1]))
                b <- which(Network$gene==y[j,2])
                c <- which(Network$secondarymiR==x[i])
                d <- intersect(a,b)
                e <- intersect(d,c)
                if (length(e)==1)
                    next
                print(i)
                if (z[k,1] == x[i]) #don't want to record an miRNA's interaction with itself
                    next
                if (any(which(x ==z[k,1]))){
                Network$primarymiR <- c(Network$primarymiR, x[i])
                Network$gene <- c(Network$gene, as.character(y[j,2]))
                Network$secondarymiR <- c(Network$secondarymiR, as.character(z[k,1]))
            }
            }
            }
        }

    Network <- data.frame(sapply(Network,c))
    return(Network)
}
condinfonetwork <- function(network, genes, miRNA){#create dataframe with gene interactions and computed conditional informaiton for each interaction and one thousand expected conditional mutual information from permutating miRNA expression
    network["observeddeltaMI"] <- NA
    a <- rownames(miRNA)
    b <- rownames(genes)
    for (i in 1:length(network[,1])){
        c <- which(a == network[i,1])
        d <- which(b == network[i,2])
        e <- which(a == network[i,3])
        f <- miRNA[c,]
        g <- genes[d,]
        h <- miRNA[e,]
        discretef<-discretize(f, disc="equalwidth")
        discreteg <- discretize(g, disc="equalwidth")
        discreteh <- discretize(h, disc="equalwidth")
        network$observeddeltaMI[i] <- condinformation(discretef,discreteh,discreteg, method="emp")
    }
    for (i in 1:1000){
        network[sprintf("expecteddeltaMI%d",i)] <- NA
        x <- length(genes[1,])
        y <- shuffle(n=x, control=CTRL)
        print(i)
        for (j in 1:length(network[,1])){
            c <- which(a == network[j,1])
            d <- which(b == network[j,2])
            e <- which(a == network[j,3])
            f <- miRNA[c,]
            g <- genes[d,]
            g <- g[y]
            h <- miRNA[e,]
            discretef<-discretize(f, disc="equalwidth")
            discreteg <- discretize(g, disc="equalwidth")
            discreteh <- discretize(h, disc="equalwidth")
            network[j,(i+4)] <- condinformation(discretef,discreteh,discreteg, method="emp")
        }
    }
    return(network)
}

condinfonetwork2 <- function(network, genes, miRNA){#compute observed conditional mutual information and then compute the first 200 of expected
    a <- rownames(miRNA)
    b <- rownames(genes)
    if (length(network[1,] == 3)){
    network["observeddeltaMI"] <- NA
    for (i in 1:length(network[,1])){
        c <- which(a == network[i,1])
        d <- which(b == network[i,2])
        e <- which(a == network[i,3])
        f <- miRNA[c,]
        g <- genes[d,]
        h <- miRNA[e,]
        discretef<-discretize(f, disc="equalwidth")
        discreteg <- discretize(g, disc="equalwidth")
        discreteh <- discretize(h, disc="equalwidth")
        network$observeddeltaMI[i] <- condinformation(discretef,discreteh,discreteg, method="emp")
    }
}
    if ((length(network[1,]))<200){
        for (i in 1:200){
        network[sprintf("expecteddeltaMI%d",i)] <- NA
        x <- length(genes[1,])
        y <- shuffle(n=x, control=CTRL)
        print(i)
        for (j in 1:length(network[,1])){
            c <- which(a == network[j,1])
            d <- which(b == network[j,2])
            e <- which(a == network[j,3])
            f <- miRNA[c,]
            g <- genes[d,]
            g <- g[y]
            h <- miRNA[e,]
            discretef<-discretize(f, disc="equalwidth")
            discreteg <- discretize(g, disc="equalwidth")
            discreteh <- discretize(h, disc="equalwidth")
            network[j,(i+4)] <- condinformation(discretef,discreteh,discreteg, method="emp")
        }
    }
    }
    return(network)
}
condinfonetwork3 <- function(network, genes, miRNA){#compute observed conditionalmutual information for columns 201 to 400
    a <- rownames(miRNA)
    b <- rownames(genes)

    if ((length(network[1,]))>200){
        for (i in 201:400){
        network[sprintf("expecteddeltaMI%d",i)] <- NA
        x <- length(genes[1,])
        y <- shuffle(n=x, control=CTRL)
        print(i)
        for (j in 1:length(network[,1])){
            c <- which(a == network[j,1])
            d <- which(b == network[j,2])
            e <- which(a == network[j,3])
            f <- miRNA[c,]
            g <- genes[d,]
            g <- g[y]
            h <- miRNA[e,]
            discretef<-discretize(f, disc="equalwidth")
            discreteg <- discretize(g, disc="equalwidth")
            discreteh <- discretize(h, disc="equalwidth")
            network[j,(i+4)] <- condinformation(discretef,discreteh,discreteg, method="emp")
        }
    }
    }
    return(network)
}
condinfonetwork4 <- function(network, genes, miRNA){#compute observed conditionalmutual information for columns 401 to 600
    a <- rownames(miRNA)
    b <- rownames(genes)

    if ((length(network[1,]))>400){
        for (i in 401:600){
        network[sprintf("expecteddeltaMI%d",i)] <- NA
        x <- length(genes[1,])
        y <- shuffle(n=x, control=CTRL)
        print(i)
        for (j in 1:length(network[,1])){
            c <- which(a == network[j,1])
            d <- which(b == network[j,2])
            e <- which(a == network[j,3])
            f <- miRNA[c,]
            g <- genes[d,]
            g <- g[y]
            h <- miRNA[e,]
            discretef<-discretize(f, disc="equalwidth")
            discreteg <- discretize(g, disc="equalwidth")
            discreteh <- discretize(h, disc="equalwidth")
            network[j,(i+4)] <- condinformation(discretef,discreteh,discreteg, method="emp")
        }
    }
    }
    return(network)
}
condinfonetwork5 <- function(network, genes, miRNA){#compute observed conditionalmutual information for columns 601 to 800
    a <- rownames(miRNA)
    b <- rownames(genes)

    if ((length(network[1,]))>600){
        for (i in 601:800){
        network[sprintf("expecteddeltaMI%d",i)] <- NA
        x <- length(genes[1,])
        y <- shuffle(n=x, control=CTRL)
        print(i)
        for (j in 1:length(network[,1])){
            c <- which(a == network[j,1])
            d <- which(b == network[j,2])
            e <- which(a == network[j,3])
            f <- miRNA[c,]
            g <- genes[d,]
            g <- g[y]
            h <- miRNA[e,]
            discretef<-discretize(f, disc="equalwidth")
            discreteg <- discretize(g, disc="equalwidth")
            discreteh <- discretize(h, disc="equalwidth")
            network[j,(i+4)] <- condinformation(discretef,discreteh,discreteg, method="emp")
        }
    }
    }
    return(network)
}
condinfonetwork6 <- function(network, genes, miRNA){#compute observed conditionalmutual information for columns 801 to 1000
    a <- rownames(miRNA)
    b <- rownames(genes)

    if ((length(network[1,]))>600){
        for (i in 801:1000){
        network[sprintf("expecteddeltaMI%d",i)] <- NA
        x <- length(genes[1,])
        y <- shuffle(n=x, control=CTRL)
        print(i)
        for (j in 1:length(network[,1])){
            c <- which(a == network[j,1])
            d <- which(b == network[j,2])
            e <- which(a == network[j,3])
            f <- miRNA[c,]
            g <- genes[d,]
            g <- g[y]
            h <- miRNA[e,]
            discretef<-discretize(f, disc="equalwidth")
            discreteg <- discretize(g, disc="equalwidth")
            discreteh <- discretize(h, disc="equalwidth")
            network[j,(i+4)] <- condinformation(discretef,discreteh,discreteg, method="emp")
        }
    }
    }
    return(network)
}#functions condinfonetwork2 through condinfornetwork6 are meant to be used in sequential order
organizenetwork <- function(network){#order network that hasn't been had conditional mutual information appended to network so that same miRNA interactions are listed together
    network2 <- list()
    network2$primarymiR <-list()
    network2$gene <- list()
    network2$secondarymiR <- list()
    for (i in 1:length(network[,1])){
        x <- which(network2$primarymiR == network[i,1])
        y <- which(network2$gene == network[i,2])
        z <- which(network2$secondarymiR == network[i,3])
        m <- intersect(x,y)
        n <- intersect(m,z)
        if (length(n) == 1)
            next
        b <- which(network[,1] == network[i,1])
        c <- which(network[,3] == network[i,3])
        d <- which(network[,3] == as.character(network[i,1]))
        e <- which(network[,1] == as.character(network[i,3]))
        f <- intersect(b,c)
        g <- intersect(d,e)
        h <- c(f,g)
        network3 <- network[h,]
        network2$primarymiR <- c(network2$primarymiR,as.character(network3$primarymiR))
        network2$gene <- c(network2$gene,as.character(network3$gene))
        network2$secondarymiR <- c(network2$secondarymiR,as.character(network3$secondarymiR))
    }
    network2 <- data.frame(sapply(network2,c))
    return(network2)
}
pvalueExp <- function(orderednetwork,condinfonetwork){#add pvalues to organized networks for pvalues ranging from 0 to 1 for experimentally validated interactions. This function actually uses the best method for determining the pvalue because it uses approximations to a normal distribution that have very accurate calculations for small values
    orderednetwork["pvalue"] <- NA
    for (i in 1:length(orderednetwork[,1])){
         c <- which(as.character(condinfonetwork[,1]) == as.character(orderednetwork[i,1]))
         d <- which(as.character(condinfonetwork[,2]) == as.character(orderednetwork[i,2]))
         e <- which(as.character(condinfonetwork[,3]) == as.character(orderednetwork[i,3]))
         f <- intersect(c,d)
         g <- intersect(e,f)
         x <- as.vector(t(condinfonetwork[g,5:1004]))
         xbar <- mean(x)
         s <- sd(x)
         y <- condinfonetwork[g,4]
         z <- (y - xbar)/(s)
         b <- pnorm(z)

         orderednetwork[i,4] <- b
     }
    return(orderednetwork)

}
Pvalue <- function(pvaluenetwork){#final interaction list with Pvalue for each miRNA-miRNA interaction. Combine pvalues using fisher's method
    pvaluenetwork2 <- pvaluenetwork[,-4]
    pvaluenetwork2 <- pvaluenetwork2[,-2]
    pvaluenetwork2 <- unique(pvaluenetwork2)
    pvaluenetwork2["Pvalue"] <- NA
    for (i in 1:length(pvaluenetwork2[,1])){
        x <- c()
        a <- which(pvaluenetwork[,1] == as.character(pvaluenetwork2[i,1]))
        b <- which(pvaluenetwork[,3] == as.character(pvaluenetwork2[i,2]))
        c <- intersect(a,b)
        for (j in 1:length(c)){
            d <- c[j]
            f <- pvaluenetwork[d,4]
            x <- c(x, log(f))
        }
        g <- (-2*sum(x))
        h <- (2*(length(c)))
        Pvalue <- 1-pchisq(g,df=h)
        pvaluenetwork2[i,3] <- Pvalue

    }
return(pvaluenetwork2)
}
#The two pvalue filter functions are in thesis.r, not in thesis2.r.
Pvaluefilter <- function(Pvaluenetwork){#filter out pvalues that are not less than .05
    index <- list()
    for (i in 1:length(Pvaluenetwork[,1])){
        if (Pvaluenetwork[i,3]<.05){
            index <- c(index, i)
        }
    }
    index <- as.numeric(index)
    Pvaluenetwork <- Pvaluenetwork[index,]
    Pvaluenetwork <- data.frame(lapply(Pvaluenetwork,function(x)factor(unlist(x))))
    return(Pvaluenetwork)
}

Pvaluefilter2 <- function(Pvaluenetwork){#filter out pvalues that are not less than exp(-5)
    index <- list()
    for (i in 1:length(Pvaluenetwork[,1])){
        if (Pvaluenetwork[i,3]< exp(-5)){
            index <- c(index, i)
        }
    }
    index <- as.numeric(index)
    Pvaluenetwork <- Pvaluenetwork[index,]
    Pvaluenetwork <- data.frame(lapply(Pvaluenetwork,function(x)factor(unlist(x))))
    return(Pvaluenetwork)
}
MicroRNA <- list()
setwd("C:/Users/Nadia/Documents/R/Expression-miRNA/UNC__H-miRNA_8x15K/Level_3")
txtfiles=list.files(pattern = 'unc.edu*')
for(i in 1:438){
    MicroRNA[i] <- lapply(txtfiles[i], read.table, header=TRUE, sep="\t", dec=".", quote="")
}
setwd("C:/Users/Nadia/Documents/R")
genes <- read.table("unifiedScaledFiltered.txt", header=TRUE, sep="\t", dec=".", quote="")
setwd("C:/Users/Nadia/Documents/glioblastoma")#final folder for all of this
genes2<- data.matrix(genes,rownames.force=NA)#convert data frame to numeric matrix
library(graph)
library(bitops)
library(RCurl)
library(XML)
library(BiocGenerics)
library(Biobase)
library(ALL)
library(ConsensusClusterPlus)
library(cluster)
library(XMLRPC)
results<-ConsensusClusterPlus(genes2, maxK=6, reps=1000, pItem=0.8, pFeature=1.0,clusterAlg="hc", title="consensuscluster", innerLinkage="average", finalLinkage="average", distance="pearson",seed=1363746070.17052 ,writeTable=TRUE, verbose=TRUE, plot="png") #use data matrix of gene expression with consensuscluster plus algorithm. This algorithm compares gene expression between cell lines and groups them together based on similarity.
#output of consensusclusterplus is list in which element of the list corresponds to results from kth clustering. For example, results[[2]] is the results for 2 cluster solution
#needed to add last argument for the clustering to happen without error
clusterindex<-results[[4]]$consensusClass#clusterindex is list of numbers of same length as the number of patient cell lines. For example, if first entry was 1, that means the first sample belongs to the cluster 1 in the 4 cluster solution
e <- dist(results[[4]]$consensusMatrix,method="euclidean")#create dissimilarity matrix, the function computes and returns distance matrix computed by using euclidean distance measure to compute the distances between rows of a data matrix
subtypesil <- silhouette(clusterindex,e)#silhouette function determines if sample belongs in cluster. Large silhouette value means observation is very well clustered. Small silhouette value, near 0, means observation lies between two clusters. Observations that are negative are probably placed in wrong cluster
clusterindex2<-positivesilhouette(subtypesil,clusterindex)#create index of only those samples that have positive silhouette width within their cluster
subtype1<-subtypesorter(clusterindex2,1)
subtype2<-subtypesorter(clusterindex2,2)
subtype3<-subtypesorter(clusterindex2,3)
subtype4<-subtypesorter(clusterindex2,4)#sort samples with positive silhouette values into subtypes
subtype1<-MicroRNASamplefilter(MicroRNA, subtype1)
subtype2<-MicroRNASamplefilter(MicroRNA, subtype2)
subtype3<-MicroRNASamplefilter(MicroRNA, subtype3)
subtype4<-MicroRNASamplefilter(MicroRNA, subtype4)#make sure that there is microRNA expression data for patient cell lines in each subtype, so each subtype list has only those samples that have microRNA expression data
MicroRNAid <- MicroRNA[[1]]$miRNA.id[]
MicroRNAid2 <- as.character(MicroRNAid)#list of all possible microRNAs from data
hsaindex <- grep("hsa",MicroRNAid2) #indices of human microRNAs
humanmiRNA <- MicroRNAfilter(MicroRNAid2,hsaindex) #create list of human microRNAs
Centroids <- read.table("ClaNC840_centroids2.txt", header=TRUE, sep="\t", dec=".", quote="")#dataset created from Verhaak paper, significance analysis of microarrays and receiver operating characteristic methods to identify marker genes for each subtype. ClaNC, nearest centroid-based classification algorithm, used to find signatures of each class. It is used to predict subtype validation
Centroids2 <- data.matrix(Centroids)
rownames(Centroids2) <- Centroids[,1]#rownames are genes
Centroids2 <- Centroids2[,-1:-2]#get rid of column with rownames and column with weird names, left with rownames and 4 columns of expression for each subtype
namesindex <- names(clusterindex2)#names of only those samples that have positive silhouette values
genes3 <- subtypedfmRNA(namesindex,genes3)#create dataframe, containing only gene expression information for only those samples that had positive silhouette values
centroidgeneset <- genes2
Centroids3 <- Centroids2 #retain original centroids
Class <- rep(NA, ncol(centroidgeneset))#Class is vector containing NA for number of columns in original scaled filtered gene set
    names(Class) <- colnames(centroidgeneset)#names in vector correspond sample names
    Features <- rownames(centroidgeneset)#Features is gene names
    Features <- Features[which(Features %in% rownames(Centroids3))] #find genes that are expressed both in Centroids and the data and reassign to the Featurespace
    centroidgeneset <- centroidgeneset[Features, ] #reducedataset such that only genes in Centroid and in gene data set are in the matrix
    Centroids3 <- Centroids3[Features, ]#do same thing as in last comment for centroid
    Correlation <- cor(centroidgeneset, Centroids3, use = "pairwise.complete.obs", method = "pearson")#correlation matrix correlation of each sample to certain subtype listed in centroids, or rather correlation to centroid
    for (i in 1:nrow(Correlation)) {
        Class[i] <- which.max(Correlation[i, ])
    } #for row i of the above correlation matrix, in Class[i], record which subtype(by column number had the highest value), so for each sample, record with which column has the greatest number. Class is vector of samples with corresponding, 1, 2,3 and 4 which indicate the highest correlation to a specific subtype
Proneural<-subtypesorter(Class,1) #Column 1 corresponds to proneural centroid. Thus, list created for samples that fall into proneural subtype
Neural<-subtypesorter(Class,2)#Column 2 corresponds to neural centroid. Thus, list created for samples that fall into neural subtype
Classical <- subtypesorter(Class,3)#Column 3 corresponds to classical centroid.Thus, list created for samples that fall into classical subtype
Mesenchymal <- subtypesorter(Class,4)#Column 4 corresponds to  mesenchymal centroid. Thus, list created for samples that fall into mesenchymal subtype
#determine which subtype is which by finding which subtype groups which had positive silhouettes have the greatest intersecting samples with the samples determined into subtypes by the centroid bases method
Neural<-intersect(Neural,subtype3)
Proneural<-intersect(subtype2,Proneural)
Mesenchymal<-intersect(subtype1,Mesenchymal)
Classical<-intersect(subtype4,Classical)
#determine significance of clusters after consensus clustering and having matching microRNA expression data, but before determining subtype types from centroids
subtype1genes <- subtypedfmRNA(subtype1,genes2) #create data frame where rows are genes and columns are samples
subtype2genes <- subtypedfmRNA(subtype2,genes2)
subtype3genes <- subtypedfmRNA(subtype3,genes2)
subtype4genes <- subtypedfmRNA(subtype4,genes2)
invert1 <- t(subtype1genes)
invert2 <- t(subtype2genes)
invert3 <- t(subtype3genes)
invert4 <- t(subtype4genes)
sigclustmatrix1 <- rbind(invert1,invert2)
label <- c(rep(1,length(subtype1)),rep(2,length(subtype2)))
pvalue12<-sigclust(sigclustmatrix1, nsim=1000, nrep=1, labflag=1, label=label, icovest=3)# significance between mesenchymal and proneural, "pvalnorm":6.6199e-22
sigclustmatrix2 <- rbind(invert1,invert3)
label <- c(rep(1,length(subtype1)),rep(2,length(subtype3)))
pvalue13 <- sigclust(sigclustmatrix2, nsim=1000,nrep=1, labflag=1, label=label, icovest=3)#significance between mesenchymal and neural, "pvalnorm":4.937556e-08
sigclustmatrix3 <- rbind(invert1,invert4)
label <- c(rep(1,length(subtype1)),rep(2,length(subtype4)))
pvalue14 <- sigclust(sigclustmatrix3, nsim=1000,nrep=1, labflag=1, label=label, icovest=3)#significance between mesenchymal and classical, "pvalnorm":3.499533e-10
sigclustmatrix4 <- rbind(invert2,invert3)
label <- c(rep(1,length(subtype2)),rep(2,length(subtype3)))
pvalue23 <- sigclust(sigclustmatrix4, nsim=1000,nrep=1, labflag=1, label=label, icovest=3)#significance between proneural and neural, "pvalnorm":0.05753577
sigclustmatrix5 <- rbind(invert2,invert4)
label <- c(rep(1,length(subtype2)),rep(2,length(subtype4)))
pvalue24 <- sigclust(sigclustmatrix5, nsim=1000,nrep=1, labflag=1, label=label, icovest=3)#significance between proneural and classical, "pvalnorm":9.928148e-13
sigclustmatrix6 <- rbind(invert3,invert4)
label <- c(rep(1,length(subtype3)),rep(2,length(subtype4)))
pvalue34 <- sigclust(sigclustmatrix6, nsim=1000,nrep=1, labflag=1, label=label, icovest=3)#significance between neural classical, "pvalnorm":0.0001140764
                                        #last determination of significance of clusters after consensus clustering but before centroid determination
Proneural<-MicroRNASamplefilter(MicroRNA,Proneural)
Neural<-MicroRNASamplefilter(MicroRNA,Neural)
Classical<-MicroRNASamplefilter(MicroRNA,Classical)
Mesenchymal<-MicroRNASamplefilter(MicroRNA,Mesenchymal)#making sure centroid determined subtypes samples all have microRNA expression data
Proneuralgenes <- subtypedfmRNA(Proneural,genes2)
Neuralgenes <- subtypedfmRNA(Neural,genes2)
Classicalgenes <- subtypedfmRNA(Classical,genes2)
Mesenchymalgenes <- subtypedfmRNA(Mesenchymal,genes2)#creating data frames for gene expression with samples of each subtype
#determine significance of clusters after, consensus clustering, centroid determination and then making sure they have matching microRNA expression data
invertP <- t(Proneuralgenes)
invertN <- t(Neuralgenes)
invertM <- t(Mesenchymalgenes)
invertC <- t(Classicalgenes)
sigclustCM <- rbind(invertC,invertM)
sigclustCN <- rbind(invertC,invertN)
sigclustCP <- rbind(invertC,invertP)
sigclustMN <- rbind(invertM,invertN)
sigclustMP <- rbind(invertM,invertP)
sigclustNP <- rbind(invertN,invertP)
label <- c(rep(1,length(Classical)),rep(2,length(Mesenchymal)))
pvalueCM<-sigclust(sigclustCM, nsim=1000, nrep=1, labflag=1, label=label, icovest=3)#"pvalnorm":5.311919e-15
label <- c(rep(1,length(Classical)),rep(2,length(Neural)))
pvalueCN<-sigclust(sigclustCN, nsim=1000, nrep=1, labflag=1, label=label, icovest=3)#"pvalnorm":5.311919e-15
label <- c(rep(1,length(Classical)),rep(2,length(Proneural)))
pvalueCP<-sigclust(sigclustCP, nsim=1000, nrep=1, labflag=1, label=label, icovest=3)#"pvalnorm":5.311919e-15
label <- c(rep(1,length(Mesenchymal)),rep(2,length(Neural)))
pvalueMN<-sigclust(sigclustMN, nsim=1000, nrep=1, labflag=1, label=label, icovest=3)#pvalnorm":2.655806e-10
label <- c(rep(1,length(Mesenchymal)),rep(2,length(Proneural)))
pvalueMP<-sigclust(sigclustMP, nsim=1000, nrep=1, labflag=1, label=label, icovest=3)#"pvalnorm":5.843272e-24
label <- c(rep(1,length(Neural)),rep(2,length(Proneural)))
pvalueNP<-sigclust(sigclustNP, nsim=1000, nrep=1, labflag=1, label=label, icovest=3)#"pvalnorm":0.0200992
ProneuralmiRNA <- subtypedfmiRNA(humanmiRNA,MicroRNA,Proneural)#humanmiRNA is vector of human microRNA, MicroRNA is microRNA expression data and Proneural is samples in proneural subtypes. ProneuralmiRNA is data frame where columns are samples and rows are genes. entries are expression
ClassicalmiRNA <- subtypedfmiRNA(humanmiRNA,MicroRNA,Classical)
NeuralmiRNA <- subtypedfmiRNA(humanmiRNA,MicroRNA,Neural)
MesenchymalmiRNA <- subtypedfmiRNA(humanmiRNA,MicroRNA,Mesenchymal)
rm(Centroids, Centroids2, Centroids3, Correlation, centroidgeneset, genes, genes2, genes3, invert1, invert2, invert3, invert4, invertC, invertM, invertN, invertP, sigclustCM,sigclustCN,sigclustCP,sigclustMN,sigclustMP,sigclustNP,sigclustmatrix1,sigclustmatrix2,sigclustmatrix3,sigclustmatrix4,sigclustmatrix5,sigclustmatrix6,subtype1genes,subtype2genes,subtype3genes,subtype4genes, Class, Classical, Features,MicroRNA,MicroRNAid,MicroRNAid2, Neural, Proneural,clusterindex,clusterindex2,e,hsaindex, i, label, namesindex,pvalue12,pvalue13,pvalue14,pvalue23,pvalue24,pvalue34,pvalueCM,pvalueCN,pvalueCP,pvalueMN,pvalueMP,pvalueNP,results,subtype1,subtype2,subtype3, subtype4,subtypesil,txtfiles, Mesenchymal)#remove unnecessary variables, data and values
library(multiMiR)
genenames <- rownames(Mesenchymalgenes)
mirnanames <- rownames(MesenchymalmiRNA)
ExpMirTar <- targetlist(mirnanames,1)
ExpMirTar <- miRtargetbasefilter(mirnanames, ExpMirTar, genenames)
ExpMirTar <- unique(ExpMirTar)#523
library(RMiR)
miRNA <- cbind(mirnanames,MesenchymalmiRNA[,1])
genes <- cbind(genesnames,Mesenchymalgenes[,1])
colnames(miRNA) <- "expr"
colnames(genes) <- "expr"
genes <- data.frame(genes)
miRNA <- data.frame(miRNA)
ExpMirTar2 <- RmiR(mirna=miRNA,genes=genes,annotation="hgug4112a.db",dbname="tarbase", id="alias")
ExpMirTar2 <- ExpMirTar2[,-3]
ExpMirTar2 <- ExpMirTar2[,-4:-5]
ExpMirTar2 <- ExpMirTar2[,-1]
miRecords <- read.csv("miRecords_version4.csv",sep=",", dec=".", fill=TRUE)
miRecords <- miRecords[miRecords$miRNA_species == "Homo sapiens",]
miRecords <- data.frame(miRecords$miRNA_mature_ID,miRecords$Target.gene_name)
mirTarBase <- read.csv("hsa_MTI.csv",sep=",", dec=".")
mirTarBase <- data.frame(mirTarBase$miRNA, mirTarBase$Target.Gene)
colnames(ExpMirTar2) <- c("miRNA","Targetgene")
colnames(mirTarBase) <- c("miRNA","Targetgene")
colnames(miRecords) <- c("miRNA","Targetgene")
ExpMirTar2 <- rbind(ExpMirTar2,mirTarBase,miRecords)
ExpMirTar2 <- unique(ExpMirTar2)
rownames(ExpMirTar2) <- NULL
ExpMirTar2 <- miRtargetbasefilter(mirnanames, ExpMirTar2, genenames)
ExpMirTar2<-data.frame(lapply(ExpMirTar2,function(x)factor(unlist(x))))
ExpMirTar3 <- rbind(ExpMirTar, ExpMirTar2)#ExpMirTar had 523 interactions and ExpMirTar2 had 2133. ExpMirTar is from multilMiR package. ExpMirTar2 created from accessing miRecords and miRTarBase directly while obtaining data from tarbase via RmiR
ExpMirTar3 <- unique(ExpMirTar3)#ExpMirTar3 has 2249 interactions
PredMirTar1 <- targetlist(genenames[1:50],2)
PredMirTar2 <- targetlist(genenames[51:100],2)
PredMirTar3 <- targetlist(genenames[101:150],2)
PredMirTar4 <- targetlist(genenames[151:200],2)
PredMirTar5 <- targetlist(genenames[201:250],2)
PredMirTar6 <- targetlist(genenames[251:300],2)
PredMirTar7 <- targetlist(genenames[301:350],2)
PredMirTar8 <- targetlist(genenames[351:400],2)
PredMirTar9 <- targetlist(genenames[401:450],2)
PredMirTar10 <- targetlist(genenames[451:500],2)
PredMirTar11 <- targetlist(genenames[501:550],2)
PredMirTar12 <- targetlist(genenames[551:600],2)
PredMirTar13 <- targetlist(genenames[601:650],2)
PredMirTar14 <- targetlist(genenames[651:700],2)
PredMirTar15 <- targetlist(genenames[701:750],2)
PredMirTar16 <- targetlist(genenames[751:800],2)
PredMirTar17 <- targetlist(genenames[801:850],2)
PredMirTar18 <- targetlist(genenames[851:900],2)
PredMirTar19 <- targetlist(genenames[901:950],2)
PredMirTar20 <- targetlist(genenames[951:1000],2)
PredMirTar21 <- targetlist(genenames[1001:1050],2)
PredMirTar22 <- targetlist(genenames[1051:1100],2)
PredMirTar23 <- targetlist(genenames[1101:1150],2)
PredMirTar24 <- targetlist(genenames[1151:1200],2)
PredMirTar25 <- targetlist(genenames[1201:1250],2)
PredMirTar26 <- targetlist(genenames[1251:1300],2)
PredMirTar27 <- targetlist(genenames[1301:1350],2)
PredMirTar28 <- targetlist(genenames[1351:1400],2)
PredMirTar29 <- targetlist(genenames[1401:1450],2)
PredMirTar30 <- targetlist(genenames[1451:1500],2)
PredMirTar31 <- targetlist(genenames[1501:1550],2)
PredMirTar32 <- targetlist(genenames[1551:1600],2)
PredMirTar33 <- targetlist(genenames[1601:1650],2)
PredMirTar34 <- targetlist(genenames[1651:1700],2)
PredMirTar35 <- targetlist(genenames[1701:1740],2)
PredMirTar <- rbind(PredMirTar1, PredMirTar2, PredMirTar3, PredMirTar4, PredMirTar5, PredMirTar6, PredMirTar7, PredMirTar8, PredMirTar9, PredMirTar10, PredMirTar11, PredMirTar12, PredMirTar13, PredMirTar14, PredMirTar15, PredMirTar16, PredMirTar17, PredMirTar18, PredMirTar19, PredMirTar20, PredMirTar21, PredMirTar22, PredMirTar23, PredMirTar24, PredMirTar25, PredMirTar26, PredMirTar27,PredMirTar28, PredMirTar29, PredMirTar30, PredMirTar31, PredMirTar32, PredMirTar33, PredMirTar34, PredMirTar35)#821,116 interactions, need to remove redundancy and unnecessary interactions
write.table(PredMirTar, "rawPredMirTar.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
PredMirTar <- miRtargetbasefilter(mirnanames, PredMirTar, genenames)#remove unnecessary interactions, only 188 mirna's between the mirna from our data in common with the mirna from the databases
PredMirTar <- unique(PredMirTar)#remove redundancy
rownames(PredMirTar) <- NULL
rm(PredMirTar1, PredMirTar2, PredMirTar3, PredMirTar4, PredMirTar5, PredMirTar6, PredMirTar7, PredMirTar8, PredMirTar9, PredMirTar10, PredMirTar11, PredMirTar12, PredMirTar13, PredMirTar14, PredMirTar15, PredMirTar16, PredMirTar17, PredMirTar18, PredMirTar19, PredMirTar20, PredMirTar21, PredMirTar22, PredMirTar23, PredMirTar24, PredMirTar25, PredMirTar26, PredMirTar27,PredMirTar28, PredMirTar29, PredMirTar30, PredMirTar31, PredMirTar32, PredMirTar33, PredMirTar34, PredMirTar35)
ExpNetwork <- Networkformation2(MesenchymalmiRNA, ExpMirTar3)
PredNetwork <- Networkformation2(MesenchymalmiRNA, PredMirTar)
write.table(PredMirTar, "processedPredMirTar.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(ExpMirTar, "rawExpMirTar.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(ExpMirTar, "processedExpMirTar.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(ExpNetwork, "ExpNetwork.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(PredNetwork, "PredNetwork.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(ExpMirTar2, "ExpMirTar2.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(ExpMirTar3, "ExpMirTar3.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
rm(ExpMirTar, ExpMirTar2,ExpMirTar3,PredMirTar,miRecords,mirTarBase)#remove unnecessary variables
rm(MicroRNASamplefilter,MicroRNAfilter,miRtargetbasefilter,positivecounter,positivesilhouette,subtypecounter,subtypedfmRNA,subtypedfmiRNA,subtypesorter,targetlist)#remove unnecessary functions
library(permute)
library(infotheo)

CTRL <- how(nperm = 1000) #necessary parameters for permutations
ClassicalcondinfoExpNet <- condinfonetwork2(ExpNetwork,Classicalgenes, ClassicalmiRNA)
ClassicalcondinfoExpNet <- condinfonetwork3(ClassicalcondinfoExpNet,Classicalgenes,ClassicalmiRNA)
ClassicalcondinfoExpNet <- condinfonetwork4(ClassicalcondinfoExpNet,Classicalgenes,ClassicalmiRNA)
ClassicalcondinfoExpNet <- condinfonetwork5(ClassicalcondinfoExpNet,Classicalgenes,ClassicalmiRNA)
ClassicalcondinfoExpNet <- condinfonetwork6(ClassicalcondinfoExpNet,Classicalgenes,ClassicalmiRNA)
ClassicalcondinfoExpNet <- data.frame(lapply(ClassicalcondinfoExpNet,function(x)factor(unlist(x))))
write.table(ClassicalcondinfoExpNet, "ClassicalcondinfoExpNet.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
ProneuralcondinfoExpNet <- condinfonetwork2(ExpNetwork,Proneuralgenes, ProneuralmiRNA)
ProneuralcondinfoExpNet <- condinfonetwork3(ProneuralcondinfoExpNet,Proneuralgenes,ProneuralmiRNA)
ProneuralcondinfoExpNet <- condinfonetwork4(ProneuralcondinfoExpNet,Proneuralgenes,ProneuralmiRNA)
ProneuralcondinfoExpNet <- condinfonetwork5(ProneuralcondinfoExpNet,Proneuralgenes,ProneuralmiRNA)
ProneuralcondinfoExpNet <- condinfonetwork6(ProneuralcondinfoExpNet,Proneuralgenes,ProneuralmiRNA)
ProneuralcondinfoExpNet <- data.frame(lapply(ProneuralcondinfoExpNet,function(x)factor(unlist(x))))
write.table(ProneuralcondinfoExpNet, "ProneuralcondinfoExpNet.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
NeuralcondinfoExpNet <- condinfonetwork2(ExpNetwork,Neuralgenes, NeuralmiRNA)
NeuralcondinfoExpNet <- condinfonetwork3(NeuralcondinfoExpNet,Neuralgenes,NeuralmiRNA)
NeuralcondinfoExpNet <- condinfonetwork4(NeuralcondinfoExpNet,Neuralgenes,NeuralmiRNA)
NeuralcondinfoExpNet <- condinfonetwork5(NeuralcondinfoExpNet,Neuralgenes,NeuralmiRNA)
NeuralcondinfoExpNet <- condinfonetwork6(NeuralcondinfoExpNet,Neuralgenes,NeuralmiRNA)
NeuralcondinfoExpNet <- data.frame(lapply(NeuralcondinfoExpNet,function(x)factor(unlist(x))))
write.table(NeuralcondinfoExpNet, "NeuralcondinfoExpNet.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
MesenchymalcondinfoExpNet <- condinfonetwork2(ExpNetwork,Mesenchymalgenes, MesenchymalmiRNA)
MesenchymalcondinfoExpNet <- condinfonetwork3(MesenchymalcondinfoExpNet,Mesenchymalgenes, MesenchymalmiRNA)
MesenchymalcondinfoExpNet <- condinfonetwork4(MesenchymalcondinfoExpNet,Mesenchymalgenes, MesenchymalmiRNA)
MesenchymalcondinfoExpNet <- condinfonetwork5(MesenchymalcondinfoExpNet,Mesenchymalgenes, MesenchymalmiRNA)
MesenchymalcondinfoExpNet <- condinfonetwork6(MesenchymalcondinfoExpNet,Mesenchymalgenes, MesenchymalmiRNA)
MesenchymalcondinfoExpNet <- data.frame(lapply(MesenchymalcondinfoExpNet,function(x)factor(unlist(x))))
write.table(MesenchymalcondinfoExpNet, "MesenchymalcondinfoExpNet.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
ExpNetwork <- organizenetwork(ExpNetwork)#something is wrong, there are not as many ExpNetwork interactions as there were before

PredNetwork <- organizenetwork(PredNetwork)#498345
organizedPredNetwork <- data.frame(lapply(PredNetwork,function(x)factor(unlist(x))))#definitely organized because tried to write PredNetwork to file and got EncodeElement error, meaning it was in some type of list format
write.table(organizedPredNetwork, "organizedPredNetwork.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
ClassicalExpp <- pvalueExp(ExpNetwork, ClassicalcondinfoExpNet)
ProneuralExpp <- pvalueExp(ExpNetwork, ProneuralcondinfoExpNet)
NeuralExpp <- pvalueExp(ExpNetwork, NeuralcondinfoExpNet)
MesenchymalExpp <- pvalueExp(ExpNetwork, MesenchymalcondinfoExpNet)
ProneuralExpp <- data.frame(lapply(ProneuralExpp,function(x)factor(unlist(x))))
NeuralExpp <- data.frame(lapply(NeuralExpp,function(x)factor(unlist(x))))
ClassicalExpp <- data.frame(lapply(ClassicalExpp,function(x)factor(unlist(x))))
MesenchymalExpp <- data.frame(lapply(MesenchymalExpp,function(x)factor(unlist(x))))
write.table(ClassicalExpp, "ClassicalExpp.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(ProneuralExpp, "ProneuralExpp.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(NeuralExpp, "NeuralExpp.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(MesenchymalExpp, "MesenchymalExpp.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
ClassicalExpPvalue <- Pvalue(ClassicalExpp)
ProneuralExpPvalue <- Pvalue(ProneuralExpp)
NeuralExpPvalue <- Pvalue(NeuralExpp)
MesenchymalExpPvalue <- Pvalue(MesenchymalExpp)
ProneuralExpPvalue <- data.frame(lapply(ProneuralExpPvalue,function(x)factor(unlist(x))))
NeuralExpPvalue <- data.frame(lapply(NeuralExpPvalue,function(x)factor(unlist(x))))
ClassicalExpPvalue <- data.frame(lapply(ClassicalExpPvalue,function(x)factor(unlist(x))))
MesenchymalExpPvalue <- data.frame(lapply(MesenchymalExpPvalue,function(x)factor(unlist(x))))
write.table(ClassicalExpPvalue, "ClassicalExpPvalue.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(ProneuralExpPvalue, "ProneuralExpPvalue.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(NeuralExpPvalue, "NeuralcondinfoExpPvalue.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(MesenchymalExpPvalue, "MesenchymalExpPvalue.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
ClassicalExpNet1 <- Pvaluefilter(ClassicalExpPvalue)
ProneuralExpNet1 <- Pvaluefilter(ProneuralExpPvalue)
NeuralExpNet1 <- Pvaluefilter(NeuralExpPvalue)
MesenchymalExpNet1 <- Pvaluefilter(MesenchymalExpPvalue)
write.table(ClassicalExpNet1, "ClassicalExpNet1.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(ProneuralExpNet1, "ProneuralExpNet1.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(NeuralExpNet1, "NeuralExpNet1.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(MesenchymalExpNet1, "MesenchymalExpNet1.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
ClassicalExpNet2 <- Pvaluefilter2(ClassicalExpPvalue)
ProneuralExpNet2 <- Pvaluefilter2(ProneuralExpPvalue)
NeuralExpNet2 <- Pvaluefilter2(NeuralExpPvalue)
MesenchymalExpNet2 <- Pvaluefilter2(MesenchymalExpPvalue)
write.table(ClassicalExpNet2, "ClassicalExpNet2.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(ProneuralExpNet2, "ProneuralExpNet2.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(NeuralExpNet2, "NeuralExpNet2.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)
write.table(MesenchymalExpNet2, "MesenchymalExpNet2.txt", quote=FALSE, sep="\t", dec=".", col.names=TRUE, row.names=TRUE)

rm(ClassicalExpNet1, ClassicalExpNet2, ClassicalcondinfoExpNet,ExpMirTar,ExpNetwork, MesenchymalExpNet1, MesenchymalExpNet2,MesenchymalcondinfoExpNet,NeuralExpNet1, NeuralExpNet2, NeuralcondinfoExpNet, ProneuralExpNet1, ProneuralExpNet2, ProneuralcondinfoExpNet,genenames, humanmiRNA, mirnanames,x, Networkformation2, PredNetwork,condinfonetwork,condinfonetwork2, condinfonetwork2, condinfonetwork3, condinfonetwork3, condinfonetwork4, condinfonetwork5, condinfonetwork6, organizenetwork, pvalueExp)
Setsofdata <- function(Network,n){#function creates a list where each member of the list is a data frame with n rows
    data <- list()
    x <- ceiling((length(Network[,1]))/n)
    a <- (((length(Network[,1]))/n)-floor((length(Network[,1]))/n))*n
    for(i in 1:x){
        y <- (i*n)
        z <- (((i-1)*n)+1)
        if (i == x){
            data[[i]] <- Network[z:(z+a),1:3]
        }
         data[[i]] <- Network[z:y,1:3]
    }
    return(data)
}
WriteListData <- function(listofdata){ #for output of Setsofdata
    for (i in 1:length(listofdata)){
        write.table(listofdata[[i]], sprintf("data%d.txt",i), quote=FALSE, sep="\t",dec=".", col.names=TRUE, row.names=TRUE)
    }
}
createpermuteset <- function(a){#a is number of permutations, in this case it would be 1000
    x <- shuffle(n=a, control=CTRL)
    for (i in 1:999){
        x <- rbind(x, permute(n=a,control=CTRL))
    }
    return(x)
}
condinfonetworkbyrow <- function(textfile, genes, miRNA,permuteset){#evaluating conditional mutual information by row, a better way to evaluate. Read in each data file, compute conditional mutual information, write it to file.
    network <- read.table(textfile, header=TRUE, sep="\t", dec=".", quote="")
    for(i in 1:1001){
        if(i == 1){
            network["observeddeltaMI"] <- NA
        }else{
            network[sprintf("expecteddeltaMI%d",i)] <- NA
        }
    }
    a <- rownames(miRNA)
    b <- rownames(genes)
    for (i in 1:length(network[,1])){
        c <- which(a == network[i,1])
        d <- which(b == network[i,2])
        e <- which(a == network[i,3])
        f <- miRNA[c,]
        g <- genes[d,]
        h <- miRNA[e,]
        for (j in 1:1001){
            if (j == 1){
                discretef<-discretize(f, disc="equalwidth" )
                discreteg <- discretize(g, disc="equalwidth")
                discreteh <- discretize(h, disc="equalwidth")
                network$observeddeltaMI[i] <- condinformation(discretef,discreteh,discreteg, method="emp")
            }else{
                z <- j-1
                y <- permuteset[z,]
                y <- as.numeric(y)
               # g <- t(data.matrix(g))
                g <- g[y]
                discreteg <- discretize(g, disc="equalwidth")
                network[i,(j+3)] <- condinformation(discretef,discreteh,discreteg, method="emp")
        }
    }
    }
    write.table(network, textfile, quote=FALSE,sep="\t",dec=".", col.names=TRUE, row.names=TRUE)
}
check <- function(textfiles){
    for (i in 1:length(textfiles)){
        network <- read.table(textfiles[i], header=TRUE, sep="\t", dec=".", quote="")
        if length(
parallelcondinfo <- function(cls, textfiles, genes, miRNA, permuteset){
    clusterApplyLB(cl=cls,x=textfiles, fun=condinfonetworkbyrow, genes, miRNA, permuteset)
}

permuteset <- createpermuteset(length(colnames(ClassicalmiRNA)))
write.table(permuteset, "Classicalpermuteset.txt", quote=FALSE, sep="\t", row.names=FALSE)
data <- Setsofdata(organizedPredNetwork, 250)
setwd("C:/Users/Nadia/Documents/glioblastoma/ClassicalPredNetwork")
WriteListData(data)
txtfiles=list.files(pattern = 'data*')#make sure all textfiles are accessible to console by making sure all there names are variables.

library(parallel)
library(infotheo)#packaages necessary for process
clust<-makePSOCKcluster(rep(c("localhost"), each=8))#create a cluster for parallel computing to compute conditional mutual information for predictive networks
clusterEvalQ(clust, library(stats))#sending all necessary packages and cluster environment
clusterEvalQ(clust, library(infotheo))
clusterEvalQ(clust, library(parallel))
varlist<- c("txtfiles","permuteset","Classicalgenes","ClassicalmiRNA","condinfonetworkbyrow")#create variable list for all necessary components.
clusterExport(clust,varlist,envir =.GlobalEnv)
parallelcondinfo(clust, txtfiles, Classicalgenes, ClassicalmiRNA, permuteset)
clusterApplyLB(clust, x=txtfiles[1:8], fun=condinfonetworkbyrow, Classicalgenes, ClassicalmiRNA, permuteset)
stopCluster(clust)

# install libraries instruction
# if(!requireNamespace('BiocManager', quietly=TRUE))
#   install.packages('BiocManager')
# BiocManager::install('GoSemSim')
# BiocManager::install('org.Sc.sgd.db')

library(GOSemSim)
library('org.Sc.sgd.db')
hs <- org.Sc.sgd.db
# keytypes(hs)
hsGO <- godata('org.Sc.sgd.db', ont='MF', keytype='GENENAME')

for (method in list('L3uvJoin', 'uvSpec_xySpec_xyContrib_uvJoin', 'xyContrib_dualCN_uvJoin')) {
    for (dataset in list('bioGRID')) {
        for (i in 0:9) {
            data <- read.csv(paste('./GOSemSimPrepData/', paste(paste(method, 'tenTrial', dataset, i, sep='_'), 'csv', sep='.'), sep=''), sep='\t')
            scoreV <- vector(mode='numeric', length=nrow(data))
            for (j in 1:nrow(data)) {
                res <- geneSim(data$nodeA[j], data$nodeB[j], semData=hsGO, measure='Wang', combine='BMA')
                if (class(res) == 'list') {
                    scoreV[j] <- res$geneSim
                } else {
                    scoreV[j] <- 0
                }
                # print(scoreV[j])
            }
            print("completed")
            # x <- data.frame("score"=scoreV)
            data<- cbind(data, score=scoreV)
            write.table(data
            , paste(file='./GOSemSimResData/', paste(paste(method, 'tenTrial', dataset, i, sep='_'), 'csv', sep='.'), sep='')
            , row.names=FALSE
            , sep="\t"
            , quote=FALSE)
        }
    }
}
#res <- geneSim('MDE1', 'PHO87', semData=hsGO, measure='Wang', combine='BWA')
#res$geneSim # returns a numeric value


# ppiList <- vector(mode='list', length=nrow(4*2*10))
# methods <- list('commonNeighbor', 'L3Normalizing', 'uvSpec_xySpec_xyContrib', 'xyContrib_dualCN')
# datasets <- list('bioGRID', 'STRING')
# for (methodI in 1:4) {
#     for (datasetI in 1:2) {
#         for (i in 0:9) {
#             dataset <- datasets[datasetI]
#             method <- methods[methodI]
#             data <- read.csv(paste('./csvResultData/', paste(paste(method, 'tenTrial', dataset, i, sep='_'), 'csv', sep='.'), sep=''), sep='\t')
#             ppi <- vector(mode='list', length=nrow(data))
#             for (j in 1:nrow(data)) {
#                 ppi[j] <- vector(mode='list', length=2)
#             }
#             for (j in 1:nrow(data)) {
#                 ppi[[j]][1] = data$nodeA[j]
#                 ppi[[j]][2] = data$nodeB[j]
#             }
        
#         }
#     }
# }
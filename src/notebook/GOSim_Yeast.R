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

for (method in list('commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'CRA', 'Sim', 'CH2_L3')) {
    for (dataset in list('IntAct_spoke', 'MINT', 'bioGRID', 'STRING')) {
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
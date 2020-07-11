# install libraries instruction
# if(!requireNamespace('BiocManager', quietly=TRUE))
#   install.packages('BiocManager')
# BiocManager::install('GoSemSim')
# BiocManager::install('org.Sc.sgd.db')

library(GOSemSim)
library('org.Hs.eg.db')
hs <- org.Hs.eg.db
# keytypes(hs)
hsGO <- godata('org.Hs.eg.db', ont='MF', keytype='UNIPROT')

for (method in list('commonNeighbor', 'L3uvJoin', 'xyContrib_dualCN_uvJoin', 'Sim', 'CH2_L3', 'CRA')) {
    for (dataset in list('HuRI', 'MINT_homo', 'bioGRID_homo', 'STRING_homo')) {
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
import pandas as pd
import json

uniprot_path = "../../UniProt/uniprot-taxonomy_9606.tab"
up_df = pd.read_csv(uniprot_path, sep="\t")
print(up_df.head())
geneCol = up_df.keys()[1]
up_gene = set(list(up_df[geneCol]))

df = pd.read_csv('essGene.csv')
print(df.head())
essGene = set(list(df['geneName']))

interEssGene = up_gene&essGene
print(len(up_gene), len(essGene), len(interEssGene))

with open("../essGene_homo.json", "w") as f:
    f.write(json.dumps(list(interEssGene)))
from bioservices.kegg import KEGG
import pandas as pd


class KEGGFinder:
    def __init__(self):
        self.kegg = KEGG()
        self.result

    def find(self, pathway):
        raw = self.kegg.get(pathway)
        self.result = self.kegg.parse(raw)
        return self.result

    def get_gene_dataframe(self):
        genes = self.result['GENE']
        symbols = []
        info = []
        for i in genes.values():
            symbols.append(i.split(';')[0])
            info.append(i.split(';')[1])
        tmp = {'id': genes.keys(), 'symbol': symbols, 'info': info}
        return pd.DataFrame(tmp)

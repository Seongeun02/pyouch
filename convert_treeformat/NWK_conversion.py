from ete3 import Tree
import pandas as pd

class NWK_conversion:
    '''
    Class for convert the Newick tree file into appropriate format(distance matrix or tree data for ou model)

    Methods :
    1) distance_matrix : convert tree file into distance matrix
    2) tree_data : convert tree file into tree data for ou model
    '''

    def __init__(self,
                 df,
                 species):
        '''
        Initialize class instance

        Arguments:
        df(str) : either file directory or real Newick tree string
        species(list) : list of species, check the name!! They should match with the data in df
        '''

        self.df = df
        self.species = species
        # convert into tree object & prune it 
        self.t = Tree(self.df, format=1)
        self.t.prune(self.species)
        # extract nodes in the tree
        self.node = self.t.search_nodes()

    
    def distance_matrix(self):
        an = {}
        for i in self.species:
            dist=[]
            for j in self.species:
                if i==j:
                    dist.append(0)
                else:
                    path = self.t.get_distance(i,j)
                    dist.append(path)

            an[i] = dist
        an = pd.DataFrame(an)
        an.insert(0, "species", self.species)
        return an
    
    def tree_data(self):
        # rename internal nodes
        a = 0
        for i in self.node:
            if not i.is_leaf():
                i.name = a
                a += 1

        # get dataframe for ancestor, time, and species 
        ancestor = []
        time = []
        sp = []
        for i in range(len(self.node)):
            if i==0:
                ancestor.append(None)
                time.append(0)
                sp.append(self.node[i].name)
            else:
                ancestor.append(self.node[i].up.name)
                time.append(self.node[0].get_distance(self.node[i]))
                sp.append(self.node[i].name)

        df = pd.DataFrame({"species":sp, "ancestor":ancestor, "time":time})

        # rearrange into right format
        idx = []
        num = len(self.species)
        df["species"] = df["species"].astype(str)
        for i in df["species"]:
            if i.isdigit():
                idx.append(int(i)+1)
            else:
                idx.append(num)
                num += 1 
        df.insert(0,"node",idx)
        df["node"] = df["node"].astype(int)
        df.sort_values(by="node", inplace=True)

        df["ancestor"] += 1
        df.index=list(range(len(self.species)*2-1))
        df.iloc[:len(self.species)-1,1] = None

        return df
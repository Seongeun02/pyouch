import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class time_simulation:
    def __init__(self, 
                 tree,
                 rate = 0.001,
                 regime = "univariate",
                 num_ou = 700,
                 num_bm = 300,
                 params = "random",
                 max_a = 10,
                 max_s = 10,
                 max_th = 100,
                 ou_var_threshold = 1):
        '''
        Arguments 
        tree(DataFrame) : contains nodes, species, ancestors, and time column
        rate(float) : size of time grid(=dt)
        regime(list or "univariate") : optimum regime, length should be same as the number of nodes
        num_ou(int) : number of nonalpha parameter sets
        num_bm(int) : number of zero-alpha parameter sets
        params(DataFrame or "random") : set of ouch parameters, containing alpha, sigma, and theta 
        max_a(float) : maximum alpha value for random parameter generation
        max_s(float) : maximum sigma value for random parameter generation
        max_th(float) : maximum theta value for random parameter generation
        ou_var_threshold(float) : maximum ou variance(sigma^2 / (2 alpha)) for random parameter generation
        '''

        self.tree = tree

        # delta t
        self.rate = rate
        
        # number of genes under selection(ou) or drift(bm)
        self.num_ou = num_ou
        self.num_bm = num_bm

        # regime
        if regime == "univariate":
            self.regime = [0 for i in range(tree.shape[0])]
        else:
            self.regime = self.set_regime()
        self.regime = np.array(self.regime)
        self.nreg = len(np.unique(self.regime))

        # set boundary value for random parameters
        self.max_a = max_a
        self.max_s = max_s
        self.max_th = max_th
        self.ou_var_threshold = ou_var_threshold

        # generate parameters
        if params == "random" : 
            self.params = self.rand_param()
        else:
            self.params = params


    def set_regime(self):
        '''
        generate regime list with unique numerical values
        '''
        reg = np.unique(self.regime)
        numerical_regime = [reg.index(item) for item in self.regime]

        return numerical_regime
    

    def rand_param(self):
        ## generate random parameters
        params = {}
        i = 1

        while len(params) < self.num_ou:
            #alpha = np.random.uniform(0, self.max_a)
            alpha = np.random.uniform(50, 100)
            theta = np.random.uniform(0, self.max_th, self.nreg)
            sigma = np.random.uniform(0, self.max_s)

            if sigma**2 <= alpha * 2 * self.ou_var_threshold:
                ## keep the ou variance lower than ou_var_threshold
                params[i] = [alpha, theta, sigma]
                i += 1

        for j in range(self.num_bm):
            theta = np.random.uniform(0, self.max_th, self.nreg)
            sigma = np.random.uniform(0, self.max_s)    

            params[i+j] = [0,theta, sigma]

        params = pd.DataFrame(params)
        params = params.T
        params.columns = ["alpha", "theta", "sigma"]

        return params


    def sim(self, param):
        theta = param["theta"]
        sigma = param["sigma"]
        alpha = param["alpha"]

        X = {1:theta[0]}
        for i in range(1,tree.shape[0]):
            ti = self.tree["time"][self.tree["ancestor"][i]-1]
            tf = self.tree["time"][i] 
            x = X[self.tree["ancestor"][i]]

            for j,t in enumerate(np.arange(ti, tf, self.rate)):
                x += alpha * (theta[self.regime[i]] - x) + sigma * np.random.normal(0, self.rate, 1)[0]

            X[i + 1]= x

        return X
    

    def run(self):
        data = {}
        nterm = self.tree[self.tree["species"].notna()].shape[0]

        for i in self.params.index.tolist():
            if i % 500 == 0:
                print(i)
                print(X)
            X = self.sim(self.params.loc[i,:])
            data[i] = X

        data = pd.DataFrame(data)
        data = data.loc[nterm:,:].T
        data.columns = self.tree["species"][nterm-1:]

        return self.params, data
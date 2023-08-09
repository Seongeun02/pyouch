import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class OuchSimulation:
    '''
    generate Ouch Simulated data with randomly generated parameters (sigma, alpha, theta)

    methods 
    set_tree : return epoch and lineage of each terminal nodes with tree data
    regime_spec : return beta matrix based on the regime data
    calculate_weight : return weight matrix(W) 
    calculate_exp : return expectation value of each terminal nodes (W * theta)
    diverged_time : return diverged time between each terminal nodes (S)
    calculate_var : return covariance matrix of each terminal nodes
    add_noise : add noise to the expectation value (sampling from the multivariate normal distribution)
    rand_param : generate random parameter sets
    run : generate simulated dataset using all the functions above
    '''

    def __init__(self, 
                 tree,
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
        self.num_ou = num_ou
        self.num_bm = num_bm

        # regime
        if regime == "univariate":
            self.regime = ["ns" for i in range(tree.shape[0])]
        else:
            self.regime = regime
        self.regime = np.array(self.regime)
        self.nreg = len(np.unique(self.regime))

        # set tree 
        self.epoch, self.lineage = self.set_tree()

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


    def set_tree(self):
        '''
        generate epoch and lineage data from tree
        epoch(numpy array) : branch timeline of each lineage
        lineage(dict) : nodes in each lineage
        '''
        terminal = self.tree[self.tree["species"].notna()].index.tolist()   
        lineage = {}
        epoch = np.zeros([len(terminal),self.tree.shape[0]])

        for i in range(len(terminal)):
            lineage[i] = [terminal[i]+1]
            an = int(self.tree["ancestor"][terminal[i]])
            epoch[i][terminal[i]] = self.tree["time"][terminal[i]]
            while an != 1:
                lineage[i].append(an)
                epoch[i][an-1] = self.tree["time"][an-1]
                an = int(self.tree["ancestor"][an-1])
            lineage[i].append(1)

        return epoch, lineage


    def regime_spec(self):
        '''
        generate beta matrix
        1st axis : terminal nodes
        2nd axis : each lineage in terminal nodes
        '''
        nterm = len(self.lineage)
        reg = np.unique(self.regime)

        beta = [[] for _ in range(nterm)]

        for i in range(nterm):
            p = np.array(self.lineage[i]) - 1
            beta[i] = np.zeros((len(p), self.nreg), dtype=int)
            for r in range(self.nreg):
                beta[i][:,r] = np.where(self.regime[p] == reg[r], 1, 0)

        return beta


    def calculate_weight(self, params):
        '''
        calculate weight matrix, W
        '''
        alpha = params["alpha"]
        T = np.max(self.tree["time"])      # depth of the tree

        W = np.zeros([len(self.lineage),len(np.unique(self.regime))+1])
        beta = self.regime_spec()

        for i in range(len(self.lineage)):
            W[i][0] = 1
            for k in range(len(self.lineage[i])-1):
                ti = self.tree["time"][self.lineage[i][k+1]-1]
                tf = self.tree["time"][self.lineage[i][k]-1]
                W[i][beta[i][k][0]] += (np.exp(alpha*tf) - np.exp(alpha*ti)) 

        W = W * np.exp(-alpha*T)

        return W

    
    def calculate_exp(self, params, W):
        ## calculate expectation values for given parameters, W*theta
        theta = np.insert(params["theta"], 0, params["theta"][0])
        exp = np.dot(W,theta).reshape([np.shape(W)[0], ])

        return exp


    def diverged_time(self):
        '''
        calculate the diverged time matrix S
        size : number of terminal nodes(species) x number of terminal nodes(species)
        '''
        n = np.shape(self.epoch)[0]
        time_mat = np.zeros([n,n])
        for i in range(n):
            nonzero_i, = np.nonzero(self.epoch[i])
            for j in range(n):
                nonzero_j, = np.nonzero(self.epoch[j])
                common_values = np.intersect1d(nonzero_i, nonzero_j) 
                if len(common_values) != 0:
                    max_common_value = np.max(common_values)
                    time_mat[i][j] = self.epoch[i][max_common_value]
        return time_mat
    

    def calculate_var(self, params, S):
        '''
        calculate covariance-variance matrix, V
        size : number of terminal nodes(species) x number of terminal nodes(species)
        '''
        alpha = params["alpha"]
        sigma = params["sigma"]
        T = np.max(self.tree["time"])
        n = np.shape(S)[0]

        if alpha == 0:
            ## In brownian model, Vij = sigma**2 * (diverged time)
            V = sigma**2 * S
        else:
            ## ou model
            V = np.zeros([n,n])
            for i in range(n):
                for j in range(n):
                    V[i][j] = np.exp(2*alpha*S[i][j]) - 1
            V *= (sigma**2 / (2*alpha)) * np.exp(-2*alpha*T)

        return V
    

    def add_noise(self, exp, V):
        ## add noise to the expectation values
        cf = np.linalg.cholesky(V)
        rnorm = np.random.normal(size=(np.shape(V)[0], 1))
        return np.dot(cf, rnorm).T + np.array(exp).reshape(1, np.shape(V)[0])
    

    def rand_param(self):
        ## generate random parameters
        params = {}
        i = 1

        while len(params) < self.num_ou:
            alpha = np.random.uniform(0, self.max_a)
            theta = np.random.uniform(0, self.max_th, self.nreg)
            sigma = np.random.uniform(0, self.max_s)

            if sigma**2 <= alpha * 2 * self.ou_var_threshold:
                ## keep the ou variance lower than 1
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
    

    def run(self):
        data = {}
        S = self.diverged_time()
        for i in self.params.index.tolist():
            if i % 500 == 0:
                print(i)
            p = self.params.loc[i,:]
            W = self.calculate_weight(p)
            V = self.calculate_var(p, S)
            exp = self.calculate_exp(p,W)
            data[i] = self.add_noise(exp,V)[0]

        data = pd.DataFrame(data)
        data = data.T
        data.columns = self.tree["species"][np.shape(S)[0]-1:]
        #data = np.log(data+0.01)
        return self.params, data


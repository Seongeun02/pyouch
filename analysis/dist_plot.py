import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import minimize_scalar
from scipy.optimize import curve_fit
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import statsmodels.api as sm

import warnings
warnings.filterwarnings("ignore")

import re


class dist_plot:
    '''
    Class for generating the plot, distance between species vs evolutionary time

    Methods:
    categorize : categorize genes based on their fdr value (stabilized selection for fdr < 0.05)
    residual : calculate orthogonal distance between species using pca
    power_fit : power fit to expression distance versus time (y=a*x^k)
    linear_fit : linear fit to expression distance versus time (y=m*x+n)
    calculate_mean : calculate mean squared distance from reference species. 
               First, mean across genes. Then, calculate mean & std across replicates in species
    plot : scatter plot with error bar & fitted lines
    run : run all the methods above, and get 3 plots, total genes, genes under stabilizing selection, 
          and genes under neutral evolution.
    '''

    def __init__(self,
                 X, 
                 dist_mat,
                 stat,
                 ref = "Human",
                 title = "distance plot"
                 ):
        '''
        Initialize class instance

        Arguments:
        X (dataframe) : RNA expression data, columns(samples) & rows(genes)
        dist_mat (dataframe) : evolutionary time difference from reference species
        stat (dataframe) : outcome of OU model containing gene_name, qvalues, thetas, var, and brownSigmas
        title (str) : title of plot
        '''
        
        self.X = X.copy()
        self.stat = stat
        self.ref = ref
        self.title = title

        # unique species in data
        col = [re.sub(r'\.\d+', '', col) for col in self.X.columns]
        self.species = np.unique(col)

        # extract the evolutionary time 
        self.timedata = dict(zip(dist_mat["species"],dist_mat[self.ref]))

        self.stabilized = []
        self.neutral = []
        self.res = self.residual()


    def categorize(self):
        for i, q in enumerate(self.stat["qvalues"]):
            if q < 0.05:
                self.stabilized.append(self.stat["gene_name"][i])      ## check the index of stat
            elif q > 0.05:
                self.neutral.append(self.stat["gene_name"][i])
        # stabilized selection
        self.stabilized = self.res[self.res.index.isin(self.stabilized)]
        # neutral evolution
        self.neutral = self.res[self.res.index.isin(self.neutral)]


    def residual(self):

        '''
        Calculate the mean squared distance using pca
        For given dataset, calculate the mean of reference species
        Then, perform pca for dataset composed of each sample and reference mean
        store the data projection to 2nd principal component
        '''

        # prepare reference mean
        ref_mean = [col for col in self.X.columns if self.ref in col]
        ref = self.X[ref_mean].mean(axis=1)

        ref_mean = np.log10(ref + 0.01) 
        data = np.log10(self.X + 0.01)
        data = data[(~data.isna().any(axis=1)) & (~ref_mean.isna())]

        # perform pca
        print("performing pca")
        residuals = {}
        for col in data.columns:
            other = data[col]
            other_data = other
            ref_data = ref_mean[other.index]
            #ref_data = ref_mean[(ref_mean > 0) & (other > 0)]
            #other_data = other[(ref_mean > 0) & (other > 0)]
            t = np.stack([ref_data, other_data], axis=1)
            m = PCA(n_components=2).fit(t)
            res = m.transform(t)[:, 1]
            if m.components_[1, 1] < 0:
                res *= -1
            residuals[col] = res

        residuals = pd.DataFrame(residuals, columns=data.columns, index=data.index)

        return residuals


    def f(self, t, a, k):
        # power law fitting
        return a*t**k

    def power_fit(self, mean, t):
        try:
            # curve fit
            popt, pcov = curve_fit(self.f, xdata=t, ydata=mean, p0=(1,1),maxfev=1000)
            a = popt[0]
            k = popt[1]

            # r2
            residuals = mean-self.f(t,*popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((mean-np.mean(mean))**2)
            r2 = 1-(ss_res/ss_tot)

        except RuntimeError:
            a = None
            k = None
            r2 = None
            
        return r2, a, k


    def linear_fit(self, y, x):
        # linear regression model & fit
        model = sm.OLS(y, x)
        result = model.fit()

        # coefficients
        coefficients = result.params
        r2 = result.rsquared
        m=coefficients[1]
        n=coefficients[0]

        return r2, m, n
    

    def calculate_mean(self, data):
        data = data.apply(lambda col: col.fillna(col.mean()), axis=0)
        col = data.columns

        # ref mean
        idx = [j for j in col if self.ref in j]    # index for reference species
        ref = data[idx].mean(axis=1)

        # mean & std across species
        mean = pd.DataFrame()
        mean_dist = {}
        std_dist = {}
        for i in self.species:
            idx = [j for j in col if i in j]
            for j in idx:
                mean[j] = (data[j]-ref)**2

            mean_dist[i] = np.mean(mean[idx].mean(axis=1, skipna=True))
            if len(idx) > 1:
                std_dist[i] = np.mean(mean[idx].std(axis=1, skipna=True))
            else:
                std_dist[i] = 0

        # sort mean, std, and time to get same order
        df = {}
        for s in self.species:
            df[s] = [self.timedata[s], mean_dist[s], std_dist[s]]

        df = pd.DataFrame(df)
        return df


    def plot(self, data):
        df = self.calculate_mean(data)

        #plt.errorbar(df.iloc[0,:], df.iloc[1,:], yerr=df.iloc[2,:], fmt='o', ecolor='b', capsize=1)
        plt.scatter(df.iloc[0,:], df.iloc[1,:])

        x = np.linspace(0, max(df.iloc[0,:])+0.1, 100)

        # linear fit
        r2_l, m, n = self.linear_fit(np.array(df.iloc[1,:]), sm.add_constant(np.array(df.iloc[0,:])))
        plt.plot(x, m*x+n, label=f'linear, r2={r2_l:.2f}')

        # power fit
        r2_p, a, k = self.power_fit(df.iloc[1,:], df.iloc[0,:])
        plt.plot(x, a*x**(k), label = f"power, r2={r2_p:.2f}")

        for i in self.species:
            plt.annotate(i[:4], (df.loc[0, i], df.loc[1,i]), fontsize='small')

        plt.xlabel("evolutionary time")
        plt.ylabel("Mean Squared Expression Distance")
        plt.legend()
        plt.title(f"gene : {data.shape[0]}, {self.title}")
        plt.show()


    def run(self):
        # total genes
        self.plot(data = self.res)

        self.categorize()
        # stabilizing selection
        self.title = "Stabilized Selection"
        self.plot(data = self.stabilized)

        # neutral evolution
        self.title = "Neutral Evolution"
        self.plot(data = self.neutral)

    
    def compare(self):
        gene = self.X.index
        
        result = {"R2_l":np.zeros(len(gene)),"R2_p":np.zeros(len(gene)),
                  "m":np.zeros(len(gene)), "n":np.zeros(len(gene)),
                  "a":np.zeros(len(gene)), "k":np.zeros(len(gene))}

        # ref mean
        idx = [j for j in self.X.columns if self.ref in j]    # index for reference species
        ref = self.X[idx].mean(axis=1)

        # mean & std across species
        mean = pd.DataFrame()
        mean_dist = {}
        for i in self.species:
            idx = [j for j in self.X.columns if i in j]
            for j in idx:
                mean[j] = (self.X[j]-ref)**2
            mean_dist[i] = mean[idx].mean(axis=1, skipna=True)

        for i,g in enumerate(gene):
            if i % 1000 == 0:
                print(i)
            df = pd.DataFrame({
            s: [self.timedata[s], mean_dist[s][i]] for s in self.species
            })
            
            try:
                # linear fit
                result["R2_l"][i], result["m"][i], result["n"][i] = self.linear_fit(np.array(df.iloc[1,:]), sm.add_constant(np.array(df.iloc[0,:])))
                # power fit
                result["R2_p"][i], result["a"][i], result["k"][i] = self.power_fit(df.iloc[1,:], df.iloc[0,:])
            except ValueError:
                result["R2_l"][i], result["m"][i], result["n"][i] = [None, None, None]
                result["R2_p"][i], result["a"][i], result["k"][i] = [None, None, None]

        result = pd.DataFrame(result)
        result.insert(0,"gene_name",gene)

        ## compare with fdr values
        pfit = result[result["R2_l"]<result["R2_p"]]["gene_name"]
        fdr = np.zeros(len(pfit))
        for i, s in enumerate(pfit):
            fdr[i] = self.stat[self.stat["gene_name"]==s]["qvalues"].values

        plt.hist(fdr, color = 'blue', edgecolor = 'black', bins=int(300))
        plt.xlabel("fdr")
        plt.ylabel("number of genes with better fit to power law")
        plt.xlim(0, 1)
        plt.show()

        return result
    
    
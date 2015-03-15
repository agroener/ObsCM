from __future__ import division
import numpy as np
from sklearn import mixture
import matplotlib.pyplot as plt


def GMM(mvir_norm,plotmodelcomparison=False):
    # Set up the gmm fitting
    numgaussians = range(1,31)
    mvir_norm_arr = np.reshape(np.array(mvir_norm),len(mvir_norm),1)
    weight_list = []
    mean_list = []
    cov_list = []
    aic_list = []
    bic_list = []
    model_list = []
    for i in numgaussians:
        g = mixture.GMM(n_components=i)
        g.fit(mvir_norm_arr)
        weight_list.append(g.weights_)
        mean_list.append(g.means_)
        cov_list.append(g.covars_)
        aic_list.append(g.aic(mvir_norm_arr))
        bic_list.append(g.bic(mvir_norm_arr))
        model_list.append(g)

    # Finding the best model according to AIC and BIC
    aic_bestmodel = numgaussians[aic_list.index(min(aic_list))]
    bic_bestmodel = numgaussians[bic_list.index(min(bic_list))]

    if plotmodelcomparison is True:
        # Plotting information criteria
        plt.title("Gaussian Mixture Model")
        plt.plot(numgaussians,aic_list,label='AIC',color='red')
        plt.arrow(aic_bestmodel,min(aic_list)*1.11,0,-0.1*min(aic_list),
                  head_width=0.5, head_length=30, fc='r', ec='r')
        plt.plot(numgaussians,bic_list,label='BIC',color='blue')
        plt.arrow(bic_bestmodel,min(bic_list)*1.11,0,-0.1*min(bic_list),
                  head_width=0.5, head_length=30, fc='b', ec='b')
        plt.xlabel(r'# of components',fontsize=18)
        plt.legend(loc=0)
        plt.show()

        # Plotting a comparison of the best-fit gmm model to the data
        plt.title('Comparison of the Best Fit Model to the Data')
        plt.hist(mvir_norm,bins=150,normed=True,color='blue',label='Data',
                 histtype='step', alpha=0.75)
        plt.hist(model_list[aic_bestmodel].sample(n_samples=1000000),bins=300,normed=True,
                 color='green',label='Model', histtype='step', alpha=0.75)
        plt.legend(loc=0)
        plt.xlim(0,300)
        plt.show()
    mindex_bic = bic_list.index(min(bic_list))
    return weight_list[mindex_bic],mean_list[mindex_bic],cov_list[mindex_bic]

if __name__ == '__main__':
    weights, means, covs = GMM()
    print(weights)
    print(means)
    print(covs)

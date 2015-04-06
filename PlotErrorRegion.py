import matplotlib.pyplot as plt
from numpy.random import normal
from numpy import linspace

import ipdb

def PlotErrorRegion(m,sigm,b,sigb,xmin,xmax):

    xlist = linspace(xmin,xmax,1000)
    
    num_draw = 1000
    m_list = linspace(m-sigm,m+sigm,10)
    b_list = linspace(b-sigb,b+sigb,10)

    plt.plot(xlist,m_list[0]*xlist+b_list[0],color='red') # lowest both
    plt.plot(xlist,m_list[-1]*xlist+b_list[-1],color='blue') # highest both
    plt.plot(xlist,m_list[0]*xlist+b_list[-1],color='green') # lowest m, highest b
    plt.plot(xlist,m_list[-1]*xlist+b_list[0],color='purple') # highest m, lowest b
    
    '''
    for i in range(len(m_list)):
        for j in range(len(b_list)):
            print("m={}  b={}".format(m_list[i],b_list[i]))
            plt.plot(xlist,m_list[i]*xlist+b_list[j])
    '''
    plt.show()
    upperlist = [max([model[i][j] for i in range(len(m_list))]) for j in range(num_draw)]
    lowerlist = [min([model[i][j] for i in range(len(m_list))]) for j in range(num_draw)]
    #ipdb.set_trace()

    


if __name__ == "__main__":
    m = -0.1645
    sigm = 0.17
    b = 3.3572
    sigb = 0.0110 + 0.1991
    xmin = 13.0
    xmax = 17.0
    PlotErrorRegion(m,sigm,b,sigb,xmin,xmax)

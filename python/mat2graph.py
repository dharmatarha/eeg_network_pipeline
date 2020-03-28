from scipy.io import loadmat
import igraph
import numpy as np

def mat2graph(matFile, arrayName='prunedConn', nanValue=0, expShape=(62, 62, 40, 4)):
    '''
    Loads a mat file containing adjacency matrices
    and converts matrices to weighted non-directed graphs

    Inputs:
    matFile     - Valid path to a .mat file, up tp version -v7
    arrayName   - Name of the array containing the adjacency
        matrices of interest, defaults to 'prunedConn'
    '''

    # load data, select adjacency matrices
    matDict = loadmat(matFile)  # returns a dict where keys are variable names
    adjArray = matDict[arrayName]
    # check shape
    if adjArray.shape != expShape:
        print('Adjacency array has unexpected shape:', str(adjArray.shape))
    # convert nans to zero
    adjArray[np.isnan(adjArray)] = nanValue

    # convert each adjacency matrix into a graph
    graphs = []  # init a list that will hold lists of epoch-level graphs
    for stim in range(expShape[-1]):  # loop through conditions
        epochs = []  # init a list for epoch-level graphs
        for epoch in range (expShape[-2]):  # loop through epochs
            # select adj matrix for given epoch in given condition
            adjTmp = adjArray[:, :, epoch, stim]
            # create graph object from a boolean version
            g = igraph.Graph.Adjacency((adjTmp > 0).tolist())
            g.es['weight'] = adjTmp[adjTmp.nonzero()]  # add weights from adj matrix
            g.to_undirected()  # cast into undirected graph
            epochs.append(g)
        graphs.append(epochs)

    return graphs



















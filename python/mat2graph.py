from scipy.io import loadmat
import igraph as ig
import leidenalg as la
import numpy as np


def mat2graph(matfile,
              arrayname='prunedConn',
              nanvalue=0,
              expshape=(62, 62, 40, 4),
              vertexnames=None):
    """
    Loads a mat file containing adjacency matrices
    and converts matrices to weighted non-directed graphs
    (igraph Graph objects)

    Inputs:
    matfile     -- Valid path to a .mat file, up tp version -v7
    arrayname   -- Name of the array containing the adjacency
            matrices of interest, defaults to 'prunedConn'
    nanvalue    -- NaN values are changed to "nanValue", defaults to 0
    expshape    -- Tuple containing the expected shape of the array
            containing adjacency matrices, defaults to (62, 62, 40, 4)
    vertexnames -- List of strings, names of vertices. Defaults to None

    Outputs:
    graphs      -- List of lists, where each element is a list for a condition.
        In each condition list, there is a Graph object for each epoch.
    """

    # load data, select adjacency matrices
    matdict = loadmat(matfile)  # returns a dict where keys are variable names
    adjarray = matdict[arrayname]

    # check shape
    if adjarray.shape != expshape:
        print('Adjacency array has unexpected shape:', str(adjarray.shape))

    # convert nans to zero
    adjarray[np.isnan(adjarray)] = nanvalue

    # convert each adjacency matrix into a graph
    graphs = []  # init a list that will hold lists of epoch-level graphs
    for stim in range(expshape[-1]):  # loop through conditions

        epochs = []  # init a list for epoch-level graphs
        for epoch in range(expshape[-2]):  # loop through epochs

            # select adj matrix for given epoch in given condition
            adjtmp = adjarray[:, :, epoch, stim]

            # create graph object from a boolean version
            g = ig.Graph.Adjacency((adjtmp > 0).tolist())
            g.es['weight'] = adjtmp[adjtmp.nonzero()]  # add weights from adj matrix
            g.to_undirected()  # cast into undirected graph

            # set vertex names if the input arg was supplied
            if vertexnames:
                g.vs['label'] = vertexnames

            epochs.append(g)
        graphs.append(epochs)

    return graphs


def mat2roinames(roinamesfile, roinamesvar='rois'):
    """
    Loads a mat file containing ROI (graph vertex) names in a cell array.

    Inputs:
    roinamesfile  -- Path to .mat file containing ROI names
            in a cell array
    roinamesvar   -- String, name of variable (dict key) containing
            the ROI names

    Outputs:
    roilist       -- List of strings, each string corresponds to a ROI
            (graph vertex)
    """

    # load data, select roi names variable
    matdict = loadmat(roinamesfile)
    rois = matdict[roinamesvar]

    # transform into a list of strings
    roilist = []
    tmp = rois[0].tolist()
    [roilist.append(roi.tolist()[0]) for roi in tmp]

    return roilist











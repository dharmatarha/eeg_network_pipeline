from scipy.io import loadmat
import igraph as ig
import leidenalg as la
import numpy as np
from math import pi
from random import sample


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
        In each condition list, there is a igraph Graph object for each epoch.
    """

    # load data, select adjacency matrices
    matdict = loadmat(matfile)  # returns a dict where keys are variable names
    adjarray = matdict[arrayname]

    # check shape
    if adjarray.shape != expshape:
        print('Adjacency array has unexpected shape:', str(adjarray.shape))

    # convert nans to specified value
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
            g.to_undirected()  # cast into undirected graph
            g.es['weight'] = adjtmp[adjtmp.nonzero()]  # add weights from adj matrix

            # set vertex names if the input arg was supplied
            if vertexnames:
                g.vs['label'] = vertexnames

            epochs.append(g)
        graphs.append(epochs)

    return graphs


def mat2roinames(roinamesfile, roinamesvar='roisShort'):
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


def la_partition(g):
    """
    Placeholder to remind me of important args - use the one-liner with appropriate args in any wrapper
    """

    p = la.find_partition(g, la.ModularityVertexPartition, weights=g.es['weight'], n_iterations=-1)

    return p


def plot_communities(g, partition, layout_edgeratio=0.2):
    """
    Plots igraph graph partition "partition" with a layout that groups
    vertices in each community together.
    Partition is an output from a leidenalg modularity
    detection method on the graph (leidenalg.VertexPartition).
    Works by defining a layout on the graph without cross-community edges,
    than utilizing that layout for plotting the partition.

    Inputs:
    g              -- igraph graph object
    partition      -- leidenalg.VertexPartition object

    Outputs:
    community_plot -- igraph plot object
    """

    # copy the graph
    gcopy = g.copy()
    # empty lists for non-crossing edges and colors of all edges
    edges = []
    edges_colors = []

    # set diff colors for intra- and inter-community edges
    for edge in g.es():
        if partition.membership[edge.tuple[0]] != partition.membership[edge.tuple[1]]:
            edges.append(edge)
            edges_colors.append('gray')
        else:
            edges_colors.append('black')

    # get layout for graph without (most) inter-community edges
    noToDelete = round(len(edges)*(1-layout_edgeratio))  # no of edges to delete from all of edges
    delIdx = sample([i for i in range(len(edges))], noToDelete)  # indices of edges to delete
    delIdx.sort(reverse=True)
    for i in delIdx:
        gcopy.delete_edges(edges[i])
#    [gcopy.delete_edges(edges[i]) for i in delIdx]
#    gcopy.delete_edges(edges[delIdx])
    layout = gcopy.layout_fruchterman_reingold()  # usually works well for eeg connectivity data
    # set edge colors
    g.es['color'] = edges_colors

    # create a visual style dict to be passed on to igraph.plot
    visual_style = {}
    visual_style['vertex_shape'] = 'circle'
    visual_style['vertex_size'] = 10
    visual_style['vertex_label_dist'] = 15
    visual_style['vertex_label_angle'] = 1.5*pi
    visual_style['edge_color'] = g.es['color']
    visual_style['edge_width'] = g.es['weight']*10
    visual_style['layout'] = layout
    visual_style['bbox'] = (1024, 768)
    visual_style['margin'] = (100, 100, 100, 200)

    # plot
    community_plot = ig.plot(partition, **visual_style)

    return community_plot


def partitionwrapper(matfile,
                     roinamesfile,
                     roinamesvar='roisShort',
                     arrayname='prunedConn',
                     expshape=(62, 62, 40, 4),
                     n_partitioning=100):
    """
    Wrapper for community / modularity detection on EEG data using the Leiden algorithm.
    Relies on other functions in mat2graph (mat2roinames, mat2graph) and
    calls leidenalg functions on each network (on the connectivity network of each epoch).

    matfile        -- Valid path to a .mat file, up to version -v7
    roinamesfile   -- Path to .mat file containing ROI names
            in a cell array
    roinamesvar    -- String, name of variable (dict key) containing
            the ROI names            
    arrayname      -- Name of the array containing the adjacency
            matrices of interest, defaults to 'prunedConn'
    expshape       -- Tuple containing the expected shape of the array
            containing adjacency matrices, defaults to (62, 62, 40, 4)
    n_partitioning -- No. of community detections to run on each network

    Outputs:
    partitions     -- list of lists of leidenalg.VertexPartition objects, one for each network (epoch)
    memberships    -- numpy array of membership indexes of graph vertices, its shape is
    (expshape[0], expshape[-2], expshape[-1])
    """

    # load roi names
    roiList = mat2roinames(roinamesfile, roinamesvar=roinamesvar)

    # load connectivity matrices, define graphs
    graphs = mat2graph(matfile, arrayname=arrayname, expshape=expshape, vertexnames=roiList)

    # init output variables
    partitions = []
    memberships = np.zeros((expshape[0], expshape[-2], expshape[-1]))

    # calculate modularity on each graph n_partitioning times
    for stim in range(len(graphs)):
        # user message
        print('\n\nCalculating for stim/cond no. ' + str(stim))
        # init list holding epoch-level partitions
        epochPartitions = []

        for epoch in range(len(graphs[0])):
            # user message
            print('Epoch no. ' + str(epoch))

            # select graph of epoch
            g = graphs[stim][epoch]

            # init temporary variables for partitioning runs
            modulValues = np.zeros((n_partitioning, 1))  # will hold modularity values from each partitioning run
            tmpPartitions = []  # empty list for partitioning outcomes

            # as community detection is somewhat random (and maximization is not guaranteed),
            # # more runs are used for each network
            for p in range(n_partitioning):
                tmp = la.find_partition(g,
                                        la.ModularityVertexPartition,
                                        weights=g.es['weight'],
                                        n_iterations=-1)
                tmpPartitions.append(tmp)
                modulValues[p] = tmp.modularity

            # find maximal modularity and corresponding partitioning
            idx = np.where(modulValues == modulValues.max())
            print(idx)
#            if idx[0].shape[0] > 1:
#                idx = tmpIdx[0][0]
#            else:
#                idx = tmpIdx[0]
            maxPartition = tmpPartitions[idx[0][0]]

            # store results, both the max-modularity partition and the corresponding membership array
            epochPartitions.append(maxPartition)
            memberships[:, epoch, stim] = np.asarray(maxPartition.membership)

        # store epoch-level partitions
        partitions.append(epochPartitions)

    return partitions, memberships








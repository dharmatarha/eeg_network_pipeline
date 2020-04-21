from scipy.io import loadmat
import igraph as ig
import leidenalg as la
import numpy as np
from random import randint

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

def la_partition(g):

    p = la.find_partition(g, la.ModularityVertexPartition, weights=g.es['weight'], n_iterations=-1)

    return p

def plot_communities(g, p):

    return

def _plot(g, membership=None):
    if membership is not None:
        gcopy = g.copy()
        edges = []
        edges_colors = []
        for edge in g.es():
            if membership[edge.tuple[0]] != membership[edge.tuple[1]]:
                edges.append(edge)
                edges_colors.append("gray")
            else:
                edges_colors.append("black")
        gcopy.delete_edges(edges)
        layout = gcopy.layout("kk")
        g.es["color"] = edges_colors
    else:
        layout = g.layout("kk")
        g.es["color"] = "gray"
    visual_style = {}
    visual_style["vertex_label_dist"] = 0
    visual_style["vertex_shape"] = "circle"
    visual_style["edge_color"] = g.es["color"]
    # visual_style["bbox"] = (4000, 2500)
    visual_style["vertex_size"] = 30
    visual_style["layout"] = layout
    visual_style["bbox"] = (1024, 768)
    visual_style["margin"] = 40
    visual_style["edge_label"] = g.es["weight"]
    for vertex in g.vs():
        vertex["label"] = vertex.index
    if membership is not None:
        colors = []
        for i in range(0, max(membership) + 1):
            colors.append('%06X' % randint(0, 0xFFFFFF))
        for vertex in g.vs():
            vertex["color"] = str('#') + colors[membership[vertex.index]]
        visual_style["vertex_color"] = g.vs["color"]
    ig.plot(g, **visual_style)

 #   if __name__ == "__main__":
 #       g = igraph.Nexus.get("karate")
 #       cl = g.community_fastgreedy()
 #       membership = cl.as_clustering().membership
 #       _plot(g, membership)










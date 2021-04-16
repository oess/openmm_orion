# Acknowledgement
# The code has been adapted from the MIT licensed package
# https://github.com/choderalab/freeenergyframework
# developed by Hannah Bruce MacDonald from the Chodera Lab


import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from MDOrion.FEC.RFEC.freeenergyframework import stats


def mle(g, factor='f_ij'):
    """
    Compute maximum likelihood estimate of free energies and covariance in their estimates.
    The number 'factor' is the node attribute on which the MLE will be calculated,
    where d'factor' will be used as the standard error of the factor

    We assume the free energy of node 0 is zero.

    Reference : https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00528
    Xu, Huafeng. "Optimal measurement network of pairwise differences." Journal of Chemical Information and Modeling 59.11 (2019): 4720-4728.

    Parameters
    ----------
    g : nx.Graph
        The graph for which an estimate is to be computed
        Each edge must have attributes 'f_ij' and 'df_ij' for the free energy and uncertainty estimate
        Will have 'bayesian_f_ij' and 'bayesian_df_ij' added to each edge
        and 'bayesian_f_i' and 'bayesian_df_i' added to each node.
    factor : string, default = 'f_ij'
        node attribute of nx.Graph that will be used for MLE
    Returns
    -------
    f_i : np.array with shape (n_ligands,)
        f_i[i] is the absolute free energy of ligand i in kT
        f_i[0] = 0

    C : np.array with shape (n_ligands, n_ligands)
        C[i,j] is the covariance of the free energy estimates of i and j

    """
    N = len(g.nodes)
    f_ij = stats.form_edge_matrix(g, factor, action='antisymmetrize')
    df_ij = stats.form_edge_matrix(g, f'd{factor}', action='symmetrize')

    # Form F matrix (Eq 4)
    F = np.zeros([N,N])
    for (i,j) in g.edges:
        F[i,j] = - df_ij[i,j]**(-2)
        F[j,i] = - df_ij[i,j]**(-2)
    for i in g.nodes:
        F[i,i] = - np.sum(F[i,:])

    # Form z vector (Eq 3)
    z = np.zeros([N])
    for (i,j) in g.edges:
        z[i] += f_ij[i,j] * df_ij[i,j]**(-2)
        z[j] += f_ij[j,i] * df_ij[j,i]**(-2)

    # Compute MLE estimate (Eq 2)
    Finv = np.linalg.pinv(F)
    f_i = - np.matmul(Finv, z) # NOTE: This differs in sign from Eq. 2!
    f_i[:] -= f_i[0]

    # Compute uncertainty
    C = Finv
    return f_i, C


class Result(object):
    def __init__(self, ligandA, ligandB,
                 calc_DDG, mbar_error, other_error,
                 exp_DDG=None, exp_dDDG=None):
        self.ligandA = str(ligandA).strip()
        self.ligandB = str(ligandB).strip()
        self.calc_DDG = float(calc_DDG)
        self.mbar_dDDG = float(mbar_error)
        self.other_dDDG = float(other_error)
        self.dcalc_DDG = self.mbar_dDDG+self.other_dDDG  # is this definitely always additive?
        if exp_DDG is not None:
            self.exp_DDG = float(exp_DDG)
        if exp_dDDG is not None:
            self.dexp_DDG = float(exp_dDDG)


class FEMap_with_exp(object):

    def __init__(self, results):
        self.results = results
        self.graph = nx.DiGraph()
        self.n_edges = len(results)

        self.generate_graph_from_results()

        # check the graph has minimal connectivity

    def generate_graph_from_results(self):
        self._name_to_id = {}
        id = 0
        for result in self.results:
            if result.ligandA not in self._name_to_id.keys():
                self._name_to_id[result.ligandA] = id
                id += 1
            if result.ligandB not in self._name_to_id.keys():
                self._name_to_id[result.ligandB] = id
                id += 1
            # TODO need some exp error for mle to converge for exp... this is a horrible hack
            if result.dexp_DDG == 0.0:
                result.dexp_DDG = 0.01
            if result.dcalc_DDG == 0.0:
                result.dcalc_DDG = 0.01
            self.graph.add_edge(self._name_to_id[result.ligandA], self._name_to_id[result.ligandB],
            exp_DDG=result.exp_DDG, dexp_DDG=result.dexp_DDG,
            calc_DDG=result.calc_DDG, dcalc_DDG=result.dcalc_DDG)

        self.n_ligands = self.graph.number_of_nodes()
        self.degree = self.graph.number_of_edges() / self.n_ligands

        # check the graph has minimal connectivity
        self.check_weakly_connected()
        if not self.weakly_connected:
            print('Graph is not connected enough to compute absolute values')
        else:
            self.generate_absolute_values()

    def check_weakly_connected(self):
        undirected_graph = self.graph.to_undirected()
        self.weakly_connected = nx.is_connected(undirected_graph)
        return nx.is_connected(undirected_graph)

    def generate_absolute_values(self):
        if self.weakly_connected:
            f_i_exp, C_exp = stats.mle(self.graph, factor='exp_DDG')
            variance = np.diagonal(C_exp)
            for i, (f_i, df_i) in enumerate(zip(f_i_exp, variance**0.5)):
                self.graph.nodes[i]['f_i_exp'] = f_i
                self.graph.nodes[i]['df_i_exp'] = df_i

            f_i_calc, C_calc = stats.mle(self.graph, factor='calc_DDG')
            variance = np.diagonal(C_calc)
            for i, (f_i, df_i) in enumerate(zip(f_i_calc, variance**0.5)):
                self.graph.nodes[i]['f_i_calc'] = f_i
                self.graph.nodes[i]['df_i_calc'] = df_i

    def draw_graph(self, title='', filename=None):
        plt.figure(figsize=(10, 10))
        self._id_to_name = {}
        for i, j in self._name_to_id.items():
            self._id_to_name[j] = i
        nx.draw_circular(self.graph, labels=self._id_to_name, node_color='hotpink', node_size=250)
        long_title = f'{title} \n Nedges={self.n_edges} \n Nligands={self.n_ligands} \n Degree={self.degree:.2f}'
        plt.title(long_title)
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename, bbox_inches='tight')



class FEMap(object):

    def __init__(self, results):
        self.results = results
        self.graph = nx.DiGraph()
        self.n_edges = len(results)

        self.generate_graph_from_results()

        # check the graph has minimal connectivity

    def generate_graph_from_results(self):
        self._name_to_id = {}
        id = 0
        for result in self.results:
            if result.ligandA not in self._name_to_id.keys():
                self._name_to_id[result.ligandA] = id
                id += 1
            if result.ligandB not in self._name_to_id.keys():
                self._name_to_id[result.ligandB] = id
                id += 1
            if result.dcalc_DDG == 0.0:
                result.dcalc_DDG = 0.01
            self.graph.add_edge(self._name_to_id[result.ligandA], self._name_to_id[result.ligandB],
            calc_DDG=result.calc_DDG, dcalc_DDG=result.dcalc_DDG)

        self.n_ligands = self.graph.number_of_nodes()
        self.degree = self.graph.number_of_edges() / self.n_ligands

        # check the graph has minimal connectivity
        self.check_weakly_connected()
        if not self.weakly_connected:
            print('Graph is not connected enough to compute absolute values')
        else:
            self.generate_absolute_values()

    def check_weakly_connected(self):
        undirected_graph = self.graph.to_undirected()
        self.weakly_connected = nx.is_connected(undirected_graph)
        return nx.is_connected(undirected_graph)

    def generate_absolute_values(self):
            f_i_calc, C_calc = stats.mle(self.graph, factor='calc_DDG')
            variance = np.diagonal(C_calc)
            for i, (f_i, df_i) in enumerate(zip(f_i_calc, variance**0.5)):
                self.graph.nodes[i]['f_i_calc'] = f_i
                self.graph.nodes[i]['df_i_calc'] = df_i

    def draw_graph(self, title='', filename=None):
        plt.figure(figsize=(10, 10))
        self._id_to_name = {}
        for i, j in self._name_to_id.items():
            self._id_to_name[j] = i
        nx.draw_circular(self.graph, labels=self._id_to_name, node_color='hotpink', node_size=250)
        long_title = f'{title} \n Nedges={self.n_edges} \n Nligands={self.n_ligands} \n Degree={self.degree:.2f}'
        plt.title(long_title)
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename, bbox_inches='tight')


def read_csv(filename):
    raw_results = []
    with open(filename,'r') as f:
        for line in f:
            if line[0] != '#':
                raw_results.append(Result(*line.split(',')))
    return raw_results

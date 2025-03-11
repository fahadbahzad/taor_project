# THIS MODULE CONTAINS THE CODE FOR RUNNING STOCHASTIC ALGORITHMS.
# Currently this only contains a stochastic implementation of the recursive alg.

# IMPORTS
import xpress as xp
# xp.init('C:/Program Files/xpressmp/bin/xpauth_personal.xpr') # For licensing.

import pandas as pd
import math
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import time

def ra(pairs:list, altruistic_donors:list, edges:dict, k:int=3, noisy:int=1, scenarios=None):
    """ Implements the recursive algorithm from Anderson et al. 2015.
    
    Args:
        pairs (list): donor-recipient pairs
        altruistic_donors (list): NDDs
        edges (dict): keys are tuples of length 2, each of which is an element
            of pairs or of altruistic_donors; values are the weights of the
            directed edge from the 0th entry to 1st entry of the tuple.
        k (int): max cycle length allowed in the optimal solution.
        noisy (int): determines how much to print as output during optimization:
            0 - switch off ALL printed output
            1 - switch off solver output log; report finding feasible solutions
            2 - show all printed output
        scenarios (list of dictionaries): 


    Returns:
        opt_val (float): value of the optimal solution; if scenarios==None this
            value is NOT calculated stochastically, otherwise it is computed
            as the expected value against all provided scenarios.
        solution_edges (list): a list of edges in the optimal solution
        time_taken (float): time taken for the optimization (excl.)
        VSS (float): value of the stochastic solution; None if scenarios==None
        opt_sol (dict): keys are edges, value is 1 if included in optimal
            solution and 0 if not.
    """


    nodes = pairs + altruistic_donors

    # Create Xpress Model
    # Initialize the model
    prob = xp.problem()


    # Define decision variables for each edge
    y = {e: xp.var(vartype=xp.binary, name=f"y_{e[0]}_{e[1]}") for e in edges}
    prob.addVariable(list(y.values()))

    # FIRST STAGE constraints
    for v in pairs:
        prob.addConstraint(xp.Sum(y[e] for e in edges if e[0] == v) <= xp.Sum(y[e] for e in edges if e[1] == v))
        prob.addConstraint(xp.Sum(y[e] for e in edges if e[1] == v) <= 1)

    for a in altruistic_donors:
        prob.addConstraint(xp.Sum(y[e] for e in edges if e[0] == a) <= 1)


    if scenarios:

        # First solve the problem deterministically, finding the edge selections
        # This is used in calculating the value of the stochastic solution
        _, _, _, _, y_D = ra(pairs, altruistic_donors, edges,k,noisy=0,scenarios=None)

        # Define decision variables for second stage (after node deletion)
        x = {(e,s) : xp.var(vartype=xp.binary, name=f"x_{e[0]}_{e[1]}_{s}") 
        for e in edges for s,_ in enumerate(scenarios)}

        prob.addVariable(list(x.values()))

        # SECOND STAGE constraints:
        for v in pairs:
            for s, _ in enumerate(scenarios):
                prob.addConstraint(xp.Sum(x[(e,s)] for e in edges if e[0] == v) <= xp.Sum(x[e,s] for e in edges if e[1] == v))
                prob.addConstraint(xp.Sum(x[(e,s)] for e in edges if e[1] == v) <= 1)
        
        for a in altruistic_donors:
            for s, _ in enumerate(scenarios):
                prob.addConstraint(xp.Sum(x[e,s] for e in edges if e[0] == a) <= 1)

        for s, S in enumerate(scenarios):
            for e in edges:
                prob.addConstraint(x[e,s] <= S[e]*y[e])

        objective = xp.Sum(w*x[e,s]/len(scenarios) for e, w in edges.items() for s, _ in enumerate(scenarios))

    else:
        # If we don't have scenarios, the objective is just first stage.
        objective = xp.Sum(y[e] * w for e, w in edges.items())

    # Objective: Maximize total benefit
    prob.setObjective(objective, sense=xp.maximize)



    finished = False # A flag to mark the end of the optimization.
    infeasible = False # A (currently unused) flag to mark no feasible solution.
    # TODO: Add a catch for no feasible solution. This shouldn't happen but it is
    # probably better to be robust about this.

    if noisy in [0,1]:
        prob.controls.outputlog = 0 # This just makes it quiet to run

    Gstart_time = time.time()

    while finished == False and infeasible == False:
        

        start_time = time.time()


        # Solve the model
        
        prob.solve()
        opt_sol = prob.getSolution()
        
        # print(f"Solver took {time.time()-start_time} seconds to run")
        start_time = time.time()

        # Construct the graph from the optimal solution:
        DG = nx.DiGraph()
        selected_edges = [list(edges.keys())[i] for i, e in enumerate(list(edges.keys())) if opt_sol[i]>0.5]
        DG.add_edges_from(selected_edges)

        # Check if there is a cycle length that is too long:
        cycles = list(nx.simple_cycles(DG))
        
        # print(f"Finding the cycles took {time.time()-start_time} seconds")
        start_time = time.time()

        # If ok, report done.
        # TODO: Rewrite this so that it is with the max_cycle bit...
        if cycles==[] or max(map(len,cycles))<=k:
            if not noisy == 0:
                print("")
                print("##########################################################")
                print(f"OPTIMIZATION COMPLETE: no cycles of length more than {k}.")
                print("##########################################################")
                print("")
            finished = True
            break
        
        # If not done, report that reoptimization is required:
        else:
            if not noisy in [0,1]:
                print("")
                print("#################################################################")
                print(f"REOPTIMIZATION REQUIRED: proposed solution contains long cycles.")
                print("#################################################################")
                print("")

        # Take the long cycle we found and make a note of its edges:
        max_cycle = max(cycles,key=len)
        cycle_edges = [(max_cycle[i],max_cycle[i+1]) for i in range(len(max_cycle)-1)]
        cycle_edges += [(max_cycle[-1],max_cycle[0])]

        # Add the constraint to remove this as an option: 
        prob.addConstraint(xp.Sum(y[e] for e in cycle_edges) <= len(max_cycle)-1)

    opt_sol = prob.getSolution(y)
    opt_val = prob.getObjVal()
    
    solution_edges = [e for e in edges if prob.getSolution(y[e]) > 0.05]

    time_taken = time.time()-Gstart_time

    if not noisy in [0]:
        print(f"The optimization took {time_taken} seconds to execute.")

        # Print the output
        print("")
        print("")
        print("")

        # print("Optimal Matches:")
        # for (u, v), var in y.items():
        #     if prob.getSolution(var) > 0.5:
        #         print(f"{u} donates to {v} with benefit {edges[(u,v)]}")


        print(f"Total Benefit: {opt_val}")
        print("")
        print("")
        print("")

    if scenarios:
        if not noisy in [0,1]: print("Calculating VSS...")
        prob.addConstraint(x[e,s]<=y_D[e] for e in edges for s, _ in enumerate(scenarios))
        prob.solve()
        val_det_sol = prob.getObjVal()
        VSS = opt_val - val_det_sol
        print(f"Expected value of DETERMINISTIC solution: {val_det_sol}")
        print(f"Expected value of STOCHASTIC solution: {opt_val}")
        print(f"VSS: {VSS}")
    else: VSS = None

    return opt_val, solution_edges, time_taken, VSS, opt_sol


def get_S(edges, deleted_nodes):
    ''' 
    Takes edges and deleted_nodes, returns S encoding which edges survive.
    
    S is a dictionary: the keys are the edges of the compatibility graph.
    S[e] is 1 if the edge e survives the deletion of nodes and 0 otherwise.
    '''
    
    S = {e : 0 if e[0] in deleted_nodes or e[1] in deleted_nodes else 1 for e in edges}
    return S


def generate_simple_scenarios(pairs, altruistic_donors, edges):
    """
    Generate a list of dictionaries which encode edge survival.

    In particular, in THIS function, the scenarios generated result from the
    deletion of a single edge. All such scenarios are built.

    Args:
        pairs: list
        altruistic_donors: list
        edges: dictionary (keys: 2-tuple (in/out vx); values: weight of edge)
        N: number of scenarios to generate

    Returns:
        scenarios: a list of dictionaries. Each dictionary is a scenario; in
            each dictionary, keys are edges, value is 1 if edge survives else 0.
    """

    nodes = pairs+altruistic_donors

    scenarios = []

    # Generate all scenarios in which exactly one node drops out:
    for node in nodes:
        deleted_nodes = [node]
        S = get_S(edges, deleted_nodes)
        scenarios.append(S)

    return scenarios


def generate_probabilistic_scenarios(pairs, altruistic_donors, edges, N, P_dropout_pairs, P_dropout_altruistic):
    ''' Return N scenarios based on a probabilities of nodes dropping out.
    
    Args:
        pairs: list
        altruistic_donors: list
        edges: dictionary (keys: 2-tuple (in/out vx); values: weight of edge)
        N: number of scenarios to generate
        P_dropout_pairs: (fixed) probability that a pair drops out
        P_dropout_altruistic: (fixed) probability that an NDD drops out

    Returns:
        scenarios: a list of dictionaries. Each dictionary is a scenario; in
            each dictionary, keys are edges, value is 1 if edge survives else 0.
    '''

    scenarios = []

    for i in range(N):
        deleted_nodes = [pair for pair in pairs if np.random.uniform() <= P_dropout_pairs] +\
                [NDD for NDD in altruistic_donors if np.random.uniform() <= P_dropout_altruistic]
        S = get_S(edges, deleted_nodes)
        scenarios.append(S)

    return scenarios

def generate_edge_dropout_scenarios(edges, N, P_edge_dropout):
    ''' Return N scenarios based on edge dropout probability P_edge_dropout.
    
    Args:
        edges: dictionary (keys: 2-tuple (in/out vx); values: weight of edge)
        N: number of scenarios to generate
        P_edge_dropout: (fixed) probability that an edge drops out

    Returns:
        scenarios: a list of dictionaries. Each dictionary is a scenario; in
            each dictionary, keys are edges, value is 1 if edge survives else 0.
    '''
    
    scenarios = []

    for i in range(N):
        S = {e : 0 if np.random.uniform() <= P_edge_dropout else 1 for e in edges}
        scenarios.append(S)

    return scenarios
# WIP

# IMPORTS
import xpress as xp
import pandas as pd
import math
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import time

def ra(pairs:list, altruistic_donors:list, edges:dict, noisy:int=0, stochastic:int=0):
    """ Implements the recursive algorithm from Anderson et al. 2015.
    
    Args:
        pairs (list): donor-recipient pairs
        altruistic_donors (list): NDDs
        edges (dict): keys are tuples of length 2, each of which is an element
            of pairs or of altruistic_donors; values are the weights of the
            directed edge from the 0th entry to 1st entry of the tuple.
        noisy (int): determines how much to print as output during optimization.
        stochastic (int): determines what form of stochastic modelling should be
            employed, if any. Default 0 is no stochasticity.


    Returns:
        opt_val
        solution_edges (dict):
        time_taken (float):

    """


    nodes = pairs + altruistic_donors

    # Create the loop!

    # Create Xpress Model
    # Initialize the model
    prob = xp.problem()

    # Define decision variables for each edge
    y = {e: xp.var(vartype=xp.binary, name=f"y_{e[0]}_{e[1]}") for e in edges}
    prob.addVariable(list(y.values()))

    # Objective: Maximize total benefit
    prob.setObjective(xp.Sum(y[e] * w for e, w in edges.items()), sense=xp.maximize)

    # Constraints
    for v in pairs:
        prob.addConstraint(xp.Sum(y[e] for e in edges if e[0] == v) <= xp.Sum(y[e] for e in edges if e[1] == v))
        prob.addConstraint(xp.Sum(y[e] for e in edges if e[1] == v) <= 1)

    for a in altruistic_donors:
        prob.addConstraint(xp.Sum(y[e] for e in edges if e[0] == a) <= 1)

    finished = False # A flag to mark the end of the optimization.
    infeasible = False # A (currently unused) flag to mark no feasible solution.
    # TODO: Add a catch for no feasible solution. This shouldn't happen but it is
    # probably better to be robust about this.

    Gstart_time = time.time()


    # prob.solve()

    while finished == False and infeasible == False:
        

        start_time = time.time()


        # Solve the model
        prob.controls.outputlog = 0 # This just makes it quiet to run
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
            print("")
            print("##########################################################")
            print(f"OPTIMIZATION COMPLETE: no cycles of length more than {k}.")
            print("##########################################################")
            print("")
            finished = True
            break
        
        # If not done, report that reoptimization is required:
        else:
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



    print(f"The optimization took {time.time() - Gstart_time} seconds to execute.")




    # Print the output
    print("")
    print("")
    print("")

    # print("Optimal Matches:")
    # for (u, v), var in y.items():
    #     if prob.getSolution(var) > 0.5:
    #         print(f"{u} donates to {v} with benefit {edges[(u,v)]}")

    solution_edges = [e for e in edges if prob.getSolution(y[e]) > 0.05]

    print(f"Total Benefit: {prob.getObjVal()}")
    print("")
    print("")
    print("")
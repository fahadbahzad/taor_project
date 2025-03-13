import xpress as xp
import networkx as nx
import pandas as pd
import numpy as np
import time


def separation_algorithm(y, f_i):
    """
    Implements the separation algorithm for cut-set constraints based on the supplementary information.
    
    >>> Inputs:
    y:                    Dictionary where keys are edge tuples (u, v) and values are their solution values i.e., y[e] gotten via current_sol[0] from above
    f_i:                  Dictionary where keys are nodes and values are the solver solution for f_i[v] gotten via current_sol[2] from above

    >>> Outputs:
    delta_minus_S:       List of edges i.e., (u,v) that are in the named set as per the paper
    v:                   Node at which the constraint is violated for
    """
    
    # Get the global Graph object from Networkx of the current problem.
    # Note: Although this is not best practice, it is more efficient rather 
    #       than creating the graph every single iteration. 
    G = globals().get("G", None)
    
    # Add edges with weights from y (solution values)
    # Note: as per the documentation adding an edge that already exists updates the edge data.
    for (u,v), sol_value in y.items():
        G.add_edge(u, v, weight=sol_value)
    
    # Solve the max-flow min-cut problem for each node with f_i[v] > 0
    for v, flow_in in f_i.items():
        if flow_in > 0:
            
            # Compute min-cut using the function from networkx
            cut_value, (Not_S, S) = nx.minimum_cut(G, "SUPER", v, capacity="weight")
            
            # If the cut weight is less than flow_in, add a violated constraint
            if cut_value < flow_in:

                # Define delta_minus_S in this case
                delta_minus_S = [(u, v) for u, v in G.edges(Not_S) if v in S]

                # Return a list that contains the violated constraint's data needed to add the cut into 
                # Xpress using the call back.
                return(delta_minus_S, v)


    # Note: This version returns the first violated constraint found; A possible extension is 
    #       to try if its better to return all violated constraints in a single potential solution. 
    return()



def separation_cbpreintsol(prob, data, soltype, cutoff):
    """
    Implements the separation callback when an integer solution is found by heuristics
    or during the branch and bound search, but before it is accepted by the optimizer.
    
    >>> Inputs:
    prob:                    Xpress model object to be passed to the callback function.
    data:                    Data we want to pass into the callback. 
    soltype:                 Type of MIP solution that has been found from 0-3 
    cutoff:                  Current cutoff value for the solution.

    >>> Outputs:
    ifreject:                Binary; 1 if the solution should be rejected and 0 otherwise
    newcutoff:               New cutoff value, to be used by the optimizer if the solution is accepted. 

    Note: The inputs and outputs are required and predetermined in this structure based on Xpress. 
    """
    
    # Initialisze outputs:
    ifreject = 0
    newcutoff = None # Note: this is not being used currently, we will potentially explore it in extensions
    
    # Get LP relaxation solution
    lp_solution = []
    prob.getlpsol(lp_solution, None, None, None)

    # Populate values from lp_solution using iteration which should be faster as it is just an O(1) process.
    # Note: we could stop at f_i_temp if it gets too slow. 

    solution_iter = iter(lp_solution)

    for key in y_temp:
        y_temp[key] = next(solution_iter)
    
    for key in z_temp:
        z_temp[key] = next(solution_iter)
    
    for key in f_i_temp:
        f_i_temp[key] = next(solution_iter)
    
    # for key in f_o_temp:
    #     f_o_temp[key] = next(solution_iter)
    
    # Get the violated constraint data using the separation algorithm from above
    violated_constraint_data = separation_algorithm(y_temp, f_i_temp)

    # Check if there is a violated constraint
    if violated_constraint_data:

        # Reject the solution
        ifreject = 1
        
        # There is a violated constraint so return 1/True i.e., reject this solution
        # Note: From the research I've found the purpose of this callback is to just reject solutions
        #       and not to add constraints to the problem dynamically. This is particularly relevant
        #       for Andrew's part because he'll probably have to do it recursively outside the solver. 
            
    # Otherwise, return 0/False i.e., this solution does not violate a constraint so its good!
    return(ifreject, newcutoff)






def separation_cboptnode(prob, data):
    """
    Implements the separation callback when called during the branch and bound search, after the LP
    relaxation has been solved for the current node, and after any internal cuts and heuristics have been
    applied, but before the optimizer checks if the current node should be branched. 
        
    >>> Inputs:
    prob:                    Xpress model object to be passed to the callback function.
    data:                    Data we want to pass into the callback. 
    
    >>> Outputs:
    feas:                    Feasibility status; If set to a nonzero value the current node is set to infeasible

    Note: The inputs and outputs are required and predetermined in this structure based on Xpress. 
    """
    
    # Initialise output 
    feas = 0
    
    # Get LP relaxation solution
    lp_solution = []
    prob.getlpsol(lp_solution, None, None, None)
    
    # Populate values from lp_solution using iteration which should be faster as it is just an O(1) process.
    # Note: we could stop at f_i_temp if it gets too slow. 
    solution_iter = iter(lp_solution)
    
    for key in y_temp:
        y_temp[key] = next(solution_iter)
    
    for key in z_temp:
        z_temp[key] = next(solution_iter)
    
    for key in f_i_temp:
        f_i_temp[key] = next(solution_iter)
    
    # for key in f_o_temp:
    #     f_o_temp[key] = next(solution_iter)

    
    # Get the violated constraint data using the separation algorithm from above
    violated_constraint_data = separation_algorithm(y_temp, f_i_temp)
    
    # If there is a violated constraint add it as a cut
    if violated_constraint_data:
        
        # Unpack the tuple to get the data we want
        delta_minus_S, v = violated_constraint_data  

        # Translate the cut into the presolve index
        colind, rowcoef = [], []
        
        drhsp, status = prob.presolverow(      rowtype = "G", # 'G' for >= constraints
                                               origcolind = [prob.getIndex(y[e]) for e in delta_minus_S] + [prob.getIndex(f_i[v])], # Index values of original variables
                                               origrowcoef = [1] * len(delta_minus_S) + [-1], # Coefficients for original variables
                                               origrhs = 0, # Right-hand side of constraint in original variables
                                               maxcoefs=prob.attributes.cols, # what does this do?
                                               colind=colind, rowcoef=rowcoef # where to output the new ones
                                        ) 
    
        # Now we add the translated cut
        prob.addcuts(
                            cuttype=[1],  # General cut
                            rowtype=['G'],  # Presolved row type
                            rhs=[drhsp],  # Presolved RHS
                            start=[0, len(colind)],  # Start indices
                            colind=colind,  # Presolved column indices
                            cutcoef=rowcoef  # Presolved coefficients
                    )

    # Continue solving
    return(feas)



def pctsp(graph:object, pairs:list, altruistic_donors:list, nodes:list, edges:dict, all_cycles:dict, noisy:int=1):
    """
    Solves the kidney exchange problem using the prize collecting travelling salesman problem approach.
    
    >>> Inputs:
    graph:                 The graph object representing the network structure.
    pairs:                 A list of recipient-donor pairs.
    altruistic_donors:     A list of altruistic donors.
    nodes:                 A list of all nodes in the graph i.e., pairs + alturistic_donors.
    edges:                 A dictionary mapping of edges with their weights.
    all_cycles:            A dictionary containing all allowed cycles.
    noisy:                 Binary to control the output printing of Xpress. 

    >>> Outputs:
    opt_val:               The optimal value of the exchanges in the solution.
    selected_edges:        A list of selected edges that are chosen in the solution.
    selected_cycles:       A list of selected cycles chosen in the solution.
    time_taken:            The total time taken to solve the problem.
    """



    global G
    G = graph

    # Define the optimization model
    model = xp.problem()

    # Decision variables
    global y, z, f_i, f_o  # Declare the variables as global to be accessible
    y = {e: xp.var(vartype=xp.binary, name=f"y_{e}") for e in edges}  # Edge selection
    z = {c: xp.var(vartype=xp.binary, name=f"z_{c}") for c in all_cycles}  # Cycle selection
    f_i = {v: xp.var(vartype=xp.binary) for v in nodes}  # Flow in decision variable
    f_o = {v: xp.var(vartype=xp.binary) for v in nodes}  # Flow out decision variable

    # Add decision variables
    model.addVariable(list(y.values()) + list(z.values())+ list(f_i.values()) + list(f_o.values()))



    # Xpress uses indexing when in callback thats why we need to create a dictionary for the ids. I suspect when you run this on actual data
    # you would want to do this in another code cell just so that it isn't repeated every time you solve the model for debugging. 
    # e.g., id_vars[("y", ('NDD1', 'P1') )] = 0 i.e., the decision variable to connect nodes NDD1 and P1 is the first decision variable in Xpress.

    # Initialize id_vars dictionary and counter to assign sequential values
    id_vars = {}
    counter = 0

    # Populate id_vars with ("y", key) tuples
    for key in y.keys():
        id_vars[("y", key)] = counter
        counter += 1
        

    # Populate id_vars with ("z", key) tuples
    for key in z.keys():
        id_vars[("z", key)] = counter
        counter += 1

    # Populate id_vars with ("f_i", key) tuples
    for key in f_i.keys():
        id_vars[("f_i", key)] = counter
        counter += 1

    # Populate id_vars with ("f_o", key) tuples
    for key in f_o.keys():
        id_vars[("f_o", key)] = counter
        counter += 1


    # Create temporary storage for the callback solutions to be used:
    global y_temp, z_temp, f_i_temp, f_o_temp  # Declare the variables as global to make it more
    y_temp = {e: 0 for e in edges}
    z_temp = {c: 0 for c in all_cycles} 
    f_i_temp = {v: 0 for v in nodes}  
    f_o_temp = {v: 0 for v in nodes}  


    # Define the objective function which is to maximize sum of weights of selected edges and cycles
    model.setObjective(xp.Sum(edges[e] * y[e] for e in edges) + xp.Sum(all_cycles[c] * z[c] for c in all_cycles), sense=xp.maximize)

    ### Constraints:

    # 1. Defining f_i and f_o i.e., the incoming and outgoing kidneys to a node v
    for v in nodes:
        model.addConstraint(xp.Sum([y[e] for e in edges if e[1] == v]) == f_i[v])
        model.addConstraint(xp.Sum([y[e] for e in edges if e[0] == v]) == f_o[v])

    # 2. Ensure that if a node is in an actived cycle, the associated edges that involve it are turned off
    for v in pairs:
        model.addConstraint(f_o[v] + xp.Sum([z[c] for c in all_cycles if v in c]) <= f_i[v] + xp.Sum([z[c] for c in all_cycles if v in c]))
        model.addConstraint(f_i[v] + xp.Sum([z[c] for c in all_cycles if v in c]) <= 1)

    # 3. Ensure that altruistic donors donate at most one kidney  
    for v in altruistic_donors:
        model.addConstraint(f_o[v] <= 1)

    # Add the call back function to detect and add violated constraints during LP relaxation
    model.addcbpreintsol(separation_cbpreintsol, None, 3)  
    model.addcboptnode(separation_cboptnode, None, 1)

    model.controls.outputlog = noisy # Toggle the output
    model.setControl("MIPRELSTOP", 0.1)
    model.setControl("maxtime", 600)

    # Solve the model
    start_time = time.time()
    model.solve()
    end_time = time.time()

    opt_val = model.getObjVal()
    selected_edges = [e for e in edges if model.getSolution(y[e]) > 0.05]
    selected_cycles = [c for c in all_cycles if model.getSolution(z[c]) > 0.05]
    time_taken = end_time - start_time

    # print(model.getProbStatusString())
    
    return (opt_val, selected_edges, selected_cycles, time_taken)
import networkx as nx
import pandas as pd

def get_data(file, pc_tsp):
    '''
    Inputs:
    file: a string containing the filename for the compatability information. Of the form 'Data.xlsx'.
    pc_tsp: boolean. TRUE if cycles need to be detected for the pc-tsp algorithm, FALSE otherwise.

    Outputs:
    pairs: List of donor pairs.
    altruistic_donors: List of altruistic donors.
    nodes: List containing both pairs and altruistic donors.
    edges: Dictionary of who can donate to who and the corresponding weights.
    all_cycles: Dictionary of the cycles in terms of donation and their total weights.
    '''

    # Read the Excel file, turn into dataframe
    pool = pd.read_excel(file, engine='openpyxl')

    # Turn into Dataframe without unnecessary columns
    df = pool.drop(['dage','source'], axis=1)

    # Create list of donor pairs
    pairs = list(set(df.loc[pd.isna(df['altruistic']), 'donor_id'].tolist()))

    # Create list of altruistic donors
    altruistic_donors = list(set(df.loc[df['altruistic'] == 1.0, 'donor_id'].tolist()))

    # Create list of nodes
    nodes = pairs + altruistic_donors

    # Create edges dictionary with weights
    edges = df.dropna(subset=['score']).set_index(['donor_id', 'recipient'])['score'].to_dict()

    # If pc_tsp = TRUE then create a graph to extract the cycles
    if pc_tsp == True:
        G = nx.DiGraph() 
        
        # Convert to numeric and handle NaN
        pool['donor_id'] = pd.to_numeric(pool['donor_id'], errors='coerce')
        pool['recipient'] = pd.to_numeric(pool['recipient'], errors='coerce')
        pool['score'] = pd.to_numeric(pool['score'], errors='coerce')

        # Create a list where each donor_id maps to a list of recipients
        matches = []
        for _, row in pool.iterrows():
            source = row['donor_id']
            recipient = row['recipient']
            score = row['score']

            # Ensure valid donor & recipient
            if not pd.isna(source) and not pd.isna(recipient):  
                matches.append((source, recipient, score))
        
        # Add nodes and edges to the graph
        for source, recipient, score in matches:
            G.add_node(source)
            G.add_node(recipient)
            G.add_edge(source, recipient, score=score)

        # Get cycles
        # Create cycles dictionary
        cycles = sorted(nx.simple_cycles(G, length_bound=3))

        # Very particular output to make sure we can turn it into a dictionary
        with open("output.txt", "w") as f: 
            for cycle in cycles:
                print(",".join(map(str, cycle)), file=f)

        # Read the file and convert each line into a tuple
        with open("output.txt", "r") as f:
            cycles = [tuple(map(float, line.strip().split(","))) for line in f if line.strip()]

        # Convert list of tuples into a dictionary (no weights yet)
        cycles_dict = {i: cycle for i, cycle in enumerate(cycles)}

        # Function to calculate the sum of weights for a given cycle
        def get_cycle_weight(cycle, edges):
            total_weight = 0
            for i in range(len(cycle) - 1):
                edge = (cycle[i], cycle[i + 1])
                if edge in edges:
                    total_weight += edges[edge]
                
            # Add the weight for the edge from the last node to the first node
            closing_edge = (cycle[-1], cycle[0])
            if closing_edge in edges:
                total_weight += edges[closing_edge]
                
            return total_weight

        # Create dictionary with cycles and weights
        all_cycles = {tuple(cycle): get_cycle_weight(cycle, edges) for cycle in cycles_dict.values()}

    else:
        all_cycles = []

    return(pairs, altruistic_donors, nodes, edges, all_cycles)

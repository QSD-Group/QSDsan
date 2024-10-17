# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 10:21:03 2024

@author: Junhyung Park
"""
import numpy as np
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

# Load the matrix from the Excel file
matrix_path = 'export_import_matrix.xlsx'
matrix = pd.read_excel(matrix_path, index_col=0)

# Convert the matrix to a NetworkX graph
G = nx.from_pandas_adjacency(matrix)

# Prepare an Excel writer
with pd.ExcelWriter('network_analysis_results1.xlsx', engine='xlsxwriter') as writer:
    
    # --- Task 1: Calculate #nodes, #links, and density ---
    # Calculate the number of nodes
    num_nodes = G.number_of_nodes()
    print(f"Number of nodes: {num_nodes}")

    # Calculate the number of links (edges)
    num_links = G.number_of_edges()
    print(f"Number of links: {num_links}")

    # Calculate the density of the network
    density = nx.density(G)
    print(f"Density: {density}")

    # Save task 1 results to Excel
    task_1_data = pd.DataFrame({
        'Metric': ['Number of nodes', 'Number of links', 'Density'],
        'Value': [num_nodes, num_links, density]
    })
    task_1_data.to_excel(writer, sheet_name='Task 1', index=False)

    #%% --- Task 2: Plot node degree distribution ---
    # Get degree of each node
    degrees = dict(G.degree())

    # Save the degrees to Excel
    degrees_df = pd.DataFrame(list(degrees.items()), columns=['Country', 'Degree'])
    degrees_df.to_excel(writer, sheet_name='Task 2 - Degrees', index=False)

    # Plot the degree distribution
    plt.figure(figsize=(8, 6))
    plt.hist(list(degrees.values()), bins=range(1, max(degrees.values()) + 2), edgecolor='black')
    plt.title('Node Degree Distribution', fontsize=20, fontweight='bold')  # Increased font size and bold
    plt.xlabel('Degree', fontsize=18, fontweight='bold')  # Increased font size and bold
    plt.ylabel('Frequency', fontsize=18, fontweight='bold')  # Increased font size and bold
    plt.xticks(fontsize=16, fontweight='bold')  # Increased font size and bold for x-ticks
    plt.yticks(fontsize=16, fontweight='bold')  # Increased font size and bold for y-ticks
    plt.grid(True)
    plt.show()

    #%% --- Task 3: Plot node strength distribution ---
    # Get strength of each node (weighted degree, if available)
    node_strengths = dict(G.degree(weight='weight'))

    # Save the node strengths to Excel
    strength_df = pd.DataFrame(list(node_strengths.items()), columns=['Country', 'Strength'])
    strength_df.to_excel(writer, sheet_name='Task 3 - Strengths', index=False)

    # Plot the node strength distribution
    plt.figure(figsize=(8, 6))
    plt.hist(list(node_strengths.values()), bins=30, edgecolor='black')
    plt.title('Node Strength Distribution', fontsize=20, fontweight='bold')  # Increased font size and bold
    plt.xlabel('Strength', fontsize=18, fontweight='bold')  # Increased font size and bold
    plt.ylabel('Frequency', fontsize=18, fontweight='bold')  # Increased font size and bold
    plt.xticks(fontsize=16, fontweight='bold')  # Increased font size and bold for x-ticks
    plt.yticks(fontsize=16, fontweight='bold')  # Increased font size and bold for y-ticks
    plt.grid(True)
    plt.show()

    #%% --- Task 4: Rank top 10 countries in terms of clustering and betweenness centrality ---
    # Clustering coefficient for each node
    clustering = nx.clustering(G, weight='weight')
    top_10_clustering = sorted(clustering.items(), key=lambda x: x[1], reverse=True)[:10]

    print("Top 10 countries by clustering:")
    for country, score in top_10_clustering:
        print(f"{country}: {score}")

    # Save clustering results to Excel
    clustering_df = pd.DataFrame(top_10_clustering, columns=['Country', 'Clustering Coefficient'])
    clustering_df.to_excel(writer, sheet_name='Task 4 - Clustering', index=False)

    # Betweenness centrality for each node
    betweenness = nx.betweenness_centrality(G, weight='weight')
    top_10_betweenness = sorted(betweenness.items(), key=lambda x: x[1], reverse=True)[:10]

    print("\nTop 10 countries by betweenness centrality:")
    for country, score in top_10_betweenness:
        print(f"{country}: {score}")

    # Save betweenness results to Excel
    betweenness_df = pd.DataFrame(top_10_betweenness, columns=['Country', 'Betweenness Centrality'])
    betweenness_df.to_excel(writer, sheet_name='Task 4 - Betweenness', index=False)

    #%% --- Task 5: Visualize the network ---
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, k=0.15)
    nx.draw_networkx_nodes(G, pos, node_size=500, node_color='blue', alpha=0.7)
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size=14, font_color='black')  # Increased font size for labels
    plt.title("Network Visualization of Trade Relationships", fontsize=20, fontweight='bold')  # Larger and bold title
    plt.show()

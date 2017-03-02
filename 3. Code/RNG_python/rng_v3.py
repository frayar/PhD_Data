# -*- coding: cp1252 -*-
# http://commons.wikimedia.org/wiki/File:Relative_neighborhood_graph.svg

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pyx import canvas,path,color
from networkx.readwrite import json_graph

# GLOBALS
G = nx.Graph()
points = []
edgelist = []


# Euclidian distance (in 2D)
# We can avoid Graph access by working directly on points
def dist2(n1,n2):
    return (G.node[n1]['x']-G.node[n2]['x'])**2 + (G.node[n1]['y']-G.node[n2]['y'])**2
 

# Check if p and are relative neighbors 
def neighbors(n1,n2):
    for n in G:
        if max(dist2(n1,n),dist2(n2,n)) < dist2(n1,n2):
            return False
    return True

# Display graph wrt coordinate (x,y)
def display():
    c = canvas.canvas()
    radius = 0.05
    
    # Draw points as red filled circles 
    for p in points:
        c.fill(path.circle(p[0],p[1],radius),[color.rgb.red])

    # Draw edges
    for e in edgelist:
        c.stroke(path.line(G.node[e[0]]['x'], G.node[e[0]]['y'], G.node[e[1]]['x'], G.node[e[1]]['y']),[color.rgb.black])

    # Export canvas 
    c.writePDFfile("Relative_neighborhood_graph")

# Export edges
def export():
    # Export format
    # http://networkx.lanl.gov/reference/readwrite.html
    
    # Write in txt file
    fp = open(".\\output.txt", "a")
    json_graph.dump(G, fp, indent=3)
    fp.close()

# Main 
def main():
    global points
    global edgelist
    
    # Load points without class label
    data =  pd.read_table('.\\zahn.txt', sep='\t', header=None)
    points = data.ix[:,0:1].values

    # Add points as nodes
    count = 0
    for p in points:
        G.add_node(count)
        G.node[count]['x'] = p[0]
        G.node[count]['y'] = p[1]
        count = count + 1

    # Compute edges of RNG(points)
    for n1 in G:
        for n2 in G:
            if n1 < n2 and neighbors(n1,n2):
                G.add_edge(n1,n2)

    # Get computed RNG edges
    edgelist = G.edges()

    print(len(edgelist))
                
    # Draw
    #nx.draw(G)
    #plt.show()

    # Display wrt coordinate (x,y) for edges validation
    #display()

    # Export
    #export()

    
if __name__ == "__main__":
	main()

% Assuming we have adjacency matrices for the fractures and their intersections, we can create 
% graph objects as follows

function = plot_fracture_graph(Floe)

    fracture_graph = graph(Adjacency_Fractures)
    intersections_graph = graph(Adjacency_Interactions)

    plot(fracture_graph)
    plot(intersections_graph)

end
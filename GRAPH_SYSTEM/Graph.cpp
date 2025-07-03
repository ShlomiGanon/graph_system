#include "Graph.h"


// Constructor: initializes graph with given 'both_ways' flag and reserves space for vertices
Graph::Graph(bool both_ways, int initial_size) : both_ways_graph(both_ways)
{
	m_vertices.reserve(initial_size); // Reserve space to avoid reallocations
}

// Constructor: initializes graph with a list of vertex names and 'both_ways' flag
Graph::Graph(const std::vector<std::string>& vertecies_names, bool both_ways) : both_ways_graph(both_ways)
{
	m_vertices.reserve(vertecies_names.size()); // Reserve space for all vertices

	// Create a new Vertex for each name and add to m_vertices
	for (auto name : vertecies_names)
	{
		Vertex* v = new Vertex(name);
		m_vertices.push_back(v);
	}
}

// Copy constructor: creates a deep copy of another graph
Graph::Graph(const Graph& other) : both_ways_graph(other.both_ways_graph)
{
	m_vertices.reserve(other.m_vertices.size()); // Reserve space for vertices

	// Map from vertex name to new vertex pointer (for edge reconstruction)
	std::unordered_map<std::string, Vertex*> name_to_ptr(other.m_vertices.size());

	// Copy vertices from 'other'
	for (const Vertex* v : other.m_vertices)
	{
		Vertex* vertex = new Vertex(v->name);

		// Handle memory allocation failure
		if (vertex == nullptr)
		{
			// Clean up already allocated vertices
			for (int i = 0; i < m_vertices.size(); i++)
			{
				delete m_vertices[i];
			}
			throw Graph_Exception("Memory allocation fail!");
		}
		else
		{
			name_to_ptr[vertex->name] = vertex; // Save mapping
			m_vertices.push_back(vertex);
		}
	}

	// Copy edges between new vertices using the mapping
	for (const Vertex* v : other.m_vertices)
	{
		Vertex* vertex = name_to_ptr[v->name];
		for (auto v_of_out : v->Out_Edges)
		{
			// Add edge to corresponding new vertex with the same weight
			vertex->AddEdgeTo(*name_to_ptr[v_of_out.first->name], v_of_out.second);
		}
	}
}

// Move constructor: transfers ownership of vertices and flags from 'other'
Graph::Graph(Graph&& other) : m_vertices(std::move(other.m_vertices)), both_ways_graph(other.both_ways_graph)
{
	other.both_ways_graph = false; // Mark 'other' as cleared or inactive
	// The moved-from graph will no longer hold vertices
}

// Destructor: deletes all allocated vertices to free memory
Graph::~Graph()
{
	std::cout << "~Graph" << std::endl;

	// Delete each vertex pointer in m_vertices
	for (auto v : m_vertices)
	{
		delete v;
	}
}


Vertex* Graph::AddVertex(const std::string& vertex_Name)
{
	// Try to find an existing vertex with the given name
	Vertex* vertex = GetVertex(vertex_Name);

	// If vertex already exists, return it
	if (vertex != nullptr)
		return vertex;

	// Otherwise, create a new vertex with the given name
	vertex = new Vertex(vertex_Name);

	// Add the new vertex to the graph's vertex list
	m_vertices.push_back(vertex);

	// Return pointer to the newly created vertex
	return vertex;
}


Vertex* Graph::GetVertex(const std::string& vertex_name) const
{
	// Iterate over all vertices in the graph
	for (auto v : m_vertices)
	{
		// If vertex name matches the requested name, return pointer to vertex
		if (v->name == vertex_name)
			return v;
	}

	// If no matching vertex is found, return nullptr
	return nullptr;
}


bool Graph::ConnectEdge(const std::string& from_vertex, const std::string& to_vertex, double weight)
{
	// Add or get the 'from' vertex (creates it if it doesn't exist)
	Vertex* from_v = AddVertex(from_vertex);

	// Add or get the 'to' vertex (creates it if it doesn't exist)
	Vertex* to_v = AddVertex(to_vertex);

	// Avoid connecting a vertex to itself
	if (from_v != to_v)
	{
		// Add an edge from 'from_v' to 'to_v' with given weight
		// If the graph is both ways, add the reverse edge as well
		from_v->AddEdgeTo(*to_v, weight, both_ways_graph);
		return true; // Connection succeeded
	}

	return false; // No connection made (same vertex)
}

bool Graph::ConnectEdges(const std::string& from_vertex, std::vector<std::string> to_vertecies, double weight)
{
	// Return false if there are no target vertices to connect
	if (to_vertecies.size() == 0)
		return false;

	// Try to connect each target vertex to the source vertex
	for (const auto& name : to_vertecies)
	{
		// If any connection fails, return false immediately
		if (!ConnectEdge(from_vertex, name, weight))
			return false;
	}

	return true; // All connections succeeded
}


bool Graph::RemoveVertex(const std::string& vertex_name)
{
	// Check if the vertex name is not empty
	if (vertex_name.length() > 0)
	{
		// Iterate over all vertices to find the matching vertex
		for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
		{
			auto v = *it;

			// If the vertex name matches
			if (v->name == vertex_name)
			{
				delete v;           // Free the memory of the vertex
				m_vertices.erase(it); // Remove vertex from the vector
				return true;        // Indicate successful removal
			}
		}
	}
	return false; // Vertex not found or empty name given
}


std::vector<std::string> Graph::getVertexNames() const
{
	std::vector<std::string> names;
	names.reserve(m_vertices.size()); // Reserve space for all vertex names

	// Loop over all vertices and collect their names
	for (auto v : m_vertices)
	{
		names.push_back(v->name); // Add vertex name to the list
	}

	return names; // Return the vector of vertex names
}

std::size_t Graph::size() const
{
	return m_vertices.size(); // Return the number of vertices in the graph
}


std::vector<Edge> Graph::getEdges() const
{
	std::vector<Edge> Edges;
	Edges.reserve(m_vertices.size()); // Reserve space, rough estimate (could be up to V^2 - V)

	// Iterate over all vertices
	for (auto v : m_vertices)
	{
		// For each outgoing edge of the vertex
		for (auto out_v : v->Out_Edges)
		{
			// Add an Edge object from v to the neighbor with the given weight
			Edges.emplace_back(v, out_v.first, out_v.second);
		}
	}

	return Edges; // Return all edges in the graph
}

const std::vector<Vertex*>& Graph::getVertices() const
{
	return m_vertices; // Return reference to the vector of all vertices
}

bool Graph::IsBothWays() const
{
	return both_ways_graph; // Return true if graph is undirected (both ways), false otherwise
}


std::unordered_map<Vertex*, int> Graph::BFS_distances(const std::string& source)
{
	// Get the start vertex by name
	Vertex* start_vertex = GetVertex(source);

	// Throw exception if the vertex doesn't exist
	if (start_vertex == nullptr)
		throw Graph_Exception("vertex not found: " + source);

	std::unordered_map<Vertex*, int> dists; // Map to store distances from source
	dists.reserve(m_vertices.size());

	// Initialize all distances to infinite
	for (auto v : m_vertices)
	{
		dists[v] = INFINITE_DISTANCE;
	}

	dists[start_vertex] = 0; // Distance to source is zero

	std::unordered_map<Vertex*, bool> visited; // Track visited vertices
	std::vector<Vertex*> vertecies; // Queue substitute for BFS traversal (using vector as queue)
	vertecies.reserve(m_vertices.size());

	// Insert start vertex into queue
	vertecies.insert(vertecies.begin(), start_vertex);

	// BFS main loop
	while (!vertecies.empty())
	{
		Vertex* current = vertecies.back(); // Get the next vertex from the back of vector (acting as queue)
		vertecies.pop_back();

		// If current vertex is not visited yet
		if (visited.find(current) == visited.end())
		{
			visited[current] = true; // Mark as visited

			// Iterate over all neighbors
			for (auto out_pair : current->Out_Edges)
			{
				auto neighbor = out_pair.first;

				// If neighbor not visited yet
				if (visited.find(neighbor) == visited.end())
				{
					vertecies.push_back(neighbor);       // Add neighbor to queue
					dists[neighbor] = dists[current] + 1; // Update distance to neighbor
				}
			}
		}
	}

	return dists; // Return map of distances from source vertex
}


std::vector<Vertex*> Graph::DFS(const std::string& source)
{
	// Get the start vertex by name
	Vertex* start_vertex = GetVertex(source);

	// If the start vertex doesn't exist, throw an exception
	if (start_vertex == nullptr)
		throw Graph_Exception("vertex not found: " + source);

	std::vector<Vertex*> vertecies; // Stores the visited vertices in DFS order
	vertecies.reserve(m_vertices.size());

	std::unordered_map<Vertex*, bool> visited; // Tracks visited vertices
	visited.reserve(m_vertices.size());

	std::stack<Vertex*> vertecies_stack; // Stack used for DFS traversal

	vertecies_stack.push(start_vertex); // Start from the source vertex

	// Main DFS loop
	while (!vertecies_stack.empty())
	{
		Vertex* node = vertecies_stack.top(); // Get the top node from the stack
		vertecies_stack.pop();

		// If this node has not been visited yet
		if (visited.find(node) == visited.end())
		{
			vertecies.push_back(node); // Add to result list
			std::cout << node->name << " , "; // Print the visited node
			visited[node] = true; // Mark as visited

			// Push all unvisited neighbors to the stack
			for (auto out_pair : node->Out_Edges)
			{
				auto neighbor = out_pair.first;
				if (visited.find(neighbor) == visited.end())
					vertecies_stack.push(neighbor);
			}
		}
	}

	std::cout << std::endl;

	return vertecies; // Return the list of visited vertices in DFS order
}


void Graph::DFS_recursive(const std::string start)
{
	// Get the starting vertex by name
	Vertex* start_vertex = GetVertex(start);

	// Throw an exception if the vertex was not found
	if (start_vertex == nullptr)
		throw Graph_Exception("vertex not found: " + start);

	std::unordered_map<Vertex*, bool> visited; // Map to track visited vertices
	std::stack<Vertex*> recStack; // Stack to track the current recursion path (for cycle detection)

	std::cout << "DFS-recursive: ";

	// Start the recursive DFS traversal
	DFS_rec(start_vertex, visited, recStack);

	std::cout << std::endl;
}


void Graph::DFS_rec(Vertex* current, std::unordered_map<Vertex*, bool>& visited, std::stack<Vertex*>& recStack)
{
	std::stack<Vertex*> helper_stack; // Temporary stack to help check for cycles

	// Check if 'current' already exists in the recursion stack (cycle detection)
	while (!recStack.empty())
	{
		Vertex* v = recStack.top();
		recStack.pop();

		if (v == current)
		{
			std::cout << "Cycle detected at " << v->name << std::endl;
			return; // Cycle found, exit the recursion
		}

		helper_stack.push(v); // Save popped elements temporarily
	}

	// Restore elements back to recStack after checking
	while (!helper_stack.empty())
	{
		Vertex* v = helper_stack.top();
		helper_stack.pop();
		recStack.push(v);
	}

	// If already visited, no need to explore again
	if (visited.find(current) != visited.end())
		return;

	std::cout << current->name << " "; // Print the current vertex

	visited[current] = true; // Mark current as visited

	recStack.push(current); // Add current to the recursion stack

	// Recursively visit all neighbors
	for (auto out_pair : current->Out_Edges)
	{
		auto neighbor = out_pair.first;
		DFS_rec(neighbor, visited, recStack);
	}

	// Remove 'current' from recStack after all neighbors are processed
	while (!recStack.empty())
	{
		Vertex* v = recStack.top();
		recStack.pop();

		if (v != current)
		{
			helper_stack.push(v); // Store other elements
		}
	}

	// Restore the stack without 'current'
	while (!helper_stack.empty())
	{
		Vertex* v = helper_stack.top();
		helper_stack.pop();
		recStack.push(v);
	}
}


std::vector<Vertex*> Graph::TopologicalSort_DFS()
{
	std::unordered_map<Vertex*, bool> visited; // Tracks visited vertices
	std::stack<Vertex*> order; // Stack to store the topological order in reverse

	// Perform DFS from each unvisited vertex
	for (auto node : m_vertices)
	{
		if (visited.find(node) == visited.end()) // If the node was not visited yet
		{
			TopologicalSort_DFS_func(node, visited, order); // Explore the component
		}
	}

	std::vector<Vertex*> result;
	result.reserve(order.size()); // Reserve space for the result vector

	// Extract vertices from the stack to form the correct topological order
	while (!order.empty())
	{
		result.push_back(order.top());
		order.pop();
	}

	return result; // Return the topological ordering of the vertices
}


int Graph::getConnectedComponents() const
{
	if (both_ways_graph == false)
		return -1; // Cannot calculate connected components on directed graph

	int result = 0; // Counter for connected components

	std::unordered_map<Vertex*, bool> visited; // Tracks visited vertices

	std::stack<Vertex*> order; // Not needed for counting, but passed to DFS helper

	// Go through all vertices
	for (auto node : m_vertices)
	{
		// If the node was not visited yet, it belongs to a new component
		if (visited.find(node) == visited.end())
		{
			// Explore the entire component using DFS
			TopologicalSort_DFS_func(node, visited, order);
			result++; // One more connected component found
		}
	}

	return result; // Return total number of connected components
}


void Graph::TopologicalSort_DFS_func(Vertex* node, std::unordered_map<Vertex*, bool>& visited, std::stack<Vertex*>& order) const
{
	visited[node] = true; // Mark the current node as visited

	// Recursively visit all unvisited neighbors
	for (auto out_pair : node->Out_Edges)
	{
		auto neighbor = out_pair.first;
		if (!visited[neighbor]) // If the neighbor has not been visited
		{
			TopologicalSort_DFS_func(neighbor, visited, order); // Visit neighbor
		}
	}

	// After all neighbors are processed, push the current node onto the stack
	order.push(node); // This ensures that dependencies are pushed before the node itself
}


std::list<Vertex*> Graph::TopologicalSort()
{
	std::list<Vertex*> order; // Final topological order of vertices

	std::unordered_map<Vertex*, int> inDegree; // Map to store in-degree for each vertex
	inDegree.reserve(m_vertices.size());

	// Step 1: Initialize all in-degrees to 0
	for (auto v : m_vertices)
	{
		inDegree[v] = 0;
	}

	// Step 2: Count in-degrees of each vertex based on incoming edges
	for (auto node : m_vertices)
	{
		for (auto out_pair : node->Out_Edges)
		{
			auto neighbor = out_pair.first;
			inDegree[neighbor] += 1; // Increase in-degree for neighbor
		}
	}

	std::queue<Vertex*> Q; // Queue for processing vertices with in-degree 0

	// Step 3: Enqueue all vertices with in-degree 0
	for (auto v : m_vertices)
	{
		if (inDegree.find(v) != inDegree.end() && inDegree[v] == 0)
		{
			Q.push(v); // Add vertex to queue if it has no incoming edges
		}
	}

	// Step 4: Process the queue
	while (!Q.empty())
	{
		Vertex* current = Q.front();
		Q.pop();
		order.push_back(current); // Add current vertex to topological order

		// Decrease in-degree of all neighbors
		for (auto out_pair : current->Out_Edges)
		{
			auto neighbor = out_pair.first;
			inDegree[neighbor] -= 1;

			// If in-degree becomes 0, add to queue
			if (inDegree[neighbor] == 0)
			{
				Q.push(neighbor);
			}
		}
	}

	// Step 5: Check for cycle (if not all vertices were processed)
	if (order.size() < m_vertices.size())
		std::cout << "Cycle detected" << std::endl; // Graph is not a DAG

	return order; // Return the topological order (may be partial if cycle exists)
}


std::vector<std::pair<Vertex*, std::vector<std::pair<Vertex*, double>>::iterator>> Graph::GetPath_BFS(const std::string& FROM, const std::string& TO) const
{
	Vertex* from_vertex = GetVertex(FROM); // Get pointer to source vertex
	Vertex* to_vertex = GetVertex(TO);     // Get pointer to destination vertex

	// Check if either vertex is not found
	if (from_vertex == nullptr || to_vertex == nullptr)
	{
		std::stringstream ss;
		if (from_vertex == nullptr)
		{
			ss << "Source vertex not found [FROM: " << FROM << "]";
		}
		if (to_vertex == nullptr)
		{
			if (from_vertex == nullptr) ss << " and ";
			ss << "Destination vertex not found [TO: " << TO << "]";
		}
		throw Graph_Exception(ss.str()); // Throw error message
		return {};
	}

	// Special case: source and destination are the same
	if (from_vertex == to_vertex)
	{
		for (auto it = from_vertex->Out_Edges.begin(); it != from_vertex->Out_Edges.end(); ++it)
		{
			if (it->first == from_vertex)
				return { { from_vertex , it } }; // Return loop to self if such edge exists
		}
		return {}; // No self-loop found
	}

	std::queue<Vertex*> Q; // Queue for BFS
	std::unordered_map<Vertex*, bool> visited; // Track visited vertices
	std::unordered_map<Vertex*, std::pair<Vertex*, std::vector<std::pair<Vertex*, double>>::iterator>> parent;
	// Stores parent and edge iterator for each visited vertex

	Q.push(from_vertex);           // Start BFS from source
	visited[from_vertex] = true;   // Mark source as visited
	bool found = false;            // Flag to mark if target was reached

	// Standard BFS loop
	while (!Q.empty() && !found)
	{
		Vertex* current = Q.front(); // Get next vertex in queue
		Q.pop();

		// Explore all outgoing edges from current vertex
		for (auto it = current->Out_Edges.begin(); it != current->Out_Edges.end(); it++)
		{
			Vertex* neigthbor = (*it).first;
			if (!visited[neigthbor])
			{
				parent[neigthbor] = std::make_pair(current, it); // Store parent and edge used to reach neighbor
				if (neigthbor == to_vertex)
				{
					found = true; // Destination reached
					break;
				}
				Q.push(neigthbor);            // Enqueue neighbor
				visited[neigthbor] = true;    // Mark as visited
			}
		}
	}

	// If destination was not found, return empty path
	if (!found)
	{
		return {};
	}

	// Build path by backtracking from destination to source using parent map
	std::stack<Vertex*> reverse_path;
	Vertex* temp = to_vertex;
	while (temp != from_vertex)
	{
		reverse_path.push(temp);
		temp = parent[temp].first;
	}

	int vertecis = reverse_path.size();
	std::vector<std::pair<Vertex*, std::vector<std::pair<Vertex*, double>>::iterator>> result(vertecis);

	// Reconstruct the path in correct order and store parent + iterator
	for (int i = 0; i < vertecis; i++)
	{
		Vertex* v = reverse_path.top();
		Vertex* v_perant = parent[v].first;
		std::vector<std::pair<Vertex*, double>>::iterator it = parent[v].second;
		result[i] = { v_perant , it };
		reverse_path.pop();
	}

	return result; // Return the reconstructed path from source to destination
}




std::vector<double> Graph::Dijastra(const std::string& source) const
{
	std::unordered_map<const Vertex*, double> dist;     // Stores the shortest distance from source to each vertex
	std::unordered_map<const Vertex*, bool> visited;    // Tracks whether a vertex was already visited

	Vertex* source_vertex = GetVertex(source);          // Get pointer to the source vertex by name
	if (source_vertex == nullptr)
		throw Graph_Exception("cant find the vertex: " + source); // Throw an exception if the source vertex is not found

	dist.reserve(m_vertices.size()); // Reserve space for all vertices in the distance map

	// Initialize distances and visited map
	for (auto v : m_vertices)
	{
		dist[v] = (v == source_vertex) ? 0 : INFINITE_DISTANCE; // Distance is 0 for source, infinite for others
		visited[v] = false; // Mark all vertices as unvisited at start
	}

	while (true) // Main loop of Dijkstra's algorithm
	{
		const Vertex* u = nullptr;

		// Step 1: Find the unvisited vertex with the smallest known distance
		for (auto dist_pair : dist)
		{
			const Vertex* v = dist_pair.first;
			if (!visited[v])
			{
				if (u == nullptr || dist_pair.second < dist[u])
					u = v;
			}
		}

		// Step 2: If there is no such vertex, or the smallest distance is infinite, stop the loop
		if (u == nullptr || dist[u] == INFINITE_DISTANCE)
			break;

		visited[u] = true; // Mark current vertex as visited

		// Step 3: For each neighbor v of u, try to relax the edge u -> v
		for (auto u_out_pair : u->Out_Edges)
		{
			const Vertex* v = u_out_pair.first;
			double newDist = dist[u] + u_out_pair.second;

			// If the current vertex u is unreachable, skip updates
			if (dist[u] == INFINITE_DISTANCE)
				newDist = INFINITE_DISTANCE;

			// If the new calculated distance is shorter, update it
			if (newDist < dist[v])
				dist[v] = newDist;
		}
	}

	// Step 4: Convert the distance map to a vector in the same order as m_vertices
	std::vector<double> dist_vec;
	dist_vec.reserve(size());
	for (auto v : m_vertices)
	{
		dist_vec.push_back(dist[v]);
	}

	return dist_vec;
}


std::vector<double> Graph::Bellman_Ford(const std::string& source) const
{
	Vertex* source_vertex = GetVertex(source); // Get pointer to the source vertex by name
	if (source_vertex == nullptr) throw Graph_Exception("cant find the vertex: " + source); // If vertex not found, throw exception

	std::unordered_map<const Vertex*, double> dist; // Map to store distances from source to each vertex
	auto Edges = this->getEdges(); // Get list of all edges in the graph

	// Step 1: Initialize all distances to infinity
	for (auto v : m_vertices)
	{
		dist[v] = INFINITE_DISTANCE;
	}
	dist[source_vertex] = 0; // Distance from source to itself is 0

	// Step 2: Relax all edges (V - 1) times
	// This loop ensures that the shortest paths are calculated by allowing paths with up to (V - 1) edges
	int size_V = m_vertices.size();
	for (int i = 1; i < size_V; i++) // Repeat the process V-1 times
	{
		const Vertex* u, * v;
		double w, newDist;
		for (const auto& e : Edges) // Go through each edge in the graph
		{
			u = e.from;
			v = e.to;
			w = e.weight;
			newDist = dist[u] + w;

			// If current distance to u is not infinite and going through u gives a shorter path to v
			if (dist[u] != INFINITE_DISTANCE && dist[v] > newDist)
			{
				dist[v] = newDist; // Update the distance to vertex v
			}
		}
	}

	// Step 3: Check for negative weight cycles
	// If we can still relax an edge, there is a negative cycle
	for (const auto& e : Edges)
	{
		const Vertex* u = e.from;
		const Vertex* v = e.to;
		double w = e.weight;
		double newDist = dist[u] + w;

		// If shorter path still exists, then a negative cycle is detected
		if (dist[u] != INFINITE_DISTANCE && dist[v] > newDist)
		{
			// Return {0.0} to indicate a negative cycle was found
			return { 0.0 };
		}
	}

	std::vector<double> dist_vec;
	dist_vec.reserve(size_V);

	// Convert distance map to a vector in the order of m_vertices
	for (auto v : m_vertices)
	{
		dist_vec.push_back(dist[v]);
	}

	return dist_vec;
}


std::vector<std::vector<double>> Graph::Floyd_Warshall() const
{
	std::unordered_map<Vertex*, int> vToIndex;
	int V_size = m_vertices.size();

	// Create a vector filled with infinite values representing unreachable distances
	std::vector<double> INF_VEC(V_size, INFINITE_DISTANCE);

	// Initialize the distance matrix with infinite values
	std::vector<std::vector<double>> dist(V_size, INF_VEC);

	// Assign an index to each vertex for easier access in the matrix
	for (int i = 0; i < V_size; i++)
	{
		vToIndex[m_vertices[i]] = i;
	}

	// Step 1: Set initial distances based on direct connections between vertices
	for (int i = 0; i < V_size; i++)
	{
		for (auto& i_out_pair : m_vertices[i]->Out_Edges)
		{
			Vertex* u = i_out_pair.first;
			dist[i][vToIndex[u]] = i_out_pair.second; // Set direct edge weight from i to u
		}
		dist[i][i] = 0; // Distance from a vertex to itself is zero
	}

	// Step 2: Update the distance matrix using the Floyd-Warshall algorithm
	// Try every vertex k as an intermediate point between all pairs of vertices (i, j)
	for (int k = 0; k < V_size; k++) // Loop over all possible intermediate vertices
	{
		// For each source vertex i
		for (int i = 0; i < V_size; i++)
		{
			// For each destination vertex j
			for (int j = 0; j < V_size; j++)
			{
				// Check if path i -> k and k -> j exists (not infinite)
				if (dist[i][k] < INFINITE_DISTANCE && dist[k][j] < INFINITE_DISTANCE)
				{
					// If going through k gives a shorter path from i to j, update it
					dist[i][j] = std::min(dist[i][j], dist[i][k] + dist[k][j]);
				}
			}
		}
	}

	// Note: This algorithm does not detect negative weight cycles
	return dist;
}


std::vector<std::vector<double>> Graph::Johnson() const
{

	Graph copy_of_this(*this);
	std::unordered_map<Vertex*, int> VertexToIndex;
	int Vertecies_size = copy_of_this.m_vertices.size();
	std::vector<std::vector<double>> dist;
	Vertex* Vertex_s = copy_of_this.AddVertex("_Johnson_Source_");
	Vertecies_size++;
	for (int i = 0; i < Vertecies_size; i++)
	{
		VertexToIndex[copy_of_this.m_vertices[i]] = i;
		if (copy_of_this.m_vertices[i] != Vertex_s)
		{
			//Connect s to any other vertex with the weigth 0
			copy_of_this.ConnectEdge(Vertex_s->name, copy_of_this.m_vertices[i]->name, 0);
		}
	}

	std::vector<double> h = copy_of_this.Bellman_Ford(Vertex_s->name);
	if (h.size() <= 1)//cycle detected!
	{
		return {};
	}
	
	int idx_s = VertexToIndex[Vertex_s];  
	copy_of_this.RemoveVertex(Vertex_s->name);
	h.erase(h.begin() + idx_s);

	VertexToIndex.clear();
	for (int i = 0; i < copy_of_this.m_vertices.size(); i++)
	{
		VertexToIndex[copy_of_this.m_vertices[i]] = i;
	}


	Vertecies_size--;
	for (Vertex* u_vertex : copy_of_this.m_vertices)
	{
		for (auto& u_out_pair : u_vertex->Out_Edges)
		{
			int u = VertexToIndex[u_vertex];
			int v = VertexToIndex[u_out_pair.first];
			u_out_pair.second = u_out_pair.second + h[u] - h[v];
			//| w'(u, v) = w(u, v) + h(u) - h(v) | -> [ new weigth ]
		}
	}
	dist.reserve(Vertecies_size);

	for (int u = 0; u < Vertecies_size; u++)
	{
		Vertex* u_vertex = copy_of_this.m_vertices[u];
		dist.push_back(copy_of_this.Dijastra(u_vertex->name));
		for (int v = 0; v < Vertecies_size; v++)
		{
			if (dist[u][v] == INFINITE_DISTANCE)continue;
			dist[u][v] = dist[u][v] - h[u] + h[v];
			//| d(u, v) = d'(u, v) - h(u) + h(v) | -> [ jump to the old dist ]
		}
	}

	
	return dist;
}



std::ostream& operator<<(std::ostream& out, const Graph& g)
{
	out << "Graph adjacency list:\n";
	for (auto v : g.m_vertices)
	{
		out << "  " << v->name << " -> ";
		if (v->Out_Edges.empty())
		{
			out << "(no outgoing edges)";
		}
		else
		{
			bool first = true;
			for (const auto& out_pair : v->Out_Edges)
			{
				if (!first) out << " , ";
				first = false;
				out << out_pair.first->name << "(" << out_pair.second << ")";
			}
		}
		out << '\n';
	}
	return out;
}
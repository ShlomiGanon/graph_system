#include "Graph.h"


Graph::Graph(bool both_ways, int initial_size) : both_ways_graph(both_ways)
{ 
	m_vertices.reserve(initial_size);
}

Graph::Graph(const std::vector<std::string>& vertecies_names , bool both_ways) : both_ways_graph(both_ways)
{
	m_vertices.reserve(vertecies_names.size());
	for (auto name : vertecies_names)
	{
		Vertex* v = new Vertex(name);
		m_vertices.push_back(v);
	}
}

Graph::Graph(const Graph& other) : both_ways_graph(other.both_ways_graph)
{
	m_vertices.reserve(other.m_vertices.size());
	std::unordered_map<std::string, Vertex*> name_to_ptr(other.m_vertices.size());
	for (const Vertex* v : other.m_vertices)
	{
		Vertex* vertex = new Vertex(v->name);
		if (vertex == nullptr)
		{
			for (int i = 0; i < m_vertices.size(); i++)
			{
				delete m_vertices[i];
			}
			throw Graph_Exception("Memory allocation fail!");
		}
		else
		{
			name_to_ptr[vertex->name] = vertex;
			m_vertices.push_back(vertex);
		}
	}

	for (const Vertex* v : other.m_vertices)
	{
		Vertex* vertex = name_to_ptr[v->name];
		for (auto v_of_out : v->Out_Edges)
		{
			vertex->AddEdgeTo(*name_to_ptr[v_of_out.first->name], v_of_out.second);
		}
	}
}

Graph::Graph(Graph&& other) : m_vertices(std::move(other.m_vertices)), both_ways_graph(other.both_ways_graph)
{
	other.both_ways_graph = false;
	//other will be cleared
}

Graph::~Graph()
{
	std::cout << "~Graph" << std::endl;
	for (auto v : m_vertices)
	{
		delete v;
	}
}

Vertex* Graph::AddVertex(const std::string& vertex_Name)
{
	Vertex* vertex = GetVertex(vertex_Name);
	if (vertex != nullptr)return vertex;
	vertex = new Vertex(vertex_Name);
	m_vertices.push_back(vertex);
	return vertex;
}

Vertex* Graph::GetVertex(const std::string& vertex_name)const
{
	for (auto v : m_vertices)
	{
		if (v->name == vertex_name)return v;
	}
	return nullptr;
}

bool Graph::ConnectEdge(const std::string& from_vertex, const std::string& to_vertex, double weight)
{
	Vertex* from_v = AddVertex(from_vertex);
	Vertex* to_v = AddVertex(to_vertex);
	if (from_v != to_v)
	{
		from_v->AddEdgeTo(*to_v, weight, both_ways_graph);
		return true;
	}
	return false;
}

bool Graph::ConnectEdges(const std::string& from_vertex, std::vector<std::string> to_vertecies, double weigth)
{
	if (to_vertecies.size() == 0)return false;
	for (const auto& name : to_vertecies)
	{
		if (!ConnectEdge(from_vertex, name, weigth))return false;
	}
	return true;
}

bool Graph::RemoveVertex(const std::string& vertex_name)
{
	if (vertex_name.length() > 0)
	{
		for (auto it = m_vertices.begin(); it != m_vertices.end(); ++it)
		{
			auto v = *it;
			if (v->name == vertex_name)
			{
				delete v;
				m_vertices.erase(it);
				return true;
			}
		}
	}
	return false;
}

std::vector<std::string> Graph::getVertexNames() const
{
	std::vector<std::string> names;
	names.reserve(m_vertices.size());
	for (auto v : m_vertices)
	{
		names.push_back(v->name);
	}
	return names;
}

std::size_t Graph::size() const
{
	return m_vertices.size();
}

std::vector<Edge> Graph::getEdges() const
{
	std::vector<Edge> Edges;
	Edges.reserve(m_vertices.size()); // can bee (v^2 - v)
	for (auto v : m_vertices)
	{
		for (auto out_v : v->Out_Edges)
		{
			Edges.emplace_back(v, out_v.first, out_v.second);
		}
	}
	return Edges;
}

const std::vector<Vertex*>& Graph::getVertices() const
{
	return m_vertices;
}

bool Graph::IsBothWays() const
{
	return both_ways_graph;
}

std::unordered_map<Vertex*,int> Graph::BFS_distances(const std::string& source)
{
	Vertex* start_vertex = GetVertex(source);
	if (start_vertex == nullptr)throw Graph_Exception("vertex not found: " + source);
	std::unordered_map<Vertex*,int> dists;
	dists.reserve(m_vertices.size());
	for (auto v : m_vertices)
	{
		dists[v] = INFINITE_DISTANCE;
	}
	dists[start_vertex] = 0;

	std::unordered_map<Vertex*,bool> visited;
	std::vector<Vertex*> vertecies;
	vertecies.reserve(m_vertices.size());
	visited.reserve(m_vertices.size());
	vertecies.insert(vertecies.begin(), start_vertex);
	while (!vertecies.empty())
	{
		Vertex* current = vertecies.back();
		vertecies.pop_back();
		if (visited.find(current) == visited.end())//not have been visited!
		{
			//std::cout << current->name << " - dist " << dists[current] << std::endl;
			visited[current] = true;//set as visited

			for (auto out_pair : current->Out_Edges)
			{
				auto neighbor = out_pair.first;
				if (visited.find(neighbor) == visited.end())//neighbor is not have been visited
				{
					vertecies.push_back(neighbor);
					dists[neighbor] = dists[current] + 1;
				}
			}
		}
	}
	return dists;
}

std::vector<Vertex*> Graph::DFS(const std::string& source)
{
	Vertex* start_vertex = GetVertex(source);
	if (start_vertex == nullptr)throw Graph_Exception("vertex not found: " + source);
	std::vector<Vertex*> vertecies;
	vertecies.reserve(m_vertices.size());

	std::unordered_map<Vertex*, bool> visited;
	visited.reserve(m_vertices.size());

	std::stack<Vertex*> vertecies_stack;

	vertecies_stack.push(start_vertex);

	while (!vertecies_stack.empty())
	{
		Vertex* node = vertecies_stack.top();
		vertecies_stack.pop();
		if (visited.find(node) == visited.end())//if node is not has been visited
		{
			vertecies.push_back(node);
			std::cout << node->name << " , ";
			visited[node] = true;//mark node as visited

			for (auto out_pair : node->Out_Edges)
			{
				auto neighbor = out_pair.first;
				if(visited.find(neighbor) == visited.end())	vertecies_stack.push(neighbor);
			}
		}
	}
	std::cout << std::endl;

	return vertecies;
}

void Graph::DFS_recursive(const std::string start)
{
	Vertex* start_vertex = GetVertex(start);
	if (start_vertex == nullptr)throw Graph_Exception("vertex not found: " + start);
	std::unordered_map<Vertex*, bool> visited;
	std::stack<Vertex*> recStack;
	std::cout << "DFS-recursive: ";
	DFS_rec(start_vertex, visited, recStack);
	std::cout << std::endl;
}

void Graph::DFS_rec(Vertex* current, std::unordered_map<Vertex*, bool>& visited, std::stack<Vertex*>& recStack)
{
	std::stack<Vertex*> helper_stack;

	//check if current is in the stack
	while (!recStack.empty())
	{
		Vertex* v = recStack.top();
		recStack.pop();
		if (v == current)
		{
			std::cout << "Cycle detected at " << v->name << std::endl;
			return;
		}
		helper_stack.push(v);
	}
	while (!helper_stack.empty())
	{
		Vertex* v = helper_stack.top();
		helper_stack.pop();
		recStack.push(v);
	}

	if (visited.find(current) != visited.end())return;//check if was visited alredy

	std::cout << current->name << " ";

	visited[current] = true; // mark as visited

	recStack.push(current);

	for (auto out_pair : current->Out_Edges)
	{
		auto neighbor = out_pair.first;
		DFS_rec(neighbor, visited, recStack);
	}

	//remove current from recStack
	while (!recStack.empty())
	{
		Vertex* v = recStack.top();
		recStack.pop();
		if (v != current)
		{
			helper_stack.push(v);
		}
	}
	while (!helper_stack.empty())
	{
		Vertex* v = helper_stack.top();
		helper_stack.pop();
		recStack.push(v);
	}
}

std::vector<Vertex*> Graph::TopologicalSort_DFS()
{
	
	std::unordered_map<Vertex*, bool> visited;
	std::stack<Vertex*> order;
	for (auto node : m_vertices)
	{
		if (visited.find(node) == visited.end())//if node not have been visited
		{
			TopologicalSort_DFS_func(node, visited, order);
		}
	}
	std::vector<Vertex*> result;
	result.reserve(order.size());
	while (!order.empty())
	{
		result.push_back(order.top());
		order.pop();
	}
	return result;
}

int Graph::getConnectedComponents() const
{
	if (both_ways_graph == false)return -1;//Cannot calculate connected components on directed graph!
	int result = 0;
	std::unordered_map<Vertex*, bool> visited;
	std::stack<Vertex*> order;
	for (auto node : m_vertices)
	{
		if (visited.find(node) == visited.end())//if node not have been visited
		{
			TopologicalSort_DFS_func(node, visited, order);
			result++;
		}
	}
	return result;
}

void Graph::TopologicalSort_DFS_func(Vertex* node, std::unordered_map<Vertex*, bool>& visited, std::stack<Vertex*>& order)const
{
	visited[node] = true; // mark as visited

	for (auto out_pair : node->Out_Edges)
	{
		auto neighbor = out_pair.first;
		if (!visited[neighbor]) // neighbor was not visited
		{
			TopologicalSort_DFS_func(neighbor, visited, order);
		}
	}
	order.push(node);
}

std::list<Vertex*> Graph::TopologicalSort()
{
	std::list<Vertex*> order;
	std::unordered_map<Vertex*, int> inDegree;
	inDegree.reserve(m_vertices.size());
	for (auto v : m_vertices)
	{
		inDegree[v] = 0;
	}

	for (auto node : m_vertices)
	{
		for (auto out_pair : node->Out_Edges)
		{
			auto neighbor = out_pair.first;
			inDegree[neighbor] += 1;
		}
	}

	std::queue<Vertex*> Q;

	for (auto v : m_vertices)
	{
		if (inDegree.find(v) != inDegree.end() && inDegree[v] == 0)
		{
			Q.push(v);
			//Push (to end of Q) start vertex
		}
	}

	while (!Q.empty())
	{
		Vertex* current = Q.front();
		Q.pop();
		order.push_back(current);

		for (auto out_pair : current->Out_Edges)
		{
			auto neighbor = out_pair.first;
			inDegree[neighbor] -= 1;
			if (inDegree[neighbor] == 0)
			{
				Q.push(neighbor);
			}

		}

	}

	if (order.size() < m_vertices.size())std::cout << "Cycle detected" << std::endl;

	return order;
}

std::vector<std::pair<Vertex*, std::vector<std::pair<Vertex*, double>>::iterator>> Graph::GetPath_BFS(const std::string& FROM, const std::string& TO) const
{
	Vertex* from_vertex = GetVertex(FROM);
	Vertex* to_vertex = GetVertex(TO);
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
		throw Graph_Exception(ss.str());
		return {};
	}

	if (from_vertex == to_vertex)
	{
		for (auto it = from_vertex->Out_Edges.begin(); it != from_vertex->Out_Edges.end(); ++it)
		{
			if (it->first == from_vertex)
				return { { from_vertex , it } };
		}
		return {};
	}

	std::queue<Vertex*> Q;
	std::unordered_map<Vertex*, bool> visited;
	std::unordered_map < Vertex*, std::pair<Vertex*, std::vector<std::pair<Vertex*, double>>::iterator>> parent;
	Q.push(from_vertex);
	visited[from_vertex] = true;
	bool found = false;
	while (!Q.empty() && !found)
	{
		Vertex* current = Q.front();
		Q.pop();
		for (auto it = current->Out_Edges.begin(); it != current->Out_Edges.end(); it++)
		{
			Vertex* neigthbor = (*it).first;
			if (!visited[neigthbor])
			{
				parent[neigthbor] = std::make_pair(current, it);
				if (neigthbor == to_vertex)
				{ 
					found = true;
					break;
				}
				Q.push(neigthbor);
				visited[neigthbor] = true;
			}
		}
	}

	if (!found)
	{
		//cant find any path to the vertecies
		return {};
	}
	
	std::stack<Vertex*> reverse_path;
	Vertex* temp = to_vertex;
	while (temp != from_vertex)
	{
		reverse_path.push(temp);
		temp = parent[temp].first;
	}

	int vertecis = reverse_path.size();
	std::vector<std::pair<Vertex*, std::vector<std::pair<Vertex*, double>>::iterator>> result(vertecis);

	for (int i = 0; i < vertecis; i++)
	{
		Vertex* v = reverse_path.top() , *v_perant = parent[v].first;
		std::vector<std::pair<Vertex*, double>>::iterator it = parent[v].second;
		result[i] = { v_perant  , it};
		reverse_path.pop();
	}
	return result;
}



std::vector<double> Graph::Dijastra(const std::string& source)const
{
	std::unordered_map<const Vertex*, double> dist;     // Distance from source to each vertex
	std::unordered_map<const Vertex*, bool> visited;    // Marks if a vertex was visited

	Vertex* source_vertex = GetVertex(source);          // Get pointer to source vertex
	if (source_vertex == nullptr)
		throw Graph_Exception("cant find the vertex: " + source); // Throw if source not found

	dist.reserve(m_vertices.size());                    // Reserve space in the distance map
	for (auto v : m_vertices)
	{
		dist[v] = (v == source_vertex) ? 0 : INFINITE_DISTANCE; // 0 for source, ∞ for others
		visited[v] = false;                                   // Mark all vertices as unvisited
	}

	while (true) // Main loop
	{
		const Vertex* u = nullptr;

		// Find the unvisited vertex with the smallest known distance
		for (auto dist_pair : dist)
		{
			const Vertex* v = dist_pair.first;
			if (!visited[v])
			{
				if (u == nullptr || dist_pair.second < dist[u])
					u = v;
			}
		}

		// If no such vertex exists or it's unreachable, stop
		if (u == nullptr || dist[u] == INFINITE_DISTANCE)
			break;

		visited[u] = true; // Mark the vertex as visited

		// Check all outgoing edges from u
		for (auto u_out_pair : u->Out_Edges)
		{
			const Vertex* v = u_out_pair.first;
			double newDist = dist[u] + u_out_pair.second;

			// If u is unreachable, skip updating
			if (dist[u] == INFINITE_DISTANCE)
				newDist = INFINITE_DISTANCE;

			// Update distance if a shorter path is found
			if (newDist < dist[v])
				dist[v] = newDist;
		}
	}

	// Convert distance map to vector
	std::vector<double> dist_vec;
	dist_vec.reserve(size());
	for (auto v : m_vertices)
	{
		dist_vec.push_back(dist[v]);
	}

	return dist_vec;
}

std::vector<double> Graph::Bellman_Ford(const std::string& source)const
{
	Vertex* source_vertex = GetVertex(source);          // Get pointer to source vertex
	if (source_vertex == nullptr)throw Graph_Exception("cant find the vertex: " + source); // Throw if source not found

	std::unordered_map<const Vertex*, double> dist;
	auto Edges = this->getEdges();

	//Step 1 : initialization
	for (auto v : m_vertices)
	{
		dist[v] = INFINITE_DISTANCE;
	}
	dist[source_vertex] = 0;


	//Step 2: relaxetion on any edge (times the vertecies)
	int size_V = m_vertices.size();
	for (int i = 1; i < size_V; i++)
	{
		const Vertex* u, * v;
		double w , newDist;
		for (const auto& e : Edges)
		{
			//edge from u to v (with weight w)
			u = e.from;
			v = e.to;
			w = e.weight;
			newDist = dist[u] + w;
			if (dist[u] != INFINITE_DISTANCE && dist[v] > newDist)
			{
				dist[v] = newDist;
			}
		}
	}

	//Step 3: negative cycle check
	for (const auto& e : Edges)
	{
		//edge from u to v (with weight w)
		const Vertex* u = e.from;
		const Vertex* v = e.to;
		double w = e.weight;
		double newDist = dist[u] + w;
		if (dist[u] != INFINITE_DISTANCE && dist[v] > newDist)
		{
			//std::cout << "negative cycle -> from: " << u->name << ", to: " << v->name << std::endl;
			return { 0.0 };
		}
	}
	std::vector<double> dist_vec;
	dist_vec.reserve(size_V);

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
	std::vector<double> INF_VEC(V_size, INFINITE_DISTANCE);
	std::vector<std::vector<double>> dist(V_size, INF_VEC);;
	for (int i = 0; i < V_size; i++)
	{
		vToIndex[m_vertices[i]] = i;
	}
	//Step 1: initialization
	for (int i = 0; i < V_size; i++)
	{
		for (auto& i_out_pair : m_vertices[i]->Out_Edges)
		{
			Vertex* u = i_out_pair.first;
			dist[i][vToIndex[u]] = i_out_pair.second;
		}
		dist[i][i] = 0;
	}

	//Step 2
	for (int k = 0; k < V_size; k++)
	{
		for (int i = 0; i < V_size; i++)
		{
			for (int j = 0; j < V_size; j++)
			{
				if (dist[i][k] < INFINITE_DISTANCE && dist[k][j] < INFINITE_DISTANCE)
				{
					dist[i][j] = std::min(dist[i][j], dist[i][k] + dist[k][j]);
				}
			}
		}
	}
	//No check negative cycle!
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
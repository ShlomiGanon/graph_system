#pragma once
#include "Vertecie.h"
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <queue>
#include <list>


struct Graph_Exception : public std::runtime_error
{
	Graph_Exception(const char* msg) : runtime_error(msg)
	{

	}

	Graph_Exception(const std::string& msg) : runtime_error(msg)
	{

	}
};

class Graph
{
private:
	std::vector<Vertex*> m_vertices;
	bool both_ways_graph;

public:
	static constexpr int INFINITE_DISTANCE = std::numeric_limits<int>::max();
	//static constexpr int INFINITE_DISTANCE_int = std::numeric_limits<int>::max();
	Graph(bool both_ways = false ,int initial_size = 0);
	Graph(const std::vector<std::string>& vertices_names , bool both_ways = false);
	Graph(const Graph& other);
	Graph(Graph&& other);
	virtual ~Graph();

	
	Vertex* AddVertex(const std::string& vertex_name);
	Vertex* GetVertex(const std::string& vertex_name)const;
	bool ConnectEdge(const std::string& from_vertex, const std::string& to_vertex, double weigth);
	bool ConnectEdges(const std::string& from_vertex, std::vector<std::string> to_vertecies, double weigth);
	bool RemoveVertex(const std::string& vertex_name);
	//bool RemoveVertecie( std::list<Vertecie>::iterator it );
	std::vector <std::string> getVertexNames()const;
	std::size_t size()const;
	std::vector <Edge> getEdges()const;
	const std::vector <Vertex*>& getVertices()const;
	bool IsBothWays()const;

	std::unordered_map<Vertex*,int> BFS_distances(const std::string& source);
	std::vector<Vertex*> DFS(const std::string& source);
	void DFS_recursive(const std::string start);
	void DFS_rec(Vertex* current, std::unordered_map<Vertex*, bool>& visited, std::stack<Vertex*>& recStack);
	std::vector<Vertex*> TopologicalSort_DFS();
	int getConnectedComponents()const;
	void TopologicalSort_DFS_func(Vertex* node, std::unordered_map<Vertex*, bool>& visited, std::stack<Vertex*>& order)const;

	std::list<Vertex*> TopologicalSort();

	std::vector<std::pair<Vertex*, std::vector<std::pair<Vertex*, double>>::iterator>> GetPath_BFS(const std::string& FROM, const std::string& TO) const;

	std::vector<double> Dijastra(const std::string& source)const;
	std::vector<double> Bellman_Ford(const std::string& source)const;
	std::vector<std::vector<double>> Floyd_Warshall()const;
	std::vector< std::vector<double>> Johnson()const;

	friend std::ostream& operator<<(std::ostream& out, const Graph& g);
};


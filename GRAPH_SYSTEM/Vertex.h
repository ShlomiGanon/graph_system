#pragma once
#include "ProjectIncludes.h"
class Vertex;
struct Edge
{
	const Vertex* from;
	const Vertex* to;
	double weight;
	Edge(const Vertex* FROM, const Vertex* TO, double WEIGHT) :from(FROM), to(TO), weight(WEIGHT) {}
};


class Vertex
{
public:
	std::string name;
	std::vector<std::pair<Vertex*,double>> Out_Edges, In_Edges;
public:
	Vertex(const std::string& name, int initial_in_edges = 0, int initial_out_edges = 0);
	virtual ~Vertex();
	bool AddEdgeTo(Vertex& v, double weight , bool both_directions = false);
	bool RemoveEdgeTo(Vertex& v, bool both_directions = false);
};


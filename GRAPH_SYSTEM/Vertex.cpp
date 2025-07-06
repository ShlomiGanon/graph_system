#include "Vertex.h"
Vertex::Vertex(const std::string& name , int initial_in_edges , int initial_out_edges) : name(name)
{
	In_Edges.reserve(initial_in_edges);
	Out_Edges.reserve(initial_out_edges);
}
Vertex::~Vertex()
{
	for (const auto& pair : Out_Edges)
	{
		Vertex* v = pair.first;
		for (auto it = v->In_Edges.begin(); it != v->In_Edges.end(); it++)
		{
			if (it->first == this)
			{
				v->In_Edges.erase(it);
				break;
			}
		}
	}
	for (const auto& pair : In_Edges)
	{
		Vertex* v = pair.first;
		for (auto it = v->Out_Edges.begin(); it != v->Out_Edges.end(); it++)
		{
			if (it->first == this)
			{
				v->Out_Edges.erase(it);
				break;
			}
		}
	}
}
bool Vertex::AddEdgeTo(Vertex& v, double weight, bool both_directions)
{
	bool out_edge_exists = false;
	// Check if edge already exists in Out_Edges and update weight
	for (auto& out_edge_of_this : Out_Edges)
	{
		if (out_edge_of_this.first == &v)
		{
			//out_edge_of_this.weight = weight; // update weight if exists
			out_edge_exists = true;
			break;
		}
	}
	if (!out_edge_exists)
		Out_Edges.emplace_back(&v, weight); // add new outgoing edge

	bool in_edge_exists = false;
	// Check if incoming edge exists in v.In_Edges and update weight
	for (auto& in_edge_of_v : v.In_Edges)
	{
		if (in_edge_of_v.first == this)
		{
			//in_edge_of_v.weight = weight; // update weight if exists
			in_edge_exists = true;
			break;
		}
	}
	if (!in_edge_exists)
		v.In_Edges.emplace_back(this, weight); // add new incoming edge

	// If bidirectional, connect back without infinite recursion
	if (both_directions)v.AddEdgeTo(*this, weight, false);

	return out_edge_exists || in_edge_exists;
}


bool Vertex::RemoveEdgeTo(Vertex& v, bool both_directions)
{
	bool Found = false;
	// Find and remove edge to v
	for (auto it = this->Out_Edges.begin(); it != this->Out_Edges.end(); it++)
	{
		if (it->first == &v)
		{
			this->Out_Edges.erase(it);

			// Remove corresponding incoming edge from v
			for (auto v_it = v.In_Edges.begin(); v_it != v.In_Edges.end(); v_it++)
			{
				if (v_it->first == this)
				{
					v.In_Edges.erase(v_it);
					Found = true;
					break;
				}
			}
			break;
		}
	}
	// If bidirectional, also disconnect reverse edge and combine results
	if (both_directions)
		return v.RemoveEdgeTo(*this, false) || Found;

	// Return true if edge removed in this direction
	return Found;
}


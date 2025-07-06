#include "Graph.h"
#include <iostream>

static int JampArray(const std::vector<int>& arr);

int main()
{
	
	Graph g;
	g.ConnectEdge("1", "2", 6);
	g.ConnectEdge("1", "3", 5);
	g.ConnectEdge("1", "4", 5);
	g.ConnectEdge("4", "3", 2);
	g.ConnectEdge("4", "6", -1);
	g.ConnectEdge("3", "2", -2);
	g.ConnectEdge("3", "5", 1);
	g.ConnectEdge("2", "5", -1);
	g.ConnectEdge("5", "7", 3);
	g.ConnectEdge("6", "7", 3);
	std::cout << g << std::endl;

	std::string FROM, TO;
	std::cout << "please enter source vertex: ";
	std::cin >> FROM;
	std::cout << "please enter destanation vertex: ";
	std::cin >> TO;
	std::vector<std::pair<Vertex*, std::vector<std::pair<Vertex*, double>>::iterator>> v;

	try
	{
		v = g.GetPath_BFS(FROM, TO);
	}
	catch (Graph_Exception& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		return 1;
	}

	if (v.empty())
	{
		std::cout << "No path found." << std::endl;
	}
	else
	{
		std::cout << "Path from " << FROM << " to " << TO << ":\n";
		std::cout << v.front().first->name;
		int sum = 0;
		for (const auto& step : v)
		{
			int weight = (*step.second).second;
			std::cout << " ==(" << weight << ")==> " << (*step.second).first->name;
			sum += weight;
		}
		std::cout << "\nTotal weight: " << sum << "\n";
	}

	std::cout << "Dijastra: " << std::endl;
	for (auto u : g.Dijastra("1"))
	{
		if(u == Graph::INFINITE_DISTANCE)std::cout <<  "inf, ";
		else std::cout << u << ", ";
	}
	std::cout << std::endl;
	return 0;
}


static int JampArray(const std::vector<int>& arr)
{
	int size = arr.size();
	std::vector<int> BestJumps(size,0);

	for (int i = size - 2; i >= 0; i--)// O(n^2)
	{
		if (arr[i] <= 0)
		{
			BestJumps[i] = 999999;//invalid [only positive numbers]
			continue;
		}
		int jump_to = i + arr[i];
		if (jump_to > size)BestJumps[i] = 1;
		else
		{
			int minimal_index = i + 1;
			for (int j = i + 2; j <= jump_to && j < size; j++)
			{
				if (BestJumps[minimal_index] > BestJumps[j])minimal_index = j;
			}
			BestJumps[i] = BestJumps[minimal_index] + 1;
		}
	}
	return BestJumps[0];
	
}
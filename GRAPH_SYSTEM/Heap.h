#pragma once
#include <vector>
#include <unordered_map>
#include <utility> 

enum HeapType { MIN_HEAP, MAX_HEAP };

template <typename KEY_TYPE, typename VALUE_TYPE>
struct Heap_Node
{
	KEY_TYPE key;        
	VALUE_TYPE value;   
};

template <typename KEY_TYPE, typename VALUE_TYPE = KEY_TYPE>
class Heap
{
private:
	HeapType heap_type;
	std::vector<Heap_Node<KEY_TYPE, VALUE_TYPE>> heap;
	std::unordered_map<KEY_TYPE, int> positions;

	bool compareByHeapRule(const VALUE_TYPE& a, const VALUE_TYPE& b)const
	{
		switch (heap_type)
		{
			case MIN_HEAP: return a < b;
			case MAX_HEAP: return a > b;
		}
		return false;// not reachable
	}

	// WARNING: Assumes index > 0.
	// No bounds check is performed to avoid unnecessary overhead.
	// Calling this with index == 0 will result in undefined behavior.
	static constexpr int parent(int index)
	{
		return (index - 1) / 2;
	}

	// WARNING: This function does not check if the returned child index is within bounds.
	// Caller must ensure index is valid and within the heap size.
	static constexpr int left_child(int index)
	{
		return 2 * index + 1;
	}

	// WARNING: This function does not check if the returned child index is within bounds.
	// Caller must ensure index is valid and within the heap size.
	static constexpr int right_child(int index)
	{
		return 2 * index + 2;
	}

	// WARNING: Assumes both indices are within bounds of the heap vector.
	// No index validation is done for performance reasons.
	// Caller must ensure both indices are valid before calling.
	void swap(int index1, int index2)
	{
		positions[heap[index1].key] = index2;
		positions[heap[index2].key] = index1;

		Heap_Node<KEY_TYPE,VALUE_TYPE> temp = heap[index1];
		heap[index1] = heap[index2];
		heap[index2] = temp;
	}
	
	// Moves the element at the given index up the heap to restore heap order.
	// WARNING: Assumes that 'index' is a valid index within the heap vector.
	// No boundary checks are performed for performance reasons.
	// This function should only be called from trusted internal methods.
	void heapifyUp(int index)
	{
		if (index == 0) return;  // reached root, no parent
		int p = parent(index);
		if (compareByHeapRule(heap[index].value, heap[p].value))
		{
			swap(index, p);
			heapifyUp(p);
		}
	}

	// Moves the element at the given index down the heap to restore heap order.
	// WARNING: Assumes that 'index' is a valid index within the heap vector.
	// This function relies on 'compareByHeapRule' and 'swap' being safe for the current heap state.
	// No checks are done for out-of-bounds access — internal use only.
	void heapifyDown(int index)
	{
		int size_ = size();
		int currentIndex = index;

		while (true)
		{
			int leftIndex = left_child(currentIndex);
			int rightIndex = right_child(currentIndex);
			int bestIndex = currentIndex;

			if (leftIndex < size_ && compareByHeapRule(heap[leftIndex].value, heap[bestIndex].value))
				bestIndex = leftIndex;
			if (rightIndex < size_ && compareByHeapRule(heap[rightIndex].value, heap[bestIndex].value))
				bestIndex = rightIndex;

			if (bestIndex == currentIndex)
				break;

			swap(currentIndex, bestIndex);
			currentIndex = bestIndex;
		}
	}

	void heapify()
	{
		for (int i = parent(size() - 1); i >= 0; i--)
		{
			heapifyDown(i);
		}
	}

public:
	Heap(HeapType heap_type, size_t initial_space = 1) : heap_type(heap_type)
	{
		heap.reserve(initial_space);
		positions.reserve(initial_space);
	}

	Heap(const Heap& other) : heap_type(other.heap_type), heap(other.heap), positions(other.positions)
	{

	}

	Heap(Heap&& other) noexcept
		: heap_type(other.heap_type), heap(std::move(other.heap)), positions(std::move(other.positions))
	{

	}

	bool empty() const
	{
		return heap.empty();
	}

	size_t size()const
	{
		return heap.size();
	}

	// Inserts a new element with both key and value.
	// Returns false if the key already exists in the heap.
	bool insert(const KEY_TYPE& key, const VALUE_TYPE& value)
	{
		if (positions.find(key) == positions.end()) // key is not in the heap
		{
			heap.push_back({ key, value });
			int index = static_cast<int>(heap.size() - 1);
			positions[key] = index;
			heapifyUp(index);
			return true;
		}
		return false; // key already exists
	}

	// Inserts a new element with only the key.
	// Assumes VALUE_TYPE has a constructor that accepts KEY_TYPE.
	// Calls the main insert method, creating the value from the key.
	bool insert(const KEY_TYPE& key)
	{
		return insert(key, VALUE_TYPE(key));
	}

	bool changeValue(const KEY_TYPE& key, const VALUE_TYPE& new_value)
	{
		auto it = positions.find(key);
		if (it != positions.end())
		{
			int index = it->second;
			heap[index].value = new_value;
			if (index > 0 && compareByHeapRule(new_value, heap[parent(index)].value))
				heapifyUp(index);
			else
				heapifyDown(index);
			return true;
		}
		return false;
	}

	bool pop_Root()
	{
		if (empty())return false;

		auto positions_it = positions.find(heap[0].key);// Find the iterator to the root's key in the positions map

		swap(0, size() - 1);// Swap the root node with the last node in the heap

		positions.erase(positions_it);// Remove the root's key from the positions map

		heap.pop_back();// Remove the last element (which is the original root)

		if (!empty())heapifyDown(0);// Restore heap property by heapifying down from the root

		return true;
	}

	std::pair<KEY_TYPE, VALUE_TYPE> peek_Root()
	{
		if(empty())throw std::runtime_error("Heap is empty, cannot peek root.");
		return { heap[0].key , heap[0].value };
	}

	Heap cloneAs(HeapType new_type) const
	{
		Heap result(*this);
		if (new_type != result.heap_type)
		{
			result.heap_type = new_type;
			result.heapify();
		}
		return result;
	}
};


/*
 * collision_detection.cpp
 *
 *  Created on: Nov 22, 2022
 *      Author: norms
 */

#include "collision_detection.h"
#include "vectors.h"
#include "rigid_bodies.h"
#include <cmath>
#include <string.h>

// Functions dealing with the broad phase of collision detection

// Create a fat box given a circular rigid body
AABB create_fatbox(Circular_Rigid_Body *body)
{
	// Enlarge box by one unit of velocity in each direction. v * 1s
	AABB fatbox;
	// Create the tight fitting box around the circle
	fatbox.min = {body->pos.x - body->radius, body->pos.y - body->radius};
	fatbox.max = {body->pos.x + body->radius, body->pos.y + body->radius};
	// Enlarge the box only in the direction of the velocity of the object
	// X direction
	if (body->linear_vel.x > 0)
	{
		// If positive add it to the maximum
		fatbox.max.x += body->linear_vel.x;
	}
	else
	{
		// If negative, add it to the minimum
		fatbox.min.x += body->linear_vel.x;
	}
	// Y direction
	if (body->linear_vel.y > 0)
	{
		fatbox.max.y += body->linear_vel.y;
	}
	else
	{
		fatbox.min.y += body->linear_vel.y;
	}
	return fatbox;
}
// Union of two AABBs
AABB union_of_aabbs(AABB a, AABB b)
{
	// Create a single AABB surrounding two AABBs
	AABB result;

	// Use minimum and maximum vector functions.
	result.min = minimum_vector(a.min, b.min);
	result.max = maximum_vector(a.max, b.max);

	return result;
}

// Surface Area (Area) of AABB
float AABB::area(void)
{
	Vector d = {max.x - min.x, max.y - min.y};
	return d.x * d.y;
}

// Node functions
bool Node::check_if_leaf(void)
{
	if (child0 == -1)
	{
		return true;
	}
	return false;
}

// AABB Tree functions

// Instantiates the tree:
void Tree::instantiate_tree(void)
{
	node_capacity = 16;
	node_list = (Node *)malloc(node_capacity * sizeof(Node));
	root_index = -1;
}

// Returns the index to a new node peeled off the array
int Tree::allocate_node(void)
{
	// Check if node array needs to be expanded
	if (node_count == node_capacity)
	{
		Node *old_list = node_list;
		// Double the capacity and reallocate to new array, copying over old data
		node_capacity *= 2;
		node_list = (Node *)malloc(node_capacity * sizeof(Node));
		memcpy(node_list, old_list, node_count * sizeof(Node));
		// Free the old list
		free(old_list);
	}

	// Take a node off the array at the index node_count
	int index = node_count;
	node_count++;
	return index;
}

// Allocate new leaf node and add to tree node list
int Tree::allocate_leaf_node(int object, AABB fatbox)
{
	// Get index of newly allocated node
	int leaf_index = allocate_node();

	// Set node as a leaf
	// Leaves have an AABB surrounding the object
	node_list[leaf_index].box = fatbox;
	// Link this node to the object in the scene
	node_list[leaf_index].object_index = object;
	// Leaves have no children
	node_list[leaf_index].child0 = -1;
	node_list[leaf_index].child1 = -1;
	// Initialize height to 0
	node_list[leaf_index].height = 0;
	// It is a leaf
	node_list[leaf_index].is_leaf = true;

	return leaf_index;
}

int Tree::allocate_branch_node(int child0, int child1)
{
	// Allocate space
	int branch_index = allocate_node();

	// Set as branch
	// Branches have two children
	node_list[branch_index].child0 = child0;
	node_list[branch_index].child1 = child1;

	// Create the union AABB
	node_list[branch_index].box = union_of_aabbs(node_list[child0].box, node_list[child1].box);
	// Set the branch as the parent of the two children
	node_list[child0].parent = branch_index;
	node_list[child1].parent = branch_index;
	// Not a leaf
	node_list[branch_index].is_leaf = false;

	return branch_index;
}

void Tree::insert_leaf(int object, AABB fatbox)
{
	// Allocate the leaf node
	int new_leaf = allocate_leaf_node(object, fatbox);
	// If tree has no root, assign this leaf to be the root. Root has parent index of null (-1)
	if (root_index == -1)
	{
		node_list[new_leaf].parent = -1;
		root_index = new_leaf;
		return;
	}

	// Stage 1: Look for the best sibling

	// Start at the root.
	int partner = root_index;
	// Initialize the best cost
	AABB leaf_aabb = node_list[new_leaf].box;

	while (!node_list[partner].is_leaf)
	{
		// Calculate the direct and inheritance costs
		float area = node_list[partner].box.area();

		AABB union_aabb = union_of_aabbs(leaf_aabb, node_list[partner].box);
		float union_area = union_aabb.area();

		// Direct cost of creating a parent node for both the partner and new leaf
		float cost = 2.0 * union_area;

		// Minimum cost of moving leaf down tree
		float inheritance_cost = 2.0 * (union_area - area);

		// Calculate costs of descending into children
		int child0 = node_list[partner].child0;
		int child1 = node_list[partner].child1;

		float cost0;
		if (node_list[child0].is_leaf)
		{
			AABB combined_aabb = union_of_aabbs(leaf_aabb, node_list[child0].box);
			cost0 = combined_aabb.area() + inheritance_cost;
		}
		else
		{
			float old_area = node_list[child0].box.area();
			AABB combined_aabb = union_of_aabbs(leaf_aabb, node_list[child0].box);
			float new_area = combined_aabb.area();
			cost0 = (new_area - old_area) + inheritance_cost;
		}

		float cost1;
		if (node_list[child1].is_leaf)
		{
			AABB combined_aabb = union_of_aabbs(leaf_aabb, node_list[child1].box);
			cost1 = combined_aabb.area() + inheritance_cost;
		}
		else
		{
			float old_area = node_list[child1].box.area();
			AABB combined_aabb = union_of_aabbs(leaf_aabb, node_list[child1].box);
			float new_area = combined_aabb.area();
			cost1 = (new_area - old_area) + inheritance_cost;
		}

		// If the cost to descend is greater than the current cost, then break out
		if ((cost < cost0) && (cost < cost1))
		{
			break;
		}

		// Descend into the child with the cheapest cost
		if (cost0 < cost1)
		{
			partner = child0;
		}
		else
		{
			partner = child1;
		}
	}

	int best_sibling = partner;

	// Stage 2: Insert the new leaf next to its sibling, and create the new parent.
	int old_parent = node_list[best_sibling].parent;
	int new_parent = allocate_branch_node(new_leaf, best_sibling);
	node_list[new_parent].parent = old_parent;
	// The only possible height contribution to the new parent would be the sibling.
	// Leaves have height 0.
	node_list[new_parent].height = node_list[best_sibling].height + 1;

	// Sibling is the root
	if (best_sibling == root_index)
	{
		// Set the parent to be the new root
		root_index = new_parent;
		node_list[new_parent].parent = -1;
	}
	// Sibling is not the root
	else
	{
		// If child0 is the sibling
		if (best_sibling == node_list[old_parent].child0)
		{
			node_list[old_parent].child0 = new_parent;
		}
		else
		{
			node_list[old_parent].child1 = new_parent;
		}
	}

	// Go back up the tree and refit the ancestor aabb's
	int index = node_list[best_sibling].parent;
	while (index != -1)
	{
		// Rotate the tree if imbalanced
		// Refit the aabb of the node
		int child0 = node_list[index].child0;
		int child1 = node_list[index].child1;
		node_list[index].box = union_of_aabbs(node_list[child0].box, node_list[child1].box);

	}
}

int Tree::balance_the_tree(int curr_node)
{
	// If the current node is a leaf. Or if both children are leaves, then nothing happens
	if (node_list[curr_node].is_leaf || node_list[curr_node].height < 2)
	{
		return curr_node;
	}

	// Determine the balance of the current_node
	int child0 = node_list[curr_node].child0;
	int child1 = node_list[curr_node].child1;

	int balance = node_list[child0].height - node_list[child1].height;

	// If it is left heavy, rotate right
	// Define left as child0, right as child1
	if (balance > 1)
	{
		// Get the children of child0
		int child00 = node_list[child0].child0;
		int child01 = node_list[child0].child1;
		// Get old parent of current node
		int old_parent = node_list[curr_node].parent;

		// Swap heavy child0 position with current node
		node_list[child0].child0 = curr_node;
		node_list[child0].parent = old_parent;
		node_list[curr_node].parent = child0;
		// Make the old parent point towards child0 as new child.
		// Current node wasnt old root
		if (old_parent != -1)
		{
			if (curr_node == node_list[old_parent].child0)
			{
				node_list[old_parent].child0 = child0;
			}
			else
			{
				node_list[old_parent].child1 = child0;
			}
		}
		else
		{
			root_index = child0;
		}

		// Perform the rotation
		// Determine which grandchild is the heaviest
		// If the left grandchild is the heaviest
		if (node_list[child00].height > node_list[child01].height)
		{
			// Do the swaps
			node_list[child0].child1 = child00;
			node_list[curr_node].child0 = child01;
			node_list[child01].parent = curr_node;
			// Calculate new aabb's for curr_node and child0
			node_list[curr_node].box = union_of_aabbs(node_list[child01].box, node_list[child1].box);
			node_list[child0].box = union_of_aabbs(node_list[curr_node].box, node_list[child00].box);
			// Calculate new heights


		}
		// Right grandchild is the heaviest
		else
		{
			node_list[child0].child1 = child01;
			node_list[curr_node].child0 = child00;
			node_list[child00].parent = curr_node;

			node_list[curr_node].box = union_of_aabbs(node_list[child00].box, node_list[child1].box);
			node_list[child0].box = union_of_aabbs(node_list[curr_node].box, node_list[child01].box);
		}
		return child0;
	}
	// Right heavy, need a left rotation
	if (balance < -1)
	{
		// Get the children of child1
		int child10 = node_list[child1].child0;
		int child11 = node_list[child1].child1;
		// Get old parent
		int old_parent = node_list[curr_node].parent;

		// Swap heavy child1 with current_node
		node_list[child1].child1 = curr_node;
		node_list[child1].parent = old_parent;
		node_list[curr_node].parent = child1;

		// Set old parent to point to child1 as child now, consider case of root
		if (old_parent != -1)
		{
			if (curr_node == node_list[old_parent].child0)
			{
				node_list[old_parent].child0 = child1;
			}
			else
			{
				node_list[old_parent].child1 = child1;
			}
		}
		else
		{
			root_index = child1;
		}

		// If the left grandchild is the heaviest
		//if node_list[child10].height
		// TODO



	}

}

// Traverses the tree recursively and finds the best partner for the new node
// Best partner minimizes the extra surface area added to the tree
//Node *Tree::find_best_partner(Node *new_leaf, Node *partner, float best_cost)
//{
//	// Write the end condition.
//	// If we have either reached a leaf, or found it is not worth exploring further, we return the pointer to the node.
//	if (partner->is_leaf)
//	{
//		return partner;
//	}
//	// If the partner is not a leaf, calculate if it is worth exploring further
//	// If the lower bound is higher than the best_cost, then stop
//	float lower_bound;
//
//}






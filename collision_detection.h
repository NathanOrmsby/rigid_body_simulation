/*
 * collision_detection.h
 *
 *  Created on: Nov 22, 2022
 *      Author: norms
 */

#ifndef COLLISION_DETECTION_H_
#define COLLISION_DETECTION_H_

#include "vectors.h"
#include <vector>

// Forward declarations
class Circular_Rigid_Body;

// Axis Aligned Bounding Boxes (AABBs). Rectangles that bound objects
class AABB
{
public:
	// Bottom left point
	Vector min;
	// Top right point
	Vector max;

	// AABB methods
	float area(void);
};

// Classes for the Binary collision tree
// Node on the binary tree. Made up of nodes and leaves
class Node
{
public:
	AABB box;
	// Index for the mass object in the mass list if this node is a leaf
	int object_index;
	// Parent index in the node list
	int parent;
	// Child indices in the node list.
	int child0;
	int child1;
	// Whether this node is a leaf or not
	bool is_leaf;
	// Height on the tree
	int height;

	// Class functions

	// Checks if node is a leaf. If false, then its a branch
	bool check_if_leaf(void);

};

// The tree containing all of the nodes
class Tree
{
public:
	// Array of nodes. These will be used to build the tree
	Node *node_list;
	// Initial capacity will be 16 nodes
	int node_capacity;
	int node_count;
	int root_index;

	// Create the tree
	void instantiate_tree(void);
	// Get index of usable node off list
	int allocate_node(void);
	// Allocate space for a leaf node and add it to the node list
	int allocate_leaf_node(int object_index, AABB box);
	// Allocate space for a branch node and add it to the node list
	int allocate_branch_node(int child0, int child1);

	// Inserting a leaf into the tree structure (HARD)
	void insert_leaf(int object_index, AABB box);

	// Balance the avl tree, rotations
	int balance_the_tree(int curr_node);
};

// Other functions relating to AABBs
AABB create_fatbox(Circular_Rigid_Body *body);
AABB union_of_aabbs(AABB a, AABB b);


#endif /* COLLISION_DETECTION_H_ */

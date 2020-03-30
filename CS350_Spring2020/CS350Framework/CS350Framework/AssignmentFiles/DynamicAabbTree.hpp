///////////////////////////////////////////////////////////////////////////////
///
/// Authors: Joshua Davis
/// Copyright 2015, DigiPen Institute of Technology
///
///////////////////////////////////////////////////////////////////////////////
#pragma once

#include "SpatialPartition.hpp"
#include "Shapes.hpp"
#include <unordered_map>
#include <stack>

/******Student:Assignment3******/
/// You must implement a dynamic aabb tree as we discussed in class.
class DynamicAabbTree : public SpatialPartition
{
public:
  DynamicAabbTree();
  ~DynamicAabbTree();

  // Spatial Partition Interface
  void InsertData(SpatialPartitionKey& key, SpatialPartitionData& data) override;
  void UpdateData(SpatialPartitionKey& key, SpatialPartitionData& data) override;
  void RemoveData(SpatialPartitionKey& key) override;

  void DebugDraw(int level, const Math::Matrix4& transform, const Vector4& color = Vector4(1), int bitMask = 0) override;

  void CastRay(const Ray& ray, CastResults& results) override;
  void CastFrustum(const Frustum& frustum, CastResults& results) override;

  void SelfQuery(QueryResults& results) override;

  void FilloutData(std::vector<SpatialPartitionQueryData>& results) const override;

  static const float mFatteningFactor;

  // Add your implementation here
private:
    class Node
    {
    public:
        //default ctor for Node
        Node ()
            :mAabb{ Aabb{} }, mClientData {nullptr}, 
            mLeft{ nullptr }, mRight{ nullptr }, mParent{nullptr}, 
            mHeight{ 0 }, mLastAxis{0}
        {} 
        Aabb mAabb;
        void* mClientData;
        Node* mLeft;
        Node* mRight;
        Node* mParent;
        int mHeight;
        size_t mLastAxis;
    };

    //helper functions =================
    // SelectNode selects which node to insert into, using Surface Area Heuristics 
    Node* SelectNode(Node* insertingNode, Node* n0, Node* n1);
    //Recursively updates the aabb and height of the parents from the leaf_node 
    void UpdateTreeParents(Node* leaf_node);
    //balances the tree
    void BalanceTree();
    //void SplitNode(Node* oldNode, Node* newNode);
    //PreOrder Traversal recursive function
    void PreOrderTraverse(Node* curr_node, std::vector<SpatialPartitionQueryData>& results, int depth) const;
    //Used to recurse down the tree in preorder traversal and calls DrawAabb to do the debug draw
    void PreOrderDraw(Node* curr_node, int level, const Math::Matrix4& transform, const Vector4& color, int bitMask, int depth);
    //Used to recurse down the tree to cast the ray 
    void CastRay_Rec(Node* curr_node, const Ray& ray, CastResults& results);
    //Used to recurse down the tree to cast the frustum
    void CastFrustum_Rec(Node* curr_node, const Frustum& frustum, CastResults& results);
    //Used in CastFrustum to recursively add all leaves node to the result
    void CastFrustum_ReturnLeavesNode_Rec(Node* curr_node, CastResults& results);
    //Used at the root node to recurse down the 3 cases, splits the recursive function into 2 functions
    void SelfQuery(Node* node, QueryResults& results);
    // second recursive function 
    void SelfQuery(Node* node0, Node* node1, QueryResults& results);
    //original split node function to test all children 
    void SplitNodes(Node* node0, Node* node1, QueryResults& results);



    //Members/Data: 
    std::unordered_map<unsigned int, Node*> mKeyNode_map; // used to map the key to the node
    std::vector<Node*> mMiddleNodes; // used to keep track all the middle ndoes of the tree
    std::stack<Node*> visted_nodes_stack; //used to balance the tree
    Node* mRoot = nullptr; // the root node of the tree
    unsigned int num_inserted = 0 ;
};

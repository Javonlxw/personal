///////////////////////////////////////////////////////////////////////////////
///
/// Authors: Joshua Davis
/// Copyright 2015, DigiPen Institute of Technology
///
///////////////////////////////////////////////////////////////////////////////
#include "Precompiled.hpp"
#include <iostream>

const float DynamicAabbTree::mFatteningFactor = 1.1f;
/******Student:Assignment3******/
DynamicAabbTree::DynamicAabbTree()
{
  mType = SpatialPartitionTypes::AabbTree;
}
/******Student:Assignment3******/
DynamicAabbTree::~DynamicAabbTree()
{
}

/*
    When inserting a node into the tree use the surface area heuristic described in class.
    After inserting you must re-balance the tree using avl-rotations. When you split a leaf node put
    the old leaf node as the left child and the new node (the one you are inserting) as the right child.
    Also make sure to fatten the new aabb’s half extent by the provided mFatteningFactor coefficient
    before inserting.
 */
void DynamicAabbTree::InsertData(SpatialPartitionKey& key, SpatialPartitionData& data)
{
  //Warn("Assignment3: Required function un-implemented");
    /******Student:Assignment3******/
    visted_nodes_stack.empty();

  //Fatten the data's aabb 

  Math::Vector3 fattenHalfExt = data.mAabb.GetHalfSize() * mFatteningFactor;
  Math::Vector3 center = data.mAabb.GetCenter();
  Aabb fat_aabb{ center - fattenHalfExt, center + fattenHalfExt  };

  //Create the right Node with the data
  Node* new_node = new Node();
  //Assign fatten aabb and client's data to the new node 
  new_node->mAabb = fat_aabb;
  new_node->mClientData = data.mClientData;
  //Store the node pointer and key to the key-node map in the spatial partition key 
  key.mUIntKey = num_inserted++;// O(1)
  //key.mVoidKey = (void*)new_node;
  //Map the newly created node with the spatial partition key 
  mKeyNode_map.emplace(std::make_pair(key.mUIntKey, new_node)); //allows for O(1) search time

  //Case 1: Empty Tree
  if (mRoot == nullptr)
  {
      mRoot = new_node;
      return;
  }
      

  //For Case 2 and 3, finding the old_node:
  Node* traversed_node = mRoot;
  Node* old_node = mRoot;
  visted_nodes_stack.push(traversed_node);
  traversed_node = SelectNode(new_node, traversed_node->mLeft, traversed_node->mRight);
  //traverse down the tree till the root node to find the 'old_node'
  while (traversed_node != nullptr)
  {
      visted_nodes_stack.push(traversed_node);
      old_node = traversed_node;
      traversed_node = SelectNode(new_node, traversed_node->mLeft, traversed_node->mRight);
  }
  //Create a new middle_node w/o clientData
  //Set left child to be old_node   
  //Set right child to be new_node
  Node* middle_node = new Node{};
  middle_node->mLeft = old_node;
  middle_node->mRight = new_node;
  //Update the vector of middle_nodes
  mMiddleNodes.push_back(middle_node);

  //Link back the old_node's parent to the middle node now
  Node* old_node_parent = old_node->mParent;

  //Check for case 2, where old node is the root node
  if (old_node_parent == nullptr)
  {
      //Make the root node the middle node
      mRoot = middle_node;
      mRoot->mParent = nullptr;
  }
  else // Case 3
  {
      if (old_node_parent->mLeft == old_node)
      {
          old_node_parent->mLeft = middle_node;
      }
      else if (old_node_parent->mRight == old_node)
      {
          old_node_parent->mRight = middle_node;
      }
      middle_node->mParent = old_node_parent;
  }

  //Link the old node with the middle node and make the old node a leaf node now 
  old_node->mParent = middle_node;
  old_node->mLeft = nullptr;
  old_node->mRight = nullptr;
  old_node->mHeight = 0;

  new_node->mParent = middle_node;
  //Update the parent of the tree
  UpdateTreeParents(new_node); 

  //Balance the Tree 
  BalanceTree();
  return;
}


/*
    If a node fully contains the new aabb then do not alter the tree structure. Otherwise
    remove then re-insert the updated node
*/
void DynamicAabbTree::UpdateData(SpatialPartitionKey& key, SpatialPartitionData& data)
{
  //Warn("Assignment3: Required function un-implemented");
  //Search for it first, so that subscript operator will not create a new pair with empty node
    /******Student:Assignment3******/
  if (mKeyNode_map.find(key.mUIntKey) == mKeyNode_map.end())
      return;
  Aabb new_aabb = data.mAabb;
  Aabb old_aabb = mKeyNode_map[key.mUIntKey]->mAabb;
  ////Update the AABB if the old AABB does NOT fully contain the new AABB
  if (old_aabb.Contains(new_aabb) == false)
  {
      //update by removing the old and inserting the new
      RemoveData(key);
      InsertData(key, data);
  }
  return;
}
/*
    After removing the node from the tree you must apply a re-balance to the tree at the
    sibling of the node that was removed as well as all parents going up the tree. Don’t forget that
    more than one balance can happen!
*/
void DynamicAabbTree::RemoveData(SpatialPartitionKey& key)
{
  //Warn("Assignment3: Required function un-implemented");
    /******Student:Assignment3******/
  //Search for it first, so that subscript operator will not create a new pair with empty node
  if (mKeyNode_map.find(key.mUIntKey) == mKeyNode_map.end())
      return;

  Node* ToDelete_node = mKeyNode_map[key.mUIntKey];
  if (ToDelete_node == nullptr)
      return; //prevent seg fault

  //Check if the ToDelete_node is the root node, and the only node in the tree
  if (ToDelete_node->mParent == nullptr)
  {
      //set the root node to null. The Tree will now be empty 
      mRoot = nullptr;
      //delete the memory allocated for it, and it's data in the map
      delete ToDelete_node;
      mKeyNode_map.erase(key.mUIntKey);
      return;
  }

  visted_nodes_stack.empty(); // empty out the stack 

  //Identify the sibling node 
  Node* sibling_node = nullptr;
  if (ToDelete_node->mParent->mLeft == ToDelete_node)
      sibling_node = ToDelete_node->mParent->mRight;
  else
      sibling_node = ToDelete_node->mParent->mLeft;
  //Identify the grandparent node 
  Node* grandparent_node = ToDelete_node->mParent->mParent;
  //Check if the tree only has ToDelete and sibling
  if (grandparent_node == nullptr)
  {
      mRoot = sibling_node;
      //delete the middle node data, i.e the old root node 
      if (mRoot->mParent)
      {
          delete mRoot->mParent;
          mRoot->mParent = nullptr;
      }
  }
  else // make the grand parent the parent of the sibling node 
  {
      if (grandparent_node->mLeft == ToDelete_node->mParent)
      {
          grandparent_node->mLeft = sibling_node;
      }
      else
          grandparent_node->mRight = sibling_node;
      sibling_node->mParent = grandparent_node;
  }

  //Update the parents aabb and height
  UpdateTreeParents(sibling_node);
  //Need to re-balance on insert and remove 
  Node* traversed_node = sibling_node;
  while (traversed_node)
  {
      visted_nodes_stack.push(traversed_node);
      traversed_node = traversed_node->mParent;
  }
  BalanceTree();
  return;
}

void DynamicAabbTree::DebugDraw(int level, const Math::Matrix4& transform, const Vector4& color, int bitMask)
{
 // Warn("Assignment3: Required function un-implemented");
    /******Student:Assignment3******/
    int depth = 0;
    PreOrderDraw(mRoot, level, transform, color, bitMask, depth);
}

void DynamicAabbTree::CastRay(const Ray& ray, CastResults& results)
{
 // Warn("Assignment3: Required function un-implemented");
    /******Student:Assignment3******/
    CastRay_Rec(mRoot, ray, results);
}

void DynamicAabbTree::CastFrustum(const Frustum& frustum, CastResults& results)
{
  //Warn("Assignment3: Required function un-implemented");
    /******Student:Assignment3******/
    CastFrustum_Rec(mRoot, frustum, results);
}

void DynamicAabbTree::SelfQuery(QueryResults& results)
{
  //Warn("Assignment3: Required function un-implemented");
    /******Student:Assignment3******/
    SelfQuery(mRoot, results);
}

void DynamicAabbTree::FilloutData(std::vector<SpatialPartitionQueryData>& results) const
{
  //Warn("Assignment3: Required function un-implemented");
    /******Student:Assignment3******/
    int depth = 0;
   PreOrderTraverse(mRoot, results, depth);
}

DynamicAabbTree::Node* DynamicAabbTree::SelectNode(Node* insertingNode, Node* n0, Node* n1)
{
    //check if n0 or n1 is nullptr 
    /******Student:Assignment3******/
    if (n0 == nullptr  || n1 == nullptr)
    {
        return nullptr;
    }

    //Get the Aabb of left and right nodes
    Aabb n0_aabb = n0->mAabb;
    Aabb n1_aabb = n1->mAabb;
    Aabb insert_aabb = insertingNode->mAabb;
    Aabb c0 = Aabb::Combine(insert_aabb, n0_aabb);
    Aabb c1 = Aabb::Combine(insert_aabb, n1_aabb);

    float a0 = c0.GetSurfaceArea() - n0_aabb.GetSurfaceArea();
    float a1 = c1.GetSurfaceArea() - n1_aabb.GetSurfaceArea();
    if (a0 < a1)
        return n0;
    return n1;

}

void DynamicAabbTree::UpdateTreeParents(Node* leaf_node)
{
    /******Student:Assignment3******/
    //Get the parent of the lead node 
    Node* parent_node = leaf_node->mParent;
    //Recursively move upwards while updating the parents 
    while (parent_node)
    {
        //Update the height of the parent node 
        if (parent_node->mLeft->mHeight > parent_node->mRight->mHeight)
        {
            //set the parent node's height to that + 1
            parent_node->mHeight = parent_node->mLeft->mHeight + 1;
        }
        else
        {
            parent_node->mHeight = parent_node->mRight->mHeight + 1;
        }
        //Update the parent's aabb 
        parent_node->mAabb = parent_node->mAabb.Combine(parent_node->mLeft->mAabb, parent_node->mRight->mAabb);
        //move up to the next parent 
        parent_node = parent_node->mParent;
    }

    return;
}

void DynamicAabbTree::BalanceTree()
{
    /******Student:Assignment3******/
    //iterate through the stack of visited nodes 
    while (visted_nodes_stack.empty() == false)
    {
        Node* old_parent_node = visted_nodes_stack.top();
        visted_nodes_stack.pop(); //remove the node from the stack, to move to the next node
        //Check if the node is a null ptr
        if (old_parent_node == nullptr)
            continue;
        //Check if the node is a leaf node, if it is skip it, don't need to balance it 
        if (old_parent_node->mHeight == 0)
            continue;

        //Check if the current node is unbalanced 
        if (Math::Abs(old_parent_node->mLeft->mHeight - old_parent_node->mRight->mHeight) > 1)
        {
            //identify the pivot and small child node. don't need the large child as we do not touch it
            Node* pivot_node = nullptr;
            Node* small_child = nullptr;
            bool smallChild_wasLeft = false; //used to determine where to attach the small child 
            bool pivot_wasLeft = false;
            //identify the pivot node base using the old parent node
            if (old_parent_node->mLeft->mHeight > old_parent_node->mRight->mHeight)
            {
                pivot_node = old_parent_node->mLeft;
                pivot_wasLeft = true;
            }
            else
            {
                pivot_node = old_parent_node->mRight;
                pivot_wasLeft = false;
            }

            //identify the small child node using the pivot node, note whether it was left side or not 
            if (pivot_node->mLeft->mHeight >= pivot_node->mRight->mHeight)
            {
                small_child = pivot_node->mRight;
                smallChild_wasLeft = false;
            }
            else
            {
                small_child = pivot_node->mLeft;
                smallChild_wasLeft = true;
            }

            //Replace the grand-parent connection of old parent and pivot
            pivot_node->mParent = old_parent_node->mParent; //swap their parents first 

            if (old_parent_node->mParent == nullptr) // check if old parent is the root node 
            {
                //if(pivot_node->mParent != nullptr)
                mRoot = pivot_node;
                mRoot->mParent = nullptr;
            }
            else
            {
                //Check whether the parent node is the left/right child of it's own parent
                if (old_parent_node->mParent->mLeft == old_parent_node)
                {
                    //attach the pivot node to where the old parent was
                    old_parent_node->mParent->mLeft = pivot_node;
                }
                else
                {
                    old_parent_node->mParent->mRight = pivot_node;
                }
            }
            //Attach where the old parent was under the pivot where small child was 
            old_parent_node->mParent = pivot_node;

            if (smallChild_wasLeft)
                pivot_node->mLeft = old_parent_node;
            else
                pivot_node->mRight = old_parent_node;

            //Attach the small child's parent to be the old parent node
            small_child->mParent = old_parent_node;

            //Attach small child to old parent where pivot was 
            if (pivot_wasLeft)
                old_parent_node->mLeft = small_child;
            else
                old_parent_node->mRight = small_child;

            //update the aabbs and heights 
            UpdateTreeParents(small_child);

             //Check from small child again if there's a need for further balancing
            while (small_child)
            {
                visted_nodes_stack.push(small_child);
                small_child = small_child->mParent;
            }
            
        }

        
    }
    return;
}

void DynamicAabbTree::PreOrderTraverse(Node* curr_node, std::vector<SpatialPartitionQueryData>& results, int depth) const
{
    /******Student:Assignment3******/
    if (curr_node == nullptr)
        return;
    SpatialPartitionQueryData mData{};
    mData.mAabb = curr_node->mAabb;
    mData.mClientData = curr_node->mClientData;
    mData.mDepth = depth;
    results.push_back(mData);
    PreOrderTraverse(curr_node->mLeft, results, depth + 1);
    PreOrderTraverse(curr_node->mRight, results, depth + 1);
}

void DynamicAabbTree::PreOrderDraw(Node* curr_node, int level, const Math::Matrix4& transform, const Vector4& color, int bitMask, int depth)
{
    /******Student:Assignment3******/
    if (curr_node == nullptr)
        return;
    if (level == -1 || depth == level )
    {
        DebugShape ToDraw = curr_node->mAabb.DebugDraw();
        ToDraw.SetMaskBit(bitMask);
        ToDraw.Color(color);
        ToDraw.SetTransform(transform);
    }
    PreOrderDraw(curr_node->mLeft, level, transform, color, bitMask, depth + 1);
    PreOrderDraw(curr_node->mRight, level, transform, color, bitMask, depth + 1);
}


void DynamicAabbTree::CastRay_Rec(Node* curr_node, const Ray& ray, CastResults& results)
{
    /******Student:Assignment3******/
    if (curr_node == nullptr)
        return;
    CastResult m_result;
    m_result.mClientData = curr_node->mClientData;
    if (RayAabb(ray.mStart, ray.mDirection, curr_node->mAabb.GetMin(), curr_node->mAabb.GetMax(), m_result.mTime) )
    {
        //if it's a leaf node 
        if (curr_node->mLeft == nullptr || curr_node->mRight == nullptr)
        {
            //return the result if there's a hit 
            results.AddResult(m_result);
        }
        else
        {
            //continue recursing down 
            CastRay_Rec(curr_node->mLeft, ray, results);
            CastRay_Rec(curr_node->mRight, ray, results);
        }
    }
}

void DynamicAabbTree::CastFrustum_Rec(Node* curr_node, const Frustum& frustum, CastResults& results)
{
    /******Student:Assignment3******/
    if (curr_node == nullptr)
        return;
    IntersectionType::Type intersection_result = IntersectionType::NotImplemented;
    //mLastAxis is cache in the node
    intersection_result = FrustumAabb(frustum.GetPlanes(), curr_node->mAabb.GetMin(), curr_node->mAabb.GetMax(), curr_node->mLastAxis);

    //if the aabb is in the frustum, return all children 
    if (intersection_result == IntersectionType::Outside)
        return; // stop recursion 
    if (intersection_result == IntersectionType::Inside )
    {
        //add the ALL the leaves node of this subtree to the result, no need for further frustumaabb tests 
        CastFrustum_ReturnLeavesNode_Rec(curr_node, results);
    }
    if (intersection_result == IntersectionType::Overlaps)
    {
        //return the result if there's a hit 
        if (curr_node->mLeft == nullptr || curr_node->mRight == nullptr)
        {
            CastResult m_result;
            m_result.mClientData = curr_node->mClientData;
            //return the result if there's a hit 
            results.AddResult(m_result);
        }
        else
        {
            //continue recursing down
            CastFrustum_Rec(curr_node->mLeft, frustum, results);
            CastFrustum_Rec(curr_node->mRight, frustum, results);
        }
    }
        
    
    
}

void DynamicAabbTree::CastFrustum_ReturnLeavesNode_Rec(Node* curr_node, CastResults& results)
{
    /******Student:Assignment3******/
    if (curr_node == nullptr)
        return;
    if (curr_node->mLeft == nullptr || curr_node->mRight == nullptr)
    {
        CastResult m_result;
        m_result.mClientData = curr_node->mClientData;
        //return the result if there's a hit 
        results.AddResult(m_result);
    }
    else
    {
        //continue recursing down
        CastFrustum_ReturnLeavesNode_Rec(curr_node->mLeft, results);
        CastFrustum_ReturnLeavesNode_Rec(curr_node->mRight, results);
    }
}

void DynamicAabbTree::SelfQuery(Node* node, QueryResults& results)
{
    /******Student:Assignment3******/
    if (node->mHeight == 0)
    {
        return;
    }

    SelfQuery(node->mLeft, node->mRight, results);
    SelfQuery(node->mLeft, results);
    SelfQuery(node->mRight, results);
}

void DynamicAabbTree::SelfQuery(Node* nodeA, Node* nodeB, QueryResults& results)
{
    /******Student:Assignment3******/
    // Not overlapping
    if (!AabbAabb(nodeA->mAabb.GetMin(), nodeA->mAabb.GetMax(), nodeB->mAabb.GetMin(), nodeB->mAabb.GetMax()))
    {
        return;
    }

    // Case 1 :Both leaf
    if (nodeA->mHeight == 0 && nodeB->mHeight == 0)
    {
        void* clientdata_A = nodeA->mClientData;
        void* clientdata_B = nodeB->mClientData;

        results.AddResult(QueryResult(clientdata_A, clientdata_B));
    }
    else if (nodeA->mHeight == 0) //Case 2 :one internal, one leaf
    {
        SelfQuery(nodeA, nodeB->mLeft, results);
        SelfQuery(nodeA, nodeB->mRight, results);
    }
    else if (nodeB->mHeight == 0) //Case 2 :one internal, one leaf
    {
        SelfQuery(nodeA->mLeft, nodeB, results);
        SelfQuery(nodeA->mRight, nodeB, results);
    }
    else // Case 3: Both internal split nodes
        SplitNodes(nodeA, nodeB, results);
}

void DynamicAabbTree::SplitNodes(Node* nodeA, Node* nodeB, QueryResults& results)
{
    /******Student:Assignment3******/
        // Both internal nodes
    if (nodeA->mAabb.GetVolume() > nodeB->mAabb.GetVolume()) // Use larger volume as heuristic, as stated
    {
        SelfQuery(nodeA->mLeft, nodeB, results);
        SelfQuery(nodeA->mRight, nodeB, results);
    }
    else
    {
        SelfQuery(nodeA, nodeB->mLeft, results);
        SelfQuery(nodeA, nodeB->mRight, results);
    }
}
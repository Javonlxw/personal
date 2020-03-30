///////////////////////////////////////////////////////////////////////////////
///
/// Authors: Joshua Davis
/// Copyright 2015, DigiPen Institute of Technology
///
///////////////////////////////////////////////////////////////////////////////
#include "Precompiled.hpp"

//-----------------------------------------------------------------------------NSquaredSpatialPartition
NSquaredSpatialPartition::NSquaredSpatialPartition()
{
  mType = SpatialPartitionTypes::NSquared;
}

void NSquaredSpatialPartition::InsertData(SpatialPartitionKey& key, SpatialPartitionData& data)
{
  // Doing this lazily (and bad, but it's n-squared...).
  // Just store as the key what the client data is so we can look it up later.
  key.mVoidKey = data.mClientData;
  mData.push_back(data.mClientData);
}

void NSquaredSpatialPartition::UpdateData(SpatialPartitionKey& key, SpatialPartitionData& data)
{
  // Nothing to do here, update doesn't do anything
}

void NSquaredSpatialPartition::RemoveData(SpatialPartitionKey& key)
{
  // Find the key data and remove it
  for(size_t i = 0; i < mData.size(); ++i)
  {
    if(mData[i] == key.mVoidKey)
    {
      mData[i] = mData.back();
      mData.pop_back();
      break;
    }
  }
}

void NSquaredSpatialPartition::DebugDraw(int level, const Math::Matrix4& transform, const Vector4& color, int bitMask)
{
  // Nothing to debug draw
}

void NSquaredSpatialPartition::CastRay(const Ray& ray, CastResults& results)
{
  // Add everything
  for(size_t i = 0; i < mData.size(); ++i)
  {
    CastResult result;
    result.mClientData = mData[i];
    results.AddResult(result);
  }
}

void NSquaredSpatialPartition::CastFrustum(const Frustum& frustum, CastResults& results)
{
  // Add everything
  for(size_t i = 0; i < mData.size(); ++i)
  {
    CastResult result;
    result.mClientData = mData[i];
    results.AddResult(result);
  }
}

void NSquaredSpatialPartition::SelfQuery(QueryResults& results)
{
  // Add everything
  for(size_t i = 0; i < mData.size(); ++i)
  {
    for(size_t j = i + 1; j < mData.size(); ++j)
    {
      results.AddResult(QueryResult(mData[i], mData[j]));
    }
  }
}

void NSquaredSpatialPartition::GetDataFromKey(const SpatialPartitionKey& key, SpatialPartitionData& data) const
{
  data.mClientData = key.mVoidKey;
}

void NSquaredSpatialPartition::FilloutData(std::vector<SpatialPartitionQueryData>& results) const
{
  for(size_t i = 0; i < mData.size(); ++i)
  {
    SpatialPartitionQueryData data;
    data.mClientData = mData[i];
    results.push_back(data);
  }
}

//-----------------------------------------------------------------------------BoundingSphereSpatialPartition
BoundingSphereSpatialPartition::BoundingSphereSpatialPartition()
{
  mType = SpatialPartitionTypes::NSquaredSphere;
}

void BoundingSphereSpatialPartition::InsertData(SpatialPartitionKey& key, SpatialPartitionData& data)
{
  //Warn("Assignment2: Required function un-implemented");
    // Add the data into our BSP map, linked to the key 
    //mBSSPmap.insert(std::pair<void*, Sphere>(key.mVoidKey, data.mBoundingSphere));
    key.mUIntKey = key.mUIntKey++;
    key.mVoidKey = data.mClientData;
    mBSSPmap.emplace(key.mVoidKey, data.mBoundingSphere);
}

void BoundingSphereSpatialPartition::UpdateData(SpatialPartitionKey& key, SpatialPartitionData& data)
{
    //Warn("Assignment2: Required function un-implemented");
    mBSSPmap[key.mVoidKey] = data.mBoundingSphere;
}

void BoundingSphereSpatialPartition::RemoveData(SpatialPartitionKey& key)
{
  //Warn("Assignment2: Required function un-implemented");
    //make use of std::map erase to erase the data from our BSP map
    mBSSPmap.erase(key.mVoidKey);
}

void BoundingSphereSpatialPartition::DebugDraw(int level, const Math::Matrix4& transform, const Vector4& color, int bitMask)
{
  //Warn("Assignment2: Required function un-implemented");
}

void BoundingSphereSpatialPartition::CastRay(const Ray& ray, CastResults& results)
{
  //Warn("Assignment2: Required function un-implemented");
    //Iterate through the BSSP map and check against all the spheres with the ray 
    for (auto& e : mBSSPmap)
    {
        float t = 0.f;
        Sphere s = e.second;
        //Make use of our Ray-Sphere intersection to test if it collides
        if (RaySphere(ray.mStart, ray.mDirection, s.mCenter, s.mRadius, t))
        {
            //if it does, Cast the result and add it to the list of results 
            CastResult m_result(e.first, t);
            results.AddResult(m_result);
        }
    }
}

void BoundingSphereSpatialPartition::CastFrustum(const Frustum& frustum, CastResults& results)
{
 // Warn("Assignment2: Required function un-implemented");
    //Iterate through the BSP map and check against all the spheres with the ray 
    for (auto& e : mBSSPmap)
    {
        //Sphere s = e.second;
        size_t axis = 0;
        // use our FrustumSphere instersection test to see if there's collision
        IntersectionType::Type res = FrustumSphere(frustum.GetPlanes(), e.second.mCenter, e.second.mRadius, axis);
        if ( res== IntersectionType::Inside || res == IntersectionType::Overlaps)
        {
            //if it does, Cast the result and add it to the list of results
            CastResult m_result(e.first);
            results.AddResult(m_result);
        }
    }

}

void BoundingSphereSpatialPartition::SelfQuery(QueryResults& results)
{
  //Warn("Assignment2: Required function un-implemented");
    typedef std::map<void*, Sphere>::iterator BSSP_iterator;

    for (BSSP_iterator s1_it = mBSSPmap.begin(); s1_it != mBSSPmap.end(); ++s1_it)
    {
        BSSP_iterator s2_it = s1_it;
        for (++s2_it ; s2_it != mBSSPmap.end(); ++s2_it)
        {
            Sphere s1 = s1_it->second, s2 = s2_it->second;
            if (SphereSphere(s1.GetCenter(), s1.GetRadius(), s2.GetCenter(), s2.GetRadius()))
            {
                //Add the query result with the void key of both spheres
                results.AddResult(QueryResult(s1_it->first, s2_it->first));
            }
        }
    }
}

void BoundingSphereSpatialPartition::FilloutData(std::vector<SpatialPartitionQueryData>& results) const
{
  //Warn("Assignment2: Required function un-implemented");
    //Iterate through the map and push all the void key data into the results 
    for (auto& e : mBSSPmap)
    {
        SpatialPartitionQueryData mData; 
        mData.mClientData = e.first;
        mData.mBoundingSphere = e.second;
        results.push_back(mData);
    }

}

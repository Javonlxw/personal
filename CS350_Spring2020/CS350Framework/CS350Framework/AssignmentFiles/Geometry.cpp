///////////////////////////////////////////////////////////////////////////////
///
/// Authors: Joshua Davis
/// Copyright 2015, DigiPen Institute of Technology
///
///////////////////////////////////////////////////////////////////////////////
#include "Precompiled.hpp"


Vector3 ProjectPointOnPlane(const Vector3& point, const Vector3& normal, float planeDistance)
{
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
    //ensures that the normal is normalized 
    Vector3 n = normal.Normalized();
    Vector3 p0 = planeDistance * n; //since n.p0 = planeDistance, p0 = <dx,dy,dz> ,where n = <x,y,z>
    float w = (point - p0).Dot(n);  //let w be the perpendicular distance of the point to the plane, w = (p-p0).n
    return point - (n*w); // p' = p - w*n 
}

bool BarycentricCoordinates(const Vector3& point, const Vector3& a, const Vector3& b,
                            float& u, float& v, float epsilon)
{
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
  // FOR A LINE
    //Check if the line is degenerated: 
    if ((a - b).LengthSq() < epsilon)
    {
        u = 0.f;
        v = 0.f;
        return false;
    }
    //Using the analytic method:
    u = ((point - b).Dot(a - b)) / (a - b).LengthSq(); //note: (a-b).(a-b) = (a - b).LengthSq
    v = 1.0f - u;
    //Check if u is within the range of [-epsilon,1+epsilon]. If it is, so will be v, vice versa 
    if (((u <= -epsilon) || (u >= (1.0f + epsilon))) )
        return false; // the point is outside the line
    else
        return true;  // the point is inside the line 
}

bool BarycentricCoordinates(const Vector3& point, const Vector3& a, const Vector3& b, const Vector3& c,
                            float& u, float& v, float& w, float epsilon)
{
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
  //FOR A TRIANGLE
    //sub w = 1 - u -v into P = uA + vB +wC, we get v0 = u*v1 + v*v2, where v0,v1,v2 is: 
    Vector3 v0 = point - c;
    Vector3 v1 = a - c;
    Vector3 v2 = b - c;
    // Solve using Cramer's rule: a1*u + b1*v = e1 and c1*u + d1*v = f1, where a,b,c,d,e,f is: 
    float a1 = v1.Dot(v1);
    float b1 = v2.Dot(v1);
    float c1 = v1.Dot(v2);
    float d1 = v2.Dot(v2);
    float e1 = v0.Dot(v1);
    float f1 = v0.Dot(v2);
    // Using Cramer's rule, we can deduce that u and v is: 
    if ((a1 * d1 - b1 * c1))
    {
        u = (e1 * d1 - b1 * f1) / (a1 * d1 - b1 * c1);
        v = (a1 * f1 - e1 * c1) / (a1 * d1 - b1 * c1);
        w = 1.0f - u - v;
    }
    else
    {
        // the signed areas will be 0 as the 3 points are on a straight line 
        u = 0;
        v = 0;
        w = 0;
        return false; //i.e the triangle is degenerated 
    }
    
    // Check if the passed in point is outside the triangle. if it is then:
    // either u,v or w will be not within the range of [-epsilon, 1+epsilon], i.e u,v,w + epsilon < 0
    // it is possible for u and v to be in range, but w to be out of range, e.g u = 0.4, v = 0.8, w = -0.2
    if (v + epsilon < 0 ||
        w + epsilon < 0 ||
        u + epsilon < 0)
        return false; // the point is outside the triangle
    else
        return true;  // the point is inside the triangle 
}

IntersectionType::Type PointPlane(const Vector3& point, const Vector4& plane, float epsilon)
{
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
  Vector4 p{ point.x, point.y, point.z, -1 };
  float w = plane.Dot(p); // since plane = <nx,ny,nz,d>, and <nx,ny,nz,d>.<px,py,pz,-1> = w, from theory 
  if (w > epsilon)
      return IntersectionType::Inside;
  else if (w < -epsilon)
      return IntersectionType::Outside;
  else //if neither greater or lesser, it's in-between epsilon values and thus co-planar
      return IntersectionType::Coplanar;
}

bool PointSphere(const Vector3& point, const Vector3& sphereCenter, float sphereRadius)
{
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
    if ((point - sphereCenter).LengthSq() <= sphereRadius * sphereRadius)
        return true; //Surface of sphere is considered to be part of the sphere
    else
        return false;
}

bool PointAabb(const Vector3& point, const Vector3& aabbMin, const Vector3& aabbMax)
{
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
    if (((point.x >= aabbMin.x) && (point.x <= aabbMax.x)) &&
        ((point.y >= aabbMin.y) && (point.y <= aabbMax.y)) &&
        ((point.z >= aabbMin.z) && (point.z <= aabbMax.z)))
        return true;
    else
        return false;
}

bool RayPlane(const Vector3& rayStart, const Vector3& rayDir,
              const Vector4& plane, float& t, float epsilon)
{
  ++Application::mStatistics.mRayPlaneTests;
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
  //compute the norrmal of the plane 
  Vector3 n{ plane.x, plane.y, plane.z };
  //normalized the normal to be safe
  n.Normalize();
  //compute p0
  Vector3 p0{ plane.z * plane.x,plane.z * plane.y,plane.z * plane.z }; // as p0 = <da,db,dc>, where n = <a,b,c>
  //check first if n.rayDir is 0, to prevent zero division, if it is 0, it means the ray is parallel to the plane and will never intersect
  // also we use Epsilon to check if the ray is close enough to parallel to the plane, if it's between epsilon value, return false.
  if ((n.Dot(rayDir) <= epsilon && n.Dot(rayDir) >= -epsilon) || n.Dot(rayDir) == 0)
      return false; //do not fill t if ray does not hit.
  // compute the value of 't', which is the time of intersection if the ray intersects with the plane
  t = (plane.w - n.Dot(rayStart)) / (n.Dot(rayDir));
  if (t > 0)
      return true;
  else
      return false; //the ray is behind the plane 
}

bool RayTriangle(const Vector3& rayStart, const Vector3& rayDir,
                 const Vector3& triP0, const Vector3& triP1, const Vector3& triP2,
                 float& t, float triExpansionEpsilon)
{
  ++Application::mStatistics.mRayTriangleTests;
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
  //define a plane eqn that the triangle lies on 
  Plane m_plane;
  m_plane.Set(triP0, triP1, triP2);
  // call RayPlane to obtain 't', which will be used to get point_intersection 
  if (RayPlane(rayStart, rayDir, m_plane.mData, t, triExpansionEpsilon))
  {
      //find p_intersect where p_intersect = rayStart + (t*rayDir)
      Vector3 p_intersect = rayStart + (t * rayDir);
      //use BaryCentric coordinates to check if the p_intersects falls inside/on the triangle 
      float u, v, w = 0.f;
      return BarycentricCoordinates(p_intersect, triP0, triP1, triP2, u, v, w, triExpansionEpsilon);
  }
  else
      return false;
}

bool RaySphere(const Vector3& rayStart, const Vector3& rayDir,
               const Vector3& sphereCenter, float sphereRadius,
               float& t)
{
  ++Application::mStatistics.mRaySphereTests;
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
  // t will be 0 if rayStart is inside the sphere, use PointSphere to check if it's true 
  if (PointSphere(rayStart, sphereCenter, sphereRadius))
  {
      t = 0; //set as required
      return true;
  }
  // by theory, t^2(dir^2) + t(-2m.dir) + (m^2 - r^2), where m = sphereCenter - rayStart
  //to find t, we can use the quadratic formula, where a,b,c is:
  float a = rayDir.LengthSq(); // again, lengthSq is also a dot product of itself, i.e rayDir.rayDir
  if (a == 0)
      return false; //the ray is not a ray, but a point. since rayStart is not in the sphere, there'll be no intersection for sure 
  Vector3 m = sphereCenter - rayStart;
  float b = -2.f * (m.Dot(rayDir)); // -2m.dir
  float c = m.LengthSq() - sphereRadius * sphereRadius; // c = m^2 - r^2
  //hardcoded way to find t 1 and t 2 values, via quadratic formula
  float t1 = (-b + Math::Sqrt(b * b - 4.f * a * c)) / (2 * a);
  float t2 = (-b - Math::Sqrt(b * b - 4.f * a * c)) / (2 * a);
  if (t1 <= t2 && t1 >= 0)
  {
      t = t1; //set t to be t1 since t1 is smaller 
      return true;
  }
  else if (t2 >= 0)
  {
      t = t2; //set t to be t2 since t2 is smaller than t1
      return true;
  }
  else
      return false;// at this point, it means both t1 and t2 is negative, i.e ray is behind the sphere 
}

bool RayAabb(const Vector3& rayStart, const Vector3& rayDir,
             const Vector3& aabbMin, const Vector3& aabbMax, float& t)
{
  ++Application::mStatistics.mRayAabbTests;
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
  //return false;
    //Check if the starting pt for Ray is inside aabb
    if (PointAabb(rayStart, aabbMin, aabbMax))
    {
        t = 0;
        return true;
    }
    // By subbing n = <1,0,0> or <0,1,0> or <0,0,1> into the formula n.(sr+t.dr-p0) = 0, we get :
    // t = (p0.i - rayStart.i)/rayDir.i, where i = x,y or z and p0 = aabbMin or aabbMax 
    // with this there'll be a total of max 6 solutions for t, and we would want to find the smallest t, that's greater than 0 
    float time_min[3]; //let this be the minimum time for x,y,z
    float time_max[3]; //let this be the max time for x,y,z

    float t1 = 0; //smallest possible time at the start
    float t2 = std::numeric_limits<float>::max(); //largest possible time at the start 

    //use a for loop to iterate through x,y,z-axis
    for (int i = 0; i < 3; ++i)
    {
        //check in each axis, is the rayDir 0
        if (rayDir[i] == 0) //note that the Vector3 is stored as a float array
        {
            //check if the rayStart is OUTSIDE the aabb in that axis
            if (rayStart[i] < aabbMin[i] || rayStart[i] > aabbMax[i])
                return false;
            else
                continue;//we need to ignore it to prevent zero-division
        }
        time_min[i] = (aabbMin[i] - rayStart[i]) / rayDir[i];
        time_max[i] = (aabbMax[i] - rayStart[i]) / rayDir[i];
        //Check for the ray's direction, if negative, swap min and max values
        if (rayDir[i] < 0)
        {
            float tmp = time_max[i];
            time_max[i] = time_min[i];
            time_min[i] = tmp;
        }
        //finally, get the correct t1 and t2 by comparing the min and max for each axis 
        t1 = Math::Max(t1, time_min[i]);
        t2 = Math::Min(t2, time_max[i]);
    }

    //Now we need to check t1 and t2 
    if (t1 > t2)
        return false; //no intersection
    else
    {
        //t1 should always be smaller than t2, and will never be negative as it 
        //starts with 0, and takes the max(0,t_min)
        t = t1;
        return true;
    }
}

IntersectionType::Type PlaneTriangle(const Vector4& plane, 
                                     const Vector3& triP0, const Vector3& triP1, const Vector3& triP2,
                                     float epsilon)
{
  ++Application::mStatistics.mPlaneTriangleTests;
  /******Student:Assignment1******/
  //Warn("Assignment1: Required function un-implemented");
  //return IntersectionType::NotImplemented;
  // Construct 3 rays base on the triangle points/edges
  Ray ray_p0_p1{ triP0, (triP1 - triP0).Normalized() };
  Ray ray_p0_p2{ triP0, (triP2 - triP0).Normalized() };
  Ray ray_p1_p2{ triP1, (triP2 - triP1).Normalized() };
  //initialize the 't' values of each ray to be -1
  float t01 = -1;
  float t02 = -1;
  float t12 = -1;
  // Run RayPlane Intersection test 
  RayPlane(ray_p0_p1.mStart, ray_p0_p1.mDirection, plane, t01, epsilon);
  RayPlane(ray_p0_p2.mStart, ray_p0_p2.mDirection, plane, t02, epsilon);
  RayPlane(ray_p1_p2.mStart, ray_p1_p2.mDirection, plane, t12, epsilon);
  // Run PointPlane Intersection test for all 3 Triangle points
  IntersectionType::Type p0_type = PointPlane(triP0, plane, epsilon);
  IntersectionType::Type p1_type = PointPlane(triP1, plane, epsilon);
  IntersectionType::Type p2_type = PointPlane(triP2, plane, epsilon);
  //Triangle can only be "coplanar" if all 3 points are coplanar 
  if (p0_type == IntersectionType::Coplanar &&
      p1_type == IntersectionType::Coplanar &&
      p2_type == IntersectionType::Coplanar)
      return p0_type;
  //Check for valid t-values: Test if triangle Overlaps
  if (t01 > 0 && t01 * t01 <= (triP1 - triP0).LengthSq()) //t-value obtained from ray_p0_p1
  {
      //Check if p0 and p1 are coplanar 
      if (p0_type != IntersectionType::Coplanar && p1_type != IntersectionType::Coplanar)
      {
          //that means there's an overlap
          return IntersectionType::Overlaps;
      }
  }
  if (t02 > 0 && t02 * t02 <= (triP2 - triP0).LengthSq()) // t-value obtained from ray_p0_p2
  {
      //Check if p0 and p1 are coplanar 
      if (p0_type != IntersectionType::Coplanar && p2_type != IntersectionType::Coplanar)
      {
          //that means there's an overlap
          return IntersectionType::Overlaps;
      }
  }
  if (t12 > 0 && t12 * t12 <= (triP2 - triP1).LengthSq()) // t-value obtained from ray_p1_p2
  {
      //Check if p0 and p1 are coplanar 
      if (p2_type != IntersectionType::Coplanar && p1_type != IntersectionType::Coplanar)
      {
          //that means there's an overlap
          return IntersectionType::Overlaps;
      }
  }
  //if t-values are invalid, it means the triangle is either Inside or Outside, at this point.
  //i.e it'll be the same IntersectionType as one of the non-coplanar points 
  if (p0_type != IntersectionType::Coplanar)
      return p0_type;
  else if (p1_type != IntersectionType::Coplanar)
      return p1_type;
  else
      return p2_type; 
}

IntersectionType::Type PlaneSphere(const Vector4& plane,
                                   const Vector3& sphereCenter, float sphereRadius)
{
  ++Application::mStatistics.mPlaneSphereTests;
  /******Student:Assignment1******/
    //get the distance,w, of the center of sphere to the plane
    // (Cs-P0).n = Cs.n -P0.n = Cs.n - d.n.n = Cs.n - d, where n is:
    Vector3 n{ plane.x, plane.y, plane.z };
    if (n.LengthSq() == 0)
        return IntersectionType::Outside;
    n.Normalize(); //make sure the normal is normalized
    float w = sphereCenter.Dot(n) - plane.w;
    if (abs(w) <= sphereRadius)
        return IntersectionType::Overlaps;
    else if (w < 0)
        return IntersectionType::Outside;
    else
        return IntersectionType::Inside;
}

IntersectionType::Type PlaneAabb(const Vector4& plane,
                                 const Vector3& aabbMin, const Vector3& aabbMax)
{
  ++Application::mStatistics.mPlaneAabbTests;
  /******Student:Assignment1******/
  //Construct a psuedo sphere that represents aabb
  //Find the center of the aabb
  Vector3 aabbCenter = (aabbMax + aabbMin) / 2;
  //Find the aabb half extent
  Vector3 aabbHExtent = aabbMax - aabbCenter;
  //Calculate the psuedo sphere's radius by taking the HE.dot(normal), where normal is:
  Vector3 normal{ abs(plane.x), abs(plane.y), abs(plane.z )}; //take abs to ensure it's all positive 
  float radius = aabbHExtent.Dot(normal);
  // Call PlaneSphere to handle it
  return PlaneSphere(plane, aabbCenter, radius);
}

IntersectionType::Type FrustumTriangle(const Vector4 planes[6],
                                       const Vector3& triP0, const Vector3& triP1, const Vector3& triP2,
                                       float epsilon)
{
  ++Application::mStatistics.mFrustumTriangleTests;
  /******Student:Assignment1******/
  // Do 6 PlaneTriangle tests, if all are inside, will it be inside
  int counter_inside = 0;
  for (int i = 0; i < 6; ++i)
  {
      IntersectionType::Type m_type = PlaneTriangle(planes[i], triP0, triP1, triP2, epsilon);
      if (m_type == IntersectionType::Outside)
          return m_type; // if one is outside, the whole triangle will be 'outside' 
      else if (m_type == IntersectionType::Inside)
          ++counter_inside;
  }
  //Check if all PlaneTriangle tests are inside, i.e the counter will be 6, if true, it's inside
  if (counter_inside == 6)
      return IntersectionType::Inside;
  else
      return IntersectionType::Overlaps; // as some either overlap/coplanar
  
}

IntersectionType::Type FrustumSphere(const Vector4 planes[6],
                                     const Vector3& sphereCenter, float sphereRadius, size_t& lastAxis)
{
  ++Application::mStatistics.mFrustumSphereTests;
  /******Student:Assignment1******/
  //same method as FrustumTriangle, call PlaneSphere instead 
  int counter_inside = 0;
  for (size_t i = 0; i < 6; ++i)
  {
      IntersectionType::Type m_type = PlaneSphere(planes[i], sphereCenter, sphereRadius);
      if (m_type == IntersectionType::Outside)
      {
          lastAxis = i;
          return m_type; // if one is outside, the whole triangle will be 'outside' 
      }
      else if (m_type == IntersectionType::Inside)
          ++counter_inside;
  }
  //Check if all PlaneTriangle tests are inside, i.e the counter will be 6, if true, it's inside
  if (counter_inside == 6)
      return IntersectionType::Inside;
  else
      return IntersectionType::Overlaps; // as some either overlap/coplanar
}

IntersectionType::Type FrustumAabb(const Vector4 planes[6],
                                   const Vector3& aabbMin, const Vector3& aabbMax, size_t& lastAxis)
{
  ++Application::mStatistics.mFrustumAabbTests;
  /******Student:Assignment1******/
  //same method as FrustumTriangle/FrustumSphere, call PlaneAabb instead
  int counter_inside = 0;
  for (size_t i = lastAxis; i < 6; ++i)
  {
      IntersectionType::Type m_type = PlaneAabb(planes[i], aabbMin, aabbMax);
      if (m_type == IntersectionType::Outside)
      {
          lastAxis = i;
          return m_type; // if one is outside, the whole triangle will be 'outside' 
      }
      else if (m_type == IntersectionType::Inside)
          ++counter_inside;
  }
  //Check if all PlaneTriangle tests are inside, i.e the counter will be 6, if true, it's inside
  if (counter_inside == 6)
      return IntersectionType::Inside;
  else
      return IntersectionType::Overlaps; // as some either overlap/coplanar

}

bool SphereSphere(const Vector3& sphereCenter0, float sphereRadius0,
                  const Vector3& sphereCenter1, float sphereRadius1)
{
  ++Application::mStatistics.mSphereSphereTests;
  /******Student:Assignment1******/
  Vector3 dist = sphereCenter1 - sphereCenter0;
  float totalRadius = sphereRadius0 + sphereRadius1;
  if (dist.LengthSq() <= totalRadius * totalRadius)
      return true;
  else
      return false;
}

bool AabbAabb(const Vector3& aabbMin0, const Vector3& aabbMax0,
              const Vector3& aabbMin1, const Vector3& aabbMax1)
{
  ++Application::mStatistics.mAabbAabbTests;
  /******Student:Assignment1******/
  //check if it's OUTSIDE
  if (aabbMin0.x > aabbMax1.x || aabbMin1.x > aabbMax0.x ||
      aabbMin0.y > aabbMax1.y || aabbMin1.y > aabbMax0.y ||
      aabbMin0.z > aabbMax1.z || aabbMin1.z > aabbMax0.z)
      return false;
  else
      return true;
}

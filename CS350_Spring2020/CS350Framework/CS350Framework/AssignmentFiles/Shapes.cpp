///////////////////////////////////////////////////////////////////////////////
///
/// Authors: Joshua Davis
/// Copyright 2015, DigiPen Institute of Technology
///
///////////////////////////////////////////////////////////////////////////////
#include "Precompiled.hpp"

//-----------------------------------------------------------------------------LineSegment
LineSegment::LineSegment()
{
  mStart = mEnd = Vector3::cZero;
}

LineSegment::LineSegment(Math::Vec3Param start, Math::Vec3Param end)
{
  mStart = start;
  mEnd = end;
}

DebugShape& LineSegment::DebugDraw() const
{
  return gDebugDrawer->DrawLine(*this);
}

//-----------------------------------------------------------------------------Ray
Ray::Ray()
{
  mStart = mDirection = Vector3::cZero;
}

Ray::Ray(Math::Vec3Param start, Math::Vec3Param dir)
{
  mStart = start;
  mDirection = dir;
}

Ray Ray::Transform(const Math::Matrix4& transform) const
{
  Ray transformedRay;
  transformedRay.mStart = Math::TransformPoint(transform, mStart);
  transformedRay.mDirection = Math::TransformNormal(transform, mDirection);
  return transformedRay;
}

Vector3 Ray::GetPoint(float t) const
{
  return mStart + mDirection * t;
}

DebugShape& Ray::DebugDraw(float t) const
{
  return gDebugDrawer->DrawRay(*this, t);
}

//-----------------------------------------------------------------------------PCA Helpers
Matrix3 ComputeCovarianceMatrix(const std::vector<Vector3>& points)
{
  /******Student:Assignment2******/
  //Warn("Assignment2: Required function un-implemented");
  return Matrix3::cIdentity;
}

Matrix3 ComputeJacobiRotation(const Matrix3& matrix)
{
  /******Student:Assignment2******/
  // Compute the jacobi rotation matrix that will turn the largest (magnitude) off-diagonal element of the input
  // matrix into zero. Note: the input matrix should always be (near) symmetric.
  //Warn("Assignment2: Required function un-implemented");
  return Matrix3::cIdentity;
}

void ComputeEigenValuesAndVectors(const Matrix3& covariance, Vector3& eigenValues, Matrix3& eigenVectors, int maxIterations)
{
  /******Student:Assignment2******/
  // Iteratively rotate off the largest off-diagonal elements until the resultant matrix is diagonal or maxIterations.
  //Warn("Assignment2: Required function un-implemented");
}


//-----------------------------------------------------------------------------Sphere
Sphere::Sphere()
{
  mCenter = Vector3::cZero;
  mRadius = 0;
}

Sphere::Sphere(const Vector3& center, float radius)
{
  mCenter = center;
  mRadius = radius;
}

void Sphere::ComputeCentroid(const std::vector<Vector3>& points)
{
  /******Student:Assignment2******/
  // The centroid method is roughly describe as: find the centroid (not mean) of all
  // points and then find the furthest away point from the centroid.
  //Warn("Assignment2: Required function un-implemented");
    Vector3 min = points.front();
    Vector3 max = points.front();
    //Find the min and max points from the cluster of points 
    for (const Vector3& p : points)
    {
        if (min.x > p.x)
            min.x = p.x;
        if (min.y > p.y)
            min.y = p.y;
        if (min.z > p.z)
            min.z = p.z;

        if (max.x < p.x)
            max.x = p.x;
        if (max.y < p.y)
            max.y = p.y;
        if (max.z < p.z)
            max.z = p.z;
    }
    //use the min and max points to compute the center point for the sphere
    mCenter = 0.5f * (max + min);
    float rad_sq = (points.front() - mCenter).LengthSq();
    for (const Vector3& p : points)
    {
        if (rad_sq < (p - mCenter).LengthSq() )
        {
            rad_sq = (p - mCenter).LengthSq();
        }
    }
    mRadius = Math::Sqrt(rad_sq);
    return;
}

void Sphere::ComputeRitter(const std::vector<Vector3>& points)
{
  /******Student:Assignment2******/
  // The ritter method:
  // Find the largest spread on each axis.
  // Find which axis' pair of points are the furthest (euclidean distance) apart.
  // Choose the center of this line as the sphere center. Now incrementally expand the sphere.
 // Warn("Assignment2: Required function un-implemented");

    //Find the 3 pair of points with the min_axis, max_axis values
    Vector3 min_x{ points.front() }, min_y{ points.front() }, min_z{ points.front() };
    Vector3 max_x{ points.front() }, max_y{ points.front() }, max_z{ points.front() };
    //iterate through the vector of points, if their x,y,z values are greater/lesser then make them the new min/max 
    for (const Vector3& p : points)
    {
        if (min_x.x > p.x)
            min_x = p;
        if (min_y.y > p.y)
            min_y = p;
        if (min_z.z > p.z)
            min_z = p;
        if (max_x.x < p.x)
            max_x = p;
        if (max_y.y < p.y)
            max_y = p;
        if (max_z.z < p.z)
            max_z = p;
    }


    //Compute the 'span/length' of each pair, and take the longest span
    float span_x = (max_x - min_x).Length();
    float span_y = (max_y - min_y).Length();
    float span_z = (max_z - min_z).Length();
    float max_span = span_x;
    if (max_span < span_y)
    {
        if (span_y < span_z)
        { 
            //span_z is the max span , so Ritter's center will be in between max_z and min_z 
            max_span = span_z;
            mCenter = 0.5f * (max_z + min_z);
        }
        else
        {
            //span_y is the max span, so Ritter's center will be in between max_y and min_y
            max_span = span_y;
            mCenter = 0.5f * (max_y + min_y);
        }
    }
    else if (max_span < span_z)
    {
        //span_z is the max span 
        max_span = span_z;
        mCenter = 0.5f * (max_z + min_z);
    }
    else
    {
        //span_x is the max span
        mCenter = 0.5f * (max_x + min_x);
    }
    //Use the half length of the chosen pair to be the radius
    mRadius = 0.5f * max_span;
    //run the sphere expansion algorithm 
    ExpandSphere(points);
    return;
}

void Sphere::ComputePCA(const std::vector<Vector3>& points)
{
  // The PCA method:
  // Compute the eigen values and vectors. Take the largest eigen vector as the axis of largest spread.
  // Compute the sphere center as the center of this axis then expand by all points.
  /******Student:Assignment2******/
  //Warn("Assignment2: Required function un-implemented");
}


bool Sphere::ContainsPoint(const Vector3& point)
{
  return PointSphere(point, mCenter, mRadius);
}

Vector3 Sphere::GetCenter() const
{
  return mCenter;
}

float Sphere::GetRadius() const
{
  return mRadius;
}

bool Sphere::Compare(const Sphere& rhs, float epsilon) const
{
  float posDiff = Math::Length(mCenter - rhs.mCenter);
  float radiusDiff = Math::Abs(mRadius - rhs.mRadius);

  return posDiff < epsilon && radiusDiff < epsilon;
}

DebugShape& Sphere::DebugDraw() const
{
  return gDebugDrawer->DrawSphere(*this);
}

// my own helper function to expand the sphere base on a vector of points
void Sphere::ExpandSphere(const std::vector<Vector3>& points)
{
    for (const Vector3& p : points)
    {
        if ((p - mCenter).LengthSq() > mRadius * mRadius)
        {
            //expand the sphere 
            Vector3 b = mCenter - mRadius * (p - mCenter).Normalized();
            mCenter = 0.5f * (p + b);
            mRadius = (p - mCenter).Length();
        }
    }
    return;
}

//-----------------------------------------------------------------------------Aabb
Aabb::Aabb()
{
  //set the aabb to an initial bad value (where the min is smaller than the max)
  mMin.Splat(Math::PositiveMax());
  mMax.Splat(-Math::PositiveMax());
}

Aabb::Aabb(const Vector3& min, const Vector3& max)
{
  mMin = min;
  mMax = max;
}

Aabb Aabb::BuildFromCenterAndHalfExtents(const Vector3& center, const Vector3& halfExtents)
{
  return Aabb(center - halfExtents, center + halfExtents);
}

float Aabb::GetVolume() const
{
  /******Student:Assignment2******/
  // Return the aabb's volume
  //Warn("Assignment2: Required function un-implemented");
    Vector3 delta = mMax - mMin;
    //volume should always be positive
    return Math::Abs(delta.x * delta.y * delta.z); 
}

float Aabb::GetSurfaceArea() const
{
  /******Student:Assignment2******/
  // Return the aabb's surface area
  //Warn("Assignment2: Required function un-implemented");
    Vector3 delta = mMax - mMin;
    return (delta.x * delta.y * 2.f) + (delta.x * delta.z * 2.f) + (delta.y * delta.z * 2.f);
}

bool Aabb::Contains(const Aabb& aabb) const
{
  /******Student:Assignment2******/
  // Return if aabb is completely contained in this
  //Warn("Assignment2: Required function un-implemented");
    //Check if the rhs aabb's min point is greater while max point is smaller than this aabb, then it will be inside this aabb
    if (aabb.mMin.x > mMin.x && aabb.mMax.x < mMax.x &&
        aabb.mMin.y > mMin.y && aabb.mMax.y < mMax.y &&
        aabb.mMin.z > mMin.z && aabb.mMax.z < mMax.z)
        return true;
  return false;
}

void Aabb::Expand(const Vector3& point)
{
  for(size_t i = 0; i < 3; ++i)
  {
    mMin[i] = Math::Min(mMin[i], point[i]);
    mMax[i] = Math::Max(mMax[i], point[i]);
  }
}

Aabb Aabb::Combine(const Aabb& lhs, const Aabb& rhs)
{
  Aabb result;
  for(size_t i = 0; i < 3; ++i)
  {
    result.mMin[i] = Math::Min(lhs.mMin[i], rhs.mMin[i]);
    result.mMax[i] = Math::Max(lhs.mMax[i], rhs.mMax[i]);
  }
  return result;
}

bool Aabb::Compare(const Aabb& rhs, float epsilon) const
{
  float pos1Diff = Math::Length(mMin - rhs.mMin);
  float pos2Diff = Math::Length(mMax - rhs.mMax);

  return pos1Diff < epsilon && pos2Diff < epsilon;
}

void Aabb::Transform(const Vector3& scale, const Matrix3& rotation, const Vector3& translation)
{
  /******Student:Assignment2******/
  // Compute aabb of the this aabb after it is transformed.
  // You should use the optimize method discussed in class (not transforming all 8 points).
  //Warn("Assignment2: Required function un-implemented");
    Vector3 new_center{GetCenter().x * scale.x ,
                       GetCenter().y* scale.y ,
                       GetCenter().z* scale.z };
    new_center = Math::Transform(rotation, new_center) + translation;
    Vector3 e = GetHalfSize(); //let 'e' be the old half-extent
    //Calculate the new halfextent with the optimized method
    Vector3 new_halfextent{ Math::Abs(rotation.m00 * e.x * scale.x) + Math::Abs(rotation.m01 * e.y* scale.y) + Math::Abs(rotation.m02 * e.z* scale.z) ,
                            Math::Abs(rotation.m10 * e.x * scale.x) + Math::Abs(rotation.m11 * e.y* scale.y) + Math::Abs(rotation.m12 * e.z* scale.z) ,
                            Math::Abs(rotation.m20 * e.x * scale.x) + Math::Abs(rotation.m21 * e.y* scale.y) + Math::Abs(rotation.m22 * e.z* scale.z) };
    
    //Update the new max and min values
    mMax = new_center + new_halfextent;
    mMin = new_center - new_halfextent;
}

Vector3 Aabb::GetMin() const
{
  return mMin;
}

Vector3 Aabb::GetMax() const
{
  return mMax;
}

Vector3 Aabb::GetCenter() const
{
  return (mMin + mMax) * 0.5f;
}

Vector3 Aabb::GetHalfSize() const
{
  return (mMax - mMin) * 0.5f;
}

DebugShape& Aabb::DebugDraw() const
{
  return gDebugDrawer->DrawAabb(*this);
}

//-----------------------------------------------------------------------------Triangle
Triangle::Triangle()
{
  mPoints[0] = mPoints[1] = mPoints[2] = Vector3::cZero;
}

Triangle::Triangle(const Vector3& p0, const Vector3& p1, const Vector3& p2)
{
  mPoints[0] = p0;
  mPoints[1] = p1;
  mPoints[2] = p2;
}

DebugShape& Triangle::DebugDraw() const
{
  return gDebugDrawer->DrawTriangle(*this);
}

//-----------------------------------------------------------------------------Plane
Plane::Plane()
{
  mData = Vector4::cZero;
}

Plane::Plane(const Vector3& p0, const Vector3& p1, const Vector3& p2)
{
  Set(p0, p1, p2);
}

Plane::Plane(const Vector3& normal, const Vector3& point)
{
  Set(normal, point);
}

void Plane::Set(const Vector3& p0, const Vector3& p1, const Vector3& p2)
{
  /******Student:Assignment1******/
  // Set mData from the 3 points. Note: You should most likely normalize the plane normal.
  //Warn("Assignment1: Required function un-implemented");
    //find the normal of the plane by crossing (P1-P0)x(P2-P0)
    Vector3 v01 = p1 - p0;
    Vector3 v02 = p2 - p0;
    Vector3 n = v01.Cross(v02);
    n.Normalize(); //ensure the plane's normal is normalized 
    // set the plane's normal
    mData.x = n.x;
    mData.y = n.y;
    mData.z = n.z;
    //Find the 'd' value (the distance from the origin) which can be computed as Dot(origin - pointOnPlane, normal)
    // d = n.p0
    mData.w = n.Dot(p0);
}

void Plane::Set(const Vector3& normal, const Vector3& point)
{
  /******Student:Assignment1******/
  // Set mData from the normal and point. Note: You should most likely normalize the plane normal.
  //Warn("Assignment1: Required function un-implemented");
    //ensure the normal is normalized
    Vector3 n = normal.Normalized();
    // set the plane's normal
    mData.x = n.x;
    mData.y = n.y;
    mData.z = n.z;
    //Find the 'd' value (the distance from the origin) which can be computed as Dot(origin - pointOnPlane, normal)
    // d = n.p0
    mData.w = n.Dot(point);
}

Vector3 Plane::GetNormal() const
{
  return Vector3(mData.x, mData.y, mData.z);
}

float Plane::GetDistance() const
{
  return mData.w;
}

DebugShape& Plane::DebugDraw(float size) const
{
  return DebugDraw(size, size);
}

DebugShape& Plane::DebugDraw(float sizeX, float sizeY) const
{
  return gDebugDrawer->DrawPlane(*this, sizeX, sizeY);
}

//-----------------------------------------------------------------------------Frustum
void Frustum::Set(const Vector3& lbn, const Vector3& rbn, const Vector3& rtn, const Vector3& ltn,
                  const Vector3& lbf, const Vector3& rbf, const Vector3& rtf, const Vector3& ltf)
{
  mPoints[0] = lbn;
  mPoints[1] = rbn;
  mPoints[2] = rtn;
  mPoints[3] = ltn;
  mPoints[4] = lbf;
  mPoints[5] = rbf;
  mPoints[6] = rtf;
  mPoints[7] = ltf;

  //left
  mPlanes[0].Set(lbf, ltf, lbn);
  //right
  mPlanes[1].Set(rbn, rtf, rbf);
  //top
  mPlanes[2].Set(ltn, ltf, rtn);
  //bot
  mPlanes[3].Set(rbn, lbf, lbn);
  //near
  mPlanes[4].Set(lbn, ltn, rbn);
  //far
  mPlanes[5].Set(rbf, rtf, lbf);
}

Math::Vector4* Frustum::GetPlanes() const
{
  return (Vector4*)mPlanes;
}

DebugShape& Frustum::DebugDraw() const
{
  return gDebugDrawer->DrawFrustum(*this);
}

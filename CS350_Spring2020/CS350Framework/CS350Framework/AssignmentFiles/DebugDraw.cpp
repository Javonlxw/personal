///////////////////////////////////////////////////////////////////////////////
///
/// Authors: Joshua Davis
/// Copyright 2015, DigiPen Institute of Technology
///
///////////////////////////////////////////////////////////////////////////////
#include "Precompiled.hpp"

#define ShowDebugDrawWarnings true

DebugDrawer* gDebugDrawer = new DebugDrawer();

//-----------------------------------------------------------------------------DebugShape
DebugShape::DebugShape()
{
  mColor = Vector4(.6f);
  mMask = (unsigned int)-1;
  mTimer = 0;
  mOnTop = false;
  mTransform.SetIdentity();
}

DebugShape& DebugShape::Color(const Vector4& color)
{
  mColor = color;
  return *this;
}

DebugShape& DebugShape::OnTop(bool state)
{
  mOnTop = state;
  return *this;
}

DebugShape& DebugShape::Time(float time)
{
  mTimer = time;
  return *this;
}

DebugShape& DebugShape::SetMaskBit(int bitIndex)
{
  mMask = 1 << bitIndex;
  return *this;
}

DebugShape& DebugShape::SetTransform(const Matrix4& transform)
{
  mTransform = transform;
  return *this;
}

//-----------------------------------------------------------------------------DebugDrawer
DebugDrawer::DebugDrawer()
{
  mActiveMask = (unsigned int)-1;
  mApplication = NULL;
}

void DebugDrawer::Update(float dt)
{
  std::vector<DebugShape> newShapes;
  for(size_t i = 0; i < mShapes.size(); ++i)
  {
    DebugShape& shape = mShapes[i];
    shape.mTimer -= dt;

    // If the shape still has time left then add it to the list of shapes to keep drawing,
    // anything that has a timer that ran out will not be in the new list
    if(shape.mTimer >= 0)
      newShapes.push_back(shape);
  }

  mShapes.swap(newShapes);
}

void DebugDrawer::Draw()
{
  for(size_t i = 0; i < mShapes.size(); ++i)
  {
    DebugShape& shape = mShapes[i];

    // If the shape doesn't have one of the active mask bits set then don't draw it
    if((shape.mMask & mActiveMask) == 0)
      continue;
    
    // If this shape always draws on top then disable depth testing
    if(shape.mOnTop)
      glDisable(GL_DEPTH_TEST);


    // Decompose the matrix to set the gl transform (too lazy to properly transform the matrix between formats)
    float radians;
    Vector3 scale, translation, axis;
    Matrix3 rotationMat;
    shape.mTransform.Decompose(&scale, &rotationMat, &translation);
    Math::ToAxisAngle(Math::ToQuaternion(rotationMat), &axis, &radians);
    glPushMatrix();
    // Set the transform
    glTranslatef(translation.x, translation.y, translation.z);
    glRotatef(Math::RadToDeg(radians), axis.x, axis.y, axis.z);
    glScalef(scale.x, scale.y, scale.z);

    glBegin(GL_LINES);
    glColor3fv(shape.mColor.array);

    // Draw all of the line segments of this shape
    for(size_t j = 0; j < shape.mSegments.size(); ++j)
    {
      LineSegment& segment = shape.mSegments[j];

      glVertex3fv(segment.mStart.array);
      glVertex3fv(segment.mEnd.array);
    }

    glEnd();
    glPopMatrix();

    // Make sure to re-enable depth testing
    if(shape.mOnTop)
      glEnable(GL_DEPTH_TEST);
  }
}

DebugShape& DebugDrawer::GetNewShape()
{
  mShapes.push_back(DebugShape());
  return mShapes.back();
}

DebugShape& DebugDrawer::DrawPoint(const Vector3& point)
{
  return DrawSphere(Sphere(point, 0.1f));
}

DebugShape& DebugDrawer::DrawLine(const LineSegment& line)
{
  /******Student:Assignment2******/
  // Draw a simple line
  DebugShape& shape = GetNewShape();
  //WarnIf(ShowDebugDrawWarnings, "Assignment2: Required function un-implemented");
  shape.mSegments.push_back(line);
  return shape;
}

DebugShape& DebugDrawer::DrawRay(const Ray& ray, float t)
{
  /******Student:Assignment2******/
  // Draw a ray to a given t-length. The ray must have an arrow head for visualization
  //WarnIf(ShowDebugDrawWarnings, "Assignment2: Required function un-implemented");
  DebugShape& shape = GetNewShape();
  CalculateRay(ray, t, shape.mSegments);
  return shape;
}

DebugShape& DebugDrawer::DrawSphere(const Sphere& sphere)
{
  /******Student:Assignment2******/
  // Draw a sphere with 4 rings: x-axis, y-axis, z-axis, and the horizon disc.
  // Note: To access the camera's position for the horizon disc calculation use mApplication->mCamera.mTranslation
  DebugShape& shape = GetNewShape();
  //Safety check to prevent seg fault 
  if (!mApplication)
      return shape;

  //Get the horizontal disc's center and radius
  Vector3 HD_center;
  float HD_rad;
  Vector3 x_unit{ 1.f, 0.f, 0.f }, y_unit{ 0.f,1.f,0.f }, z_unit{ 0.f,0.f,1.f };
  //Calculate the line segments to draw on the z-axis, therefore we use the z_unit vector as the axis vector 
  CalculateDisc(sphere.mCenter, sphere.mRadius, z_unit, shape.mSegments);
  //Calculate the line segments to draw on the x-axis
  CalculateDisc(sphere.mCenter, sphere.mRadius, x_unit, shape.mSegments);
  //Calculate the line segments to draw on the y-axis
  CalculateDisc(sphere.mCenter, sphere.mRadius, y_unit, shape.mSegments);
  
  //Compute the horizontal disc's center and radius
  GetHorizontalDisc(sphere, HD_center, HD_rad);
  //Calculate the 'axis' that the horizontal disc will be on.
  Vector3 HD_axis = sphere.mCenter - mApplication->mCamera.mTranslation;
  // Lastly, calculate the line segments to draw the horizontal disc
  CalculateDisc(HD_center, HD_rad, HD_axis.Normalized(), shape.mSegments);
  //WarnIf(ShowDebugDrawWarnings, "Assignment2: Required function un-implemented");
  return shape;
}

DebugShape& DebugDrawer::DrawAabb(const Aabb& aabb)
{
  /******Student:Assignment2******/
  // Draw all edges of an aabb. Make sure to not mis-match edges!
  DebugShape& shape = GetNewShape();  
  // Use the aabb's Min and Max to get the 6 other points 
  // Note that Max is the front_top_right point, while Min is back_bottom_left point 
  Vector3 delta = aabb.GetMax() - aabb.GetMin();
  Vector3 front_top_left    = aabb.mMax - Vector3{ delta.x,0.f,0.f };
  Vector3 front_btm_left    = front_top_left    -  Vector3{ 0.f, delta.y, 0.f };
  Vector3 front_btm_right   = aabb.GetMax()     -  Vector3{ 0.f, delta.y, 0.f };
  Vector3 back_top_left     = aabb.GetMin()     +  Vector3{ 0.f, delta.y, 0.f };
  Vector3 back_top_right    = back_top_left     +  Vector3{ delta.x, 0.f,0.f };
  Vector3 back_btm_right    = back_top_right    -  Vector3{ 0.f, delta.y, 0.f };

  //Calculate the lines for the front face (quad)
  CalculateQuad(front_top_left, aabb.GetMax(), front_btm_right, front_btm_left, shape.mSegments);
  //Calcualte the lines for the back face (quad)
  CalculateQuad(back_top_left, back_top_right, back_btm_right, aabb.GetMin(), shape.mSegments);
 
  //Push the right-side top line
  shape.mSegments.push_back(LineSegment{ aabb.GetMax(), back_top_right });
  //Push the right-side bottom line
  shape.mSegments.push_back(LineSegment{ front_btm_right, back_btm_right});
  //Push the left-side top line
  shape.mSegments.push_back(LineSegment{ front_top_left, back_top_left });
  //Push the left-side bottom line
  shape.mSegments.push_back(LineSegment{ front_btm_left, aabb.GetMin() });
  //WarnIf(ShowDebugDrawWarnings, "Assignment2: Required function un-implemented");
  return shape;
}

DebugShape& DebugDrawer::DrawTriangle(const Triangle& triangle)
{
  /******Student:Assignment2******/
  // Draw the 3 edges of a triangles
  DebugShape& shape = GetNewShape();
  //WarnIf(ShowDebugDrawWarnings, "Assignment2: Required function un-implemented");
  for (int i = 0; i < 3; ++i)
  {
      if (i != 2)
          shape.mSegments.push_back(LineSegment{ triangle.mPoints[i], triangle.mPoints[i + 1] });
      else
          shape.mSegments.push_back(LineSegment{ triangle.mPoints[0], triangle.mPoints[i] });
  }
  return shape;
}

DebugShape& DebugDrawer::DrawPlane(const Plane& plane, float sizeX, float sizeY)
{
  /******Student:Assignment2******/
  // Draw a quad with a normal at the plane's center.
  DebugShape& shape = GetNewShape();
  //Draw the Normal using DrawRay
  Vector3 p0 = plane.GetDistance() * plane.GetNormal(); // calculate mStart of the ray
  Ray m_ray{ p0, plane.GetNormal() };
  //DrawRay(m_ray, 1.0f);
  CalculateRay(m_ray, 2.0f, shape.mSegments);
  //Calculate p1,2,3,4 that'll be used to draw the Quad that represents the plane 
  Vector3 v, w;
  CalculateOrthogonalBasis(m_ray.mDirection, v, w);
  Vector3 p1 = p0 + sizeX * v + sizeY * w;
  Vector3 p2 = p0 - sizeX * v + sizeY * w;
  Vector3 p3 = p0 - sizeX * v - sizeY * w;
  Vector3 p4 = p0 + sizeX * v - sizeY * w;
  //Calculate the lines needed to draw the quad that represents the plane
  CalculateQuad(p1, p2, p3, p4, shape.mSegments);
  //WarnIf(ShowDebugDrawWarnings, "Assignment2: Required function un-implemented");
  return shape;
}

DebugShape& DebugDrawer::DrawQuad(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3)
{
  /******Student:Assignment2******/
  // Draw the4 edges of a quad. Make sure to look at this and make sure the quad is not bow-tied.
  DebugShape& shape = GetNewShape();
  CalculateQuad(p0,p1,p2,p3,shape.mSegments);
//  WarnIf(ShowDebugDrawWarnings, "Assignment2: Required function un-implemented");
  return shape;
}

DebugShape& DebugDrawer::DrawFrustum(const Frustum& frustum)
{
  /******Student:Assignment2******/
  // Draw the 6 faces of the frustum using the 8 frustum points.
  // See Frustum.Set for the point order. For example, Points[4] is left-bottom-front.
  DebugShape& shape = GetNewShape();
  //For easy reference
  /*
  mPoints[0] = Left-Bottom-Near
  mPoints[1] = Right-Bottom_Near
  mPoints[2] = Right-Top-Near
  mPoints[3] = Left-Top-Near

  mPoints[4] = Left-Bottom-Far
  mPoints[5] = Right-Bottom-Far
  mPoints[6] = Right-Top-Far
  mPointsp7[ = Left-Top-Far
  */
  //Calculate the lines to draw the near-face (CCW)
  CalculateQuad(frustum.mPoints[0], frustum.mPoints[1], frustum.mPoints[2], frustum.mPoints[3], shape.mSegments);
  //Calculate the lines to draw the far-face (CCW)
  CalculateQuad(frustum.mPoints[4], frustum.mPoints[5], frustum.mPoints[6], frustum.mPoints[7], shape.mSegments);

  //Calculate the lines to draw the right-side
  shape.mSegments.push_back(LineSegment{ frustum.mPoints[2], frustum.mPoints[6] });
  shape.mSegments.push_back(LineSegment{ frustum.mPoints[1], frustum.mPoints[5] });
  //Calculate the lines to draw the left-side
  shape.mSegments.push_back(LineSegment{ frustum.mPoints[3], frustum.mPoints[7] });
  shape.mSegments.push_back(LineSegment{ frustum.mPoints[0], frustum.mPoints[4] });
  //WarnIf(ShowDebugDrawWarnings, "Assignment2: Required function un-implemented");
  return shape;
}

void DebugDrawer::CalculateOrthogonalBasis(const Vector3& u, Vector3& v, Vector3& w)
{
    v = u.Cross(Vector3{ 0.f,1.f,0.f });
    //edge case where the ray is parallel to the global up axis
    if (v.LengthSq() == 0.f)
    {
        v = Vector3{ 1.0f, 0.f,0.f };// set v to be parallel to the global right
        w = Vector3{ 0.f, 0.f, 1.f };// set w to be prallel to the global forward 
        return;
    }
    w = v.Cross(u);

    //Normalize v and w
    v = v.Normalized();
    w = w.Normalized();

    
    return;
}

void DebugDrawer::CalculateRay(const Ray& ray, float t, std::vector<LineSegment>& shape_lines)
{
    //Let p0 be the end of the ray 
    Vector3 p0 = ray.mStart + ray.mDirection * t;
    LineSegment ray_line{ ray.mStart, p0 };
    shape_lines.push_back(ray_line);
    //Calculate the lines needed to draw the arrow-head //////////////////////
    float b = 0.5f; //set my own 'b' value to calcualte p0_prime 
    Vector3 p0_prime = p0 - (b * ray.mDirection); // Let p0' be the center of the disc for the arrow head 
    float r = 0.25f; // set my own radius for the arrow head's disc 
    int n = 100; // 'n' refers to the level of detail for the arrow-head, or the number of line segments used to draw the arrow disc
    float theta = 0.0f; // theta is the angle used to find the next point to draw for the arrow disc 
    float theta_delta = (Math::cTwoPi) / (float)n;
    Vector3 v, w;
    CalculateOrthogonalBasis(ray.mDirection, v, w);
    // Calculating the line segments that are used to draw the arrow disc 
    for (int i = 0; i < n; ++i)
    {
        Vector3 p_first = p0_prime + r * Math::Cos(theta) * v + r * Math::Sin(theta) * w;
        shape_lines.push_back(LineSegment{ p_first, p0 }); //line to create the 'slanted surface' of the arrowhead 
        theta += theta_delta;
        Vector3 p_second = p0_prime + r * Math::Cos(theta) * v + r * Math::Sin(theta) * w;
        // shape.mSegments.push_back(LineSegment{ p_second,p0 }); //line to create the 'cone' of the arrowhead 
        shape_lines.push_back(LineSegment{ p_first,p_second }); // line to draw the 'disc/base' of the arrowhead
    }

    return;
}

void DebugDrawer::CalculateQuad(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3, std::vector<LineSegment>& shape_lines)
{
    // Calculate line p1-p0
    shape_lines.push_back(LineSegment{ p0,p1 });
    // Calculate line p2-p1
    shape_lines.push_back(LineSegment{ p1, p2 });
    // Calculate line p3-p2
    shape_lines.push_back(LineSegment{ p2, p3 });
    // Calculate line p0-p3
    shape_lines.push_back(LineSegment{ p3, p0 });
}

void DebugDrawer::GetHorizontalDisc(const Sphere& sphere, Vector3& HD_center, float& HD_radius)
{
    //Get the camera pos, aka 'e', the eye of the camera
    Vector3 e = mApplication->mCamera.mTranslation;
    //Calculate the vector from e to the center of the sphere
    Vector3 CE = sphere.GetCenter() - e;
    //Calculate 'd_sq' which is the distance squared of the center of the sphere to the camera
    //float d_sq = CE.LengthSq();
    ////Get r^2
    //float r_sq = sphere.GetRadius() * sphere.GetRadius();
    ////Calculate 'l_sq' which is the hypo squared of the triangle d-r-l
    //float l_sq = d_sq - r_sq; 
    ////From r' = (r*l)/d, we know r'*r' = (r*r*l*l)/d*d 
    //float r_prime_sq = (r_sq * l_sq) / d_sq;
    //HD_radius = Math::Sqrt(r_prime_sq);
    ////calculate 'z' which is the distance between the sphere and the horizontal sphere's center.
    //float z = Math::Sqrt(r_sq - r_prime_sq );
    ////Finally calculate the center of the horizontal disc by: c' = c - z*CE
    //HD_center = sphere.GetCenter() - z * CE.Normalized();

    float d = CE.Length();
    float l = Math::Sqrt(d * d - sphere.mRadius * sphere.mRadius);
    HD_radius = (sphere.mRadius * l) / d;
    float z = Math::Sqrt(sphere.mRadius * sphere.mRadius - HD_radius * HD_radius);
    HD_center = sphere.mCenter - z * CE.Normalized();

    return;
}

void DebugDrawer::CalculateDisc(const Vector3& center, const float& radius, const Vector3& axis, std::vector<LineSegment>& shape_lines)
{
    //Calculate the orthogonal vectors to the axis
    Vector3 v, w;
    CalculateOrthogonalBasis(axis, v, w);
    int n = 100; // 'n' refers to the level of detail for the arrow-head, or the number of line segments used to draw the arrow disc
    float theta = 0.0f; // theta is the angle used to find the next point to draw for the arrow disc 
    float theta_delta = (Math::cTwoPi) / (float)n;
    // Calculating the line segments that are used to draw the arrow disc 
    for (int i = 0; i < n; ++i)
    {
        Vector3 p_first = center + radius * Math::Cos(theta) * v + radius * Math::Sin(theta) * w;
        theta += theta_delta;
        Vector3 p_second = center + radius * Math::Cos(theta) * v + radius * Math::Sin(theta) * w;
        shape_lines.push_back(LineSegment{ p_first,p_second }); // line to draw the 'disc/base' of the arrowhead
    }
}
/*******************************************************************************
 *
 * Descritpion: Rectangle class
 *
 ******************************************************************************/
#include "VolumeRectangle.h"

#include <math.h>

using namespace ParticlesFillaar;

//------------------------------------------------------------------------------
// Constructors / Destructor
VolumeRectangle::VolumeRectangle(const Vec3d& minCoords, const Vec3d& maxCoords)
: _minCoords(minCoords),
	_maxCoords(maxCoords)
{
	//_maxCoords.z = _minCoords.z;	// We assume the rectangle is aligned on the XY plane
}

VolumeRectangle::VolumeRectangle(const Vec3d& center, double dimX, double dimY)
{
	double hDimX = dimX / 2.0;
	double hDimY = dimY / 2.0;
	
	_minCoords.x = center.x - hDimX;
	_maxCoords.x = center.x + hDimX;
	_minCoords.y = center.y - hDimY;
	_maxCoords.y = center.y + hDimY;
	
	_maxCoords.z = _minCoords.z;	// We assume the rectangle is aligned on the XY plane
}

VolumeRectangle::~VolumeRectangle()
{
}

//------------------------------------------------------------------------------
// Public functions
double VolumeRectangle::projectPointOntoSurface(const Vec3d& p, Vec3d& projection, Vec3d& normal) const
{
	Vec3d center = _minCoords;
	center += _maxCoords;
	center /= 2.0;
	
	projection = p;
	projection.z = _minCoords.z;
	normal.z = 0.0;
	
	// Project onto the rectangle toward (or away from) its center
	// Project into lines x=minX, x=maxX, y=minY and y=maxY using the parametric equation of the line:
	// P' = P + (P-C)*t, where P is the point and C is the center of the rectangle
	bool intersectionFound = false;
	double tmin = 0.0;

	// x=minX
	double t = (_minCoords.x - p.x) / (p.x-center.x);
	double lineIntersection = p.y + (p.y - center.y) * t;
	if ((!intersectionFound || fabs(t)<tmin) && (lineIntersection >= _minCoords.y) && (lineIntersection <= _maxCoords.y))
	{
		intersectionFound = true;
		tmin = fabs(t);
		projection.x = _minCoords.x;
		projection.y = lineIntersection;
		normal.x = -1.0;
		normal.y = 0.0;
	}

	// x=maxX
	t = (_maxCoords.x - p.x) / (p.x-center.x);
	lineIntersection = p.y + (p.y - center.y) * t;
	if ((!intersectionFound || fabs(t)<tmin) && (lineIntersection >= _minCoords.y) && (lineIntersection <= _maxCoords.y))
	{
		intersectionFound = true;
		tmin = fabs(t);
		projection.x = _maxCoords.x;
		projection.y = lineIntersection;
		normal.x = 1.0;
		normal.z = 0.0;
	}

	// y=minY
	t = (_minCoords.y - p.y) / (p.y-center.y);
	lineIntersection = p.x + (p.x - center.x) * t;
	if ((!intersectionFound || fabs(t)<tmin) && (lineIntersection >= _minCoords.x) && (lineIntersection <= _maxCoords.x))
	{
		intersectionFound = true;
		tmin = fabs(t);
		projection.x = lineIntersection;
		projection.y = _minCoords.y;
		normal.x = 0.0;
		normal.y = -1.0;
	}

	// y =maxY
	t = (_maxCoords.y - p.y) / (p.y-center.y);
	lineIntersection = p.x + (p.x - center.x) * t;
	if ((!intersectionFound || fabs(t)<tmin) && (lineIntersection >= _minCoords.x) && (lineIntersection <= _maxCoords.x))
	{
		intersectionFound = true;
		tmin = fabs(t);
		projection.x = lineIntersection;
		projection.y = _maxCoords.y;
		normal.x = 0.0;
		normal.y = 1.0;
	}

	// Return the length of the projection
	Vec3d delta = p;
	delta -= projection;
	double sign = pointIsInside(p) ? -1.0: 1.0;	// Inside will be negative distance, outside is positive
	
	return sign*delta.length();
}

bool VolumeRectangle::pointIsInside(const Vec3d& p) const
{
	bool isInside = false;
	
	if ((p.x >= _minCoords.x) && (p.x <= _maxCoords.x) &&
			(p.y >= _minCoords.y) && (p.y <= _maxCoords.y))
	{
		isInside = true;
	}
	
	return isInside;
}

void VolumeRectangle::getAABB(Vec3d& aabbMin, Vec3d& aabbMax) const
{
	aabbMin = _minCoords;
	aabbMax = _maxCoords;
}

Vec3d VolumeRectangle::getmin() const
{
	return _minCoords;
}

/*******************************************************************************
 *
 * Description: Fillar class
 *
 ******************************************************************************/
#include "ParticlesFillaar.h"
#include "Volume.h"

using namespace ParticlesFillaar;

//------------------------------------------------------------------------------
// Constructor / Destructor
Fillaar::Fillaar()
{}

Fillaar::~Fillaar()
{
	clear();	
}

//------------------------------------------------------------------------------
// Public functions
std::deque< Vec3d >& Fillaar::fill(double resolution, int xyz)
{
	// Pre-filling operations
	preFill(resolution,xyz);
	
	// Fill volumes
	std::vector<Volume*>::iterator it = _volumes.begin();
	std::vector<Volume*>::iterator itEnd = _volumes.end();
	for (; it != itEnd; ++it)
	{
		fillVolume(*it);
	}
	
	// Post-filling operations
	postFill();

	return _points;
}

void Fillaar::addVolume(Volume* v)
{
	_volumes.push_back(v);
}

void Fillaar::removeVolume(Volume* v)
{
	std::vector<Volume*>::iterator it = _volumes.begin();
	std::vector<Volume*>::iterator itEnd = _volumes.end();
	for (; it != itEnd; ++it)
	{
		if ((*it) == v)
		{
			_volumes.erase(it);
			break;
		}
	}
}

void Fillaar::deleteVolume(Volume* v)
{
	removeVolume(v);
	delete v;
}

void Fillaar::getVolumesAABB(Vec3d& aabbMin, Vec3d& aabbMax) const
{
	Vec3d volumeAABBMin;
	Vec3d volumeAABBMax;
	
	aabbMin = Vec3d(0,0,0);
	aabbMax = Vec3d(0,0,0);
	
	std::vector<Volume*>::const_iterator it = _volumes.begin();
	std::vector<Volume*>::const_iterator itEnd = _volumes.end();
	bool isFirst = true;
	for (; it != itEnd; ++it)
	{
		Volume* v = *it;
		
		if (isFirst)
		{
			isFirst = false;
			v->getAABB(aabbMin, aabbMax);
		}
		else
		{
			v->getAABB(volumeAABBMin, volumeAABBMax);
			
			if (volumeAABBMin.x < aabbMin.x)
			{
				aabbMin.x = volumeAABBMin.x;
			}
			if (volumeAABBMax.x > aabbMax.x)
			{
				aabbMax.x = volumeAABBMax.x;
			}
			
			if (volumeAABBMin.y < aabbMin.y)
			{
				aabbMin.y = volumeAABBMin.y;
			}
			if (volumeAABBMax.y > aabbMax.y)
			{
				aabbMax.y = volumeAABBMax.y;
			}
			
			if (volumeAABBMin.z < aabbMin.z)
			{
				aabbMin.z = volumeAABBMin.z;
			}
			if (volumeAABBMax.z > aabbMax.z)
			{
				aabbMax.z = volumeAABBMax.z;
			}
		}
	}
}

void Fillaar::clear()
{
	clearVolumes();
	clearPoints();
}

void Fillaar::clearVolumes()
{
	std::vector<Volume*>::iterator it = _volumes.begin();
	std::vector<Volume*>::iterator itEnd = _volumes.end();
	for (; it != itEnd; ++it)
	{
		Volume* v = *it;
		delete v;
	}
	
	_volumes.clear();
}

void Fillaar::clearPoints()
{
	_points.clear();
}

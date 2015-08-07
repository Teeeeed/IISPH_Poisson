/*******************************************************************************
 *
 * Descritpion: Rectangle class
 *
 ******************************************************************************/
#ifndef PARTICLESFILLAAR_VOLUMERECTANGLE_H
#define PARTICLESFILLAAR_VOLUMERECTANGLE_H

#include "Volume.h"


namespace ParticlesFillaar {

class VolumeRectangle : public Volume
{
public:
    VolumeRectangle(const Vec3d& minCoords, const Vec3d& maxCoords);
    VolumeRectangle(const Vec3d& center, double dimX, double dimY);
    virtual ~VolumeRectangle();

    virtual double projectPointOntoSurface(const Vec3d& p, Vec3d& projection, Vec3d& normal) const;
    virtual bool pointIsInside(const Vec3d& p) const;
    virtual void getAABB(Vec3d& aabbMin, Vec3d& aabbMax) const;
	virtual Vec3d getmin() const;

private:
    Vec3d	_minCoords;
    Vec3d _maxCoords;
};

}

#endif // PARTICLESFILLAAR_VOLUMERECTANGLE_H

/*******************************************************************************
 *
 * Descritpion: Volume class
 *
 ******************************************************************************/
#ifndef VOLUME_H
#define VOLUME_H

#include "ParticlesFillaar.h"
#include "Vec3.h"

namespace ParticlesFillaar
{
class Volume
{

public:
        Volume() {}
        virtual ~Volume() {}

        virtual void getAABB(Vec3d& aabbMin, Vec3d& aabbMax) const = 0;
        virtual bool pointIsInside(const Vec3d& p) const = 0;
        virtual double projectPointOntoSurface(const Vec3d& p, Vec3d& projection, Vec3d& normal) const = 0;
		virtual Vec3d getmin() const = 0;
};
}

#endif // VOLUME_H

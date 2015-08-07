/*******************************************************************************
 *
 * Description: Fillar class
 *
 ******************************************************************************/
#ifndef PARTICLESFILLAAR_H
#define PARTICLESFILLAAR_H

#include "Vec3.h"

#include <vector>
#include <deque>

namespace ParticlesFillaar
{
    class Volume;

    class Fillaar
    {
    public:
        Fillaar();
        virtual ~Fillaar();

        virtual std::deque<Vec3d>& fill(double resolution,int xyz);

        void addVolume(Volume* v);
        void removeVolume(Volume* v);
        void deleteVolume(Volume* v);

        void getVolumesAABB(Vec3d& aabbMin, Vec3d& aabbMax) const;

        void clear();
        void clearVolumes();
        void clearPoints();

        std::deque<Vec3d>& getPoints() { return _points; }
        const std::deque<Vec3d>& getPoints() const { return _points; }

        const std::vector<Volume*>& getVolumes() const { return _volumes; }

    protected:
        virtual void preFill(double resolution,int xyz) = 0;
        virtual void fillVolume(const Volume* v) = 0;
        virtual void postFill() = 0;

    private:
        std::deque<Vec3d>		_points;
        std::vector<Volume*>	_volumes;
    };
}

#endif	// PARTICLESFILLAAR_H

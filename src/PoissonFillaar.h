/*******************************************************************************
 *
 * Description:	Possion sampling
 *
 ******************************************************************************/
#ifndef PARTICLESFILLAAR_POISSONFILLAAR_H
#define PARTICLESFILLAAR_POISSONFILLAAR_H

#include "ParticlesFillaar.h"
#include <vector>

namespace ParticlesFillaar
{
    class PoissonFillaar : public ParticlesFillaar::Fillaar
    {
    public:
    PoissonFillaar();
    virtual ~PoissonFillaar();

    protected:
    virtual void postFill();
    virtual void fillVolume(const Volume* v);
    virtual void preFill(double resolution,int xyz);

    private:
        inline int getGridIndex(const Vec3d& p) const;
        inline int getGridIndex(const int xIndex, const int yIndex) const;
        inline int getGridXIndex(const int index) const;
        inline int getGridXIndex(const double x) const;
        inline int getGridYIndex(const int index) const;
        inline int getGridYIndex(const double y) const;

        inline void addPoint(const Vec3d& p);

        inline bool satisfyPoissonDiskCriterion(const Vec3d& p) const;

        void surfaceSampling(const Volume* v, std::vector<int>& surfaceSamples);
        void volumeSampling(const Volume*v, const std::vector<int>& initialSeeds, std::vector<int>& volumeSamples);
        void relaxation(const Volume* v, const std::vector<int>& pointsToRelax, int nbSweeps, bool alwaysProjectToSurface = false);

    private:
        double						_resolution;
        double						_cellSize;
        std::vector<int>	_grid;
        Vec3d						_gridMin;
        int								_gridResX;
        int								_gridResY;
        int								_gridSize;
		int								_xyz;
    };
}

#endif // PARTICLESFILLAAR_POISSONFILLAAR_H

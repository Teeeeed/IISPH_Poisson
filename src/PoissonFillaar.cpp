/*******************************************************************************
 *
 * Description:	Possion sampling based on the paper "Fast Poisson Disk 
 *				Sampling in Arbitrary Dimensions" (Bridson 2007). Only
 *				xOy-dimension poisson disk is implemented.
 *
 ******************************************************************************/
#include "PoissonFillaar.h"
#include "Volume.h"

#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace ParticlesFillaar;

//------------------------------------------------------------------------------
// Constructor / Destructor
PoissonFillaar::PoissonFillaar()
: _resolution(0.1),
	_grid(0x0),
	_xyz(0)
{
}

PoissonFillaar::~PoissonFillaar()
{
}

//------------------------------------------------------------------------------
// Overloaded functions from Fillaar
void PoissonFillaar::preFill(double resolution, int xyz)
{
	_resolution = 0.92*resolution;	// Value used in the article
	_cellSize = _resolution / sqrt(2.0);	// TODO: Make a 3D implementation and rename this class so that we know it's 2D or make this class able to do 2D AND 3D... 
	
	// Initialize grid
	Vec3d gridDim;
	Vec3d gridBorder(_cellSize, _cellSize,0.0);
	getVolumesAABB(_gridMin, gridDim);
	_gridMin -= gridBorder;
	gridDim += gridBorder;
	gridDim -= _gridMin;
	_gridResX = static_cast<int>(ceil(gridDim.x/ _cellSize));
	_gridResY = static_cast<int>(ceil(gridDim.y/ _cellSize));
	_gridSize = _gridResX * _gridResY;
	_grid.resize(_gridSize);
	
	for (int i=0; i<_gridSize; ++i)
	{
		_grid[i] = -1;
	}
	//
	//std::cout << "CellSize: " << _cellSize << std::endl;
	//std::cout << "GridMin: (" << _gridMin.x << ", " << _gridMin.y << ")" << std::endl;
	//std::cout << "GridRes: (" << _gridResX << ", " << _gridResY << ")" << std::endl;
	//std::cout << "GridDim: (" << gridDim.x << ", " << gridDim.y << ")" << std::endl;
}

void PoissonFillaar::postFill()
{
	// Clear the grid and free the memory! (If we simply call clear(), the memory won't be freed)
	std::vector<int>().swap(_grid);
}

void PoissonFillaar::fillVolume(const ParticlesFillaar::Volume* v)
{
	std::vector<int> generatedSurfaceSamples;
	std::vector<int> generatedVolumeSamples;
	
	surfaceSampling(v, generatedSurfaceSamples);
    relaxation(v, generatedSurfaceSamples, 5, true);
	
	volumeSampling(v, generatedSurfaceSamples, generatedVolumeSamples);
    relaxation(v, generatedVolumeSamples, 5, false);
}

//------------------------------------------------------------------------------
// Private functions
int PoissonFillaar::getGridIndex(const Vec3d& p) const
{
	// BIG WARNING!!!: We assume that the point lies inside the grid and do not check
	//                 if it's out of bound. It's ok since the grid is only used
	//								 internally and the grid has been constructed to surround
	//								 everything that will be filled.
	int ix = static_cast<int>(floor((p.x-_gridMin.x)/_cellSize));
	int iy = static_cast<int>(floor((p.y-_gridMin.y)/_cellSize));
	
	return (iy*_gridResX + ix);
}

int PoissonFillaar::getGridIndex(const int xIndex, const int yIndex) const
{
	// BIG WARNING!!!: We assume that the point lies inside the grid and do not check
	//                 if it's out of bound. It's ok since the grid is only used
	//								 internally and the grid has been constructed to surround
	//								 everything that will be filled.
	return (yIndex*_gridResX + xIndex);
}


int PoissonFillaar::getGridXIndex(const int index) const
{
	return index%_gridResX;
}

int PoissonFillaar::getGridXIndex(const double x) const
{
	return static_cast<int>(floor((x-_gridMin.x)/_cellSize));
}

int PoissonFillaar::getGridYIndex(const int index) const
{
	return index/_gridResX;
}

int PoissonFillaar::getGridYIndex(const double y) const
{
	return static_cast<int>(floor((y-_gridMin.y)/_cellSize));
}

bool PoissonFillaar::satisfyPoissonDiskCriterion(const Vec3d& p) const
{
	// Declarations
	double r = _resolution;
	double r2 = r*r;
	
	const std::deque<Vec3d>& points = getPoints();
	
	// Test cell containing point
	int pIndex = getGridIndex(p);
	bool doesNotSatisfyCriterion = false;
	if (_grid[pIndex] < 0)
	{
		// Test neighbor cells
		int xIndex = getGridXIndex(pIndex);
		int yIndex = getGridYIndex(pIndex);
		int startX = getGridXIndex(p.x-_resolution); if (startX < 0) { startX = 0; }
		int startY = getGridYIndex(p.y-_resolution); if (startY < 0) { startY = 0; }
		int endX = getGridXIndex(p.x+_resolution); if (endX >(_gridResX-1)) { endX = _gridResX-1; }
		int endY = getGridYIndex(p.y+_resolution); if (endY >(_gridResY-1)) { endY = _gridResY-1; }
		for (int ix=startX; !doesNotSatisfyCriterion && (ix<=endX); ++ix)
		{
			for (int iy=startY; !doesNotSatisfyCriterion && (iy<=endY); ++iy)
			{
				if (ix!=xIndex || iy!=yIndex)
				{
					int neighborIndex = getGridIndex(ix, iy);
					int neighborID = _grid[neighborIndex];
					if (neighborID > -1)
					{
						Vec3d delta = p;
						delta -= points[neighborID];
						
						if (delta.length2() < r2)
						{
							doesNotSatisfyCriterion = true;
						}
					}
				}
			}
		}
	}
	else
	{
		doesNotSatisfyCriterion = true;
	}
	
	return !doesNotSatisfyCriterion;
}

//------------------------------------------------------------------------------
// Sampling/Relaxation functions
void PoissonFillaar::surfaceSampling(const Volume* v, std::vector< int >& surfaceSamples)
{
	// Declarations
	double randMax = static_cast<double>(RAND_MAX);
	double epsilon = 1.085;
	double r = _resolution;
	
	std::deque<Vec3d>& points = getPoints();
	
	// First, fill the surface (See Algo 1 in Schechter2012)
	// 1) For each cells try to generate and project a random point
	//    that satisfy the poisson disk creterion
	for (int cell=0; cell<_gridSize; ++cell)
	{
		// Determine if its an interface cell
		// (Check if there is a random point inside and a 
		//  random point outside the volume inside the same cell)
		bool isSurfaceCell = false;
		/*bool hasOnePointInside = false;
		bool hasOnePointOutside = false;
		for (int k=0; !isSurfaceCell && (k<30); ++k)
		{
			Vec3d p;
			p.x = (static_cast<double>(rand())/randMax) * _cellSize + _gridMin.x + getGridXIndex(cell)*_cellSize;
			p.y = (static_cast<double>(rand())/randMax) * _cellSize + _gridMin.y + getGridYIndex(cell)*_cellSize;
			
			if (v->pointIsInside(p))
			{
				hasOnePointInside = true;
			}
			else
			{
				hasOnePointOutside = true;
			}
			
			if (hasOnePointInside && hasOnePointOutside)
			{
				isSurfaceCell = true;
			}
		}
		*/
		int xIndex = getGridXIndex(cell);
		int yIndex = getGridYIndex(cell);
		int startX = (xIndex>0) ? xIndex-1 : 0;
		int startY = (yIndex>0) ? yIndex-1 : 0;
		int endX = (xIndex<(_gridResX-1)) ? xIndex+1 : (_gridResX-1);
		int endY = (yIndex<(_gridResY-1)) ? yIndex+1 : (_gridResY-1);
		bool hasOneCellOutside = false;
		bool hasOneCellInside = false;
		for (int ix=startX; !isSurfaceCell && (ix<=endX); ++ix)
		{
			for (int iy=startY; !isSurfaceCell && (iy<=endY); ++iy)
			{
				Vec3d middlePoint;
				middlePoint.x = _gridMin.x + static_cast<double>(ix)*_cellSize + 0.5*_cellSize;
				middlePoint.y = _gridMin.y + static_cast<double>(iy)*_cellSize + 0.5*_cellSize;
				
				if (v->pointIsInside(middlePoint))
				{
					hasOneCellInside = true;
				}
				else
				{
					hasOneCellOutside = true;
				}
				
				if (hasOneCellInside && hasOneCellOutside)
				{
					isSurfaceCell = true;
				}
			}
		}
		
		// Process only surface cells
		if (!isSurfaceCell)
		{
			continue;
		}
		
		// Try random points until poisson criterion or stopping criterion is met
		Vec3d p;
		Vec3d normal;
		Vec3d projP;
		Vec3d randPoint;
		Vec3d temp;
		bool pointFound = false;
		p.z = 0.0;
		for (int t=0; t<30; ++t)
		{
			// Create a random point and project it onto the surface
			randPoint.x = (static_cast<double>(rand())/randMax) * _cellSize + _gridMin.x + getGridXIndex(cell)*_cellSize;
			randPoint.y = (static_cast<double>(rand())/randMax) * _cellSize + _gridMin.y + getGridYIndex(cell)*_cellSize;
			v->projectPointOntoSurface(randPoint, projP, normal);
			
			// Verify poisson disk criterion
			if (satisfyPoissonDiskCriterion(projP))
			{
				// Add the point!
				p = projP;
				p.z = v->getmin().z;
				points.push_back(p);
				int index = points.size()-1;
				_grid[getGridIndex(p)] = index;
				surfaceSamples.push_back(index);
				pointFound = true;
				break;
			}
		}
		
		// If a point has been found, expand from there!!!! Mouahahah, conquer the whole surface!!!
		while (pointFound)
		{
			pointFound = false;
			
			// Generate a random tangential direction
			Vec3d dir;
			dir.x = normal.y;
			dir.y = normal.x;
			dir.z = 0.0;
			if ((static_cast<double>(rand())/randMax) >= 0.5)
			{
				dir.x *= -1.0;
			}
			else
			{
				dir.y *= -1.0;
			}
			
			// Generate a point in the direction of this direction
			dir *= r*epsilon;
			Vec3d q = p;
			q += dir;
			
			// Project onto surface and check to see if Poisson disk criterion is met
			v->projectPointOntoSurface(q, projP, normal);
			if (satisfyPoissonDiskCriterion(projP))
			{
				p = projP;
				p.z = v->getmin().z;
				points.push_back(p);
				int index = points.size()-1;
				_grid[getGridIndex(p)] = index;
				surfaceSamples.push_back(index);
				pointFound = true;
			}
		}
	}
}

void PoissonFillaar::relaxation(const Volume* v, const std::vector<int>& pointsToRelax, int nbSweeps, bool alwaysProjectToSurface)
{
	// Declarations
	double r = _resolution;
	double twoR = 2.0*r;
    double twoR2 = twoR*twoR;
	double randMax = static_cast<double>(RAND_MAX);
	double PI = 3.14159265;
	double twoPI = 2.0*PI;
	
	std::deque<Vec3d>& points = getPoints();
	std::vector<int> neighbors;

	// Relax those tense points ;)
	for (int k=0; k<nbSweeps; ++k)
	{
		std::vector<int>::const_iterator it = pointsToRelax.begin();
		std::vector<int>::const_iterator itEnd = pointsToRelax.end();
		for (; it!=itEnd; ++it)
		{
			Vec3d *p = &points[*it];
			
			neighbors.clear();
			
			// Find neighbors within 2.0*r of p and compute
			// the smallest distance between p and its neighbors
            double dmin2 = 100000000000000000.0;
			int startX = getGridXIndex(p->x - twoR); if (startX<0) startX = 0;
			int startY = getGridYIndex(p->y - twoR); if (startY<0) startY = 0;
			int endX = getGridXIndex(p->x + twoR); if (endX>=_gridResX) endX = _gridResX-1;
			int endY = getGridYIndex(p->y + twoR); if (endY>=_gridResY) endY = _gridResY-1;
			bool hasNeighbors = false;
			for (int ix=startX; ix<=endX; ++ix)
			{
				for (int iy=startY; iy<=endY; ++iy)
				{
					int id = _grid[getGridIndex(ix, iy)];
					if ((id < 0) || (id==(*it)))
						continue;
					
					Vec3d delta = points[id];
					delta -= *p;
                    double dist2 = delta.length2();
                    if (dist2 <= twoR2)
					{
						hasNeighbors = true;
						neighbors.push_back(id);
                        if (dist2<dmin2)
						{
                            dmin2 = dist2;
						}
					}
				}
			}
			
			if (!hasNeighbors)
			{
				continue;
			}
			
			// Try to find a better point (a point further from its nearest neighbors)
			Vec3d pnew = *p;
			double tmax = 100;
			for (int t=0; t<tmax; ++t)
			{
				double tau = static_cast<double>(tmax-t)/static_cast<double>(tmax);
				
				// Generate a random point around p
				double randTheta = (static_cast<double>(rand())/randMax) * twoPI;
				Vec3d pcand(cos(randTheta),sin(randTheta),0.0);
				pcand *= r*tau;
				pcand += *p;
				
				// Project point onto surface if necessary
				if (alwaysProjectToSurface || !v->pointIsInside(pcand))
				{
					Vec3d proj;
					Vec3d normal;
					v->projectPointOntoSurface(pcand, proj, normal);
					pcand = proj;
				}
				
				// Get the smallest distance between pcand and its neighbors
				int nbNeighbors = neighbors.size();
                double dminCand2 = 100000000000000000.0;
				for (int i=0; i<nbNeighbors; ++i)
				{
					Vec3d delta = points[neighbors[i]];
					delta -= pcand;
                    double dist2 = delta.length();
                    if (dist2<dminCand2)
					{
                        dminCand2 = dist2;
					}
				}
				
                if (dminCand2 > dmin2)
				{
					pnew = pcand;
                    dmin2 = dminCand2;
				}
			}
			
			// Update p!
			*p = pnew;
		}
	}
}

void PoissonFillaar::volumeSampling(const Volume* v, const std::vector<int>& initialSeeds, std::vector< int >& volumeSamples)
{
	// Declarations
	const double r = _resolution;
	const double r2 = r*r;
	const int k = 30;
	const double randMax = static_cast<double>(RAND_MAX);
	const double twoPI = 2.0*3.1415926535897932384626433832795;
	
	Vec3d temp;
	std::deque<Vec3d>& points = getPoints();
	std::vector<int> activePoints(initialSeeds);
	
	// Add particles inside the volume using poisson disk (see Schechter2012 & Bridson2007)
	int i=0;
	while (activePoints.size() > 0)
	{
		++i;
		// Process random point from the active list
		int index = static_cast<int>(floor((static_cast<double>(rand())/randMax)*(activePoints.size()-1)+0.5));
		Vec3d *p = &points[activePoints[index]];
		
		bool pointFound = false;
		for (int t=0; !pointFound && (t<k); ++t)
		{
			// Generate a random point in the spherical annulus between radius r and 2r of p
			double theta = (static_cast<double>(rand())/randMax)*twoPI;
			double pr = (static_cast<double>(rand())/randMax)*r + r;
			
			Vec3d pos = *p;
			pos.x += pr*cos(theta);
			pos.y += pr*sin(theta);
			
			// If point is outside volume, project it inside
			if (!v->pointIsInside(pos))
			{
				Vec3d proj;
				Vec3d normal;
				v->projectPointOntoSurface(pos, proj, normal);
				
				pos = proj;
			}
			
			// Add point if poisson criterion is satisfied
			if (satisfyPoissonDiskCriterion(pos))
			{
				pointFound = true;
				
				pos.z = v->getmin().z;
				points.push_back(pos);
				int newIndex = points.size()-1;
				_grid[getGridIndex(pos)] = newIndex;
				activePoints.push_back(newIndex);
				volumeSamples.push_back(newIndex);
			}
		}
		
		// Remove point from active list if no new neighbor has been found
		if (!pointFound)
		{
			activePoints.erase(activePoints.begin()+index);
		}
	}
}




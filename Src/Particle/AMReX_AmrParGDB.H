#ifndef _AmrParGDB_H
#define _AmrParGDB_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_AmrCore.H>

namespace amrex {

class AmrParGDB
{
public:

    AmrParGDB (AmrCore* amr)
	: m_amrcore(amr),
	  m_dmap(amr->maxLevel()+1),
	  m_ba(amr->maxLevel()+1)
	{ }

    ~AmrParGDB () {;}
    
    const Geometry& Geom (int level) const;
    const DistributionMapping& ParticleDistributionMap (int level) const;
    const DistributionMapping&         DistributionMap (int level) const;
    const BoxArray&            ParticleBoxArray (int level) const;
    const BoxArray&                    boxArray (int level) const;

    void SetParticleBoxArray (int level, const BoxArray& new_ba);
    void SetParticleDistributionMap (int level,	const DistributionMapping& new_dm);

    void ClearParticleBoxArray (int level);
    void ClearParticleDistributionMap (int level);

    bool LevelDefined (int level) const;
    int finestLevel () const;
    int maxLevel () const;

    IntVect refRatio (int level) const;
    int MaxRefRatio (int level) const;

    bool OnSameGrids (int level, const MultiFab& mf) const;

protected:
    AmrCore* m_amrcore;
    Array<DistributionMapping> m_dmap;
    Array<BoxArray>            m_ba;
};

inline 
const Geometry& 
AmrParGDB::Geom (int level) const
{
    return m_amrcore->Geom(level);
}

inline 
const DistributionMapping& 
AmrParGDB::ParticleDistributionMap (int level) const
{
    if (m_dmap[level].empty()) {
	return m_amrcore->DistributionMap(level);
    } else {
	return m_dmap[level];
    }
}

inline 
const DistributionMapping& 
AmrParGDB::DistributionMap (int level) const
{
    return m_amrcore->DistributionMap(level);
}

inline 
const BoxArray&
AmrParGDB::ParticleBoxArray (int level) const
{
    if (m_ba[level].empty()) {
	return m_amrcore->boxArray(level);
    } else {
	return m_ba[level];
    }
}

inline
const BoxArray&
AmrParGDB::boxArray (int level) const
{
    return m_amrcore->boxArray(level);
}

inline
void AmrParGDB::SetParticleBoxArray (int level, const BoxArray& new_ba)
{
    m_ba[level] = new_ba;
}

inline
void AmrParGDB::SetParticleDistributionMap (int level, const DistributionMapping& new_dmap)
{
    m_dmap[level] = new_dmap;
}

inline
void AmrParGDB::ClearParticleBoxArray (int level)
{
    m_ba[level] = BoxArray();
}

inline
void AmrParGDB::ClearParticleDistributionMap (int level)
{
    m_dmap[level] = DistributionMapping();
}

inline 
bool 
AmrParGDB::LevelDefined (int level) const
{
    return m_amrcore->LevelDefined(level);
}

inline 
int 
AmrParGDB::finestLevel () const
{
    return m_amrcore->finestLevel();
}

inline 
int 
AmrParGDB::maxLevel () const
{
    return m_amrcore->maxLevel();
}
 
inline 
IntVect 
AmrParGDB::refRatio (int level) const
{
    return m_amrcore->refRatio(level);
}

inline 
int 
AmrParGDB::MaxRefRatio (int level) const
{
    return m_amrcore->MaxRefRatio(level);
}

inline
bool
AmrParGDB::OnSameGrids (int level, const MultiFab& mf) const
{
    return (mf.DistributionMap() == ParticleDistributionMap(level) &&
	    mf.boxArray().CellEqual(ParticleBoxArray(level)));
}

}

#endif

/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(LUDING,luding,5)
#else
#ifndef NORMAL_MODEL_LUDING_H_
#define NORMAL_MODEL_LUDING_H_
#include "contact_models.h"
#include <math.h>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "global_properties.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<LUDING> : public NormalModel<HOOKE>
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_SURFACES_INTERSECT | CM_SURFACES_CLOSE;

    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *c) :
        NormalModel<HOOKE>(lmp, hsetup,c),
        k2MAX(NULL),
        k1(NULL),
        kC(NULL),
        phiF(NULL)
    {
      history_offset = hsetup->add_history_value("deltaMax", "0");
      
    }

    inline void registerSettings(Settings & settings){
      NormalModel<HOOKE>::registerSettings(settings);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    inline void connectToProperties(PropertyRegistry & registry) {
      NormalModel<HOOKE>::connectToProperties(registry);

      //registry.registerProperty("kn2Max", &MODEL_PARAMS::createCoeffMaxElasticStiffness);
      //registry.registerProperty("kn2kc", &MODEL_PARAMS::createCoeffAdhesionStiffness);
      //registry.registerProperty("phiF", &MODEL_PARAMS::createCoeffPlasticityDepth);

      //registry.connect("kn2kcMax", kn2k2Max,"model hooke/hysteresis");
      //registry.connect("kn2kc", kn2kc,"model hooke/hysteresis");
      //registry.connect("phiF", phiF,"model hooke/hysteresis");
      //---- Mohammad

      registry.registerProperty("k2MAX", &MODEL_PARAMS::createUnloadingStiffness);
      registry.registerProperty("k1", &MODEL_PARAMS::createLoadingStiffness);
      registry.registerProperty("kC", &MODEL_PARAMS::createCohesionStrength);
      registry.registerProperty("phiF", &MODEL_PARAMS::createCoeffPlasticityDepth);

      registry.connect("k2MAX", k2MAX,"model luding");
      registry.connect("k1", k1,"model luding");
      registry.connect("kC", kC,"model luding");
      registry.connect("phiF", phiF,"model luding");

      registry.registerProperty("gammaN", &MODEL_PARAMS::createGamman);
      registry.connect("gammaN", gammaN,"model luding");

      //--------------------
      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model luding");
	//error->cg(FLERR,"model hooke/hysteresis");
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.;
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      // use these values from HOOKE implementation
      bool & viscous = NormalModel<HOOKE>::viscous;
      double ** & Yeff = NormalModel<HOOKE>::Yeff;
      double & charVel = NormalModel<HOOKE>::charVel;
      bool & tangential_damping = NormalModel<HOOKE>::tangential_damping;
      Force * & force = NormalModel<HOOKE>::force;

      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double deltan = sidata.deltan;
      const double radi = sidata.radi;
      const double radj = sidata.radj;
      double reff=sidata.is_wall ? radi : (radi*radj/(radi+radj));
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical) {
        if(sidata.is_wall)
          reff = MathExtraLiggghtsNonspherical::get_effective_radius_wall(sidata, atom->roundness[sidata.i], error);
        else
          reff = MathExtraLiggghtsNonspherical::get_effective_radius(sidata, atom->roundness[sidata.i], atom->roundness[sidata.j], error);
      }
#endif
      double meff=sidata.meff;
      double coeffRestLogChosen;

      if (viscous)  {
        double ** & coeffMu = NormalModel<HOOKE>::coeffMu;
        double ** & coeffRestMax = NormalModel<HOOKE>::coeffRestMax;
        double ** & coeffStc = NormalModel<HOOKE>::coeffStc;
        // Stokes Number from MW Schmeeckle (2001)
        const double stokes=sidata.meff*sidata.vn/(6.0*M_PI*coeffMu[itype][jtype]*reff*reff);
        // Empirical from Legendre (2006)
        coeffRestLogChosen=log(coeffRestMax[itype][jtype])+coeffStc[itype][jtype]/stokes;
      } else {
        double ** & coeffRestLog = NormalModel<HOOKE>::coeffRestLog;
        coeffRestLogChosen=coeffRestLog[itype][jtype];
      }

      const double sqrtval = sqrt(reff);
      
     // double kn = 16./15.*sqrtval*(Yeff[itype][jtype])*pow(15.*meff*charVel*charVel/(16.*sqrtval*Yeff[itype][jtype]),0.2);
     // double kt = kn;
      double kn = k1[itype][jtype];
      double kt = kn;
      const double gamman = gammaN[itype][jtype];
      const double gammat = tangential_damping ? gamman : 0.0;

      // convert Kn and Kt from pressure units to force/distance^2
      //kn /= force->nktv2p;
      //kt /= force->nktv2p;

      // coefficients
      
      //const double k2Max = kn * kn2k2Max[itype][jtype]; 
      //const double kc = kn * kn2kc[itype][jtype]; 
      //--- Mohammad
        const double k2Max = k2MAX[itype][jtype]; 
        const double kc = kC[itype][jtype];
      //-----------
      // get the history value -- maximal overlap
      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
      double * const history = &sidata.contact_history[history_offset];
      if (deltan > history[0]) {
          history[0] = deltan;
      }


      // k2 dependent on the maximum overlap
      // this accounts for an increasing stiffness with deformation
      //const double deltaMaxLim =(k2Max/(k2Max-kn))*phiF[itype][jtype]*2*reff;
      const double deltaPmax = (k2Max/(k2Max-kn))*phiF[itype][jtype]*2.0*reff;

      double k2, fHys;
      const bool update_history = sidata.computeflag && sidata.shearupdate;
      double deltaMax = std::min(deltaPmax,std::max(history[0],deltan)); 

      // k2(deltamax)
      if (deltaMax >= deltaPmax) 
	k2 = k2Max ;
      else 
        k2 = kn+(k2Max-kn)*deltaMax/deltaPmax;

      const double deltaZero = (1.-kn/k2)*deltaMax ;
      // Large overlap
      if (deltan > deltaPmax) 
      {
            fHys = k2*(deltan-deltaZero);//k2*(deltan-delta0);
      } 
      // Loading part
      else if (k2*(deltan-deltaZero) >= kn*deltan)
      { 
            fHys = kn*deltan;
      }
      // un/re-loading part
      else if (k2*(deltan-deltaZero) >= -kc*deltan and k2*(deltan-deltaZero) < kn*deltan)
      { 
            fHys = k2*(deltan-deltaZero);
      }
      // Cohesion part
      else 
      { 
            fHys = -kc*deltan;
            const double newDeltaMax = ((k2+kc)/(k2-kn))*deltan;
           // if (update_history)
                  history[0] = newDeltaMax;
      }

      const double Fn_damping = -gamman*sidata.vn;
      const double Fn = fHys + Fn_damping;

      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      sidata.gamman = gamman;
      sidata.gammat = gammat;

      #ifdef NONSPHERICAL_ACTIVE_FLAG
          double Fn_i[3] = { Fn * sidata.en[0], Fn * sidata.en[1], Fn * sidata.en[2]};
          double torque_i[3] = {0.0, 0.0, 0.0}; //initialized here with zeros to avoid compiler warnings
          if(sidata.is_non_spherical) {
            double xci[3];
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Fn_i, torque_i);
          }
      #endif
      // apply normal force
      if(sidata.is_wall) {
        const double Fn_ = Fn * sidata.area_ratio;
        i_forces.delta_F[0] += Fn_ * sidata.en[0];
        i_forces.delta_F[1] += Fn_ * sidata.en[1];
        i_forces.delta_F[2] += Fn_ * sidata.en[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
                if(sidata.is_non_spherical) {
                  //for non-spherical particles normal force can produce torque!
                  i_forces.delta_torque[0] += torque_i[0];
                  i_forces.delta_torque[1] += torque_i[1];
                  i_forces.delta_torque[2] += torque_i[2];
                }
        #endif
      } else {
        i_forces.delta_F[0] += sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] += sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] += sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] += -i_forces.delta_F[0];
        j_forces.delta_F[1] += -i_forces.delta_F[1];
        j_forces.delta_F[2] += -i_forces.delta_F[2];
        #ifdef NONSPHERICAL_ACTIVE_FLAG
                if(sidata.is_non_spherical) {
                  //for non-spherical particles normal force can produce torque!
                  double xcj[3], torque_j[3];
                  double Fn_j[3] = { -Fn_i[0], -Fn_i[1], -Fn_i[2]};
                  vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
                  vectorCross3D(xcj, Fn_j, torque_j);

                  i_forces.delta_torque[0] += torque_i[0];
                  i_forces.delta_torque[1] += torque_i[1];
                  i_forces.delta_torque[2] += torque_i[2];

                  j_forces.delta_torque[0] += torque_j[0];
                  j_forces.delta_torque[1] += torque_j[1];
                  j_forces.delta_torque[2] += torque_j[2];
                }
        #endif
      }
    }

    inline void surfacesClose(SurfacesCloseData & scdata, ForceData&, ForceData&)
    {
      if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_NORMAL_MODEL;
      double * const history = &scdata.contact_history[history_offset];
      history[0] = 0.0;
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double **k2MAX;
    double **k1;
    double **kC;
    double **phiF;
    double **gammaN;
    int history_offset;
  };
}
}
#endif // NORMAL_MODEL_LUDING_H_
#endif

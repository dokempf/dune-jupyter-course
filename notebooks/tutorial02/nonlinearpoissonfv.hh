#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/type.hh>

#include<dune/geometry/referenceelements.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

/** a local operator for solving the nonlinear Poisson equation with cell-centered finite volume method
 *
 * \f{align*}{
 *   -\Delta u(x) + q(u(x)) &=& f(x) x\in\Omega,  \\
 *                     u(x) &=& g(x) x\in\partial\Omega_D \\
 *  -\nabla u(x) \cdot n(x) &=& j(x) x\in\partial\Omega_N \\
 * \f}
 *
 */
template<typename Param>
class NonlinearPoissonFV :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::FullSkeletonPattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
  Param& param;        // parameter functions

public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  enum { doLambdaVolume = true };
  enum { doLambdaBoundary = true };
  enum { doAlphaVolume = true };
  enum { doAlphaSkeleton  = true };
  enum { doAlphaBoundary  = true };

  //! constructor stores a copy of the parameter object
  NonlinearPoissonFV (Param& param_)
    : param(param_)
  {}

  //! residual contribution of volume source term
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv,
                      R& r) const
  {
    // center of reference element
    auto cellgeo = eg.geometry();
    auto cellcenterlocal =
      referenceElement(cellgeo).position(0,0);

    // accumulate residual
    auto f = param.f(eg.entity(),cellcenterlocal);
    r.accumulate(lfsv,0,-f*cellgeo.volume());
  }

  //! residual contribution of boundary integral (Neumann boundary condition)
  template<typename IG, typename LFSV, typename R>
  void lambda_boundary (const IG& ig, const LFSV& lfsv_i,
                        R& r_i) const
  {
    // face volume for integration
    auto facegeo = ig.geometry();
    auto facecenterlocal =
      referenceElement(facegeo).position(0,0);

    // evaluate boundary condition and quit on Dirichlet
    bool isdirichlet =
      param.b(ig.intersection(),facecenterlocal);

    if (isdirichlet)
      {
        // inside cell center
        auto insidecenterglobal=ig.inside().geometry().center();

        // face center in global coordinates
        auto facecenterglobal = facegeo.center();

        // compute distance of these two points
        insidecenterglobal -= facecenterglobal;
        auto distance = insidecenterglobal.two_norm();

        // face center in local coordinates of the element
        auto facecenterinelement=ig.geometryInInside().center();

        // evaluate Dirichlet condition
        auto g = param.g(ig.inside(),facecenterinelement);

        // face volume for integration
        auto face_volume = facegeo.volume();

        // contribution to residual
        r_i.accumulate(lfsv_i,0,-g/distance*face_volume);
      }
    else
      {
        // contribution to residual from Neumann boundary
        auto j = param.j(ig.intersection(),facecenterlocal);
        r_i.accumulate(lfsv_i,0,j*facegeo.volume());
      }
  }

  //! residual contribution of volume integral (reaction term)
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    // get cell value
    auto u = x(lfsu,0);

    // evaluate reaction term
    auto q = param.q(u);

    // and accumulate
    r.accumulate(lfsv,0,q*eg.geometry().volume());
  }

  //! jacobian contribution of volume term (reaction term)
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x,
                        const LFSV& lfsv, M& mat) const
  {
    // evaluate derivative reaction term
    auto u = x(lfsu,0);
    auto qprime = param.qprime(u);

    // and accumulate
    mat.accumulate(lfsv,0,lfsu,0,qprime*eg.geometry().volume());
  }

  //! apply local jacobian of the volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z,
                              const LFSV& lfsv, R& r) const
  {
    // evaluate derivative reaction term
    auto u = x(lfsu,0);
    auto qprime = param.qprime(u);

    // and accumulate
    r.accumulate(lfsv,0,qprime*z(lfsu,0)*eg.geometry().volume());
  }

  //! residual contribution from skeleton terms
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
         const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
         const LFSU& lfsu_o, const X& x_o, const LFSV& lfsv_o,
         R& r_i, R& r_o) const
  {
    // inside and outside cells
    auto cell_inside = ig.inside();
    auto cell_outside = ig.outside();

    // inside and outside cell geometries
    auto insidegeo = cell_inside.geometry();
    auto outsidegeo = cell_outside.geometry();

    // cell centers in global coordinates
    auto inside_global = insidegeo.center();
    auto outside_global = outsidegeo.center();

    // distance between the two cell centers
    inside_global -= outside_global;
    auto distance = inside_global.two_norm();

    // face volume for integration
    auto facegeo = ig.geometry();
    auto face_volume = facegeo.volume();

    // contribution to residual on inside and outside elements
    auto dudn = (x_o(lfsu_o,0)-x_i(lfsu_i,0))/distance;
    r_i.accumulate(lfsv_i,0,-dudn*face_volume);
    r_o.accumulate(lfsv_o,0, dudn*face_volume);
  }

  //! Jacobian contribution from skeleton terms
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_skeleton (const IG& ig,
         const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
         const LFSU& lfsu_o, const X& x_o, const LFSV& lfsv_o,
         M& mat_ii, M& mat_io,
         M& mat_oi, M& mat_oo) const
  {
    // inside and outside cells
    auto cell_inside = ig.inside();
    auto cell_outside = ig.outside();

    // inside and outside cell geometries
    auto insidegeo = cell_inside.geometry();
    auto outsidegeo = cell_outside.geometry();

    // cell centers in global coordinates
    auto inside_global = insidegeo.center();
    auto outside_global = outsidegeo.center();

    // distance between the two cell centers
    inside_global -= outside_global;
    auto distance = inside_global.two_norm();

    // face volume for integration
    auto facegeo = ig.geometry();
    auto face_volume = facegeo.volume();

    // contribution to jacobian entries
    mat_ii.accumulate(lfsv_i,0,lfsv_i,0, face_volume/distance);
    mat_io.accumulate(lfsv_i,0,lfsv_o,0,-face_volume/distance);
    mat_oi.accumulate(lfsv_o,0,lfsv_i,0,-face_volume/distance);
    mat_oo.accumulate(lfsv_o,0,lfsv_o,0, face_volume/distance);
  }

  //! apply local jacobian of the skeleton term
  template<typename IG, typename LFSU, typename X, typename LFSV,
           typename Y>
  void jacobian_apply_skeleton
  ( const IG& ig,
    const LFSU& lfsu_i, const X& x_i, const X& z_i, const LFSV& lfsv_i,
    const LFSU& lfsu_o, const X& x_o, const X& z_o, const LFSV& lfsv_o,
    Y& y_i, Y& y_o) const
  {
    // reuse alpha_boundary because it is linear
    alpha_skeleton(ig,lfsu_i,z_i,lfsv_i,lfsu_o,z_o,lfsv_o,y_i,y_o);
  }

  //! residual contribution of boundary integral (Dirichlet condition)
  // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                       const LFSU& lfsu_i, const X& x_i,
                       const LFSV& lfsv_i, R& r_i) const
  {
    // check for Dirichlet boundary condition
    auto facegeo = ig.geometry();
    auto facecenterlocal =
      referenceElement(facegeo).position(0,0);
    bool isdirichlet = param.b(ig.intersection(),facecenterlocal);
    if (!isdirichlet) return;

    // inside cell center
    auto insidecenterglobal = ig.inside().geometry().center();

    // face center in global coordinates
    auto facecenterglobal = facegeo.center();

    // compute distance of these two points
    insidecenterglobal -= facecenterglobal;
    auto distance = insidecenterglobal.two_norm();

    // face volume for integration
    auto face_volume = facegeo.volume();

    // contribution to residual
    r_i.accumulate(lfsv_i,0,x_i(lfsu_i,0)/distance*face_volume);
  }

  //! Jacobian contribution from boundary integral (Dirichlet condition)
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_boundary (const IG& ig,
                          const LFSU& lfsu_i, const X& x_i,
                          const LFSV& lfsv_i, M& mat_ii) const
  {
    // face volume for integration
    auto facegeo = ig.geometry();
    auto facecenterlocal = referenceElement(facegeo).position(0,0);

    // evaluate boundary condition and quit on NOT Dirichlet
    bool isdirichlet = param.b(ig.intersection(),facecenterlocal);
    if (!isdirichlet) return;

    // inside cell center
    auto insidecenterglobal = ig.inside().geometry().center();

    // face center in global coordinates
    auto facecenterglobal = facegeo.center();

    // compute distance of these two points
    insidecenterglobal -= facecenterglobal;
    auto distance = insidecenterglobal.two_norm();

    // face volume for integration
    auto face_volume = facegeo.volume();

    // contribution to matrix
    mat_ii.accumulate(lfsv_i,0,lfsv_i,0,face_volume/distance);
  }

  //! apply local jacobian of the boundaryterm
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename Y>
  void jacobian_apply_boundary
  ( const IG& ig,
    const LFSU& lfsu_i, const X& x_i, const X& z_i,
    const LFSV& lfsv_i, Y& y_i) const
  {
    // reuse alpha_boundary because it is linear
    alpha_boundary(ig,lfsu_i,z_i,lfsv_i,y_i);
  }
};

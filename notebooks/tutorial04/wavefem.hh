#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/type.hh>

#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

/** a local operator for the spatial part of the wave equation as a first order system in time
 *
 * \f{align*}{
 *   -c^2\Delta u_1(x) &=& 0  x\in\Omega, \\
 *             -u_0(x) &=& 0  x\in\Omega, \\
 *              u_1(x) &=& 0  x\in\partial\Omega \\
 *  no boundary condition for u_1
 * \f}
 *
 * We assume that both components use the same finite element space.
 * \tparam FEM is a finite element map of the space used for both components
 */
template<typename FEM>
class WaveFEM :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  // types
  using LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  using RF = typename LocalBasis::Traits::RangeFieldType;

  // private data members
  Dune::PDELab::LocalBasisCache<LocalBasis> cache; // a cache for local basis evaluations
  RF c; // the wave speed

public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  //! constructor stores a copy of the parameter object
  WaveFEM (RF c_) : c(c_) {}

  //! volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // select the two components (but assume Galerkin scheme U=V)
    using namespace Dune::TypeTree::Indices;
    auto lfsu0 = lfsu.child(_0);
    auto lfsu1 = lfsu.child(_1);

    // types & dimension
    const int dim = EG::Entity::dimension;

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = 2*lfsu0.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
        // evaluate basis functions; Assume basis is the same for both components
        auto& phihat = cache.evaluateFunction(ip.position(),lfsu0.finiteElement().localBasis());

        // evaluate u1
        RF u1=0.0;
        for (std::size_t i=0; i<lfsu0.size(); i++) u1 += x(lfsu1,i)*phihat[i];

        // evaluate gradient of shape functions
        auto& gradphihat = cache.evaluateJacobian(ip.position(),lfsu0.finiteElement().localBasis());

        // transform gradients of shape functions to real element
        const auto S = geo.jacobianInverseTransposed(ip.position());
        auto gradphi = makeJacobianContainer(lfsu0);
        for (std::size_t i=0; i<lfsu0.size(); i++)
          S.mv(gradphihat[i][0],gradphi[i][0]);

        // compute gradient of u0
        Dune::FieldVector<RF,dim> gradu0(0.0);
        for (std::size_t i=0; i<lfsu0.size(); i++)
          gradu0.axpy(x(lfsu0,i),gradphi[i][0]);

        // integrate both equations
        RF factor = ip.weight() * geo.integrationElement(ip.position());
        for (std::size_t i=0; i<lfsu0.size(); i++) {
          r.accumulate(lfsu0,i,c*c*(gradu0*gradphi[i][0])*factor);
          r.accumulate(lfsu1,i,-u1*phihat[i]*factor);
        }
      }
  }

    //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                        M& mat) const
  {
    // select the two components (assume Galerkin scheme U=V)
    using namespace Dune::TypeTree::Indices;
    auto lfsu0 = lfsu.child(_0);
    auto lfsu1 = lfsu.child(_1);

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = 2*lfsu0.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
        // evaluate basis functions; Assume basis is the same for both components
        auto& phihat = cache.evaluateFunction(ip.position(),lfsu0.finiteElement().localBasis());

        // evaluate gradient of shape functions
        auto& gradphihat = cache.evaluateJacobian(ip.position(),lfsu0.finiteElement().localBasis());

        // transform gradients of shape functions to real element
        const auto S = geo.jacobianInverseTransposed(ip.position());
        auto gradphi = makeJacobianContainer(lfsu0);
        for (std::size_t i=0; i<lfsu0.size(); i++)
          S.mv(gradphihat[i][0],gradphi[i][0]);

        // integrate both equations
        RF factor = ip.weight() * geo.integrationElement(ip.position());
        for (std::size_t j=0; j<lfsu0.size(); j++)
          for (std::size_t i=0; i<lfsu0.size(); i++) {
            mat.accumulate(lfsu0,i,lfsu0,j,c*c*(gradphi[j][0]*gradphi[i][0])*factor);
            mat.accumulate(lfsu1,i,lfsu1,j,-phihat[j]*phihat[i]*factor);
          }
      }
  }

  //! apply local jacobian of the volume term -> nonlinear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,x,lfsv,r);
  }
};


/** a local operator for the temporal operator in the wave equation assuming identical components
 *
 * \f{align*}{
 \int_\Omega uv dx
 * \f}
 * \tparam FEM      Type of a finite element map
 */
template<typename FEM>
class WaveL2
  : public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  // types
  using LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  using RF = typename LocalBasis::Traits::RangeFieldType;

  // private data members
  Dune::PDELab::LocalBasisCache<LocalBasis> cache; // a cache for local basis evaluations

public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // select the two components (assume Galerkin scheme U=V)
    using namespace Dune::TypeTree::Indices;
    auto lfsu0 = lfsu.child(_0);
    auto lfsu1 = lfsu.child(_1);

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = 2*lfsu0.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
        // evaluate basis functions at first child
        auto& phihat = cache.evaluateFunction(ip.position(),lfsu0.finiteElement().localBasis());

        // evaluate u0
        RF u0=0.0, u1=0.0;
        for (std::size_t i=0; i<lfsu0.size(); i++) {
          u0 += x(lfsu0,i)*phihat[i];
          u1 += x(lfsu1,i)*phihat[i];
        }

        // integration factor
        RF factor = ip.weight() * geo.integrationElement(ip.position());

        // integrate u*phi_i
        for (std::size_t i=0; i<lfsu0.size(); i++) {
          r.accumulate(lfsu0,i,u1*phihat[i]*factor);
          r.accumulate(lfsu1,i,u0*phihat[i]*factor);
        }
      }
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                        M& mat) const
  {
    // get first child, assuming PowerGridFunctionSpace
    using namespace Dune::TypeTree::Indices;
    auto lfsu0 = lfsu.child(_0);
    auto lfsu1 = lfsu.child(_1);

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = 2*lfsu0.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
        // evaluate basis functions at first child
        auto& phihat = cache.evaluateFunction(ip.position(),lfsu0.finiteElement().localBasis());

        // integration factor
        RF factor = ip.weight() * geo.integrationElement(ip.position());

        // loop over all components
        for (std::size_t j=0; j<lfsu0.size(); j++)
          for (std::size_t i=0; i<lfsu0.size(); i++) {
            mat.accumulate(lfsu0,i,lfsu1,j,phihat[j]*phihat[i]*factor);
            mat.accumulate(lfsu1,i,lfsu0,j,phihat[j]*phihat[i]*factor);
          }
      }
  }

  //! apply local jacobian of the volume term -> nonlinear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,x,lfsv,r);
  }
};

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/type.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>


/** a local operator for residual-based error estimation for the problem
 *
 * \f{align*}{
 *   -\Delta u(x) + q(u(x)) &=& f(x) x\in\Omega,  \\
 *                     u(x) &=& g(x) x\in\partial\Omega_D \\
 *  -\nabla u(x) \cdot n(x) &=& j(x) x\in\partial\Omega_N \\
 * \f}
 *
 * A call to residual() of a grid operator space will assemble
 * the quantity \f$\eta_T^2\f$ for each cell. Note that the squares
 * of the cell indicator \f$\eta_T\f$ is stored. To compute the global
 * error estimate sum up all values and take the square root.
 *
 * Assumptions and limitations:
 * - Assumes that LFSU is \f$P_k\f$/\f$Q_k\f$ finite element space
 *   and LFSV is a \f$P_0\f$ finite element space (one value per cell).
 * - However, the second order derivatives are ignored!
 *
 */
template<typename Param, typename FEM>
class NonlinearPoissonFEMEstimator
  : public Dune::PDELab::LocalOperatorDefaultFlags
{
  // a cache for local basis evaluations
  typedef typename FEM::Traits::FiniteElementType
    ::Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
  Param& param; // parameter functions
  int incrementorder; // additional increment for integration order

  // a function to compute the diameter of an entity
  template<class GEO>
  typename GEO::ctype diameter (const GEO& geo) const
  {
    typedef typename GEO::ctype DF;
    DF hmax = -1.0E00;
    for (int i=0; i<geo.corners(); i++)
      {
        auto xi = geo.corner(i);
        for (int j=i+1; j<geo.corners(); j++)
          {
            auto xj = geo.corner(j);
            xj -= xi;
            hmax = std::max(hmax,xj.two_norm());
          }
      }
    return hmax;
  }

public:
  // pattern assembly flags
  enum { doPatternVolume = false };
  enum { doPatternSkeleton = false };

  // residual assembly flags
  enum { doAlphaVolume  = true };
  enum { doAlphaSkeleton  = true };
  enum { doAlphaBoundary  = true };

  //! constructor: pass parameter object
  NonlinearPoissonFEMEstimator (Param& param_, int incrementorder_=0)
    : param(param_), incrementorder(incrementorder_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
    typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu,
      const X& x, const LFSV& lfsv, R& r) const
  {
    // types & dimension
    typedef decltype(makeZeroBasisFieldValue(lfsu)) RF;

    // select quadrature rule
    auto geo = eg.geometry();
    const int order =
      incrementorder+2*lfsu.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    RF sum(0.0);
    for (const auto& ip : rule)
      {
        // evaluate basis functions
        auto& phihat = cache.evaluateFunction(
            ip.position(),lfsu.finiteElement().localBasis());

        // evaluate u
        RF u=0.0;
        for (size_t i=0; i<lfsu.size(); i++) u += x(lfsu,i)*phihat[i];

        // evaluate reaction term
        auto q = param.q(u);

        // evaluate right hand side parameter function
        auto f = param.f(eg.entity(),ip.position());

        // integrate f^2
        RF factor =
          ip.weight() * geo.integrationElement(ip.position());
        sum += (f-q)*(f-q)*factor;
      }

    // accumulate cell indicator
    auto h_T = diameter(eg.geometry());
    r.accumulate(lfsv,0,h_T*h_T*sum);
  }

  // skeleton integral depending on test and ansatz functions
  // each face is only visited ONCE!
  template<typename IG, typename LFSU, typename X,
    typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
      const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
      const LFSU& lfsu_o, const X& x_o, const LFSV& lfsv_o,
      R& r_i, R& r_o) const
  {
    // geometries in local coordinates of the elements
    auto insidegeo = ig.geometryInInside();
    auto outsidegeo = ig.geometryInOutside();

    // inside and outside cells
    auto cell_inside = ig.inside();
    auto cell_outside = ig.outside();

    // geometries from local to global in elements
    auto geo_i = cell_inside.geometry();
    auto geo_o = cell_outside.geometry();

    // dimensions
    const int dim = IG::Entity::dimension;

    // select quadrature rule
    auto globalgeo = ig.geometry();
    const int order =
      incrementorder+2*lfsu_i.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(globalgeo,order);

    // loop over quadrature points and integrate normal flux
    typedef decltype(makeZeroBasisFieldValue(lfsu_i)) RF;
    RF sum(0.0);
    for (const auto& ip : rule)
      {
        // position of quadrature point in local coordinates of elements
        auto iplocal_i = insidegeo.global(ip.position());
        auto iplocal_o = outsidegeo.global(ip.position());

        // unit outer normal direction
        auto n_F = ig.unitOuterNormal(ip.position());

        // gradient in normal direction in self
        auto& gradphihat_i = cache.evaluateJacobian(
            iplocal_i,lfsu_i.finiteElement().localBasis());
        const auto S_i = geo_i.jacobianInverseTransposed(iplocal_i);
        RF gradun_i = 0.0;
        for (size_t i=0; i<lfsu_i.size(); i++)
          {
            Dune::FieldVector<RF,dim> v;
            S_i.mv(gradphihat_i[i][0],v);
            gradun_i += x_i(lfsu_i,i)*(v*n_F);
          }

        // gradient in normal direction in neighbor
        auto& gradphihat_o = cache.evaluateJacobian(
            iplocal_o,lfsu_o.finiteElement().localBasis());
        const auto S_o = geo_o.jacobianInverseTransposed(iplocal_o);
        RF gradun_o = 0.0;
        for (size_t i=0; i<lfsu_o.size(); i++)
          {
            Dune::FieldVector<RF,dim> v;
            S_o.mv(gradphihat_o[i][0],v);
            gradun_o += x_o(lfsu_o,i)*(v*n_F);
          }

        // integrate
        RF factor =
          ip.weight()*globalgeo.integrationElement(ip.position());
        RF jump = gradun_i-gradun_o;
        sum += jump*jump*factor;
      }

    // accumulate indicator
    auto h_T = diameter(globalgeo);
    r_i.accumulate(lfsv_i,0,0.5*h_T*sum);
    r_o.accumulate(lfsv_o,0,0.5*h_T*sum);
  }

  // boundary integral depending on test and ansatz functions
  // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                       const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
                       R& r_i) const
  {
    // geometries in local coordinates of the elements
    auto insidegeo = ig.geometryInInside();

    // inside and outside cells
    auto cell_inside = ig.inside();

    // geometries from local to global in elements
    auto geo_i = cell_inside.geometry();

    // dimensions
    const int dim = IG::Entity::dimension;

    // select quadrature rule
    auto globalgeo = ig.geometry();
    const int order = incrementorder+2*lfsu_i.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(globalgeo,order);

    // loop over quadrature points and integrate normal flux
    typedef decltype(makeZeroBasisFieldValue(lfsu_i)) RF;
    RF sum(0.0);
    for (const auto& ip : rule)
      {
        // skip body if we are on Dirichlet boundary
        if (param.b(ig.intersection(),ip.position())) continue;

        // position of quadrature point in local coordinates of elements
        auto iplocal_i = insidegeo.global(ip.position());

        // unit outer normal direction
        auto n_F = ig.unitOuterNormal(ip.position());

        // gradient in normal direction in self
        auto& gradphihat_i = cache.evaluateJacobian(
            iplocal_i,lfsu_i.finiteElement().localBasis());
        const auto S_i = geo_i.jacobianInverseTransposed(iplocal_i);
        RF gradun_i = 0.0;
        for (size_t i=0; i<lfsu_i.size(); i++)
          {
            Dune::FieldVector<RF,dim> v;
            S_i.mv(gradphihat_i[i][0],v);
            gradun_i += x_i(lfsu_i,i)*(v*n_F);
          }

        // Neumann boundary condition value
        auto j = param.j(ig.intersection(),ip.position());

        // integrate
        RF factor =
          ip.weight()*globalgeo.integrationElement(ip.position());
        RF jump = gradun_i+j;
        sum += jump*jump*factor;
      }

    // accumulate indicator
    auto h_T = diameter(globalgeo);
    r_i.accumulate(lfsv_i,0,h_T*sum);
  }
};

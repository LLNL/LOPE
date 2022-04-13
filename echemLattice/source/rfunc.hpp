#pragma once

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

/**
 * This header file defines a series of methods to generate simple lattice level sets and expressions to use with an automatic-differentiation framework.
 *
 * The assumption in this implementation, is that AD will be done by reference w.r.t to a particular parameter. This has different ramifications and considerations. The functions themselves are not dependent on the byref assumption.
 * The byref assumption means, that any primitive that depends on it's input types should be handled by references. 
 */

/** templated point types **/
template <typename T, int dim>
using PointType = Eigen::Matrix<T, dim, 1>;

/** explicit point types **/
template <typename T>
using Point2D = PointType<T, 2>;

template <typename T>
using Point3D = PointType<T, 3>;

/** math functions **/
namespace math {
  template <typename T, int p>
  struct Pow {
    constexpr T operator()(T value) { return value * Pow<T, p-1>()(value); }
  };

  template <typename T>
  struct Pow<T, 1> {
    constexpr T operator()(T value) { return value; }
  };

  template <typename T, int p>
  struct Norm {
    constexpr T operator()(T f1, T f2) { return pow(Pow<T,p>()(f1) + Pow<T,p>()(f2), 1./p); }
  };

  template <typename T>
  struct Norm<T, 2> {
    constexpr T operator()(T f1, T f2) { return sqrt(Pow<T,2>()(f1) + Pow<T,2>()(f2)); }
  };

}
/** end math functions **/

template <typename FuncType>
auto GenerateFunction(FuncType & f) {
  return [&] (auto x) { return f(x); };
}

template <typename T, int dim>
struct Function {
  static constexpr int Dimension = dim;
  using type = T;
  virtual T operator()(PointType<T, dim> p) {return 0;}
  virtual ~Function() {};
};

template <typename FuncType, typename GradType>
struct Normalize : Function<typename FuncType::type, FuncType::Dimension> {
  using T = typename FuncType::type;
  using P = PointType<T, FuncType::Dimension>;
  FuncType func;
  GradType grad; // copy by value
  
  Normalize(FuncType f, GradType g) : func(f), grad(g) {}

  T operator()(P & point) {
    T f = func(point); // f @ point
    return f / sqrt(f * f + grad.dot(grad));
  }
};

template <typename FuncType>
struct SmoothHeaviside : Function<typename FuncType::type, FuncType::Dimension> {
  using T = typename FuncType::type;
  using P = PointType<T, FuncType::Dimension>;
  FuncType func;  
  T & k;
  
  SmoothHeaviside (T& penalty, FuncType f) : func(f), k(penalty) {}

  T operator()(P point) {
    return 0.5 * (1. + tanh(k*func(point)));
    //return tanh(k*func(point));
  }
};

template <typename FuncType1, typename FuncType2>
struct Mult : Function<typename FuncType1::type, FuncType1::Dimension> {
  using T = typename FuncType1::type;
  using P = PointType<T, FuncType1::Dimension>;
  FuncType1 func1;
  FuncType2 func2;

  Mult(FuncType1 f1, FuncType2 f2) : func1(f1), func2(f2) {}
  T operator()(P point) override{
    return func1(point) * func2(point);
  }
    
};

template <typename FuncType1, typename ... FuncTypes>
struct Average : Function<typename FuncType1::type, FuncType1::Dimension> {
  using T = typename FuncType1::type;
  using P = PointType<T, FuncType1::Dimension>;
  FuncType1 func1;
  std::tuple<FuncTypes...> funcs;

  Average(FuncType1 f1, FuncTypes... types) : func1(f1), funcs ( types...) {};

  T operator()(P point) override {
    T sum = func1(point);
    int n = 1;
    std::apply([&](auto &... fs) { ((sum += fs(point), n++), ...); }, funcs);
    return sum/n;
  };
  
};

/// R-function namespace
namespace R{
  template <typename FuncType1, typename FuncType2>
  struct And : Function<typename FuncType1::type, FuncType1::Dimension> {
    using T = typename FuncType1::type;
    using P = PointType<T, FuncType1::Dimension>;
    FuncType1 func1;
    FuncType2 func2;

    And(FuncType1 f1, FuncType2 f2) : func1(f1), func2(f2) {}

    T operator()(P point) override {      
      T f1_eval = func1(point);
      T f2_eval = func2(point);
      return f1_eval + f2_eval - math::Norm<T, 2>()(f1_eval, f2_eval);
    };
  };

  template <typename FuncType1, typename FuncType2>
  struct Or : Function<typename FuncType1::type, FuncType1::Dimension> {
    using T = typename FuncType1::type;
    using P = PointType<T, FuncType1::Dimension>;
    FuncType1 func1;
    FuncType2 func2;

    Or(FuncType1 f1, FuncType2 f2) : func1(f1), func2(f2) {}

    T operator()(P point) override {
      T f1_eval = func1(point);
      T f2_eval = func2(point);      
      return f1_eval + f2_eval + math::Norm<T, 2>()(f1_eval, f2_eval);
    };
  };

  template <typename FuncType1, typename FuncType2>
  struct Minus : Function<typename FuncType1::type, FuncType1::Dimension> {
    using T = typename FuncType1::type;
    using P = PointType<T, FuncType1::Dimension>;
    FuncType1 func1;
    FuncType2 func2;

    Minus(FuncType1 f1, FuncType2 f2) : func1(f1), func2(f2) {}

    T operator()(P point) override{
      T f1_eval = func1(point);
      T f2_eval = func2(point);            
      return f1_eval - f2_eval - math::Norm<T, 2>()(f1_eval, f2_eval);
    };
  };

} // Namespace R


/// all level set primitivies
namespace primitive {

  /**
   * @brief Creates a levelset circle
   */
  template <typename T>
  struct Circle : Function <T, 2>{
    struct Params {
      T & radius;
      Point2D<T> & center;
    };

    Params params;
    Circle(T & r, Point2D<T> & c) : params{r, c} {};
  
    T operator()(Point2D<T> point) { 
      return  1. - math::Pow<T,2>()((point - params.center).norm()) / math::Pow<T,2>()(params.radius);
    }
  };

  /**
   * @brief Creates a level set of a sphere
   */
  template <typename T>
  struct Sphere : Function <T, 3>{
    struct Params {
      T & radius;
      Point3D<T> & center;
    };

    Params params;
    Sphere(T & r, Point3D<T> & c) : params{r, c} {};
  
    T operator()(Point3D<T> point) { 
      //return  1. - math::Pow<T,2>()((point - params.center).norm()) / math::Pow<T,2>()(params.radius);
      // return  1. - (math::Pow<T,2>()(point[0] - params.center[0])
      //              +math::Pow<T,2>()(point[1] - params.center[1])
      //              +math::Pow<T,2>()(point[2] - params.center[2]))
      //               / math::Pow<T,2>()(params.radius);
      return  math::Pow<T,4>()(params.radius) - (math::Pow<T,4>()(point[0] - params.center[0])
                   +math::Pow<T,4>()(point[1] - params.center[1])
                   +math::Pow<T,4>()(point[2] - params.center[2]))
                    ;
      
    }
  };

  
  /**
   * @brief Creates a level set of a line with a radius in 2D or 3D
   *
   * @param[in] r radius
   * @param[in] point1 first point on the line
   * @param[in] point2 second point on the line
   */
  template <typename T, int dim>
  struct Line : Function<T, dim> {
    using Point = PointType<T, dim>;
    struct Params {
      Point &p1;
      Point &p2;
      T &radius;
    };

    Point v_line;
    Params params;
    Line(T & r, Point & point1, Point & point2) : params { point1, point2, r} {
    }

    T operator()(Point p) {
      // find the vector from p1 to p2
      v_line = params.p2 - params.p1;
      v_line = v_line.normalized();
      
      // find the vector from p1 to p
      auto v_p = p - params.p1;
      // v_p project onto v_line
      auto projV = v_line.dot(v_p);
      // find closest projection point on v_line
      auto closest = projV * v_line;
      // calculate the distance
      auto d = (v_p - closest).norm();
      //      return 1. - math::Pow<T,2>()(d/params.radius);
      return math::Pow<T,2>()(params.radius) - math::Pow<T,2>()(d);
    }
  };

  // Creates a band level set that is positive perpendicular to the line from p1 to p2
  template <typename T, int dim>
  struct Band : Function<T, dim> {
    using Point = PointType<T, dim>;
    struct Params {
      Point &p1;
      Point &p2;
    };

    Point v_line;
    T center;
    Params params;
    Band(Point &point1, Point &point2) : params {point1, point2} {
    }

    T operator()(Point p) {
      // find the vector from p1 to p2
      v_line = params.p2 - params.p1;
      // find the center of the line along v_line
      center = v_line.norm() * 0.5;     
      v_line = v_line.normalized();
      
      // construct v_p from p1 to p
      auto v_p = p - params.p1;
      // find the projection of the point on the line
      auto proj = v_p.dot(v_line);

      //find the center of v_line to calculate the distance
      auto d = 1. -math::Pow<T,2>()((proj - center)/center);
      return d;
    }
  };

  // Constructs a rod from p1 to p2 with radius, r
  template <typename T, int dim>
  struct Rod : Function<T, dim> {
    using Point = PointType<T, dim>;
    std::shared_ptr<Function<T, dim>> combined;

    Rod(Band<T, dim> band, Line<T, dim> line) {
      combined = std::make_shared<R::And<decltype(band), decltype(line)>>(band, line);
    }
    Rod(T &r,  Point & point1, Point & point2) {
      combined = std::make_shared<R::And<Band<T, dim>, Line<T, dim>>>(Band(point1, point2), Line(r, point1, point2));
      // combined = std::make_shared<Mult<Band<T, dim>, Line<T, dim>>>(Band(point1, point2), Line(r, point1, point2));
    }

    T operator()(Point p) {
      return (*combined)(p);
    }    
  };

  // Constructs a WindowBox aligned with x, y, and z
  template <typename T, int dim>
  struct Window : Function <T, dim> {
    using Point = PointType<T, dim>;
    
    std::shared_ptr<Function<T, dim>> combined;
    
    Window(Point p1, Point p2) {
      if constexpr(dim == 1) {
	  // band is the same as window
	  combined = std::make_shared<Band<T, dim>>(p1, p2);
	  
	} else if constexpr(dim == 2) {
	  Point p1_x; p1_x << p1[0] , 0;
	  Point p2_x; p2_x << p2[0] , 0;
	  auto band_x = Band<T, dim>(p1_x, p2_x);
	  Point p1_y; p1_y << 0, p1[1];
	  Point p2_y; p2_y << 0, p2[1];
	  auto band_y = Band<T, dim>(p1_y, p2_y);
	  combined = std::make_shared<R::And>(band_x, band_y);
	  
	} else if constexpr(dim==3) {
	  Point p1_x; p1_x << p1[0] , 0 , 0;
	  Point p2_x; p2_x << p2[0] , 0 , 0;
	  auto band_x = Band<T, dim>(p1_x, p2_x);
	  Point p1_y; p1_y << 0, p1[1], 0;
	  Point p2_y; p2_y << 0, p2[1], 0;
	  auto band_y = Band<T, dim>(p1_y, p2_y);
	  Point p1_z; p1_z << 0, 0, p1[2];
	  Point p2_z; p2_z << 0, 0, p2[2];
	  auto band_z = Band<T, dim>(p1_z, p2_z);
	  
	  combined = std::make_shared<R::And>(R::And(band_x, band_y), band_z);	  
	} else {
	// TODO: error
      }
    }

    T operator()(Point p) {
      return (*combined)(p);
    }    

  };
  
} // primitive namespace

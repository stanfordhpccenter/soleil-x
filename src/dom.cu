#include <assert.h>
#include <math.h>
#include <fstream>
#include <stdio.h>

#include <thrust/adjacent_difference.h>
#include <thrust/copy.h>
#include <thrust/device_delete.h>
#include <thrust/device_new.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/iterator/iterator_facade.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/scan.h>
#include <thrust/tuple.h>

typedef std::size_t size_t;
typedef thrust::tuple<double,double> Double2;
typedef thrust::tuple<double,double,double,double> Double4;
using thrust::get;

//=============================================================================
// HELPER CLASSES
//=============================================================================

class DeviceArray;

class HostArray {
  friend class DeviceArray;
public:
  HostArray(size_t nx, size_t ny, size_t nz)
    : nx_(nx),
      ny_(ny),
      nz_(nz),
      data_(new double[nx * ny * nz]) {}
  ~HostArray() {
    delete [] data_;
  }
public:
  HostArray(const HostArray&) = delete;
  HostArray& operator=(const HostArray&) = delete;
  void copy_from(const DeviceArray& src);
  void read_from(const char* filename) {
    std::fstream fin(filename, std::fstream::in);
    for (size_t i = 0; i < nx_; ++i) {
      for (size_t j = 0; j < ny_; ++j) {
        for (size_t k = 0; k < nz_; ++k) {
          fin >> (*this)(i,j,k);
        }
      }
    }
  }
  void write_to(const char* filename) {
    FILE* fout = fopen(filename, "w");
    if (fout == NULL) {
      perror("Error opening file");
      exit(1);
    }
    for (size_t i = 0; i < nx_; ++i) {
      for (size_t j = 0; j < ny_; ++j) {
        for (size_t k = 0; k < nz_; ++k) {
          fprintf(fout, "%.17e\n", (*this)(i,j,k));
        }
      }
    }
    fclose(fout);
  }
public:
  double& operator()(size_t x, size_t y, size_t z) {
    return *(data_ + x*ny_*nz_ + y*nz_ + z);
  }
  const double& operator()(size_t x, size_t y, size_t z) const {
    return *(data_ + x*ny_*nz_ + y*nz_ + z);
  }
  size_t nx() const { return nx_; }
  size_t ny() const { return ny_; }
  size_t nz() const { return nz_; }
private:
  size_t nx_;
  size_t ny_;
  size_t nz_;
  double* data_;
};

class DeviceArray {
  friend class HostArray;
public:
  DeviceArray(size_t nx, size_t ny, size_t nz)
    : nx_(nx),
      ny_(ny),
      nz_(nz),
      data_(thrust::device_new<double>(nx * ny * nz)) {}
public:
  DeviceArray(const DeviceArray&) = delete;
  DeviceArray& operator=(const DeviceArray&) = delete;
  void copy_from(const HostArray& src) {
    assert(src.nx() == nx_ && src.ny() == ny_ && src.nz() == nz_);
    thrust::copy_n(src.data_, nx_ * ny_ * nz_, data_);
  }
public:
  thrust::device_reference<double>
  operator()(size_t x, size_t y, size_t z) {
    return *(data_ + x*ny_*nz_ + y*nz_ + z);
  }
  thrust::device_reference<const double>
  operator()(size_t x, size_t y, size_t z) const {
    return *(data_ + x*ny_*nz_ + y*nz_ + z);
  }
  thrust::device_ptr<double> z_row(size_t x, size_t y) {
    return data_ + x*ny_*nz_ + y*nz_;
  }
  thrust::device_ptr<const double> z_row(size_t x, size_t y) const {
    return data_ + x*ny_*nz_ + y*nz_;
  }
  size_t nx() const { return nx_; }
  size_t ny() const { return ny_; }
  size_t nz() const { return nz_; }
private:
  size_t nx_;
  size_t ny_;
  size_t nz_;
  thrust::device_ptr<double> data_;
};

void HostArray::copy_from(const DeviceArray& src) {
  assert(src.nx() == nx_ && src.ny() == ny_ && src.nz() == nz_);
  thrust::copy_n(src.data_, nx_ * ny_ * nz_, data_);
}

// template<typename ValueT, typename IterT> class with_header_iterator
//   : public thrust::iterator_facade<with_header_iterator<ValueT,IterT>,ValueT> {
// public:
//   typedef with_header_iterator<ValueT,IterT> SelfT;
//   typedef thrust::iterator_facade<with_header_iterator<ValueT,IterT>,
//                                   ValueT> SuperT;
//   typedef SuperT::difference_type DiffT;
//   friend class thrust::iterator_core_access;
// public:
//   __host__ __device__
//   with_header_iterator(ValueT header, IterT iter)
//     : header_(std::move(header)),
//       pos_(0),
//       iter_(std::move(iter)) {}
// private:
//   __host__ __device__
//   const ValueT& dereference() const {
//     if (pos_ == 0) {
//       return header_;
//     }
//     return *iter_;
//   }
//   __host__ __device__
//   bool equal(const SelfT& other) const {
//     return
//       // header_ == other.header_ &&
//       // iter_ == other.iter_ &&
//       pos_ == other.pos_;
//   }
//   __host__ __device__
//   void increment() {
//     if (++pos_ > 1) {
//       ++iter_;
//     }
//   }
//   __host__ __device__
//   void decrement() {
//     if (--pos_ >= 1) {
//       --iter_;
//     }
//   }
//   __host__ __device__
//   void advance(DiffT n) {
//     if (pos_ == 0) {
//       pos_++;
//       n--;
//     }
//     iter_ += n;
//   }
//   __host__ __device__
//   void distance_to(const SelfT& other) const {
//     return other.pos_ - pos_;
//   }
// private:
//   const ValueT header_;
//   DiffT pos_;
//   IterT iter_;
// };

// template<typename ValueT, typename IterT>
// add_header(ValueT header, IterT iter) {
//   return with_header_iterator(header, iter);
// }

//=============================================================================
// ALGORITHM
//=============================================================================

// Start with general TWD scheme:
// cell_int[i,j,k] = (cell_source[i,j,k] * dV
//                    + fabs(xi)  * dAx / x_gamma * x_face_int[i,j,k]
//                    + fabs(eta) * dAy / y_gamma * y_face_int[i,j,k]
//                    + fabs(mu)  * dAz / z_gamma * z_face_int[i,j,k])
//                 / (cell_sigma[i,j,k] * dV
//                    + fabs(xi)  * dAx / x_gamma
//                    + fabs(eta) * dAy / y_gamma
//                    + fabs(mu)  * dAz / z_gamma)
// x_face_int[i+1,j,k] = cell_int[i,j,k] / x_gamma
//                     - (1-x_gamma)/x_gamma * x_face_int[i,j,k]
// y_face_int[i,j+1,k] = cell_int[i,j,k] / y_gamma
//                     - (1-y_gamma)/y_gamma * y_face_int[i,j,k]
// z_face_int[i,j,k+1] = cell_int[i,j,k] / z_gamma
//                     - (1-z_gamma)/z_gamma * z_face_int[i,j,k]

// Select a specific z-row (i,j),
// treat x_face_int[i,j,*], y_face_int[i,j,*], z_face_int[i,j,0] as constants,
// reformulate in terms of z_face_int:
// a(k) = fabs(mu) * dAz / z_gamma
//      / (cell_sigma[i,j,k] * dV
//         + fabs(xi)  * dAx / x_gamma
//         + fabs(eta) * dAy / y_gamma
//         + fabs(mu)  * dAz / z_gamma)
//      / z_gamma
//      - (1-z_gamma)/z_gamma
// b(k) = (cell_source[i,j,k] * dV
//         + fabs(xi)  * dAx / x_gamma * x_face_int[i,j,k]
//         + fabs(eta) * dAy / y_gamma * y_face_int[i,j,k])
//      / (cell_sigma[i,j,k] * dV
//         + fabs(xi)  * dAx / x_gamma
//         + fabs(eta) * dAy / y_gamma
//         + fabs(mu)  * dAz / z_gamma)
//      / z_gamma
// z_face_int[i,j,k+1] = a(k) * z_face_int[i,j,k] + b(k)

// Other quantities can be computed based on z_face_int[i,j,k+1]:
// cell_int[i,j,k] = (1-z_gamma) * z_face_int[i,j,k]
//                 + z_gamma * z_face_int[i,j,k+1]
// x_face_int[i+1,j,k] = cell_int[i,j,k] / x_gamma
//                     - (1-x_gamma)/x_gamma * x_face_int[i,j,k]
// y_face_int[i,j+1,k] = cell_int[i,j,k] / y_gamma
//                     - (1-y_gamma)/y_gamma * y_face_int[i,j,k]

// To solve the recurrence relation for z_face_int:
// Create a vector of pairs:
// (any,z_face_int[i,j,0]) (a(0),b(0)) ... (a(NZ-1),b(NZ-1))
// Do a prefix sum using operator:
// (a,b) X (c,d) = ( a*c, (b*c)+d )
// Take the 2nd element of each pair.

struct A : public thrust::unary_function<double,double> {
  __host__ __device__
  double operator()(double sigma) const {
    double dAx = dy_*dz_;
    double dAy = dx_*dz_;
    double dAz = dx_*dy_;
    double dV = dx_*dy_*dz_;
    return
      fabs(mu_) * dAz / z_gamma_
      / (sigma * dV
         + fabs(xi_)  * dAx / x_gamma_
         + fabs(eta_) * dAy / y_gamma_
         + fabs(mu_)  * dAz / z_gamma_)
      / z_gamma_
      - (1-z_gamma_)/z_gamma_;
  }
  __host__ __device__
  A(double dx, double dy, double dz,
    double xi, double eta, double mu,
    double x_gamma, double y_gamma, double z_gamma)
    : dx_(dx), dy_(dy), dz_(dz),
      xi_(xi), eta_(eta), mu_(mu),
      x_gamma_(x_gamma), y_gamma_(y_gamma), z_gamma_(z_gamma) {}
  double dx_, dy_, dz_;
  double xi_, eta_, mu_;
  double x_gamma_, y_gamma_, z_gamma_;
};

struct B : public thrust::unary_function<Double4,double> {
  __host__ __device__
  double operator()(const Double4& args) const {
    double cell_source = get<0>(args);
    double x_face_int = get<1>(args);
    double y_face_int = get<2>(args);
    double cell_sigma = get<3>(args);
    double dAx = dy_*dz_;
    double dAy = dx_*dz_;
    double dAz = dx_*dy_;
    double dV = dx_*dy_*dz_;
    return
      (cell_source * dV
       + fabs(xi_)  * dAx / x_gamma_ * x_face_int
       + fabs(eta_) * dAy / y_gamma_ * y_face_int)
      / (cell_sigma * dV
         + fabs(xi_)  * dAx / x_gamma_
         + fabs(eta_) * dAy / y_gamma_
         + fabs(mu_)  * dAz / z_gamma_)
      / z_gamma_;
  }
  __host__ __device__
  B(double dx, double dy, double dz,
    double xi, double eta, double mu,
    double x_gamma, double y_gamma, double z_gamma)
    : dx_(dx), dy_(dy), dz_(dz),
      xi_(xi), eta_(eta), mu_(mu),
      x_gamma_(x_gamma), y_gamma_(y_gamma), z_gamma_(z_gamma) {}
  double dx_, dy_, dz_;
  double xi_, eta_, mu_;
  double x_gamma_, y_gamma_, z_gamma_;
};

struct X : public thrust::binary_function<Double2,Double2,Double2> {
  __host__ __device__
  Double2 operator()(const Double2& lhs, const Double2& rhs) const {
    return thrust::make_tuple(get<0>(lhs) * get<0>(rhs),
                              get<1>(lhs) * get<0>(rhs) + get<1>(rhs));
  }
};

struct FaceInt2CellInt : public thrust::unary_function<Double2,double> {
  __host__ __device__
  double operator()(const Double2& args) const {
    double face_int_prev = get<0>(args);
    double face_int_next = get<1>(args);
    return (1-gamma_) * face_int_prev + gamma_ * face_int_next;
  }
  __host__ __device__
  FaceInt2CellInt(double gamma) : gamma_(gamma) {}
  double gamma_;
};

struct CellInt2FaceInt : public thrust::unary_function<Double2,double> {
  __host__ __device__
  double operator()(const Double2& args) const {
    double cell_int = get<0>(args);
    double face_int_prev = get<1>(args);
    return cell_int / gamma_ - (1-gamma_)/gamma_ * face_int_prev;
  }
  __host__ __device__
  CellInt2FaceInt(double gamma) : gamma_(gamma) {}
  double gamma_;
};

struct Project2nd : public thrust::unary_function<Double2,double> {
  __host__ __device__
  double operator()(const Double2& tup) const {
    return get<1>(tup);
  }
};

void sweep(DeviceArray/*(NX,NY,NZ)*/& cell_int,
           const DeviceArray/*(NX,NY,NZ)*/& cell_source,
           const DeviceArray/*(NX,NY,NZ)*/& cell_sigma,
           DeviceArray/*(NX+1,NY,NZ)*/& x_face_int,
           DeviceArray/*(NX,NY+1,NZ)*/& y_face_int,
           DeviceArray/*(NX,NY,NZ+1)*/& z_face_int,
           double dx, double dy, double dz,
           double xi, double eta, double mu,
           double x_gamma, double y_gamma, double z_gamma) {
  size_t NX = cell_int.nx();
  size_t NY = cell_int.ny();
  size_t NZ = cell_int.nz();
  // Construct function objects
  A a(dx, dy, dz, xi, eta, mu, x_gamma, y_gamma, z_gamma);
  B b(dx, dy, dz, xi, eta, mu, x_gamma, y_gamma, z_gamma);
  X x;
  FaceInt2CellInt z2c(z_gamma);
  CellInt2FaceInt c2x(x_gamma);
  CellInt2FaceInt c2y(y_gamma);
  // Iterate over all z-rows, in row-major order
  for (size_t i = 0; i < NX; ++i) {
    for (size_t j = 0; j < NY; ++j) {
      // Fill in pairs
      thrust::device_vector<Double2> pairs(NZ+1);
      pairs[0] = thrust::make_tuple(-1.0, z_face_int(i,j,0));
      thrust::copy_n(
        thrust::make_zip_iterator(thrust::make_tuple(
          thrust::make_transform_iterator(
            cell_sigma.z_row(i,j),
            a),
          thrust::make_transform_iterator(
            thrust::make_zip_iterator(thrust::make_tuple(
              cell_source.z_row(i,j),
              x_face_int.z_row(i,j),
              y_face_int.z_row(i,j),
              cell_sigma.z_row(i,j))),
            b))),
        NZ,
        pairs.begin() + 1);
      // Run prefix sum over pairs
      thrust::inclusive_scan(pairs.begin(), pairs.end(), pairs.begin(), x);
      // Retrieve z-face intensity values
      thrust::copy_n(
        thrust::make_transform_iterator(
          pairs.begin() + 1,
          Project2nd()),
        NZ,
        z_face_int.z_row(i,j) + 1);
      // Calculate cell intensity values
      thrust::copy_n(
        thrust::make_transform_iterator(
          thrust::make_zip_iterator(thrust::make_tuple(
            z_face_int.z_row(i,j),
            z_face_int.z_row(i,j) + 1)),
          z2c),
        NZ,
        cell_int.z_row(i,j));
      // Calculate downstream x-face intensity values
      thrust::copy_n(
        thrust::make_transform_iterator(
          thrust::make_zip_iterator(thrust::make_tuple(
            cell_int.z_row(i,j),
            x_face_int.z_row(i,j))),
          c2x),
        NZ,
        x_face_int.z_row(i+1,j));
      // Calculate downstream y-face intensity values
      thrust::copy_n(
        thrust::make_transform_iterator(
          thrust::make_zip_iterator(thrust::make_tuple(
            cell_int.z_row(i,j),
            y_face_int.z_row(i,j))),
          c2y),
        NZ,
        y_face_int.z_row(i,j+1));
    }
  }
}

//=============================================================================
// MAIN
//=============================================================================

int main() {
  size_t NX, NY, NZ;
  double x_gamma, y_gamma, z_gamma;
  double dx, dy, dz;
  double xi, eta, mu;

  {
    std::fstream setup("setup.dat", std::fstream::in);
    setup >> NX >> NY >> NZ;
    setup >> x_gamma >> y_gamma >> z_gamma;
    setup >> dx >> dy >> dz;
    setup >> xi >> eta >> mu;
  }

  HostArray h_cell_int(NX,NY,NZ);
  HostArray h_cell_source(NX,NY,NZ);
  HostArray h_cell_sigma(NX,NY,NZ);
  HostArray h_x_face_int(NX+1,NY,NZ);
  HostArray h_y_face_int(NX,NY+1,NZ);
  HostArray h_z_face_int(NX,NY,NZ+1);

  DeviceArray d_cell_int(NX,NY,NZ);
  DeviceArray d_cell_source(NX,NY,NZ);
  DeviceArray d_cell_sigma(NX,NY,NZ);
  DeviceArray d_x_face_int(NX+1,NY,NZ);
  DeviceArray d_y_face_int(NX,NY+1,NZ);
  DeviceArray d_z_face_int(NX,NY,NZ+1);

  h_cell_int.read_from("cell_int_prev.dat");
  h_cell_source.read_from("cell_source.dat");
  h_cell_sigma.read_from("cell_sigma.dat");
  h_x_face_int.read_from("x_face_int_prev.dat");
  h_y_face_int.read_from("y_face_int_prev.dat");
  h_z_face_int.read_from("z_face_int_prev.dat");

  d_cell_int.copy_from(h_cell_int);
  d_cell_source.copy_from(h_cell_source);
  d_cell_sigma.copy_from(h_cell_sigma);
  d_x_face_int.copy_from(h_x_face_int);
  d_y_face_int.copy_from(h_y_face_int);
  d_z_face_int.copy_from(h_z_face_int);

  sweep(d_cell_int,
        d_cell_source,
        d_cell_sigma,
        d_x_face_int,
        d_y_face_int,
        d_z_face_int,
        dx, dy, dz,
        xi, eta, mu,
        x_gamma, y_gamma, z_gamma);

  h_cell_int.copy_from(d_cell_int);
  h_x_face_int.copy_from(d_x_face_int);
  h_y_face_int.copy_from(d_y_face_int);
  h_z_face_int.copy_from(d_z_face_int);

  h_cell_int.write_to("cell_int_cu.dat");
  h_x_face_int.write_to("x_face_int_cu.dat");
  h_y_face_int.write_to("y_face_int_cu.dat");
  h_z_face_int.write_to("z_face_int_cu.dat");
}

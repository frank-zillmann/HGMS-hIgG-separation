#include <Eigen/Dense>
#include <stdexcept>

template <typename VectorType>
class LinearInterpolator {
   public:
    LinearInterpolator(const VectorType& x, const VectorType& y) : x(x), y(y) {
        if (x.size() != y.size() || x.size() < 2) {
            throw std::invalid_argument("Vectors must have the same size and at least two elements.");
        }
    }

    typename VectorType::Scalar operator()(typename VectorType::Scalar x_query) const {
        // interpolate within bounds
        for (Eigen::Index i = 0; i < x.size() - 1; ++i) {
            if (x_query >= x(i) && x_query <= x(i + 1)) {
                auto t = (x_query - x(i)) / (x(i + 1) - x(i));
                return y(i) + t * (y(i + 1) - y(i));
            }
        }
        // return lowest value below bounds
        if (x_query < x(0)) {
            return y(0);
        }
        // return highest value above bounds
        if (x_query > x(x.size() - 1)) {
            return y(y.size() - 1);
        }
        throw std::runtime_error("Interpolation failed: could not find interval.");
    }

   private:
    VectorType x, y;
};

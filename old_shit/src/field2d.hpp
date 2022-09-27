#ifndef FLUID_SOLVER_RAW_FIELD2D_HPP
#define FLUID_SOLVER_RAW_FIELD2D_HPP

#include <functional>
#include <Eigen/Core>
#include <cfloat>

#include "macros.hpp"

#define INTEGRAL_IN_RANGE(V, MIN, MAX){ (V) < (MIN) ? (MIN) : ((V) > (MAX) ? (MAX) : (V))}
#define COORD_TO_INDEX(ARR, X, Y, YS) ARR[((X) * (YS)) + (Y)]

using namespace Eigen;
using namespace std;
using namespace std::placeholders;

namespace field_traits {

    /**
     * Type dependant value for zero
     * @tparam T
     */
    template<typename T>
    struct default_zero_value;

    template<>
    struct default_zero_value<float> {
        static const float value() { return 0.0f; }
    };

    template<>
    struct default_zero_value<double> {
        static const float value() { return 0.0; }
    };


    template<>
    struct default_zero_value<int> {
        static const float value() { return 0; }
    };


    template<>
    struct default_zero_value<Eigen::Vector2f> {
        static const Eigen::Vector2f value() { return Eigen::Vector2f::Zero(); }
    };

    template<>
    struct default_zero_value<Eigen::Vector2d> {
        static const Eigen::Vector2d value() { return Eigen::Vector2d::Zero(); }
    };


    template<class T>
    T abs(const T v) {
        return v == std::numeric_limits<T>::min() ? v :
               (v < field_traits::default_zero_value<T>::value() ? -v : v);
    }

    /**
     * Enable type dependant value to be appened to std::ostream
     * @tparam T
     */

    template<typename T>
    struct append_to_ostream;

    template<>
    struct append_to_ostream<float> {
        static std::ostream &value(std::ostream &stream, float value) {
            return stream << value;
        }
    };

    template<>
    struct append_to_ostream<double> {
        static std::ostream &value(std::ostream &stream, double value) {
            return stream << value;
        }
    };

    template<>
    struct append_to_ostream<int> {
        static std::ostream &value(std::ostream &stream, int value) {
            using std::operator<<;
            return stream << value;
        }
    };

    template<>
    struct append_to_ostream<Eigen::Vector2f> {
        static std::ostream &value(std::ostream &stream, Eigen::Vector2f &value) {
            using std::operator<<;
            return stream << "(" << value.x() << "," << value.y() << ")";
        }
    };


    template<>
    struct append_to_ostream<Eigen::Vector2d> {
        static std::ostream &value(std::ostream &stream, Eigen::Vector2d &value) {
            using std::operator<<;
            return stream << "(" << value.x() << "," << value.y() << ")";
        }
    };

}

template<typename T>
class Field2D {
private:
    T *q;
    int q_size;
protected:
    int length_x;
    int length_y;
    int size;
public:

    Field2D(int length_x, int length_y) {
        this->q_size = sizeof(T);
        this->length_x = length_x;
        this->length_y = length_y;
        this->size = length_x * length_y;
        this->q = new T[this->length_x * this->length_y];
        this->forEach([](auto &cell, int x, int y) {
            cell = field_traits::default_zero_value<T>::value();
        });
    }

    virtual ~Field2D() {
        delete[]this->q;
    }

    T inline &operator()(const Eigen::Vector2i &p) const {
        return this->operator()(p.x(), p.y());
    }

    T inline &operator()(const Eigen::Vector2f &p) const {
        return this->operator()(p.x(), p.y());
    }

    T inline &operator()(const Eigen::Vector2d &p) const {
        return this->operator()(p.x(), p.y());
    }

    T inline &operator()(double x, double y) const {
        return this->operator()(static_cast<int>(std::round(x)), static_cast<int>(std::round(y)));
    }

    T inline &operator()(int x, int y) const {
        int _x = std::clamp(x, 0, this->length_x - 1);
        int _y = std::clamp(y, 0, this->length_y - 1);
        return this->q[_x + length_x * _y];
    }

    void inline set(const Eigen::Vector2i &p, const T &v) {
        this->set(p.x(), p.y(), v);
    }

    void inline set(const Eigen::Vector2f &p, const T &v) {
        return this->set(std::nearbyint(p.x()), std::nearbyint(p.y()), v);
    }

    void inline set(const Eigen::Vector2d &p, const T &v) {
        this->set(p.x(), p.y(), v);
    }

    void inline set(double x, double y, const T &v) {
        return this->set(std::nearbyint(x), std::nearbyint(y), v);
    }

    void inline set(int x, int y, const T &v) {
        int _x = std::clamp(x, 0, this->length_x - 1);
        int _y = std::clamp(y, 0, this->length_y - 1);
        this->q[_x + length_x * _y] = v;
    };

    int inline getLengthX() const { return length_x; }

    int inline getLengthY() const { return length_y; }

    int inline getSize() const { return size; }

    friend inline std::ostream &operator<<(std::ostream &output, const Field2D<T> &v) {
        auto cols = v.length_x;
        auto rows = v.length_y;
        auto cells = ((v.q_size * v.length_x * v.length_y) / v.q_size);
        auto cell_size = v.q_size;
        auto field_size = v.q_size * v.length_x * v.length_y;

        output << "Field2d(cols:" << cols << ", rows:" << rows << ", cells:" << cells << ", cell_size:" << cell_size
               << "b, field_size:" << field_size << "b) = {"
               << std::endl;

        for (int x = 0; x < v.length_x; x++) {
            for (int y = 0; y < v.length_y; y++) {
                auto cell = v(x, y);
                output << "\t";
                field_traits::append_to_ostream<T>::value(output, cell) << " ";
            }
            output << endl;
        }
        output << "}" << std::endl;
        return output;
    }

    friend inline std::ostream &operator<<(std::ostream &output, Field2D<T> *&v) {
        output << "// cell size in byte " << v->q_size << endl;
        output << "// field size in byte " << v->q_size * v->length_x * v->length_y << endl;
        output << "// number of cols " << v->length_y << endl;
        output << "// number of rows " << v->length_x << endl;
        output << "// number of cells " << (v->q_size * v->length_x * v->length_y) / v->q_size << endl;
        output << "F=[" << std::endl;
        for (int x = 0; x < v->length_x; x++) {
            for (int y = 0; y < v->length_y; y++) {
                auto cell = v->operator()(x, y);
                field_traits::append_to_ostream<T>::value(output, cell) << " ";
            }
            output << endl;
        }
        output << "]" << std::endl;
        return output;
    }

    void inline copyTo(Field2D<T> *&b) {

        assert(b->q_size == this->q_size && "Size of data-type for quantity must be equal!");
        assert(b->length_x == this->length_x && "Length in x-direction must be equal!");
        assert(b->length_y == this->length_y && "Length in y-direction must be equal!");

        assert(b->q_size > 0 && this->q_size > 0 && "Size of data-type for quantity must not be 0!");
        assert(b->length_x > 0 && this->length_x > 0 && "Length in x-direction must not be 0!");
        assert(b->length_y > 0 && this->length_y > 0 && "Length in y-direction must not be 0!");

        memcpy(b->q, this->q, this->size * this->q_size);

    }

    void inline print() {
        std::cout << (*this) << std::endl;
    }

    void inline
    forEach(const std::function<void(T &/* current cell */, int /* x coordinate */, int /* y coordinate */)> &f) {
// #pragma omp parallel for
        for (int x = 0; x < this->length_x; ++x) {
// #pragma omp parallel for
            for (int y = 0; y < this->length_y; ++y) {
                f(this->operator()(x, y), x, y);
            }
        }
    }

};

#endif //FLUID_SOLVER_RAW_FIELD2D_HPP

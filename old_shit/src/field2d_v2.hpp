#ifndef FLUID_SOLVER_RAW_FIELD2D_HPP
#define FLUID_SOLVER_RAW_FIELD2D_HPP

#include <functional>
#include <Eigen/Core>
#include <iostream>

#include "macros.hpp"
#include "math.hpp"

#define XY2I(x, y, xs) ((x)+(xs)*(y))

template <typename T> T CLAMP(const T& value, const T& low, const T& high)
{
    return value < low ? low : (value > high ? high : value);
}

using namespace Eigen;
using namespace std;
using namespace std::placeholders;
using namespace fluid::util::math::interpolation;

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
class Field2DV2 {
private:
    int width;
    int height;
    int size;
    int element_size;
protected:
public:

    T *_u;
    T *_v;

    Field2DV2(int width, int height) {
        this->width = width;
        this->height = height;
        this->size = width * height;
        this->element_size = sizeof(T);
        this->_u = new T[this->width * this->height];
        this->_v = new T[this->width * this->height];

        this->forEach([](auto &u_xy, auto &v_xy, int x, int y) {
            u_xy = field_traits::default_zero_value<T>::value();
            v_xy = field_traits::default_zero_value<T>::value();
        });
    }

    virtual ~Field2DV2() {
        delete[]this->_u;
        delete[]this->_v;
    }

    /* u : get / set */

    T inline u(const Eigen::Vector2i &p) const {
        return this->u(p.x(), p.y());
    }

    T inline &u(const Eigen::Vector2i &p) {
        return this->u(p.x(), p.y());
    }

    T inline u(const Eigen::Vector2f &p) const {
        return this->u(p.x(), p.y());
    }

    T inline &u(const Eigen::Vector2f &p) {
        return this->u(p.x(), p.y());
    }

    T inline u(const Eigen::Vector2d &p) const {
        return this->u(p.x(), p.y());
    }

    T inline &u(const Eigen::Vector2d &p) {
        return this->u(p.x(), p.y());
    }


    T inline u(double x, double y) const {
        return bridson_u(x,y);
    }

    T inline u(double x, double y) {
        return bridson_u(x,y);
    }

    T inline u(int x, int y) const {
        auto _x = CLAMP(x, 0, this->width-1);
        auto _y = CLAMP(y, 0, this->height-1);
        return this->_u[XY2I(_x,_y, this->width)];
    }

    T inline &u(int x, int y) {
        auto _x = CLAMP(x, 0, this->width-1);
        auto _y = CLAMP(y, 0, this->height-1);
        return this->_u[XY2I(_x,_y, this->width)];
    }

    T inline u(int i) const {
        auto _i = CLAMP(i, 0, this->size-1);
        return this->_u[_i];
    }

    T inline &u(int i)  {
        auto _i = CLAMP(i, 0, this->size-1);
        return this->_u[_i];
    }
    /* v : get / set */
    T inline v(const Eigen::Vector2i &p) const {
        return this->_v(p.x(), p.y());
    }

    T inline &v(const Eigen::Vector2i &p) {
        return this->_v(p.x(), p.y());
    }

    T inline v(const Eigen::Vector2f &p) const {
        return this->_v(p.x(), p.y());
    }

    T inline &v(const Eigen::Vector2f &p) {
        return this->_v(p.x(), p.y());
    }

    T inline v(const Eigen::Vector2d &p) const {
        return this->_v(p.x(), p.y());
    }

    T inline &v(const Eigen::Vector2d &p) {
        return this->_v(p.x(), p.y());
    }

    T inline v(double x, double y) const {
        return bridson_v(x,y);
    }

    T inline v(double x, double y) {
        return bridson_v(x,y);
    }

    T inline v(int x, int y) const {
        auto _x = CLAMP(x, 0, this->width-1);
        auto _y = CLAMP(y, 0, this->height-1);
        return this->_v[XY2I(_x,_y, this->width)];
    }

    T inline &v(int x, int y) {
        auto _x = CLAMP(x, 0, this->width-1);
        auto _y = CLAMP(y, 0, this->height-1);
        return this->_v[XY2I(_x,_y, this->width)];
    }

    T inline v(int i) const {
        auto _i = CLAMP(i, 0, this->size-1);
        return this->_v[_i];
    }

    T inline &v(int i) {
        auto _i = CLAMP(i, 0, this->size-1);
        return this->_v[_i];
    }

    /*
    void inline _v(const Eigen::Vector2i &p, const T &value) {
        this->_v(p.x(), p.y(), value);
    }

    void inline _v(const Eigen::Vector2f &p, const T &value) {
        return this->_v(std::nearbyint(p.x()), std::nearbyint(p.y()), value);
    }

    void inline _v(const Eigen::Vector2d &p, const T &value) {
        this->_v(p.x(), p.y(), value);
    }

    void inline _v(double x, double y, const T &value) {
        return this->_v(std::nearbyint(x), std::nearbyint(y), value);
    }

    void inline _v(int x, int y, const T &value) {
        auto _x = CLAMP(x, 0, this->width-1);
        auto _y = CLAMP(y, 0, this->height-1);
        this->v[XY2I(_x,_y,this->width)] = value;
    };

    void inline _v(int i, const T &value) {
        auto _i = CLAMP(i, 0, this->size-1);
        this->v[_i] = value;
    };
    */

    int inline getWidth() const { return width; }

    int inline getHeight() const { return height; }

    int inline getSize() const { return size; }

    friend inline std::ostream &operator<<(std::ostream &output, const Field2DV2<T> &field) {
        auto cols = field.width;
        auto rows = field.height;
        auto cells = ((field.element_size * field.width * field.height) / field.element_size);
        auto cell_size = field.element_size;
        auto field_size = field.element_size * field.width * field.height;

        output << "Field2d(cols:" << cols << ", rows:" << rows << ", cells:" << cells << ", cell_size:" << cell_size
               << "b, field_size:" << field_size << "b) = {"
               << std::endl;

        for (int x = 0; x < field.width; x++) {
            for (int y = 0; y < field.height; y++) {
                output << "\t";
                output << "(";
                field_traits::append_to_ostream<T>::value(output, field.u(x, y));
                output << ",";
                field_traits::append_to_ostream<T>::value(output, field.v(x, y));
                output << ")";
            }
            output << endl;
        }
        output << "}" << std::endl;
        return output;
    }

    friend inline std::ostream &operator<<(std::ostream &output, Field2DV2<T> *&field) {
        output << "// cell size in byte " << field->element_size << endl;
        output << "// field size in byte " << field->element_size * field->width * field->height << endl;
        output << "// number of cols " << field->height << endl;
        output << "// number of rows " << field->width << endl;
        output << "// number of cells " << (field->element_size * field->width * field->height) / field->element_size << endl;
        output << "F=[" << std::endl;
        for ( int x = 0; x < field->width; x++) {
            for ( int y = 0; y < field->height; y++) {
                output << "\t";
                output << "(";
                field_traits::append_to_ostream<T>::value(output, field->u(x, y));
                output << ",";
                field_traits::append_to_ostream<T>::value(output, field->v(x, y));
                output << ")";
            }
            output << endl;
        }
        output << "]" << std::endl;
        return output;
    }

    void inline copyTo(Field2DV2<T> *&b) {

        assert(b->element_size == this->element_size && "Size of data-type for quantity must be equal!");
        assert(b->width == this->width && "Width must be equal!");
        assert(b->height == this->height && "Height must be equal!");

        assert(b->element_size > 0 && this->element_size > 0 && "Size of data-type for quantity must not be 0!");
        assert(b->width > 0 && this->width > 0 && "Width must not be 0!");
        assert(b->height > 0 && this->height > 0 && "Height must not be 0!");

        memcpy(b->_u, this->_u, this->size * this->element_size);
        memcpy(b->_v, this->_v, this->size * this->element_size);

    }

    double cerp_u(double x, double y) const {
        const double offset_x = 0.0;
        const double offset_y = 0.0;
        x = min(max(x - offset_x, 0.0), width - 1.0);
        y = min(max(y - offset_y, 0.0), height - 1.0);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;

        int x0 = std::max(ix - 1, 0), x1 = ix, x2 = ix + 1, x3 = std::min(ix + 2, width - 1);
        int y0 = std::max(iy - 1, 0), y1 = iy, y2 = iy + 1, y3 = std::min(iy + 2, height - 1);

        double q0 = cerp(this->u(x0, y0), this->u(x1, y0), this->u(x2, y0), this->u(x3, y0), x);
        double q1 = cerp(this->u(x0, y1), this->u(x1, y1), this->u(x2, y1), this->u(x3, y1), x);
        double q2 = cerp(this->u(x0, y2), this->u(x1, y2), this->u(x2, y2), this->u(x3, y2), x);
        double q3 = cerp(this->u(x0, y3), this->u(x1, y3), this->u(x2, y3), this->u(x3, y3), x);

        return cerp(q0, q1, q2, q3, y);
    }

    double cerp_v(double x, double y) const {
        const double offset_x = 0.0;
        const double offset_y = 0.0;
        x = min(max(x - offset_x, 0.0), width - 1.0);
        y = min(max(y - offset_y, 0.0), height - 1.0);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;

        int x0 = std::max(ix - 1, 0), x1 = ix, x2 = ix + 1, x3 = std::min(ix + 2, width - 1);
        int y0 = std::max(iy - 1, 0), y1 = iy, y2 = iy + 1, y3 = std::min(iy + 2, height - 1);

        double q0 = cerp(this->v(x0, y0), this->v(x1, y0), this->v(x2, y0), this->v(x3, y0), x);
        double q1 = cerp(this->v(x0, y1), this->v(x1, y1), this->v(x2, y1), this->v(x3, y1), x);
        double q2 = cerp(this->v(x0, y2), this->v(x1, y2), this->v(x2, y2), this->v(x3, y2), x);
        double q3 = cerp(this->v(x0, y3), this->v(x1, y3), this->v(x2, y3), this->v(x3, y3), x);

        return cerp(q0, q1, q2, q3, y);
    }

    double bridson_u(double x, double y) const {
        x = min(max(x, 0.0), width - 1.001);
        y = min(max(y, 0.0), height - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;

        const int p0_offset = -1;
        const int p1_offset = 0;
        const int p2_offset = 1;
        const int p3_offset = 2;
        /* q_{j-1} */
        auto q_j_m1 = fluid::util::math::interpolation::w_m1(x) * this->u(ix + p0_offset, iy + p0_offset) +
                fluid::util::math::interpolation::w_0(x) * this->u(ix + p1_offset, iy + p0_offset) +
                fluid::util::math::interpolation::w_1(x) * this->u(ix + p2_offset, iy + p0_offset) +
                fluid::util::math::interpolation::w_2(x) * this->u(ix + p3_offset, iy + p0_offset);
        /* q_{j} */
        auto q_j = fluid::util::math::interpolation::w_m1(x) * this->u(ix + p0_offset, iy + p1_offset) +
                fluid::util::math::interpolation::w_0(x) * this->u(ix + p1_offset, iy + p1_offset) +
                fluid::util::math::interpolation::w_1(x) * this->u(ix + p2_offset, iy + p1_offset) +
                fluid::util::math::interpolation::w_2(x) * this->u(ix + p3_offset, iy + p1_offset);

        /* q_{j+1} */
        auto q_j_p1 = fluid::util::math::interpolation::w_m1(x) * this->u(ix + p0_offset, iy + p2_offset) +
                fluid::util::math::interpolation::w_0(x) * this->u(ix + p1_offset, iy + p2_offset) +
                fluid::util::math::interpolation::w_1(x) * this->u(ix + p2_offset, iy + p2_offset) +
                fluid::util::math::interpolation::w_2(x) * this->u(ix + p3_offset, iy + p2_offset);
        /* q_{j+2} */
        auto q_j_p2 = fluid::util::math::interpolation::w_m1(x) * this->u(ix + p0_offset, iy + p3_offset) +
                fluid::util::math::interpolation::w_0(x) * this->u(ix + p1_offset, iy + p3_offset) +
                fluid::util::math::interpolation::w_1(x) * this->u(ix + p2_offset, iy + p3_offset) +
                fluid::util::math::interpolation::w_2(x) * this->u(ix + p3_offset, iy + p3_offset);

        auto interp = fluid::util::math::interpolation::w_m1(y) * q_j_m1 +
                fluid::util::math::interpolation::w_0(y) * q_j +
                fluid::util::math::interpolation::w_1(y) * q_j_p1 +
                fluid::util::math::interpolation::w_2(y) * q_j_p2;

        return interp;
        // return interp < 0.000000001 ? 0.0 : interp;
    }


    double bridson_v(double x, double y) const {
        x = min(max(x, 0.0), width - 1.001);
        y = min(max(y, 0.0), height - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;

        const int p0_offset = -1;
        const int p1_offset = 0;
        const int p2_offset = 1;
        const int p3_offset = 2;
        /* q_{j-1} */
        auto q_j_m1 = fluid::util::math::interpolation::w_m1(x) * this->v(ix + p0_offset, iy + p0_offset) +
                      fluid::util::math::interpolation::w_0(x) * this->v(ix + p1_offset, iy + p0_offset) +
                      fluid::util::math::interpolation::w_1(x) * this->v(ix + p2_offset, iy + p0_offset) +
                      fluid::util::math::interpolation::w_2(x) * this->v(ix + p3_offset, iy + p0_offset);
        /* q_{j} */
        auto q_j = fluid::util::math::interpolation::w_m1(x) * this->v(ix + p0_offset, iy + p1_offset) +
                   fluid::util::math::interpolation::w_0(x) * this->v(ix + p1_offset, iy + p1_offset) +
                   fluid::util::math::interpolation::w_1(x) * this->v(ix + p2_offset, iy + p1_offset) +
                   fluid::util::math::interpolation::w_2(x) * this->v(ix + p3_offset, iy + p1_offset);

        /* q_{j+1} */
        auto q_j_p1 = fluid::util::math::interpolation::w_m1(x) * this->v(ix + p0_offset, iy + p2_offset) +
                      fluid::util::math::interpolation::w_0(x) * this->v(ix + p1_offset, iy + p2_offset) +
                      fluid::util::math::interpolation::w_1(x) * this->v(ix + p2_offset, iy + p2_offset) +
                      fluid::util::math::interpolation::w_2(x) * this->v(ix + p3_offset, iy + p2_offset);
        /* q_{j+2} */
        auto q_j_p2 = fluid::util::math::interpolation::w_m1(x) * this->v(ix + p0_offset, iy + p3_offset) +
                      fluid::util::math::interpolation::w_0(x) * this->v(ix + p1_offset, iy + p3_offset) +
                      fluid::util::math::interpolation::w_1(x) * this->v(ix + p2_offset, iy + p3_offset) +
                      fluid::util::math::interpolation::w_2(x) * this->v(ix + p3_offset, iy + p3_offset);

        auto interp = fluid::util::math::interpolation::w_m1(y) * q_j_m1 +
                      fluid::util::math::interpolation::w_0(y) * q_j +
                      fluid::util::math::interpolation::w_1(y) * q_j_p1 +
                      fluid::util::math::interpolation::w_2(y) * q_j_p2;

        return interp;
        //return interp < 0.000000001 ? 0.0 : interp;
    }

    void inline print() {
        std::cout << (*this) << std::endl;
    }

    void inline
    forEach(const std::function<void(T &/* current cell */, T &/* current cell */, int /* x coordinate */, int /* y coordinate */)> &f) {
// #pragma omp parallel for
        for (int x = 0; x < this->width; ++x) {
// #pragma omp parallel for
            for (int y = 0; y < this->height; ++y) {
                f(this->u(x, y), this->v(x, y), x, y);
            }
        }
    }

};

#endif //FLUID_SOLVER_RAW_FIELD2D_HPP

//
// Created by David Baldin on 2019-04-28.
//

#ifndef FLUID_SOLVER_RAW_MACGRID_HPP
#define FLUID_SOLVER_RAW_MACGRID_HPP

#include "field2d.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>

// Half index access like `u_{i+1/2,j}` is translated like `u[i+1,j]`

class MACGrid {
private:

    /* dimensions of velocity field */
    int n_x = 100, n_y = 100;

    /* cell pressure */
    Field2D<float> *pressure_field;

    /* velocity in u direction */
    Field2D<float> *u_velocity_field;

    /* velocity in v direction */
    Field2D<float> *v_velocity_field;

public:

    MACGrid(int n_x, int n_y) {
        this->n_x = n_x;
        this->n_y = n_y;

        this->pressure_field = new Field2D<float>(n_x, n_y);
        this->u_velocity_field = new Field2D<float>(n_x + 1, n_y);
        this->v_velocity_field = new Field2D<float>(n_x + 1, n_y);
    }

    virtual ~MACGrid() {
        delete this->pressure_field;
        delete this->u_velocity_field;
        delete this->v_velocity_field;
    }

    /**
     * Pressure value at (i,j)
     */
    float p_value(int i, int j) {
        return (*pressure_field)(i, j);
        // return p[i * n_x + j];
    }

    void p_value(int i, int j, float value) {
        (*pressure_field)(i, j) = value;
        //p[i * n_x + j] = value;
    }

    /**
     * Velocity in u direction at (i,j)
    */
    float u_value(int i, int j) {
        return (*u_velocity_field)(i, j);
        // return u[i * (n_x + 1) + j];
    }

    void u_value(int i, int j, float value) {
        (*u_velocity_field)(i, j) = value;
        // u[i * (n_x + 1) + j] = value;
    }

    /**
     * Velocity in v direction at (i,j)
    */
    float v_value(int i, int j) {
        return (*v_velocity_field)(i, j);
        // return v[i * n_x + j];
    }

    void v_value(int i, int j, float value) {
        (*v_velocity_field)(i, j) = value;
        // v[i * n_x + j] = value;
    }

    /**
     * Velocity vector at (i,j)
     *
     * Note: Half indices alà (i+0.5) are for mathematical reference and will be automatically casted to the  according integers
     */
    Eigen::Vector2f u_vec2(int i, int j) {
        return Eigen::Vector2f((u_value(i - 0.5, j) + u_value(i + 0.5, j)) * 0.5,
                               (v_value(i, j - 0.5) + v_value(i, j + 0.5)) * 0.5);
    }

    Eigen::Vector2f u_vec2(Eigen::Vector2i x) {
        return Eigen::Vector2f((u_value(x.x() - 0.5, x.y()) + u_value(x.x() + 0.5, x.y())) * 0.5,
                               (v_value(x.x(), x.y() - 0.5) + v_value(x.x(), x.y() + 0.5)) * 0.5);
    }

/* get u at i+0.5,j */
/**
 * Velocity vector at (i,j)
 * @param i
 * @param j
 * @return
 */
    Eigen::Vector2f u_vec2_i_plus_point_five(int i, int j) {
        return Eigen::Vector2f(u_value(i + 0.5, j),
                               (v_value(i, j - 0.5) + v_value(i, j + 0.5) +
                                v_value(i + 1, j - 0.5) + v_value(i + 1, j + 0.5)) * 0.25);
    }

/* get u at i,j+0.5 */
    Eigen::Vector2f u_vec2_j_plus_point_five(int i, int j) {
        return Eigen::Vector2f((u_value(i - 0.5, j) + u_value(i + 0.5, j) +
                                u_value(i - 0.5, j + 1) + u_value(i + 0.5, j + 1)) * 0.25,
                               v_value(i, j + 0.5));
    }

    void operator+=(Eigen::Vector2f v) {

    };


    friend std::ostream &operator<<(std::ostream &output, const MACGrid &v) {
        output << "// number of cols " << v.n_x << std::endl;
        output << "// number of rows " << v.n_y << std::endl;
        output << "// number of cells " << (v.n_x * v.n_y) << std::endl;
        output << "Grid={" << std::endl;
        output << "\tp : " << std::endl << (v.pressure_field) << std::endl;
        output << "\tu : " << std::endl << (v.u_velocity_field) << std::endl;
        output << "\tv : " << std::endl << (v.v_velocity_field) << std::endl;
        output << "}" << std::endl;
        return output;
    }

    void print() {
        std::cout << (*this) << std::endl;
    }

};


#endif //FLUID_SOLVER_RAW_MACGRID_HPP

//
// Created by David Baldin on 2019-06-11.
//

#ifndef FLUID_SOLVER_RAW_MACROS_HPP
#define FLUID_SOLVER_RAW_MACROS_HPP

/*
 * common macros
 */

/**
 * Helper to access 1D array with 2D coordinates
 */
#define COORDINATES_TO_INDEX_2D(X, Y, YS) ((X) * (YS)) + (Y)

/**
 * Helper to quickly wrap statements being compiled for debug purpose only
 */
#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

#endif //FLUID_SOLVER_RAW_MACROS_HPP

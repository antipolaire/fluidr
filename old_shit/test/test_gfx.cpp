/**
 * Just a simple copy'n'paste template class to quickly add tests
 */

#include <unistd.h>
#include <string>
#include <algorithm>
#include <functional>
#include <iostream>

#include "../src/gfx.h"
#include "../src/math.hpp"

/**
 * Check if file under "destination_filename" exists
 * @param destination_filename
 * @return true if file exists, otherwise false
 */
bool file_exists(const char *destination_filename) {
    return access(destination_filename, F_OK) != -1;
}

using namespace std;

using namespace cfd::util::gfx;
using namespace fluid::util::math;

void test_pgm_io() {

    unsigned int width;
    unsigned int height;

    const char *source_filename = "../../resource/assets/images/lena.pgm";
    const char *destination_filename = "../../resource/assets/images/lena2.pgm";

    uint8_t *pgm = cfd::util::gfx::loadPGM(source_filename, width, height);

    unsigned int i, j;
    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            pgm[((i) * (height)) + (j)] /= 2;
        }
    }

    cfd::util::gfx::saveAsPGM(destination_filename, pgm, width, height);

    free(pgm);

    assert(file_exists(destination_filename));

    std::remove(destination_filename);

    assert(!file_exists(destination_filename));
}


void test_pgm_interpolation() {

    unsigned int width;
    unsigned int height;

    const char *source_filename = "../../resource/assets/images/simple.pgm";
    const char *destination_filename = "../../resource/assets/images/simple_processed.pgm";

    if (file_exists(destination_filename)) {
        std::remove(destination_filename);
    }

    uint8_t *pgm = cfd::util::gfx::loadPGM(source_filename, width, height);

    std::function<double(double)> kernel = fluid::util::math::interpolation::general::kernel_bilinear;

    std::function<double(unsigned int, unsigned int)> f = [&pgm, &height](unsigned int x, unsigned int y) {
        unsigned int _x = std::clamp(x, (unsigned int) 0, (unsigned int) 513);
        unsigned int _y = std::clamp(y, (unsigned int) 0, (unsigned int) 513);
        uint8_t color = pgm[((_x) * (height)) + (_y)];
        // std::cout << "f(" << _x << "," << _y <<  ") = "<< (int)color << std::endl;
        return color / 255.0;
    };

    uint8_t *pgm_out = new uint8_t[width * height];

    unsigned int i, j;
    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {

            auto interpolated_color = fluid::util::math::interpolation::general::bi_apply_4th_order(
                    -1.0,
                    0.0,
                    1.0,
                    2.0,
                    i,
                    j,
                    f,
                    kernel,
                    0.0,
                    -0.5,
                    1.0);

            pgm_out[((i) * (height)) + (j)] = interpolated_color * 255;
        }
    }

    cfd::util::gfx::saveAsPGM(destination_filename, pgm_out, width, height);

    free(pgm);
    free(pgm_out);

    assert(file_exists(destination_filename));

    //std::remove(destination_filename);
    //assert(!file_exists(destination_filename));
}

int main(int argc, const char *argv[]) {
    // test_pgm_io();
    test_pgm_interpolation();
    return 0;
}
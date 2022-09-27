#ifndef GFX_UTILS_H
#define GFX_UTILS_H

#include <string>
#include <fstream>
#include <sstream>
#include <memory.h>
#include <assert.h>

using namespace std;

/**
 * Subprograms for graphics and bitmaps
 */
namespace cfd::util::gfx {

//---------------------------------------------------------------------------------------
/// @brief  Load a grayscale image from a PGM file
///
/// @param      strFileName     Name of the PGM file
/// @param[out] width           Width the image
/// @param[out] height          Height the image
/// @return                     A buffer containing the image, 0 in case of failure
///
/// @remark The application must release the buffer when done with it
//---------------------------------------------------------------------------------------
    uint8_t *loadPGM(const std::string &strFileName, unsigned int &width, unsigned int &height) {
        // Assertions
        assert(!strFileName.empty());

        // Constants
        const size_t BUFFER_SIZE = 50;

        // Declarations
        uint8_t *pBuffer = 0;
        char tmp[BUFFER_SIZE];
        std::ifstream stream;

        width = 0;
        height = 0;

        stream.open(strFileName.c_str(), std::ios::in | std::ios::binary);
        if (!stream.is_open())
            return 0;

        stream.getline(tmp, BUFFER_SIZE);
        if (strcmp(tmp, "P5") != 0) {
            stream.close();
            return 0;
        }

        do {
            stream.getline(tmp, BUFFER_SIZE);
        } while (tmp[0] == '#');

        std::istringstream str(tmp);
        str >> width >> height;

        if ((width <= 0) || (height <= 0)) {
            stream.close();
            return 0;
        }

        do {
            stream.getline(tmp, BUFFER_SIZE);
        } while (tmp[0] == '#');

        if (strcmp(tmp, "255") != 0) {
            stream.close();
            return 0;
        }

        pBuffer = new uint8_t[width * height];

        stream.read((char *) pBuffer, width * height * sizeof(uint8_t));

        stream.close();

        return pBuffer;
    }


//---------------------------------------------------------------------------------------
/// @brief  Save a grayscale image in a PGM file
///
/// @param  strFileName     Name of the PGM file
/// @param  pImage          Buffer containing the image
/// @param  width           Width the image
/// @param  height          Height the image
//---------------------------------------------------------------------------------------
    void saveAsPGM(const std::string &strFileName, unsigned char *pImage, unsigned int width,
                   unsigned int height) {
        // Assertions
        assert(!strFileName.empty());
        assert(pImage);
        assert(width > 0);
        assert(height > 0);

        // Declarations
        std::ofstream stream;

        stream.open(strFileName.c_str(), std::ios::out | std::ios::binary);
        if (stream.is_open()) {
            std::ostringstream str;
            str << "P5" << std::endl << width << " " << height << std::endl << 255 << std::endl;

            stream.write(str.str().c_str(), (std::streamsize) str.str().length());
            stream.write((char *) pImage, width * height * sizeof(uint8_t));

            stream.close();
        }
    }

};


#endif //GFX_UTILS_H

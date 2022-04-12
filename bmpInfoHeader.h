
#ifndef _BMP_INFO_HEADER_H_
#define _BMP_INFO_HEADER_H_

#include <cstdint> // for specific size integers
#include <fstream>

using namespace std;

class BmpInfoHeader
{
    public:
    void setImageDimensions(uint32_t mwidth, uint32_t mheight)
    {
        width = static_cast<int32_t>(mwidth);
        height = static_cast<int32_t>(mheight);
    }

    void save_on_file(std::ofstream& fout) {
        fout.write((char*)&this->sizeOfThisHeader, sizeof(uint32_t));
        fout.write((char*)&this->width, sizeof(int32_t));
        fout.write((char*)&this->height, sizeof(int32_t));
        fout.write((char*)&this->numberOfColorPlanes, sizeof(uint16_t));
        fout.write((char*)&this->colorDepth, sizeof(uint16_t));
        fout.write((char*)&this->compressionMethod, sizeof(uint32_t));
        fout.write((char*)&this->rawBitmapDataSize, sizeof(uint32_t));
        fout.write((char*)&this->horizontalResolution, sizeof(int32_t));
        fout.write((char*)&this->verticalResolution, sizeof(int32_t));
        fout.write((char*)&this->colorTableEntries, sizeof(uint32_t));
        fout.write((char*)&this->importantColors, sizeof(uint32_t));
    }

    public:
        uint32_t sizeOfThisHeader = 40;
        int32_t width = 0; // in pixels
        int32_t height = 0; // in pixels
        uint16_t numberOfColorPlanes = 1; // must be 1
        uint16_t colorDepth = 24;
        uint32_t compressionMethod = 0;
        uint32_t rawBitmapDataSize = 0; // generally ignored
        int32_t horizontalResolution = 0; // in pixel per meter
        int32_t verticalResolution = 0; // in pixel per meter
        uint32_t colorTableEntries = 0;
        uint32_t importantColors = 0;
};

#endif // _BMP_INFO_HEADER_H_
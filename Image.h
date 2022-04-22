#pragma once

#include <cstdint>
#include <vector>
#include <iostream>
#include <fstream>

#include "bmpHeader.h"
#include "bmpInfoHeader.h"

struct Pixel
{
    float blue;
    float green;
    float red;

    Pixel(){}
    Pixel(float b, float g, float r) : 
                        blue(b), 
                        green(g), 
                        red(r)
    {}
    ~Pixel(){}

    inline void save_on_file(std::ofstream& fout) const {
        
        unsigned char b = static_cast<unsigned char>(blue  * 255.f);
        unsigned char g = static_cast<unsigned char>(green * 255.f);
        unsigned char r = static_cast<unsigned char>(red   * 255.f);

        fout.write(reinterpret_cast<char*>(&b), sizeof(uint8_t));
        fout.write(reinterpret_cast<char*>(&g), sizeof(uint8_t));
        fout.write(reinterpret_cast<char*>(&r), sizeof(uint8_t));
    }
};

class Image
{
    public:
    Image(){}
    Image(uint32_t width, uint32_t height) : 
                                    m_width(width), 
                                    m_height(height),
                                    m_pixels(std::vector<Pixel>(width * height))
    {

        // Each row in the pixel data must be 4 byte aligned.
        // Calculated how much padding is needed per row
        paddingAmount = ((4 - (m_width * 3) % 4) % 4);
    }

    ~Image() {}

    Pixel GetPixel(uint32_t col, uint32_t row) const
    {
        return m_pixels[row * m_width + col];
    }

    void SetPixel(Pixel val, uint32_t col, uint32_t row)
    {
        m_pixels[row * m_width + col] = val;
    }

    void Export(const char* path) const
    {
        std::ofstream fout;
        fout.open(path, std::ios::out | std::ios::binary);

        if(!fout.is_open())
        {
            std::cout << "File could not be opened" << std::endl;
            return;
        }

        // Set dimensions on BMP header
        BmpHeader bmpheader;
        bmpheader.setImageDimensions(m_width, m_height, paddingAmount);
        BmpInfoHeader bmpinfoheader;
        bmpinfoheader.setImageDimensions(m_width, m_height);

        bmpheader.save_on_file(fout);
        bmpinfoheader.save_on_file(fout);

        for(uint32_t row = 0; row < m_height; row++)
        {
            for(uint32_t col = 0; col < m_width; col++)
            {
                Pixel p = GetPixel(col, row);
                p.save_on_file(fout);
            }

            uint8_t* padding_const = const_cast<uint8_t*>(padding); // remove const
            fout.write(reinterpret_cast<char*>(*padding_const), paddingAmount); //padding
        }
        fout.close();
    }

    private:
    uint32_t m_width;
    uint32_t m_height;
    uint8_t padding[3] = {0};
    uint32_t paddingAmount;
    std::vector<Pixel> m_pixels;
};
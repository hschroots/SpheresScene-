
#ifndef _BMP_HEADER_H_
#define _BMP_HEADER_H_

#include <cstdint> // for specific size integers
#include <fstream>

using namespace std;

class BmpHeader
{
	public:

	void setImageDimensions(uint32_t width, uint32_t height, uint32_t paddingAmount)
	{
		sizeOfBitmapFile = 54 + (width * height * 3) + (paddingAmount * height);
	}

	void save_on_file(std::ofstream& fout) {
		fout.write(this->bitmapSignatureBytes, 2);
		fout.write((char*)&this->sizeOfBitmapFile, sizeof(uint32_t));
		fout.write((char*)&this->reservedBytes, sizeof(uint32_t));
		fout.write((char*)&this->pixelDataOffset, sizeof(uint32_t));
	}

private:
	char bitmapSignatureBytes[2] = {'B', 'M'};
	uint32_t sizeOfBitmapFile = 54; // total size of bitmap file
	uint32_t reservedBytes = 0;
	uint32_t pixelDataOffset = 54;
};

#endif //_BMP_HEADER_H_

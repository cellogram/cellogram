////////////////////////////////////////////////////////////////////////////////
#include "tif_reader.h"

#include <tinytiffreader.h>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	bool read_tif_image(const std::string &path, Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> &img)
	{
		bool ok = false;
		TinyTIFFReaderFile* tiffr=NULL;
		tiffr=TinyTIFFReader_open(path.c_str());
		if (!tiffr) {
			std::cerr<<"    ERROR reading (not existent, not accessible or no TIFF file)\n";
		} else {
			uint32_t width=TinyTIFFReader_getWidth(tiffr);
			uint32_t height=TinyTIFFReader_getHeight(tiffr);
			uint16_t* image=(uint16_t*)calloc(width*height, sizeof(uint16_t));
			TinyTIFFReader_getSampleData(tiffr, image, 0);

			//image in (ROW-MAJOR!)
			img = Eigen::Map<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>>(image, height, width).transpose();

			free(image);

			ok = true;
		}
		TinyTIFFReader_close(tiffr);
		return ok;
	}

} // namespace cellogram

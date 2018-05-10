////////////////////////////////////////////////////////////////////////////////
#include "image_reader.h"

#include "StringUtils.h"
#include <tinytiffreader.h>
#include <igl/png/readPNG.h>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	namespace
	{
		template<typename T>
		bool read_mem_to_eigen(TinyTIFFReaderFile* tiffr, Eigen::MatrixXd &img)
		{
			uint32_t width = TinyTIFFReader_getWidth(tiffr);
			uint32_t height = TinyTIFFReader_getHeight(tiffr);

			T* image = (T*)calloc(width*height, sizeof(T));
			bool ok = TinyTIFFReader_getSampleData(tiffr, image, 0);


			//image in (ROW-MAJOR!)
			if(ok)
			{
				Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp;
				tmp = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(image, width, height);

				//normalize img
				double min = tmp.minCoeff();
				double max = tmp.maxCoeff();

				img = ((tmp.template cast<double>().array() - min) / (max - min)).rowwise().reverse();
			}

			free(image);

			return ok;
		}
	}

	bool read_tif_image(const std::string &path, Eigen::MatrixXd &img)
	{
		bool ok = false;
		TinyTIFFReaderFile* tiffr=NULL;
		tiffr=TinyTIFFReader_open(path.c_str());
		if (!tiffr) {
			std::cerr<<"    ERROR reading (not existent, not accessible or no TIFF file)"<< std::endl;
		} else {
			uint16_t samples = TinyTIFFReader_getBitsPerSample(tiffr);


			ok = true;

			if (samples == 8)
			{
				ok = read_mem_to_eigen<uint8_t>(tiffr, img);
			}
			else if (samples == 16)
			{
				ok = read_mem_to_eigen<uint16_t>(tiffr, img);
			}
			else if (samples == 32)
			{
				ok = read_mem_to_eigen<uint32_t>(tiffr, img);
			}
			else
			{
				std::cerr << "    ERROR wrong format" << std::endl;
				ok = false;
			}


		}
		TinyTIFFReader_close(tiffr);
		return ok;
	}

	bool read_png_image(const std::string &path, Eigen::MatrixXd &img)
	{
		Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R, G, B, A; // Image
		bool ok = igl::png::readPNG(path, R, G, B, A);

		if(ok)
		{
		//normalize img
			double min = R.minCoeff();
			double max = R.maxCoeff();

			img = (R.cast<double>().array() - min) / (max - min);
		}

		return ok;
	}

	bool read_image(const std::string & path, Eigen::MatrixXd & img)
	{
		if (StringUtils::endswidth(path, ".png"))
			return read_png_image(path, img);
		if (StringUtils::endswidth(path, ".tif") || StringUtils::endswidth(path, ".tiff"))
			return read_tif_image(path, img);

		std::cout << "Unsupported image format " << path << std::endl;
		return false;
	}

} // namespace cellogram

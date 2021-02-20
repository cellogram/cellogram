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
		template <typename T>
		bool read_mem_to_eigen(TinyTIFFReaderFile *tiffr, T *image, Eigen::MatrixXd &img) {
			uint32_t width = TinyTIFFReader_getWidth(tiffr);
			uint32_t height = TinyTIFFReader_getHeight(tiffr);

			if(!image)
				image = (T *)calloc(width * height, sizeof(T));
			bool ok = TinyTIFFReader_getSampleData(tiffr, image, 0);

			//image in (ROW-MAJOR!)
			if (ok) {
				Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp;
				tmp = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >(image, width, height);

				// normalize img
				// double min = tmp.minCoeff();
				// double max = tmp.maxCoeff();

				// img = ((tmp.template cast<double>().array() - min) / (max - min)).rowwise().reverse();
				tmp.transposeInPlace();
				img = tmp.template cast<double>();
			}

			return ok;
		}
	}


	bool read_tif_image(const std::string &path, std::vector<Eigen::MatrixXd> & img3D)
	{

		bool ok = false;
		TinyTIFFReaderFile *tiffr = NULL;
		tiffr = TinyTIFFReader_open(path.c_str());

		if (!tiffr) {

			std::cerr << "ERROR reading (not existent, not accessible or no TIFF file)" << std::endl;
		} else {

			// logger().info(path);
			// logger().info(TinyTIFFReader_getImageDescription(tiffr));

			img3D.clear();
			int currentImgCount = 0;
			uint32_t countSlice = 0;

			Eigen::MatrixXd sliceMat;

			uint8_t *buffer8 = nullptr;
			uint16_t *buffer16 = nullptr;
			uint32_t *buffer32 = nullptr;

			do {
				// prepare to read this slice
				ok = true;
				++countSlice;

				// read data to "sliceMat"
				uint16_t samples = TinyTIFFReader_getBitsPerSample(tiffr);
				if (samples == 8) {
					ok = read_mem_to_eigen<uint8_t>(tiffr, buffer8, sliceMat);
					// scale to 0~1 for visualization
					sliceMat /= double((1 << 8) - 1);
				}
				else if (samples == 16) {
					ok = read_mem_to_eigen<uint16_t>(tiffr, buffer16, sliceMat);
					// scale to 0~1 for visualization
					sliceMat /= double((1 << 16) - 1);
				}
				else if (samples == 32) {
					ok = read_mem_to_eigen<uint32_t>(tiffr, buffer32, sliceMat);
					// scale to 0~1 for visualization
					sliceMat /= double((1ll << 32) - 1);
				} 
				else {
					std::cerr << "ERROR: TinyTIFFReader_getBitsPerSample wrong format" << std::endl;
					ok = false;
					break;
				}

				// push a slice to "img"
				img3D.push_back(sliceMat);
			} while (TinyTIFFReader_readNext(tiffr));

			free(buffer8);
			free(buffer16);
			free(buffer32);
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


	bool read_image(const std::string & path, std::vector<Eigen::MatrixXd> & img3D)
	{
		// Only read TIFF
		// if (StringUtils::endswidth(path, ".png"))
		//  	return read_png_image(path, img);
		if (StringUtils::endswidth(path, ".tif") || StringUtils::endswidth(path, ".tiff"))
			return read_tif_image(path, img3D);

		std::cout << "Unsupported image format " << path << std::endl;
		return false;
	}

} // namespace cellogram

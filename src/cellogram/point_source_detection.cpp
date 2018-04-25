////////////////////////////////////////////////////////////////////////////////
#include "convex_hull.h"
#include "delaunay.h"
#include "navigation.h"
#include "fitGaussian2D.h"
#include <igl/edges.h>
#include <igl/boundary_loop.h>
#include <igl/triangle/cdt.h>
#include <algorithm>
#include <numeric>
#include <stack>
#include <geogram/basic/geometry.h>
#undef IGL_STATIC_LIBRARY
#include <igl/edge_lengths.h>
#include <igl/colon.h>
#include "tcdf/tcdf.h"
#include <queue>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------
	typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

namespace {

	//put here all helpers functions
	void padarray_symmetric(const Eigen::MatrixXd &img, int padSize, Eigen::MatrixXd &imgXT)
	{
		int ydim = img.rows();
		int xdim = img.cols();

		imgXT.resize(ydim + 2 * padSize, xdim + 2 * padSize);
		imgXT.block(padSize, padSize, ydim, xdim) = img;

		for (int i = 1; i <= padSize; i++)
		{
			imgXT.row(padSize - i) = imgXT.row(padSize + i);
			imgXT.row(ydim + padSize + i - 1) = imgXT.row(ydim + padSize - i - 1);

			imgXT.col(padSize - i) = imgXT.col(padSize + i);
			imgXT.col(xdim + padSize + i - 1) = imgXT.col(xdim + padSize - i - 1);
		}

	}

	Eigen::MatrixXd conv2(const Eigen::VectorXd &k1, const Eigen::VectorXd &k2,const Eigen::MatrixXd &imgXT)
	{
		int padSizeY = (k1.size() - 1) / 2;
		int padSizeX = (k2.size() - 1) / 2;
		int ydim = imgXT.rows() - 2* padSizeY;
		int xdim = imgXT.cols() - 2* padSizeX;

		Eigen::MatrixXd img_filtered(ydim, xdim + 2* padSizeX);

		//cols
		for (int j = 0; j < ydim; j++)
		{
			for (int i = 0; i < xdim + 2*padSizeX; i++)
			{
				double c_value = 0;
				for (int k = 0; k < k1.rows(); k++)
				{
					c_value += imgXT(j + k, i)*k1(k);
				}
				img_filtered(j, i) = c_value;
			}
		}

		// rows
		Eigen::MatrixXd tmp = imgXT;
		tmp.block(padSizeY, 0, ydim, xdim + 2* padSizeX) = img_filtered;
		img_filtered.resize(ydim, xdim);
		//std::cout << tmp << std::endl;
		for (int i = 0; i < xdim; i++)
		{
			for (int j = 0; j < ydim; j++)
			{
				double c_value = 0;
				for (int k = 0; k < k2.rows(); k++)
				{
					double asd = tmp(j + padSizeY, i + k);
					c_value += tmp(j + padSizeY, i + k)*k2(k);
				}
				img_filtered(j, i) = c_value;
			}
		}

		return img_filtered;
	}
	
	Eigen::MatrixXd ordfilt2(Eigen::MatrixXd &img, int order, Eigen::MatrixXi mask)
	{
		// analogous to matlabs ordfilt2

		// pad image
		int padSize = (mask.cols() - 1) / 2;
		int ydim = img.rows();
		int xdim = img.cols();
		Eigen::MatrixXd filtered_img(ydim,xdim);

		Eigen::MatrixXd imgPad;
		imgPad.resize(ydim + 2 * padSize, xdim + 2 * padSize);
		imgPad.setZero();
		imgPad.block(padSize, padSize, ydim, xdim) = img;

		// scan through image and replace each pixel with the highest value of it's domain
		for (int i = padSize; i < ydim + padSize; i++)
		{
			for (int j = padSize; j < xdim + padSize; j++)
			{
				// Loop through domain and find highest value
				std::vector<double> data; data.reserve(mask.size());
				for (int k = -padSize; k <= padSize; k++)
				{
					for (int l = -padSize; l <= padSize; l++)
					{
							data.push_back(imgPad(i + k, j + l));
					}
				}
				assert(data.size() == mask.size());
				std::sort(data.begin(), data.end());
				filtered_img(i - padSize, j - padSize) = data[order];
			}
		}
		return filtered_img;
	}


	Eigen::MatrixXd locmax2d(Eigen::MatrixXd &img, int maskSize)
	{
		int rows = maskSize;
		int padSize = (maskSize - 1) / 2;
		int cols = maskSize;
		Eigen::MatrixXi mask;
		mask.setZero(maskSize,maskSize);
		int numEl = rows * cols;

		Eigen::MatrixXd fImg = ordfilt2(img, numEl-1, mask);
		Eigen::MatrixXd fImg2 = ordfilt2(img, numEl - 2, mask);

		// take only those positions where the max filter and the original image value
		// are equal -> this is a local maximum
		for (int i = 0; i < img.rows(); i++)
		{
			for (int j = 0; j < img.cols(); j++)
			{
				if (fImg2(i, j) == fImg(i, j))
					fImg(i, j) = 0;
				if (fImg(i, j) != img(i, j))
					fImg(i, j) = 0;
			}
		}

		// set image border to zero
		Eigen::MatrixXd fImgFinal(img.rows(),img.cols());
		fImgFinal.setZero(img.rows(), img.cols());

		fImgFinal.block(padSize, padSize, img.rows() - 2 * padSize, img.cols() - 2 * padSize) = fImg.block(padSize, padSize, img.rows() - 2 * padSize, img.cols() - 2 * padSize);

		return fImgFinal;
	}

	Eigen::MatrixXi bwlabel(const MatrixXb &mask)
	{

		typedef Eigen::Vector2d vec2;
		MatrixXb tmp = mask;
		Eigen::MatrixXi labels(mask.rows(), mask.cols());
		labels.setZero();


		Eigen::MatrixXi visited(mask.rows(), mask.cols());
		int label = 1;
		// loop through entire mask
		for (int i = 0; i < mask.rows(); i++)
		{
			for (int j = 0; j < mask.cols(); j++)
			{
				if (tmp(i, j))
				{
					std::queue<vec2> q;
					q.emplace(i, j);
					while (!q.empty())
					{
						auto index = q.front();
						q.pop();
						if (!tmp(index(0), index(1)))
							continue;
						labels(index(0), index(1)) = label;
						tmp(index(0), index(1)) = false;

						if (index(0) > 0)
							q.emplace(index(0) - 1, index(1));
						if (index(1) > 0)
							q.emplace(index(0), index(1) - 1);
						if (index(0) < mask.rows() - 1)
							q.emplace(index(0) + 1, index(1));
						if (index(1) < mask.cols() - 1)
							q.emplace(index(0), index(1) + 1);
					}
					label++;
				}
			}
		}
		return labels;
	}

	void fitGaussians2D(const Eigen::MatrixXd &img, const Eigen::MatrixXi &xy, const Eigen::VectorXd &A, const Eigen::VectorXd &sigma,
		const Eigen::VectorXd &c, const MatrixXb &mask, Eigen::MatrixXd &V, Eigen::MatrixXd &V_std, Eigen::VectorXd &pval_Ar)
	{
		int np = xy.rows();
		std::string mode = "xy";
		Eigen::MatrixXi labels = bwlabel(mask);
		//std::cout << "labels" << labels << std::endl;
		int ny = img.rows();
		int nx = img.cols();
		double kLevel = 1.959963984540054;

		Eigen::Vector2d iRange(img.minCoeff(), img.maxCoeff());

		Eigen::Vector2i estIdx(1, 2);
		double sigma_max = sigma.maxCoeff();
		int w2 = ceil(2 * sigma_max);
		int w4 = ceil(4 * sigma_max);

		Eigen::VectorXd g(2 * w4 + 1);
		int k = 0;
		for (int i = -w4; i <= w4; i++)
		{
			g(k) = std::exp(-i * i / (2 * sigma_max *sigma_max));
			k++;
		}
		Eigen::MatrixXd g2 = g*g.transpose();
		Eigen::Map<Eigen::RowVectorXd> gv(g2.data(), g2.size());

		for (int p = 0; p < np; p++)
		{
			// ignore points in border
			if (xy(p, 0) < w4 || xy(p, 0) >= nx - w4 || xy(p, 1) < w4 || xy(p, 1) >= ny - w4)
				continue;

			// label mask
			Eigen::MatrixXi maskWindow;
			Eigen::MatrixXd window;
			maskWindow = labels.block(xy(p, 1) - w4, xy(p, 0) - w4, 2 * w4 + 1, 2 * w4 + 1);
			
			maskWindow = (maskWindow.array() == maskWindow(w4 + 1, w4 + 1)).select(0, maskWindow);
			//maskWindow(maskWindow == maskWindow(w4 + 1, w4 + 1)) = 0;
			window = img.block(xy(p, 1) - w4, xy(p, 0) - w4, 2 * w4 + 1, 2 * w4 + 1);

			double c_init = c(p);

			std::cout << "\n\n\n-------------\n" << window << "\n\n\n" << std::endl;

			// set any other components to NaN
			int npx = 0;
			for (int i = 0; i < window.rows(); i++)
			{
				for (int j = 0; j < window.cols(); j++)
				{
					if (maskWindow(i,j) != 0)
						window(i, j) = std::nan("");
					else
						npx++;
				}
			}

			if (npx < 10) // only perform fit if window contains sufficient data points
				continue;

			double A_init = A(p);
			fitGaussian2D(window, 0, 0, A_init, sigma(p), c_init);
		}

			

	}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////



void point_source_detection(const Eigen::MatrixXd &img, const double sigma, Eigen::MatrixXd &V, Eigen::MatrixXd &V_std, Eigen::VectorXd &pval_Ar)
{
	assert(img.minCoeff()>=0);
	assert(img.maxCoeff()<=1);

	// Gaussian kernel
	const int w = std::ceil(4 * sigma);
	const int kernel_size = 2 * w + 1;
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(kernel_size, -w, w);
	Eigen::VectorXd u;
	u.setOnes(kernel_size);
	Eigen::VectorXd g(kernel_size);
	for (int i = 0; i < kernel_size; i++)
	{
		g(i) = std::exp( (-1) * x(i) * x(i) / (2 * sigma * sigma));
	}
	
	// Convolutions
	Eigen::MatrixXd imgXT;
	padarray_symmetric(img, w, imgXT);

	Eigen::MatrixXd fg = conv2(g, g, imgXT);
	Eigen::MatrixXd fu = conv2(u, u, imgXT);
	Eigen::MatrixXd fu2 = conv2(u, u, imgXT.array().square());
	
	// Laplacian of Gaussian
	Eigen::MatrixXd gx2 = g.array()*x.array().square();
	double sigma2 = sigma * sigma;
	Eigen::MatrixXd imgLoG = 2 * fg / sigma2 - (conv2(g, gx2, imgXT) + conv2(gx2, g, imgXT)) / (sigma2*sigma2);
	imgLoG = imgLoG / (2 * 3.1415 *sigma2);

	
	
	// 2 - D kernel
	Eigen::MatrixXd g2 = g*g.transpose();
	double n = g2.size();
	double gsum = g2.sum();
	double g2sum = g2.array().square().sum();

	

	// solution to linear system
	Eigen::MatrixXd A_est = (fg - gsum * fu / n) / (g2sum - gsum*gsum / n);
	Eigen::MatrixXd c_est = (fu - A_est * gsum) / n;

	

	// Prefilter
	Eigen::Map<Eigen::VectorXd> g2_vector(g2.data(), g2.size());
	Eigen::VectorXd ones = Eigen::VectorXd::Ones(g2.size());
	Eigen::MatrixXd J(g2.size(), 2);
	J.col(0) = g2_vector;
	J.col(1) = ones;
		
	Eigen::MatrixXd C;
	C = J.transpose()*J;
	C = C.inverse();

	Eigen::MatrixXd	f_c = fu2 - 2 * c_est.cwiseProduct(fu) + n * c_est.array().square().matrix();
	Eigen::MatrixXd	RSS;
	{
		Eigen::MatrixXd tmp = fg - c_est * gsum;
		RSS = A_est.array().square().matrix() * g2sum - 2 * A_est.cwiseProduct(tmp) + f_c;
	}
	RSS = (RSS.array().abs() < 1e-10).select(0, RSS); // negative numbers may result from machine epsilon / roundoff precision



	Eigen::MatrixXd	sigma_e2 = RSS / (n - 3);
	Eigen::MatrixXd sigma_A = sigma_e2*C(0, 0);
	sigma_A = sigma_A.cwiseSqrt();


	// standard deviation of residuals
	Eigen::MatrixXd	sigma_res = RSS / (n - 1);
	sigma_res = sigma_res.cwiseSqrt();

	double kLevel = 1.959963984540054;

	Eigen::MatrixXd	SE_sigma_c = sigma_res / sqrt(2 * (n - 1)) * kLevel; //checked
	Eigen::MatrixXd	df2;
	Eigen::MatrixXd scomb;
	{
		//df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
		Eigen::MatrixXd tmp = sigma_A.array().square() + SE_sigma_c.array().square();
		Eigen::MatrixXd tmp2 = tmp.array().square();
		Eigen::MatrixXd tmp3 = sigma_A.array().square().square() + SE_sigma_c.array().square().square(); //checked
		df2 = (n - 1) * tmp2.array() / tmp3.array(); //checked
		scomb = tmp/n;
	}
		
	scomb = scomb.cwiseSqrt(); //checked
	Eigen::MatrixXd	 T = A_est - sigma_res * kLevel;
	T = T.cwiseQuotient(scomb); //checked
	Eigen::MatrixXd pval;
	matlabautogen::tcdf(-T, df2, pval);
	
	//mask of admissible positions for local maxima
	MatrixXb mask = pval.array() < 0.05;

	// all local max
	Eigen::MatrixXd allMax = locmax2d(imgLoG, 2 * std::ceil(sigma) + 1); // checked

	// local maxima above threshold in image domain
	Eigen::MatrixXd imgLM = allMax.array() * mask.cast<double>().array(); // checked

	if (imgLM.cwiseAbs().sum() - 1e-8 < 0)
		return;

	//->set threshold in LoG domain
	Eigen::MatrixXd tmp;
	tmp = (imgLM.array().abs() < 1e-10).select(std::numeric_limits<double>::max(), imgLoG);
	double logThreshold = tmp.minCoeff();
	MatrixXb logMask = (imgLoG.array() >= logThreshold);
	
	// combine masks
	mask = mask.array() || logMask.array();
	std::cout << imgLM << std::endl;

	// re - select local maxima
	imgLM = allMax.array()*mask.cast<double>().array();

	Eigen::MatrixXi lm(imgLM.size(),2);
	Eigen::VectorXd A_est_idx(imgLM.size()), c_est_idx(imgLM.size()), s_est_idx;
	int k = 0;
	for (int i = 0; i < imgLM.cols(); i++)
	{
		for (int j = 0; j < imgLM.rows(); j++)
		{
			if (std::abs(imgLM(j, i)) > 1e-8)
			{
				lm(k,0) = i;
				lm(k,1) = j;
				A_est_idx(k) = A_est(j, i);
				c_est_idx(k) = c_est(j, i);
				k++;
			}
		}
	}

	lm.conservativeResize(k, 2);
	A_est_idx.conservativeResize(k);
	c_est_idx.conservativeResize(k);
	s_est_idx = Eigen::VectorXd::Constant(k, sigma);

	std::cout << "lm:\n" << lm << std::endl;
	std::cout << "A_est_idx:\n" << A_est_idx << std::endl;
	std::cout << "c_est_idx:\n" << c_est_idx << std::endl;
	std::cout << "mask:\n" << mask << std::endl;

	fitGaussians2D(img, lm, A_est_idx, s_est_idx, c_est_idx, mask, V, V_std, pval_Ar); //inputs checked

	//% remove NaN values
	//	idx = ~isnan([pstruct.x]);
}

} // namespace cellogram

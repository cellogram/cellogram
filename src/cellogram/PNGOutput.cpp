#include "PNGOutput.h"

#include <cellogram/State.h>
#include <cellogram/laplace_energy.h>

namespace cellogram
{
	PNGOutput::PNGOutput(const double scale)
	: scale_(scale)
	{
	}

	void PNGOutput::save()
	{
#ifdef CELLOGRAM_WITH_PNG
		writer_.close();
#endif
	}

	void PNGOutput::init(const std::string &name, const int w, const int h)
	{
#ifdef CELLOGRAM_WITH_PNG
		w_ = w;
		h_ = h;
		writer_.clear();
		writer_.pngwriter_rename(name.c_str());
		writer_.resize(w*scale_, h*scale_ + text_offset);

		writer_.filledsquare(0, 0, w*scale_, h*scale_ + text_offset, 1., 1., 1.);
#endif
	}

	void PNGOutput::draw_image()
	{
#ifdef CELLOGRAM_WITH_PNG
		const auto &state = cellogram::State::state();
		const auto &img = state.img;

		for(long i = 0; i < img.rows(); ++i)
		{
			for(long j = 0; j < img.cols(); ++j)
			{
				writer_.filledsquare(i*scale_, j*scale_, i*scale_ + scale_, j*scale_ + scale_, img(i,j), img(i,j), img(i,j));
			}
		}
#endif
	}


	void PNGOutput::draw_detection()
	{
#ifdef CELLOGRAM_WITH_PNG
		const int circle_radius = 500./std::max(w_, h_)*scale_;
		static const double red = 1;
		static const double green = 1;
		static const double blue = 0;
		static const double opacity = 0.6;


		const auto &state = cellogram::State::state();
		const auto &pts = state.mesh.detected;

		for(int i = 0; i < pts.rows(); ++i)
		{
			writer_.filledcircle_blend(pts(i, 0)*scale_, (h_ - pts(i, 1))*scale_, circle_radius, opacity, red, green, blue);
		}
#endif
	}


	void PNGOutput::draw_untangle()
	{
#ifdef CELLOGRAM_WITH_PNG
		const int circle_radius = 500./std::max(w_, h_)*scale_;
		static const double tri_red = 1;
		static const double tri_green = 1;
		static const double tri_blue = 1;
		static const double bad_region_opacity = 0.6;

		static const std::string font_dir = DATA_DIR;
		static const std::string font = font_dir + "font.ttf";
		static const int font_size = 20;

		const auto &state = cellogram::State::state();
		const auto &mesh = state.mesh;
		const auto &pts = mesh.detected;
		const auto &tris = mesh.triangles;

		for(int i = 0; i < tris.rows(); ++i)
		{
			const int p1 = tris(i, 0);
			const int p2 = tris(i, 1);
			const int p3 = tris(i, 2);

			writer_.triangle(pts(p1, 0)*scale_, (h_ - pts(p1, 1))*scale_, pts(p2, 0)*scale_, (h_ - pts(p2, 1))*scale_, pts(p3, 0)*scale_, (h_ - pts(p3, 1))*scale_, tri_red, tri_green, tri_blue);
		}


		const auto &added = mesh.added_by_untangler;
		const auto &deleted = mesh.deleted_by_untangler;

		for(int i = 0; i < added.size(); ++i)
		{
			writer_.filledcircle(pts(added[i], 0)*scale_, (h_ - pts(added[i], 1))*scale_, circle_radius, 1., 0., 0.);
		}

		for(int i = 0; i < deleted.rows(); ++i)
		{
			writer_.filledcircle(deleted(i, 0)*scale_, (h_ - deleted(i, 1))*scale_, circle_radius, 1., 0., 0.);
		}



		const auto &regions = state.regions;
		int n_bad = 0;

		for(const auto &r : regions)
		{
			if(r.status == Region::OK || r.status != Region::NOT_CHECKED)
				continue;

			n_bad ++;

			for(long i = 0; i < r.region_triangles.size(); ++i)
			{
				const int p1 = tris(r.region_triangles(i), 0);
				const int p2 = tris(r.region_triangles(i), 1);
				const int p3 = tris(r.region_triangles(i), 2);

				writer_.filledtriangle_blend(
					pts(p1, 0)*scale_, (h_ - pts(p1, 1))*scale_,
					pts(p2, 0)*scale_, (h_ - pts(p2, 1))*scale_,
					pts(p3, 0)*scale_, (h_ - pts(p3, 1))*scale_,
					bad_region_opacity, 1., 0., 0.);
			}
		}


		Eigen::VectorXd energy;
		laplace_energy(pts, tris, energy);


		const std::string added_msg = "Added vertices: " + std::to_string(added.size());
		const std::string deleted_msg = "Deleted vertices: " + std::to_string(deleted.rows());
		const std::string bad_msg = "Bad regions: " + std::to_string(n_bad);

		const std::string max_msg = "Max energy: " + std::to_string(energy.maxCoeff());
		const std::string min_msg = "Min energy: " + std::to_string(energy.minCoeff());

		std::vector<char> cadded_msg(added_msg.c_str(), added_msg.c_str() + added_msg.size() + 1);
		std::vector<char> cdeleted_msg(deleted_msg.c_str(), deleted_msg.c_str() + deleted_msg.size() + 1);

		std::vector<char> cbad_msg(bad_msg.c_str(), bad_msg.c_str() + bad_msg.size() + 1);
		std::vector<char> cfont(font.c_str(), font.c_str() + font.size() + 1);

		std::vector<char> cmax_msg(max_msg.c_str(), max_msg.c_str() + max_msg.size() + 1);
		std::vector<char> cmin_msg(min_msg.c_str(), min_msg.c_str() + min_msg.size() + 1);

		writer_.plot_text(&cfont[0], font_size, 0, h_*scale_ + 5, 0, &cadded_msg[0],	added.size() == 0 ? 0. : 1., 0., 0.);
		writer_.plot_text(&cfont[0], font_size, 0, h_*scale_ + 35, 0, &cdeleted_msg[0], deleted.rows() == 0 ? 0. : 1., 0., 0.);

		writer_.plot_text(&cfont[0], font_size, 0, h_*scale_ + 75, 0, &cbad_msg[0], n_bad == 0 ? 0. : 1., 0., 0.);

		writer_.plot_text(&cfont[0], font_size, 0, h_*scale_ + 115, 0, &cmax_msg[0], 0., 0., 0.);
		writer_.plot_text(&cfont[0], font_size, 0, h_*scale_ + 145, 0, &cmin_msg[0], 0., 0., 0.);
#endif
	}
}
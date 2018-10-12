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
#else
		writer_ << "</svg>\n";
		std::ofstream out(name_.c_str());
		out << writer_.str() << std::endl;
		out.close();
#endif
	}

	void PNGOutput::init(const std::string &name, const int w, const int h)
	{
		w_ = w;
		h_ = h;
		name_ = name;
#ifdef CELLOGRAM_WITH_PNG
		writer_.clear();
		writer_.pngwriter_rename(name.c_str());
		writer_.resize(w*scale_, h*scale_ + text_offset);
		writer_.filledsquare(0, 0, w*scale_, h*scale_ + text_offset, 1., 1., 1.);
#else
		writer_ <<"<?xml version=\"1.0\" standalone=\"no\"?>\n";
		writer_ <<"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
		writer_ <<"<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" ";
		writer_ << "width=\""<< w*scale_ <<"\" height=\""<< h*scale_ + text_offset <<"\" ";
		writer_ << "x=\""<< 0 <<"\" y=\""<< 0 <<"\">\n";
#endif
	}

	void PNGOutput::draw_image()
	{
		const auto &state = cellogram::State::state();
		const auto &img = state.img;

#ifndef CELLOGRAM_WITH_PNG
		writer_ << "<rect x=\""<<0<<"\" y=\""<<0<<"\" ";
		writer_ << "width=\""<<img.rows()*scale_<<"\" height=\""<<img.cols()*scale_<<"\" ";
		writer_ << "fill=\"rgb("<<0<<","<<0<<","<<0<<")\" stroke=\"none\"";
		writer_ << "/>\n";
#endif

		for(long i = 0; i < img.rows(); ++i)
		{
			for(long j = 0; j < img.cols(); ++j)
			{
				const int x = i*scale_;
				const int y = j*scale_;
				const double c = img(i,j);
#ifdef CELLOGRAM_WITH_PNG
				const int w = i*scale_ + scale_;
				const int h = j*scale_ + scale_;

				writer_.filledsquare(x, y, w, h, c, c, c);
#else
				const int ci = c*255;
				if(ci < 50)
					continue;
				writer_ << "<rect x=\""<<x<<"\" y=\""<<y<<"\" ";
				writer_ << "width=\""<<scale_<<"\" height=\""<<scale_<<"\" ";
				writer_ << "fill=\"rgb("<<ci<<","<<ci<<","<<ci<<")\" stroke=\"none\"";
				writer_ << "/>\n";
#endif
			}
		}

	}


	void PNGOutput::draw_detection()
	{

		const int circle_radius = 500./std::max(w_, h_)*scale_;
		static const double red = 1;
		static const double green = 1;
		static const double blue = 0;
		static const double opacity = 0.6;


		const auto &state = cellogram::State::state();
		const auto &pts = state.mesh.detected;

		for(int i = 0; i < pts.rows(); ++i)
		{
			const int cx = pts(i, 0)*scale_;
			const int cy = (h_ - pts(i, 1))*scale_;
#ifdef CELLOGRAM_WITH_PNG
			writer_.filledcircle_blend(cx, cy, circle_radius, opacity, red, green, blue);
#else
			writer_ << "<circle cx=\""<<cx<<"\" cy=\""<<cy<<"\" ";
			writer_ << "r=\""<< circle_radius <<"\" ";
			writer_ << "fill=\"rgba("<<((int)(red*255))<<","<<((int)(green*255))<<","<<((int)(blue*255))<<","<<((int)(opacity*255))<<")\" stroke=\"none\"";
			writer_ << "/>\n";
#endif
		}

	}


	void PNGOutput::draw_untangle()
	{

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

			const int p1x = pts(p1, 0)*scale_;
			const int p1y = (h_ - pts(p1, 1))*scale_;

			const int p2x = pts(p2, 0)*scale_;
			const int p2y = (h_ - pts(p2, 1))*scale_;

			const int p3x = pts(p3, 0)*scale_;
			const int p3y = (h_ - pts(p3, 1))*scale_;
#ifdef CELLOGRAM_WITH_PNG
			writer_.triangle(p1x, p1y, p2x, p2y, p3x, p3y, tri_red, tri_green, tri_blue);
#else
			writer_ << "<polygon fill=\"none\" stroke=\"rgb("<<((int)(tri_red*255)) << "," << ((int)(tri_green*255)) << "," << ((int)(tri_blue*255)) <<")\" points=\"";
			writer_ <<p1x<<", "<<p1y<<" ";
			writer_ <<p2x<<", "<<p2y<<" ";
			writer_ <<p3x<<", "<<p3y;
			writer_ << "\"";
			writer_ << "/>\n";
#endif
		}



		const auto &added = mesh.added_by_untangler;
		const auto &deleted = mesh.deleted_by_untangler;

		for(int i = 0; i < added.size(); ++i)
		{
			const int cx = pts(added[i], 0)*scale_;
			const int cy = (h_ - pts(added[i], 1))*scale_;
#ifdef CELLOGRAM_WITH_PNG
			writer_.filledcircle(cx, cy, circle_radius, 1., 0., 0.);
#else
			writer_ << "<circle cx=\""<<cx<<"\" cy=\""<<cy<<"\" ";
			writer_ << "r=\""<< circle_radius <<"\" ";
			writer_ << "fill=\"rgb(255,0,0)\" stroke=\"none\"";
			writer_ << "/>\n";
#endif
		}

		for(int i = 0; i < deleted.rows(); ++i)
		{
			const int cx = deleted(i, 0)*scale_;
			const int cy = (h_ - deleted(i, 1))*scale_;
#ifdef CELLOGRAM_WITH_PNG
			writer_.filledcircle(cx, cy, circle_radius, 1., 0., 0.);
#else
			writer_ << "<circle cx=\""<<cx<<"\" cy=\""<<cy<<"\" ";
			writer_ << "r=\""<< circle_radius <<"\" ";
			writer_ << "fill=\"rgb(255,0,0)\" stroke=\"none\"";
			writer_ << "/>\n";
#endif
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

				const int p1x = pts(p1, 0)*scale_;
				const int p1y = (h_ - pts(p1, 1))*scale_;

				const int p2x = pts(p2, 0)*scale_;
				const int p2y = (h_ - pts(p2, 1))*scale_;

				const int p3x = pts(p3, 0)*scale_;
				const int p3y = (h_ - pts(p3, 1))*scale_;
#ifdef CELLOGRAM_WITH_PNG
				writer_.filledtriangle_blend(p1x, p1y, p2x, p2y, p3x, p3y, bad_region_opacity, 1., 0., 0.);
#else
				writer_ << "<polygon stroke=\"none\" fill=\"rgba(255, 0, 0, "<<((int)(bad_region_opacity*255)) <<")\" points=\"";
				writer_ <<p1x<<", "<<p1y<<" ";
				writer_ <<p2x<<", "<<p2y<<" ";
				writer_ <<p3x<<", "<<p3y;
				writer_ << "\"";
				writer_ << "/>\n";
#endif
			}
		}

		Eigen::VectorXd energy;
		laplace_energy(pts, tris, energy);


		const std::string added_msg = "Added vertices: " + std::to_string(added.size());
		const std::string deleted_msg = "Deleted vertices: " + std::to_string(deleted.rows());
		const std::string bad_msg = "Bad regions: " + std::to_string(n_bad);

		const std::string max_msg = "Max energy: " + std::to_string(energy.maxCoeff());
		const std::string min_msg = "Min energy: " + std::to_string(energy.minCoeff());


#ifdef CELLOGRAM_WITH_PNG
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
#else
		writer_ << "<text x=\""<<0<<"\" y=\""<<h_*scale_ + 5 + 20<<"\" font-family=\"Helvetica\" font-size=\""<<font_size<<"\" ";
		writer_ << "fill=\"rgb("<<(added.size() == 0 ? 0 : 255)<<", 0, 0)\" stroke=\"none\">";
		writer_ << added_msg <<"</text>\n";

		writer_ << "<text x=\""<<0<<"\" y=\""<<h_*scale_ + 35 + 20<<"\" font-family=\"Helvetica\" font-size=\""<<font_size<<"\" ";
		writer_ << "fill=\"rgb("<<(deleted.rows() == 0 ? 0 : 255)<<", 0, 0)\" stroke=\"none\">";
		writer_ << deleted_msg <<"</text>\n";


		writer_ << "<text x=\""<<0<<"\" y=\""<<h_*scale_ + 75 + 20<<"\" font-family=\"Helvetica\" font-size=\""<<font_size<<"\" ";
		writer_ << "fill=\"rgb("<<(n_bad == 0 ? 0 : 255)<<", 0, 0)\" stroke=\"none\">";
		writer_ << bad_msg <<"</text>\n";


		writer_ << "<text x=\""<<0<<"\" y=\""<<h_*scale_ + 115 + 20<<"\" font-family=\"Helvetica\" font-size=\""<<font_size<<"\" ";
		writer_ << "fill=\"rgb(0, 0, 0)\" stroke=\"none\">";
		writer_ << max_msg <<"</text>\n";

		writer_ << "<text x=\""<<0<<"\" y=\""<<h_*scale_ + 145 + 20<<"\" font-family=\"Helvetica\" font-size=\""<<font_size<<"\" ";
		writer_ << "fill=\"rgb(0, 0, 0)\" stroke=\"none\">";
		writer_ << min_msg <<"</text>\n";
#endif
	}
}
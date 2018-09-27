////////////////////////////////////////////////////////////////////////////////
#pragma once

#include <pngwriter.h>

#include <string>

namespace cellogram {

    class PNGOutput
    {
    public:
        PNGOutput(const double scale = 5);

        void init(const std::string &name, const int w, const int h);

        void draw_image();
        void draw_detection();
        void draw_untangle();


        void save() {
            writer_.close();
        }

    private:
        pngwriter writer_;
        double scale_;
        int h_;
        int w_;

        static const int text_offset = 200;
    };

} // namespace cellogram

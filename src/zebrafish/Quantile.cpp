#include <zebrafish/Quantile.h>

#include <queue>
#include <math.h>

namespace zebrafish {

double QuantileImage(const zebrafish::image_t &image, double q, int layerBegin, int layerEnd) {

    int i, j, depth;
    int N = image.size();
    int M = image[0].size();  // rows() * cols()
    int num = floor((1.0-q) * N * M);  // target number for quantile(q)
    std::priority_queue<double, std::vector<double>, std::greater<double> > heap;

    if (layerBegin == -1 || layerEnd == -1) {
        layerBegin = 0;
        layerEnd = N-1;
    }

    for (depth=layerBegin; depth<=layerEnd; depth++) {

        const Eigen::MatrixXd &layer = image[depth];  // alias
        for (i=0; i<layer.size(); i++) {

            if (heap.size() < num)
                heap.push(layer(i));
            else if (layer(i) > heap.top()) {
                heap.pop();
                heap.push(layer(i));
            }
        }
    }

    // DEBUG only
    // printf("Quantile %f = %f", q, heap.top());
    return heap.top();
}

}  // namespace zebrafish

/** Utilities
 */

#ifndef WIRECELLAUX_UTIL
#define WIRECELLAUX_UTIL

#include "WireCellIface/ITensorSet.h"
#include "WireCellIface/ITensor.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellUtil/Exceptions.h"

#include <Eigen/Core>

#include <iostream>

namespace {
}
namespace WireCell {
    namespace Aux {

        // TODO this function seems too complicated...
        template <typename ElementType>
        inline std::string dump(const WireCell::ITensor::pointer iten, std::vector<size_t> loc, const size_t limit = 10)
        {
            std::stringstream ss;
            auto shape = iten->shape();
            if (shape.size() > loc.size()) {
                loc.push_back(0);
                ss << "[";
                bool skip = false;
                if (shape[loc.size() - 1] > limit) {
                    skip = true;
                }
                for (size_t i = 0; i < shape[loc.size() - 1]; ++i) {
                    if (i >= limit / 2 && skip == true) {
                        ss << " ... ";
                        if (loc.size() != shape.size()) {
                            ss << "\n";
                        }
                        i = shape[loc.size() - 1] - limit / 2;
                        skip = false;
                    }
                    loc.back() = i;
                    ss << dump<ElementType>(iten, loc, limit);
                }
                ss << "]";
                if (loc.size() == 1 && loc[0] == shape[0] - 1) {
                    return ss.str();
                }
                if (loc[loc.size() - 2] != shape[loc.size() - 2] - 1) {
                    for (int i = 0; i < (int) shape.size() - (int) loc.size() + 1; ++i) {
                        ss << "\n";
                    }
                }
            }
            else if (shape.size() == loc.size()) {
                auto data = (ElementType *) iten->data();
                auto ind = loc.back();
                for (size_t d = 0; d < loc.size() - 1; ++d) {
                    size_t plus = loc[d];
                    for (size_t i = d + 1; i < shape.size(); ++i) {
                        plus *= shape[i];
                    }
                    ind += plus;
                }
                ss << data[ind] << " ";
            }
            else {
                THROW(WireCell::ValueError() << WireCell::errmsg{"shape.size() < loc.size()"});
            }
            return ss.str();
        }

        template <typename ElementType>
        inline std::string dump(const ITensor::pointer iten, const size_t limit = 10)
        {
            std::stringstream ss;
            ss << "ITensor: ";
            Json::FastWriter jwriter;
            ss << "shape: [";
            for (auto l : iten->shape()) {
                ss << l << " ";
            }
            ss << "]\n";
            ss << dump<ElementType>(iten, {}, limit);
            ss << "\n";

            return ss.str();
        }

        inline std::string dump(const ITensorSet::pointer itens)
        {
            std::stringstream ss;
            ss << "ITensorSet: ";
            Json::FastWriter jwriter;
            ss << itens->ident() << ", " << jwriter.write(itens->metadata());
            for (auto iten : *itens->tensors()) {
                const auto &md = iten->metadata();
                ss << "tag: " << md["tag"];
                ss << ", type: " << md["type"];
                ss << ", shape: [";
                for (auto l : iten->shape()) {
                    ss << l << " ";
                }
                ss << "]\n";
            }
            return ss.str();
        }

        template <typename ElementType>
        inline Eigen::Array<ElementType, Eigen::Dynamic, Eigen::Dynamic> itensor_to_eigen_array(
            const ITensor::pointer iten)
        {
            if (!iten) {
                THROW(ValueError() << errmsg{"(iten == null"});
            }

            if (iten->shape().size() != 2) {
                THROW(ValueError() << errmsg{"iten->shape().size()!=2"});
            }

            // Eigen storage is col major by default
            auto nrows = iten->shape()[0];
            auto ncols = iten->shape()[1];
            Eigen::Map<const Eigen::Array<ElementType, Eigen::Dynamic, Eigen::Dynamic> > arr(
                (const ElementType *) iten->data(), nrows, ncols);

            return arr;
        }

        template <typename ElementType>
        inline SimpleTensor *eigen_array_to_simple_tensor(
            const Eigen::Array<ElementType, Eigen::Dynamic, Eigen::Dynamic> &arr)
        {
            std::vector<size_t> shape = {(size_t) arr.rows(), (size_t) arr.cols()};
            SimpleTensor *st = new SimpleTensor(shape, (ElementType*)nullptr);
            auto dst = (ElementType *) st->data();
            auto src = (ElementType *) arr.data();
            size_t size = sizeof(ElementType) * arr.rows() * arr.cols();
            memcpy(dst, src, size);

            return st;
        }

        template <typename ElementType>
        inline ITensor::pointer eigen_array_to_itensor(
            const Eigen::Array<ElementType, Eigen::Dynamic, Eigen::Dynamic> &arr)
        {
            auto st = eigen_array_to_simple_tensor<ElementType>(arr);
            return ITensor::pointer(st);
        }

    };  // namespace Aux
}  // namespace WireCell

#endif  // WIRECELLSIG_UTIL

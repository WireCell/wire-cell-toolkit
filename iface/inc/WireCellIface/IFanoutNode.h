#ifndef WIRECELL_IFANOUTNODE
#define WIRECELL_IFANOUTNODE

#include "WireCellIface/INode.h"

#include <boost/any.hpp>
#include <vector>

namespace WireCell {

    /** A node which accepts one object of the given InputType and
     *  produces N of the given OutputType.
     */
    class IFanoutNodeBase : public INode {
       public:
        typedef std::shared_ptr<IFanoutNodeBase> pointer;

        virtual ~IFanoutNodeBase();

        typedef std::vector<boost::any> any_vector;

        /// The calling signature:
        virtual bool operator()(const boost::any& anyin, any_vector& anyout) = 0;

        virtual NodeCategory category() { return fanoutNode; }
    };

    // This converts between any and typed.
    template <typename InputType, typename OutputType, int FanoutMultiplicity = 3>
    class IFanoutNode : public IFanoutNodeBase {
       public:
        typedef InputType input_type;
        typedef OutputType output_type;
        typedef std::shared_ptr<const InputType> input_pointer;
        typedef std::shared_ptr<const OutputType> output_pointer;
        typedef std::vector<output_pointer> output_vector;

        virtual ~IFanoutNode() {}

        virtual bool operator()(const boost::any& anyin, any_vector& anyv)
        {
            const input_pointer& in = boost::any_cast<const input_pointer&>(anyin);
            output_vector outv;
            bool ok = (*this)(in, outv);
            if (!ok) return false;
            const size_t mult = output_types().size();  // don't use FanoutMultiplicity
            anyv.resize(mult);
            for (size_t ind = 0; ind < mult; ++ind) {
                anyv[ind] = outv[ind];
            }
            return true;
        }

        // The typed interface
        virtual bool operator()(const input_pointer& in, output_vector& outv) = 0;

        // Return the names of the types this node takes as input.
        virtual std::vector<std::string> input_types() { return std::vector<std::string>{typeid(input_type).name()}; }
        // Return the names of the types this node produces as output.
        // Note: if subclass wants to supply FaninMultiplicity at
        // construction time, this needs to be overridden.
        virtual std::vector<std::string> output_types()
        {
            std::vector<std::string> ret(FanoutMultiplicity, std::string(typeid(output_type).name()));
            return ret;
        }
    };


}  // namespace WireCell

#endif

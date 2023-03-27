/**
   WireSchema holds "wire geometry" data structures and functions.

   The top-most structure is a "store" holding a hierarchy of
   detectors of anodes of faces of planes of wires of points.

   The functions:

   - generate() :: parameterized wire store generator

   - load() :: loader of WCT "wires file" following the "wire data
   schema" (typically as compressed JSON) producing a "store" object.

   - dump() :: save a store to file.

   - validate() :: validate a "store" object against WCT conventions. 

   Content validity conventions are described in

       util/docs/wire-schema.org.

   See also the Python module wirecell.util.wires.schema and family of
   CLI commands: "wirecell-util convert-*".
 */
#ifndef WIRECELLUTIL_WIREPLANES
#define WIRECELLUTIL_WIREPLANES

#include "WireCellUtil/Point.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/BoundingBox.h"
#include <memory>

namespace WireCell {

    namespace WireSchema {

        // IWire
        struct Wire {
            int ident;
            int channel;
            int segment;
            Point tail, head;  // end points, direction of signal to channel
        };

        struct Plane {
            int ident;
            std::vector<int> wires;
        };

        struct Face {
            int ident;
            std::vector<int> planes;
        };

        struct Anode {
            int ident;
            std::vector<int> faces;
        };

        struct Detector {
            int ident;
            std::vector<int> anodes;
        };

        struct StoreDB {
            std::vector<Detector> detectors;
            std::vector<Anode> anodes;
            std::vector<Face> faces;
            std::vector<Plane> planes;
            std::vector<Wire> wires;
        };

        /// Return reference to the plane reached by the hierarchy
        /// path of idents.  Any missing objects are created in the
        /// store.
        Plane& get_append(StoreDB& store, int pid, int fid, int aid, int did);

        /// Return reference to the object in the array "arr" with the
        /// ident, appending a new one if needed.
        template<typename Type>
        Type& get_append(std::vector<Type>& arr, int ident) {
            auto it = std::find_if(arr.begin(), arr.end(),
                                   [&](const Type& obj) { return obj.ident == ident; });
            if (it == arr.end()) {
                arr.emplace_back(Type{ident});
                return arr.back();
            }
            return *it;
        }

        /// Return reference to the object in the array "arr" at one
        /// of the given indices in "inds" with the given "ident",
        /// appending a new object and its index if needed.
        template<typename Type>
        Type& get_append(std::vector<Type>& arr, std::vector<int>& inds, int ident)
        {
            auto it = std::find_if(inds.begin(), inds.end(),
                                   [&](int ind) { return arr[ind].ident == ident; });
            if (it == inds.end()) {
                inds.push_back(arr.size());
                arr.emplace_back(Type{ident});
                return arr.back();
            }
            return arr[*it];
        }


        // Access store via shared pointer to allow for caching of underlying data.
        typedef std::shared_ptr<const StoreDB> StoreDBPtr;

        // Bolt on some const functions to the underlying and shared store.
        class Store {
            StoreDBPtr m_db;

          public:
            Store();  // underlying store will be null!
            Store(StoreDBPtr db);
            Store(const Store& other);  // copy ctro
            Store& operator=(const Store& other);

            // Access underlying data store as shared pointer.
            StoreDBPtr db() const;

            const std::vector<Detector>& detectors() const;
            const std::vector<Anode>& anodes() const;
            const std::vector<Face>& faces() const;
            const std::vector<Plane>& planes() const;
            const std::vector<Wire>& wires() const;

            const Anode& anode(int ident) const;

            std::vector<Anode> anodes(const Detector& detector) const;
            std::vector<Face> faces(const Anode& anode) const;
            std::vector<Plane> planes(const Face& face) const;

            // Return ordered wires.
            std::vector<Wire> wires(const Plane& plane) const;

            // Return the smallest axis-aligned box bounding the
            // object.  If region is true, the box bounds the wire
            // regions and not just the wire centers.  Wire-center
            // bounds are returned by default.
            BoundingBox bounding_box(const Anode& anode, bool region=false) const;
            BoundingBox bounding_box(const Face& face, bool region=false) const;
            BoundingBox bounding_box(const Plane& plane, bool region=false) const;

            // Return wire and pitch direction unit vectors.  This
            // tries to be robust in the face of imprecise wire data.
            Ray wire_pitch(const Plane& plane) const;

            // Return pitch distance vector averaged over wires in plane.
            Vector mean_pitch(const Plane& plane) const;

            // Return wire distance vector averaged over wires in plane.
            Vector mean_wire(const Plane& plane) const;

            std::vector<int> channels(const Plane& plane) const;
        };

        /**
           Load WCT "wires file".

           A number of cumulative corrections can be applied.

           - laod=1 :: The store represents data as-is from the file.

           - order=2 :: Correct order of wire-in-plane (WIP) and wire
             endpoints.

             Normally, this orders by ascending wire center Z value
             but in the case wires are nearly parallel with the Z axis
             this orders by ascending wire center Y value.

             Normally, this orders (tail,head) endpoints by ascending
             Y value but in the case of Z-parallel wires (tail,head)
             is ordered by DEscending Z.

             These ordering rules assure that the diretion of
             increasing WIP is same as increasing pitch vector and
             that X x W = P.

           - direction=3 :: Rotate wires about their centers so that
             they are all parallel.  The common wire direction is
             taken as the average over the original wires, rotated
             into the Y-Z plane.  Wire length and centers are held
             fixed.

           - pitch=4 :: Translate wires along a common pitch direction
             so that they become uniformly distributed.  The common
             pitch is taken as the average over all wires rotated into
             the Y-Z plane.  The center Y/Z of the central wire at WIP
             = nwires/2 is kept fixed and X is set to the average of
             all center X.
           
        */
        enum struct Correction { empty, load, order, direction, pitch, ncorrections };
        inline Correction& operator++(Correction& c) {
            const int in = static_cast<int>(c) + 1;
            return c = static_cast<Correction>( in );
        }
        inline Correction& operator--(Correction& c) {
            const int in = static_cast<int>(c) - 1;
            return c = static_cast<Correction>( in );
        }

        /// Generate a plane of wires adding them to the store and plane.
        ///
        /// - store :: A Store to fill
        ///
        /// - plane :: A Plane to fill, see get_append() for an easy
        ///   way to provide a Plane.
        ///
        /// - pitch :: A ray parallel to the Y-Z plane with tail point
        ///   at the center of wire zero, head on the neighboring wire
        ///   and direction perpendicular to both wire directions.
        ///
        /// - bounds :: A ray from opposite corners of a rectangular
        ///   box on which wire endpoints terminate.
        ///
        /// - wid :: The ident number of wire zero.  Subsequent
        ///   wire.ident are counted from this number in order of
        ///   wire-in-plane.
        ///
        /// Return the number of wires generated and appended to the
        /// store.wires array and their indices appended to the
        /// plane.wires attribute of the Plane.  Order of both appends
        /// are identical and in order of increasing wire-in-plane.
        ///
        /// The wire.channel and wire.segment attributes of new wires
        /// are zero.  The user is free to modify the store contents
        /// after calling.
        ///
        /// ValueError is thrown when input requirements are not met.
        int generate(StoreDB& store, Plane& plane,
                     const Ray& pitch, const Ray& bounds, int wid=0);

        /// Load file into store performing correction
        Store load(const char* filename, Correction correction = Correction::pitch);

        /// Dump store to file
        void dump(const char* filename, const Store& store);


        /// Return only if store is considered valid else throw ValueError.
        ///
        /// Logging will print out descriptions of invalid state at "error" level.
        ///
        /// The "repsilon" sets a unitless relative threshold used to
        /// identify imprecision.
        void validate(const Store& store, double repsilon = 1e-6, bool fail_fast=false);

        /**
           Return vector of ray pairs that categorize each layer
           of a face as suitable for use in constructing a
           RayGrid::Coordinates.
            
           First two elements are horizontal and vertical bounds
           and subsequent elements bound wire0 region for each
           wire plane in the face.
        */
        ray_pair_vector_t ray_pairs(const Store& store, const Face& face, bool region=true);

        /**
           Return ray pairs bounding horiz/vert active area.
        */
        ray_pair_vector_t ray_pairs_active(const Store& store, const Face& face, bool region=true);


    }  // namespace WireSchema

}  // namespace WireCell
#endif

// Caution: this file is generated.  See the .org file of the same name.

#include "WireCellUtil/PointTree.h"

#include "WireCellUtil/PointTesting.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree; // for "Points"
using namespace WireCell::PointTesting;     // make_janky_track()
using namespace spdlog;                     // for debug() etc.

  TEST_CASE("point tree example simple point cloud")
  {
      Dataset pc({
          {"x", Array({1.0, 1.0, 1.0})},
          {"y", Array({2.0, 1.0, 3.0})},
          {"z", Array({1.0, 4.0, 1.0})}});
      
      // Each array is size 3 and thus the PC has that major axis size
      CHECK( pc.size_major() == 3 );
      
      // Accessing a single element in the PC takes two steps:
      // get the array by name, get the element by index.
      auto arr = pc.get("x");
      CHECK( arr );
      
      // We must supply a C++ type in order to extract a value.
      CHECK( arr->element<double>(0) == 1.0 );

    Dataset::selection_t sel = pc.selection({"x","y","z"});

    // The selection must be the same size as the list of names.
    CHECK( sel.size() == 3 );

    // Get first array ("x") and the value at its index 0.
    CHECK( sel[0]->element<double>(0) == 1.0 );

      // Defaults to point elements of type double
      coordinate_array ca(sel);
      CHECK( ca.ndims() == sel.size() );
  
      // Iterate over the selection, point by point.
      size_t count = 0;
      for (const auto& cpt : ca) {
          // Each cpt is a "column" down the rows of selected arrays.
          CHECK( cpt.size() == sel.size() );
          CHECK( cpt[0] == sel[0]->element(count) );
          ++count;
      }
  }

using node_ptr = std::unique_ptr<Points::node_t>;

  TEST_CASE("point tree example nodeless points")
  {
      Points p( { {"3d", Dataset({
              {"x", Array({1.0, 1.0, 1.0})},
              {"y", Array({2.0, 1.0, 3.0})},
              {"z", Array({1.0, 4.0, 1.0})}}) } });
  
      // Normally, we can get the node that holds this value
      // but in this example that node does not exist.
      CHECK( p.node() == nullptr );
  
      // But we can access the local point clouds.
      auto& lpcs = p.local_pcs();
      CHECK( lpcs.size() == 1 );
      CHECK( lpcs.find("3d") != lpcs.end() );
  }

  TEST_CASE("point tree example single node")
  {
      // Normally, we would not create a node on the stack
      // as it could never be a child because as a child
      // must be held by unique_ptr.  
      Points::node_t n(Points( { {"3d", Dataset({
              {"x", Array({1.0, 1.0, 1.0})},
              {"y", Array({2.0, 1.0, 3.0})},
              {"z", Array({1.0, 4.0, 1.0})}}) } }));
  
      // The node directly exposes its value as a data member  
      Points& p = n.value;
  
      // And we can go full circle to get a pointer to the value's node.
      Points::node_t* nptr = p.node();
  
      // This time, that node must exist because we just made it
      CHECK( nptr == &n );
  }

  static
  Points::node_ptr make_simple_pctree()
  {
      // We will return this unique pointer to node
      Points::node_ptr root = std::make_unique<Points::node_t>();
  
      // Insert first child with a set of named points clouds 
      // containing one point cloud build from a track.
      auto* n1 = root->insert(Points({ {"3d", make_janky_track()} }));
      REQUIRE(n1 != nullptr);
  
      // Insert a second child with a point cloud
      // from a different track.
      auto* n2 = root->insert(Points({ {"3d", make_janky_track(
                          Ray(Point(-1, 2, 3), Point(1, -2, -3)))} }));
      REQUIRE(n2 != nullptr);
  
      return root;
  }

// Loop over children and dump info about them.
static void dump_children(const node_ptr& node)
{
    // Loop over children node "Points" values 
    for (const auto& cval : node->child_values()) {
        // The named point clouds held by this node value.
        const auto& pcs = cval.local_pcs();

        debug("child node at {} with {} local point clouds",
              (void*)cval.node(), pcs.size());

        // loop over the set of point clouds
        for (const auto& [name, pc] : pcs) {
            debug("\tchild has pc named \"{}\" with {} points",
                  name, pc.size_major());
        }
    }
}

TEST_CASE("point tree example simple tree operations")
{
    auto root = make_simple_pctree();
    dump_children(root);

    // Find iterator to first child in child list
    CHECK( root->nchildren() == 2 );
    auto cptr = root->children().front();

    // Remove it as a child and we get ownership as unique_ptr.
    auto cuptr = root->remove(cptr);
    CHECK( root->nchildren() == 1 );

    // We can add that orphaned child back.  This transfers ownership.
    auto* cptr2 = root->insert(std::move(cuptr));
    CHECK( root->nchildren() == 2 );

    // But, we caught the return and now have a loaned bare pointer
    CHECK( cptr2 );

    // We should now see the reverse order as above dump.
    dump_children(root);
}

  TEST_CASE("point tree example scope")
  {
      // A default scope can be constructed.
      // It is a null scope as it would match no local 
      // point cloud arrays despite having unlimited depth.
      Scope s;
      CHECK( s.pcname == "" );
      CHECK( s.coords.empty() );
      CHECK( s.depth == 0 );
  
      // Some non-empty scopes.  Note the case sensitivity.
      Scope s0{ "pcname", {"x","y","z"}, 0};
      Scope s1{ "pcname", {"x","y","z"}, 1};
      Scope s2{ "PCNAME", {"x","y","z"}, 0};
      Scope s3{ "pcname", {"X","Y","Z"}, 0};
      Scope sc{ "pcname", {"x","y","z"}, 0};
  
      // Scopes can be compared for equality
      CHECK( s0 == s0 );
      CHECK( s0 == sc );
      CHECK( s0 != s1 );
      CHECK( s0 != s2 );
      CHECK( s0 != s3 );
  
      // A scope has a std::hash().
      CHECK( s0.hash() == sc.hash() );
      CHECK( s0.hash() != s1.hash() );
      CHECK( s0.hash() != s2.hash() );
      CHECK( s0.hash() != s3.hash() );
  
      // So, it may be used as a key.
      std::unordered_map<Scope, size_t> m;
      m[s1] = 1;
      m[s2] = 2;
      CHECK( m[s0] == 0 );          // size_t default value
      CHECK( m[s1] == 1 );
      CHECK( m[s2] == 2 );
  
      // One can also print a scope
      debug("Here is a scope: {}", s0);
  }

  TEST_CASE("point tree example scoped point cloud")
  {
      auto root = make_simple_pctree();
  
      // Specify point cloud name and the arrays to select.  
      // We take the default depth number of 0.
      Scope scope{ "3d", {"x","y","z"} };
  
      // The corresponding scoped view. 
      auto const & sv = root->value.scoped_view(scope);
  
      CHECK(sv.npoints() > 0);
  }

  TEST_CASE("point tree example scoped nodes")
  {
      auto root = make_simple_pctree();
  
      Scope scope{ "3d", {"x","y","z"} };
      auto const & sv = root->value.scoped_view(scope);
  
      // Vector of pointers to n-ary tree nodes
      auto const& snodes = sv.nodes();
      CHECK(snodes.size() > 0);
      for (const auto& node : snodes) {
          const auto& lpcs = node->value.local_pcs();
          REQUIRE(node);
          CHECK(lpcs.find("3d") != lpcs.end());
      }
  }

    TEST_CASE("point tree example scoped nodes")
    {
        auto root = make_simple_pctree();
  
        Scope scope{ "3d", {"x","y","z"} };
        auto const & sv = root->value.scoped_view(scope);
  
        auto const& snodes = sv.nodes();
        auto const& spcs = sv.pcs();
        CHECK(spcs.size() == snodes.size());
  
        // A scoped point cloud is a vector of
        // references to node-local point clouds
        for (const Dataset& pc : spcs) {
            debug("pc {} arrays and {} points",
                  pc.size(), pc.size_major());
        }
  }

  #include "WireCellUtil/DisjointRange.h"
  TEST_CASE("point tree example scoped point cloud disjoint")
  {
      // As above.
      auto root = make_simple_pctree();
      auto spc = root->value.scoped_view(Scope{ "3d", {"x","y","z"} }).pcs();
  
      size_t npoints = 0;
      for (const Dataset& pc : spc) {
          npoints += pc.size_major();
      }
  
      // This part is a bit verbose, see text.
      // First we select a subset of arrays from each point cloud.
      std::vector<Dataset::selection_t> sels;
      for (Dataset& pc : spc) {
          sels.push_back(pc.selection({"x","y","z"}));
      }
      // Then we wrap them as coordinate points.
      using point_array = coordinate_array<double>;
      std::vector<point_array> pas;
      for (auto& sel : sels) {
          pas.emplace_back(sel);
      }
      // Finally we join them together as a disjoint range
      using points_t = disjoint_range<point_array>;
      using point_type = points_t::value_type;
  
      points_t points;
      for (auto& pa : pas) {
          points.append(pa);
      }
  
      // Finally, we can perform a simple iteration
      // as if we had a single contiguous selection.
      CHECK( points.size() == npoints );
  
      for (points_t::iterator pit = points.begin();
          pit != points.end(); ++pit) {
          point_type& pt = *pit;
          CHECK( pt.size() == 3);
          for (size_t dim=0; dim<3; ++dim) {
              debug("pt[{}][{}]={}", std::distance(points.begin(), pit),
                    dim, pt.at(dim));
          }
      }
  
      // We know the scoped PC has two local PCs.
      // Pick an index in either
      const size_t in0 = 0.25*npoints; // in local PC at index 0
      const size_t in1 = 0.75*npoints; // in local PC at index 1
  
      // Iterator to the first point
      auto beg = points.begin();
  
      // Iterators to points inside the flat disjoint_range.  
      // Eg, as we will see returned later from a k-d tree query.
      points_t::iterator pit0 = beg + in0;
      CHECK( pit0 == beg + in0 );
      points_t::iterator pit1 = beg + in1;
      CHECK( pit0 != pit1 );
  
      for (size_t dim=0; dim<3; ++dim) {
          debug("dim={}: {} {}", dim, pit0->at(dim), pit1->at(dim));
      }
  
      // Full circle, find that these points are provided 
      // by the first and second major range.
      const auto djind0 = pit0.index();
      const auto djind1 = pit1.index();
  
      CHECK( 0 == djind0.first );
      CHECK( 1 == djind1.first );
  
      // Can also use major/minor indices to get elements
      // with checked and unchecked lookups.
      {
          auto& qt0 = points.at(djind0);
          auto& qt1 = points.at(djind1);
          for (size_t dim=0; dim<3; ++dim) {
              CHECK( pit0->at(dim) == qt0.at(dim) );
              CHECK( pit1->at(dim) == qt1.at(dim) );
          }
      }
  }

  TEST_CASE("point tree example scoped k-d tree")
  {
      auto root = make_simple_pctree();
  
      // Form a k-d tree query over a scoped point cloud.
      Scope scope = { "3d", {"x","y","z"} };
      const auto& sv = root->value.scoped_view(scope);
      const auto& skd = sv.kd();
      const auto& points = skd.points();

      // Some query point.
      const std::vector<double> origin = {0,0,0};
  
      // Find three nearest neighbors.
      const size_t nnn = 10;    
      REQUIRE(points[0].size() >= nnn);

      auto knn = skd.knn(nnn, origin);
      CHECK( knn.size() == nnn );
  
      // Get the underlying scoped PC if you want to go back to the full dataset.
      // const auto& spc = sv.pcs();
  
      // Loop over results and refrence back to original 
      // scoped PC at both major and minor range level.
      for (const auto& [point_index, metric] : knn) {
          const auto node_index = skd.major_index(point_index);
          // point-in-node
          const auto pin_index = skd.minor_index(point_index);
          /// if you want to get back to the full dataset
          // const Dataset& pc = spc[node_index];
          const size_t ndim = scope.coords.size();
          // Iterate over the point's dimensions and compare 
          // what is in the PC and what is in the iterator.
          for (size_t dim=0; dim<ndim; ++dim) {
              std::string const& name = scope.coords[dim];
              debug("\tdim={} is named {} and holds value {} at point_index {} from node {} local point {}",
                    dim, name, points[dim][point_index], point_index, node_index, pin_index);
          }
      }
  }

  #include "WireCellUtil/Point.h"
  TEST_CASE("point tree example point means")
  {
      Dataset pc({
          {"q", Array({1.0, 2.0, 5.0})},
          {"x", Array({1.0, 1.0, 1.0})},
          {"y", Array({2.0, 1.0, 3.0})},
          {"z", Array({1.0, 4.0, 1.0})}});
      
      auto xyz = pc.selection({"x","y","z"});
      coordinate_range coords(xyz);

      // Find the mean point as a Point type
      Point mu;
      mean_point(mu, coords);
      debug("mean point: ({},{},{})", mu.x(), mu.y(), mu.z());
      CHECK(mu.x() == 1.0);
      CHECK(mu.y() == 2.0);
      CHECK(mu.z() == 2.0);

      // Find the weighted mean point as a std::vector type.  We must pre-size
      // it.
      std::vector<double> wmu(3);
      mean_point(wmu, coords, pc.get("q")->elements());
      debug("weighted mean point: ({},{},{})", wmu[0], wmu[1], wmu[2]);
      CHECK(wmu[0] == 1.0);
      CHECK(wmu[1] == 2.375);
      CHECK(wmu[2] == 1.75);
  }

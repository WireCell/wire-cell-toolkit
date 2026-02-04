#include "WireCellHio/Tensors.h"
#include "WireCellHio/HIO.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Configuration.h"

#include <vector>
#include <cstdint>
#include <boost/filesystem.hpp>

using namespace WireCell;
using namespace WireCell::Hio;
namespace fs = boost::filesystem;

TEST_SUITE("hio_tensors") {

    TEST_CASE("write and read ITensor with int16") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_int16.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create a SimpleTensor
        std::vector<int16_t> data = {1, 2, 3, 4, 5, 6};
        ITensor::shape_t shape = {2, 3};
        auto tensor = std::make_shared<Aux::SimpleTensor>(shape, data.data());

        // Write tensor
        CHECK_NOTHROW(write_itensor(file_id, tensor, "/tensor/int16"));

        // Read tensor back
        ITensor::pointer read_tensor;
        CHECK_NOTHROW(read_tensor = read_itensor(file_id, "/tensor/int16"));

        // Verify properties
        REQUIRE(read_tensor);
        CHECK(read_tensor->shape() == shape);
        CHECK(read_tensor->dtype() == "i2");
        CHECK(read_tensor->size() == data.size() * sizeof(int16_t));

        // Verify data
        const int16_t* read_data = reinterpret_cast<const int16_t*>(read_tensor->data());
        for (size_t i = 0; i < data.size(); ++i) {
            CHECK(read_data[i] == data[i]);
        }

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read ITensor with float") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_float.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create a SimpleTensor
        std::vector<float> data = {1.5f, 2.5f, 3.5f, 4.5f};
        ITensor::shape_t shape = {2, 2};
        auto tensor = std::make_shared<Aux::SimpleTensor>(shape, data.data());

        // Write tensor
        write_itensor(file_id, tensor, "/tensor/float");

        // Read tensor back
        auto read_tensor = read_itensor(file_id, "/tensor/float");

        // Verify properties
        REQUIRE(read_tensor);
        CHECK(read_tensor->shape() == shape);
        CHECK(read_tensor->dtype() == "f4");

        // Verify data
        const float* read_data = reinterpret_cast<const float*>(read_tensor->data());
        for (size_t i = 0; i < data.size(); ++i) {
            CHECK(read_data[i] == doctest::Approx(data[i]));
        }

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read ITensor with double") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_double.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create a SimpleTensor
        std::vector<double> data = {1.123456789, 2.987654321, 3.14159265};
        ITensor::shape_t shape = {3};
        auto tensor = std::make_shared<Aux::SimpleTensor>(shape, data.data());

        // Write tensor
        write_itensor(file_id, tensor, "/tensor/double");

        // Read tensor back
        auto read_tensor = read_itensor(file_id, "/tensor/double");

        // Verify properties
        REQUIRE(read_tensor);
        CHECK(read_tensor->shape() == shape);
        CHECK(read_tensor->dtype() == "f8");

        // Verify data
        const double* read_data = reinterpret_cast<const double*>(read_tensor->data());
        for (size_t i = 0; i < data.size(); ++i) {
            CHECK(read_data[i] == doctest::Approx(data[i]));
        }

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read ITensor with int32") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_int32.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create a SimpleTensor
        std::vector<int32_t> data = {100, 200, 300, 400};
        ITensor::shape_t shape = {4};
        auto tensor = std::make_shared<Aux::SimpleTensor>(shape, data.data());

        // Write tensor
        write_itensor(file_id, tensor, "/tensor/int32");

        // Read tensor back
        auto read_tensor = read_itensor(file_id, "/tensor/int32");

        // Verify properties
        REQUIRE(read_tensor);
        CHECK(read_tensor->shape() == shape);
        CHECK(read_tensor->dtype() == "i4");

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read ITensor with int64") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_int64.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create a SimpleTensor
        std::vector<int64_t> data = {1000000000LL, 2000000000LL};
        ITensor::shape_t shape = {2};
        auto tensor = std::make_shared<Aux::SimpleTensor>(shape, data.data());

        // Write tensor
        write_itensor(file_id, tensor, "/tensor/int64");

        // Read tensor back
        auto read_tensor = read_itensor(file_id, "/tensor/int64");

        // Verify properties
        REQUIRE(read_tensor);
        CHECK(read_tensor->shape() == shape);
        CHECK(read_tensor->dtype() == "i8");

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read ITensor with metadata") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_metadata.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create a SimpleTensor with metadata
        std::vector<float> data = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
        ITensor::shape_t shape = {2, 3};

        Configuration metadata;
        metadata["description"] = "Test tensor";
        metadata["version"] = 1;
        metadata["scale"] = 2.5;
        metadata["active"] = true;
        metadata["tags"] = Json::arrayValue;
        metadata["tags"].append("test");
        metadata["tags"].append("example");

        auto tensor = std::make_shared<Aux::SimpleTensor>(shape, data.data(), metadata);

        // Write tensor
        write_itensor(file_id, tensor, "/tensor/with_metadata");

        // Read tensor back
        auto read_tensor = read_itensor(file_id, "/tensor/with_metadata");

        // Verify metadata
        REQUIRE(read_tensor);
        Configuration read_meta = read_tensor->metadata();

        CHECK(read_meta["description"].asString() == "Test tensor");
        CHECK(read_meta["version"].asInt() == 1);
        CHECK(read_meta["scale"].asDouble() == doctest::Approx(2.5));
        CHECK(read_meta["active"].asBool() == true);
        REQUIRE(read_meta["tags"].isArray());
        CHECK(read_meta["tags"].size() == 2);
        CHECK(read_meta["tags"][0].asString() == "test");
        CHECK(read_meta["tags"][1].asString() == "example");

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read multi-dimensional ITensor") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_3d.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create 3D tensor
        std::vector<double> data;
        for (int i = 0; i < 2 * 3 * 4; ++i) {
            data.push_back(static_cast<double>(i));
        }
        ITensor::shape_t shape = {2, 3, 4};
        auto tensor = std::make_shared<Aux::SimpleTensor>(shape, data.data());

        // Write tensor
        write_itensor(file_id, tensor, "/tensor/3d");

        // Read tensor back
        auto read_tensor = read_itensor(file_id, "/tensor/3d");

        // Verify
        REQUIRE(read_tensor);
        CHECK(read_tensor->shape() == shape);
        CHECK(read_tensor->size() == data.size() * sizeof(double));

        // Verify data
        const double* read_data = reinterpret_cast<const double*>(read_tensor->data());
        for (size_t i = 0; i < data.size(); ++i) {
            CHECK(read_data[i] == doctest::Approx(data[i]));
        }

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write ITensor with null pointer throws error") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_null.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        ITensor::pointer null_tensor;
        CHECK_THROWS_AS(write_itensor(file_id, null_tensor, "/tensor/null"), IOError);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("multiple ITensors in same file") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_multiple.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create and write multiple tensors
        std::vector<int16_t> data1 = {1, 2, 3};
        auto tensor1 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{3}, data1.data());
        write_itensor(file_id, tensor1, "/tensors/t1");

        std::vector<float> data2 = {1.1f, 2.2f, 3.3f, 4.4f};
        auto tensor2 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{2, 2}, data2.data());
        write_itensor(file_id, tensor2, "/tensors/t2");

        std::vector<double> data3 = {10.0, 20.0};
        auto tensor3 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{2}, data3.data());
        write_itensor(file_id, tensor3, "/tensors/t3");

        // Read all back
        auto read1 = read_itensor(file_id, "/tensors/t1");
        auto read2 = read_itensor(file_id, "/tensors/t2");
        auto read3 = read_itensor(file_id, "/tensors/t3");

        // Verify dtypes
        CHECK(read1->dtype() == "i2");
        CHECK(read2->dtype() == "f4");
        CHECK(read3->dtype() == "f8");

        // Verify shapes
        CHECK(read1->shape() == ITensor::shape_t{3});
        CHECK(read2->shape() == ITensor::shape_t{2, 2});
        CHECK(read3->shape() == ITensor::shape_t{2});

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("read non-existent ITensor path throws error") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_notfound.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        show_errors(false);
        CHECK_THROWS_AS(read_itensor(file_id, "/does/not/exist"), IOError);
        show_errors(true);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("ITensor roundtrip preserves all properties") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensor_roundtrip.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create tensor with all properties
        std::vector<float> data = {1.1f, 2.2f, 3.3f, 4.4f, 5.5f, 6.6f, 7.7f, 8.8f};
        ITensor::shape_t shape = {2, 2, 2};

        Configuration metadata;
        metadata["units"] = "meters";
        metadata["timestamp"] = 12345678;

        auto original = std::make_shared<Aux::SimpleTensor>(shape, data.data(), metadata);

        // Write and read
        write_itensor(file_id, original, "/roundtrip");
        auto restored = read_itensor(file_id, "/roundtrip");

        // Verify everything matches
        CHECK(restored->shape() == original->shape());
        CHECK(restored->dtype() == original->dtype());
        CHECK(restored->size() == original->size());
        CHECK(restored->element_size() == original->element_size());

        // Verify metadata
        CHECK(restored->metadata()["units"].asString() == "meters");
        CHECK(restored->metadata()["timestamp"].asInt() == 12345678);

        // Verify data byte-by-byte
        const std::byte* orig_data = original->data();
        const std::byte* rest_data = restored->data();
        for (size_t i = 0; i < original->size(); ++i) {
            CHECK(rest_data[i] == orig_data[i]);
        }

        close(file_id);
        fs::remove(tmpfile);
    }

    // ITensorSet tests

    TEST_CASE("write and read ITensorSet with multiple tensors") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensorset_basic.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create tensors
        std::vector<int16_t> data1 = {1, 2, 3, 4};
        auto tensor1 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{2, 2}, data1.data());

        std::vector<float> data2 = {1.1f, 2.2f, 3.3f};
        auto tensor2 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{3}, data2.data());

        std::vector<double> data3 = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
        auto tensor3 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{2, 3}, data3.data());

        // Create tensor set
        auto tensors = std::make_shared<ITensor::vector>();
        tensors->push_back(tensor1);
        tensors->push_back(tensor2);
        tensors->push_back(tensor3);

        Configuration metadata;
        metadata["description"] = "Test tensor set";

        auto tensorset = std::make_shared<Aux::SimpleTensorSet>(42, metadata, tensors);

        // Write tensor set
        CHECK_NOTHROW(write_itensorset(file_id, tensorset, "/tensorset/basic"));

        // Read tensor set back
        ITensorSet::pointer read_set;
        CHECK_NOTHROW(read_set = read_itensorset(file_id, "/tensorset/basic"));

        // Verify properties
        REQUIRE(read_set);
        CHECK(read_set->ident() == 42);
        CHECK(read_set->metadata()["description"].asString() == "Test tensor set");

        // Verify tensors
        auto read_tensors = read_set->tensors();
        REQUIRE(read_tensors);
        CHECK(read_tensors->size() == 3);

        // Verify first tensor
        CHECK((*read_tensors)[0]->dtype() == "i2");
        CHECK((*read_tensors)[0]->shape() == ITensor::shape_t{2, 2});

        // Verify second tensor
        CHECK((*read_tensors)[1]->dtype() == "f4");
        CHECK((*read_tensors)[1]->shape() == ITensor::shape_t{3});

        // Verify third tensor
        CHECK((*read_tensors)[2]->dtype() == "f8");
        CHECK((*read_tensors)[2]->shape() == ITensor::shape_t{2, 3});

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read ITensorSet with metadata") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensorset_metadata.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create tensor
        std::vector<int32_t> data = {100, 200, 300};
        auto tensor = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{3}, data.data());

        auto tensors = std::make_shared<ITensor::vector>();
        tensors->push_back(tensor);

        // Create metadata
        Configuration metadata;
        metadata["experiment"] = "test";
        metadata["run_number"] = 12345;
        metadata["temperature"] = 273.15;
        metadata["active"] = true;
        metadata["tags"] = Json::arrayValue;
        metadata["tags"].append("cosmic");
        metadata["tags"].append("beam");

        auto tensorset = std::make_shared<Aux::SimpleTensorSet>(999, metadata, tensors);

        // Write and read
        write_itensorset(file_id, tensorset, "/tensorset/with_metadata");
        auto read_set = read_itensorset(file_id, "/tensorset/with_metadata");

        // Verify metadata
        REQUIRE(read_set);
        CHECK(read_set->ident() == 999);
        auto read_meta = read_set->metadata();
        CHECK(read_meta["experiment"].asString() == "test");
        CHECK(read_meta["run_number"].asInt() == 12345);
        CHECK(read_meta["temperature"].asDouble() == doctest::Approx(273.15));
        CHECK(read_meta["active"].asBool() == true);
        REQUIRE(read_meta["tags"].isArray());
        CHECK(read_meta["tags"].size() == 2);
        CHECK(read_meta["tags"][0].asString() == "cosmic");
        CHECK(read_meta["tags"][1].asString() == "beam");

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read empty ITensorSet") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensorset_empty.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create empty tensor set
        auto tensors = std::make_shared<ITensor::vector>();
        auto tensorset = std::make_shared<Aux::SimpleTensorSet>(0, Configuration(), tensors);

        // Write and read
        write_itensorset(file_id, tensorset, "/tensorset/empty");
        auto read_set = read_itensorset(file_id, "/tensorset/empty");

        // Verify
        REQUIRE(read_set);
        CHECK(read_set->ident() == 0);
        auto read_tensors = read_set->tensors();
        REQUIRE(read_tensors);
        CHECK(read_tensors->size() == 0);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write ITensorSet with null pointer throws error") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensorset_null.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        ITensorSet::pointer null_set;
        CHECK_THROWS_AS(write_itensorset(file_id, null_set, "/tensorset/null"), IOError);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("multiple ITensorSets in same file") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensorset_multiple.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create first tensor set
        std::vector<int16_t> data1 = {1, 2, 3};
        auto tensor1 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{3}, data1.data());
        auto tensors1 = std::make_shared<ITensor::vector>();
        tensors1->push_back(tensor1);
        auto set1 = std::make_shared<Aux::SimpleTensorSet>(1, Configuration(), tensors1);

        // Create second tensor set
        std::vector<float> data2 = {4.0f, 5.0f};
        auto tensor2 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{2}, data2.data());
        auto tensors2 = std::make_shared<ITensor::vector>();
        tensors2->push_back(tensor2);
        auto set2 = std::make_shared<Aux::SimpleTensorSet>(2, Configuration(), tensors2);

        // Write both
        write_itensorset(file_id, set1, "/sets/set1");
        write_itensorset(file_id, set2, "/sets/set2");

        // Read both back
        auto read1 = read_itensorset(file_id, "/sets/set1");
        auto read2 = read_itensorset(file_id, "/sets/set2");

        // Verify
        CHECK(read1->ident() == 1);
        CHECK(read2->ident() == 2);
        CHECK(read1->tensors()->size() == 1);
        CHECK(read2->tensors()->size() == 1);
        CHECK((*read1->tensors())[0]->dtype() == "i2");
        CHECK((*read2->tensors())[0]->dtype() == "f4");

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("ITensorSet roundtrip preserves all properties") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensorset_roundtrip.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create tensor set with various properties
        std::vector<double> data1 = {1.1, 2.2, 3.3, 4.4};
        Configuration meta1;
        meta1["name"] = "tensor1";
        auto tensor1 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{2, 2}, data1.data(), meta1);

        std::vector<int64_t> data2 = {100, 200, 300};
        Configuration meta2;
        meta2["name"] = "tensor2";
        auto tensor2 = std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{3}, data2.data(), meta2);

        auto tensors = std::make_shared<ITensor::vector>();
        tensors->push_back(tensor1);
        tensors->push_back(tensor2);

        Configuration set_metadata;
        set_metadata["experiment"] = "roundtrip_test";
        set_metadata["version"] = 123;

        auto original = std::make_shared<Aux::SimpleTensorSet>(777, set_metadata, tensors);

        // Write and read
        write_itensorset(file_id, original, "/roundtrip_set");
        auto restored = read_itensorset(file_id, "/roundtrip_set");

        // Verify set properties
        CHECK(restored->ident() == original->ident());
        CHECK(restored->metadata()["experiment"].asString() == "roundtrip_test");
        CHECK(restored->metadata()["version"].asInt() == 123);

        // Verify tensors
        auto orig_tensors = original->tensors();
        auto rest_tensors = restored->tensors();
        REQUIRE(rest_tensors->size() == orig_tensors->size());

        // Verify first tensor
        CHECK((*rest_tensors)[0]->shape() == (*orig_tensors)[0]->shape());
        CHECK((*rest_tensors)[0]->dtype() == (*orig_tensors)[0]->dtype());
        CHECK((*rest_tensors)[0]->metadata()["name"].asString() == "tensor1");

        // Verify second tensor
        CHECK((*rest_tensors)[1]->shape() == (*orig_tensors)[1]->shape());
        CHECK((*rest_tensors)[1]->dtype() == (*orig_tensors)[1]->dtype());
        CHECK((*rest_tensors)[1]->metadata()["name"].asString() == "tensor2");

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("read non-existent ITensorSet path throws error") {
        fs::path tmpfile = fs::temp_directory_path() / "test_itensorset_notfound.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        show_errors(false);
        CHECK_THROWS_AS(read_itensorset(file_id, "/does/not/exist"), IOError);
        show_errors(true);

        close(file_id);
        fs::remove(tmpfile);
    }

}  // TEST_SUITE

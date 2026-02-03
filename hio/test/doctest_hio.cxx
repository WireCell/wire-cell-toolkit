#include "WireCellHio/HIO.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Configuration.h"

#include <vector>
#include <cstdint>
#include <cstddef>
#include <string>
#include <boost/filesystem.hpp>

using namespace WireCell;
using namespace WireCell::Hio;
namespace fs = boost::filesystem;

TEST_SUITE("hio") {

    TEST_CASE("file open and close") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_open.h5";

        // Create a new file with truncate mode
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);
        CHECK(file_id >= 0);

        // Close the file
        CHECK_NOTHROW(close(file_id));

        // Open existing file in read-only mode
        hid_t file_id2 = open(tmpfile.native(), FileMode::rdonly);
        CHECK(file_id2 >= 0);
        close(file_id2);

        // Cleanup
        fs::remove(tmpfile);
    }

    TEST_CASE("file open with excl mode") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_excl.h5";

        // Create file first
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);
        close(file_id);

        // Try to open with excl mode - should fail
        CHECK_THROWS_AS(open(tmpfile.native(), FileMode::excl), IOError);

        // Cleanup
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read dataset with int16") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_int16.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Prepare test data
        std::vector<int16_t> write_data = {1, 2, 3, 4, 5, 6};
        std::vector<int64_t> shape = {2, 3};

        // Write dataset
        CHECK_NOTHROW(write_dataset(file_id, write_data, shape, "/test/data"));

        // Read dataset back
        std::vector<int16_t> read_data;
        std::vector<int64_t> read_shape;
        CHECK_NOTHROW(read_dataset(file_id, read_data, read_shape, "/test/data"));

        // Verify shape
        REQUIRE(read_shape.size() == shape.size());
        CHECK(read_shape[0] == shape[0]);
        CHECK(read_shape[1] == shape[1]);

        // Verify data
        REQUIRE(read_data.size() == write_data.size());
        for (size_t i = 0; i < write_data.size(); ++i) {
            CHECK(read_data[i] == write_data[i]);
        }

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read dataset with int32") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_int32.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        std::vector<int32_t> write_data = {100, 200, 300, 400};
        std::vector<int64_t> shape = {2, 2};

        write_dataset(file_id, write_data, shape, "/data32");

        std::vector<int32_t> read_data;
        std::vector<int64_t> read_shape;
        read_dataset(file_id, read_data, read_shape, "/data32");

        REQUIRE(read_data.size() == write_data.size());
        CHECK(read_data == write_data);
        CHECK(read_shape == shape);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read dataset with int64") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_int64.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        std::vector<int64_t> write_data = {1000000, 2000000, 3000000};
        std::vector<int64_t> shape = {3};

        write_dataset(file_id, write_data, shape, "/data64");

        std::vector<int64_t> read_data;
        std::vector<int64_t> read_shape;
        read_dataset(file_id, read_data, read_shape, "/data64");

        CHECK(read_data == write_data);
        CHECK(read_shape == shape);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read dataset with float") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_float.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        std::vector<float> write_data = {1.5f, 2.5f, 3.5f, 4.5f};
        std::vector<int64_t> shape = {2, 2};

        write_dataset(file_id, write_data, shape, "/float_data");

        std::vector<float> read_data;
        std::vector<int64_t> read_shape;
        read_dataset(file_id, read_data, read_shape, "/float_data");

        REQUIRE(read_data.size() == write_data.size());
        for (size_t i = 0; i < write_data.size(); ++i) {
            CHECK(read_data[i] == doctest::Approx(write_data[i]));
        }
        CHECK(read_shape == shape);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read dataset with double") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_double.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        std::vector<double> write_data = {1.123456789, 2.987654321, 3.14159265};
        std::vector<int64_t> shape = {3};

        write_dataset(file_id, write_data, shape, "/double_data");

        std::vector<double> read_data;
        std::vector<int64_t> read_shape;
        read_dataset(file_id, read_data, read_shape, "/double_data");

        REQUIRE(read_data.size() == write_data.size());
        for (size_t i = 0; i < write_data.size(); ++i) {
            CHECK(read_data[i] == doctest::Approx(write_data[i]));
        }

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write dataset with raw bytes") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_bytes.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create raw byte data (e.g., int32 values)
        std::vector<int32_t> values = {10, 20, 30, 40};
        const std::byte* byte_data = reinterpret_cast<const std::byte*>(values.data());
        std::vector<int64_t> shape = {2, 2};

        write_dataset(file_id, byte_data, shape, DataType::int32, "/raw_data");

        // Read back as typed data
        std::vector<int32_t> read_data;
        std::vector<int64_t> read_shape;
        read_dataset(file_id, read_data, read_shape, "/raw_data");

        CHECK(read_data == values);
        CHECK(read_shape == shape);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read metadata on dataset") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_metadata.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create a dataset first
        std::vector<int32_t> data = {1, 2, 3, 4};
        std::vector<int64_t> shape = {2, 2};
        write_dataset(file_id, data, shape, "/data_with_meta");

        // Write metadata
        Configuration metadata;
        metadata["description"] = "Test dataset";
        metadata["version"] = 1;
        metadata["temperature"] = 273.15;
        metadata["active"] = true;

        CHECK_NOTHROW(write_metadata(file_id, metadata, "/data_with_meta"));

        // Read metadata back
        Configuration read_meta = read_metadata(file_id, "/data_with_meta");

        CHECK(read_meta["description"].asString() == "Test dataset");
        CHECK(read_meta["version"].asInt() == 1);
        CHECK(read_meta["temperature"].asDouble() == doctest::Approx(273.15));
        CHECK(read_meta["active"].asBool() == true);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("write and read metadata on group") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_group_meta.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create a dataset to establish the group path
        std::vector<int32_t> data = {1, 2};
        std::vector<int64_t> shape = {2};
        write_dataset(file_id, data, shape, "/mygroup/data");

        // Write metadata on the group
        Configuration metadata;
        metadata["group_name"] = "MyGroup";
        metadata["purpose"] = "Testing group metadata";

        CHECK_NOTHROW(write_metadata(file_id, metadata, "/mygroup"));

        // Read metadata back
        Configuration read_meta = read_metadata(file_id, "/mygroup");

        CHECK(read_meta["group_name"].asString() == "MyGroup");
        CHECK(read_meta["purpose"].asString() == "Testing group metadata");

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("read metadata from non-existent path") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_no_meta.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Try to read metadata from non-existent path
        Configuration meta = read_metadata(file_id, "/does/not/exist");

        // Should return empty Configuration (null or has no members)
        bool is_empty = meta.isNull() || meta.empty();
        CHECK(is_empty);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("metadata with complex types") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_complex_meta.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create dataset
        std::vector<float> data = {1.0f, 2.0f, 3.0f};
        std::vector<int64_t> shape = {3};
        write_dataset(file_id, data, shape, "/complex_meta_data");

        // Write metadata with arrays and objects
        Configuration metadata;
        metadata["tags"] = Json::arrayValue;
        metadata["tags"].append("tag1");
        metadata["tags"].append("tag2");
        metadata["tags"].append("tag3");

        metadata["params"]["alpha"] = 0.5;
        metadata["params"]["beta"] = 1.5;

        write_metadata(file_id, metadata, "/complex_meta_data");

        // Read back
        Configuration read_meta = read_metadata(file_id, "/complex_meta_data");

        REQUIRE(read_meta["tags"].isArray());
        CHECK(read_meta["tags"].size() == 3);
        CHECK(read_meta["tags"][0].asString() == "tag1");

        REQUIRE(read_meta["params"].isObject());
        CHECK(read_meta["params"]["alpha"].asDouble() == doctest::Approx(0.5));
        CHECK(read_meta["params"]["beta"].asDouble() == doctest::Approx(1.5));

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("multi-dimensional dataset") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_multidim.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Create 3D dataset
        std::vector<double> data;
        for (int i = 0; i < 2 * 3 * 4; ++i) {
            data.push_back(static_cast<double>(i));
        }
        std::vector<int64_t> shape = {2, 3, 4};

        write_dataset(file_id, data, shape, "/multidim");

        std::vector<double> read_data;
        std::vector<int64_t> read_shape;
        read_dataset(file_id, read_data, read_shape, "/multidim");

        CHECK(read_shape == shape);
        CHECK(read_data == data);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("shape mismatch error") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_shape_error.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Data size doesn't match shape
        std::vector<int32_t> data = {1, 2, 3, 4, 5};  // 5 elements
        std::vector<int64_t> shape = {2, 2};           // 4 elements expected

        // Should throw IOError
        CHECK_THROWS_AS(write_dataset(file_id, data, shape, "/bad_shape"), IOError);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("empty shape error") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_empty_shape.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        std::vector<int32_t> data = {1, 2, 3};
        std::vector<int64_t> shape = {};  // Empty shape

        // Should throw IOError
        CHECK_THROWS_AS(write_dataset(file_id, data, shape, "/empty_shape"), IOError);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("nested group paths") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_nested.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Write to deeply nested path
        std::vector<float> data = {1.0f, 2.0f, 3.0f};
        std::vector<int64_t> shape = {3};

        CHECK_NOTHROW(write_dataset(file_id, data, shape, "/level1/level2/level3/data"));

        // Read back
        std::vector<float> read_data;
        std::vector<int64_t> read_shape;
        read_dataset(file_id, read_data, read_shape, "/level1/level2/level3/data");

        CHECK(read_data == data);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("multiple datasets in same file") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_multiple.h5";
        hid_t file_id = open(tmpfile.native(), FileMode::trunc);

        // Write multiple datasets
        std::vector<int16_t> data1 = {1, 2, 3};
        std::vector<int32_t> data2 = {100, 200};
        std::vector<double> data3 = {1.1, 2.2, 3.3, 4.4};

        write_dataset(file_id, data1, {3}, "/set1");
        write_dataset(file_id, data2, {2}, "/set2");
        write_dataset(file_id, data3, {2, 2}, "/set3");

        // Read all back
        std::vector<int16_t> read1;
        std::vector<int32_t> read2;
        std::vector<double> read3;
        std::vector<int64_t> shape1, shape2, shape3;

        read_dataset(file_id, read1, shape1, "/set1");
        read_dataset(file_id, read2, shape2, "/set2");
        read_dataset(file_id, read3, shape3, "/set3");

        CHECK(read1 == data1);
        CHECK(read2 == data2);
        CHECK(read3 == data3);

        close(file_id);
        fs::remove(tmpfile);
    }

    TEST_CASE("read from read-write mode file") {
        fs::path tmpfile = fs::temp_directory_path() / "test_hio_rdrw.h5";

        // Create file and write data
        {
            hid_t file_id = open(tmpfile.native(), FileMode::trunc);
            std::vector<int32_t> data = {42, 43, 44};
            write_dataset(file_id, data, {3}, "/data");
            close(file_id);
        }

        // Open in read-write mode and modify
        {
            hid_t file_id = open(tmpfile.native(), FileMode::rdrw);

            // Read existing data
            std::vector<int32_t> read_data;
            std::vector<int64_t> shape;
            read_dataset(file_id, read_data, shape, "/data");
            CHECK(read_data[0] == 42);

            // Write new dataset
            std::vector<int32_t> new_data = {99, 98};
            write_dataset(file_id, new_data, {2}, "/new_data");

            close(file_id);
        }

        // Verify both datasets exist
        {
            hid_t file_id = open(tmpfile.native(), FileMode::rdonly);

            std::vector<int32_t> data1, data2;
            std::vector<int64_t> shape1, shape2;

            read_dataset(file_id, data1, shape1, "/data");
            read_dataset(file_id, data2, shape2, "/new_data");

            CHECK(data1[0] == 42);
            CHECK(data2[0] == 99);

            close(file_id);
        }

        fs::remove(tmpfile);
    }

}  // TEST_SUITE

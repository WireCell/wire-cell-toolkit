#include <iostream>
#include <fstream>
#include <iomanip>
#include <json/json.h>

int main()
{
    Json::Value value;
    value["pi"] = 3.14159265358979323846;
    value["x"] = 42;
    value["y"] = 42.12345;
    std::ofstream file("test.json");
    if (file.is_open()) {
        Json::StreamWriterBuilder builder;
        builder.settings_["precision"] = 6;
        // builder.settings_["precisionType"] = "decimal"; // option available in jsoncpp 1.9.0
        std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
        writer->write(value, &file);
        file.close();
    }
    else {
        std::cout << "Failed to open file" << std::endl;
    }

    return 0;
}

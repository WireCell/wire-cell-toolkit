// A main program to test regex.


#include <regex>
#include <iostream>


int main(int argc, char* argv[])
{
    // 0:prog 1:text 2:pattern 3:format
    if (argc < 3 || argc > 4) {
        std::cerr << "usage: check-regex \"string\" \"regex\" [\"format\"]\n";
        std::cerr << "2-arg will apply std::regex_match to regex to string\n";
        std::cerr << "3-arg will apply std::regex_replace to regex and format to string\n";
        return 1;
    }
    const std::string text = argv[1];
    const std::string pattern = argv[2];
    const std::regex match(pattern);
    if (argc == 4) {
        const std::string format = argv[3];

        // std::cerr << "\"" << text << "\" -> s/\"" << pattern << "\"/\"" << format << "\"/\n";
        const std::string result = std::regex_replace(text, match, format);
        std::cerr << result << "\n";
        return 0;
    }

    // std::cerr << "\"" << text << "\" -> s/\"" << pattern << "\"/\n";
    if (std::regex_match(text, match)) {
        // std::cerr << "matches\n";
        return 0;
    }
    // std::cerr << "misses\n";
    return 1;
}


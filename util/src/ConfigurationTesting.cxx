#include "WireCellUtil/ConfigurationTesting.h"


const std::string& WireCell::ConfigurationTesting::tracks_json()
{
    static std::string ret = R"(
[
   {
      "data": {
         "step_size": 1,
         "tracks": [
            {
               "charge": -1,
               "ray": {
                  "head": {
                     "x": 100,
                     "y": 10,
                     "z": 10
                  },
                  "tail": {
                     "x": 10,
                     "y": 0,
                     "z": 0
                  }
               },
               "time": 10
            },
            {
               "charge": -2,
               "ray": {
                  "head": {
                     "x": 2,
                     "y": -100,
                     "z": 0
                  },
                  "tail": {
                     "x": 1,
                     "y": 0,
                     "z": 0
                  }
               },
               "time": 120
            },
            {
               "charge": -3,
               "ray": {
                  "head": {
                     "x": 11,
                     "y": -50,
                     "z": -30
                  },
                  "tail": {
                     "x": 130,
                     "y": 50,
                     "z": 50
                  }
               },
               "time": 99
            }
         ]
      },
      "type": "TrackDepos"
   }
]
)";

    return ret;
}

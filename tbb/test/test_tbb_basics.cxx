#include <tbb/tbb.h>
#include <iostream>

int main(int argc, char* argv[]) {
    int numThreads = 1;
    if (argc > 1) {
        numThreads = std::stoi(argv[1]);
    }

    // Set the number of threads using tbb::global_control
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, numThreads);

    tbb::flow::graph graph;

    // Define nodes
    tbb::flow::function_node<int, int> node1(graph, tbb::flow::unlimited, [](const int& input) {
        std::cout << "node 1 input: " << input << std::endl;
        // Perform computation for node1
        return input + 1;
    });

    tbb::flow::function_node<int, int> node2(graph, tbb::flow::unlimited, [](const int& input) {
        std::cout << "node 2 input: " << input << std::endl;
        // Perform computation for node2
        return input + 1;
    });

    // Define edges
    tbb::flow::make_edge(node1, node2);

    tbb::flow::sequencer_node<int> sequencer(graph, [](const int& input) {
        std::cout << "sequencer input: " << input << std::endl;
        // Perform computation for sequencer
        return input;
    });
    tbb::flow::make_edge(node2, sequencer);

    tbb::flow::function_node<int, int> node3(graph, tbb::flow::unlimited, [](const int& input) {
        std::cout << "node 3 input: " << input << std::endl;
        // Perform computation for node3
        return input + 1;
    });
    tbb::flow::make_edge(sequencer, node3);

    // Add input to the graph
    const int limit = 2;
    int count = 0;
    tbb::flow::input_node<int> src(graph, [&](tbb::flow_control& fc) -> int {
        std::cout << "src input" << std::endl;
        if (count++ < limit) {
            return count * 10;
        }
        fc.stop();
        return {};
    });
    tbb::flow::make_edge(src, node1);

    std::cout << "Running graph" << std::endl;
    // Wait for all nodes to finish
    src.activate();
    graph.wait_for_all();

    return 0;
}

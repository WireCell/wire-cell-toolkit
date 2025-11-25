
// this is for the convert_wires_to_channels_extended()
/*
int main() {
    // R=2 rows, W=5 wires
    // The MP value (0-3) for each wire in each row
    torch::Tensor input = torch::tensor({
        {1, 3, 0, 2, 1}, // Row 0
        {2, 1, 3, 0, 0}  // Row 1
    }, torch::kInt32);

    // W=5 wires. Maps each wire to a channel index (0, 1, or 2)
    torch::Tensor indices = torch::tensor({0, 1, 0, 2, 1}, torch::kInt64);

    // W=5 wires. Maps each wire to a segment (0, 1, or 2)
    torch::Tensor segment = torch::tensor({0, 1, 2, 0, 1}, torch::kInt32);
    
    // Expected output shape: (2, 3) -> 2 rows, 3 channels (0, 1, 2)

    // --- Row 0 Analysis ---
    // Wire:    0   1   2   3   4
    // Channel: 0   1   0   2   1
    // Segment: 0   1   2   0   1
    // MP Val:  1   3   0   2   1
    // Bit Pos: 1  (3+4)=7 (0+8)=8 (2+0)=2 (1+4)=5
    // Value:   2^1=2 2^7=128 2^8=256 2^2=4 2^5=32
    // ---------------------------------------------------------------------
    // Channel 0 (Wires 0, 2): 2 + 256 = 258 (Binary: 100000010)
    // Channel 1 (Wires 1, 4): 128 + 32 = 160 (Binary: 10100000)
    // Channel 2 (Wire 3): 4
    // Row 0 Expected: {258, 160, 4}

    // --- Row 1 Analysis ---
    // Wire:    0   1   2   3   4
    // Channel: 0   1   0   2   1
    // Segment: 0   1   2   0   1
    // MP Val:  2   1   3   0   0
    // Bit Pos: 2  (1+4)=5 (3+8)=11 (0+0)=0 (0+4)=4
    // Value:   2^2=4 2^5=32 2^11=2048 2^0=1 2^4=16
    // ---------------------------------------------------------------------
    // Channel 0 (Wires 0, 2): 4 + 2048 = 2052 (Binary: 100000000100)
    // Channel 1 (Wires 1, 4): 32 + 16 = 48 (Binary: 110000)
    // Channel 2 (Wire 3): 1
    // Row 1 Expected: {2052, 48, 1}

    // torch::Tensor result = convert_wires_to_channels_extended(input, indices, segment);
    // std::cout << "Result:\n" << result << std::endl;
    // Expected Output:
    // {{ 258,  160,    4},
    //  {2052,   48,    1}}
    
    return 0;
}
*/

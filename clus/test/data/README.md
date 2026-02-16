## Test Data Directory

This directory contains configuration and input files used for testing the Wire-Cell Toolkit clustering and track fitting modules.

---

### `uboone-mabc_config.json`

**Description:**  
Configuration file from the QLport directory.

**Purpose:**  
Defines detector configuration parameters for the MicroBooNE MABC (Magnetically Assisted Bubble Chamber) system.

**Contents include:**
- Detector geometry definitions
- Electronics and readout configuration
- Processing and reconstruction parameters

These settings are consumed by the reconstruction pipeline to ensure consistent detector modeling and data processing.

---

### `init_first_segment_input.json`

**Description:**  
Debug dump of the initial first-segment input data.

**Generation condition:**  
This file is produced when the environment variable `WCT_DUMP_INIT_FIRST_SEGMENT` is set.

**Initialization workflow captured in this file:**
1. Create a Pattern Recognition (PR) graph to represent Wire-Cell clustering structures
2. Initialize the first track segment from the main cluster using:
   - Pattern-matching algorithms
   - Track-fitting information
3. Register the PR graph with the track fitter to enable multi-track reconstruction
4. Execute the full multi-tracking algorithm
5. Store track-fitting results in the grouping object for downstream processing

**Use case:**  
This file is primarily intended for debugging and validation, allowing developers to:
- Inspect pattern-recognition initialization
- Verify first-segment generation
- Validate track-fitting setup within the Wire-Cell reconstruction pipeline
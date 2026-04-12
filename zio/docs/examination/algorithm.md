# Algorithm and Architecture - zio Package

This document explains the overall architecture, data flow, and algorithmic details of the WireCellZio package.

---

## 1. Overview

The WireCellZio package bridges the Wire-Cell Toolkit (WCT) data processing framework with the ZIO (ZeroMQ I/O) messaging library. Its primary purpose is to enable WCT tensor data to be streamed over a network using ZIO's flow protocol, supporting use cases such as:

- **Remote inference:** Sending detector data to a remote PyTorch model server for deep neural network inference (e.g., ROI finding)
- **Data persistence:** Streaming tensor data to an HDF5 file server via ZIO's flow-file-server
- **Distributed processing:** Connecting WCT pipeline components across processes or machines

The package consists of three functional layers:
1. **Serialization** (`TensUtil`): Convert between WCT `ITensorSet` and ZIO `Message`
2. **Flow transport** (`FlowConfigurable`, `TensorSetSink`, `TensorSetSource`): Manage ZIO flow connections for streaming
3. **Application integration** (`ZioTorchScript`): PyTorch TorchScript inference via ZIO DOMO client

---

## 2. Core Data Model

### ITensorSet (WCT side)
WCT represents multi-dimensional numerical data as `ITensorSet`, which contains:
- An **ident** (integer sequence number)
- A **metadata** object (JSON)
- A vector of **ITensor** objects, each with:
  - Shape (vector of dimension sizes)
  - Element type and size
  - Raw data buffer
  - Per-tensor metadata (JSON)

### zio::Message (ZIO side)
ZIO represents data as `zio::Message`, which contains:
- A **label** (JSON string in the header)
- A **payload** (ZMQ multipart message with binary frames)
- Protocol fields: origin, granule, seqno, level, form

The "TENS" form is ZIO's convention for tensor data, storing tensor metadata in the label and tensor data in payload frames.

---

## 3. Serialization Algorithm (TensUtil)

### pack(): ITensorSet -> zio::Message

```
Input: ITensorSet with N tensors

1. Create a zio::Message with form "TENS"
2. Extract ITensorSet metadata and store in label["TENS"]["metadata"]
3. Serialize the label as JSON string and set as the message label
4. For each tensor i in [0, N):
   a. Convert tensor metadata from jsoncpp to nlohmann::json (via string roundtrip)
   b. Copy tensor raw data into a zio::message_t (ZMQ frame)
   c. Call zio::tens::append() which:
      - Adds the data frame to the message payload
      - Records in label["TENS"]["tensors"][i]:
        - "shape": dimension array
        - "word": element size in bytes
        - "dtype": type character ('f'=float, 'i'=signed int, 'u'=unsigned int)
        - "part": payload frame index
        - "metadata": per-tensor metadata
5. Return the message
```

### unpack(): zio::Message -> ITensorSet

```
Input: zio::Message in TENS form

1. Parse the message label as JSON
2. For each tensor descriptor in label["TENS"]["tensors"]:
   a. Extract shape, part index, word size, dtype
   b. Retrieve the raw data frame from payload at the part index
   c. Type-dispatch based on (dtype, word):
      - 'f',4 -> float    |  'f',8 -> double
      - 'i',1 -> int8_t   |  'i',2 -> int16_t  | 'i',4 -> int32_t  | 'i',8 -> int64_t
      - 'u',1 -> uint8_t  |  'u',2 -> uint16_t | 'u',4 -> uint32_t | 'u',8 -> uint64_t
   d. Allocate a SimpleTensor with the shape
   e. memcpy raw data from the ZMQ frame into the tensor's storage
   f. Convert per-tensor metadata from nlohmann::json to jsoncpp
3. Extract seqno from message as the ITensorSet ident
   (Note: this is acknowledged as a layer violation in the code)
4. Parse set-level metadata from label
5. Return a SimpleTensorSet wrapping the tensors
```

### JSON Format Conversion

Because WCT uses jsoncpp (`Json::Value`) and ZIO uses nlohmann::json (`zio::json`), the code must convert between them. This is done via string serialization:

```
jsoncpp -> string -> nlohmann::json    (wct_to_zio)
nlohmann::json -> string -> jsoncpp    (zio_to_wct)
```

---

## 4. Flow Transport Architecture

### FlowConfigurable Base Class

`FlowConfigurable` provides the common infrastructure for any component that sends or receives data over ZIO flow. It manages:

1. **ZIO Node**: A Zyre-based network node with a unique nickname and origin ID
2. **Port**: A single ZMQ port (CLIENT or SERVER socket type)
3. **Flow**: A credit-based flow control wrapper around the port

#### Configuration Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| timeout | 1000ms | ZMQ operation timeout |
| credit | 10 | Flow credit (max unacknowledged messages) |
| stype | ZMQ_CLIENT | Socket type (CLIENT or SERVER only) |
| nodename | (auto) | Zyre node nickname for discovery |
| portname | "flow" | ZIO port name |
| binds | [] | Explicit bind addresses |
| connects | [] | Connect targets (address or node/port name pairs) |
| headers | {} | Zyre headers for discovery |
| attributes | {} | Flow BOT attributes |
| verbose | (unset) | Enable Zyre verbose logging |

#### Lifecycle

```
1. Constructor(direction, nodename)
   - Sets flow direction: "inject" (receive) or "extract" (send)
   - Creates ZIO node

2. configure(cfg)
   a. Call user_configure() for subclass-specific config
   b. Parse timeout, credit, node parameters
   c. Create port with specified socket type
   d. Bind or connect the port (ephemeral if no addresses given)
   e. Parse Zyre headers
   f. Create zio::Flow object with direction, credit, timeout
   g. Go online (starts Zyre discovery)
   h. Call user_online()
   i. Call post_configure()
   j. Perform BOT (beginning of transmission) handshake
      - BOT carries flow attributes as JSON label

3. operator() [in subclass]
   - Use m_flow->put() or m_flow->get() for data transfer

4. finalize()
   - Send EOT (end of transmission)
   - Go offline
```

### TensorSetSink (WCT -> ZIO)

Direction: "extract" (sends data out from WCT)

```
operator()(ITensorSet in):
  if in == nullptr:    // end-of-stream signal
    if already_had_eos:
      finalize()       // second null = shutdown
      return false
    else:
      send empty "FLOW" message as EOS marker
      set m_had_eos = true
      return true
  else:
    reset m_had_eos
    pack(in) -> zio::Message
    m_flow->put(msg)   // send via flow with credit control
    return true

Exception handling:
  - end_of_transmission: acknowledge EOT, null out flow
  - timeout: finalize and return false
```

### TensorSetSource (ZIO -> WCT)

Direction: "inject" (receives data into WCT)

```
operator()(ITensorSet& out):
  out = nullptr
  msg = m_flow->get()     // receive via flow
  if payload is empty:    // EOS
    return true (out stays nullptr = EOS signal)
  else:
    out = unpack(msg)
    return true

Exception handling:
  - end_of_transmission: acknowledge EOT, return true (as EOS)
  - timeout: finalize and return false
```

---

## 5. ZIO Flow Protocol

The flow protocol provides credit-based flow control on top of ZMQ:

### Message Types
- **BOT** (Beginning of Transmission): Handshake to establish flow. Carries attributes.
- **DAT**: Data message. Carries tensor payload.
- **PAY**: Credit payment. Grants permission to send more DAT messages.
- **EOT** (End of Transmission): Graceful shutdown signal.

### Credit-Based Flow Control
The flow protocol uses a credit system to prevent fast producers from overwhelming slow consumers:

1. The receiver (injector) grants an initial credit of N (default 10) to the sender
2. Each DAT message sent consumes 1 credit
3. When the sender runs out of credit, `put()` blocks until more PAY messages arrive
4. The receiver sends PAY messages as it processes data

This ensures backpressure: if the consumer falls behind, the producer automatically slows down.

### Direction Semantics
- **extract**: The component sends data (is a producer). Uses `m_flow->put()`.
- **inject**: The component receives data (is a consumer). Uses `m_flow->get()`.

Note the naming is from the component's perspective relative to the ZIO network - "extract" means extracting data from WCT to push into ZIO, "inject" means injecting data from ZIO into WCT.

---

## 6. ZioTorchScript: Remote DNN Inference

`ZioTorchScript` implements `ITensorSetFilter`, taking a tensor set in and producing a tensor set out, with the actual computation happening on a remote server.

### Architecture

```
WCT Pipeline -> ZioTorchScript -> [ZMQ] -> DOMO Broker -> TorchScript Server
                                                              |
WCT Pipeline <- ZioTorchScript <- [ZMQ] <- DOMO Broker <-----+
```

### Protocol
Uses ZIO's DOMO (Distributed Object Management Organization) pattern, based on the Majordomo Protocol:
- **Client** (ZioTorchScript): Sends requests, receives replies
- **Broker**: Routes requests to available workers by service name
- **Worker** (external TorchScript server): Processes inference requests

### Per-Call Flow

```
1. Create ZMQ context + socket + DOMO client (per call)
2. Pack input ITensorSet into zio::Message
3. Convert message to ZMQ multipart format
4. Send to DOMO broker with service name (e.g., "torch:dnnroi")
5. Wait for reply (blocking)
6. Convert reply back to zio::Message
7. Unpack to output ITensorSet
8. On failure: sleep wait_time ms, retry (infinite loop)
```

### Configuration
| Parameter | Default | Description |
|-----------|---------|-------------|
| address | tcp://localhost:5555 | DOMO broker address |
| service | torch:dnnroi | Service name for routing |
| model | model.ts | TorchScript model file (unused in client) |
| gpu | true | GPU flag (unused in client, used by server) |
| wait_time | 500ms | Retry delay on failure |

---

## 7. Network Discovery

ZIO uses Zyre (ZeroMQ Realtime Exchange) for automatic peer discovery:

1. Each node goes online with a nickname and optional headers
2. Zyre uses UDP beacons (or gossip) to discover peers on the local network
3. Ports can connect by (nodename, portname) pair instead of explicit addresses
4. This enables dynamic wiring: a TensorSetSink can find a flow-file-server by name without hardcoded addresses

The `binds` and `connects` configuration arrays support three modes:
- **Empty**: Ephemeral bind (random port, discoverable via Zyre)
- **String address**: Explicit bind/connect to a TCP or IPC address
- **Object {nodename, portname}**: Connect via Zyre discovery

---

## 8. Test Architecture

The test suite validates the system at multiple levels:

### Unit Tests
- **test_tens**: Raw ZIO tensor message creation and readback
- **test_tensutil**: Pack/unpack roundtrip for all supported data types (float, double, int8-64, uint8-64)
- **test_flowconfigurable**: IFrame -> ITensorSet -> zio::Message -> ITensorSet pipeline

### Integration Tests
- **test_zioflow**: Three-actor test (giver -> middleman -> taker) with in-process ZMQ transport
- **check_zioflow_{give,take}**: Two-process tests requiring an external ZIO flow-file-server
- **check_client_torch**: End-to-end DNN inference test (requires DOMO broker + TorchScript server)
- **check_gdm**: DOMO client test with echo service

### Test Data Flow (test_zioflow)

```
giver (CLIENT) ──DAT──> middleman (SERVER/SERVER) ──DAT──> taker (CLIENT)
      <──────PAY──────            <──────PAY──────
```

The middleman binds two SERVER ports, the giver and taker connect as CLIENTs. The middleman forwards DAT messages and manages flow credit. The test creates 10 dummy integer tensor sets (2x3x4), sends them through the pipeline, and verifies receipt.

---

## 9. Configuration Example (issue54.jsonnet)

The `issue54.jsonnet` demonstrates a complete WCT+ZIO pipeline:

```
Particle tracks -> Drifter -> Bagger -> [Signal+Noise sim per APA] -> TaggedFrameTensorSet -> ZioTensorSetSink
                                                                                                    |
                                                                                                    v
                                                                                          ZIO flow-file-server
                                                                                                    |
                                                                                                    v
                                                                                              issue54.hdf
```

Key configuration points:
- Multiple parallel APA pipelines fan out from a single depo source
- Each pipeline has its own `ZioTensorSetSink` with a unique stream name (`simsn0`, `simsn1`, ...)
- All sinks connect to a single `issue54server` node by Zyre discovery
- The flow-file-server stores each stream in a separate HDF5 group

The flow-file-server configuration (`check_zioflow.jsonnet`) defines rules matching direction to read/write mode and mapping streams to HDF5 groups.

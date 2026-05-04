# Efficiency Analysis - zio Package

This document examines performance and efficiency concerns in the WireCellZio package.

---

## 1. ZMQ Context + Socket + Client Created Per Inference Call

**File:** `src/ZioTorchScript.cxx:88-90`

```cpp
zmq::context_t ctx;
zmq::socket_t sock(ctx, ZMQ_CLIENT);
zio::domo::Client m_client(sock, get<std::string>(m_cfg, "address", "tcp://localhost:5555"));
```

Every call to `operator()` creates a new ZMQ context (which spawns I/O threads), a new socket, and a new DOMO client. These are then destroyed at the end of the call. ZMQ context creation/destruction is expensive:
- Context creation spawns 1+ I/O threads
- Socket creation involves kernel resource allocation
- DOMO client establishment involves protocol handshake

**Impact:** For a pipeline processing thousands of tensor sets, this overhead is paid on every single invocation. The context, socket, and client should be created once during `configure()` and reused across calls.

**Estimated overhead:** ~1-10ms per call for context/socket lifecycle, depending on OS. For a 10,000-frame job, this adds 10-100 seconds of pure overhead.

---

## 2. JSON Serialization Roundtrip via String for Format Conversion

**File:** `src/TensUtil.cxx:10-20`

```cpp
void Zio::zio_to_wct(const zio::json& zioj, Json::Value& wctj) {
    std::stringstream ss(zioj.dump());
    ss >> wctj;
}
void Zio::wct_to_zio(const Json::Value& wctj, zio::json& zioj) {
    std::stringstream ss;
    ss << wctj;
    zioj = zio::json::parse(ss.str());
}
```

Converting between nlohmann::json (zio) and jsoncpp (WCT) is done by serializing to a string and re-parsing. This involves:
1. Full JSON serialization to string (`dump()` or `operator<<`)
2. String stream allocation and copying
3. Full JSON parsing from string (`parse()` or `operator>>`)

This is done once per tensor in `fill_tensor()` (for metadata) and once per tensor set in both `pack()` and `unpack()`.

**Impact:** For tensor sets with many tensors, the repeated serialize-parse cycle adds up. A direct tree-walking converter would be more efficient, avoiding string allocation and parsing overhead entirely.

---

## 3. Logger Obtained by Name on Every operator() Call

**File:** `src/TensorSetSink.cxx:28`, `src/TensorSetSource.cxx:21`

```cpp
auto log = Log::logger("wctzio");
```

The spdlog logger is looked up by name via string on every call to `operator()`. While spdlog's `get()` uses a mutex-protected map lookup, this is unnecessary per-call overhead. The logger should be obtained once (in the constructor or `configure()`) and stored as a member, as `ZioTorchScript` does with its `l` member.

**Impact:** Minor per-call overhead (mutex acquisition + map lookup). Negligible for low-throughput use but adds up in tight loops.

---

## 4. Full Tensor Data Copy in `pack()`

**File:** `src/TensUtil.cxx:38`

```cpp
zio::tens::append(msg, zio::message_t(ten->data(), ten->size()), ...);
```

The `zio::message_t` constructor copies the tensor data into a new ZMQ message. For large tensors (e.g., a 2400-channel x 6000-tick frame at 4 bytes/float = ~57 MB), this is a full memcpy. The data then gets copied again when sent over ZMQ (into the kernel send buffer).

This is somewhat inherent to the ZMQ message model, but worth noting that the pack path involves at least one full copy of every tensor's data. If zero-copy semantics are available in the ZIO/ZMQ layer, they could be explored.

**Impact:** For large detector frames, this can be significant. A single APA frame copy at ~57 MB takes ~10ms on modern hardware.

---

## 5. Full Tensor Data Copy in `unpack()` via `fill_tensor()`

**File:** `src/TensUtil.cxx:48-50`

```cpp
Aux::SimpleTensor *st = new Aux::SimpleTensor(shape, (TYPE*)nullptr);
auto& store = st->store();
memcpy(store.data(), one.data(), one.size());
```

The unpack path allocates a new tensor storage buffer and copies the entire payload into it. Combined with the copy in `pack()`, a round-trip (pack -> send -> receive -> unpack) involves at minimum 2 full data copies on each side.

**Impact:** Same as #4. The double-copy overhead is inherent to the current architecture but could potentially be reduced with zero-copy message handling.

---

## 6. Polling Spin Loop in `flow_middleman()`

**File:** `src/TestHelpers.cxx:91-99`

```cpp
const auto wait = std::chrono::milliseconds{100};
// ...
while (keep_going) {
    int nhits = poller.wait_all(events, wait);
```

The middleman actor uses a 100ms poll timeout in a loop. When there is no activity, this wakes up 10 times per second to check for events. While this is test code, if the pattern were adopted in production, it would waste CPU cycles on idle polling.

A longer timeout (or blocking poll with a dedicated interrupt mechanism) would be more efficient. The current design is acknowledged in the comment as "a slow spin" (line 89).

**Impact:** Low for tests. Would matter if this pattern were used in production.

---

## 7. String Formatting with `fmt::format` in Hot Path

**File:** `src/TestHelpers.cxx:187-188`

```cpp
const std::string prefix =
    fmt::format("[{} {}]", get<std::string>(cfg, "nodename"), get<std::string>(cfg, "portname"));
```

This is test code, so the impact is minimal. But in the production path, `TensorSetSink::operator()` constructs log messages with `m_node.nick()` and `m_portname` even when the log level would suppress them. The spdlog library does have lazy formatting, but the argument evaluation (including `nick()` calls) still happens.

**Impact:** Negligible in practice. spdlog's macro-based level check avoids formatting when the level is disabled.

---

## 8. `Json::FastWriter` Created Per-Call in `pack()`

**File:** `src/TensUtil.cxx:31`

```cpp
Json::FastWriter jwriter;
msg.set_label(jwriter.write(label));
```

A `Json::FastWriter` is created on the stack every time `pack()` is called. While `FastWriter` is lightweight, a `static` or member-held writer would avoid repeated construction. Similarly, `Json::Reader` is created in `unpack()` (line 132).

Note: `Json::FastWriter` and `Json::Reader` are deprecated in newer jsoncpp in favor of `StreamWriterBuilder`/`CharReaderBuilder`.

**Impact:** Negligible. The construction cost is trivial compared to the actual JSON writing.

---

## 9. Label Object Parsed Twice in `unpack()`

**File:** `src/TensUtil.cxx:64,132-133`

```cpp
auto label = zmsg.label_object();          // parse 1
// ... loop over tensors ...
reader.parse(label[...]["metadata"].dump(), md);  // serialize + parse again
```

The label is parsed from the message once (line 64), then the metadata sub-object is serialized back to string (`dump()`) and re-parsed into a different JSON library format (line 133). This is a specialized case of issue #2 above but worth noting as an additional roundtrip.

**Impact:** Minor. One extra serialize-parse cycle per unpack call.

---

## Summary Table

| Issue | Location | Severity | Fix Complexity |
|-------|----------|----------|----------------|
| ZMQ context per call | ZioTorchScript.cxx | High | Low - move to configure() |
| JSON string roundtrip | TensUtil.cxx | Medium | Medium - write tree walker |
| Logger lookup per call | Sink/Source .cxx | Low | Low - cache as member |
| Tensor data copy (pack) | TensUtil.cxx | Medium | High - requires ZMQ zero-copy |
| Tensor data copy (unpack) | TensUtil.cxx | Medium | High - requires ZMQ zero-copy |
| Polling spin loop | TestHelpers.cxx | Low | Low - test only |
| fmt::format in hot path | TestHelpers.cxx | Negligible | N/A - test only |
| FastWriter per call | TensUtil.cxx | Negligible | Low |
| Double label parse | TensUtil.cxx | Low | Low |

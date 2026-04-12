# Potential Bugs - zio Package

This document examines potential bugs in the WireCellZio package source code.

---

## 1. Infinite Retry Loop with No Exit in ZioTorchScript::operator()

**File:** `src/ZioTorchScript.cxx:94-109`

The retry loop catches *all* exceptions and retries indefinitely:

```cpp
while (!success) {
    try {
        // ... pack, send, recv, unpack ...
        success = true;
    }
    catch (...) {
        std::this_thread::sleep_for(std::chrono::milliseconds(wait_time));
        thread_wait_time += wait_time;
    }
}
```

**Problem:** If the remote service is permanently down, misconfigured, or the data itself causes a repeated error (e.g., serialization bug), this loop runs forever. There is no maximum retry count, no escalating backoff, and no way to detect non-transient failures (e.g., a serialization error will never succeed on retry). The `catch(...)` also silently swallows all exception types including programming errors.

**Severity:** High. This can cause a processing pipeline to hang indefinitely with no diagnostic output.

---

## 2. Uninitialized Member `m_had_eos` in TensorSetSource

**File:** `inc/WireCellZio/TensorSetSource.h:22`

```cpp
bool m_had_eos;
```

In `TensorSetSource`, `m_had_eos` is declared but never initialized in the constructor (`src/TensorSetSource.cxx:12-14`). Compare with `TensorSetSink` which explicitly initializes `m_had_eos(false)`.

**Problem:** The member will have an indeterminate value. While `m_had_eos` does not appear to be read in the current `TensorSetSource` implementation (it's declared but unused), this is itself a secondary issue: the member exists suggesting it was intended to be used for EOS tracking but the logic was never completed.

**Severity:** Low (currently unused), but indicates incomplete implementation.

---

## 3. Raw `new` Without Exception Safety in `unpack()`

**File:** `src/TensUtil.cxx:60`

```cpp
ITensor::vector* itv = new ITensor::vector;
```

This raw pointer is used throughout the function and only wrapped in a shared_ptr at the very end (line 135). If any exception is thrown between lines 60-135 (e.g., during JSON parsing at line 70, `get<>()` calls, or `fill_tensor` calls), the `itv` pointer leaks.

Similarly in `fill_tensor()` (line 48):
```cpp
Aux::SimpleTensor *st = new Aux::SimpleTensor(shape, (TYPE*)nullptr);
```
If `memcpy` or `zio_to_wct` throws, `st` leaks.

**Severity:** Medium. Memory leaks on error paths.

---

## 4. Silent Data Corruption on Unknown dtype/word in `unpack()`

**File:** `src/TensUtil.cxx:77-124`

When `unpack()` encounters an unknown `dtype` or `word` size, it only logs an error but does *not* skip the tensor or throw:

```cpp
else {
    log->error("unknown floating point size: {}", word);
}
```

After logging, execution continues and the tensor is simply omitted from the output vector. However, the tensor indexing is silently shifted - the output `ITensorSet` will have fewer tensors than expected with no indication of which tensor was dropped. Downstream code expecting tensor N at index N will get the wrong tensor.

**Severity:** Medium. Silent data misalignment is worse than a crash.

---

## 5. `pre_flow()` No Longer Guards Against Repeated BOT

**File:** `src/FlowConfigurable.cxx:167`

The header comment (line 48-51) says:
> Subclass must call this before any actual flow. It is safe to call at the top of each execution. If it returns false, no flow is possible.

But the implementation is just:
```cpp
bool Zio::FlowConfigurable::pre_flow() { return true; }
```

The original guard `m_did_bot` (declared at line 65) is never set or checked. The `pre_flow()` method was historically supposed to handle BOT handshake once and guard against re-entry. After refactoring (BOT was moved into `configure()`), `pre_flow()` became a no-op but the comments and `m_did_bot` member were left behind. In `TensorSetSink::operator()` (line 30), `pre_flow()` is called every invocation but does nothing.

**Severity:** Low. The dead code is confusing but not currently harmful since BOT is now done in `configure()`.

---

## 6. Exception Caught by Value Instead of Reference

**File:** `src/TensorSetSink.cxx:49`, `src/TensorSetSink.cxx:75`, `src/TensorSetSource.cxx:35`

```cpp
catch (zio::flow::end_of_transmission) {
```

These catch exception objects by value rather than by reference (`const zio::flow::end_of_transmission&`). This causes object slicing if the exception type has derived classes and involves unnecessary copying.

**Severity:** Low. Functionally works but is non-idiomatic C++ and risks slicing.

---

## 7. `check_client_torch.cxx` Calls Non-Existent Static Methods

**File:** `test/check_client_torch.cxx:135,151`

```cpp
auto msg = Zio::FlowConfigurable::pack(iitens);
// ...
auto oitens = Zio::FlowConfigurable::unpack(msg);
```

`pack()` and `unpack()` are free functions in `WireCell::Zio`, not static methods of `FlowConfigurable`. This test file would not compile as written. The correct calls should be `Zio::pack()` and `Zio::unpack()`.

**Severity:** High (for this test file). The test is broken and cannot be compiled.

---

## 8. ZioTorchScript Creates ZMQ Context Per Call

**File:** `src/ZioTorchScript.cxx:88-90`

```cpp
zmq::context_t ctx;
zmq::socket_t sock(ctx, ZMQ_CLIENT);
zio::domo::Client m_client(sock, ...);
```

A new ZMQ context, socket, and DOMO client are created on *every* invocation of `operator()`. ZMQ contexts are heavyweight (they manage I/O threads). Creating and destroying them per-call is:
1. A performance issue (see efficiency.md)
2. A potential resource exhaustion bug if calls are rapid - ZMQ contexts may not release OS resources immediately upon destruction

**Severity:** Medium. Could lead to file descriptor / thread exhaustion under load.

---

## 9. TensorSetSink: EOS Handling and Double-Call Protocol

**File:** `src/TensorSetSink.cxx:36-63`

The EOS protocol requires the caller to send `nullptr` twice:
- First `nullptr`: sends an empty message as EOS marker, sets `m_had_eos = true`, returns `true`
- Second `nullptr`: calls `finalize()`, returns `false`

This is fragile:
1. If a non-null tensor arrives between the two EOS calls (`m_had_eos` gets reset to `false` at line 66), the EOS state machine resets silently.
2. The caller (pipeline framework) must know to send exactly two nullptrs. This contract is not enforced and not documented in the header.

**Severity:** Medium. The protocol is fragile and can lead to incomplete shutdown if the calling sequence varies.

---

## 10. Potential Data Size Mismatch in `fill_tensor()`

**File:** `src/TensUtil.cxx:48-50`

```cpp
Aux::SimpleTensor *st = new Aux::SimpleTensor(shape, (TYPE*)nullptr);
auto& store = st->store();
memcpy(store.data(), one.data(), one.size());
```

There is no validation that `one.size()` matches `store.size()`. If the ZIO message payload size doesn't match the expected tensor size (shape * sizeof(TYPE)), the `memcpy` will either:
- Read past the source buffer (if `one.size() < store.size()`) - undefined behavior
- Overflow the destination buffer (if `one.size() > store.size()`) - buffer overflow

**Severity:** High. Buffer overflow / out-of-bounds read on malformed input.

---

## 11. Attributes Values Forced to String

**File:** `src/FlowConfigurable.cxx:155-158`

```cpp
for (auto key : cfg["attributes"].getMemberNames()) {
    auto val = cfg["attributes"][key];
    fobj[key] = val.asString();
}
```

All attribute values are converted with `asString()`. If a configuration provides numeric or boolean attributes, `asString()` on a non-string Json::Value will either throw or produce unexpected results (e.g., integer 42 becomes "42"). The `issue54.jsonnet` config uses string values (`"simsn%d"`), but this code would break for non-string attribute values.

**Severity:** Low. Works for current usage but is fragile.

---

## 12. Unused `m_tensors` Member in TensorSetSource

**File:** `inc/WireCellZio/TensorSetSource.h:23`

```cpp
ITensorSet::vector m_tensors;  // current set of depos
```

This member is declared but never used anywhere in `TensorSetSource.cxx`. The comment "current set of depos" suggests it was part of an earlier design. Dead members increase object size and are misleading.

**Severity:** Low. Dead code.

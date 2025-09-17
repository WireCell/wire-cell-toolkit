#include "WireCellTriton/TritonService.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellIface/ITensor.h"         
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

#include <grpcpp/grpcpp.h>
#include <cstring>
#include <numeric>
#include <sstream>

WIRECELL_FACTORY(TritonService,
                 WireCell::Triton::TritonService,
                 WireCell::ITensorForward,
                 WireCell::IConfigurable)

using namespace WireCell;

namespace {
  
  static std::string shape_str(const std::vector<int64_t>& s) {
    std::ostringstream os; os << "[";
    for (size_t i = 0; i < s.size(); ++i) { if (i) os << ","; os << s[i]; }
    os << "]";
    return os.str();
  }
  
  // product helpers
  static size_t numel(const std::vector<size_t>& s) {
    return std::accumulate(s.begin(), s.end(), size_t{1}, std::multiplies<size_t>());
  }
  
  // make a zero ITensorSet with arbitrary [N,C,H,W]
  static WireCell::ITensorSet::pointer make_zero_tensorset(
							   size_t ident, const std::vector<size_t>& shape)
  {
    // default to [1,1,800,600] if shape is odd
    std::vector<size_t> shp = shape;
    if (shp.size() != 4) shp = {1,1,800,600};
    std::vector<float> zeros(numel(shp), 0.0f);
    auto ten = std::make_shared<WireCell::Aux::SimpleTensor>(shp, zeros.data());
    WireCell::ITensor::shared_vector tv(new WireCell::ITensor::vector{ten});
    WireCell::Configuration meta;
    return std::make_shared<WireCell::Aux::SimpleTensorSet>(ident, meta, tv);
  }

}


namespace WireCell::Triton {

  TritonService::TritonService()
    : Aux::Logger("TritonService", "triton")
  {}

  TritonService::~TritonService() {}

  WireCell::Configuration TritonService::default_configuration() const {
    Configuration cfg;
    cfg["url"]         = "localhost:8001";
    cfg["model"]       = "dnn";
    cfg["input_name"]  = "INPUT__0";
    cfg["output_name"] = "OUTPUT__0";
    cfg["soft_fail"]   = true;

    // Optional gRPC channel knobs
    cfg["grpc_max_metadata_bytes"]   = 32*1024*1024;
    cfg["grpc_keepalive_time_ms"]    = 120000;
    cfg["grpc_keepalive_timeout_ms"] = 20000;
    cfg["grpc_http2_bdp_probe"]      = true;
    return cfg;
  }

  void TritonService::configure(const Configuration& cfg)
  {
    m_url         = get<std::string>(cfg, "url", "localhost:8001");
    m_model_name  = get<std::string>(cfg, "model", "dnn");
    m_input_name  = get<std::string>(cfg, "input_name", "INPUT__0");
    m_output_name = get<std::string>(cfg, "output_name", "OUTPUT__0");
    m_soft_fail   = get<bool>(cfg, "soft_fail", true);

    if (m_client) return;

    // channel args to avoid metadata/keepalive issues on large payloads
    const int  max_meta      = get(cfg, "grpc_max_metadata_bytes",   32*1024*1024);
    const int  ka_time_ms    = get(cfg, "grpc_keepalive_time_ms",    120000);
    const int  ka_timeout_ms = get(cfg, "grpc_keepalive_timeout_ms", 20000);
    const bool bdp_probe     = get(cfg, "grpc_http2_bdp_probe",      true);

    grpc::ChannelArguments ch_args;
    ch_args.SetMaxReceiveMessageSize(-1);
    ch_args.SetMaxSendMessageSize(-1);
    ch_args.SetInt(GRPC_ARG_MAX_METADATA_SIZE, max_meta);
    ch_args.SetInt(GRPC_ARG_KEEPALIVE_TIME_MS, ka_time_ms);
    ch_args.SetInt(GRPC_ARG_KEEPALIVE_TIMEOUT_MS, ka_timeout_ms);
    ch_args.SetInt(GRPC_ARG_KEEPALIVE_PERMIT_WITHOUT_CALLS, 1);
    ch_args.SetInt(GRPC_ARG_HTTP2_MAX_PINGS_WITHOUT_DATA, 0);
    ch_args.SetInt(GRPC_ARG_HTTP2_BDP_PROBE, bdp_probe ? 1 : 0);

    triton::client::SslOptions ssl; // no SSL
    auto err = triton::client::InferenceServerGrpcClient::Create(
								 &m_client, m_url, ch_args, /*verbose*/false,
								 /*use_ssl*/false, ssl, /*use_cached_channel*/false);

    if (!err.IsOk()) {
      log->critical("Failed to create Triton gRPC client for {}: {}", m_url, err.Message());
      THROW(RuntimeError() << errmsg{err.Message()});
    }
    log->info("Connected to Triton at {}", m_url);
  }


  WireCell::ITensorSet::pointer TritonService::forward(const ITensorSet::pointer& input) const
  {
    // Derive a reasonable fallback shape early (used by soft_fail exits).
    auto fallback_from_wc = [&](const WireCell::ITensor::pointer& t)->std::vector<size_t> {
      if (!t) return {1,1,800,600};
      const auto& s = t->shape();            // [N,C,H,W] expected
      if (s.size() == 4) return {1,1, (size_t)s[2], (size_t)s[3]};
      return {1,1,800,600};
    };

    if (!input || !input->tensors() || input->tensors()->empty()) {
      log->error("no input tensor");
      return m_soft_fail ? make_zero_tensorset(0, {1,1,800,600}) : nullptr;
    }

    // Expect first tensor, float32, rank-4 (N,C,H,W)
    auto tin = input->tensors()->front();
    const auto& wc_shape = tin->shape();
    if (wc_shape.size() != 4) {
      log->error("incoming ITensor has rank {}, expected 4", wc_shape.size());
      return m_soft_fail ? make_zero_tensorset(input->ident(), fallback_from_wc(tin)) : nullptr;
    }

    // Validate input byte size against shape
    const size_t want_in_bytes =
      (size_t)wc_shape[0] * (size_t)wc_shape[1] * (size_t)wc_shape[2] * (size_t)wc_shape[3] * sizeof(float);
    const size_t in_bytes = tin->size();
    if (in_bytes != want_in_bytes) {
      log->warn("input byte size {} differs from shape product {} bytes (N={},C={},H={},W={})",
		in_bytes, want_in_bytes, wc_shape[0], wc_shape[1], wc_shape[2], wc_shape[3]);
    }
    const uint8_t* in_raw = reinterpret_cast<const uint8_t*>(tin->data());

    
    std::vector<int64_t> triton_in_shape{
      (int64_t)wc_shape[0], (int64_t)wc_shape[1], (int64_t)wc_shape[2], (int64_t)wc_shape[3]
	};

    triton::client::InferInput* triton_input_raw = nullptr;
    auto err = triton::client::InferInput::Create(&triton_input_raw, m_input_name, triton_in_shape, "FP32");
    if (!err.IsOk()) {
      log->error("InferInput::Create failed: {}", err.Message());
      return m_soft_fail ? make_zero_tensorset(input->ident(), fallback_from_wc(tin)) : nullptr;
    }
    std::unique_ptr<triton::client::InferInput> triton_input(triton_input_raw);

    err = triton_input->AppendRaw(in_raw, in_bytes);
    if (!err.IsOk()) {
      log->error("AppendRaw failed: {}", err.Message());
      return m_soft_fail ? make_zero_tensorset(input->ident(), fallback_from_wc(tin)) : nullptr;
    }

    triton::client::InferRequestedOutput* triton_output_raw = nullptr;
    err = triton::client::InferRequestedOutput::Create(&triton_output_raw, m_output_name);
    if (!err.IsOk()) {
      log->error("InferRequestedOutput::Create failed: {}", err.Message());
      return m_soft_fail ? make_zero_tensorset(input->ident(), fallback_from_wc(tin)) : nullptr;
    }
    std::unique_ptr<triton::client::InferRequestedOutput> triton_output(triton_output_raw);

    triton::client::InferOptions options(m_model_name);
    std::vector<triton::client::InferInput*> inputs{ triton_input.get() };
    std::vector<const triton::client::InferRequestedOutput*> outputs{ triton_output.get() };
    triton::client::Headers headers;

    // Measure inference time
    auto t0 = std::chrono::steady_clock::now();
    triton::client::InferResult* result_raw = nullptr;
    err = m_client->Infer(&result_raw, options, inputs, outputs, headers);
    auto t1 = std::chrono::steady_clock::now();
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    if (!err.IsOk() || !result_raw) {
      log->error("Infer() failed: {}", err.Message());
      return m_soft_fail ? make_zero_tensorset(input->ident(), fallback_from_wc(tin)) : nullptr;
    }
    std::unique_ptr<triton::client::InferResult> result(result_raw);

    // Output shape
    std::vector<int64_t> got_shape64;
    err = result->Shape(m_output_name, &got_shape64);
    if (!err.IsOk()) {
      log->error("Shape() failed: {}", err.Message());
      return m_soft_fail ? make_zero_tensorset(input->ident(), fallback_from_wc(tin)) : nullptr;
    }
    // Accept [N,1,H,W] or [N,H,W] → insert C=1
    if (got_shape64.size() == 3) got_shape64.insert(got_shape64.begin() + 1, 1);
    if (got_shape64.size() != 4) {
      log->error("unexpected Triton output rank {} (want 4), shape={}", got_shape64.size(), shape_str(got_shape64));
      return m_soft_fail ? make_zero_tensorset(input->ident(), fallback_from_wc(tin)) : nullptr;
    }
    if (got_shape64[1] != 1) {
      log->warn("output C dimension is {}, expected 1 (continuing)", got_shape64[1]);
    }
    log->debug("Infer took {} ms, in_shape={}, out_shape={}",
               ms, shape_str(triton_in_shape), shape_str(got_shape64));

    // Pull bytes
    const uint8_t* out_raw = nullptr;
    size_t out_bytes = 0;
    err = result->RawData(m_output_name, &out_raw, &out_bytes);
    if (!err.IsOk() || !out_raw) {
      log->error("RawData failed: {}", err.Message());
      return m_soft_fail ? make_zero_tensorset(input->ident(), fallback_from_wc(tin)) : nullptr;
    }

    // Build WCT tensor (deep copy)
    std::vector<size_t> out_wc_shape{
      (size_t)got_shape64[0], (size_t)got_shape64[1],
        (size_t)got_shape64[2], (size_t)got_shape64[3]
	};
    const size_t need_floats = numel(out_wc_shape);
    const size_t got_floats  = out_bytes / sizeof(float);
    if (got_floats < need_floats) {
      log->error("output byte size too small: {} floats, need {}", got_floats, need_floats);
      return m_soft_fail ? make_zero_tensorset(input->ident(), fallback_from_wc(tin)) : nullptr;
    }

    std::vector<float> out_copy(need_floats);
    std::memcpy(out_copy.data(), out_raw, need_floats * sizeof(float));

    auto out_tensor = std::make_shared<WireCell::Aux::SimpleTensor>(out_wc_shape, out_copy.data());
    WireCell::ITensor::shared_vector otv(new WireCell::ITensor::vector{out_tensor});
    WireCell::Configuration meta;
    return std::make_shared<WireCell::Aux::SimpleTensorSet>(input->ident(), meta, otv);
  }

} // namespace WireCell::Triton

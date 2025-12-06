#ifndef WIRECELLTRITON_TRITONSERVICE
#define WIRECELLTRITON_TRITONSERVICE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorForward.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"

// Triton client includes:
#include "grpc_client.h"

namespace WireCell::Triton {
  class TritonService : public Aux::Logger,
    public ITensorForward,
                          public IConfigurable
			  {
			  public:
			    TritonService();
			    virtual ~TritonService();

			    // IConfigurable
			    virtual void configure(const WireCell::Configuration& config) override;
			    virtual WireCell::Configuration default_configuration() const override;

			    // ITensorForward interface
			    virtual ITensorSet::pointer forward(const ITensorSet::pointer& input) const override;

			  private:
			    // Triton server info
			    std::string m_url;
			    std::string m_model_name;
			    std::string m_input_name;
			    std::string m_output_name;

			    // shape, etc, if needed
			    std::vector<int64_t> m_input_shape;
			    std::vector<size_t>  m_output_shape;
			    bool m_soft_fail{true};
			    // mutable client to be able to reuse in a const method
			    mutable std::unique_ptr<triton::client::InferenceServerGrpcClient> m_client;
			  };
} // namespace WireCell::Triton

#endif // WIRECELLTRITON_TRITONSERVICE

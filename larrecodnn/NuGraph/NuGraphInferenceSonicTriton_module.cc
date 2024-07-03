////////////////////////////////////////////////////////////////////////
// Class:       NuGraphInferenceSonicTriton
// Plugin Type: producer (Unknown Unknown)
// File:        NuGraphInferenceSonicTriton_module.cc
//
// Generated at Tue Nov 14 14:41:30 2023 by Giuseppe Cerati using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <array>
#include <limits>
#include <memory>

#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h" //this creates a conflict with torch script if included before it...

class NuGraphInferenceSonicTriton;

using anab::FeatureVector;
using anab::MVADescription;
using recob::Hit;
using recob::SpacePoint;
using std::array;
using std::vector;

namespace {
  template <typename T, typename A>
  int arg_max(std::vector<T, A> const& vec)
  {
    return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
  }

  template <typename T, size_t N>
  void softmax(std::array<T, N>& arr)
  {
    T m = -std::numeric_limits<T>::max();
    for (size_t i = 0; i < arr.size(); i++) {
      if (arr[i] > m) { m = arr[i]; }
    }
    T sum = 0.0;
    for (size_t i = 0; i < arr.size(); i++) {
      sum += expf(arr[i] - m);
    }
    T offset = m + logf(sum);
    for (size_t i = 0; i < arr.size(); i++) {
      arr[i] = expf(arr[i] - offset);
    }
    return;
  }
}

class NuGraphInferenceSonicTriton : public art::EDProducer {
public:
  explicit NuGraphInferenceSonicTriton(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  NuGraphInferenceSonicTriton(NuGraphInferenceSonicTriton const&) = delete;
  NuGraphInferenceSonicTriton(NuGraphInferenceSonicTriton&&) = delete;
  NuGraphInferenceSonicTriton& operator=(NuGraphInferenceSonicTriton const&) = delete;
  NuGraphInferenceSonicTriton& operator=(NuGraphInferenceSonicTriton&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  vector<std::string> planes;
  art::InputTag hitInput;
  art::InputTag spsInput;
  size_t minHits;
  bool debug;
  // vector<vector<float>> avgs;
  // vector<vector<float>> devs;
  bool filterDecoder;
  bool semanticDecoder;
  bool vertexDecoder;
};

NuGraphInferenceSonicTriton::NuGraphInferenceSonicTriton(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , planes(p.get<vector<std::string>>("planes"))
  , hitInput(p.get<art::InputTag>("hitInput"))
  , spsInput(p.get<art::InputTag>("spsInput"))
  , minHits(p.get<size_t>("minHits"))
  , debug(p.get<bool>("debug"))
  , filterDecoder(p.get<bool>("filterDecoder"))
  , semanticDecoder(p.get<bool>("semanticDecoder"))
  , vertexDecoder(p.get<bool>("vertexDecoder"))
{

  // for (size_t ip = 0; ip < planes.size(); ++ip) {
  //   avgs.push_back(p.get<vector<float>>("avgs_" + planes[ip]));
  //   devs.push_back(p.get<vector<float>>("devs_" + planes[ip]));
  // }

  if (filterDecoder) { produces<vector<FeatureVector<1>>>("filter"); }
  //
  if (semanticDecoder) {
    produces<vector<FeatureVector<5>>>("semantic");
    produces<MVADescription<5>>("semantic");
  }
  //
  if (vertexDecoder) { produces<vector<recob::Vertex>>("vertex"); }
}

void NuGraphInferenceSonicTriton::produce(art::Event& e)
{

  art::Handle<vector<Hit>> hitListHandle;
  vector<art::Ptr<Hit>> hitlist;
  if (e.getByLabel(hitInput, hitListHandle)) { art::fill_ptr_vector(hitlist, hitListHandle); }

  std::unique_ptr<vector<FeatureVector<1>>> filtcol(
    new vector<FeatureVector<1>>(hitlist.size(), FeatureVector<1>(std::array<float, 1>({-1.}))));

  std::unique_ptr<vector<FeatureVector<5>>> semtcol(new vector<FeatureVector<5>>(
    hitlist.size(), FeatureVector<5>(std::array<float, 5>({-1., -1., -1., -1., -1.}))));
  std::unique_ptr<MVADescription<5>> semtdes(
    new MVADescription<5>(hitListHandle.provenance()->moduleLabel(),
                          "semantic",
                          {"MIP", "HIP", "shower", "michel", "diffuse"}));

  std::unique_ptr<vector<recob::Vertex>> vertcol(new vector<recob::Vertex>());

  if (debug) std::cout << "Hits size=" << hitlist.size() << std::endl;
  if (hitlist.size() < minHits) {
    if (filterDecoder) { e.put(std::move(filtcol), "filter"); }
    if (semanticDecoder) {
      e.put(std::move(semtcol), "semantic");
      e.put(std::move(semtdes), "semantic");
    }
    if (vertexDecoder) { e.put(std::move(vertcol), "vertex"); }
    return;
  }

  // event id
  int run = e.id().run();
  int subrun = e.id().subRun();
  int event = e.id().event();

  array<int, 3> evtID;
  evtID[0] = run;
  evtID[1] = subrun;
  evtID[2] = event;  

  // hit table
  vector<size_t> hit_table_hit_id;
  vector<size_t> hit_table_local_plane;
  vector<float>  hit_table_local_time;
  vector<size_t> hit_table_local_wire;
  vector<float>  hit_table_integral;
  vector<float>  hit_table_rms;
  for (auto h : hitlist) {
    hit_table_hit_id.push_back(h.key());
    hit_table_local_plane.push_back(h->View());
    hit_table_local_time.push_back(h->PeakTime());
    hit_table_local_wire.push_back(h->WireID().Wire);
    hit_table_integral.push_back(h->Integral());
    hit_table_rms.push_back(h->RMS());
  }


  // Get spacepoints from the event record
  art::Handle<vector<SpacePoint>> spListHandle;
  vector<art::Ptr<SpacePoint>> splist;
  if (e.getByLabel(spsInput, spListHandle)) { art::fill_ptr_vector(splist, spListHandle); }
  // Get assocations from spacepoints to hits
  vector<vector<art::Ptr<Hit>>> sp2Hit(splist.size());
  if (splist.size() > 0) {
    art::FindManyP<Hit> fmp(spListHandle, e, "sps");
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    }
  }

  // space point table
  vector<size_t> spacepoint_table_spacepoint_id;
  vector<int> spacepoint_table_hit_id_u;
  vector<int> spacepoint_table_hit_id_v;
  vector<int> spacepoint_table_hit_id_y;
  for (size_t i = 0; i < splist.size(); ++i) {
    spacepoint_table_spacepoint_id.push_back(i);
    spacepoint_table_hit_id_u.push_back(-1);
    spacepoint_table_hit_id_v.push_back(-1);
    spacepoint_table_hit_id_y.push_back(-1);
    for (size_t j = 0; j < sp2Hit[i].size(); ++j) {
      if (sp2Hit[i][j]->View()==0) spacepoint_table_hit_id_u.back() = sp2Hit[i][j].key();
      if (sp2Hit[i][j]->View()==1) spacepoint_table_hit_id_v.back() = sp2Hit[i][j].key();
      if (sp2Hit[i][j]->View()==2) spacepoint_table_hit_id_y.back() = sp2Hit[i][j].key();
    }
  }

  //Here the input should be sent to Triton

  //This needs to be replaced with the Triton output
  /*
  auto outputs = model.forward(inputs).toGenericDict();

  if (debug) std::cout << "output =" << outputs << std::endl;
  if (semanticDecoder) {
    for (size_t p = 0; p < planes.size(); p++) {
      torch::Tensor s = outputs.at("x_semantic").toGenericDict().at(planes[p]).toTensor();
      for (int i = 0; i < s.sizes()[0]; ++i) {
        size_t idx = idsmap[p][i];
        std::array<float, 5> input({s[i][0].item<float>(),
                                    s[i][1].item<float>(),
                                    s[i][2].item<float>(),
                                    s[i][3].item<float>(),
                                    s[i][4].item<float>()});
        softmax(input);
        FeatureVector<5> semt = FeatureVector<5>(input);
        (*semtcol)[idx] = semt;
      }
      if (debug) {
        for (int j = 0; j < 5; j++) {
          std::cout << "x_semantic category=" << j << " : ";
          for (size_t p = 0; p < planes.size(); p++) {
            torch::Tensor s = outputs.at("x_semantic").toGenericDict().at(planes[p]).toTensor();
            for (int i = 0; i < s.sizes()[0]; ++i)
              std::cout << s[i][j].item<float>() << ", ";
          }
          std::cout << std::endl;
        }
      }
    }
  }
  if (filterDecoder) {
    for (size_t p = 0; p < planes.size(); p++) {
      torch::Tensor f = outputs.at("x_filter").toGenericDict().at(planes[p]).toTensor();
      for (int i = 0; i < f.numel(); ++i) {
        size_t idx = idsmap[p][i];
        std::array<float, 1> input({f[i].item<float>()});
        (*filtcol)[idx] = FeatureVector<1>(input);
      }
    }
    if (debug) {
      std::cout << "x_filter : ";
      for (size_t p = 0; p < planes.size(); p++) {
        torch::Tensor f = outputs.at("x_filter").toGenericDict().at(planes[p]).toTensor();
        for (int i = 0; i < f.numel(); ++i)
          std::cout << f[i].item<float>() << ", ";
      }
      std::cout << std::endl;
    }
  }
  if (vertexDecoder) {
    torch::Tensor v = outputs.at("x_vtx").toGenericDict().at("evt").toTensor()[0];
    double vpos[3];
    vpos[0] = v[0].item<float>();
    vpos[1] = v[1].item<float>();
    vpos[2] = v[2].item<float>();
    vertcol->push_back(recob::Vertex(vpos));
  }
  */

  if (filterDecoder) { e.put(std::move(filtcol), "filter"); }
  if (semanticDecoder) {
    e.put(std::move(semtcol), "semantic");
    e.put(std::move(semtdes), "semantic");
  }
  if (vertexDecoder) { e.put(std::move(vertcol), "vertex"); }
}

DEFINE_ART_MODULE(NuGraphInferenceSonicTriton)

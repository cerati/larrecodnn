////////////////////////////////////////////////////////////////////////
// Class:       NuSliceHitsProducer
// Plugin Type: producer (art v3_06_03)
// File:        NuSliceHitsProducer_module.cc
//
// Generated at Tue May 25 10:39:19 2021 by Giuseppe Cerati using cetskelgen
// from cetlib version v3_11_01.
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

#include <memory>

#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCParticle.h"

class NuSliceHitsProducer;

using HitParticleAssociations =
  art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>;

class NuSliceHitsProducer : public art::EDProducer {
public:
  explicit NuSliceHitsProducer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuSliceHitsProducer(NuSliceHitsProducer const&) = delete;
  NuSliceHitsProducer(NuSliceHitsProducer&&) = delete;
  NuSliceHitsProducer& operator=(NuSliceHitsProducer const&) = delete;
  NuSliceHitsProducer& operator=(NuSliceHitsProducer&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  // Declare member data here.
  std::string fPfpLabel;
  std::string fSliceLabel;
  std::string fHitLabel;
  std::string fHitTruthLabel;
  bool fRecoverHighestNuScoreSlice;
  bool fRecover2ndShower;
  float fVtxDistCut;
  int fMaxHitCut;

  void AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,
                    const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h,
                    std::vector<art::Ptr<recob::PFParticle> > &pfp_v,
                    std::map<unsigned int, unsigned int>& pfpmap);
};

NuSliceHitsProducer::NuSliceHitsProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fPfpLabel(p.get<std::string>("PfpLabel", "pandora"))
  , fSliceLabel(p.get<std::string>("SliceLabel", "pandora"))
  , fHitLabel(p.get<std::string>("HitLabel", "gaushit"))
  , fHitTruthLabel(p.get<std::string>("HitTruthLabel", ""))
  , fRecoverHighestNuScoreSlice(p.get<bool>("RecoverHighestNuScoreSlice"))
  , fRecover2ndShower(p.get<bool>("Recover2ndShower"))
  , fVtxDistCut(p.get<float>("VtxDistCut"))
  , fMaxHitCut(p.get<int>("MaxHitCut"))
// More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>>();
  if (!fHitTruthLabel.empty()) produces<HitParticleAssociations>();

  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void NuSliceHitsProducer::produce(art::Event& e)
{
  // Implementation of required member function here.

  auto outputHits = std::make_unique<std::vector<recob::Hit>>();
  auto outputHitPartAssns = std::make_unique<HitParticleAssociations>();
  art::PtrMaker<recob::Hit> hitPtrMaker(e);

  art::ValidHandle<std::vector<recob::PFParticle>> inputPfp =
    e.getValidHandle<std::vector<recob::PFParticle>>(fPfpLabel);
  auto assocPfpSlice = std::unique_ptr<art::FindManyP<recob::Slice>>(
    new art::FindManyP<recob::Slice>(inputPfp, e, fPfpLabel));
  auto assocPfpMetadata = e.getValidHandle<std::vector<larpandoraobj::PFParticleMetadata>>(fPfpLabel);
  auto assocPfpVertex = std::unique_ptr<art::FindManyP<recob::Vertex>>(
    new art::FindManyP<recob::Vertex>(inputPfp, e, fPfpLabel));
  auto assocPfpCluster = std::unique_ptr<art::FindManyP<recob::Cluster>>(
    new art::FindManyP<recob::Cluster>(inputPfp, e, fPfpLabel));

  art::ValidHandle<std::vector<recob::Slice>> inputSlice =
    e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
  auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(
    new art::FindManyP<recob::Hit>(inputSlice, e, fSliceLabel));

  auto const& cluster_h = e.getValidHandle<std::vector<recob::Cluster> >(fPfpLabel);
  art::FindManyP<recob::Hit> cluster_hit_assn_v(cluster_h, e, fPfpLabel);

  art::Handle<std::vector<recob::Hit>> hitListHandle;
  e.getByLabel(fHitLabel, hitListHandle);
  std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> hittruth;
  if (!fHitTruthLabel.empty()) {
    hittruth = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
      hitListHandle, e, fHitTruthLabel);
  }

  //select the neutrino slice, or the one with largest NuScore if there is no neutino slice
  size_t primaryIdx = inputPfp->size();
  float maxNuScore = std::numeric_limits<float>::lowest();
  bool foundNuSlice = false;
  recob::Vertex::Point_t nuvtx;
  for (size_t ipfp = 0; ipfp < inputPfp->size(); ipfp++) {
    art::Ptr<recob::PFParticle> pfp(inputPfp, ipfp);
    if (pfp->IsPrimary() == false) continue;
    auto PDG = fabs(pfp->PdgCode());
    if (PDG == 12 || PDG == 14) {
      primaryIdx = ipfp;
      foundNuSlice = true;
      nuvtx = assocPfpVertex->at(pfp.key()).at(0)->position();
      break;
    }
    //
    if (fRecoverHighestNuScoreSlice==false) continue;
    //
    auto pfParticleMetadata = assocPfpMetadata->at(ipfp);
    auto pfParticlePropertiesMap = pfParticleMetadata.GetPropertiesMap();
    if (!pfParticlePropertiesMap.empty()) {
      auto it = pfParticlePropertiesMap.begin();
      while (it != pfParticlePropertiesMap.end()) {
	//std::cout << "ipfp=" << ipfp << " primary=" << pfp->IsPrimary() << " meta=" << it->first << " " << it->second << std::endl;
	if (it->first == "NuScore") {
	  if (pfParticlePropertiesMap.at(it->first)>maxNuScore) {
	    primaryIdx = ipfp;
	    maxNuScore = pfParticlePropertiesMap.at(it->first);
	    nuvtx = assocPfpVertex->at(pfp.key()).at(0)->position();
	  }
	}
	it++;
      }
    } // if PFP metadata exists!
  }

  //now save the hits for this slice
  if (primaryIdx < inputPfp->size()) {
    art::Ptr<recob::PFParticle> pfp(inputPfp, primaryIdx);

    auto assocSlice = assocPfpSlice->at(pfp.key());
    auto sliceHits = assocSliceHit->at(assocSlice[0].key());

    for (size_t ihit = 0; ihit < sliceHits.size(); ++ihit) {
      auto hit = sliceHits.at(ihit);
      outputHits->emplace_back(*hit);

      if (!hittruth) continue;
      std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(hit.key());
      const art::Ptr<recob::Hit> ahp = hitPtrMaker(outputHits->size() - 1);
      for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
        outputHitPartAssns->addSingle(particle_vec[i_p], ahp, *match_vec[i_p]);
      }
    }
  }

  std::map<unsigned int, unsigned int> pfpmap;
  for (unsigned int p=0; p < inputPfp->size(); p++) pfpmap[inputPfp->at(p).Self()] = p;

  //now try to recover hits from potential 2nd showers in other slices
  if (fRecover2ndShower) {
    for (size_t ipfp = 0; ipfp < inputPfp->size(); ipfp++) {
      art::Ptr<recob::PFParticle> pfp(inputPfp, ipfp);
      if (pfp->IsPrimary() == false) continue;
      if (ipfp == primaryIdx) continue;

      // skip clear cosmics
      bool isCC = false;
      auto pfParticleMetadata = assocPfpMetadata->at(ipfp);
      auto pfParticlePropertiesMap = pfParticleMetadata.GetPropertiesMap();
      for (auto it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); it++) {
	//std::cout << "ipfp=" << ipfp << " " << it->first << " " << it->second << std::endl;
	if (it->first == "IsClearCosmic") {
	  isCC = true;
	  break;
	}
      }
      if (isCC) continue;

      std::vector<art::Ptr<recob::PFParticle> > pfp_ptr_v;
      AddDaughters(pfp, inputPfp, pfp_ptr_v,pfpmap);

      for (unsigned int q=0; q < pfp_ptr_v.size(); q++) {

	// only pfps within 1m of the neutrino vertex
	auto pfvtx = assocPfpVertex->at(pfp_ptr_v[q].key()).at(0)->position();
	std::cout << "pfp vtx dist=" << std::sqrt( (nuvtx-pfvtx).Mag2() ) << std::endl;
	if ( (nuvtx-pfvtx).Mag2() > fVtxDistCut*fVtxDistCut ) continue;

	const std::vector< art::Ptr<recob::Cluster> > this_cluster_ptr_v = assocPfpCluster->at( pfp_ptr_v[q].key() );

	int nhits = 0;
	for (auto cluster_ptr : this_cluster_ptr_v) nhits += cluster_hit_assn_v.at( cluster_ptr.key() ).size();

	//consider only showers or tracks with a limited number of hits consistent with a pi0 shower
	std::cout << "pfp pdg=" << pfp_ptr_v[q]->PdgCode() << " nhits=" << nhits << std::endl;
	if (pfp_ptr_v[q]->PdgCode()==13 && nhits>fMaxHitCut) continue;

	std::cout << "adding hits from this pfp" << std::endl;
	for (auto cluster_ptr : this_cluster_ptr_v) {
	  const std::vector< art::Ptr<recob::Hit> > this_hit_ptr_v = cluster_hit_assn_v.at( cluster_ptr.key() );
	  for (auto hit : this_hit_ptr_v) {
	    outputHits->emplace_back(*hit);

	    if (!hittruth) continue;
	    std::vector<art::Ptr<simb::MCParticle>> particle_vec = hittruth->at(hit.key());
	    std::vector<anab::BackTrackerHitMatchingData const*> match_vec = hittruth->data(hit.key());
	    const art::Ptr<recob::Hit> ahp = hitPtrMaker(outputHits->size() - 1);
	    for (size_t i_p = 0; i_p < particle_vec.size(); ++i_p) {
	      outputHitPartAssns->addSingle(particle_vec[i_p], ahp, *match_vec[i_p]);
	    }
	  }
	}
      }
    }
  }

  std::cout << "NuSliceHitProducer nhits=" << outputHits->size() << " assns=" << outputHitPartAssns->size() << " foundNuSlice=" << foundNuSlice << std::endl;
  e.put(std::move(outputHits));
  if (!fHitTruthLabel.empty()) e.put(std::move(outputHitPartAssns));
}

void NuSliceHitsProducer::AddDaughters(const art::Ptr<recob::PFParticle>& pfp_ptr,  const art::ValidHandle<std::vector<recob::PFParticle> >& pfp_h, std::vector<art::Ptr<recob::PFParticle> > &pfp_v,std::map<unsigned int, unsigned int>& pfpmap) {

  auto daughters = pfp_ptr->Daughters();

  pfp_v.push_back(pfp_ptr);

  for(auto const& daughterid : daughters) {

    if (pfpmap.find(daughterid) == pfpmap.end()) {
      std::cout << "Did not find DAUGHTERID in map! error"<< std::endl;
      continue;
    }

    const art::Ptr<recob::PFParticle> pfp_ptr(pfp_h, pfpmap.at(daughterid) );

    AddDaughters(pfp_ptr, pfp_h, pfp_v, pfpmap);

  }// for all daughters

  return;
}


DEFINE_ART_MODULE(NuSliceHitsProducer)

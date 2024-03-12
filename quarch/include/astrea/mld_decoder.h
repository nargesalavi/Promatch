/*
 *  author: Suhas Vittal
 *  date:   30 March 2023
 * */

#ifndef ASTREA_MLD_DECODER_h
#define ASTREA_MLD_DECODER_h

#include "mwpm_decoder.h"

#include <array>

namespace qrc {
namespace astrea {

class MLDDecoder : public MWPMDecoder {
public:
    MLDDecoder(const stim::Circuit&);

    std::string name(void) override;
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    fp_t get_obs_prob(uint);
    std::map<uint, uint> get_representative(uint);
protected:
    void explore(const std::map<uint, uint>& base, fp_t weight);

    std::vector<uint> curr_detector_list;
    std::array<std::map<uint, uint>, 2> representatives;
    std::array<fp_t, 2> obs_prob;
};

}   // astrea
}   // qrc

#endif

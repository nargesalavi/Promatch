/*
 *  author: Suhas Vittal
 *  date:   1 December 2022
 * */

#ifndef TMR_DECODER_h
#define TMR_DECODER_h

#include "mwpm_decoder.h"
#include "defs.h"

#include <algorithm>

#define TMR_DEBUG

namespace qrc {

class TMRDecoder : public MWPMDecoder {
public:
    TMRDecoder(const stim::Circuit&, Decoder*, uint detectors_per_round);

    std::string name(void) override;
    bool is_software(void) override;
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
private:
    void unxor(std::vector<uint8_t>&);
    std::vector<uint8_t> majority(
            const std::vector<uint8_t>&, 
            const std::vector<uint8_t>&,
            const std::vector<uint8_t>&);

    Decoder * baseline;
    uint detectors_per_round;
};


}   // qrc

#endif

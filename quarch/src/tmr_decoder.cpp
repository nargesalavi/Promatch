/*
 *  author: Suhas Vittal
 *  date:   1 December 2022
 * */

#include "tmr_decoder.h"

namespace qrc {

TMRDecoder::TMRDecoder(
        const stim::Circuit& circuit,
        Decoder * baseline,
        uint detectors_per_round)
    :MWPMDecoder(circuit),
    baseline(baseline),
    detectors_per_round(detectors_per_round)
{}

std::string
TMRDecoder::name() {
    return (baseline->name() + "[TMR]");
}

bool
TMRDecoder::is_software() {
    return baseline->is_software();
}

DecoderShotResult
TMRDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();
    
    uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
    if (hw <= 10) {
        return MWPMDecoder::decode_error(syndrome);
    }

    uint tmr_detectors = (n_detectors - detectors_per_round) / 3;
    std::vector<uint8_t> tmr_syndrome(tmr_detectors
                                      + detectors_per_round
                                      + n_observables);

    std::vector<uint8_t> ssyndrome(syndrome);
    unxor(ssyndrome);
    // Convert unxor'd syndrome to TMR syndrome.
    std::vector<uint8_t> last_vector(detectors_per_round, 0);
    for (uint i = 0; i < tmr_detectors; i += detectors_per_round) {
        uint k1 = 3*i,
             k2 = 3*i + detectors_per_round,
             k3 = 3*i + 2*detectors_per_round,
             k4 = 3*i + 3*detectors_per_round;
        std::vector<uint8_t> r1(ssyndrome.begin() + k1, ssyndrome.begin() + k2);
        std::vector<uint8_t> r2(ssyndrome.begin() + k2, ssyndrome.begin() + k3);
        std::vector<uint8_t> r3(ssyndrome.begin() + k3, ssyndrome.begin() + k4);

        std::vector<uint8_t> res = majority(r1, r2, r3);
        for (uint j = 0; j < detectors_per_round; j++) {
            tmr_syndrome[i+j] = last_vector[j] ^ res[j];
            last_vector[j] = res[j];
        }
    }
    // Add last round as is.
    for (uint i = 0; i < detectors_per_round; i++) {
        tmr_syndrome[tmr_detectors + i] =
                ssyndrome[n_detectors - detectors_per_round + i] ^ last_vector[i];
    }
    // Also add observables (unchanged).
    for (uint i = 0; i < n_observables; i++) {
        tmr_syndrome[tmr_detectors + detectors_per_round + i] =
                        ssyndrome[n_detectors + i];
    }

    auto res = baseline->decode_error(tmr_syndrome);
#ifdef TMR_DEBUG
    if (res.is_logical_error) {
        uint rounds = n_detectors / detectors_per_round;
        for (uint r = 0; r < rounds; r++) {
            std::string orig, tmr;
            for (uint i = 0; i < detectors_per_round; i++) {
                uint k1 = r*detectors_per_round + i;
                if (syndrome[k1]) {
                    orig.push_back('!');
                } else {
                    orig.push_back('.');
                }

                uint k2 = (r/3)*detectors_per_round + i;
                if (tmr_syndrome[k2]) {
                    tmr.push_back('!');
                } else {
                    tmr.push_back('.');
                }
            }
            std::cout << orig;
            if (r % 3 == 0) {
                std::cout << "\t|\t" << tmr;
            }
            std::cout << "\n";
        }
        std::cout << "=======================\n";
    }
#endif
    return res;
}

void
TMRDecoder::unxor(std::vector<uint8_t>& syndrome) {
    uint n_observables = circuit.count_observables();
    for (uint i = detectors_per_round; i < syndrome.size() - n_observables; i++) 
    {
        syndrome[i] ^= syndrome[i-detectors_per_round];
    }
}

std::vector<uint8_t> 
TMRDecoder::majority(
        const std::vector<uint8_t>& r1,
        const std::vector<uint8_t>& r2,
        const std::vector<uint8_t>& r3)
{
    std::vector<uint8_t> maj(r1.size());
    for (uint i = 0; i < maj.size(); i++) {
        maj[i] = (r1[i] + r2[i] + r3[i]) >= 2;
    }
    return maj;
}

}   // qrc

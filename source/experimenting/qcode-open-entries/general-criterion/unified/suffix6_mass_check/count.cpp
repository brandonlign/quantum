// Count all one-dimensional extensions of binary additive [9,5,6] codes.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

int main() {
    std::vector<unsigned> low;
    for (unsigned w = 0; w < (1u << 18); ++w)
        if (__builtin_popcount((w & 511u) | (w >> 9)) < 6) low.push_back(w);
    unsigned count;
    if (!(std::cin >> count)) return 1;
    for (unsigned row = 0; row < count; ++row) {
        unsigned id, g;
        if (!(std::cin >> id)) return 2;
        std::array<unsigned, 18> pivots{};
        unsigned rank = 0;
        for (unsigned i = 0; i < 5; ++i) {
            if (!(std::cin >> g) || g >= (1u << 18)) return 3;
            for (int b = 17; b >= 0; --b) {
                if (!(g & (1u << b))) continue;
                if (pivots[b]) g ^= pivots[b];
                else { pivots[b] = g; ++rank; break; }
            }
        }
        if (rank != 5) return 4;
        std::array<unsigned, 18> label{};
        for (unsigned i = 0; i < 18; ++i) {
            unsigned v = 1u << i;
            for (int b = 17; b >= 0; --b)
                if (pivots[b] && (v & (1u << b))) v ^= pivots[b];
            unsigned k = 0;
            for (unsigned b = 0; b < 18; ++b)
                if (!pivots[b]) label[i] |= ((v >> b) & 1u) << k++;
            if (k != 13) return 5;
        }
        std::array<unsigned, 512> lo{}, hi{};
        for (unsigned w = 1; w < 512; ++w) {
            unsigned b = __builtin_ctz(w), previous = w & (w - 1);
            lo[w] = lo[previous] ^ label[b];
            hi[w] = hi[previous] ^ label[b + 9];
        }
        std::array<unsigned char, 8192> reached{};
        for (unsigned w : low) reached[lo[w & 511u] ^ hi[w >> 9]] = 1;
        unsigned good = 0;
        uint64_t checksum = 14695981039346656037ULL;
        for (unsigned i = 0; i < 8192; ++i) {
            if (!reached[i]) ++good;
            checksum = (checksum ^ reached[i]) * 1099511628211ULL;
        }
        std::cout << id << ' ' << good << ' ' << checksum << '\n';
    }
    return std::cin.eof() || std::cin.good() ? 0 : 6;
}

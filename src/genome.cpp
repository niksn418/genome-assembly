#include "genome.h"

#include <algorithm>
#include <array>
#include <deque>
#include <unordered_map>
#include <vector>

namespace genome {

namespace {

struct genome_string_hasher
{
    static constexpr size_t p = 5;

    const size_t k;
    std::array<size_t, 4> max_hash_values;
    size_t last_hash;

    genome_string_hasher(const size_t k)
        : k(k)
    {
        const size_t max_power = pow(p, k);
        max_hash_values[0] = 0;
        max_hash_values[1] = max_power;
        for (size_t i = 2; i < 4; ++i) {
            max_hash_values[i] = max_hash_values[i - 1] + max_power;
        }
    }

    static constexpr size_t encode_char(const char c)
    {
        return (c & 6) >> 1;
    }

    size_t hash(const std::string_view read, const size_t pos)
    {
        return last_hash = pos != 0
                ? from_previous_hash(last_hash, read[pos - 1], read[k + pos - 1])
                : hash(read.substr(0, k));
    }

    static constexpr size_t pow(size_t x, size_t power)
    {
        size_t result = 1;
        while (power) {
            if (power & 1) {
                result *= x;
            }
            x *= x;
            power >>= 1;
        }
        return result;
    }

    static constexpr size_t hash(const std::string_view & sv)
    {
        size_t hash = 0;
        for (const char c : sv) {
            hash = hash * p + encode_char(c);
        }
        return hash;
    }

    size_t from_previous_hash(const size_t previous_hash, const char last_symbol, const char new_symbol) const
    {
        return previous_hash * p - max_hash_values[encode_char(last_symbol)] + encode_char(new_symbol);
    }
};

struct genome_string
{
    const std::string_view m_value;
    const size_t m_hash;

    bool operator==(const genome_string & other) const
    {
        return m_value == other.m_value;
    }

    struct hash
    {
        size_t operator()(const genome_string & s) const
        {
            return s.m_hash;
        }
    };
};

struct genome_graph
{
    struct vertex
    {
        std::vector<vertex *> neighbours;
        const char value;
        std::size_t in_degree = 0;

        vertex(char value)
            : value(value)
        {
        }
    };

    std::unordered_map<genome_string, vertex, genome_string::hash> vertices;

    void add_edge(vertex * const u, vertex * const v)
    {
        u->neighbours.push_back(v);
        ++v->in_degree;
    }

    vertex * add_vertex(genome_string && key)
    {
        auto res = vertices.emplace(key, key.m_value.back());
        return &res.first->second;
    }
};

} // anonymous namespace

std::string assembly(const size_t k, const std::vector<std::string> & reads)
{
    if (k == 0 || reads.empty()) {
        return "";
    }
    const size_t d = reads[0].size();
    genome_graph graph;
    genome_string_hasher hasher(k);
    using vertex = genome_graph::vertex;
    graph.vertices.reserve(reads.size() * (d - k) + 1);
    for (const std::string_view read : reads) {
        vertex * u = graph.add_vertex({read.substr(0, k), hasher.hash(read, 0)});
        for (size_t i = 1; i <= d - k; ++i) {
            vertex * const v = graph.add_vertex({read.substr(i, k), hasher.hash(read, i)});
            graph.add_edge(u, v);
            u = v;
        }
    }
    std::deque<vertex *> s;
    std::string result;
    result.resize(k + reads.size() * (d - k));
    for (auto it = graph.vertices.begin(), end = graph.vertices.end(); it != end; ++it) {
        auto & [key, v] = *it;
        if (v.neighbours.size() > v.in_degree) {
            s.push_front(&v);
            std::copy(key.m_value.begin(), key.m_value.end(), result.begin());
            break;
        }
    }
    size_t index = result.size();
    while (!s.empty()) {
        vertex & u = *s.front();
        if (u.neighbours.empty()) {
            result[--index] = u.value;
            s.pop_front();
        }
        else {
            s.push_front(u.neighbours.back());
            u.neighbours.pop_back();
        }
    }
    return result;
}

} // namespace genome

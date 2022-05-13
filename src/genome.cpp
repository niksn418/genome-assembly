#include "genome.h"

#include <algorithm>
#include <deque>
#include <unordered_map>
#include <vector>

namespace genome {

namespace {

struct vertex
{
    using edge = std::pair<vertex * const, const std::string_view>;

    size_t in_degree;
    std::vector<edge> edges;
};

using string_graph = std::unordered_map<std::string_view, vertex>;

void add_edge(vertex * const u, vertex * const v, const std::string_view edge_value)
{
    u->edges.emplace_back(v, edge_value);
    ++(v->in_degree);
}

vertex * add_vertex(const std::string_view key, string_graph & graph)
{
    const auto res = graph.try_emplace(key);
    return &res.first->second;
}

} // anonymous namespace

std::string assembly(const size_t k, const std::vector<std::string> & reads)
{
    if (k == 0 || reads.empty()) {
        return "";
    }
    const size_t d = reads[0].size();
    string_graph graph;
    graph.reserve(reads.size() + 1);
    for (const std::string_view read : reads) {
        const auto u = add_vertex(read.substr(0, k), graph);
        const auto v = add_vertex(read.substr(d - k, k), graph);
        add_edge(u, v, read);
    }
    const auto start_it = std::find_if(graph.begin(), graph.end(), [](const auto & it) {
        return it.second.edges.size() > it.second.in_degree;
    });
    std::string result;
    result.resize(k + reads.size() * (d - k));
    std::copy(start_it->first.cbegin(), start_it->first.cend(), result.begin());
    auto result_it = result.end();
    std::deque<vertex::edge> stack;
    stack.emplace_front(&start_it->second, start_it->first);
    while (!stack.empty()) {
        const auto & [v, edge_value] = stack.front();
        std::vector<vertex::edge> & edges = v->edges;
        if (edges.empty()) {
            result_it = std::copy_backward(edge_value.cbegin() + k, edge_value.cend(), result_it);
            stack.pop_front();
        }
        else {
            stack.push_front(edges.back());
            edges.pop_back();
        }
    }
    return result;
}

} // namespace genome

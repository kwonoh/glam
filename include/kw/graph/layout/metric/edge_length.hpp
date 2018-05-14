// Copyright 2017 Oh-Hyun Kwon (kwonoh.net)
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef __KW_GRAPH_LAYOUT_METRIC_EDGE_LENGTH_HPP__
#define __KW_GRAPH_LAYOUT_METRIC_EDGE_LENGTH_HPP__

#include <cmath>

#include <tbb/tbb.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <boost/foreach.hpp>
#include <boost/geometry.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/range/combine.hpp>

#include <kw/algorithm/normalize.hpp>
#include <kw/iterator_range.hpp>

namespace kw {

// The coefficient of variance of the edge length
template <typename Graph, typename PositionMap>
double
edge_length_cv(Graph const& g, PositionMap pos)
{
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_desc;
    auto const num_edges = boost::num_edges(g);
    std::vector<double> lengths(num_edges);

    tbb::parallel_for_each(
        boost::combine(kw::make_iterator_range(boost::edges(g)), lengths),
        [&](auto itr) {
            auto const e = boost::get<0>(itr);
            auto& l = boost::get<1>(itr);

            auto const p = boost::get(pos, boost::source(e, g));
            auto const q = boost::get(pos, boost::target(e, g));
            l = boost::geometry::distance(p, q);
        });

    bool const same = !kw::minmax_normalize(lengths);
    if (same) return 0.0;

    auto const sum = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, num_edges), double(0),
        [&](tbb::blocked_range<size_t> const& r, double s) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                s += lengths[i];
            }
            return s;
        },
        [](double const a, double const b) { return a + b; });

    auto const mean = sum / num_edges;

    std::vector<double> sq_diffs(num_edges);
    tbb::parallel_for(std::size_t(0), num_edges, std::size_t(1),
                      [&](std::size_t const i) {
                          sq_diffs[i] = boost::math::pow<2>(lengths[i] - mean);
                      });

    auto const sq_diff_sum = tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_edges), double(0),
        [&](tbb::blocked_range<std::size_t> const& r, double s) {
            for (std::size_t i = r.begin(); i != r.end(); ++i) {
                s += sq_diffs[i];
            }
            return s;
        },
        [](double const a, double const b) { return a + b; });

    auto const stddev = std::sqrt(sq_diff_sum / num_edges);

    return stddev / mean;
}

// The normalized coefficient of variance of the edge length
// The upper bound of the coefficient of variation of n values is sqrt(n - 1).
template <typename Graph, typename PositionMap>
double
edge_length_normalized_cv(Graph const& g, PositionMap pos, double& cv)
{
    cv = edge_length_cv(g, pos);
    return cv / std::sqrt(boost::num_edges(g) - 1);
}

template <typename Graph, typename PositionMap>
double
edge_length_normalized_cv(Graph const& g, PositionMap pos)
{
    double cv;
    return edge_length_normalized_cv(g, pos, cv);
}

}  // namespace kw

#endif  // __KW_GRAPH_LAYOUT_METRIC_EDGE_LENGTH_HPP__

// Copyright 2017 Oh-Hyun Kwon (kwonoh.net)
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef VIDI_GRAPH_LAYOUT_METRIC_SHAPE_BASED_HPP
#define VIDI_GRAPH_LAYOUT_METRIC_SHAPE_BASED_HPP

#include <cmath>
#include <iostream>

#include <tbb/tbb.h>
#include <boost/foreach.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/graph/adjacency_list.hpp>

#ifdef KW_USE_CGAL
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#endif

namespace kw {

#ifdef KW_USE_CGAL
// template <typename Points, typename DelaunayGraph>
// void
// delaunay_triangulation(Points const& points, DelaunayGraph& dg)
// {
//     namespace bg = boost::geometry;
//
//     typedef typename boost::graph_traits<DelaunayGraph>::vertex_descriptor
//         vertex_descriptor;
//
//     typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//     typedef CGAL::Triangulation_vertex_base_with_info_2<vertex_descriptor, K>
//         Vb;
//     typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
//     typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
//     typedef typename Delaunay::Point DPoint;
//
//     Delaunay dt;
//     for (auto const& p : points) {
//         dt.insert(DPoint(bg::get<0>(p), bg::get<1>(p)))->info() =
//             boost::add_vertex(dg);
//     }
//
//     for (auto e_itr = dt.finite_edges_begin(), e_end = dt.finite_edges_end();
//          e_itr != e_end; ++e_itr) {
//         auto const f = e_itr->first;
//         auto const i = e_itr->second;
//         vertex_descriptor const s = f->vertex(dt.cw(i))->info();
//         vertex_descriptor const t = f->vertex(dt.ccw(i))->info();
//         boost::add_edge(s, t, dg);
//     }
// }
template <typename DelaunayGraph, typename PositionMap>
void
delaunay_triangulation(DelaunayGraph& g, PositionMap pos)
{
    namespace bg = boost::geometry;

    typedef typename boost::graph_traits<DelaunayGraph>::vertex_descriptor
        vertex_descriptor;

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<vertex_descriptor, K>
        Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
    typedef typename Delaunay::Point DPoint;

    Delaunay dt;
    BOOST_FOREACH (auto u, boost::vertices(g)) {
        auto const& p = boost::get(pos, u);
        dt.insert(DPoint(bg::get<0>(p), bg::get<1>(p)))->info() = u;
    }

    for (auto e_itr = dt.finite_edges_begin(), e_end = dt.finite_edges_end();
         e_itr != e_end; ++e_itr) {
        auto const f = e_itr->first;
        auto const i = e_itr->second;
        vertex_descriptor const s = f->vertex(dt.cw(i))->info();
        vertex_descriptor const t = f->vertex(dt.ccw(i))->info();
        boost::add_edge(s, t, g);
    }
}
#endif

template <typename DelaunayGraph, typename PositionMap>
void
convert_delaunay_triangulation_to_gabriel_graph(DelaunayGraph& g,
                                                PositionMap vp_pos)
{
    namespace bg = boost::geometry;
    namespace bgi = boost::geometry::index;

    typedef typename boost::graph_traits<DelaunayGraph>::vertices_size_type
        vertices_size_type;
    typedef
        typename boost::graph_traits<DelaunayGraph>::edge_descriptor edge_desc;

    typedef typename boost::property_traits<PositionMap>::value_type point_t;
    typedef std::pair<point_t, vertices_size_type> indexed_point_t;

    auto const vp_idx = boost::get(boost::vertex_index, g);
    std::vector<indexed_point_t> indexed_points(boost::num_vertices(g));
    BOOST_FOREACH (auto u, boost::vertices(g)) {
        auto const v_idx = vp_idx[u];
        indexed_points[v_idx] = std::make_pair(vp_pos[u], v_idx);
    }

    bgi::rtree<indexed_point_t, bgi::rstar<8>> rtree(indexed_points.begin(),
                                                     indexed_points.end());

    std::vector<edge_desc> not_gabriel_edges;
    BOOST_FOREACH (auto e, boost::edges(g)) {
        auto const u = boost::source(e, g);
        auto const v = boost::target(e, g);
        auto const p = vp_pos[u];
        auto const q = vp_pos[v];

        auto center = p;
        bg::add_point(center, q);
        bg::divide_value(center, 2);

        for (auto it = rtree.qbegin(bgi::nearest(center, 1));
             it != rtree.qend(); ++it) {
            if (it->second != vp_idx[u] && it->second != vp_idx[v]) {
                not_gabriel_edges.push_back(e);
            }
        }
    }

    for (auto e : not_gabriel_edges) {
        boost::remove_edge(e, g);
    }
}

template <typename Graph, typename ShapeGraph>
double
shape_based_metric(Graph const& graph, ShapeGraph const& shape)
{
    typedef typename boost::graph_traits<Graph>::vertex_descriptor g_v_desc;
    typedef typename boost::graph_traits<Graph>::vertices_size_type
        vertices_size_type;
    typedef
        typename boost::graph_traits<ShapeGraph>::vertex_descriptor s_v_desc;

    auto g_vp_idx = boost::get(boost::vertex_index, graph);
    auto s_vp_idx = boost::get(boost::vertex_index, shape);

    auto s_v_itrs = boost::vertices(shape);
    std::vector<s_v_desc> s_v_descs(s_v_itrs.first, s_v_itrs.second);

    std::vector<double> jaccard_similarities(boost::num_vertices(graph));

    tbb::task_group task_group;
    BOOST_FOREACH (auto g_u, boost::vertices(graph)) {
        task_group.run([&, g_u] {
            auto const g_u_idx = g_vp_idx[g_u];
            auto const s_u = s_v_descs[g_u_idx];
            assert(s_vp_idx[s_u] == g_u_idx);

            std::vector<vertices_size_type> g_adjs;
            BOOST_FOREACH (auto g_v, boost::adjacent_vertices(g_u, graph)) {
                g_adjs.push_back(g_vp_idx[g_v]);
            }
            std::sort(g_adjs.begin(), g_adjs.end());

            std::vector<vertices_size_type> s_adjs;
            BOOST_FOREACH (auto s_v, boost::adjacent_vertices(s_u, shape)) {
                s_adjs.push_back(s_vp_idx[s_v]);
            }
            std::sort(s_adjs.begin(), s_adjs.end());

            std::vector<vertices_size_type> adjs_union;
            std::set_union(g_adjs.begin(), g_adjs.end(), s_adjs.begin(),
                           s_adjs.end(), std::back_inserter(adjs_union));

            std::vector<vertices_size_type> adjs_intersection;
            std::set_intersection(g_adjs.begin(), g_adjs.end(), s_adjs.begin(),
                                  s_adjs.end(),
                                  std::back_inserter(adjs_intersection));

            jaccard_similarities[g_u_idx] =
                double(adjs_intersection.size()) / double(adjs_union.size());
        });
    }
    task_group.wait();

    return std::accumulate(jaccard_similarities.begin(),
                           jaccard_similarities.end(), 0.0) /
           boost::num_vertices(graph);
}

}  // namespace kw

#endif  // VIDI_GRAPH_LAYOUT_METRIC_SHAPE_BASED_HPP

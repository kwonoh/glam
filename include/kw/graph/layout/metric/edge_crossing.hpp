// Copyright 2017 Oh-Hyun Kwon (kwonoh.net)
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef __KW_GRAPH_LAYOUT_METRIC_EDGE_CROSSING_HPP__
#define __KW_GRAPH_LAYOUT_METRIC_EDGE_CROSSING_HPP__

#include <cmath>
#include <mutex>

#include <tbb/tbb.h>
#include <boost/compute.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/make_shared.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace kw {

namespace detail {

struct edge_crossings_kernel {
    typedef boost::compute::float_ d_vertex_pos_coord_t;
    typedef boost::compute::float2_ d_vertex_pos_t;
    typedef boost::compute::ulong_ d_vertex_idx_t;
    typedef boost::compute::ulong2_ d_edge_t;
    typedef boost::compute::ulong_ d_edge_size_t;

    boost::compute::device device;
    boost::compute::context context;
    boost::compute::command_queue queue;
    boost::compute::program program;
    boost::compute::kernel kernel;

    boost::compute::ulong2_ local_work_size;
    boost::compute::ulong2_ global_work_size;

    boost::shared_ptr<boost::compute::vector<d_vertex_pos_t>> d_pos;
    boost::shared_ptr<boost::compute::vector<d_edge_t>> d_edges;
    boost::shared_ptr<boost::compute::vector<d_edge_size_t>> d_crossing_count;

    d_edge_size_t n_crossings;

    edge_crossings_kernel() = default;

    edge_crossings_kernel(boost::compute::device const& d)
        : device(d),
          context(device),
          queue(context, device),
          program(),
          kernel(),
          local_work_size(16, 16),
          global_work_size(),
          d_pos(),
          d_edges(),
          d_crossing_count()
    {
        build_kernel();
    }

    void
    build_kernel()
    {
        std::stringstream program_source_ss;
        program_source_ss
            << "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n";
        program_source_ss << BOOST_COMPUTE_STRINGIZE_SOURCE(

            inline int orient(float2 p, float2 q, float2 r) {
                const float2 a = q - p;
                const float2 b = r - q;
                const float o = a.y * b.x - a.x * b.y;
                if (o == 0.0) {
                    return 0;
                }
                return copysign(1, o);
            }

            inline bool on_seg(float2 p, float2 q, float2 r) {
                return all(min(p, r) <= q && q <= max(p, r));
            }

            __kernel void count_edge_crossings(
                __global const float2* pos,     // 0
                __global const ulong2* edges,   // 1
                const ulong n_edges,            // 2
                const ulong work_offset_i,      // 3
                const ulong work_offset_j,      // 4
                const ulong work_size_1d,       // 5
                __global ulong* crossing_count  // 6
            ) {
                const size_t global_i = get_global_id(0) + work_offset_i;
                const size_t global_j = get_global_id(1) + work_offset_j;

                if (global_j <= global_i ||  //
                    global_i >= n_edges || global_j >= n_edges ||
                    global_i >= work_offset_i + work_size_1d ||
                    global_j >= work_offset_j + work_size_1d) {
                    return;
                }

                const ulong2 e_i = edges[global_i];
                const ulong2 e_j = edges[global_j];

                const ulong u_i = e_i.s0;
                const ulong v_i = e_i.s1;
                const ulong u_j = e_j.s0;
                const ulong v_j = e_j.s1;

                if (u_i == u_j || u_i == v_j || v_i == u_j || v_i == v_j) {
                    return;
                }

                const float2 p_i = pos[u_i];
                const float2 q_i = pos[v_i];
                const float2 p_j = pos[u_j];
                const float2 q_j = pos[v_j];

                const int o_0 = orient(p_i, q_i, p_j);
                const int o_1 = orient(p_i, q_i, q_j);
                const int o_2 = orient(p_j, q_j, p_i);
                const int o_3 = orient(p_j, q_j, q_i);

                if (o_0 != o_1 && o_2 != o_3) {
                    crossing_count[get_global_id(1) * get_global_size(0) +
                                   get_global_id(0)] += 1;
                }
                else if ((o_0 == 0 && on_seg(p_i, p_j, q_i)) ||
                         (o_1 == 0 && on_seg(p_i, q_j, q_i)) ||
                         (o_2 == 0 && on_seg(p_j, p_i, q_j)) ||
                         (o_3 == 0 && on_seg(p_j, q_i, q_j))) {
                    crossing_count[get_global_id(1) * get_global_size(0) +
                                   get_global_id(0)] += 1;
                }
            }

        );

        program = boost::compute::program::build_with_source(
            program_source_ss.str(), context);
        kernel = boost::compute::kernel(program, "count_edge_crossings");
    }

    template <typename Graph, typename PositionMap>
    void
    setup(Graph const& g, PositionMap pos, d_edge_size_t const work_size_1d)
    {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor
            vertex_descriptor;
        typedef typename boost::graph_traits<Graph>::edge_descriptor
            edge_descriptor;
        typedef typename boost::graph_traits<Graph>::edges_size_type
            edges_size_type;

        auto const n_vertices = boost::num_vertices(g);
        d_edge_size_t const n_edges = boost::num_edges(g);
        n_crossings = 0;

        auto v_idx = boost::get(boost::vertex_index, g);
        auto v_itr = boost::vertices(g);
        auto e_itr = boost::edges(g);

        d_pos = boost::make_shared<boost::compute::vector<d_vertex_pos_t>>(
            n_vertices, context);
        d_edges = boost::make_shared<boost::compute::vector<d_edge_t>>(n_edges,
                                                                       context);

        d_crossing_count =
            boost::make_shared<boost::compute::vector<d_edge_size_t>>(
                work_size_1d * work_size_1d, context);
        boost::compute::fill(d_crossing_count->begin(), d_crossing_count->end(),
                             d_edge_size_t(0), queue);

        {
            // TODO: is it possible to bypass h_pos?
            std::vector<d_vertex_pos_t> h_pos(n_vertices);
            std::transform(v_itr.first, v_itr.second, h_pos.begin(),
                           [&](vertex_descriptor u) {
                               auto const& p = boost::get(pos, u);
                               return d_vertex_pos_t(
                                   boost::numeric_cast<d_vertex_pos_coord_t>(
                                       boost::geometry::get<0>(p)),
                                   boost::numeric_cast<d_vertex_pos_coord_t>(
                                       boost::geometry::get<1>(p)));
                           });
            boost::compute::copy(h_pos.begin(), h_pos.end(), d_pos->begin(),
                                 queue);
        }

        {
            // TODO: is it possible to bypass h_edges?
            std::vector<d_edge_t> h_edges(n_edges);
            std::transform(e_itr.first, e_itr.second, h_edges.begin(),
                           [&](edge_descriptor e) {
                               return d_edge_t(
                                   boost::numeric_cast<d_vertex_idx_t>(
                                       boost::get(v_idx, boost::source(e, g))),
                                   boost::numeric_cast<d_vertex_idx_t>(
                                       boost::get(v_idx, boost::target(e, g))));
                           });
            boost::compute::copy(h_edges.begin(), h_edges.end(),
                                 d_edges->begin(), queue);
        }

        kernel.set_arg(0, *d_pos);
        kernel.set_arg(1, *d_edges);
        kernel.set_arg(2, n_edges);
        // kernel.set_arg(3, work_offset_i);
        // kernel.set_arg(4, work_offset_j);
        kernel.set_arg(5, work_size_1d);
        kernel.set_arg(6, *d_crossing_count);

        global_work_size[0] = work_size_1d;
        global_work_size[1] = work_size_1d;
    }

    void
    compute(d_edge_size_t const work_offset_i,
            d_edge_size_t const work_offset_j)
    {
        kernel.set_arg(3, work_offset_i);
        kernel.set_arg(4, work_offset_j);

        queue.enqueue_nd_range_kernel(
            kernel, boost::compute::dim(0, 0),
            boost::compute::dim(global_work_size[0], global_work_size[1]),
            boost::compute::dim(local_work_size[0], local_work_size[1]));
    }

    void
    finish()
    {
        boost::compute::reduce(d_crossing_count->begin(),
                               d_crossing_count->end(), &n_crossings,
                               boost::compute::plus<d_edge_size_t>(), queue);
    }
};

}  // namespace detail

template <typename Graph, typename PositionMap>
typename boost::graph_traits<Graph>::edges_size_type
num_edge_crossings(Graph const& g, PositionMap pos)
{
    static std::map<cl_device_id, detail::edge_crossings_kernel> kernels;
    std::vector<cl_device_id> devices;
    for (auto const& platform : boost::compute::system::platforms()) {
        for (auto& device : platform.devices(CL_DEVICE_TYPE_GPU |
                                             CL_DEVICE_TYPE_ACCELERATOR)) {
            devices.push_back(device.id());
            auto kernel_itr = kernels.find(device.id());
            if (kernel_itr == kernels.end()) {
                kernels.emplace(device.id(),
                                detail::edge_crossings_kernel(device));
            }
        }
    }

    std::size_t const n_edges = boost::num_edges(g);

    std::size_t const work_size_1d =
        std::min(((n_edges / 16) + 1) * 16, std::size_t(4096));
    std::size_t const n_groups_1d = (n_edges / work_size_1d) + 1;

    std::vector<std::pair<std::size_t, std::size_t>> groups_2d;
    for (std::size_t group_i = 0; group_i < n_groups_1d; group_i++) {
        for (std::size_t group_j = group_i; group_j < n_groups_1d; group_j++) {
            groups_2d.emplace_back(group_i, group_j);
        }
    }

    std::size_t const n_devices = devices.size();
    std::uint32_t current_group_idx = 0;
    std::mutex group_idx_mutex;

    tbb::task_group task_group;
    for (std::size_t device_idx = 0; device_idx < n_devices; device_idx++) {
        task_group.run([&, device_idx] {
            auto& kernel = kernels[devices[device_idx]];
            kernel.setup(g, pos, work_size_1d);

            while (current_group_idx < groups_2d.size()) {
                group_idx_mutex.lock();
                std::uint32_t const group_idx = current_group_idx++;
                group_idx_mutex.unlock();

                auto const& group = groups_2d[group_idx];
                kernel.compute(group.first * work_size_1d,
                               group.second * work_size_1d);
            }
        });
    }
    task_group.wait();

    typename boost::graph_traits<Graph>::edges_size_type n_crossings = 0;
    for (auto& device : devices) {
        auto& kernel = kernels[device];
        kernel.finish();
        n_crossings += kernel.n_crossings;
    }

    return n_crossings;
}

namespace detail {

template <typename Graph>
double
crosslessness_impl(
    Graph const& g,
    typename boost::graph_traits<Graph>::edges_size_type const n_crossings)
{
    typedef
        typename boost::graph_traits<Graph>::edges_size_type edges_size_type;

    edges_size_type const n_edges = boost::num_edges(g);
    edges_size_type const n_all = n_edges * (n_edges - 1) / 2;
    edges_size_type n_impossible = 0;
    BOOST_FOREACH (auto u, boost::vertices(g)) {
        auto const d = boost::degree(u, g);
        n_impossible += d * (d - 1);
    }
    n_impossible /= 2;
    auto const n_max = n_all - n_impossible;
    return n_max > 0 ? 1.0 - (double(n_crossings) / n_max) : 1.0;
}

}  // namespace detail

template <typename Graph, typename PositionMap>
double
crosslessness(Graph const& g,
              PositionMap pos,
              typename boost::graph_traits<Graph>::edges_size_type& n_crossings)
{
    n_crossings = num_edge_crossings(g, pos);
    return detail::crosslessness_impl(g, n_crossings);
}

template <typename Graph, typename PositionMap>
double
crosslessness(Graph const& g, PositionMap pos)
{
    return detail::crosslessness_impl(g, num_edge_crossings(g, pos));
}

}  // namespace kw

#endif  // __KW_GRAPH_LAYOUT_METRIC_EDGE_CROSSING_HPP__

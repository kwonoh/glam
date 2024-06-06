// Copyright 2017 Oh-Hyun Kwon (kwonoh.net)
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef KW_GRAPH_LAYOUT_METRIC_EDGE_CROSSING_HPP
#define KW_GRAPH_LAYOUT_METRIC_EDGE_CROSSING_HPP

#include <tbb/tbb.h>

#include <algorithm>
#include <atomic>
#include <boost/compute.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace kw {

namespace detail {

boost::compute::program
make_edge_crossings_program(boost::compute::context const& context) {
  std::stringstream kernel_source;

  // We can't embed a #pragma directive within macro arguments
  kernel_source
    << "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n";

  // OpenCL kernel code for counting edge crossings based on line segment
  // intersection
  kernel_source << BOOST_COMPUTE_STRINGIZE_SOURCE(
    inline int orient(const float2 p, const float2 q, const float2 r) {
      const float2 a = q - p;
      const float2 b = r - q;
      return convert_int(sign(a.y * b.x - a.x * b.y));
    }

    inline bool on_seg(const float2 p, const float2 q, const float2 r) {
      return all(min(p, r) <= q && q <= max(p, r));
    }

    __kernel void count_edge_crossings(
      __global const float2* pos,     // 0
      __global const ulong2* edges,   // 1
      const ulong n_edges,            // 2
      const ulong work_offset_i,      // 3
      const ulong work_offset_j,      // 4
      const ulong work_size_1d,       // 5
      __global ulong* crossing_count  // 6 (Output)
    ) {
      const size_t global_i = get_global_id(0);
      const size_t global_j = get_global_id(1);

      const size_t edge_i = global_i + work_offset_i;
      const size_t edge_j = global_j + work_offset_j;

      // Skip if the current edge indices are out of range
      // or in one triangular section of the pairwise matrix
      if (edge_j <= edge_i || edge_i >= n_edges || edge_j >= n_edges ||
          edge_i >= work_offset_i + work_size_1d ||
          edge_j >= work_offset_j + work_size_1d) {
        return;
      }

      const ulong2 e_i = edges[edge_i];
      const ulong2 e_j = edges[edge_j];

      const ulong u_i = e_i.s0;
      const ulong v_i = e_i.s1;
      const ulong u_j = e_j.s0;
      const ulong v_j = e_j.s1;

      // Skip if the two edges share a vertex
      if (u_i == u_j || u_i == v_j || v_i == u_j || v_i == v_j) {
        return;
      }

      //         p_j
      //         /
      // p_i ---/--- q_i (e_i)
      //       /
      //     q_j (e_j)

      const float2 p_i = pos[u_i];
      const float2 q_i = pos[v_i];
      const float2 p_j = pos[u_j];
      const float2 q_j = pos[v_j];

      const int o_0 = orient(p_i, q_i, p_j);
      const int o_1 = orient(p_i, q_i, q_j);
      const int o_2 = orient(p_j, q_j, p_i);
      const int o_3 = orient(p_j, q_j, q_i);

      const size_t output_idx = global_i * get_global_size(0) + global_j;

      // Two edges are crossed
      if (o_0 != o_1 && o_2 != o_3) {
        crossing_count[output_idx] += 1;
      }
      // One end of an edge is on the other edge
      else if ((o_0 == 0 && on_seg(p_i, p_j, q_i)) ||
               (o_1 == 0 && on_seg(p_i, q_j, q_i)) ||
               (o_2 == 0 && on_seg(p_j, p_i, q_j)) ||
               (o_3 == 0 && on_seg(p_j, q_i, q_j))) {
        crossing_count[output_idx] += 1;
      }
    }

  );

  return boost::compute::program::build_with_source(kernel_source.str(),
                                                    context);
}

class edge_crossings_kernel {
protected:
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

  boost::compute::vector<d_vertex_pos_t> d_pos;
  boost::compute::vector<d_edge_t> d_edges;
  boost::compute::vector<d_edge_size_t> d_crossing_count;

public:
  d_edge_size_t n_crossings;

  edge_crossings_kernel() = delete;

  edge_crossings_kernel(boost::compute::device const& d)
    : device(d),
      context(device),
      queue(context, device),
      program(make_edge_crossings_program(context)),
      kernel(boost::compute::kernel(program, "count_edge_crossings")),
      local_work_size(),
      global_work_size(),
      d_pos(),
      d_edges(),
      d_crossing_count() {
  }

protected:
  template <typename Graph, typename PositionMap>
  void
  setup_device_vertex_pos(Graph const& g, PositionMap pos) {
    d_pos =
      boost::compute::vector<d_vertex_pos_t>(boost::num_vertices(g), context);

    auto make_device_vertex_pos =
      [&](typename boost::graph_traits<Graph>::vertex_descriptor const u) {
        auto const& p = boost::get(pos, u);
        auto const x = boost::geometry::get<0>(p);
        auto const y = boost::geometry::get<1>(p);

        return d_vertex_pos_t(boost::numeric_cast<d_vertex_pos_coord_t>(x),
                              boost::numeric_cast<d_vertex_pos_coord_t>(y));
      };

    auto v_itr = boost::vertices(g);
    boost::compute::copy(
      boost::make_transform_iterator(v_itr.first, make_device_vertex_pos),
      boost::make_transform_iterator(v_itr.second, make_device_vertex_pos),
      d_pos.begin(), queue);
  }

  template <typename Graph>
  void
  setup_device_edges(Graph const& g) {
    auto v_idx = boost::get(boost::vertex_index, g);

    d_edges = boost::compute::vector<d_edge_t>(boost::num_edges(g), context);

    auto make_device_edge =
      [&](typename boost::graph_traits<Graph>::edge_descriptor const e) {
        auto const u = boost::get(v_idx, boost::source(e, g));
        auto const v = boost::get(v_idx, boost::target(e, g));

        return d_edge_t(boost::numeric_cast<d_vertex_idx_t>(u),
                        boost::numeric_cast<d_vertex_idx_t>(v));
      };

    auto e_itr = boost::edges(g);
    boost::compute::copy(
      boost::make_transform_iterator(e_itr.first, make_device_edge),
      boost::make_transform_iterator(e_itr.second, make_device_edge),
      d_edges.begin(), queue);
  }

public:
  template <typename Graph, typename PositionMap>
  void
  setup(Graph const& g,
        PositionMap pos,
        d_edge_size_t const global_work_size_1d,
        d_edge_size_t const local_work_size_1d = 16) {
    n_crossings = 0;

    d_crossing_count = boost::compute::vector<d_edge_size_t>(
      global_work_size_1d * global_work_size_1d, context);
    boost::compute::fill(d_crossing_count.begin(), d_crossing_count.end(),
                         d_edge_size_t(0), queue);

    setup_device_vertex_pos(g, pos);

    setup_device_edges(g);

    // work_offset_i (3) and work_offset_j (4) are set in
    // edge_crossings_kernel::compute_batch
    kernel.set_arg(0, d_pos);
    kernel.set_arg(1, d_edges);
    kernel.set_arg(2, boost::num_edges(g));
    kernel.set_arg(5, global_work_size_1d);
    kernel.set_arg(6, d_crossing_count);

    global_work_size = {global_work_size_1d, global_work_size_1d};
    local_work_size = {local_work_size_1d, local_work_size_1d};
  }

  void
  compute_batch(d_edge_size_t const work_offset_i,
                d_edge_size_t const work_offset_j) {
    kernel.set_arg(3, work_offset_i);
    kernel.set_arg(4, work_offset_j);

    queue.enqueue_nd_range_kernel(
      kernel, boost::compute::dim(0, 0),
      boost::compute::dim(global_work_size[0], global_work_size[1]),
      boost::compute::dim(local_work_size[0], local_work_size[1]));
  }

  void
  finish() {
    boost::compute::reduce(d_crossing_count.begin(), d_crossing_count.end(),
                           &n_crossings, boost::compute::plus<d_edge_size_t>(),
                           queue);
  }
};

std::vector<boost::compute::device>
get_compute_devices() {
  std::vector<boost::compute::device> devices;
  for (auto const& platform : boost::compute::system::platforms()) {
    // edge_crossings_kernel is designed for GPUs and Accelerators.
    // For CPUs, sweep line algorithms should be used for checking line
    // segment intersections.
    for (auto& device :
         platform.devices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_ACCELERATOR)) {
      devices.emplace_back(device);
    }
  }
  return devices;
}

}  // namespace detail

/**
 * Computes the number of edge crossings of a graph layout.
 *
 * This implementation is primarily designed for handling large graphs, with
 * testing conducted on datasets up to several hundred million edges.
 *
 * Counting edge crossings equates to counting line segment intersections. While
 * conventional O(n log n) sweep line algorithms exist for this purpose, I
 * found that employing a massively parallelized pairwise segment intersection
 * approach, as demonstrated in this implementation with a complexity of O(n^2),
 * often outperforms serial O(n log n) implementations for large graphs.
 *
 * Leveraging OpenCL (via Boost.Compute), this implementation efficiently
 * computes pairwise segment intersections in parallel. Moreover, it supports
 * multiple GPUs or accelerator, further boosting its performance capabilities.
 *
 * @param g Graph
 * @param pos A vertex property map for vertex positions
 * @param max_global_work_size_1d The max global work size per batch
 *                                (default: 4096)
 * @param local_work_size_1d The local work size (default: 16)
 * @return the number of edge crossings
 */
template <
  typename Graph,
  typename PositionMap,
  typename edge_size_t = typename boost::graph_traits<Graph>::edges_size_type>
edge_size_t
num_edge_crossings(Graph const& g,
                   PositionMap pos,
                   edge_size_t const max_global_work_size_1d = 4096,
                   edge_size_t const local_work_size_1d = 16) {
  static std::map<cl_device_id, detail::edge_crossings_kernel> kernels;
  std::vector<boost::compute::device> devices{detail::get_compute_devices()};

  if (devices.size() == 0) {
    throw std::runtime_error("Unable to find OpenCL devices");
  }

  // Init a kernel for each device. Reuse a kernel that was initialized in a
  // previous execution.
  for (auto const& device : devices) {
    auto kernel_itr = kernels.find(device.id());
    if (kernel_itr == kernels.end()) {
      kernels.emplace(device.id(), detail::edge_crossings_kernel(device));
    }
  }

  edge_size_t const n_edges = boost::num_edges(g);

  // Determine global work size per batch
  edge_size_t const global_work_size_1d =
    std::min(((n_edges / local_work_size_1d) + 1) * local_work_size_1d,
             max_global_work_size_1d);

  edge_size_t const n_groups_1d = (n_edges / global_work_size_1d) + 1;

  // Construct batches
  std::vector<std::pair<edge_size_t, edge_size_t>> batches;
  for (edge_size_t group_i = 0; group_i < n_groups_1d; group_i++) {
    for (edge_size_t group_j = group_i; group_j < n_groups_1d; group_j++) {
      batches.emplace_back(group_i, group_j);
    }
  }

  std::size_t const n_devices = devices.size();
  std::atomic_size_t current_batch_idx(0);

  tbb::task_group task_group;
  for (std::size_t device_idx = 0; device_idx < n_devices; device_idx++) {
    task_group.run([&, device_idx] {
      auto kernel_itr = kernels.find(devices[device_idx].id());
      if (kernel_itr == kernels.end()) {
        return;
      }

      auto& kernel = kernel_itr->second;
      kernel.setup(g, pos, global_work_size_1d, local_work_size_1d);

      // Each device processes next batch as soon as it completes the previous
      // one, continue until all batches are processed.
      while (current_batch_idx < batches.size()) {
        auto const& batch = batches[current_batch_idx++];
        kernel.compute_batch(batch.first * global_work_size_1d,
                             batch.second * global_work_size_1d);
      }

      kernel.finish();
    });
  }
  task_group.wait();

  // Sum n_crossings from all devices
  edge_size_t n_crossings = 0;
  for (auto const& device : devices) {
    auto kernel_itr = kernels.find(device.id());
    if (kernel_itr != kernels.end()) {
      n_crossings += kernel_itr->second.n_crossings;
    }
  }

  return n_crossings;
}

namespace detail {

template <
  typename Graph,
  typename edge_size_t = typename boost::graph_traits<Graph>::edges_size_type>
double
crosslessness_impl(Graph const& g, edge_size_t const n_crossings) {
  edge_size_t const n_edges = boost::num_edges(g);
  edge_size_t const n_all = n_edges * (n_edges - 1) / 2;
  edge_size_t n_impossible = 0;
  BOOST_FOREACH (auto u, boost::vertices(g)) {
    auto const d = boost::degree(u, g);
    n_impossible += d * (d - 1);
  }
  n_impossible /= 2;
  auto const n_max = n_all - n_impossible;
  return n_max > 0 ? 1.0 - (double(n_crossings) / n_max) : 1.0;
}

}  // namespace detail

/**
 * Computes the crosslessness metric of a graph layout.
 *
 * The definition of crosslessness metric can be found in
 * https://arxiv.org/pdf/1710.04328 (page 5)
 *
 * @param g Graph
 * @param pos A vertex property map for vertex positions
 * @param max_global_work_size_1d The max global work size per batch
 *                                (default: 4096)
 * @param local_work_size_1d The local work size (default: 16)
 * @return std::pair<crosslessness, the number of edge crossings>
 */
template <
  typename Graph,
  typename PositionMap,
  typename edge_size_t = typename boost::graph_traits<Graph>::edges_size_type>
std::pair<double, edge_size_t>
crosslessness(Graph const& g,
              PositionMap pos,
              edge_size_t const max_global_work_size_1d = 4096,
              edge_size_t const local_work_size_1d = 16) {
  edge_size_t const n_crossings =
    num_edge_crossings(g, pos, max_global_work_size_1d, local_work_size_1d);

  return std::pair<double, edge_size_t>(
    detail::crosslessness_impl(g, n_crossings), n_crossings);
}

}  // namespace kw

#endif  // KW_GRAPH_LAYOUT_METRIC_EDGE_CROSSING_HPP

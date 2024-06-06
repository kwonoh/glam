// Copyright 2017 Oh-Hyun Kwon (kwonoh.net)
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef KW_GRAPH_LAYOUT_METRIC_MINIMUM_ANGLE_HPP
#define KW_GRAPH_LAYOUT_METRIC_MINIMUM_ANGLE_HPP

#include <tbb/tbb.h>

#include <boost/algorithm/clamp.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <cmath>
#include <kw/iterator_range.hpp>

namespace kw {

namespace detail {
double
get_min_angle(double min_angle, std::vector<double> const& angles) {
  if (angles.size() > 1) {
    for (int i = 0, j = 1, n = angles.size(); j < n; i = j, j++) {
      min_angle = std::min(min_angle, angles[j] - angles[i]);
    }
  } else if (angles.size() == 1) {
    min_angle = std::min(min_angle, std::abs(angles[0]));
  }
  return min_angle;
}
}  // namespace detail

template <typename Graph, typename PositionMap>
double
min_angle_metric(Graph const& g, PositionMap pos) {
  typedef typename boost::property_traits<PositionMap>::value_type point_t;

  auto vi = boost::get(boost::vertex_index, g);
  auto const num_vertices = boost::num_vertices(g);
  std::vector<double> avg_dev_min_angles(num_vertices, 0.0);

  tbb::parallel_for_each(
    kw::make_iterator_range(boost::vertices(g)), [&](auto const u) {
      auto const u_degree = boost::out_degree(u, g);
      if (u_degree > 1) {
        double const desired_angle =
          boost::math::constants::two_pi<double>() / u_degree;

        std::vector<point_t> directions;
        directions.reserve(u_degree);

        std::vector<double> negative_angles;
        std::vector<double> positive_angles;
        negative_angles.reserve(u_degree);
        positive_angles.reserve(u_degree);

        auto const& center = pos[u];
        BOOST_FOREACH (auto v, boost::adjacent_vertices(u, g)) {
          auto const& target = pos[v];
          auto direction = target;
          boost::geometry::subtract_point(direction, center);
          directions.push_back(direction);
        }

        auto const& dir_0 = directions[0];
        auto const x_0 = boost::geometry::get<0>(dir_0);
        auto const y_0 = boost::geometry::get<1>(dir_0);
        for (int i = 1; i < u_degree; i++) {
          auto const dir_i = directions[i];
          auto const x_i = boost::geometry::get<0>(dir_i);
          auto const y_i = boost::geometry::get<1>(dir_i);
          auto const theta =
            std::atan2(x_0 * y_i - y_0 * x_i, x_0 * x_i + y_0 * y_i);

          if (theta < 0) {
            negative_angles.push_back(theta);
          } else {
            positive_angles.push_back(theta);
          }
        }

        std::sort(negative_angles.begin(), negative_angles.end());
        std::sort(positive_angles.begin(), positive_angles.end());

        double min_angle = boost::math::constants::pi<double>();
        if (negative_angles.size() > 0 && positive_angles.size() > 0) {
          double min_angle = boost::math::constants::two_pi<double>() +
                             negative_angles[0] - positive_angles.back();
        }

        min_angle = detail::get_min_angle(min_angle, negative_angles);
        min_angle = detail::get_min_angle(min_angle, positive_angles);

        avg_dev_min_angles[vi[u]] =
          std::abs((desired_angle - min_angle) / desired_angle);
      }
    });

  return boost::algorithm::clamp(
    1.0 - std::accumulate(avg_dev_min_angles.begin(), avg_dev_min_angles.end(),
                          0.0) /
            num_vertices,
    0.0, 1.0);
}

}  // namespace kw

#endif  // KW_GRAPH_LAYOUT_METRIC_MINIMUM_ANGLE_HPP
